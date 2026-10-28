#!/usr/bin/env python3
"""Host tests for the shared Teensy device lock (tools/device_lock.sh and
tools/device_lock_guard.py).

Run:  python -m unittest discover -s tools/profile_tests
"""

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import isolated_env  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
LOCK_SH = REPO / "tools" / "device_lock.sh"
LOCK_GUARD = REPO / "tools" / "device_lock_guard.py"

GRACE = 120  # HS_DEVICE_STALE_GRACE default


class LockGuardUsageTests(unittest.TestCase):
    def test_unicode_claim_roundtrips_and_releases(self):
        with tempfile.TemporaryDirectory() as directory:
            claim = Path(directory) / "claim"
            token = "claim-雪"
            environment = {**os.environ, "PYTHONIOENCODING": "utf-8"}
            result = subprocess.run(
                [sys.executable, str(LOCK_GUARD), "claim", str(claim)],
                input=f"token={token}\nsession=雪\n", text=True, encoding="utf-8",
                capture_output=True, env=environment, timeout=5)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((claim / "info").read_text(encoding="utf-8"),
                             f"token={token}\nsession=雪\n")
            result = subprocess.run(
                [sys.executable, str(LOCK_GUARD), "break", str(claim), token],
                capture_output=True, text=True, env=environment, timeout=5)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(claim.exists())

    def test_invalid_arguments_report_usage(self):
        for args in ([], ["claim"], ["break", "directory"],
                     ["unknown", "directory", "token"], ["claim", "directory", "extra"]):
            with self.subTest(args=args):
                result = subprocess.run([sys.executable, str(LOCK_GUARD), *args],
                                        capture_output=True, text=True, timeout=5)
                self.assertEqual(result.returncode, 2)
                self.assertIn("usage:", result.stderr)
                self.assertNotIn("Traceback", result.stderr)


class CheckoutBuildLockTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        tools = self.root / "tools"
        tools.mkdir()
        for name in ("device_lock.sh", "device_lock_guard.py", "find_python.sh", "teensy_cold_build.sh"):
            shutil.copyfile(REPO / "tools" / name, tools / name)
        subprocess.run(["git", "-C", str(self.root), "init", "--quiet"], check=True, env=isolated_env())
        self.env = isolated_env()

    def command(self, *args):
        return subprocess.run(["bash", *map(str, args)], cwd=self.root,
                              env=self.env, capture_output=True, text=True, timeout=10)

    def test_busy_build_tree_blocks_cleanup_and_size_builds(self):
        image = self.root / ".pio" / "build" / "firmware.hex"
        image.parent.mkdir(parents=True)
        image.write_text("existing image", encoding="utf-8")
        holder, pid = _live_bash()
        self.addCleanup(_stop_bash, holder, pid)
        lock = self.root / ".profile-lock"
        lock.mkdir()
        (lock / "info").write_text(f"token=holder\npid={pid}\n", encoding="utf-8")
        for script, args in (("teensy_cold_build.sh", []),
                             ("device_lock.sh", ["tree", "bash", "-c", "touch ran"])):
            with self.subTest(script=script):
                result = self.command(self.root / "tools" / script, *args)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("already claimed", result.stderr)
                self.assertEqual(image.read_text(encoding="utf-8"), "existing image")
                self.assertFalse((self.root / "ran").exists())
                self.assertTrue((lock / "info").exists())

    def test_failed_checkout_command_releases_its_claim(self):
        result = self.command(self.root / "tools" / "device_lock.sh", "tree",
                              "bash", "-c", "test -f .profile-lock/info || exit 9; exit 3")
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertFalse((self.root / ".profile-lock").exists())

    def test_failed_cold_build_releases_its_claim(self):
        pio = self.root / "pio"
        pio.write_text("#!/bin/bash\ntest -f .profile-lock/info || exit 9\nexit 3\n",
                       encoding="utf-8", newline="\n")
        pio.chmod(0o755)
        self.env["PATH"] = str(self.root) + os.pathsep + os.environ["PATH"]
        result = self.command(self.root / "tools" / "teensy_cold_build.sh")
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertFalse((self.root / ".profile-lock").exists())


def clean_device_env():
    return {key: value for key, value in os.environ.items()
            if not key.startswith(("HS_DEVICE_", "HS_TEENSY_"))}


def is_stale(lock_dir):
    """Run _hs_lock_is_stale against lock_dir; True = breakable."""
    script = f'. "{LOCK_SH}"; _hs_lock_is_stale "{lock_dir}"'
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                       env=clean_device_env())
    return r.returncode == 0


def break_lock(lock_dir, expected="stale", prelude=""):
    """Run _hs_break_lock against lock_dir; True = we won the right to evict."""
    script = (f'. "{LOCK_SH}"; {prelude} '
              f'_hs_break_lock "{lock_dir}" "{expected}"')
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                       env=clean_device_env())
    return r.returncode == 0


def _live_bash():
    """A running bash and its own $$, as the lock stores it.

    The pid must come from bash: acquire writes `pid=$$`, and under Git Bash
    that is an MSYS pid its `kill -0` can resolve, while a native Windows pid
    (os.getpid()) reads back as dead and would fake a stale lock here.
    """
    p = subprocess.Popen(["bash", "-c", "echo $$; exec sleep 30"],
                         stdout=subprocess.PIPE, text=True)
    return p, int(p.stdout.readline())


def _stop_bash(process, pid):
    """Terminate and reap the fixture's reported MSYS/POSIX owner."""
    subprocess.run(["bash", "-c", 'kill -TERM "$1" 2>/dev/null || :',
                    "fixture-stop", str(pid)], check=True, timeout=5)
    process.communicate(timeout=5)
    subprocess.run([
        "bash", "-c",
        'for ((i=0; i<100; ++i)); do '
        'kill -0 "$1" 2>/dev/null || exit 0; sleep 0.01; done; exit 1',
        "fixture-wait", str(pid)], check=True, timeout=5)


def _live_pid(test_case):
    p, pid = _live_bash()
    test_case.addCleanup(_stop_bash, p, pid)
    return pid


def _dead_pid():
    p, pid = _live_bash()
    _stop_bash(p, pid)
    return pid


class LockStaleness(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp()) / "lock.d"
        self.d.mkdir(parents=True)
        self.addCleanup(shutil.rmtree, self.d.parent, ignore_errors=True)

    def _write_info(self, pid, started, deadline):
        (self.d / "info").write_text(
            f"pid={pid}\nstarted={started}\ndeadline={deadline}\n")

    def _age_dir(self, seconds):
        old = time.time() - seconds
        os.utime(self.d, (old, old))

    def test_claim_being_written_is_not_stale(self):
        # acquire() mkdirs the lock, then writes info. A peer reading in that
        # gap sees no deadline; calling it stale hands two sessions the board.
        self.assertFalse(is_stale(self.d))

    def test_unwritten_claim_past_grace_is_stale(self):
        # Nobody ever wrote info: a crash between mkdir and write, not a race.
        self._age_dir(GRACE + 60)
        self.assertTrue(is_stale(self.d))

    def test_live_holder_within_eta_is_not_stale(self):
        now = int(time.time())
        self._write_info(_live_pid(self), now, now + 600)
        self.assertFalse(is_stale(self.d))

    def test_live_holder_past_eta_and_grace_is_not_stale(self):
        now = int(time.time())
        self._write_info(_live_pid(self), now - 900, now - GRACE - 60)
        self.assertFalse(is_stale(self.d))

    def test_incomplete_pidless_claim_past_eta_within_grace_is_not_stale(self):
        now = int(time.time())
        (self.d / "info").write_text(
            f"started={now - 900}\ndeadline={now - 10}\n", encoding="utf-8")
        self.assertFalse(is_stale(self.d))

    def test_dead_holder_is_stale(self):
        now = int(time.time())
        self._write_info(_dead_pid(), now - 300, now + 600)
        self.assertTrue(is_stale(self.d))

    def test_dead_holder_claimed_seconds_ago_is_not_stale(self):
        # A dead holder with a complete claim stays within its 60 s start grace.
        now = int(time.time())
        self._write_info(_dead_pid(), now, now + 600)
        self.assertFalse(is_stale(self.d))


class LockBreak(unittest.TestCase):
    """Evicting a stale claim. Two peers can judge one claim stale at the same
    moment, so the eviction itself has to pick a single winner."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.d = self.root / "lock.d"

    def _claim(self, token):
        self.d.mkdir(parents=True, exist_ok=True)
        (self.d / "info").write_text(f"token={token}\n")

    def test_breaks_the_claim_it_judged_stale(self):
        self._claim("stale")
        self.assertTrue(break_lock(self.d))
        self.assertFalse(self.d.exists())

    def test_loser_of_the_race_breaks_nothing(self):
        # A peer already evicted this claim. Falling through to delete whatever
        # sits at the path would take out the winner's own fresh lock.
        self.assertFalse(break_lock(self.d))
        self.assertFalse(self.d.exists())

    def test_claim_retaken_since_the_judgement_survives(self):
        self._claim("fresh")
        self.assertFalse(break_lock(self.d, expected="stale"))
        self.assertEqual((self.d / "info").read_text(), "token=fresh\n")

    def test_three_sessions_cannot_claim_during_stale_retirement(self):
        self._claim("stale")
        script = """
import sys
sys.path.insert(0, sys.argv[1])
import device_lock_guard as lock
original = lock.read_token
def paused_read(directory):
    token = original(directory)
    print("checked", flush=True)
    input()
    return token
lock.read_token = paused_read
sys.exit(0 if lock.update_claim(sys.argv[2], "break", "stale") else 1)
"""
        retire = subprocess.Popen(
            [sys.executable, "-c", script, str(LOCK_SH.parent), str(self.d)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
        peers = []
        wrapper = self.root / "guard_wrapper.py"
        wrapper.write_text(
            "import runpy, sys\n"
            'print("guarding", file=sys.stderr, flush=True)\n'
            f"runpy.run_path({str(LOCK_GUARD)!r}, run_name='__main__')\n",
            encoding="utf-8")
        try:
            self.assertEqual(retire.stdout.readline().strip(), "checked")
            for name in ("B", "C"):
                recovery = '_hs_break_lock "$D" stale || :; ' if name == "B" else ""
                body = (f'. "{LOCK_SH}"; _HS_LOCK_HELPER={shlex.quote(wrapper.as_posix())}; '
                        f'D="{self.d}"; echo started; '
                        f'{recovery}_hs_try_claim "$D" COM3 {name} test 60; '
                        'rc=$?; echo "CLAIM=$rc"; exit "$rc"')
                peer = subprocess.Popen(
                    ["bash", "-c", body], stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                peers.append(peer)
                self.assertEqual(peer.stdout.readline().strip(), "started")
            for peer in peers:
                self.assertEqual(peer.stderr.readline().strip(), "guarding")
            time.sleep(0.2)
            self.assertEqual((self.d / "info").read_text(), "token=stale\n")
            self.assertTrue(all(peer.poll() is None for peer in peers))
            retire.communicate("continue\n", timeout=10)
            self.assertEqual(retire.returncode, 0)
            results = [peer.stdout.readline().strip() for peer in peers]
            self.assertEqual(sorted(results), ["CLAIM=0", "CLAIM=1"])
            winner = "B" if results[0] == "CLAIM=0" else "C"
            self.assertIn(f"effect={winner}\n", (self.d / "info").read_text())
        finally:
            if retire.poll() is None:
                retire.kill()
            retire.communicate(timeout=10)
            for peer in peers:
                peer.communicate(timeout=10)

    def test_break_leaves_no_scratch_directory_behind(self):
        self._claim("stale")
        self.assertTrue(break_lock(self.d))
        self.assertEqual([p.name for p in self.root.iterdir()], ["lock.d.guard"])

    def test_guard_is_released_when_its_process_dies(self):
        self._claim("stale")
        script = """
import sys
sys.path.insert(0, sys.argv[1])
from device_lock_guard import guard
with guard(sys.argv[2]):
    print("locked", flush=True)
    input()
"""
        holder = subprocess.Popen(
            [sys.executable, "-c", script, str(LOCK_SH.parent), str(self.d)],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True)
        try:
            self.assertEqual(holder.stdout.readline().strip(), "locked")
        finally:
            holder.kill()
            holder.communicate(timeout=10)
        self.assertTrue(break_lock(self.d))


def run_lock(script, lock_base, ports=("COM3", "COM4"), env=None, pnp=()):
    """Run a snippet against device_lock.sh with a stubbed board list.

    hs_device_ports is stubbed after sourcing; the real enumerator needs
    attached boards. Get-PnpDevice always is: pnp holds the live
    "<InstanceId>|<FriendlyName>" rows it reports.
    """
    stub = "_hs_pnp_teensys() { printf '%%s\\n' %s; };" % (
        " ".join(shlex.quote(row) for row in pnp) or "''")
    if ports is not None:
        body = "; ".join(f"echo {p}" for p in ports) or ":"
        stub += " hs_device_ports() { %s; };" % body
    full = f'. "{LOCK_SH}"; {stub} {script}'
    e = dict(clean_device_env(), HS_DEVICE_LOCK=str(lock_base))
    e.update(env or {})
    return subprocess.run(["bash", "-c", full], capture_output=True, text=True,
                          env=e)


class BoardSelection(unittest.TestCase):
    """Several boards attached: one lock each, and acquire finds a free one."""

    def setUp(self):
        self.base = Path(tempfile.mkdtemp()) / "lock"
        self.addCleanup(shutil.rmtree, self.base.parent, ignore_errors=True)

    def test_force_claims_first_busy_board(self):
        self.hold("COM3")
        self.hold("COM4")
        result = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                          self.base, env={"HS_DEVICE_FORCE": "1"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PORT=COM3", result.stdout)
        self.assertNotIn("token=peer", (self.lock_dir("COM3") / "info").read_text())

    def test_force_preserves_young_tokenless_claim(self):
        self.lock_dir("COM3").mkdir()
        self.hold("COM4")
        result = run_lock("hs_device_acquire E profile 60", self.base,
                          env={"HS_DEVICE_FORCE": "1"})
        self.assertEqual(result.returncode, 1)
        self.assertTrue(self.lock_dir("COM3").is_dir())
        self.assertFalse((self.lock_dir("COM3") / "info").exists())

    def test_force_preserves_mid_break_retake(self):
        self.hold("COM3")
        self.hold("COM4")
        script = ('_hs_break_stale() { echo "token=peerB" > "$1/info"; '
                  '_hs_break_lock "$1" "$2"; }; hs_device_acquire E profile 60')
        result = run_lock(script, self.base, env={"HS_DEVICE_FORCE": "1"})
        self.assertEqual(result.returncode, 1)
        self.assertEqual((self.lock_dir("COM3") / "info").read_text().strip(),
                         "token=peerB")

    def test_no_enumerated_board_fails_without_a_claim(self):
        result = run_lock("hs_device_acquire E profile 60", self.base, ports=())
        self.assertEqual(result.returncode, 1)
        self.assertIn("no Teensy is enumerated", result.stderr)
        self.assertEqual(list(self.base.parent.glob("lock*.d")), [])

    def test_wait_reenumerates_a_board_after_it_appears(self):
        ready = self.base.parent / "ready"
        script = (f'hs_device_ports() {{ test ! -f "{ready}" || echo COM3; }}; '
                  f'sleep() {{ touch "{ready}"; }}; '
                  'hs_device_acquire E profile 60; rc=$?; '
                  'echo "port=$HS_DEVICE_PORT"; hs_device_release; exit "$rc"')
        result = run_lock(script, self.base, ports=(), env={"HS_DEVICE_WAIT": "5"})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("port=COM3", result.stdout)
        self.assertEqual(list(self.base.parent.glob("lock*.d")), [])

    def lock_dir(self, port):
        return Path(f"{self.base}-{port}.d")

    def hold(self, port, deadline_in=600, pid=None):
        d = self.lock_dir(port)
        d.mkdir(parents=True)
        now = int(time.time())
        (d / "info").write_text(
            f"token=peer\nsession=peer\npid={pid or _live_pid(self)}\n"
            f"port={port}\neffect=Peer\nenv=profile\nstarted={now}\n"
            f"deadline={now + deadline_in}\n")

    def test_acquires_a_board_and_pins_the_port(self):
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM3", r.stdout)
        self.assertTrue(self.lock_dir("COM3").is_dir())
        self.assertFalse(self.lock_dir("COM4").exists())

    def test_busy_board_is_skipped_for_a_free_one(self):
        self.hold("COM3")
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM4", r.stdout)

    def test_all_boards_busy_fails_fast(self):
        self.hold("COM3")
        self.hold("COM4")
        r = run_lock("hs_device_acquire E profile 60", self.base)
        self.assertEqual(r.returncode, 1)
        self.assertIn("ALL DEVICES BUSY", r.stderr)

    def test_free_board_is_preferred_over_breaking_a_stale_lock(self):
        # A stale lock is never broken while another board sits free.
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=_dead_pid())
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM4", r.stdout)
        self.assertNotIn("breaking", r.stderr)

    def test_stale_lock_is_broken_when_no_board_is_free(self):
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=_dead_pid())
        self.hold("COM4")
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM3", r.stdout)
        self.assertIn("stale", r.stderr)

    def test_aged_unwritten_claim_is_recovered(self):
        directory = self.lock_dir("COM3")
        directory.mkdir()
        old = time.time() - GRACE - 60
        os.utime(directory, (old, old))
        self.hold("COM4")
        result = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                          self.base)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("PORT=COM3", result.stdout)
        self.assertTrue((directory / "info").is_file())

    def test_a_claim_retaken_mid_break_is_not_evicted(self):
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=_dead_pid())
        self.hold("COM4")
        script = (
            '_hs_break_stale() { echo "token=peerB" > "$1/info"; '
            '_hs_break_lock "$1" "$2"; }; '
            'hs_device_acquire E profile 60')
        r = run_lock(script, self.base)
        self.assertEqual(r.returncode, 1)
        self.assertIn("token=peerB",
                      (self.lock_dir("COM3") / "info").read_text())

    def test_pinned_port_does_not_wander_to_another_board(self):
        # HS_TEENSY_PORT names the board under test; falling back to a free
        # peer board would profile the wrong hardware silently.
        self.hold("COM3")
        tools = self.base.parent / "stub-tools"
        tools.mkdir()
        enumerator = tools / "teensy_ports.exe"
        enumerator.write_text(
            "#!/bin/bash\nprintf '%s\\n' 'usb-a COM3 (Teensy 4.0)' 'usb-b COM4 (Teensy 4.0)'\n",
            encoding="utf-8", newline="\n")
        enumerator.chmod(0o755)
        r = run_lock("hs_device_acquire E profile 60", self.base, ports=None,
                     env={"HS_TEENSY_PORT": "COM3", "HS_TEENSY_TOOLS": str(tools)})
        self.assertEqual(r.returncode, 1)

        self.assertIn("ALL DEVICES BUSY", r.stderr)
        self.assertFalse(self.lock_dir("COM4").exists())

    def test_claim_whose_info_never_landed_is_not_handed_out(self):
        # The stubbed info write reads back with no token, a claim
        # hs_device_release could never match.
        script = ('_hs_lock_field() { :; }; hs_device_acquire E profile 60; '
                  'echo "RC=$?"; echo "TOKEN=[$_HS_TOKEN]"; '
                  'echo "PIN=[$HS_TEENSY_PORT]"')
        r = run_lock(script, self.base)
        self.assertIn("RC=1", r.stdout)
        self.assertIn("TOKEN=[]", r.stdout)
        self.assertIn("PIN=[]", r.stdout)
        self.assertIn("cannot record the claim", r.stderr)
        self.assertFalse(self.lock_dir("COM3").exists())
        self.assertFalse(self.lock_dir("COM4").exists())

    def test_release_frees_only_our_own_claim(self):
        r = run_lock("hs_device_acquire E profile 60; hs_device_release", self.base)
        self.assertEqual(r.returncode, 0)
        self.assertFalse(self.lock_dir("COM3").exists())

    def test_release_drops_the_pin_so_the_next_acquire_can_roam(self):
        # acquire -> release -> acquire in one shell: a kept pin would steer the
        # second claim onto the freed board.
        script = ('hs_device_ports() { if [ -n "${HS_TEENSY_PORT:-}" ]; then '
                  'echo "$HS_TEENSY_PORT"; else echo COM3; echo COM4; fi; }; '
                  'hs_device_acquire E profile 60; hs_device_release; '
                  'echo "PIN=[${HS_TEENSY_PORT:-}]"; '
                  f'mkdir "{self.base}-COM3.d"; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"')
        r = run_lock(script, self.base, ports=None)
        self.assertIn("PIN=[]", r.stdout)
        self.assertIn("PORT=COM4", r.stdout)

    def test_release_keeps_a_pin_the_caller_set(self):
        script = ('hs_device_acquire E profile 60; hs_device_release; '
                  'echo "PIN=[${HS_TEENSY_PORT:-}]"')
        r = run_lock(script, self.base, ports=("COM3",),
                     env={"HS_TEENSY_PORT": "COM3"})
        self.assertIn("PIN=[COM3]", r.stdout)

    def test_release_restores_a_peers_claim_it_declined_to_free(self):
        # A declined release preserves the peer's claim and leaves no scratch directory.
        script = ('hs_device_acquire E profile 60; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_release; echo "RC=$?"')
        r = run_lock(script, self.base)
        self.assertIn("RC=0", r.stdout)
        self.assertEqual((self.lock_dir("COM3") / "info").read_text().strip(),
                         "token=peer")
        self.assertEqual(sorted(x.name for x in self.base.parent.iterdir()),
                         ["lock-COM3.d", "lock-COM3.d.guard"])

    def test_wait_queues_under_errexit(self):
        # Under set -e, an unclaimable board still reaches the wait loop and guidance.
        self.hold("COM3")
        self.hold("COM4")
        script = ('set -e; sleep() { :; }; hs_device_acquire E profile 60; '
                  'echo UNREACHED')
        r = run_lock(script, self.base, env={"HS_DEVICE_WAIT": "5"})
        self.assertEqual(r.returncode, 1)
        self.assertIn("waiting up to 5s", r.stderr)
        self.assertIn("Every attached Teensy is in use", r.stderr)
        self.assertNotIn("UNREACHED", r.stdout)

    def test_fail_fast_reports_the_guidance_under_errexit(self):
        self.hold("COM3")
        self.hold("COM4")
        script = 'set -e; hs_device_acquire E profile 60; echo UNREACHED'
        r = run_lock(script, self.base)
        self.assertEqual(r.returncode, 1)
        self.assertIn("Every attached Teensy is in use", r.stderr)
        self.assertNotIn("UNREACHED", r.stdout)

    def test_release_survives_errexit(self):
        # Under a caller's `set -e`, a declined release must not abort
        # mid-teardown.
        script = ('set -e; hs_device_acquire E profile 60; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_release; echo DONE')
        r = run_lock(script, self.base)
        self.assertIn("DONE", r.stdout)
        self.assertEqual(r.returncode, 0)


class MissingLockRoot(unittest.TestCase):
    """A lock base whose parent directory does not exist.

    A missing parent must be reported as a claim failure, not a busy device.
    """

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.base = self.root / "absent" / "lock"

    def test_the_guard_names_the_missing_root(self):
        r = subprocess.run(
            [sys.executable, str(LOCK_GUARD), "claim", f"{self.base}-COM3.d"],
            input="token=t\n", capture_output=True, text=True,
            encoding="utf-8")
        self.assertEqual(r.returncode, 2)
        self.assertIn("lock root", r.stderr)

    def test_acquire_reports_the_path_not_a_busy_bench(self):
        r = run_lock("hs_device_acquire E profile 60", self.base)
        self.assertEqual(r.returncode, 2)
        self.assertNotIn("ALL DEVICES BUSY", r.stderr)
        self.assertIn("lock root", r.stderr)

    def test_tree_acquire_reports_guard_failure(self):
        result = run_lock(f'TREE="{self.root / "absent"}"; acquire_tree_lock', self.base)
        self.assertEqual(result.returncode, 2)
        self.assertIn("lock guard could not run", result.stderr)
        self.assertNotIn("already claimed", result.stderr)

    def test_an_existing_root_still_acquires(self):
        self.base.parent.mkdir(parents=True)
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM3", r.stdout)


class PinnedPortEnumeration(unittest.TestCase):
    """A caller-set HS_TEENSY_PORT is checked against the loader's board list.

    A board replugged onto a new COM name leaves the old pin naming nothing.
    """

    def setUp(self):
        self.base = Path(tempfile.mkdtemp()) / "lock"
        self.addCleanup(shutil.rmtree, self.base.parent, ignore_errors=True)
        self.tools = self.base.parent / "tool-teensy"
        self.tools.mkdir(parents=True)
        loader = self.tools / "teensy_ports.exe"
        loader.write_text("#!/bin/bash\necho '1 COM7 Teensy'\necho '2 COM9 Teensy'\n")
        loader.chmod(0o755)

    def run_pinned(self, script, pin):
        return run_lock(script, self.base, ports=None,
                        env={"HS_TEENSY_TOOLS": str(self.tools),
                             "HS_TEENSY_PORT": pin})

    def test_empty_enumeration_rejects_a_pinned_port(self):
        (self.tools / "teensy_ports.exe").write_text("#!/bin/bash\nexit 0\n")
        result = self.run_pinned("hs_device_ports", "COM3")
        self.assertEqual(result.returncode, 1)
        self.assertIn("not attached", result.stderr)

    def test_attached_pin_is_returned(self):
        r = self.run_pinned("hs_device_ports", "COM9")
        self.assertEqual(r.returncode, 0)
        self.assertEqual(r.stdout.split(), ["COM9"])

    def test_unattached_pin_is_reported(self):
        r = self.run_pinned("hs_device_ports", "COM3")
        self.assertEqual(r.returncode, 1)
        self.assertIn("not attached", r.stderr)

    def hold(self, port):
        d = Path(f"{self.base}-{port}.d")
        d.mkdir()
        now = int(time.time())
        (d / "info").write_text(
            f"token=peer\nsession=peer\npid={_live_pid(self)}\n"
            f"port={port}\neffect=Peer\nenv=profile\nstarted={now}\n"
            f"deadline={now + 600}\n")

    def test_status_distinguishes_an_unattached_pin_from_busy(self):
        result = self.run_pinned("hs_device_status", "COM3")
        self.assertEqual(result.returncode, 2)
        self.hold("COM7")
        result = self.run_unpinned("hs_device_status")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("COM7 BUSY", result.stdout)
        self.assertIn("COM9 free", result.stdout)
        self.hold("COM9")
        result = self.run_unpinned("hs_device_status")
        self.assertEqual(result.returncode, 1, result.stdout)
        self.assertIn("COM9 BUSY", result.stdout)

    def test_acquire_refuses_an_unattached_pin_without_locking(self):
        r = self.run_pinned("hs_device_acquire E profile 60", "COM3")
        self.assertEqual(r.returncode, 1)
        self.assertEqual(list(self.base.parent.glob("lock-*.d")), [])

    def fail_loader(self):
        """A loader that is present but whose enumeration fails."""
        loader = self.tools / "teensy_ports.exe"
        loader.write_text("#!/bin/bash\necho 'boom' >&2\nexit 3\n")
        loader.chmod(0o755)

    def run_unpinned(self, script):
        return run_lock(script, self.base, ports=None,
                        env={"HS_TEENSY_TOOLS": str(self.tools)})

    def test_missing_loader_reports_configuration_error(self):
        (self.tools / "teensy_ports.exe").unlink()
        r = self.run_unpinned("hs_device_ports")
        self.assertEqual(r.returncode, 2)
        self.assertIn("not found (set HS_TEENSY_TOOLS)", r.stderr)

    def test_failed_enumeration_reports_loader_failure(self):
        self.fail_loader()
        r = self.run_unpinned("hs_device_ports")
        self.assertEqual(r.returncode, 2)
        self.assertEqual(r.stdout.strip(), "")
        self.assertIn("failed", r.stderr)

    def test_acquire_refuses_a_failed_enumeration_without_locking(self):
        self.fail_loader()
        r = self.run_unpinned("hs_device_acquire E profile 60")
        self.assertEqual(r.returncode, 2)
        self.assertEqual(list(self.base.parent.glob("lock-*.d")), [])

    def test_status_reports_a_failed_enumeration(self):
        self.fail_loader()
        r = self.run_unpinned("hs_device_status")
        self.assertEqual(r.returncode, 2)


# Get-PnpDevice rows for a live Teensy in HalfKay and for a running one.
HALFKAY = r"USB\VID_16C0&PID_0478\0009B551|USB Input Device"


def serial_row(port):
    return rf"USB\VID_16C0&PID_0483&MI_00\6&15795FF9&0&0000|USB Serial Device ({port})"


class PinnedPortReenumeration(unittest.TestCase):
    """A queued waiter's pinned board dropping off the port list."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root, ignore_errors=True)
        self.base = self.root / "lock"
        self.tools = self.root / "tool-teensy"
        self.tools.mkdir()
        self.peer = Path(f"{self.base}-COM4.d")

    def enumerate_com4_except(self, absent_calls):
        """teensy_ports.exe omits COM4 on the 1-based calls in absent_calls."""
        calls = (self.root / "calls").as_posix()
        absent = " ".join(map(str, absent_calls)) if absent_calls != "all" else "*"
        loader = self.tools / "teensy_ports.exe"
        loader.write_text(
            "#!/bin/bash\n"
            f'n=$(( $(cat "{calls}" 2>/dev/null || echo 0) + 1 )); echo "$n" > "{calls}"\n'
            "echo 'usb-a COM3 (Teensy 4.0)'\n"
            + ('' if absent == "*" else
               f'case " {absent} " in *" $n "*) ;; *) echo "usb-b COM4 (Teensy 4.0)";; esac\n'),
            encoding="utf-8", newline="\n")
        loader.chmod(0o755)

    def hold_com4(self):
        self.peer.mkdir()
        now = int(time.time())
        (self.peer / "info").write_text(
            f"token=peer\nsession=peer\npid={_live_pid(self)}\nport=COM4\n"
            f"effect=Peer\nenv=profile\nstarted={now}\ndeadline={now + 600}\n")

    def acquire(self, release_after=None, pnp=(), **env):
        """Acquire pinned to COM4; the peer releases on sleep #release_after."""
        release = (f'[ "$s" -lt {release_after} ] || rm -rf "{self.peer.as_posix()}"; '
                   if release_after else "")
        script = (f's=0; sleep() {{ s=$((s + 1)); {release}}}; '
                  'hs_device_acquire E profile 60; rc=$?; '
                  'echo "RC=$rc SLEEPS=$s PORT=$HS_DEVICE_PORT"; hs_device_release; exit "$rc"')
        return run_lock(script, self.base, ports=None, pnp=pnp,
                        env={"HS_TEENSY_TOOLS": str(self.tools),
                             "HS_TEENSY_PORT": "COM4", **env})

    def test_waiter_outlasts_a_board_reenumerating_mid_wait(self):
        self.hold_com4()
        self.enumerate_com4_except((3, 4, 5))
        r = self.acquire(release_after=4, HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=4 PORT=COM4", r.stdout)
        self.assertIn("dropped off the port list", r.stderr)
        self.assertNotIn("not attached", r.stderr)

    def test_waiter_queued_while_a_claimed_board_is_off_the_list(self):
        self.hold_com4()
        self.enumerate_com4_except((1, 2))
        r = self.acquire(release_after=3, HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=3 PORT=COM4", r.stdout)

    def test_board_seen_by_the_waiter_may_drop_after_its_claim_is_released(self):
        self.hold_com4()
        self.enumerate_com4_except((3, 4))
        r = self.acquire(release_after=1, HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=3 PORT=COM4", r.stdout)

    def test_grace_counts_one_continuous_absence(self):
        self.hold_com4()
        self.enumerate_com4_except((1, 2, 4, 5))
        r = self.acquire(release_after=5, HS_DEVICE_WAIT="3600",
                         HS_DEVICE_ABSENT_GRACE="10")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=5 PORT=COM4", r.stdout)

    def test_known_board_that_never_returns_fails_after_the_grace(self):
        self.hold_com4()
        self.enumerate_com4_except("all")
        r = self.acquire(HS_DEVICE_WAIT="3600", HS_DEVICE_ABSENT_GRACE="15")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=3 ", r.stdout)
        self.assertIn("stayed off the port list for 15s", r.stderr)
        self.assertIn("token=peer", (self.peer / "info").read_text())

    def test_absent_board_fails_when_the_wait_runs_out_first(self):
        self.hold_com4()
        self.enumerate_com4_except("all")
        r = self.acquire(HS_DEVICE_WAIT="10")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=2 ", r.stdout)

    def test_never_seen_port_fails_fast_despite_a_wait(self):
        self.enumerate_com4_except("all")
        r = self.acquire(HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=0 ", r.stdout)
        self.assertIn("COM4 is not attached", r.stderr)
        self.assertEqual(list(self.root.glob("lock-*.d")), [])

    def test_unseen_pin_is_waited_for_while_a_teensy_is_in_halfkay(self):
        self.enumerate_com4_except((1, 2, 3))
        r = self.acquire(pnp=(HALFKAY,), HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=3 PORT=COM4", r.stdout)
        self.assertIn("shows a Teensy in HalfKay", r.stderr)
        self.assertNotIn("not attached", r.stderr)

    def test_unseen_pin_is_waited_for_while_windows_shows_its_port_live(self):
        self.enumerate_com4_except((1,))
        r = self.acquire(pnp=(serial_row("COM4"),), HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("RC=0 SLEEPS=1 PORT=COM4", r.stdout)
        self.assertIn("shows it live", r.stderr)

    def test_board_stuck_in_halfkay_fails_after_the_grace(self):
        self.enumerate_com4_except("all")
        r = self.acquire(pnp=(HALFKAY,), HS_DEVICE_WAIT="3600",
                         HS_DEVICE_ABSENT_GRACE="15")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=3 ", r.stdout)
        self.assertIn("stayed off the port list for 15s", r.stderr)
        self.assertEqual(list(self.root.glob("lock-*.d")), [])

    def test_halfkay_wait_is_bounded_by_a_shorter_device_wait(self):
        self.enumerate_com4_except("all")
        r = self.acquire(pnp=(HALFKAY,), HS_DEVICE_WAIT="10")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=2 ", r.stdout)

    def test_other_live_ports_are_no_reason_to_wait(self):
        self.enumerate_com4_except("all")
        r = self.acquire(pnp=(serial_row("COM3"), serial_row("COM44")),
                         HS_DEVICE_WAIT="3600")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=0 ", r.stdout)
        self.assertIn("COM4 is not attached", r.stderr)

    def test_halfkay_is_not_consulted_without_a_wait(self):
        self.enumerate_com4_except("all")
        r = self.acquire(pnp=(HALFKAY,))
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=0 ", r.stdout)
        self.assertIn("COM4 is not attached", r.stderr)

    def test_claimed_absent_board_fails_fast_without_a_wait(self):
        self.hold_com4()
        self.enumerate_com4_except("all")
        r = self.acquire()
        self.assertEqual(r.returncode, 1)
        self.assertIn("RC=1 SLEEPS=0 ", r.stdout)
        self.assertIn("COM4 is not attached", r.stderr)

    def test_absent_budget_is_zero_without_a_wait_and_capped_by_the_grace(self):
        for env, expected in (({}, "0"), ({"HS_DEVICE_WAIT": "10"}, "10"),
                              ({"HS_DEVICE_WAIT": "3600"}, "120"),
                              ({"HS_DEVICE_WAIT": "3600", "HS_DEVICE_ABSENT_GRACE": "45"}, "45")):
            with self.subTest(env=env):
                r = run_lock("hs_device_absent_budget", self.base, ports=None, env=env)
                self.assertEqual(r.stdout.strip(), expected, r.stderr)


if __name__ == "__main__":
    unittest.main()
