#!/usr/bin/env python3
"""Host tests for the shared Teensy device lock (tools/device_lock.sh).

One bench board is shared by every concurrent session, and the failure mode is
silent both ways: an evicted holder's capture gets spliced across two firmware
images, and the evictor can capture the peer's firmware under its own effect
name. So staleness must never fire on a claim whose owner is alive -- these
drive _hs_lock_is_stale through bash directly.

Run:  python -m unittest discover -s tools/profile_tests
"""

import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LOCK_SH = REPO / "tools" / "device_lock.sh"
LOCK_GUARD = REPO / "tools" / "device_lock_guard.py"

GRACE = 120  # HS_DEVICE_STALE_GRACE default


class LockGuardUsageTests(unittest.TestCase):
    def test_invalid_arguments_report_usage(self):
        for args in ([], ["claim"], ["break", "directory"],
                     ["unknown", "directory", "token"], ["claim", "directory", "extra"]):
            with self.subTest(args=args):
                result = subprocess.run([sys.executable, str(LOCK_GUARD), *args],
                                        capture_output=True, text=True, timeout=5)
                self.assertEqual(result.returncode, 2)
                self.assertIn("usage:", result.stderr)
                self.assertNotIn("Traceback", result.stderr)


def is_stale(lock_dir):
    """Run _hs_lock_is_stale against lock_dir; True = breakable."""
    script = f'. "{LOCK_SH}"; _hs_lock_is_stale "{lock_dir}"'
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    return r.returncode == 0


def break_lock(lock_dir, expected="stale", prelude=""):
    """Run _hs_break_lock against lock_dir; True = we won the right to evict."""
    script = (f'. "{LOCK_SH}"; {prelude} '
              f'_hs_break_lock "{lock_dir}" "{expected}"')
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
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


class LockStaleness(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp()) / "lock.d"
        self.d.mkdir(parents=True)
        self.addCleanup(shutil.rmtree, self.d.parent, ignore_errors=True)

    def _write_info(self, pid, started, deadline):
        (self.d / "info").write_text(
            f"pid={pid}\nstarted={started}\ndeadline={deadline}\n")

    def _live_pid(self):
        p, pid = _live_bash()
        self.addCleanup(_stop_bash, p, pid)
        return pid

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
        self._write_info(self._live_pid(), now, now + 600)
        self.assertFalse(is_stale(self.d))

    def test_holder_past_eta_and_grace_is_stale(self):
        now = int(time.time())
        self._write_info(self._live_pid(), now - 900, now - GRACE - 60)
        self.assertTrue(is_stale(self.d))

    def test_holder_just_past_eta_within_grace_is_not_stale(self):
        # A long capture that overruns its own estimate still owns the board.
        now = int(time.time())
        self._write_info(self._live_pid(), now - 900, now - 10)
        self.assertFalse(is_stale(self.d))

    def test_dead_holder_is_stale(self):
        now = int(time.time())
        self._write_info(self._dead_pid(), now - 300, now + 600)
        self.assertTrue(is_stale(self.d))

    def test_dead_holder_claimed_seconds_ago_is_not_stale(self):
        # The PID check waits out a 60 s window so a peer that has claimed the
        # lock but not yet forked its build is not mistaken for a corpse.
        now = int(time.time())
        self._write_info(self._dead_pid(), now, now + 600)
        self.assertFalse(is_stale(self.d))

    @staticmethod
    def _dead_pid():
        p, pid = _live_bash()
        _stop_bash(p, pid)
        return pid


class LockBreak(unittest.TestCase):
    """Evicting a stale claim. Two peers can judge one claim stale at the same
    moment, so the eviction itself has to pick a single winner -- otherwise the
    loser deletes the winner's fresh lock and both flash the same board."""

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
        try:
            self.assertEqual(retire.stdout.readline().strip(), "checked")
            for name in ("B", "C"):
                recovery = '_hs_break_lock "$D" stale || :; ' if name == "B" else ""
                body = (f'. "{LOCK_SH}"; D="{self.d}"; echo started; '
                        f'{recovery}_hs_try_claim "$D" COM3 {name} test 60; '
                        'rc=$?; echo "CLAIM=$rc"; read -r done; exit "$rc"')
                peer = subprocess.Popen(
                    ["bash", "-c", body], stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                peers.append(peer)
                self.assertEqual(peer.stdout.readline().strip(), "started")
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
                peer.communicate("done\n", timeout=10)

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


def run_lock(script, lock_base, ports=("COM3", "COM4"), env=None):
    """Run a snippet against device_lock.sh with a stubbed board list.

    hs_device_ports is overridden after sourcing rather than faking
    teensy_ports.exe: the selection logic is what these test, and the real
    enumerator needs boards physically attached.
    """
    stub = ""
    if ports is not None:
        body = "; ".join(f"echo {p}" for p in ports) or ":"
        stub = "hs_device_ports() { %s; };" % body
    full = f'. "{LOCK_SH}"; {stub} {script}'
    e = dict(os.environ, HS_DEVICE_LOCK=str(lock_base))
    e.pop("HS_TEENSY_PORT", None)
    e.update(env or {})
    return subprocess.run(["bash", "-c", full], capture_output=True, text=True,
                          env=e)


class BoardSelection(unittest.TestCase):
    """Several boards attached: one lock each, and acquire finds a free one."""

    def setUp(self):
        self.base = Path(tempfile.mkdtemp()) / "lock"
        self.addCleanup(shutil.rmtree, self.base.parent, ignore_errors=True)

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
            f"token=peer\nsession=peer\npid={pid or self._live_pid()}\n"
            f"port={port}\neffect=Peer\nenv=profile\nstarted={now}\n"
            f"deadline={now + deadline_in}\n")

    def _live_pid(self):
        p, pid = _live_bash()
        self.addCleanup(_stop_bash, p, pid)
        return pid

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
        # Breaking a claim is a last resort; a stale lock on one board must
        # never be taken while another board sits free.
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=self._dead_pid())
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM4", r.stdout)
        self.assertNotIn("breaking", r.stderr)

    def test_stale_lock_is_broken_when_no_board_is_free(self):
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=self._dead_pid())
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
        # Two peers judge one claim stale at once: the first evicts it and takes
        # the board, and the second must not then delete that fresh lock. The
        # stub replays the interleaving by re-claiming during our own read.
        self.hold("COM3", deadline_in=-(GRACE + 60), pid=self._dead_pid())
        self.hold("COM4")
        script = (
            'MARK="%s"; _hs_lock_field() { '
            'if [ "$2" = token ] && [ ! -e "$MARK" ]; then : > "$MARK"; '
            'sed -n "s/^token=//p" "$1/info" | head -1; '
            'echo "token=peerB" > "$1/info"; return; fi; '
            'sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1; }; '
            'hs_device_acquire E profile 60' % (self.base.parent / "seen"))
        r = run_lock(script, self.base)
        self.assertEqual(r.returncode, 1)
        self.assertIn("token=peerB",
                      (self.lock_dir("COM3") / "info").read_text())

    def test_pinned_port_does_not_wander_to_another_board(self):
        # HS_TEENSY_PORT names the board under test; falling back to a free
        # peer board would profile the wrong hardware silently.
        self.hold("COM3")
        r = run_lock("hs_device_acquire E profile 60", self.base, ports=None,
                     env={"HS_TEENSY_PORT": "COM3"})
        self.assertEqual(r.returncode, 1)

    def test_claim_whose_info_never_landed_is_not_handed_out(self):
        # An unreadable claim is one hs_device_release can never match, so the
        # lock dir would hold the board until the stale grace expired. The stub
        # replays that: the info write reads back with no token.
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
        # acquire -> release -> acquire in one shell: the pin the first claim
        # exported would otherwise steer the second onto the freed board, and
        # a peer holding it by then reads as every board busy.
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

    def test_release_leaves_a_lock_reclaimed_by_a_peer(self):
        # Our claim was broken as stale and re-taken; our teardown must not
        # unlock the board out from under whoever holds it now.
        script = ('hs_device_acquire E profile 60; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_release')
        run_lock(script, self.base)
        self.assertTrue(self.lock_dir("COM3").is_dir())

    def test_release_restores_a_peers_claim_it_declined_to_free(self):
        # The declined release goes through the same rename as a break, so the
        # peer's info must be put back and no scratch directory left behind.
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
        # hs_device_status returns 1 when no board is claimable, which is the
        # case every time the wait branch is reached. Under `set -e` that
        # aborted the caller where it should have queued, so neither the loop
        # nor the guidance that follows it ran.
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
        # profile_one.sh runs under `set -e`; a declined release must not abort
        # the caller mid-teardown.
        script = ('set -e; hs_device_acquire E profile 60; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_release; echo DONE')
        r = run_lock(script, self.base)
        self.assertIn("DONE", r.stdout)
        self.assertEqual(r.returncode, 0)

    @staticmethod
    def _dead_pid():
        p, pid = _live_bash()
        _stop_bash(p, pid)
        return pid


class MissingLockRoot(unittest.TestCase):
    """A lock base whose parent directory does not exist.

    No claim can be recorded there, and a claim that cannot be recorded looks
    exactly like a board somebody else holds -- acquire once reported "ALL
    DEVICES BUSY" over a status listing every board as free.
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

    def test_an_existing_root_still_acquires(self):
        self.base.parent.mkdir(parents=True)
        r = run_lock('hs_device_acquire E profile 60 && echo "PORT=$HS_TEENSY_PORT"',
                     self.base)
        self.assertIn("PORT=COM3", r.stdout)


class PinnedPortEnumeration(unittest.TestCase):
    """A caller-set HS_TEENSY_PORT is checked against the loader's board list.

    A board replugged onto a new COM name leaves the old pin naming nothing, and
    locking a port no board answers on only surfaces after two image builds.
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

    def test_status_distinguishes_an_unattached_pin_from_busy(self):
        result = self.run_pinned("hs_device_status", "COM3")
        self.assertEqual(result.returncode, 2)

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

    def test_failed_enumeration_is_distinct_from_no_loader(self):
        """rc 2, not the rc 0 + empty output that means an enumerate-less host."""
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


if __name__ == "__main__":
    unittest.main()
