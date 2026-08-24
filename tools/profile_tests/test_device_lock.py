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
import tempfile
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LOCK_SH = REPO / "tools" / "device_lock.sh"

GRACE = 120  # HS_DEVICE_STALE_GRACE default


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
        self.addCleanup(p.kill)
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
        p.kill()
        p.wait()
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
        # The holder released and a peer re-claimed between our staleness
        # judgement and the break; that claim is live and must be put back.
        self._claim("fresh")
        mark = self.root / "seen"
        prelude = (
            'MARK="%s"; _hs_lock_field() { '
            'if [ "$2" = token ] && [ ! -e "$MARK" ]; then : > "$MARK"; '
            'echo stale; return; fi; '
            'sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1; };' % mark)
        self.assertFalse(break_lock(self.d, prelude=prelude))
        self.assertTrue(self.d.is_dir())
        self.assertIn("token=fresh", (self.d / "info").read_text())

    def test_a_peer_claiming_during_the_put_back_is_not_nested_into(self):
        # The window between "is the slot free?" and the restoring move: a peer
        # that claims it there must get its lock back intact, with no scratch
        # directory left inside it for the next debugger to find.
        self._claim("fresh")
        mark = self.root / "seen"
        prelude = (
            'MARK="%s"; D="%s"; MVN=0; '
            'mv() { MVN=$((MVN+1)); '
            '[ "$MVN" = 2 ] && { mkdir -p "$D"; echo token=peer > "$D/info"; }; '
            'command mv "$@"; }; '
            '_hs_lock_field() { '
            'if [ "$2" = token ] && [ ! -e "$MARK" ]; then : > "$MARK"; '
            'echo stale; return; fi; '
            'sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1; };'
            % (mark, self.d))
        self.assertFalse(break_lock(self.d, prelude=prelude))
        self.assertTrue(self.d.is_dir())
        self.assertEqual([p.name for p in self.d.iterdir()], ["info"])
        self.assertEqual(sorted(p.name for p in self.root.iterdir()),
                         ["lock.d", "seen"])

    def test_break_leaves_no_scratch_directory_behind(self):
        self._claim("stale")
        self.assertTrue(break_lock(self.d))
        self.assertEqual(list(self.root.iterdir()), [])


def run_lock(script, lock_base, ports=("COM3", "COM4"), env=None):
    """Run a snippet against device_lock.sh with a stubbed board list.

    hs_device_ports is overridden after sourcing rather than faking
    teensy_ports.exe: the selection logic is what these test, and the real
    enumerator needs boards physically attached.
    """
    stub = ""
    if ports is not None:
        body = "; ".join(f"echo {p}" for p in ports)
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
        self.addCleanup(p.kill)
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

    def test_release_leaves_a_lock_reclaimed_by_a_peer(self):
        # Our claim was broken as stale and re-taken; our teardown must not
        # unlock the board out from under whoever holds it now.
        script = ('hs_device_acquire E profile 60; '
                  f'echo token=peer > "{self.base}-COM3.d/info"; '
                  'hs_device_release')
        run_lock(script, self.base)
        self.assertTrue(self.lock_dir("COM3").is_dir())

    @staticmethod
    def _dead_pid():
        p, pid = _live_bash()
        p.kill()
        p.wait()
        return pid


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

    def test_attached_pin_is_returned(self):
        r = self.run_pinned("hs_device_ports", "COM9")
        self.assertEqual(r.returncode, 0)
        self.assertEqual(r.stdout.split(), ["COM9"])

    def test_unattached_pin_is_reported(self):
        r = self.run_pinned("hs_device_ports", "COM3")
        self.assertEqual(r.returncode, 1)
        self.assertIn("not attached", r.stderr)

    def test_acquire_refuses_an_unattached_pin_without_locking(self):
        r = self.run_pinned("hs_device_acquire E profile 60", "COM3")
        self.assertEqual(r.returncode, 1)
        self.assertEqual(list(self.base.parent.glob("lock-*.d")), [])


if __name__ == "__main__":
    unittest.main()
