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


if __name__ == "__main__":
    unittest.main()
