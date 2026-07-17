#!/usr/bin/env python3
"""Host tests for the profile log parser (tools/parse_profile.py).

Covers spilled_frames, whose reading the README cadence colours are defined
against: spills under 25% of a phase's frames flap 16<->8, at or above it the
phase has slipped a tier. Both only mean that if the count is frames.

Run:  python -m unittest discover -s tools/profile_tests
"""

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import parse_profile as pp   # noqa: E402

W = pp.DISPLAY_WINDOW_US     # 62_500 us = one display window at 480 RPM


def _window(renders=(), wall_sum=None, frames=None):
    """A Window carrying per-frame renders (us), or only a wall sum.

    Wall is render plus a nonzero sync idle: equal wall and render is the
    separate no-*_buffer_wait-scope case that RenderIsWall covers.
    """
    n = frames if frames is not None else len(renders)
    w = pp.Window("Fx", 288, 144, 1, n, 1)
    w.frame_rows = [(i + 1, r + 1_000, r) for i, r in enumerate(renders)]
    if wall_sum is not None:
        w.wall = (0, 0, 0, wall_sum)
    return w


class SpilledFrames(unittest.TestCase):
    def test_frame_under_one_window_does_not_spill(self):
        self.assertEqual(pp.spilled_frames(_window([W - 1, 1, W // 2])), 0)

    def test_frame_over_one_window_spills(self):
        self.assertEqual(pp.spilled_frames(_window([W + 1])), 1)

    def test_frame_at_exactly_one_window_does_not_spill(self):
        # It still makes its flip, so the boundary is exclusive.
        self.assertEqual(pp.spilled_frames(_window([W])), 0)

    def test_heavy_frame_counts_once_not_per_missed_flip(self):
        # 3.2 windows: 3 missed flips, but one spilled frame. Counting flips
        # here is what let the aggregate read 532 spilled of 512 frames.
        self.assertEqual(pp.spilled_frames(_window([W * 3 + W // 5])), 1)

    def test_never_exceeds_frame_count(self):
        renders = [W * 4] * 10
        w = _window(renders)
        self.assertEqual(pp.spilled_frames(w), 10)
        self.assertLessEqual(pp.spilled_frames(w), w.frames)

    def test_wall_fallback_is_bounded_by_frame_count(self):
        # No per-frame rows: 10 frames whose wall sum spans 40 windows. The
        # extra-window estimate (30) exceeds the frames, so it clamps.
        w = _window(wall_sum=40 * W, frames=10)
        self.assertEqual(pp.spilled_frames(w), 10)

    def test_wall_fallback_counts_extra_windows(self):
        # 10 frames, 13 windows of wall: 3 frames took a second window.
        w = _window(wall_sum=13 * W, frames=10)
        self.assertEqual(pp.spilled_frames(w), 3)

    def test_no_wall_and_no_rows_is_zero(self):
        self.assertEqual(pp.spilled_frames(_window(frames=4)), 0)


if __name__ == "__main__":
    unittest.main()


class RenderIsWall(unittest.TestCase):
    """An effect with no *_buffer_wait scope: Profile.ino's wall-minus-wait
    degenerates to wall, so nothing derived from it may be called render."""

    def _w(self, walls, wait_scope):
        n = len(walls)
        w = pp.Window("Fx", 288, 144, 1, n, 1)
        # No wait scope => render == wall on every frame (the degenerate case).
        w.frame_rows = [(i + 1, x, x if not wait_scope else x - 10_000)
                        for i, x in enumerate(walls)]
        w.render = (sum(walls) // n, max(walls))
        w.wall = (min(walls), sum(walls) // n, max(walls), sum(walls))
        w.counters = {"frame": {"us": sum(walls), "calls": n, "cyc": 0, "pct": 100}}
        if wait_scope:
            w.counters["fx_buffer_wait"] = {"us": 10_000 * n, "calls": n,
                                            "cyc": 0, "pct": 10}
        return w

    def test_detects_missing_wait_scope(self):
        self.assertTrue(self._w([125_000] * 4, wait_scope=False).render_is_wall())

    def test_wait_scope_present_is_not_wall(self):
        self.assertFalse(self._w([125_000] * 4, wait_scope=True).render_is_wall())

    def test_peak_is_not_reported_as_render(self):
        # 125 ms wall is two display windows of a ~77 ms render plus idle;
        # reporting 125 as a "peak render" is the bug this guards.
        peak, exact = self._w([125_000] * 4, wait_scope=False).peak_render_ms()
        self.assertIsNone(peak)
        self.assertFalse(exact)

    def test_render_ms_unknown_without_wait_scope(self):
        self.assertIsNone(self._w([125_000] * 4, wait_scope=False).render_ms())

    def test_render_ms_known_with_wait_scope(self):
        self.assertAlmostEqual(
            self._w([125_000] * 4, wait_scope=True).render_ms(), 115.0)

    def test_spill_uses_wall_formula_when_render_is_wall(self):
        # 4 frames, 8 windows of wall => 4 frames each took a second window.
        # The per-frame >62.5ms test would also say 4 here, but it counts
        # jitter above the boundary as a spill; the wall formula does not.
        self.assertEqual(pp.spilled_frames(self._w([125_000] * 4, wait_scope=False)), 4)

    def test_jitter_above_boundary_is_not_a_spill(self):
        # Wall quantizes to whole windows; 63 ms is one window plus jitter.
        w = self._w([63_000] * 8, wait_scope=False)
        self.assertEqual(pp.spilled_frames(w), 0)
