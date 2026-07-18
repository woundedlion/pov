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


def _synth_log(path, shapes, per_window=4):
    """A capture where each shape advances mid-window.

    `shapes` is [(name, F, [render_us, ...])]. The marker is emitted before the
    shape's first frame line, exactly as the device logs it: spawn_shape runs
    during a frame, so its serial line precedes that frame's row and the
    enclosing window's dump -- which is what made the window-level attribution
    credit the outgoing shape's frames to the incoming one.
    """
    lines = ["profile harness: effect=Fx segments=4 rpm=480 f_cpu=600000000"]
    rows, n = [], 0
    for name, F, renders in shapes:
        lines.append(f"Spawning Shape: {name} (V=1, E=1, F={F}, I=1)")
        for r in renders:
            n += 1
            lines.append(f"f {n} w={r + 5_000} r={r}")
            rows.append(r)
            if len(rows) == per_window:
                s, e = n - per_window + 1, n
                lines.append(f"=== profile Fx [288x144] frames {s}-{e} "
                             f"window=1000000 us ===")
                lines.append(f"frame wall us: min=0 avg=0 max=0 sum=0 "
                             f"({per_window} frames)")
                lines.append(f"frame render us: avg={sum(rows)//len(rows)} "
                             f"max={max(rows)}")
                lines.append(f"frame                  1 us (100%)  "
                             f"{per_window} calls  1 cyc")
                rows = []
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class StraddleWindowAttribution(unittest.TestCase):
    """A window spanning a shape advance holds frames of BOTH shapes.

    Guards the bug where a cheap shape following an expensive one inherited the
    expensive one's peak (two byte-identical 182-face solids reported 87.5 and
    98.9 ms peaks -- their predecessors' -- against a true 51.2/52.4 ms).
    """

    def _buckets(self, shapes):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = _synth_log(Path(d) / "cap.log", shapes)
            windows, _ = pp.parse(p)
        out = {}
        for w in windows:
            for f in w.frame_rows:
                out.setdefault(f[3]["name"], []).append(f[2])
        return out

    def test_cheap_shape_after_expensive_keeps_its_own_peak(self):
        # 6 frames each, windows of 4 => the window at frames 5-8 straddles.
        got = self._buckets([("expensive", 1082, [100_000] * 6),
                             ("cheap", 182, [50_000] * 6)])
        self.assertEqual(max(got["cheap"]), 50_000)
        self.assertEqual(max(got["expensive"]), 100_000)

    def test_every_frame_is_attributed_exactly_once(self):
        # 6/6/4 frames over windows of 4: the window at frames 5-8 straddles
        # a->b. Total is a whole number of windows, as a real capture is (a
        # trailing partial window never dumps, so its rows never parse).
        got = self._buckets([("a", 74, [10_000] * 6),
                             ("b", 182, [20_000] * 6),
                             ("c", 542, [30_000] * 4)])
        self.assertEqual(sum(len(v) for v in got.values()), 16)
        self.assertEqual({k: len(v) for k, v in got.items()},
                         {"a": 6, "b": 6, "c": 4})

    def test_first_frame_of_a_shape_belongs_to_that_shape(self):
        # The advance frame pays the new shape's rebuild and draws nothing;
        # it is the incoming shape's cost, not the outgoing one's.
        got = self._buckets([("first", 74, [10_000] * 4),
                             ("second", 182, [5_000] + [40_000] * 5)])
        self.assertIn(5_000, got["second"])
        self.assertNotIn(5_000, got["first"])


class ScanMetricsLines(unittest.TestCase):
    """The HS_SCAN_METRICS 'scan totals' window line (Profile.ino)."""

    SCAN = ("scan totals: tested=800 culled=200 lut=300 convex=200 "
            "sector=80 walk=20 cand=64 backstop=0")

    def _log(self, path, scan_line=SCAN, shade_calls=32):
        lines = [
            "=== profile Fx [288x144] frames 1-4 window=1000000 us ===",
            "frame wall us: min=0 avg=0 max=0 sum=0 (4 frames)",
            "frame render us: avg=0 max=0",
            "frame                  1 us (100%)  4 calls  1 cyc",
            f"  raster_shade         1 us (100%)  {shade_calls} calls  1 cyc",
        ]
        if scan_line:
            lines.append(scan_line)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _parse(self, **kw):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            windows, _ = pp.parse(self._log(Path(d) / "cap.log", **kw))
        return windows

    def test_totals_are_parsed(self):
        w = self._parse()[0]
        self.assertEqual(w.scan["tested"], 800)
        self.assertEqual(w.scan["culled"], 200)
        self.assertEqual(w.scan["walk"], 20)
        self.assertEqual(w.scan["cand"], 64)

    def test_absent_line_leaves_scan_none(self):
        self.assertIsNone(self._parse(scan_line=None)[0].scan)

    def test_scan_line_is_not_read_as_a_counter(self):
        # 'scan totals:' must not fall through to COUNTER_RE and invent a scope.
        w = self._parse()[0]
        self.assertEqual(set(w.counters), {"frame", "raster_shade"})

    def test_row_reports_per_frame_and_path_split(self):
        w = self._parse()[0]
        row = pp._metrics_row(w.scan, w.frames, 32).split()
        # 800/4 probes, 200/4 culled; 600 served splits 50/33.3/13.3/3.3.
        self.assertEqual(row[0], "200")
        self.assertEqual(row[1], "50")
        self.assertEqual(row[2], "50.0")
        self.assertEqual(row[3], "33.3")
        self.assertEqual(row[4], "13.3")
        self.assertEqual(row[5], "3.3")

    def test_probes_per_shade_uses_the_raster_shade_call_count(self):
        w = self._parse()[0]
        # 800 probes over 32 shade events.
        self.assertEqual(pp._metrics_row(w.scan, w.frames, 32).split()[-1],
                         "25.00")

    def test_zero_shade_events_do_not_divide_by_zero(self):
        w = self._parse()[0]
        self.assertEqual(pp._metrics_row(w.scan, w.frames, 0).split()[-1], "nan")

    def test_metrics_command_reports_missing_instrumentation(self):
        self.assertEqual(pp.cmd_metrics(self._parse(scan_line=None)), 2)

    def test_metrics_command_succeeds_on_an_instrumented_capture(self):
        self.assertEqual(pp.cmd_metrics(self._parse()), 0)
