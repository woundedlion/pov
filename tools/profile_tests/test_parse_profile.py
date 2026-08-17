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


class ProbeBreakdownLines(unittest.TestCase):
    """The HS_PROBE_BREAKDOWN 'probe cycles'/'probe counts' lines."""

    # 100 probes; tick = 100 read pairs at 10 cyc => read = 5 cyc/read.
    CYC = ("probe cycles: point=1500 project=3000 lut=0 convex=500 "
           "sector=2000 exact=4000 pack=1200 alpha=400 tick=1000")
    CNT = ("probe counts: probe=100 cull_cos=4 cull_r=6 lut=0 convex=20 "
           "sector=30 exact=40 alpha=50")

    def _log(self, path, cyc=CYC, cnt=CNT):
        lines = [
            "=== profile Fx [288x144] frames 1-4 window=1000000 us ===",
            "frame wall us: min=0 avg=0 max=0 sum=0 (4 frames)",
            "frame render us: avg=0 max=0",
            "frame                  1 us (100%)  4 calls  1 cyc",
        ]
        lines += [x for x in (cyc, cnt) if x]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _parse(self, **kw):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            windows, _ = pp.parse(self._log(Path(d) / "cap.log", **kw))
        return windows

    def test_both_lines_merge_into_one_dict(self):
        w = self._parse()[0]
        self.assertEqual(w.probe["exact"], 4000)
        self.assertEqual(w.probe["n_exact"], 40)
        self.assertEqual(w.probe["tick"], 1000)

    def test_absent_lines_leave_probe_none(self):
        self.assertIsNone(self._parse(cyc=None, cnt=None)[0].probe)

    def test_probe_lines_are_not_read_as_counters(self):
        self.assertEqual(set(self._parse()[0].counters), {"frame"})

    def test_command_reports_net_of_the_measured_read_cost(self):
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(pp.cmd_probe(self._parse()), 0)
        out = buf.getvalue()
        # tick 1000 over 2*100 reads = 5.0 cyc per counter read.
        self.assertIn("counter read=5.0 cyc", out)
        # exact: 4000/40 = 100.0 cyc/event, 95.0 net, 40/100 probes => 38.0.
        row = next(l for l in out.splitlines() if l.startswith("exact"))
        self.assertEqual(row.split()[1:5], ["40", "100.0", "95.0", "38.0"])

    def test_command_exits_2_without_the_flag(self):
        import contextlib
        import io
        with contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(pp.cmd_probe(self._parse(cyc=None, cnt=None)), 2)
        self.assertIn("HS_PROBE_BREAKDOWN", err.getvalue())


class PlotCountLines(unittest.TestCase):
    """The HS_PLOT_COUNTS count-only workload line."""

    PLOT = ("plot counts: r=8,e=40,p=16,g=12,d=2,o=6,t=40,c=4,"
            "s=200,y=160,u=240,a=80,n=320,h=120,x=120,k=11,b=1")

    def _log(self, path, plot_line=PLOT):
        lines = [
            "=== profile Fx [288x144] frames 1-4 window=1000000 us ===",
            "frame wall us: min=0 avg=0 max=0 sum=0 (4 frames)",
            "frame render us: avg=0 max=0",
            "frame                  1 us (100%)  4 calls  1 cyc",
        ]
        if plot_line:
            lines.append(plot_line)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _parse(self, **kw):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            windows, _ = pp.parse(self._log(Path(d) / "cap.log", **kw))
        return windows

    def test_counts_are_parsed(self):
        p = self._parse()[0].plot
        self.assertEqual(p["rings"], 8)
        self.assertEqual(p["planar_unprojects"], 240)
        self.assertEqual(p["steps_peak"], 11)

    def test_absent_line_leaves_plot_none(self):
        self.assertIsNone(self._parse(plot_line=None)[0].plot)

    def test_plot_line_is_not_read_as_a_counter(self):
        self.assertEqual(set(self._parse()[0].counters), {"frame"})

    def test_row_reports_per_frame_and_per_edge_counts(self):
        row = pp._plot_row(self._parse()[0].plot, 4).split()
        self.assertEqual(row[:2], ["2.0", "10.0"])
        self.assertEqual(row[8:13], ["5.0", "4.0", "6.0", "2.0", "8.0"])
        self.assertEqual(row[-4:], ["30.0", "30.0", "11.0", "1.0"])

    def test_command_reports_missing_instrumentation(self):
        self.assertEqual(pp.cmd_plot(self._parse(plot_line=None)), 2)

    def test_command_succeeds_on_count_capture(self):
        self.assertEqual(pp.cmd_plot(self._parse()), 0)

    def test_aggregate_uses_peak_cache_depth(self):
        import contextlib
        import io
        windows = self._parse()
        second = self._parse()[0]
        second.plot["steps_peak"] = 7
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(pp.cmd_plot(windows + [second]), 0)
        aggregate = out.getvalue().splitlines()[-1]
        self.assertTrue(aggregate.endswith("11.0      2.0  8 frames"))


class MindSplatterInstrumentationLines(unittest.TestCase):
    """Dedicated MindSplatter count and short-batch stall lines."""

    LOG = "\n".join([
        "=== profile MindSplatter [288x144] frames 1-4 window=1000000 us ===",
        "frame wall us: min=0 avg=0 max=0 sum=0 (4 frames)",
        "frame render us: avg=0 max=0",
        "frame                  1 us (100%)  4 calls  1 cyc",
        "msp counts particles: resident=400 live=360 full=200 partial=160 draining=40",
        "msp counts gate: cart_lat=10 cart_mer=20 cart_fallback=30 row=4 col=5 edge=6 visible=7 exact=8",
        "msp counts render: dot=90 long=10 adaptive=50 shader=100 pal_end=2 pal_lerp=98 hole_early=70",
        "msp counts aa: tap0=1 tap1=2 tap2=3 tap3=4 tap4=90 interior=90 boundary=10",
        "msp stall: stage=history_vertex batches=20 cyc=2000 cpi=200 lsu=100 exc=2",
        "msp stall: stage=aa_weights batches=100 cyc=4000 cpi=300 lsu=50 exc=1",
    ]) + "\n"

    def _parse(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capture.log"
            path.write_text(self.LOG, encoding="utf-8")
            return pp.parse(path)[0]

    def test_count_sections_merge_and_do_not_become_cycle_counters(self):
        window = self._parse()[0]
        self.assertEqual(window.msp_counts["resident"], 400)
        self.assertEqual(window.msp_counts["exact_fallback"], 8)
        self.assertEqual(window.msp_counts["tap4"], 90)
        self.assertEqual(set(window.counters), {"frame"})

    def test_stall_stages_parse(self):
        window = self._parse()[0]
        self.assertEqual(window.msp_stalls["history_vertex"]["cyc"], 2000)
        self.assertEqual(window.msp_stalls["aa_weights"]["batches"], 100)

    def test_commands_accept_instrumented_capture(self):
        import contextlib
        import io
        windows = self._parse()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(pp.cmd_msp_counts(windows), 0)
            self.assertEqual(pp.cmd_msp_stalls(windows), 0)


class ValidateRequiresData(unittest.TestCase):
    """`validate` is the mandatory pre-trust step: it must not certify nothing."""

    HEADERS = "\n".join(
        f"=== profile Fx [288x144] frames {1 + 10 * i}-{10 + 10 * i} "
        f"window=62500 us ===" for i in range(3))

    @staticmethod
    def _measurable(cyc, wall_sum):
        return "\n".join([
            "=== profile Fx [288x144] frames 1-10 window=62500 us ===",
            f"frame wall us: min=100 avg=100 max=100 sum={wall_sum} (10 frames)",
            f"frame                 {wall_sum} us (100%)   10 calls   {cyc} cyc",
            "=== profile Fx [288x144] frames 11-20 window=62500 us ===",
            f"frame wall us: min=100 avg=100 max=100 sum={wall_sum} (10 frames)",
            f"frame                 {wall_sum} us (100%)   10 calls   {cyc} cyc",
            "=== profile Fx [288x144] frames 21-30 window=62500 us ===",
            f"frame wall us: min=100 avg=100 max=100 sum={wall_sum} (10 frames)",
            f"frame                 {wall_sum} us (100%)   10 calls   {cyc} cyc",
        ])

    def _validate(self, text, scope="frame"):
        import contextlib
        import io
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "prof.log"
            log.write_text(text, encoding="utf-8")
            windows, effect = pp.parse(str(log))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            ok = pp.cmd_validate(windows, effect, scope)
        return ok, out.getvalue()

    def test_bare_headers_are_not_valid(self):
        ok, out = self._validate(self.HEADERS)
        self.assertFalse(ok)
        self.assertIn("INVALID", out)

    def test_measurable_capture_is_valid(self):
        # 600_000 cyc / 600 == 1000 us, matching the wall sum exactly.
        ok, out = self._validate(self._measurable(600_000, 1000))
        self.assertTrue(ok)
        self.assertIn("VALID", out)

    def test_ppm_drift_is_still_caught(self):
        ok, _ = self._validate(self._measurable(700_000, 1000))
        self.assertFalse(ok)

    def test_five_ppm_boundary(self):
        ok, _ = self._validate(self._measurable(600_002_999, 1_000_000))
        self.assertTrue(ok)
        ok, _ = self._validate(self._measurable(600_003_001, 1_000_000))
        self.assertFalse(ok)

    def test_scope_selects_the_window_to_validate(self):
        text = "\n".join([
            "=== profile Fx [288x144] frames 1-10 window=62500 us ===",
            "frame wall us: min=100 avg=100 max=100 sum=2000 (10 frames)",
            "frame                 2000 us (100%)   10 calls   1200000 cyc",
            "xx_render              100 us (5%)     10 calls   60000 cyc",
            "=== profile Fx [288x144] frames 11-20 window=62500 us ===",
            "frame wall us: min=100 avg=100 max=100 sum=1000 (10 frames)",
            "frame                 1000 us (100%)   10 calls   700000 cyc",
            "xx_render              900 us (90%)    10 calls   540000 cyc",
            "=== profile Fx [288x144] frames 21-30 window=62500 us ===",
            "frame wall us: min=100 avg=100 max=100 sum=1500 (10 frames)",
            "frame                 1500 us (100%)   10 calls   900000 cyc",
            "xx_render              200 us (13%)    10 calls   120000 cyc",
        ])
        ok, out = self._validate(text, "xx_render")
        self.assertFalse(ok)
        self.assertIn("frames 11-20", out)


class PullbackTelemetryValidation(unittest.TestCase):
    MANIFEST_PATH = (Path(__file__).resolve().parents[2] /
                     "tests/data/pullback/programs.json")

    @classmethod
    def setUpClass(cls):
        cls.manifest = pp.load_shaderball_program_manifest(cls.MANIFEST_PATH)
        cls.by_preset = {}
        for program in cls.manifest["programs"]:
            for preset in program["presets"]:
                cls.by_preset[preset] = program["id"]

    def _telemetry(self):
        events = []
        preset_count = len(self.by_preset)
        for preset in range(preset_count):
            next_preset = (preset + 1) % preset_count
            events.extend([
                {"preset": preset, "total": preset_count,
                 "pipeline": self.by_preset[preset], "endpoint": "steady"},
                {"preset": preset, "total": preset_count,
                 "pipeline": self.by_preset[preset], "endpoint": "from"},
                {"preset": next_preset, "total": preset_count,
                 "pipeline": self.by_preset[next_preset], "endpoint": "to"},
            ])
        return {
            "arms": [{"arm": "LANDED", "sha": self.manifest["capture_sha"][:12]}],
            "programs": events,
        }

    def _validate(self, telemetry, expected="LANDED", manifest=None):
        checks = []

        def check(condition, message):
            checks.append((condition, message))

        pp.validate_pullback_telemetry(
            telemetry, expected, manifest or self.manifest, check)
        return all(condition for condition, _message in checks), checks

    def test_complete_synthetic_cycle(self):
        ok, _ = self._validate(self._telemetry())
        self.assertTrue(ok)

    def test_arm_and_program_records_parse(self):
        import tempfile
        text = "\n".join([
            f"Pullback arm: LANDED sha={self.manifest['capture_sha'][:12]}",
            "Pullback program: preset=0/13 "
            f"pipeline={self.by_preset[0]} endpoint=steady",
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capture.log"
            path.write_text(text, encoding="utf-8")
            _windows, _effect, telemetry = pp.parse_capture(path)
        self.assertEqual(telemetry["arms"][0]["arm"], "LANDED")
        self.assertEqual(telemetry["programs"][0]["pipeline"],
                         self.by_preset[0])

    def test_arm_rejections(self):
        telemetry = self._telemetry()
        telemetry["arms"] = []
        self.assertFalse(self._validate(telemetry)[0])
        telemetry = self._telemetry()
        telemetry["arms"].append(dict(telemetry["arms"][0]))
        self.assertFalse(self._validate(telemetry)[0])
        telemetry = self._telemetry()
        telemetry["arms"][0]["arm"] = "CORE"
        self.assertFalse(self._validate(telemetry)[0])

    def test_sha_mismatch_is_rejected(self):
        telemetry = self._telemetry()
        telemetry["arms"][0]["sha"] = "0" * 12
        self.assertFalse(self._validate(telemetry)[0])

    def test_unknown_and_mismatched_programs_are_rejected(self):
        telemetry = self._telemetry()
        telemetry["programs"][0]["pipeline"] = "UNKNOWN_PROGRAM"
        self.assertFalse(self._validate(telemetry)[0])
        telemetry = self._telemetry()
        telemetry["programs"][0]["pipeline"] = self.by_preset[1]
        self.assertFalse(self._validate(telemetry)[0])

    def test_missing_steady_is_rejected(self):
        telemetry = self._telemetry()
        telemetry["programs"] = [
            event for event in telemetry["programs"]
            if not (event["preset"] == 5 and event["endpoint"] == "steady")]
        self.assertFalse(self._validate(telemetry)[0])

    def test_missing_transition_halves_are_rejected(self):
        for endpoint in ("from", "to"):
            telemetry = self._telemetry()
            removed = False
            kept = []
            for event in telemetry["programs"]:
                if not removed and event["endpoint"] == endpoint:
                    removed = True
                    continue
                kept.append(event)
            telemetry["programs"] = kept
            self.assertFalse(self._validate(telemetry)[0])

    def test_dynamic_none_is_rejected(self):
        telemetry = self._telemetry()
        telemetry["programs"][0]["pipeline"] = "NONE"
        self.assertFalse(self._validate(telemetry)[0])

    def test_duplicate_tuple_is_rejected(self):
        telemetry = self._telemetry()
        telemetry["programs"].insert(1, dict(telemetry["programs"][0]))
        self.assertFalse(self._validate(telemetry)[0])

    def test_missing_last_to_first_wrap_is_rejected(self):
        telemetry = self._telemetry()
        telemetry["programs"] = telemetry["programs"][:-1]
        self.assertFalse(self._validate(telemetry)[0])


class FramesMode(unittest.TestCase):
    """frame_rows carry a fourth preset-owner field; `frames` must unpack it.

    The spill column counts frames, like spilled_frames: a 190 ms frame misses
    three flips but is one spilled frame.
    """

    LOG = "\n".join([
        "f 1 w=70000 r=64000",
        "f 2 w=60000 r=55000",
        "f 3 w=200000 r=190000",
        "=== profile Fx [288x144] frames 1-3 window=62500 us ===",
        "frame wall us: min=60000 avg=110000 max=200000 sum=330000 (3 frames)",
        "frame                 330000 us (100%)   3 calls   198000000 cyc",
    ])

    def test_frames_mode_prints_a_row_per_frame(self):
        import contextlib
        import io
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "prof.log"
            log.write_text(self.LOG, encoding="utf-8")
            argv = sys.argv
            sys.argv = ["parse_profile.py", str(log), "frames"]
            try:
                with contextlib.redirect_stdout(io.StringIO()) as out:
                    self.assertEqual(pp.main(), 0)
            finally:
                sys.argv = argv
        rows = [l for l in out.getvalue().splitlines() if not l.startswith("#")]
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0].split(), ["1", "70000", "64000", "1"])
        self.assertEqual(rows[1].split(), ["2", "60000", "55000", "0"])
        self.assertEqual(rows[2].split(), ["3", "200000", "190000", "1"])
        self.assertEqual(sum(int(r.split()[3]) for r in rows),
                         pp.spilled_frames(_window([64000, 55000, 190000])))


class TruncatedTrailingWindow(unittest.TestCase):
    """A capture cut mid-dump leaves a header whose stats never arrived.

    Read as a window, its frames enter the spill denominator unmeasured and its
    absent per-frame render marks every OTHER window's exact peak a '~'
    placeholder -- one severed serial line downgrading a whole run.
    """

    def _log(self, path, windows=3, truncated=False):
        lines = []
        n = 0
        for _ in range(windows):
            for _ in range(4):
                n += 1
                lines.append(f"f {n} w=55000 r=50000")
            lines.append(f"=== profile Fx [288x144] frames {n - 3}-{n} "
                         f"window=250000 us ===")
            lines.append("frame wall us: min=55000 avg=55000 max=55000 "
                         "sum=220000 (4 frames)")
            lines.append("frame render us: avg=50000 max=50000")
            lines.append("frame             220000 us (100%)  4 calls  1 cyc")
            lines.append("  fx_buffer_wait   20000 us (9%)  4 calls  1 cyc")
        if truncated:
            for _ in range(4):
                n += 1
                lines.append(f"f {n} w=55000 r=50000")
            lines.append(f"=== profile Fx [288x144] frames {n - 3}-{n} "
                         f"window=250000 us ===")
            lines.append("frame wall us: min=55000 avg=550")  # cut mid-line
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    def _parse(self, **kw):
        import contextlib
        import io
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            path = self._log(Path(d) / "cap.log", **kw)
            with contextlib.redirect_stderr(io.StringIO()) as err:
                windows, _ = pp.parse(path)
        return windows, err.getvalue()

    def _windows_view(self, **kw):
        import contextlib
        import io
        windows, _ = self._parse(**kw)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            pp.cmd_windows(windows, "frame")
        return out.getvalue()

    def test_the_cut_window_is_dropped_and_reported(self):
        windows, err = self._parse(truncated=True)
        self.assertEqual(len(windows), 3)
        self.assertIn("13-16", err)

    def test_a_complete_capture_keeps_every_window(self):
        windows, err = self._parse()
        self.assertEqual(len(windows), 3)
        self.assertEqual(err, "")

    def test_the_run_stays_exact_and_the_denominator_is_measured_frames(self):
        self.assertEqual(self._windows_view(truncated=True),
                         self._windows_view())
        self.assertIn("peak frame render", self._windows_view(truncated=True))
        self.assertIn("spilled 0/12 frames", self._windows_view(truncated=True))

    def test_a_capture_that_never_dumps_counters_keeps_its_last_window(self):
        # Only a window that differs from its peers is a cut dump; a build
        # whose windows all carry no counter tree is not truncated.
        import contextlib
        import io
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "cap.log"
            path.write_text("\n".join(
                f"=== profile Fx [288x144] frames {i * 4 + 1}-{i * 4 + 4} "
                f"window=250000 us ===" for i in range(3)) + "\n",
                encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()) as err:
                parsed, _ = pp.parse(path)
        self.assertEqual(len(parsed), 3)
        self.assertEqual(err.getvalue(), "")


class BucketOrdering(unittest.TestCase):
    """A preset's clean-hold scope time cannot exceed its own peak render.

    The bucket holds every frame the preset was on screen for, transitions
    included, so it is a superset of the frames the clean holds average over.
    Two committed profile rows once reported the inversion.
    """

    def _run(self, render_us, shader_us_per_frame, frames=4):
        import contextlib
        import io
        import tempfile
        lines = ["Spawning Shape: solid (V=1, E=1, F=74, I=1)"]
        for n in range(1, frames + 1):
            lines.append(f"f {n} w={render_us + 5_000} r={render_us}")
        lines.append(f"=== profile Fx [288x144] frames 1-{frames} "
                     f"window=1000000 us ===")
        lines.append(f"frame wall us: min=0 avg=0 max=0 sum=0 ({frames} frames)")
        lines.append(f"frame render us: avg={render_us} max={render_us}")
        lines.append(f"frame            {render_us * frames} us (100%)  "
                     f"{frames} calls  1 cyc")
        lines.append(f"  shader         {shader_us_per_frame * frames} us (50%)  "
                     f"{frames} calls  1 cyc")
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "cap.log"
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            windows, _ = pp.parse(path)
        with contextlib.redirect_stderr(io.StringIO()) as err:
            with contextlib.redirect_stdout(io.StringIO()) as out:
                status = pp.cmd_buckets(windows, "shader", None)
        return status, out.getvalue(), err.getvalue()

    def test_a_clean_hold_above_its_peak_render_fails(self):
        status, _, err = self._run(34_940, 39_650)
        self.assertEqual(status, 1)
        self.assertIn("exceeds peak render", err)
        self.assertIn("39.65", err)
        self.assertIn("34.94", err)

    def test_an_ordered_capture_passes_and_reports_both(self):
        status, out, err = self._run(49_520, 43_470)
        self.assertEqual(status, 0)
        self.assertEqual(err, "")
        row = [l for l in out.splitlines() if not l.startswith("#")][0]
        self.assertIn("49.52", row)
        self.assertIn("43.47", row)


class CleanHoldSelection(unittest.TestCase):
    def test_transition_window_with_one_missing_call_is_excluded(self):
        marker = {"key": "preset", "idx": 12, "total": 17, "name": "12"}
        held = pp.Window("Fx", 288, 144, 1, 16, 1)
        held.marker = marker
        held.counters["shader"] = {
            "us": 16 * 50_000, "calls": 16, "cyc": 0, "pct": 80}
        transition = pp.Window("Fx", 288, 144, 17, 32, 1)
        transition.marker = marker
        transition.counters["shader"] = {
            "us": 15 * 100_000, "calls": 15, "cyc": 0, "pct": 80}

        rows = pp.clean_hold_rows([held, transition], "shader", None)

        self.assertEqual(rows, [("12", 2, 1, 50.0, 1.0, "")])


if __name__ == "__main__":
    unittest.main()
