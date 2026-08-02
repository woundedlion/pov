"""Parse an on-device HS_PROFILE capture log into per-window and per-preset views.

Companion to targets/Profile/Profile.ino + tools/profile_capture.py (see the
teensy-profile skill). Reads a capture produced by `profile`/`profile_o3` and:

  windows   per-window per-frame cost for one counter scope; the footer gives
            the pass aggregate the READMEs quote (peak render + spilled)
  presets   per-preset/-shape/-mode table, each row read from that preset's
            clean-hold windows (modal draw-call count — transition windows,
            whose call count differs, are excluded)
  buckets   per-preset cadence buckets (how many presets lock / flap / slip a
            tier) — fills a cycling effect's README cell. Unlike `presets`,
            each preset owns the transition that follows it, so its counts are
            stricter than the clean-hold view.
  metrics   per-window scan-probe counts and their path split, from a capture
            built with -D HS_SCAN_METRICS. Counts only: every probe in such a
            build pays a global increment, so its times are not comparable.
  probe     per-probe stage cycle split, from a capture built with
            -D HS_PROBE_BREAKDOWN. Ratios only: every stage boundary is a
            cycle-counter read, whose measured cost the view subtracts.
  plot      per-window Plot workload attribution, from a capture built with
            -D HS_PLOT_COUNTS. Counts only; timings from this build are perturbed.
  msp-counts MindSplatter particle/gate/raster counts from a dedicated count image.
  msp-stalls MindSplatter short-batch CYCCNT/CPICNT/LSUCNT/EXCCNT attribution.
  validate  sanity checks a cycling-effect capture: preset markers present,
            the cycle wraps back to its first index, the effect instance is
            never torn down mid-capture (frame numbers stay monotonic), and the
            root counter matches the wall clock

Counter lines are WINDOW TOTALS since the previous dump; per-frame = total /
window frames, us -> ms / 1000. A display window is 62.5 ms at 480 RPM, so a
draw_frame's wall time quantizes to whole 62.5 ms windows.

Preset markers the effects emit (one per advance), matched here:
  Preset: <i>/<N>            Presets<>::next()/prev() users + DreamBalls
  Shape: <i>/<N>             MeshFeedback solids
  Mode: <i>/<N>              SphericalHarmonics
  Spawning Shape: <name> ... IslamicStars (carries V/E/F/I)
  Loading shape: '<name>'    HankinSolids
"""

import argparse
import re
import sys
from collections import Counter

CPU_HZ = 600_000_000  # 600 MHz; cyc / 600 = us
DISPLAY_WINDOW_US = 62_500  # half-revolution at 480 RPM

HEADER_RE = re.compile(
    r"=== profile (\S+) \[(\d+)x(\d+)\] frames (\d+)-(\d+) window=(\d+) us ===")
WALL_RE = re.compile(
    r"frame wall us: min=(\d+) avg=(\d+) max=(\d+) sum=(\d+) \((\d+) frames\)")
RENDER_RE = re.compile(r"frame render us: avg=(\d+) max=(\d+)")
FRAME_RE = re.compile(r"^f (\d+) w=(\d+) r=(\d+)$")
COUNTER_RE = re.compile(
    r"^(\s*)(\S+)\s+(\d+) us \((\d+)%\)\s+(\d+) calls\s+(\d+) cyc\s*$")
# HS_SCAN_METRICS window totals (Profile.ino dump_scan_totals).
SCAN_RE = re.compile(
    r"^scan totals: tested=(\d+) culled=(\d+) lut=(\d+) convex=(\d+) "
    r"sector=(\d+) walk=(\d+) cand=(\d+) backstop=(\d+)\s*$")
SCAN_FIELDS = ("tested", "culled", "lut", "convex", "sector", "walk",
               "cand", "backstop")
# HS_PROBE_BREAKDOWN window totals (Profile.ino dump_probe_breakdown).
PROBE_CYC_RE = re.compile(
    r"^probe cycles: point=(\d+) project=(\d+) lut=(\d+) convex=(\d+) "
    r"sector=(\d+) exact=(\d+) pack=(\d+) alpha=(\d+) tick=(\d+)\s*$")
PROBE_CYC_FIELDS = ("point", "project", "lut", "convex", "sector", "exact",
                    "pack", "alpha", "tick")
PROBE_CNT_RE = re.compile(
    r"^probe counts: probe=(\d+) cull_cos=(\d+) cull_r=(\d+) lut=(\d+) "
    r"convex=(\d+) sector=(\d+) exact=(\d+) alpha=(\d+)\s*$")
PROBE_CNT_FIELDS = ("n_probe", "n_cull_cos", "n_cull_r", "n_lut", "n_convex",
                    "n_sector", "n_exact", "n_alpha")
PLOT_RE = re.compile(
    r"^plot counts: r=(\d+),e=(\d+),p=(\d+),g=(\d+),d=(\d+),o=(\d+),"
    r"t=(\d+),c=(\d+),s=(\d+),y=(\d+),u=(\d+),a=(\d+),n=(\d+),"
    r"h=(\d+),x=(\d+),k=(\d+),b=(\d+)\s*$")
PLOT_FIELDS = ("rings", "edges", "planar", "geodesic", "degenerate",
               "one_dot", "cull_tests", "culled", "sim_samples",
               "replay_samples", "planar_unprojects", "planar_arc_samples",
               "normalizations", "shader_calls", "plotted_samples",
               "steps_peak", "backstops")
MSP_PARTICLE_RE = re.compile(
    r"^msp counts particles: resident=(\d+) live=(\d+) full=(\d+) "
    r"partial=(\d+) draining=(\d+)\s*$")
MSP_PARTICLE_FIELDS = ("resident", "live", "full", "partial", "draining")
MSP_GATE_RE = re.compile(
    r"^msp counts gate: cart_lat=(\d+) cart_mer=(\d+) cart_fallback=(\d+) "
    r"row=(\d+) col=(\d+) edge=(\d+) visible=(\d+) exact=(\d+)\s*$")
MSP_GATE_FIELDS = ("cart_lat", "cart_mer", "cart_fallback", "row", "col",
                   "edge_reject", "visible", "exact_fallback")
MSP_RENDER_RE = re.compile(
    r"^msp counts render: dot=(\d+) long=(\d+) adaptive=(\d+) shader=(\d+) "
    r"pal_end=(\d+) pal_lerp=(\d+) hole_early=(\d+)\s*$")
MSP_RENDER_FIELDS = ("one_dot", "long", "adaptive", "shader", "pal_end",
                     "pal_lerp", "hole_early")
MSP_AA_RE = re.compile(
    r"^msp counts aa: tap0=(\d+) tap1=(\d+) tap2=(\d+) tap3=(\d+) "
    r"tap4=(\d+) interior=(\d+) boundary=(\d+)\s*$")
MSP_AA_FIELDS = ("tap0", "tap1", "tap2", "tap3", "tap4", "interior",
                 "boundary")
MSP_STALL_RE = re.compile(
    r"^msp stall: stage=(\S+) batches=(\d+) cyc=(\d+) cpi=(\d+) lsu=(\d+) "
    r"exc=(\d+)\s*$")

# Preset/shape/mode advance markers. `key` groups them; `idx`/`total`/`name`
# are pulled when present.
MARKER_RES = [
    ("profile_preset", re.compile(r"^Profile preset: (\d+)/(\d+)")),
    ("preset", re.compile(r"^Preset: (\d+)/(\d+)")),
    ("shape", re.compile(r"^Shape: (\d+)/(\d+)")),
    ("mode", re.compile(r"^Mode: (\d+)/(\d+)")),
    ("islamic", re.compile(
        r"^Spawning Shape: (\S+) \(V=(\d+), E=(\d+), F=(\d+), I=(\d+)\)")),
    ("hankin", re.compile(r"^Loading shape: '([^']*)'")),
]


class Window:
    """One readout window: header stats plus its counter tree (label -> node)."""

    def __init__(self, effect, w, h, f_start, f_end, window_us):
        self.effect = effect
        self.w, self.h = w, h
        self.f_start, self.f_end = f_start, f_end
        self.window_us = window_us
        self.frames = f_end - f_start + 1
        self.wall = None  # (min, avg, max, sum)
        self.render = None  # (avg, max) per-frame render us (wall - sync wait)
        # per-frame (n, wall_us, render_us, owner_marker), when logged. The
        # owner is stamped per FRAME, not per window: a window straddling a
        # preset advance holds frames of both, so only the frame carries a
        # truthful attribution.
        self.frame_rows = []
        self.counters = {}  # label -> dict(us, pct, calls, cyc, depth)
        self.marker = None  # active preset marker at this window
        self.scan = None  # HS_SCAN_METRICS window totals, when the build has them
        self.probe = None  # HS_PROBE_BREAKDOWN window buckets + counts
        self.plot = None  # HS_PLOT_COUNTS window workload
        self.msp_counts = None  # dedicated MindSplatter count image
        self.msp_stalls = {}  # stage -> short-batch DWT totals

    def per_frame_ms(self, label):
        n = self.counters.get(label)
        return None if not n else n["us"] / self.frames / 1000.0

    def calls_per_frame(self, label):
        n = self.counters.get(label)
        return None if not n else n["calls"] / self.frames

    def root_us(self):
        n = self.counters.get("frame")
        return n["us"] if n else None

    def render_ms(self):
        """Per-frame render (frame minus the display-sync wait), averaged.

        None when the effect opens no *_buffer_wait scope: the sync idle is
        then indistinguishable from work, and root/frames would be wall. Read
        such an effect through its own render scope (--scope <xx>_render).
        """
        root = self.counters.get("frame")
        if not root:
            return None
        wait = next((n["us"] for label, n in self.counters.items()
                     if label.endswith("buffer_wait")), None)
        if wait is None:
            return None
        return (root["us"] - wait) / self.frames / 1000.0

    def render_is_wall(self):
        """The telemetry's render is really wall: the effect has no wait scope.

        Profile.ino derives each frame's render as wall minus the effect's
        *_buffer_wait counter delta. An effect that never opens such a scope
        yields a zero delta, so render == wall on every frame and the "peak
        render" would silently be a peak WALL -- render plus the sync idle,
        quantized up to a whole display window.
        """
        return bool(self.frame_rows) and all(f[1] == f[2] for f in self.frame_rows)

    def peak_render_ms(self):
        """(ms, exact) worst frame's render.

        Exact from the per-frame telemetry when the capture has it; otherwise
        a placeholder: this window's MEAN render, which understates the peak
        it stands in for. Wall time is render + sync wait, so its max is not a
        peak render and is never used here -- including when the telemetry
        itself is wall for want of a *_buffer_wait scope, which falls back to
        the same placeholder rather than passing wall off as a render peak.
        """
        if self.render and not self.render_is_wall():
            return (self.render[1] / 1000.0, True)
        r = self.render_ms()
        return (None, False) if r is None else (r, False)


def parse(path):
    """Return (windows, effect_name). Markers are attached to trailing windows.

    A marker logged mid-window (the effect advanced during a frame, so frames of
    the OUTGOING preset already streamed into the open window) takes effect at
    the next window boundary. Attaching it to the window about to close would
    credit the outgoing preset's frames -- up to 15 of the 16 -- to the incoming
    one, which reads as a spurious peak on whichever preset follows an expensive
    one. A marker logged between windows applies immediately.
    """
    windows = []
    cur = None
    active_marker = None    # in force for the window being dumped
    deferred_marker = None  # seen mid-window; applies from the next boundary
    frame_owner = None      # preset owning the frames streaming right now
    effect = None
    pending_frames = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = FRAME_RE.match(line)
            if m:
                # A marker logged since the last frame line names the preset the
                # advance switched TO; this frame is its first (on an effect that
                # rebuilds on advance, the one that pays the rebuild).
                if deferred_marker is not None:
                    frame_owner = deferred_marker
                pending_frames.append(
                    tuple(int(m.group(i)) for i in (1, 2, 3)) + (frame_owner,))
                continue
            m = HEADER_RE.search(line)
            if m:
                effect = effect or m.group(1)
                nw = Window(m.group(1), int(m.group(2)), int(m.group(3)),
                            int(m.group(4)), int(m.group(5)), int(m.group(6)))
                # Frame numbers restarting = the effect was torn down and
                # reconstructed (epoch); the fresh instance is back on its
                # first preset, so the old marker no longer applies.
                if cur and nw.f_start <= cur.f_start:
                    active_marker = deferred_marker = frame_owner = None
                cur = nw
                cur.marker = active_marker
                if deferred_marker is not None:
                    active_marker, deferred_marker = deferred_marker, None
                # The per-frame lines stream before their window's dump.
                cur.frame_rows = [f for f in pending_frames
                                  if nw.f_start <= f[0] <= nw.f_end]
                pending_frames = []
                windows.append(cur)
                continue
            m = WALL_RE.search(line)
            if m and cur:
                cur.wall = tuple(int(m.group(i)) for i in range(1, 5))
                continue
            m = RENDER_RE.search(line)
            if m and cur:
                cur.render = (int(m.group(1)), int(m.group(2)))
                continue
            m = SCAN_RE.match(line)
            if m and cur:
                cur.scan = dict(zip(SCAN_FIELDS,
                                    (int(g) for g in m.groups())))
                continue
            m = PROBE_CYC_RE.match(line)
            if m and cur:
                cur.probe = dict(cur.probe or {})
                cur.probe.update(zip(PROBE_CYC_FIELDS,
                                     (int(g) for g in m.groups())))
                continue
            m = PROBE_CNT_RE.match(line)
            if m and cur:
                cur.probe = dict(cur.probe or {})
                cur.probe.update(zip(PROBE_CNT_FIELDS,
                                     (int(g) for g in m.groups())))
                continue
            m = PLOT_RE.match(line)
            if m and cur:
                cur.plot = dict(zip(PLOT_FIELDS,
                                    (int(g) for g in m.groups())))
                continue
            matched_msp_counts = False
            for regex, fields in ((MSP_PARTICLE_RE, MSP_PARTICLE_FIELDS),
                                  (MSP_GATE_RE, MSP_GATE_FIELDS),
                                  (MSP_RENDER_RE, MSP_RENDER_FIELDS),
                                  (MSP_AA_RE, MSP_AA_FIELDS)):
                m = regex.match(line)
                if m and cur:
                    cur.msp_counts = dict(cur.msp_counts or {})
                    cur.msp_counts.update(zip(fields,
                                              (int(g) for g in m.groups())))
                    matched_msp_counts = True
                    break
            if matched_msp_counts:
                continue
            m = MSP_STALL_RE.match(line)
            if m and cur:
                cur.msp_stalls[m.group(1)] = dict(
                    zip(("batches", "cyc", "cpi", "lsu", "exc"),
                        (int(g) for g in m.groups()[1:])))
                continue
            m = COUNTER_RE.match(line)
            if m and cur:
                indent, label = m.group(1), m.group(2)
                cur.counters[label] = dict(
                    us=int(m.group(3)), pct=int(m.group(4)),
                    calls=int(m.group(5)), cyc=int(m.group(6)),
                    depth=len(indent) // 2)
                continue
            for key, rgx in MARKER_RES:
                mm = rgx.match(line.strip())
                if not mm:
                    continue
                mk = _marker(key, mm)
                deferred_marker = mk
                if not pending_frames:
                    # Between windows: no outgoing frames to protect.
                    active_marker = mk
                break
    return windows, effect


def _marker(key, mm):
    if key in ("profile_preset", "preset", "shape", "mode"):
        return dict(key=key, idx=int(mm.group(1)), total=int(mm.group(2)),
                    name=str(int(mm.group(1))))
    if key == "islamic":
        return dict(key=key, name=mm.group(1), V=int(mm.group(2)),
                    E=int(mm.group(3)), F=int(mm.group(4)), I=int(mm.group(5)))
    return dict(key=key, name=mm.group(1))  # hankin


def dominant_leaf(windows):
    """Label of the costliest non-root, non-buffer_wait scope across the run."""
    tot = Counter()
    for w in windows:
        for label, n in w.counters.items():
            if label == "frame" or label.endswith("buffer_wait"):
                continue
            tot[label] += n["us"]
    return tot.most_common(1)[0][0] if tot else "frame"


def spilled_frames(w):
    """Frames of this window that took >1 display window (62.5 ms).

    Counts FRAMES, not missed flips: a 180 ms frame misses two flips but is one
    spilled frame, so spilled/frames stays a true fraction the cadence colour
    thresholds can be read against.
    """
    if w.frame_rows and not w.render_is_wall():
        # Renders start flip-aligned (the buffer_free gate opens at a flip),
        # so a frame spills exactly when its render exceeds one window.
        return sum(1 for f in w.frame_rows if f[2] > DISPLAY_WINDOW_US)
    if not w.wall:
        return 0
    # Either no per-frame rows, or they carry wall for want of a *_buffer_wait
    # scope. Wall already includes the sync idle and quantizes to whole
    # windows, so testing it against one window would count jitter above 62.5
    # as a spill; the sum-derived estimate below is the treatment wall data
    # gets either way.
    # Without per-frame rows only the window's wall sum is known, and it gives
    # the extra windows consumed (= missed flips), which bounds the spilled
    # frames from above: at most every frame spilled. Marked '~' at the callers.
    extra = max(0, round(w.wall[3] / DISPLAY_WINDOW_US) - w.frames)
    return min(extra, w.frames)


def cmd_windows(windows, scope):
    print(f"# window  frames        {scope} ms/f  calls/f  wall_ms  "
          f"peak_ms  spill  marker")
    # peak_ms is the worst frame's RENDER: exact from the per-frame telemetry,
    # else '~' = the worst window's MEAN render as a placeholder until the
    # effect is re-captured. Wall time is render + the display-sync wait, so
    # its max is not a peak render and is never substituted here.
    agg_spill, agg_frames = 0, 0
    worst_peak = 0.0
    worst_peak_w = None
    exact_run = all(w.render and not w.render_is_wall() for w in windows)
    for w in windows:
        pf = w.per_frame_ms(scope)
        cf = w.calls_per_frame(scope)
        wall = w.wall[3] / w.frames / 1000.0 if w.wall else 0.0
        peak, exact = w.peak_render_ms()
        spill = spilled_frames(w)
        mk = w.marker["name"] if w.marker else "-"
        pf_s = f"{pf:8.2f}" if pf is not None else "     n/a"
        cf_s = f"{cf:7.1f}" if cf is not None else "    n/a"
        peak_s = ("    n/a " if peak is None
                  else f"{peak:7.2f}" + (" " if exact else "~"))
        print(f"{w.f_start:6d}-{w.f_end:<6d} {pf_s}  {cf_s}  {wall:7.2f}  "
              f"{peak_s}{spill:5d}  {mk}")
        if peak is not None and peak > worst_peak:
            worst_peak, worst_peak_w = peak, w
        agg_spill += spill
        agg_frames += w.frames
    if worst_peak_w is None:
        print(f"# aggregate: peak render n/a, "
              f"spilled {agg_spill}/{agg_frames} frames")
        return
    kind = ("peak frame render" if exact_run else
            "peak WINDOW MEAN render (~placeholder, no per-frame telemetry)")
    print(f"# aggregate: {kind} {worst_peak:.2f} ms "
          f"(frames {worst_peak_w.f_start}-{worst_peak_w.f_end}), "
          f"spilled {agg_spill}/{agg_frames} frames"
          f"{'' if exact_run else ' (wall-derived, +-1)'}")


def cmd_presets(windows, scope, gate):
    """One row per preset, read from its modal-call-count (clean-hold) windows."""
    gate = gate or scope
    # Windows before the first marker belong to the initial preset: index-style
    # cyclers construct on entry 0 and only log on advance, so fold the
    # unmarked prefix into index 0's group.
    idx_key = next((w.marker["key"] for w in windows
                    if w.marker and "idx" in w.marker), None)
    groups = {}
    for w in windows:
        if w.marker is None:
            key = (idx_key, "0") if idx_key else ("start", "0")
        else:
            key = (w.marker["key"], w.marker.get("name"))
        groups.setdefault(key, []).append(w)

    print(f"# preset   holds  clean  {scope} ms/f  calls/f  meta")
    rows = []
    for key, ws in groups.items():
        gate_counts = [round(w.calls_per_frame(gate)) for w in ws
                       if w.calls_per_frame(gate) is not None]
        clean = ws
        if gate_counts:
            modal = Counter(gate_counts).most_common(1)[0][0]
            clean = [w for w in ws
                     if w.calls_per_frame(gate) is not None
                     and round(w.calls_per_frame(gate)) == modal]
        vals = [(w.per_frame_ms(scope), w) for w in clean
                if w.per_frame_ms(scope) is not None]
        if not vals:
            continue
        best_ms, best_w = max(vals, key=lambda t: t[0])
        cf = best_w.calls_per_frame(scope)
        meta = ""
        mk = best_w.marker
        if mk and mk["key"] == "islamic":
            meta = f"V={mk['V']} E={mk['E']} F={mk['F']} I={mk['I']}"
        elif mk and mk["key"] == "hankin":
            meta = mk["name"]
        rows.append((key[1], len(ws), len(clean), best_ms, cf, meta))
    for name, holds, clean, ms, cf, meta in sorted(
            rows, key=lambda r: -r[3]):
        print(f"{name:>8}  {holds:5d}  {clean:5d}  {ms:8.2f}  {cf:7.1f}  {meta}")


def cmd_buckets(windows):
    """Per-preset cadence buckets: how many presets hold 16 fps vs spill.

    A preset owns the FRAMES its marker was in force for, not whole windows: the
    window straddling an advance holds frames of both the outgoing and incoming
    preset, so crediting it to either invents a peak the preset never rendered
    (whichever preset neighbours an expensive one inherits its cost). Per-frame
    attribution needs the per-frame telemetry; without it this falls back to the
    window-level split, which carries that error.

    Colour is binary: no spill anywhere is green; any spilled frame is red.
    """
    idx_key = next((w.marker["key"] for w in windows
                    if w.marker and "idx" in w.marker), None)

    def group_key(mk):
        if mk is None:
            return (idx_key, "0") if idx_key else ("start", "0")
        return (mk["key"], mk.get("name"))

    exact_run = all(w.render and not w.render_is_wall() for w in windows)
    have_frames = exact_run and all(w.frame_rows for w in windows)

    groups = {}
    if have_frames:
        for w in windows:
            for f in w.frame_rows:
                groups.setdefault(group_key(f[3]), []).append(f)
    else:
        for w in windows:
            groups.setdefault(group_key(w.marker), []).append(w)

    rows = []
    for key, items in groups.items():
        if have_frames:
            peak = max(f[2] for f in items) / 1000.0
            spill = sum(1 for f in items if f[2] > DISPLAY_WINDOW_US)
            frames = len(items)
        else:
            peaks = [p for p, _ in (w.peak_render_ms() for w in items)
                     if p is not None]
            if not peaks:
                continue
            peak = max(peaks)
            spill = sum(spilled_frames(w) for w in items)
            frames = sum(w.frames for w in items)
        colour = "green" if spill == 0 else "red"
        rows.append((key[1], colour, peak, spill, frames))

    mark = "" if exact_run else "~"
    print(f"# preset  cadence  peak_render_ms{mark}  spilled/frames")
    for name, colour, peak, spill, frames in sorted(
            rows, key=lambda r: (-r[3], -r[2])):
        print(f"{name:>8}  {colour:>7}  {peak:11.2f}{mark}  {spill:5d}/{frames}")
    print("# buckets: " + " ".join(
        f"{c}={sum(1 for r in rows if r[1] == c)}"
        for c in ("green", "red")))
    for c in ("green", "red"):
        sel = [r for r in rows if r[1] == c]
        if not sel:
            continue
        print(f"#   {c}: {len(sel)} presets, peak render "
              f"{min(r[2] for r in sel):.1f}-{max(r[2] for r in sel):.1f} "
              f"ms{mark}, spilled {sum(r[3] for r in sel)}")
    if not exact_run:
        print("#   ~ = window-mean placeholder (understates the true peak); "
              "spilled is wall-derived, +-1. Re-capture for exact peaks.")
    return 0


def cmd_metrics(windows):
    """Per-window scan-probe counts and the Face::distance path split.

    Probes are Face::distance calls; culled are those the back-face / radius
    guards reject before any edge work. lut/convex/sector/walk partition the
    survivors. `cand` counts pixels passing the scan's d < pixel_width test and
    `shade` the ones that then survived the alpha test, read from the
    raster_shade scope's call count (the two differ by the alpha-rejected AA
    fringe). probes/shade is the figure a cycles-per-probe estimate divides by.
    """
    have = [w for w in windows if w.scan]
    if not have:
        print("no 'scan totals' lines: rebuild with -D HS_SCAN_METRICS",
              file=sys.stderr)
        return 2
    print(f"{'# window':<13} " + " ".join(
        f"{h:>{w}}" for h, w in zip(METRICS_COLS, METRICS_WIDTHS))
        + "  marker")
    agg = Counter()
    agg_frames = 0
    agg_shade = 0
    for w in have:
        s = w.scan
        shade = w.counters.get("raster_shade", {}).get("calls", 0)
        for k in SCAN_FIELDS:
            agg[k] += s[k]
        agg_frames += w.frames
        agg_shade += shade
        print(f"{w.f_start:6d}-{w.f_end:<6d} " + _metrics_row(s, w.frames, shade)
              + f"  {w.marker['name'] if w.marker else '-'}")
    print(f"{'# aggregate':<13} " + _metrics_row(agg, agg_frames, agg_shade)
          + f"  {agg_frames} frames")
    return 0


METRICS_COLS = ("probes/f", "culled/f", "lut%", "convex%", "sector%", "walk%",
                "cand/f", "shade/f", "probes/shade")
METRICS_WIDTHS = (9, 9, 7, 8, 8, 7, 7, 8, 13)


def _metrics_row(s, frames, shade):
    """One formatted metrics row: per-frame counts plus the path split."""
    served = s["lut"] + s["convex"] + s["sector"] + s["walk"]

    def pct(v):
        return 100.0 * v / served if served else 0.0

    vals = (s["tested"] / frames, s["culled"] / frames, pct(s["lut"]),
            pct(s["convex"]), pct(s["sector"]), pct(s["walk"]),
            s["cand"] / frames, shade / frames,
            s["tested"] / shade if shade else float("nan"))
    prec = (0, 0, 1, 1, 1, 1, 0, 0, 2)
    return " ".join(f"{v:{w}.{p}f}"
                    for v, w, p in zip(vals, METRICS_WIDTHS, prec))


PROBE_STAGES = (("point", "n_probe"), ("project", "n_probe"),
                ("lut", "n_lut"), ("convex", "n_convex"),
                ("sector", "n_sector"), ("exact", "n_exact"),
                ("pack", "n_probe"), ("alpha", "n_alpha"))


def cmd_probe(windows):
    """Per-probe stage cycle split from an HS_PROBE_BREAKDOWN capture.

    Each stage prints its mean cycles per event of its own denominator (`point`,
    `project` and `pack` run once per probe; the edge stages once per probe on
    that path; `alpha` once per shade candidate), the same mean with the
    self-measured counter-read cost removed, and the share of the summed
    per-probe cost it accounts for. `read` is that measured cost: `tick` sums a
    back-to-back read pair per probe, so read = tick/(2*probes), and every stage
    carries exactly one read. Ratios from this build are meaningful; absolute
    times are not.
    """
    have = [w for w in windows if w.probe and "n_probe" in w.probe
            and "point" in w.probe]
    if not have:
        print("no 'probe cycles'/'probe counts' lines: rebuild with "
              "-D HS_PROBE_BREAKDOWN", file=sys.stderr)
        return 2
    agg = Counter()
    for w in have:
        agg.update(w.probe)
    probes = agg["n_probe"]
    if not probes:
        print("probe count is zero", file=sys.stderr)
        return 2
    read = agg["tick"] / (2.0 * probes)
    print(f"probes={probes}  windows={len(have)}  "
          f"counter read={read:.1f} cyc (measured)")
    print(f"{'# stage':<10} {'events':>10} {'cyc/event':>10} "
          f"{'net':>9} {'cyc/probe':>10} {'share':>7}")
    net_per_probe = {}
    for stage, den_key in PROBE_STAGES:
        den = agg[den_key]
        if not den:
            continue
        mean = agg[stage] / den
        net = mean - read
        net_per_probe[stage] = net * den / probes
    total = sum(v for v in net_per_probe.values() if v > 0)
    for stage, den_key in PROBE_STAGES:
        if stage not in net_per_probe:
            continue
        den = agg[den_key]
        mean = agg[stage] / den
        per_probe = net_per_probe[stage]
        share = 100.0 * per_probe / total if total else 0.0
        print(f"{stage:<10} {den:10d} {mean:10.1f} {mean - read:9.1f} "
              f"{per_probe:10.1f} {share:6.1f}%")
    print(f"{'# total':<10} {'':>10} {'':>10} {'':>9} {total:10.1f} "
          f"{100.0:6.1f}%")
    return 0


PLOT_COLS = ("rings/f", "edges/f", "plan/f", "geo/f", "deg/f", "dot/f",
             "cull/f", "reject/f", "sim/e", "replay/e", "unproj/e",
             "arc/e", "norm/e", "shader/f", "plot/f", "peak", "back")


def _plot_row(p, frames):
    edges = p["edges"]
    per_frame = (p["rings"] / frames, edges / frames, p["planar"] / frames,
                 p["geodesic"] / frames, p["degenerate"] / frames,
                 p["one_dot"] / frames, p["cull_tests"] / frames,
                 p["culled"] / frames)
    per_edge = tuple(p[k] / edges if edges else 0.0 for k in
                     ("sim_samples", "replay_samples", "planar_unprojects",
                      "planar_arc_samples", "normalizations"))
    tail = (p["shader_calls"] / frames, p["plotted_samples"] / frames,
            p["steps_peak"], p["backstops"])
    return " ".join(f"{v:8.1f}" for v in per_frame + per_edge + tail)


def cmd_plot(windows):
    """Per-window count-only Plot rasterizer workload attribution."""
    have = [w for w in windows if w.plot]
    if not have:
        print("no 'plot counts' lines: rebuild with -D HS_PLOT_COUNTS",
              file=sys.stderr)
        return 2
    print(f"{'# window':<13} " + " ".join(f"{h:>8}" for h in PLOT_COLS)
          + "  marker")
    agg = Counter()
    agg_frames = 0
    peak = 0
    for w in have:
        agg.update(w.plot)
        agg_frames += w.frames
        peak = max(peak, w.plot["steps_peak"])
        marker = w.marker["name"] if w.marker else "-"
        print(f"{w.f_start:6d}-{w.f_end:<6d} " +
              _plot_row(w.plot, w.frames) + f"  {marker}")
    agg["steps_peak"] = peak
    print(f"{'# aggregate':<13} " + _plot_row(agg, agg_frames) +
          f"  {agg_frames} frames")
    return 0


def cmd_msp_counts(windows):
    """Aggregate a dedicated MindSplatter count capture."""
    have = [w for w in windows if w.msp_counts]
    if not have:
        print("no 'msp counts' lines: use HS_PROFILE_MINDSPLATTER=counts",
              file=sys.stderr)
        return 2
    agg = Counter()
    frames = 0
    for w in have:
        agg.update(w.msp_counts)
        frames += w.frames
    print(f"frames={frames} windows={len(have)}")
    print(f"{'# counter':<24} {'total':>12} {'per-frame':>12}")
    for field in (MSP_PARTICLE_FIELDS + MSP_GATE_FIELDS + MSP_RENDER_FIELDS +
                  MSP_AA_FIELDS):
        print(f"{field:<24} {agg[field]:12d} {agg[field] / frames:12.2f}")
    return 0


def cmd_msp_stalls(windows):
    """Aggregate MindSplatter short-batch DWT cycle/stall attribution."""
    have = [w for w in windows if w.msp_stalls]
    if not have:
        print("no 'msp stall' lines: use HS_PROFILE_MINDSPLATTER=stalls",
              file=sys.stderr)
        return 2
    stages = {}
    for w in have:
        for stage, values in w.msp_stalls.items():
            stages.setdefault(stage, Counter()).update(values)
    print(f"{'# stage':<22} {'batches':>10} {'cyc/batch':>10} "
          f"{'cpi/batch':>10} {'lsu/batch':>10} {'exc/batch':>10}")
    for stage, values in stages.items():
        batches = values["batches"]
        divisor = batches or 1
        print(f"{stage:<22} {batches:10d} {values['cyc'] / divisor:10.1f} "
              f"{values['cpi'] / divisor:10.2f} "
              f"{values['lsu'] / divisor:10.2f} "
              f"{values['exc'] / divisor:10.2f}")
    return 0


def cmd_validate(windows, effect, scope):
    ok = True

    def check(cond, msg):
        nonlocal ok
        ok = ok and cond
        print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")

    print(f"effect={effect}  windows={len(windows)}")
    check(len(windows) >= 3, f"captured >=3 windows ({len(windows)})")

    # Single instance: frame numbers strictly increase (no mid-capture teardown).
    resets = sum(1 for a, b in zip(windows, windows[1:])
                 if b.f_start <= a.f_start)
    check(resets == 0,
          f"frame numbers monotonic (no epoch reset) ({resets} resets)")

    # Preset markers: only cycling effects emit them. When present, require the
    # cycle to wrap back to its first index; when absent, this is a non-cycling
    # capture and the wrap checks do not apply.
    marks = [w.marker for w in windows if w.marker]
    idxs = [m["idx"] for m in marks if "idx" in m]
    names = [m["name"] for m in marks]
    if not marks:
        print("  [INFO] no preset markers - non-cycling capture")
    selected = [m for m in marks if m["key"] == "profile_preset"]
    if selected:
        choices = {(m["idx"], m.get("total")) for m in selected}
        index, total = selected[0]["idx"], selected[0].get("total")
        check(len(choices) == 1 and total is not None and index < total,
              f"fixed profile preset {index}/{total}")
    elif idxs:
        peak = max(idxs)
        wrapped = any(idxs[i] < idxs[i - 1] for i in range(1, len(idxs)))
        check(wrapped, f"cycle wraps to 0 (peak idx {peak}, "
                       f"{len(set(idxs))} distinct)")
        total = marks[0].get("total")
        if total:
            check(len(set(idxs)) >= total,
                  f"all {total} presets visited ({len(set(idxs))})")
    elif names:
        distinct = len(set(names))
        first_repeat = next((i for i in range(1, len(names))
                             if names[i] in names[:i]), None)
        check(first_repeat is not None,
              f"a shape repeats (cycle closed): {distinct} distinct")

    # Per-frame render telemetry is only render if the effect opened a
    # *_buffer_wait scope for Profile.ino to subtract; without one it is wall,
    # and every peak it feeds is a peak WALL wearing a render's name.
    # The diagnosis belongs on the failing branch only: on the passing one it
    # asserts the opposite of what was just measured, and a reader who takes it
    # at face value carries "this effect has no buffer_wait scope" away from a
    # run that proved it has one.
    wall_render = [w for w in windows if w.render_is_wall()]
    check(not wall_render,
          f"per-frame render is render, not wall ({len(wall_render)} of "
          f"{len(windows)} windows have render == wall"
          + (": the effect opens no *_buffer_wait scope)" if wall_render
             else ")"))

    # Exactness: root cyc/600 vs wall sum for the richest window.
    # Every check above is skip-on-absent, so this one carries the whole
    # certification: without it a capture holding no counters at all is VALID.
    scoped = [w for w in windows if scope in w.counters]
    check(bool(scoped), f"scope '{scope}' appears in at least one window")
    candidates = scoped or windows
    w = max(candidates, key=lambda w: (w.per_frame_ms(scope) or 0))
    root = w.counters.get("frame")
    have = bool(root) and bool(w.wall) and w.wall[3] > 0
    check(have, "richest window carries a root 'frame' counter and a wall sum"
                + ("" if have else " — nothing in this capture is measurable"))
    if have:
        model_us = root["cyc"] / 600.0
        ppm = abs(model_us - w.wall[3]) / w.wall[3] * 1e6
        check(ppm < 100,
              f"root cyc/600 vs wall sum within {ppm:.1f} ppm "
              f"(frames {w.f_start}-{w.f_end})")

    print(f"\n{'VALID' if ok else 'INVALID'}: {effect}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log")
    ap.add_argument("mode", choices=["windows", "presets", "buckets",
                                     "validate", "frames", "metrics",
                                     "probe", "plot", "msp-counts",
                                     "msp-stalls"])
    ap.add_argument("--scope", help="counter label to read (default: costliest leaf)")
    ap.add_argument("--gate", help="call-count scope gating clean holds "
                                    "(default: --scope)")
    args = ap.parse_args()

    windows, effect = parse(args.log)
    if not windows:
        print(f"no windows parsed from {args.log}", file=sys.stderr)
        return 2
    scope = args.scope or dominant_leaf(windows)

    if args.mode == "windows":
        cmd_windows(windows, scope)
    elif args.mode == "frames":
        print("# frame  wall_us  render_us  spill")
        for w in windows:
            for n, wall, render, _owner in w.frame_rows:
                print(f"{n:7d} {wall:8d} {render:9d} {render // DISPLAY_WINDOW_US:5d}")
    elif args.mode == "presets":
        cmd_presets(windows, scope, args.gate)
    elif args.mode == "buckets":
        return cmd_buckets(windows)
    elif args.mode == "metrics":
        return cmd_metrics(windows)
    elif args.mode == "probe":
        return cmd_probe(windows)
    elif args.mode == "plot":
        return cmd_plot(windows)
    elif args.mode == "msp-counts":
        return cmd_msp_counts(windows)
    elif args.mode == "msp-stalls":
        return cmd_msp_stalls(windows)
    else:
        return 0 if cmd_validate(windows, effect, scope) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
