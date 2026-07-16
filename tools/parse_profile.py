"""Parse an on-device HS_PROFILE capture log into per-window and per-preset views.

Companion to targets/Profile/Profile.ino + tools/profile_capture.py (see the
teensy-profile skill). Reads a capture produced by `profile`/`profile_o3` and:

  windows   per-window per-frame cost for one counter scope
  presets   per-preset/-shape/-mode table, each row read from that preset's
            clean-hold windows (modal draw-call count — transition windows,
            whose call count differs, are excluded)
  validate  sanity checks a cycling-effect capture: preset markers present,
            the cycle wraps back to its first index, the effect instance is
            never torn down mid-capture (frame numbers stay monotonic), and the
            root counter matches the wall clock

Counter lines are WINDOW TOTALS since the previous dump; per-frame = total /
window frames, us -> ms / 1000. A display window is 62.5 ms at 480 RPM, so a
draw_frame's wall time quantizes to whole 62.5 ms windows.

Preset markers the effects emit (one per advance), matched here:
  Preset: <i>/<N>            Presets<>::next()/prev() users + DreamBalls
  Shape: <i>/<N>             ShapeShifter, MeshFeedback solids
  Mode: <i>/<N>              SphericalHarmonics
  Spawning Shape: <name> ... IslamicStars (carries V/E/F/I)
  Loading shape: '<name>'    HankinSolids
"""

import argparse
import re
import sys
from collections import Counter

CPU_HZ = 600_000_000  # 600 MHz; cyc / 600 = us

HEADER_RE = re.compile(
    r"=== profile (\S+) \[(\d+)x(\d+)\] frames (\d+)-(\d+) window=(\d+) us ===")
WALL_RE = re.compile(
    r"frame wall us: min=(\d+) avg=(\d+) max=(\d+) sum=(\d+) \((\d+) frames\)")
COUNTER_RE = re.compile(
    r"^(\s*)(\S+)\s+(\d+) us \((\d+)%\)\s+(\d+) calls\s+(\d+) cyc\s*$")

# Preset/shape/mode advance markers. `key` groups them; `idx`/`total`/`name`
# are pulled when present.
MARKER_RES = [
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
        self.counters = {}  # label -> dict(us, pct, calls, cyc, depth)
        self.marker = None  # active preset marker at this window

    def per_frame_ms(self, label):
        n = self.counters.get(label)
        return None if not n else n["us"] / self.frames / 1000.0

    def calls_per_frame(self, label):
        n = self.counters.get(label)
        return None if not n else n["calls"] / self.frames

    def root_us(self):
        n = self.counters.get("frame")
        return n["us"] if n else None


def parse(path):
    """Return (windows, effect_name). Markers are attached to trailing windows."""
    windows = []
    cur = None
    pending_marker = None
    effect = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = HEADER_RE.search(line)
            if m:
                effect = effect or m.group(1)
                nw = Window(m.group(1), int(m.group(2)), int(m.group(3)),
                            int(m.group(4)), int(m.group(5)), int(m.group(6)))
                # Frame numbers restarting = the effect was torn down and
                # reconstructed (epoch); the fresh instance is back on its
                # first preset, so the old marker no longer applies.
                if cur and nw.f_start <= cur.f_start:
                    pending_marker = None
                cur = nw
                cur.marker = pending_marker
                windows.append(cur)
                continue
            m = WALL_RE.search(line)
            if m and cur:
                cur.wall = tuple(int(m.group(i)) for i in range(1, 5))
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
                pending_marker = _marker(key, mm)
                break
    return windows, effect


def _marker(key, mm):
    if key in ("preset", "shape", "mode"):
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


def cmd_windows(windows, scope):
    print(f"# window  frames        {scope} ms/f  calls/f  wall_ms  marker")
    for w in windows:
        pf = w.per_frame_ms(scope)
        cf = w.calls_per_frame(scope)
        wall = w.wall[3] / w.frames / 1000.0 if w.wall else 0.0
        mk = w.marker["name"] if w.marker else "-"
        pf_s = f"{pf:8.2f}" if pf is not None else "     n/a"
        cf_s = f"{cf:7.1f}" if cf is not None else "    n/a"
        print(f"{w.f_start:6d}-{w.f_end:<6d} {pf_s}  {cf_s}  {wall:7.2f}  {mk}")


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
    if idxs:
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

    # Exactness: root cyc/600 vs wall sum for the richest window.
    w = max(windows, key=lambda w: (w.root_us() or 0))
    root = w.counters.get("frame")
    if root and w.wall:
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
    ap.add_argument("mode", choices=["windows", "presets", "validate"])
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
    elif args.mode == "presets":
        cmd_presets(windows, scope, args.gate)
    else:
        return 0 if cmd_validate(windows, effect, scope) else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
