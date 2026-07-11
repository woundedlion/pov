"""Rank Quilter autoplace/autoroute candidates by signal integrity + ergonomics.

Parses each candidate's `phantasm_unplaced.kicad_pcb` (KiCad 10 s-expr) directly
-- no KiCad needed -- and scores the fast nets (the 24 MHz SPI DATA/CLK to the
strip, DATA_IN/CLK_IN from the Teensy, and the SYNC pair) plus placement quality.

Usage:
    python analyze_candidates.py [DIR ...]

With no args it globs `../Quilter_phantasm_unplaced.kicad_pcb_Candidate_*` next to
the project. Pass explicit candidate folders (or .kicad_pcb files) to override.

A DRC gate runs kicad-cli on each candidate (override path with env KICAD_CLI; the
gate is skipped if it is missing) so a geometry-clean but DRC-broken board can't win
the ranking. Errors are split: 'refill-fixable' clearance/hole errors (Quilter exports
pours without via antipads -- they clear on a KiCad zone refill) vs 'REAL FAULTS'
(shorts/crossings/opens), which disqualify a candidate from the recommended pick.

What matters here, and why:
- The board is 4-layer SIG/GND/GND/SIG with BOTH inner layers poured GND, so any
  outer-layer trace has a solid adjacent reference. A signal via only hops between
  two *same-net* GND planes (benign), but it is still a stub + impedance bump --
  so fewer fast-net vias and shorter fast nets win. DATA/CLK staying on ONE layer
  (zero vias) is ideal: one continuous reference, no return-path transition.
- Connectors J1-J4 are locked, so connector accessibility is equal; ergonomics is
  decided by how the loose parts group (decoupling near U1, terminators near the
  strip connector, the high-Z sync divider kept tight).
"""
import glob
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)

# kicad-cli for the DRC gate (override with env KICAD_CLI). DRC is run on each
# candidate so geometry-clean-but-DRC-broken boards (e.g. Quilter zone-fill
# artifacts) can't quietly win the ranking. Set to "" / missing to skip the gate.
KCLI = os.environ.get("KICAD_CLI", r"C:\Program Files\KiCad\10.0\bin\kicad-cli.exe")
# zone clearance/hole errors usually clear on a KiCad zone refill (Quilter exports
# pours without antipads around signal vias) -- flagged separately from real faults.
REFILL_FIXABLE = {"clearance", "hole_clearance"}

# Fast / critical nets for the 24 MHz SPI + sync (post-relabel source-side names
# DATA_SRC/CLK_SRC/SYNC_SRC are pre-terminator stubs -- included as fast too).
SPI = ["DATA", "CLK", "DATA_IN", "CLK_IN", "DATA_SRC", "CLK_SRC"]
SYNC = ["SYNC_BUS", "FRAME_SYNC", "SYNC_SRC"]
CRIT = SPI + SYNC
# the two pairs we care about phase/length matching on
PAIRS = [("DATA", "CLK"), ("DATA_IN", "CLK_IN")]


def blocks(text, kind):
    """Yield the text of each top-level (kind ...) block, depth-matched."""
    out = []
    i = 0
    key = "(" + kind
    while True:
        j = text.find(key, i)
        if j < 0:
            break
        if text[j + len(key)] not in " \n\t(":
            i = j + len(key)
            continue
        depth = 0
        k = j
        while k < len(text):
            c = text[k]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    out.append(text[j:k + 1])
                    break
            k += 1
        i = k + 1
    return out


def net_of(b):
    m = re.search(r'\(net\s+"([^"]*)"\)', b)
    return m.group(1) if m else None


def field(b, name):
    m = re.search(r'\(%s\s+([^\)]*)\)' % re.escape(name), b)
    return m.group(1).strip() if m else None


def xy(b, tag):
    m = re.search(r'\(%s\s+([-\d.]+)\s+([-\d.]+)\)' % tag, b)
    return (float(m.group(1)), float(m.group(2))) if m else None


def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def real_faults(by_type):
    """Count DRC entries that are NOT refill-fixable zone artifacts -- real
    shorts/crossings/opens that disqualify a candidate."""
    return sum(v for t, v in by_type.items() if t not in REFILL_FIXABLE)


def run_drc(pcb_path):
    """Run kicad-cli DRC; return dict(errors, unconnected, by_type Counter) or None
    if kicad-cli is unavailable. Errors are bucketed by rule so refill-fixable zone
    artifacts are distinguishable from real shorts/crossings."""
    if not KCLI or not os.path.exists(KCLI):
        return None
    # Unique report per call + returncode check: a shared fixed path lets a
    # kicad-cli early-exit leave a stale neighbor's report to be misattributed.
    fd, rpt = tempfile.mkstemp(suffix=".rpt", prefix="cand_drc_")
    os.close(fd)
    try:
        r = subprocess.run([KCLI, "pcb", "drc", "--severity-error", pcb_path, "-o", rpt],
                           capture_output=True, timeout=120, check=False)
        if r.returncode != 0:
            return None
        txt = open(rpt, encoding="utf-8").read()
    except (subprocess.SubprocessError, OSError):
        return None
    finally:
        try:
            os.remove(rpt)
        except OSError:
            pass
    by_type = Counter(re.findall(r"^\[(\w+)\]", txt, re.M))
    m = re.search(r"Found (\d+) unconnected", txt)
    unconnected = int(m.group(1)) if m else 0
    return dict(errors=sum(by_type.values()), unconnected=unconnected, by_type=by_type)


def analyze(path):
    text = open(path, encoding="utf-8").read()
    segs, arcs, vias = blocks(text, "segment"), blocks(text, "arc"), blocks(text, "via")

    netlen, netseg, netlayers = {}, {}, {}
    total_len = 0.0
    for b, is_arc in [(s, False) for s in segs] + [(a, True) for a in arcs]:
        st, en = xy(b, "start"), xy(b, "end")
        if not st or not en:
            continue
        dl = dist(st, en)
        if is_arc:
            mid = xy(b, "mid")
            if mid:
                dl = dist(st, mid) + dist(mid, en)
        n = net_of(b)
        ly = field(b, "layer")
        netlen[n] = netlen.get(n, 0) + dl
        netseg[n] = netseg.get(n, 0) + 1
        netlayers.setdefault(n, set()).add((ly or "?").strip('"'))
        total_len += dl

    netvias, gnd_pts, crit_via_pts = {}, [], []
    for b in vias:
        n = net_of(b)
        p = xy(b, "at")
        netvias[n] = netvias.get(n, 0) + 1
        if n == "GND" and p:
            gnd_pts.append(p)
        if n in CRIT and p:
            crit_via_pts.append(p)

    # fast-net vias with a GND stitching via within 1 mm (clean return-path hop)
    near = sum(1 for p in crit_via_pts
               if min((dist(p, g) for g in gnd_pts), default=99) < 1.0)

    # footprint positions for ergonomics
    pos = {}
    for b in blocks(text, "footprint"):
        m = re.search(r'\(property "Reference" "([^"]*)"', b)
        at = re.search(r'\(at\s+([-\d.]+)\s+([-\d.]+)', b)
        if m and at:
            pos[m.group(1)] = (float(at.group(1)), float(at.group(2)))

    def d(a, b):
        return dist(pos[a], pos[b]) if a in pos and b in pos else float("nan")

    ergo = dict(
        decap_u1=min(d("C_DEC1", "U1"), d("C_DEC2", "U1")),   # decoupling near buffer
        term_j2=(d("R_D1", "J2") + d("R_D2", "J2")) / 2,      # terminators near strip
        rs_sync=d("R_S", "J3A"),                              # sync term near daisy
        divider=(d("R1", "R2") + d("R2", "C_SYNC")) / 2,      # high-Z RC divider tight
    )

    return dict(
        nseg=len(segs) + len(arcs), nvia=len(vias), total_len=total_len,
        netlen=netlen, netseg=netseg, netlayers=netlayers, netvias=netvias,
        spi_vias=sum(netvias.get(n, 0) for n in SPI),
        sync_vias=sum(netvias.get(n, 0) for n in SYNC),
        crit_len=sum(netlen.get(n, 0) for n in CRIT),
        crit_vias=len(crit_via_pts), crit_vias_stitched=near,
        gnd_vias=len(gnd_pts), pos=pos, ergo=ergo,
    )


def score(r):
    """Cheap composite scores (0-10) for a quick ranking. Lower fast-net length,
    fewer fast-net vias, and tighter functional grouping all score higher."""
    si = 10.0
    si -= r["crit_len"] / 90.0          # total fast-net copper (reflection/EMI)
    si -= r["spi_vias"] * 0.6           # SPI-bus vias hurt most (stub + plane hop)
    si -= r["sync_vias"] * 0.4
    # (crit_vias_stitched is reported but not penalised: both inner planes are the
    #  same GND net, so a signal-via plane hop is benign vs a split-plane crossing.)
    e = r["ergo"]
    erg = 10.0
    erg -= max(0, e["decap_u1"] - 5) * 0.3
    erg -= max(0, e["divider"] - 4) * 0.25
    erg -= max(0, e["term_j2"] - 12) * 0.08
    return max(0.0, min(10.0, si)), max(0.0, min(10.0, erg))


def main(argv):
    args = argv[1:]
    if args:
        dirs = args
    else:
        dirs = sorted(glob.glob(os.path.join(
            PROJ, "Quilter_phantasm_unplaced.kicad_pcb_Candidate_*")))
    if not dirs:
        print("no candidate folders found (pass them as args)")
        return 1

    R = {}
    for d in dirs:
        f = d if d.endswith(".kicad_pcb") else os.path.join(d, "phantasm_unplaced.kicad_pcb")
        if not os.path.exists(f):
            print(f"skip {d}: no phantasm_unplaced.kicad_pcb")
            continue
        # name from the Candidate_N anywhere in the path (works whether a folder or a
        # .kicad_pcb file was passed, and regardless of the board's filename)
        m = re.search(r"Candidate_(\w+)", d)
        name = m.group(1) if m else os.path.splitext(os.path.basename(d))[0]
        R[name] = analyze(f)
        R[name]["drc"] = run_drc(f)

    drc_ran = any(R[k]["drc"] is not None for k in R)
    print("=" * 76)
    print("OVERALL ROUTING" + ("  +  DRC GATE" if drc_ran else "  (DRC skipped: no kicad-cli)"))
    print("=" * 76)
    hdr = f"{'Cand':>6} {'tracks':>7} {'vias':>5} {'totLen(mm)':>11} {'GNDvias':>8}"
    if drc_ran:
        hdr += f" {'DRCerr':>7} {'unconn':>6}  flag"
    print(hdr)
    for k in sorted(R):
        r = R[k]
        line = f"{k:>6} {r['nseg']:>7} {r['nvia']:>5} {r['total_len']:>11.1f} {r['gnd_vias']:>8}"
        if drc_ran:
            dc = r["drc"] or {}
            err, unc = dc.get("errors", 0), dc.get("unconnected", 0)
            real = real_faults(dc.get("by_type", {}))
            if err == 0 and unc == 0:
                flag = "clean"
            elif real == 0 and unc == 0:
                flag = "refill-fixable (zone)"   # only clearance/hole vs zones
            else:
                flag = "REAL FAULTS"             # shorts/crossings/unconnected
            line += f" {err:>7} {unc:>6}  {flag}"
        print(line)
    if drc_ran:
        print("  note: 'refill-fixable' = clearance/hole errors that clear on a KiCad zone"
              " refill (Quilter omits via antipads); 'REAL FAULTS' = shorts/crossings/opens.")

    print("\n" + "=" * 76)
    print("CRITICAL NETS  -  length(mm) / vias / layers")
    print("=" * 76)
    for net in CRIT:
        if not any(net in R[k]["netlen"] or net in R[k]["netvias"] for k in R):
            continue
        print(f"\n--- {net} ---")
        print(f"{'Cand':>6} {'len':>8} {'vias':>5} {'segs':>5}  layers")
        for k in sorted(R):
            r = R[k]
            ly = ",".join(sorted(r["netlayers"].get(net, set())))
            print(f"{k:>6} {r['netlen'].get(net,0):>8.1f} {r['netvias'].get(net,0):>5} "
                  f"{r['netseg'].get(net,0):>5}  {ly}")

    print("\n" + "=" * 76)
    print("FAST-NET SUMMARY  (the SI levers)")
    print("=" * 76)
    print(f"{'Cand':>6} {'SPIvias':>8} {'SYNCvias':>9} {'critLen':>9} "
          f"{'vias/stitched':>14}")
    for k in sorted(R):
        r = R[k]
        print(f"{k:>6} {r['spi_vias']:>8} {r['sync_vias']:>9} {r['crit_len']:>9.1f} "
              f"{r['crit_vias']:>7}/{r['crit_vias_stitched']:<6}")

    print("\n" + "=" * 76)
    print("PAIR SKEW (mm)")
    print("=" * 76)
    hdr = "  ".join(f"{a}-{b}" for a, b in PAIRS)
    print(f"{'Cand':>6}  {hdr}")
    for k in sorted(R):
        r = R[k]
        cells = []
        for a, b in PAIRS:
            cells.append(f"{abs(r['netlen'].get(a,0)-r['netlen'].get(b,0)):>7.1f}")
        print(f"{k:>6}  " + "   ".join(cells))

    print("\n" + "=" * 76)
    print("ERGONOMICS  (mm; smaller = better-grouped)")
    print("=" * 76)
    print(f"{'Cand':>6} {'decap->U1':>10} {'term->J2':>9} {'R_S->J3A':>9} {'divider':>8}")
    for k in sorted(R):
        e = R[k]["ergo"]
        print(f"{k:>6} {e['decap_u1']:>10.1f} {e['term_j2']:>9.1f} "
              f"{e['rs_sync']:>9.1f} {e['divider']:>8.1f}")

    print("\n" + "=" * 76)
    print("COMPOSITE SCORE (rough, 0-10)")
    print("=" * 76)
    print(f"{'Cand':>6} {'SignalInteg':>12} {'Ergonomics':>11} {'Overall':>8}  DRC")
    ranked = []
    for k in sorted(R):
        si, erg = score(R[k])
        overall = 0.6 * si + 0.4 * erg
        ranked.append((overall, k, si, erg))

    def drc_tag(k):
        dc = R[k].get("drc")
        if dc is None:
            return "?"
        real = real_faults(dc["by_type"])
        if dc["errors"] == 0 and dc["unconnected"] == 0:
            return "ok"
        return "REAL" if (real or dc["unconnected"]) else "refill"

    for overall, k, si, erg in sorted(ranked, reverse=True):
        print(f"{k:>6} {si:>12.1f} {erg:>11.1f} {overall:>8.1f}  {drc_tag(k)}")
    if ranked:
        # a candidate with REAL faults (not refill-fixable) is disqualified from the pick
        eligible = [t for t in ranked if drc_tag(t[1]) != "REAL"]
        best = max(eligible or ranked)[1]
        note = "" if (eligible and max(eligible) == max(ranked)) else \
            "  (top geometric scorer had REAL DRC faults -- skipped)"
        print(f"\n>> best by composite: Candidate {best}{note}")
        print("   verify by eye -- these are auto-routed; refill zones + DRC, then"
              " hand-polish the winner. 'refill'-tagged DRC clears on a KiCad zone refill.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
