# Raymarch optimization results — 2026-09-20

Changes 1 and 2 landed. They reduce the observed live peak from **62.725 ms to
56.071 ms**, saving **6.654 ms (10.61%)** and leaving **6.429 ms** before the
62.5 ms deadline. Change 3 was implemented and measured in four forms; all
regressed the peak, so its code remains an archived experiment.

## Peak measurements

Every row is a 110-second on-device capture at 600 MHz with the real segmented
driver and interrupts active. All shipping comparisons use COM3 except the
explicitly marked projected/forced-inline experiment. Runtime statistics use
every raw per-frame row after setup frame 1, including the final partial
counter window. No warmup frames, transitions or real spills are discarded.

| Cumulative implementation | Live peak ms | Change from previous accepted version | Spills / live frames | Decision |
|---|--:|--:|--:|---|
| Baseline | 62.725 | — | 2 / 1734 | Reference |
| 1: scalar extent checks + exact probe seed reuse | 59.332 | −3.393 ms (−5.41%) | 0 / 1736 | Landed |
| 1 + 2: twist-2 specialization + compact O3/RGB shader | 56.071 | −3.261 ms (−5.50%) | 0 / 1736 | Landed |
| 1 + 2 + 3: cylinder/slab/hole rejection | 57.255 | +1.184 ms | 0 / 1736 | Rejected |
| Bounds with forced distance inlining | 57.798 | +1.727 ms | 0 / 1736 | Rejected |
| Projected bounds, forced inlining, COM4 | 57.650 | +1.579 ms | 0 / 1736 | Rejected |
| Projected bounds, automatic inlining | 60.001 | +3.930 ms | 0 / 1736 | Rejected |

Every listed peak occurred at frame 349. A second baseline spill occurred at
frame 1345. These are observed peaks for the default seeded trajectory, not a
worst-case proof across all controls and animation durations. Means played no
part in accepting changes. The fresh 70-second baseline was 62.722 ms, within
3 microseconds of the longer baseline's peak.

Raw captures, build logs, source SHAs, ELF hashes and compiler settings are in
`build/prof/raymarch_implementation_20260920/`. Attested profile and full-roster
ELFs are in `build/prof/artifacts/`. The final shipping and global-O3 reports
are linked from the ranked profile READMEs.

## What shipped

1. **Extent checks and seed reuse.** Trace and probe compare scalar travel
   against a precomputed exit distance instead of recalculating a 3D dot
   product at every step. The background probe reuses its first distance from
   the foreground trace. Incremental sample coordinates and the chosen closest
   point retain the original arithmetic.
2. **Torus/shader specialization.** Twist 2 uses the explicit one-step harmonic
   expression, with the original generic recurrence for other controls. The
   compact surface/shader functions receive selective O3. Palette sampling
   returns RGB directly, omitting alpha work whose result was overwritten.
   World-space lighting and its per-hit quaternion rotation remain intact.

No step limit, AA width, resolution, noise-table resolution, shading model or
live control was reduced. No approximate reciprocal or square root was added.

## Quality qualification

The framebuffer harness renders 1,480 full 288×144 frames and compares 78
checkpoints: 400 frames of default evolution plus 27 control scenarios covering
all 21 placement solids, twists 0–8, fill and AA extremes, and step budgets.
That is **3,234,816 compared pixels** per candidate.

- Change 1: bit-identical to baseline. Its differential trace/probe oracle
  covers 393,721 rays with zero position, distance, coverage or probe changes.
- Final change 2: maximum channel difference **4/65,535**, RMS **0.03302**;
  133 pixels differ by more than one unit, none by more than four.
- Endpoint bounds, forced-inline endpoint bounds and projected bounds with
  automatic inlining: bit-identical to the accepted shader output in the full
  framebuffer comparisons.
- The expanded bounds oracle checks 191,484 sphere-surviving rays, including
  near-axis and near-horizontal directions. It rejects 24,111 and eliminates
  589,032 distance evaluations, with zero lost AA coverage or independent
  double-precision oracle violations. Those corpus counts do not predict the
  default animation's peak-time improvement.

Native effect, SDF, shading, scan, effect-smoke, stack/arena and test-registration
checks passed. The full native run's three findings were corrected and all ten
affected checks passed on rerun. The unrelated Clang replay-form check retains
its expected skip. All full-roster firmware size/layout gates passed for the
accepted changes; the budget-policy suite passed 214 tests.

## Codegen and ITCM

| Accepted image | Full-roster ITCM code | Padding to 196,608 B | DTCM local reserve |
|---|--:|--:|--:|
| Baseline | 192,072 B | 4,536 B | 12,800 B |
| Change 1 | 192,184 B | 4,424 B | 12,800 B |
| Changes 1 + 2 | 193,560 B | 3,048 B | 12,800 B |

The accepted changes cost **1,488 bytes of ITCM** and consume no additional
32 KiB bank. RAM2 remains 520,064 bytes used, with 4,224 bytes free. Per the
requested policy, the extra 3,072-byte software padding reserve is removed;
the stack-floor-derived bank ceiling remains enforced. A regression test
accepts exactly 196,608 code bytes and rejects 196,609.

The final global-O3 reference peaks at 56.044 ms on COM4 versus 56.071 ms for
shipping selective O3 on COM3. The 27-microsecond difference is too small to
interpret as a meaningful advantage from these separate-board captures.

The Cortex-M7 still executes hardware `vsqrt.f32` and `vdiv.f32`; the existing
radial-root common-subexpression elimination was already effective. The wins
come from eliminating repeated extent arithmetic, duplicate SDF evaluation,
twist-loop control and shader call/alpha overhead, not replacing hardware FP
with approximations.

Bounds changed GCC's full-roster inlining decisions: the primary SDF became
outlined even though the smaller single-effect profile still inlined it.
Forcing the two wrapper boundaries inline removed five static distance call
sites, reduced the full-roster pixel stack frame from 304 to 240 bytes and cost
432 ITCM bytes. **It still regressed measured peak time.** Fewer calls and a
smaller stack frame are insufficient reasons to ship a change.

The cheaper projected-cylinder formulation also failed the timing gate. It
replaces per-ray slab intersection points with homogeneous squared tests:
`u = x*dz-z*dx`, `v = (dx²+dz²)*y-dy*(x*dx+z*dz)`. The nearest and farthest
slab radii follow from `dy²*u² + (abs(v) ± h*(dx²+dz²))²`, avoiding reciprocals
and endpoint reconstruction. Despite preserving coverage and removing work,
the complete compiled renderer remained slower at its peak.

## Rejected numerical changes

Reconstructing positions from scalar `t` changed ownership between nearly tied
SDF minima. A distance change around `3e-8` could select a shading point 0.0161
world units away. The original incremental coordinates were retained, and a
regression test now pins the observed tie.

Hoisting lighting into object space changed shared quaternion/transform
codegen. A correction that fixed one case produced another color outlier in
the full corpus. The final shader preserves original world-space lighting;
its maximum difference is four U16 units across the entire corpus.

The archived bounds patch includes its tests and both formulas' capture
evidence. It is available for future work without imposing a measured
regression on the shipping effect.
