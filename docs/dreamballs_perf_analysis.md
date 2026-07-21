# DreamBalls deep perf analysis — device phase split, M7 codegen audit, and lever outcomes (2026-07-20)

Companion to the 2026-07-20 profile pair
(`docs/profiles/shipping|O3/profile_dreamballs_teensy_2026-07-20.md`), which
concluded from an -Os-vs-O3 twin that DreamBalls is coverage-bound and "the way
to 16 fps is a coverage cut, not codegen". That is true of the **opt level** and
false of the **algorithm**: a deep phase split showed ~46% of `db_mesh_plot`
never deposits a pixel. Two levers landed off that finding; five were killed or
blocked, each with a measurement.

## New instrumentation (landed `b3b6ca24`)

`HS_PROFILE_DEEP` sub-scopes in `core/render/plot.h` (the facility existed in
scan.h/filter.h; plot.h had none): `plot_edge_presample`, `plot_seg_cull`,
`plot_seg_sim`, `plot_seg_draw`. Zero cost unless `HS_PROFILE_DEEP_ENABLE`.

Capture: `HS_PROFILE_DEEP=1 HS_TEENSY_PORT=COM3 bash tools/profile_one.sh
DreamBalls profile 100 16 -D HS_PROFILE_EPOCH_REVS=1600`.

## The finding: half of the rasterizer's time is pre-pixel

Pre-change device split, worst window (`db_mesh_plot` 93.56 ms/f, sprite
overlap, 24 copies):

| phase | ms/f | share | calls/f | per call |
|---|--:|--:|--:|--:|
| `plot_seg_draw` | 42.87 | 46% | 3,979 | 10.8 µs |
| `plot_seg_cull` | 22.33 | 24% | 15,840 | **846 cyc** |
| `plot_edge_presample` | 14.81 | 16% | 1,584 | **5.6 kcyc/edge** |
| `plot_seg_sim` | 8.90 | 10% | 3,979 | 2.2 µs |

Every mesh edge was chopped into 10 sub-segments by `Line::sample(…, 10)`, and
each sub-segment then paid a full `edge_visible_in_clip` that re-derived
endpoint rows/columns (acos/atan2) its neighbour had just computed. **Only 25%
of sub-segments survive** (3,979 of 15,840), yet all of them — and the whole
presample — ran before anything was culled.

## M7 codegen audit (profile image, gcc 15.2.1, `-Os` + `HS_O3` regions)

Per **drawn pixel** the O3 rasterizer executed: `fast_sinf`+`fast_cosf` (a
`vdiv` each, Bhaskara rational), `.normalized()` (vsqrt+vdiv), a full 56-byte
`Fragment` materialized on the stack, an indirect `blx` to the fragment shader,
`bl vector_to_pixel` (**not** inlined despite `HS_O3_FN` — the optimize-attribute
inline seam), and `bl AntiAlias::plot`. Inside `vector_to_pixel` sat
**`bl __fmodf_veneer`** — a newlib FLASH call, once per plotted pixel.
`AntiAlias::plot` itself is clean (no calls, no div/sqrt, `vrintm` floors,
unrolled integer blends) — it is the real work. No D-cache thrash evidence.

## Landed

### `985bbfa4` — drop the per-pixel `fmodf` (was lever 3)

`vector_to_theta` wrapped `fast_atan2`'s output through the generic float
`wrap()`, whose `std::fmod` is an out-of-line FLASH call. `fast_atan2` is
bounded by |π|, so the modulus is the identity and the wrap is one conditional
add plus the existing half-open guard.

- **Bit-exact**: 42 M inputs per width (288/144/480), 0 bit differences,
  covering ±0.0, denormals, ±π, ±W/2 and the whole 2 M-float band below zero
  where `t + W` can round up to exactly W.
- Codegen: `bl __fmodf_veneer` 1 → 0 in `vector_to_pixel` (10 → 8 image-wide),
  body 360 → 336 B.
- **Device A/B** (same base, deep capture, worst window): `db_mesh_plot`
  93.56 → 88.77 ms/f (**−5.1%**), peak frame 96.20 → 91.10 ms, spill
  18.4% → 14.9%. Clean control: the scopes that use the projection moved
  (`plot_seg_draw` −7.2%, `plot_seg_cull` −9.2%) while those that do not were
  flat (`plot_edge_presample` +0.1%, `plot_seg_sim` −0.3%).

### `c25667ff` — gate wireframe edges before the presample (was lever 1)

Rejects the whole edge first, then hoists the per-sub-segment cull:

- **Whole-edge reject before `Line::sample`.** Every sub-segment lies on the
  u→v great circle (`Line::sample` rotates about that axis, and
  `make_geodesic_edge_span` resolves an antipodal pair to the *same* axis via
  `stable_perpendicular_axis(a)` under an identical threshold), so the edge's
  own span bounds the chopped polyline exactly — no bulge approximation needed.
  A reject now skips the presample too.
- **Hoisted bits for survivors** via the existing `gate_trail_edges`, whose
  per-edge verdicts are identical to `rasterize`'s own cull, fed through
  `rasterize`'s existing `edge_visible` parameter.
- Both paths are `if constexpr`-gated on `pipeline_hoistable_cull` (world-cull
  pipelines keep the old behavior — culling them on raw points would be wrong)
  and opt out under a vertex shader, which may move samples off the arc.
- `gate_trail_edges` moved above `struct Mesh` to be visible at the call site.
- **Test**: new `test_mesh_edge_gate_pixel_parity` (randomized orientations,
  4 clip topologies incl. a seam-adjacent wedge) asserts a clipped
  `Plot::Mesh::draw` is pixel-identical to the full render inside the band.
  Wired into `unit_plot_scan`. **Negative control**: injecting an over-cull
  (drop every 7th edge) makes it fail, so the test has teeth.
- Cost: ~2.1 KB ITCM (phantasm slack 5,240 → 2,104 B; gate still PASSes).
- **Device A/B**: see the results table below.

### `161a8188` — divide-free reciprocal square roots in `screen_step` (was lever 6)

Both reciprocal square roots go through `fast_rsqrt`; the speed floor moves
into the squared domain (sqrt is monotone, so it is the same clamp) to keep
that function's strictly-positive domain. `Plot::rasterize` drops
**28 → 16 `vdiv.f32`** and **15 → 7 `vsqrt.f32`**.

- **Device A/B** (stacked on the two above): `plot_seg_sim` — the phase
  `screen_step` dominates — 8.89 → 8.19 ms/f (**−7.9%**); worst frame
  73.95 → 73.19 ms (−1.0%); `db_mesh_plot` flat within noise (+0.3%), since
  the step simulation is only ~12% of it.
- **It gives ITCM back**: phantasm RAM1 code 194,504 → 193,000 B, restoring
  headroom 2,104 → 3,608 B.

## Killed or blocked, with the measurement

- **`screen_step` via `fast_rsqrt` (lever 6) — LANDED `161a8188`.** An earlier
  pass in this session called it dead on a byte-identical disassembly; that
  comparison was reading a **stale build** (the worktree's first `pio run`
  had failed transiently and the object was not rebuilt). A clean rebuild
  with `CCACHE_DISABLE=1` shows the real effect: `Plot::rasterize`
  28 → 16 `vdiv.f32` and 15 → 7 `vsqrt.f32` for +72 B. A standalone compile
  at the image's own flags confirms `fast_rsqrt` emits the bit-hack + two
  Newton steps with no divide or square root. **Method note:** the magic
  constant lands in a literal pool printed as `.word 1597463007`, so
  searching a disassembly for `5f3759df` (or for `movw`/`movt` immediates)
  finds nothing and proves nothing — grep the little-endian bytes
  (`df59375f`) or the decimal.
- **Rotor recurrence for `sample_geodesic` (lever 5) — NOT VIABLE.** The draw
  loop advances by *adaptive* steps, so the rotor needs `sin δ`/`cos δ` for a
  varying δ each step — exactly the two trig calls it was meant to replace.
  The only residue is dropping the per-plot `.normalized()` (~30 of ~1,100
  cyc/plot), and that renormalize exists precisely because fast_sinf's 0.17%
  error breaks c²+s²=1. Not worth the unit-length risk.
- **`always_inline` on `vector_to_pixel` (lever 4) — HELD; priors say ~1–2%.**
  It works, at **+6,320 B** ITCM in the profile image and **+1,072 B** in
  phantasm. The tempting precedent is HopfFibration `d9bd43da` (always_inline
  on O3-region leaf callees: worst frame 71.3 → 56.8 ms, 1.25×; an
  arithmetic-only subset +1,312 B for 1.16×) — but that landing's own set
  already includes `fast_acos/atan2/sinf/cosf`, `normalized`, `rotate` and
  **`vector_to_theta`**, and it records that the other region effects
  (DreamBalls included) *inherited* that leaf inlining. The expensive interior
  is therefore already inlined here; lever 4 adds only the wrapper's own
  `bl`/`ret` and its two-float struct return. The global -O3 twin — which does
  inline it, no option mismatch — buys just 1.067× overall, part of that in
  `isr_pack` and cheap holds. Revisit against the recovered 3,608 B only with
  a device A/B in hand.
- **Per-plot shader call barrier (lever 2) — ITCM-BLOCKED, ~2–5% at best.**
  Templating `rasterize` on the shader type duplicates it per (pipeline,
  shader) pair; the phantasm image already carries **18**
  `Plot::rasterize`-family instantiations (largest 5,756 B). The governing
  precedent is MindSplatter's selective-O3 region over `ParticleSystem::draw`,
  **built, measured and parked** at +6,096 B (81% of the roster's headroom) for
  ≈−6.5 ms because it crossed no cadence tier. Sizing the gain: ~18–20k
  plots/frame at ~1.3 kcyc each, so 100 cyc/plot ≈ 3.3 ms ≈ 4% — and the
  column-cull ledger's own estimate lesson (a share eyeballed at ~15 ms
  measured ~3.5) says discount structural estimates of this kind. The
  batch-shade variant avoids the duplication but needs a new seam through
  every primitive.
- **Reusing the edge span in the strategy (lever 7) — NOT WORTH IT.**
  `angle_between` + axis derivation per surviving sub-segment ≈ 90 cyc ×
  3,979/f ≈ 0.6 ms of ~89 ms (**0.7%**), and plumbing spans from the gate into
  `rasterize` needs a per-edge span array plus coupling between two modules.

## Results

Each column is a separate deep capture from a worktree at that exact commit,
same settings, `validate` PASS on all four.

| | baseline `14f4ed3e` | +fmodf kill `985bbfa4` | +edge gate `c25667ff` | +rsqrt `161a8188` |
|---|--:|--:|--:|--:|
| `db_mesh_plot` worst | 93.56 ms/f | 88.77 (−5.1%) | **69.70 (−21.5%)** | 69.88 (noise) |
| `plot_seg_cull` | 22.33 | 20.30 | **0.20** | 0.22 |
| `plot_edge_presample` | 14.81 | 14.81 | **5.44** | 5.70 |
| `plot_seg_sim` | 8.90 | 8.87 | 8.89 | **8.19 (−7.9%)** |
| `plot_seg_draw` | 42.87 | 39.80 | 39.84 | 39.77 |
| peak frame render | 96.20 ms | 91.10 | 73.95 | **73.19** |
| spilled frames | 244/1328 (18.4%) | 200/1344 (14.9%) | 100/1456 (6.9%) | 99/1456 (6.8%) |
| phantasm ITCM slack | 5,240 B | 5,240 B | 2,104 B | 3,608 B |

**Cumulative: peak frame 96.20 → 73.19 ms (−23.9%), spill 18.4% → 6.8%, for
1,632 B of net ITCM.** `plot_seg_draw` moved only with the fmodf kill and was
otherwise flat, which is the control: the edge gate removed non-drawing work
only, exactly as designed.

## The next lever, per the ledgers

MindSplatter's top-ranked and **landed** lever was neither of the two above: it
was killing a double projection — the gate already computes `rows[]`/`cols[]`
by formulas bit-identical to `vector_to_pixel`, and the raster path recomputes
them per plot. Routing the gate's values through so a dot skips
`vector_to_pixel` was worth **1.14×** at negligible ITCM. The edge gate landed
here computes exactly those values per presampled point, so the same redundancy
now exists in DreamBalls. Caveat before believing the number transfers:
MindSplatter banked it on its single-dot fast path (91% of its edges), whereas
DreamBalls' long shell edges sub-step to ~5–6 plots each, so the reuse applies
only to dot-path edges and to sub-segment endpoints. Measure the plot/dot split
first (`HS_SCAN_METRIC` counter run, as MindSplatter did).

## Still the binding constraint

The 62.5 ms window leaves ≈58.5 ms for render after the 6.35% ISR share. The
32-frame two-sprite overlap (`SPRITE_LIFE` 320 > `SPAWN_PERIOD` 288) draws two
whole shell sets and is what the worst windows measure; preset 0's own hold
(18 copies × radius 0.30) is the other red. Code work has taken real time out
of both, but closing the overlap outright is still a coverage decision —
shortening the crossfade or trimming preset 0's copy count — not a codegen one.

## Validation protocol used

Same-base worktrees per configuration (`14f4ed3e`, `985bbfa4`, `c25667ff`),
deep captures at identical settings, `parse_profile.py … validate` PASS on
each; native suite (51 tests) green per commit; pre-commit Teensy size gate
(3 envs) per commit; disassembly compared with `arm-none-eabi-objdump` from the
matching worktree's own ELF. Host timing was not used to judge any lever —
host -O2 CSEs redundant work device GCC does not.
