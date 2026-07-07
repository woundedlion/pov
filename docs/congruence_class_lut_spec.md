# Congruence-Class Canonical Distance LUTs — Design Spec

Status: FACILITY ONLY — the machinery is landed, tested, and gate-green, but
NO EFFECT USES IT. IslamicStars, the intended Phase-1 consumer, was unwired:
its meshes deform most frames (ripple + warping segues), and under
deformation the canonical premise loses every way it can be played — served
interiors mis-shade, faces switching to the exact path pop visibly, and a
widened sign guard makes probes pay the lookup AND the walk (measured slower
than no LUT during ripple). See §11-12. Available for a future consumer whose
meshes hold still between spawns.
Fourth structural lever of the 2026-07-01 rasterizer campaign, and the only
one of the last four to survive its go/no-go census. Predecessors: per-face
LUT removal (`6241a24b`), convex half-plane path (`7b4873d9`), small-face atan
drop (`bae893c4`); rejected: row-span tightening (two variants), inverse
sampling (see `inverse_sampling_spec.md` for that tombstone).

## 1. The measured fact this exploits

Every islamic mesh's faces are exact copies — within **0.25 px** Procrustes
residual at 288 wide — of a handful of canonical 2D shapes. Census over all
23 registry meshes (gnomonic projection about each face's own centroid,
alignment over cyclic vertex offsets × reflection × optimal rotation):

| Mesh (representative rows) | F | topo classes | geo classes | faces in shared geo classes |
|---|---|---|---|---|
| dodecahedron_ambo_bevel33_relax_hk66 | 362 | 6 | 8 | 100% |
| truncatedIcosidodecahedron_bevel5_relax_hk77 | 722 | 10 | 14 | 100% |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk43 | 1082 | 13 | 19 | 100% |
| truncatedOctahedron_gyro_kis_hk17 | 542 | 24 | 24 | 100% |
| icosahedron_ambo_truncate033_hankin59 | 182 | 4 | 5 | 100% |
| **range over all 23** | 74–1082 | 3–24 | **3–24** | **100% (every mesh)** |

Two supporting facts:

- The *topology* classes (`classify_faces_by_topology`: vertex count + sorted
  1°-bucketed angles + neighbor hashes) are nearly congruence classes already,
  but a topology class can contain 2–3 distinct congruent families (mirror
  pairs, separate symmetry orbits) with up to 6 px cross-family deviation —
  the clustering must be geometric, seeded per topology class.
- Gnomonic projection about a face's own centroid is position-covariant:
  congruent spherical polygons project to identical 2D polygons up to an
  in-plane rotation (+ reflection), for any mesh orientation. Congruence is
  therefore frame-invariant and can be baked once per spawn.

## 2. Why this resurrects the deleted LUT

`6241a24b` removed the per-face distance LUT because its **build** ran per
face per frame (~8 ms of IslamicStars' 21 ms WASM frame) while serving only
~38% of samples at its coarse adaptive resolution. Both defects invert at
class scope:

- **Build cost**: ~3–24 LUT builds per *spawn* (every ~160 frames) instead of
  ~360–720 per *frame*. Effectively free, regardless of resolution.
- **Hit rate**: free build affords fine grids (n = 64–96 per axis vs the old
  adaptive ~12), shrinking the sign-unsafe fallback band (one cell diagonal)
  to the true AA fringe. Target: exact-fallback share ≤ ~20% of probes
  (bake-time computable; a gate, see §8).

The payoff concentrates exactly where current cost is: the concave exact walk
(~12 flops/edge + sqrt, faces of 5–10 edges) becomes a ~20-flop bilinear
lookup. Measured concave probe shares: 15% (cube_relax…expand5) to **100%**
(truncatedIcosahedron…hankin59). Convex faces keep the half-plane path
(`7b4873d9`) — already ~as cheap as a lookup; do not LUT them in Phase 1.

## 3. Data model

### Baked once per spawned mesh (persistent arena, effect-owned)

```
MeshClassBake
  n_classes                       (≤ 32; census max 24)
  per class:
    lut          n×n signed distances, canonical gnomonic plane units
    n            per-class resolution (§6 policy)
    cx, cy       canonical bounding-box center
    Rx, Ry       half-extents (+ margin)
    inv_step_x/y, safe_dist (cell diagonal — sign-pure interpolation guard)
  per face:
    class_id     u8
    vert_offset  u8   cyclic offset aligning mesh vertex order to canonical
    reflected    u8   (bool) mirror family
```

Per-face records ride alongside the existing `MeshState::topology` array in
the effect's persistent slot. `MeshState` itself is untouched.

### Canonical shape and rep selection

Cluster greedily within each topology class (census algorithm): assign a face
to the first sub-class rep within `CONGRUENCE_EPS_PX = 0.25` px residual, else
it founds a new sub-class. Rep = the founding face's centered 2D polygon;
optionally refine to the aligned class mean (halves worst deviation — do this
if the visual gate is marginal, not before). A face whose sub-class stays
size 1 gets `class_id = NO_CLASS` and keeps today's per-face path — the
system degrades gracefully to the status quo, never depends on clustering
succeeding.

### LUT contents

Same field as the deleted `build_distance_lut` (git history at `6241a24b^`):
signed point-to-polygon distance over the canonical polygon's bounding box +
margin, sign from the crossing test, computed with the exact edge walk at
bake. Margin = `BOUNDS_MARGIN_WIDE` as before so the lookup domain covers the
`max_dist_sq` cull ring.

## 4. Per-frame work (unchanged + one small addition)

The forward scan is untouched: per frame, per face, the `Face` ctor still
builds the tangent frame, `poly_2d`, bounds, intervals, packed edges, and
half-planes from the **actual transformed vertices** — bounds and the exact
fallback therefore always reflect the true (rippled) geometry.

Added per face: the alignment rotation between the current projection and
the canonical frame. With the baked `(vert_offset, reflected)`
correspondence, this is one complex correlation over the face's vertices:

```
m  = mean(poly_2d)                      # 2D centroid of current projection
r  = Σ_k  canon[k] * conj(±(poly_2d[(k+off)%n] − m))   # conj pair per refl
rot = r / |r|                           # (cos, sin) stored on the Face
```

~4 flops/vertex, ~25 flops/face/frame — noise next to the ctor's existing
work. Guard: `|r| < kAlignMinCorr` (degenerate correlation) ⇒ treat the face
as `NO_CLASS` this frame. Rotational-symmetry ambiguity (regular n-gons admit
n equivalent rotations) is harmless: the LUT is invariant under the shape's
own symmetry group, so any group element aligns correctly.

## 5. Per-probe path (the hot loop)

`Face::distance` gains a third branch ahead of the existing two:

```
project p → (px, py)                    # unchanged: 2 dots + culls
if (class_lut) {
    q = rot⁻¹ · ((px,py) − m)           # 4 mults + 2 adds (+x-flip if reflected)
    4-tap fetch; if same-sign && min|d| > safe_dist:
        d = bilinear(q)                 # ~15 flops, LUT hit
    else:
        d = plane_dist_convex/exact()   # AA band + sign-unsafe cells: exact,
                                        # on the TRUE per-frame edges
} else existing convex/exact selection
```

Constraints carried over from this campaign, non-negotiable:

- **Both edge loops stay `noinline`** (register-spill cliff, measured +40% on
  every probe when either inlines into the scan chain — `7b4873d9`). The LUT
  branch is loop-free and may inline, but A/B it both ways; if the inline
  variant regresses the no-LUT meshes, make the whole hybrid body `noinline`.
- The `linear_dist` small-face gate (`bae893c4`) applies to the LUT output
  identically — canonical LUT distances are plane units already; large faces
  (inradius ≥ 0.2) still convert via `fast_atan2`.
- Resurrect the `lut_hits` scan-metrics counter (removed in `6241a24b`) for
  the hit-rate gate.

## 6. Resolution and memory policy

Per-class n chosen so the sign-unsafe band (cell diagonal) is a small
fraction of the face: `n = clamp(ceil(k · bb_extent / target_diag), 32, 96)`
with `target_diag ≈ 0.35 × pixel_width` in plane units. At typical face
bb ~0.15 tan units that lands n ≈ 64–96.

| Config | per class | typical mesh (5–14 classes) | worst (24 classes) |
|---|---|---|---|
| f32, n=96 | 36 KB | 180–500 KB | 865 KB |
| f32, n=64 | 16 KB | 80–230 KB | 393 KB |
| **int16, n=64** (quantize d over ±(Rx+Ry); step ≈ 1e-5) | **8 KB** | **40–115 KB** | 197 KB |
| int16, n=48 | 4.6 KB | 23–65 KB | 110 KB |

- **Host/WASM** (8 MB arena): f32 n=96 everywhere, no policy needed.
- **Device** (136 KB persistent partition, 49.5 KB measured headroom after
  `2ef8f081`; both carousel slots need bakes during a crossfade): int16 n=48
  fits typical meshes (~23–65 KB/mesh is over for the big ones ×2 slots).
  Device policy: fixed per-mesh LUT budget (e.g. 20 KB/slot); allocate
  classes by descending `faces × concave` benefit until the budget is spent;
  remaining classes run `NO_CLASS`. Graceful, measurable, no cull required.
  (Quantization error ~1e-5 plane units ≈ 5e-4 px — irrelevant.)
- Bit-identity note: int16-on-device vs f32-on-host would fork sim/device
  output. Use the **same format on both targets** (int16) unless the visual
  gate rejects it; only then consider a ledger entry.

## 7. Ripple (and morphing effects)

The division of labor keeps ripple correctness simple:

- **Edges are exact under ripple.** Bounds, intervals, packed edges, and the
  AA-band/sign-unsafe fallback all come from the per-frame transformed
  vertices, exactly as today. Silhouettes and edge AA do not change class.
- **Interior gradient reads the canonical (undeformed) shape.** The error is
  the intra-face ripple bend: field wavelength ~0.7 rad (~30 px) vs face
  ~7 px, amplitude ≤ 0.15 → sub-pixel interior deviation, same class of
  approximation as the three landed changes. Additionally bounded by the
  0.25 px congruence epsilon when unrippled.
- The alignment rotation is computed from the *rippled* `poly_2d`, so the
  canonical shape is placed at the face's current least-squares pose — the
  interior gradient follows the face's rigid motion through the ripple.

Effects whose meshes change shape per frame (HankinSolids' angle sweep,
ShapeShifter/MeshMorph) must not reuse a spawn-time clustering: they simply
don't pass a bake (null ⇒ status quo). IslamicStars is the Phase-1 consumer.
HankinSolids could rebake per morph target later (its sweep changes vertex
positions every frame, so per-frame congruence would need re-validation —
out of scope).

## 8. Gates (in order; each blocks the next)

1. **Bake telemetry** (log at spawn): classes found, % faces in shared
   classes (expect 100%), worst residual (expect < 0.25 px), predicted
   LUT-hit share from cell geometry. Gate: hit share ≥ ~75%.
2. **ctest 45/45** — the distance-oracle test gains a canonical-LUT case:
   LUT-path samples within `|d| > safe_dist` must match the exact oracle to
   the interpolation bound; the sign must always match (reuse the invariants
   of the deleted `check_face_lut`, git history).
3. **Runtime hit-rate**: `lut_hits / exact_hits` on the split bench across
   solids 0/1/4/6 — confirms the bake-time prediction on real scans.
4. **Visual delta** vs master (400-frame dump): budget = the accepted
   envelope of this campaign (worst-frame mean ≤ ~0.2% FS expected — interior
   deviation ≤ 0.25 px against palette slope); worst frame eyeballed at 1×.
5. **Interleaved same-session A/B**, native + WASM, IslamicStars +
   HankinSolids + ShapeShifter control. Never compare across sessions
   (machine-state swings 30%+, measured). Expected: native scan −15–25%,
   WASM −10–20%, concentrated on concave-heavy meshes; solid 6
   (100% concave) is the bellwether. Gate: no regression on any mesh.
6. Teensy: compile + size gates via pre-commit as usual; runtime remains
   manual (flag in the commit message, as for the previous three).

## 9. Implementation plan

1. `core/sdf.h`: resurrect `build_distance_lut` (from `6241a24b^`) as a
   free function over a centered canonical polygon at fixed n; add the
   hybrid branch + members (`const ClassLut *`, `rot`, `mean`, `reflected`)
   to `Face`; `lut_hits` metric back into `platform.h`.
2. `core/mesh.h` (or a new `core/mesh_classes.h`): the Procrustes clustering
   (port of the census code, arena-based, no std::vector), `MeshClassBake`
   build, per-face record table.
3. `core/scan.h`: `Scan::Mesh::draw` takes an optional
   `const MeshClassBake *` (default null = today's behavior); when present,
   computes the per-face alignment after the `Face` ctor and binds the LUT.
4. `effects/IslamicStars.h`: bake per slot in `spawn_shape` (after
   `classify_faces_by_topology`; rebake both slots after
   `compact_keep_front`, which resets the persistent arena the bakes live
   in); pass the slot's bake to `draw`.
5. Tests: clustering unit test (census invariants: 100% coverage, ≤32
   classes, residual < ε on a couple of registry meshes), oracle extension
   (§8.2), wired into the existing sdf/mesh_raster modules — no orphan files.

## 10. Risks and open questions

- **Spill cliff**: the third distance branch is the same codegen hazard that
  bit `7b4873d9`'s first attempt. Mitigation is baked into §5 and gate §8.5
  (per-mesh regressions, especially the no-bake path and no-LUT meshes).
- **Arena lifecycle**: bakes live in the persistent arena and die at every
  `compact_keep_front`; the rebake-both-slots step must be unconditional or
  the fading front mesh silently degrades to `NO_CLASS` (correct but slower —
  acceptable as the failure mode, assert-logged).
- **Class-count cap**: u8 ids cap classes at 255; census max is 24. Trap at
  bake if exceeded (fail-fast, per project philosophy).
- **Canonical-mean refinement** (§3) is held in reserve for the visual gate.
- **Device budget policy** (§6): fixed LUT budget per slot vs re-running the
  scratch high-water exercise for more headroom — decide when device
  measurements exist at all (all four landed/planned rasterizer changes are
  still unverified on hardware).
- Whether HankinSolids ever joins (per-morph rebake) — deferred until the
  IslamicStars numbers are in.

## 11. Implementation results (2026-07-01)

Landed as `core/mesh_classes.h` (clustering + bake) with the hybrid branch in
`SDF::Face::distance`, the per-frame alignment bound by `Scan::Mesh::draw`,
and per-slot bakes in IslamicStars (rebaked unconditionally after every
`compact_keep_front`). Deviations from the design above, all forced by
measurement:

- **§6 is wrong about WASM**: the WASM release build uses the device
  330 KB `GLOBAL_ARENA_SIZE` (the 8 MB figure is `HS_TEST_BUILD` only), so
  the "f32 n=96, no policy" host row applies nowhere that ships. The device
  policy is the only policy: **int16, n adaptive 32-64 (degrading to fit
  before dropping a class), 18 KB/mesh budget on every target** — 20 KB
  overflowed the 136 KB device persistent partition by 2 KB once the
  per-face records and both slots' bakes were counted (now pinned by
  test_solids' persistent-budget sweep: peak 137.2 KB / 139.3 KB).
- A class whose bake-predicted hit share is under 40% is not kept
  (small faces are mostly fallback band; the LUT would cost its guard on
  every probe and then walk anyway). No registry mesh currently trips this.
- **§7's fixed one-cell fallback band is not ripple-safe.** The ripple slides
  vertices tangentially with a steep Ricker slope, so a face can deviate from
  its rigidly-aligned canonical shape by well over one cell diagonal; the
  canonical field then serves wrong-SIGNED distances near the true edges and
  faces visibly separate. Since the true face is exactly the polygon of the
  displaced vertices, |d_true − d_canon| is bounded by the worst aligned
  vertex deviation — so bind_class_lut measures that deviation per face per
  frame (one extra pass over the vertices) and widens the sign-purity guard
  by it; a face bent beyond ONE cell diagonal keeps the exact path. The tight
  cap is a value-accuracy bound, not just a sign bound: a sign-safe but bent
  face still shades a visibly misplaced interior gradient (stars swim under
  ripple), so a rippling face goes exact until the wavefront passes. Pinned
  by the rippled render A/B, whose delta envelope now equals the static one
  (max 4% FS at a deliberately steep distance-encoding shader; was 43% FS as
  coverage cracks, then 26% FS as gradient shift at a 6-diagonal cap).

Gate outcomes:

1. **Bake telemetry**: every registry mesh 100% shared classes, 3-24 classes,
   worst residual ≤ 0.008 px (vs 0.25 px epsilon); predicted hit share on
   LUT-bound faces 75-91%. The runtime share runs ~2x below the prediction
   because probes concentrate near faces, where the fallback cells are — the
   predictor is a ranking/quality tool, not a calibrated forecast.
2. **ctest 45/45**, including the resurrected LUT-vs-oracle invariants
   (sign always correct, magnitude floor, interpolation bound) across
   identity / cyclic-offset+rotation / mirror bindings, census invariants on
   registry meshes, and a rendered A/B (delta mean 0.02% FS, max 1.9% FS).
3. **Runtime hit rate**: `lut_hits/(lut_hits+exact_hits)` 10-39% by mesh,
   proportional to concave share x coverage as expected.
4. **Interleaved same-session A/B** (Scan::Mesh::draw, 288x144, white
   shader), all 24 islamic meshes:
   - **WASM: −10% to −28% on every LUT-bound mesh** (bellwether
     icosahedron_ambo_truncate033_hankin59: −21 to −29% across runs); the
     only positive deltas occur on meshes with zero LUTs (all-convex
     gyro/dual recipes, where the LUT path never executes) and flip sign
     between runs — noise, not regression.
   - **Native x86: −2% to −8%**, far below the −15-25% §8 expectation: the
     out-of-order host pipelines the exact walk well enough that the lookup
     wins little. The WASM/device targets are where the walk actually hurts,
     and WASM confirms the design's prediction band.
5. **Teensy**: compile + size gates in CI; runtime on hardware still
   unverified (as for the previous three rasterizer changes). DTCM-resident
   int16 LUTs + in-order VSQRT should favor the LUT even more than WASM.

## 12. Why IslamicStars was unwired (2026-07-02)

The gates above were all run on static meshes. In the app, IslamicStars'
meshes deform most frames — ripple bursts recur every 96 frames with a broad
Ricker skirt, and the Segue carousel adds warping transitions — and under
deformation the canonical premise loses every way it can be played:

- **Serve a bent interior** and the gradient mis-shades: the star's shading
  swims against its true edges (caught by eye in daydream at a 6-diagonal
  deviation cap, ~2 px of value error).
- **Fall back per face** and the switch itself pops: canonical and exact
  fields differ by the interpolation envelope, so a face crossing the
  deviation threshold visibly jumps between frames.
- **Widen the sign guard** to stay crack-free and the band swallows a
  concave star's interior: probes pay the 4-tap lookup and then walk anyway
  — ripple frames measured ~3x quiescent, slower than no LUT at all.

A sign-only far-field mode (serve only the cull-ring skirt beyond the AA
pad, exact interior) is visually lossless and keeps most of the win, but it
still leaves the mode-transition pop on the interior and adds a third
serving regime to the hot path — for an effect that deforms most frames the
complexity is not worth the residual win. The facility remains for effects
whose meshes hold still between spawns; the deformation guard
(SDF::ALIGN_MAX_DEV_DIAGS = 0.25 cell diagonals) stays as the safety net
that keeps any future consumer's deformed frames correct rather than fast.
