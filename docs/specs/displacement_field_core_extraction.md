# DisplacementField core extraction

**Status: IMPLEMENTED (rev 4, 2026-11-12).**
The ring-stack, column-time sweep and swept-bump machinery now live in `core/`.

**Historical baseline:** the problem statement, inventory and original line
references describe master `7a22dc38c` before extraction. Section 4 records
the implemented interfaces; section 9 records qualification results.

**The bar is correctness, performance and beauty, not pixel identity.**
Changed pixels are fine. Strict pixel parity is only a tool here: it proves
that a pure move broke nothing (§6.1). Steps that change behaviour (S3b, S1b)
are judged by property tests and a visual check.

Revision history:

- **Rev 2** followed an outside review:
  - the swept-bump cover gained a validity contract;
  - the pool/alternates binding was specified;
  - the parity harness moved to isolated trees and reuses the flush fixture;
  - the device bar became measurable;
  - the S1b arithmetic direction was corrected.
- **Rev 3** records the owner's decisions on Q1–Q3 (§7):
  - the cover fix becomes the exact bound (S3b);
  - S1b lands without a gate.
- **Rev 4** records the implementation, the value-owned fused late-window
  visitor, the allocation-order callbacks and qualification evidence (§9).

## 1. Problem

`effects/DisplacementField.h` is 2,080 lines. The field primitives
(`BumpParams`, `NoiseProductTransformer`, `DominantFieldAccumulator`), the ring
SDF, the fused stack scan, the chunk culls and the hue table are already in
core. About 600 lines of rendering machinery are still in the effect, and none
of it depends on noise or falling balls:

- **The producer half of a core contract.** `Scan::DistortedRingStack`
  probes its shape range for the optional `knot_range` / `table_ready` /
  `table_built` hooks (`core/render/scan/shapes.h` (baseline lines 232-237), `:458-463`). The
  only implementation of those hooks is DF's private `RingPool::ShapeView`
  (`effects/DisplacementField.h` (baseline lines 1662-1702)). Any second ring-stack client would have to
  rebuild the pool.
- **Column-time sampling for ring stacks.** The two-window swept bake
  (`SweepCache`, `retain_alternate`, `apply_second_window`,
  `scan_swept_stack`, `mark_boundary`) is generic to any stack whose knots
  change over a frame. `mark_boundary` reads `SDF::DistortedRing` internals
  (`POLE_SIN2_FLOOR`, `search_reach_cells`).
- **Duplicates of core helpers.** `knot_fraction` (`:367`) repeats
  `math::half_turns` (`core/math/3dmath.h` (baseline lines 498)). `prepare_sweep`'s step angle
  (`:532-540`) repeats `Spacetime::history_turn`
  (`core/render/spacetime/pose_table.h` (baseline lines 129-131)).
- **A core field with no core consumer for its endpoint state.** `BumpParams`
  keeps `previous_center` / `previous_envelope` (`core/animation/params.h` (baseline lines 1211)).
  The only code that evaluates a bump between those endpoints is DF's
  `BallPath` / `ball_at` (`:323-417`). That code also under-covers fast
  motion (§4.3).

This applies the standing rule that rendering machinery belongs in core and
effects supply only geometry and fields.

### 1.1 Prior decision this supersedes

Monolith decomposition rev 2 kept the ring pool storage-only inside DF, with
shading staying in DF. Rev 8 moved "RingPool → core" to its §11 Reassess list,
and that spec landed with no step moving effect code into core. This spec is
that reassessment. The owner ruled (Q1) that the pool **and** its hue-row
shader move to core.

## 2. Goals and non-goals

**Goals**

- DF keeps only what is specific to it: the noise→balls phase machine, ball
  spawning, the palette roll and spin, `Params` and presets, composing the
  noise and ball fields into knot shifts, and its hue shading policy.
- Every moved symbol has one owner in core, with property tests, death cases
  for its `HS_CHECK`s, and doxygen coverage.
- The swept-bump cover is correct for every input (S3b).
- No device regression beyond measured variation without owner approval, and
  no demotion of other code out of ITCM (§6.2).

**Non-goals**

- No new perf levers. Nothing on the dead-lever list in the DF perf ledger is
  reopened.
- No change to sliders, presets or choreography.
- No change to `Scan::DistortedRingStack`'s behaviour. S4a moves it into
  `core/render/scan/ring_stack.h` verbatim (§4.4).

## 3. Inventory

| Symbol (DF lines) | Destination | Step |
|---|---|---|
| `prepare_sweep` quaternion step → rotation vector (532-540) | `Spacetime::rotation_vector` in `core/render/spacetime/pose_table.h` | S1a |
| `knot_fraction` body (365-374) | `math::half_turns` call | S1b |
| `KnotStep`, `RingKnotWalk` (954-988) | `SDF::KnotStep`, `SDF::KnotWalk` in new `core/render/sdf/ring_knots.h` | S2 |
| Strided Catmull-Rom fill in `bake_noise_octaves` (1223-1233) | `SDF::catmull_rom_strided` in `ring_knots.h` | S2 |
| `BallPath`, `ball_at`, path setup in `draw_swept` (574-595), `FIELD_REACH_CUT`, `FIELD_AMPLITUDE_CUT` | `Animation::BumpPath`, `Animation::BumpCover`, `make_bump_path`, `bump_field_at` in `core/animation/transformer.h` | S3, S3b |
| Cap ∩ ring azimuth-arc solve in `bake_ball_spans` (1326-1348) | `Plot::cap_ring_arc` in `core/render/plot/cull/visibility.h` | S3 |
| `Scan::DistortedRingStack`, `CandidateTable` (`core/render/scan/shapes.h` (baseline lines 125-560)) | new `core/render/scan/ring_stack.h`, verbatim | S4a |
| `RingPool`, `RingPool::ShapeView`, `RingShader` (1630-1840) | `Scan::RingStackPool`, `Scan::HueRowShader` in `core/render/scan/ring_stack.h` | S4 |
| Shift-side `SweepCache` fields, `retain_alternate`, the capacity flush in `draw_pass` (803-810), `scan_swept_stack`, the shift install in `apply_second_window`, `mark_boundary`, `SWEEP_QUARTER`, `SWEEP_WINDOWS`, `retains_endpoint`, `other_endpoint` | `Scan::KnotAlternates`, `Scan::SweptRingStack` in `core/render/scan/ring_stack.h` | S5 |
| `ScratchFootprint` (2038-2050) | `ArenaFootprint` in `core/memory/arena.h` | S6 |

**Stays in DF:**

- the phase machine, spawning, palette, `Params`
- `draw_rings`, `draw_pass`, `draw_swept` (now thin)
- both bakers, `ball_field` / `ball_field_at`, `cache_ball_frame`
- `select_hue_mode`, `HueMode`, `noise_lut_samples`
- the shading-side `SweepCache` fields (`modes`, `colors`, `alphas`)
- the oracle flags `force_window_bakes`, `force_column_reference`, `force_exact_hue`

`ball_field` and `ball_field_at` stay because they index DF's struct-of-arrays
ball cache and are named in the selective-O3 profile image.

## 4. Design

### 4.1 S1: reuse existing helpers

- **S1a.** Add `Spacetime::rotation_vector(const Quaternion &previous, const
  Quaternion &current)` to `pose_table.h`. It returns the shortest rotation
  carrying `previous` to `current`, as `prepare_sweep` computes it today,
  including the `r < 0` sign flip and the small-angle branch. Add
  `Spacetime::rotation_angle(const Quaternion &step)`, the
  `2·asin(min(|v|, 1))` term. Both `history_turn` and `rotation_vector` call
  it. This is a pure move.
- **S1b.** `sweep_fraction` (§4.5) calls `math::half_turns`, leaving one
  implementation of the azimuth-in-half-turns concept.
  - The inline code multiplies by `1/π`; `half_turns` divides by `π`. A sample
    time can move by 1 ulp. Both are equally correct.
  - It lands last, as its own commit, only so the S4+S5 move can still be
    verified by strict parity.

### 4.2 S2: `core/render/sdf/ring_knots.h`

- **`KnotStep { cos_t, sin_t, cos_d, sin_d }`.** Passed by value so it travels
  in VFP registers, as today.
- **`KnotWalk`.** Moved verbatim from `RingKnotWalk`: basis reference, the four
  trig terms, a `position()` / `advance()` angle-addition recurrence starting
  at knot 0. Keep the walker form, not a visitor: as a visitor lambda it
  outlined FastNoiseLite at `-Os` (+448 B, monolith D2).
- **`catmull_rom_strided(value, d, k, n)`.** `value(index)` is a functor
  returning the grid sample. It returns `value(k)` on grid knots and otherwise
  the uniform Catmull-Rom over the four stride-`d` neighbours, with periodic
  wrap over `n`. Use the same arithmetic order as `:1230-1233`. DF's
  alternate-aware `value` lambda stays in DF and is passed in.

### 4.3 S3 and S3b: swept bump and cap-ring arc

**S3 (pure move).**

- **`Animation::BumpPath`.** The `BallPath` fields, unchanged so
  `SCRATCH_A_PEAK_BYTES` does not move: `from`, `to`, `axis_from`, `axis_to`,
  `reach_from`, `reach_to`, `cos_reach`, `cos_cover`.
- **`Animation::BumpCover { Vector centre; float reach; }`.** The cap
  covering the whole swept footprint. `cache_ball_frame` takes it in place of
  today's `mid` and `cover` locals. It is returned separately so per-ball
  scratch does not grow.
- **`make_bump_path(const BumpParams &, const Vector &turn, const Vector
  &axis, BumpPath &path, BumpCover &cover) -> bool`.** Returns `false` when
  both endpoint reaches are at or below `BUMP_REACH_CUT`, and writes neither
  output. That is the case where DF writes `ball_reach = -1` today. On `true`,
  it fills both outputs. In S3 it keeps today's bound (`:591-595`).
- **`bump_field_at(p, bp, path, theta, t)`.** `ball_at` verbatim, `HS_HOT_INLINE`
  so it still inlines into DF's `HS_O3_FN ball_field_at`.
- **Renamed constants.** `BUMP_REACH_CUT` (was `FIELD_REACH_CUT`) and
  `FIELD_AMPLITUDE_CUT` live beside the bump profile. DF's noise-amplitude
  gates refer to the core constant.
- **`Plot::cap_ring_arc(cos_t, sin_t, axial, radial, azimuth, cos_radius,
  n, margin) -> KnotArc { first, count, empty }`.** The padded arc solve from
  `:1326-1348`, including the `b_term ≤ margin` whole-ring case and the
  `q ≥ 1` empty case. `COS_MARGIN` becomes the caller's `margin`.

**Cover contract.** For every `t ∈ [0, 1]`:

- the centre `c(t)` (the normalized lerp of `from` and `to`) lies within
  `cover.reach − reach(t)` of `cover.centre`;
- so the cap of radius `reach(t)` around `c(t)` lies inside the cover.

`cos_cover = cos(min(cover.reach, π))`.

**Why today's bound breaks the contract.** Today's code pads by half the
chord, `h = |to − from| / 2`, plus `1e-4`. The centre path runs along the
minor great arc, and its farthest point from the midpoint is half the endpoint
separation, `α/2 = asin(h)`, which is more than `h`. The `1e-4` pad covers the
gap only while `α/2 − sin(α/2) ≤ 1e-4`, that is `α ≲ 0.169 rad`.

DF exceeds that range today:

- the fall moves `π / duration` per frame, 0.157 rad at the Speed slider's
  maximum (duration 20);
- the stack walk turns up to 0.02 rad per frame (`RandomWalkOptions::speed`);
- worst case together: about 0.177 rad, short by about 1.6e-5 rad.

**S3b (behaviour change): the exact bound.**

```cpp
cover.centre = (from + to).normalized();
cover.reach  = reach + asinf(fminf(h, 1.0f)) + COVER_PAD;  // COVER_PAD = 1e-4f
```

- **Correct for every non-antipodal input.** `asin(h)` is the exact
  half-separation, and `COVER_PAD` stays as float slop for normalization and
  `asinf`.
- **Never smaller than today's bound.** Because `asin(h) ≥ h`, the cover never
  shrinks; it grows by about `h³/6`. At default sliders that is about 6e-6 rad.
- **Cost.** One `asinf` per active ball per swept frame, at most 56.
- **Near-antipodal endpoints.** When `|from + to|² < MID_EPS²` the midpoint is
  undefined. `make_bump_path` then returns the whole-sphere cover:
  `cover.reach = π`, `cos_cover = −1`, `cover.centre = to`. This is correct
  for every input, so it is not a trap. It is unreachable at today's sliders,
  because `BallDrop`'s duration ≥ 2 limits a uniform fall step to `π/2`.

The cover feeds the ring prefilter, the band, the hue-table domain and boundary
marking, so S3b changes ball-phase pixels slightly in every frame. Its bar is
the containment property test (§6.3) plus the visual check (§6.1).

### 4.4 S4a and S4: `core/render/scan/ring_stack.h`

**One file for the ring-stack concept.** `core/render/scan/ring_stack.h` holds everything
that builds, feeds and scans a same-axis `DistortedRing` stack:

- the fused scan (`DistortedRingStack`) and its `CandidateTable`;
- the pool and hue-row shader (S4);
- the column-time alternates and two-window driver (S5).

The stack's optional-hook contract and its only implementation then sit side
by side. `core/render/scan/shapes.h` keeps the single-shape scans.

**S4a: verbatim move.** Move `DistortedRingStack`, including its
`HS_O3_BEGIN` / `HS_O3_END` bracket and every placement attribute, into the
new file in its own commit.

- **Includes.** `core/render/scan/ring_stack.h` includes `core/render/scan/shapes.h` for `check_canvas_dims`,
  `swept_rows`, `window_pieces`, `walk_pieces_once`, `plot_admitted` and
  `source_clip`. The umbrella `core/render/scan.h` includes it right after
  `core/render/scan/shapes.h`, so the single TU sees nearly the same definition order.
- **Callers.** No caller changes; every user reaches the stack through
  `scan.h`.
- **Symbol maps.** Compare the phantasm symbol map with the parent's: the TU
  sits at GCC's inline-growth cap.

**S4: pool and shader.**

```cpp
template <int SLOTS, int ROW> class RingStackPool;  // ROW = W + 1
template <typename Pool> struct HueRowShader;
```

**Binding.**

- **Pool type determines the alternates type.** The pool declares
  `using Alternates = KnotAlternates<ROW - 1>;` (one alternate bit per knot,
  `ROW − 1 = W` knots).
- **Runtime sizes.** Capacity and slot count are runtime values on the
  alternates (§4.5), so they never enter the pool's type.
- **Mutable binding.** The pool binds `Alternates *`, non-const, because
  `ShapeView::table_built()` writes `table_ready`.
- **`begin_frame` checks the slot count.** With alternates bound, it
  `HS_CHECK`s `alternates->slots() ≥ SLOTS`. That runs once per frame, on a
  cold path.

**API.**

- **The pool API stays exactly as it is today:**
  - `init_storage(Arena &)`, `BYTES`
  - `begin_frame(Alternates *alternates = nullptr)`
  - `next()` (stateless pending rows)
  - `commit(ring, frag_alpha, lut_n, ring_args...)` (construct, then count++;
    with alternates bound, it sets the knot band from `knot_range`)
  - `release()`, `size`, `slot_of`, `slot_map`, `shapes`, `hue_row`,
    `lut_columns`, `frag_alpha`, `set_alpha`, `shift_row`, `mutable_hue_row`
  - `CULLED`, plus the `SLOTS ≤ INT8_MAX` static assert
- **`ShapeView`** implements the stack's three optional hooks. That is the
  only coupling between the pool and S5:
  - `knot_range(slot, knot, lo, hi)` sets `lo = hi =` the primary shift from
    pool storage, maps `slot → ring` through `ring_by_slot`, then calls
    `alternates->widen(ring, knot, lo, hi)` when alternates are bound.
  - `table_ready()` / `table_built()` read and write `alternates->table_ready`.
  - Without alternates, `knot_range` returns the primary alone,
    `table_ready()` is `false`, and `table_built()` does nothing. That matches
    today's `cache == nullptr` behaviour.
- **`HueRowShader`** is `RingShader` verbatim: quintic lerp of the slot's hue
  row at `v0 · lut_columns` × `frag_alpha` × `v2`.
- **HS_CHECK messages** change from `"DisplacementField: ring pool …"` to
  `"ring stack pool …"`. The two death cases in
  `tests/death/effects_cases.h` (baseline lines 188-215) move to `tests/death/scan_cases.h`,
  and their `DEATH_CASE_COUNT` is derived from the table.
- **Placement attributes** carry over unchanged: `HS_NOINLINE_NOCLONE
  begin_frame`, `HS_HOT_INLINE next`, `HS_HOT_FLASH_MEMBER commit` and
  `HS_FLASH_MEMBER init_storage`.

### 4.5 S5: `Scan::KnotAlternates` and `Scan::SweptRingStack`

**`KnotAlternates<KNOTS>`** is the shift-side state of `SweepCache`.

- **Template parameter.** It is templated only on knots per ring, with
  `WORDS = (KNOTS + 31) / 32` constexpr as today, so index arithmetic stays
  constant-folded.
- **Members:** `values`, `masks`, `starts`, `boundary`, `ring`, `count`,
  `table_ready`, plus runtime capacity and slot count with read-only accessors
  `capacity()` and `slots()`.
- **Allocation:** `allocate(Arena &values_arena, Arena &table_arena, int slots,
  int capacity)`. It `HS_CHECK`s `capacity ≥ KNOTS`, a cold path that replaces
  DF's `static_assert` (`:291`). It reproduces today's split: `values` from
  scratch B; `masks` (zeroed), `starts` (zeroed) and `boundary` from scratch
  A. The four-argument overload delegates to an overload taking a caller
  storage callback. That callback runs after `starts` and before `boundary`,
  so DF allocates its hue modes, colors and alphas in the original
  `draw_swept` order (`:598-605`). DF passes `RING_SLOTS` and capacity `3840`.
  Capacity is also checked against `UINT16_MAX`, the prefix-storage limit.
- **What stays in DF.** The shading fields DF allocates there move into DF's
  own `SecondWindowShading` struct.
- **`begin_ring(i)`.** Sets `ring` and seeds `starts[i·WORDS] = count`.
- **`fits(lut_n)`.** `count + lut_n ≤ capacity()`.
- **`retain(k, primary, alternate, time, row)`.** `retain_alternate` verbatim.
- **`finish_ring(i)`.** The popcount prefix pass at `:879-887`.
- **`widen(ring, knot, float &lo, float &hi) const`.** Precondition: `lo` and
  `hi` hold the knot's primary shift. If the knot's mask bit is set, it widens
  them by the retained alternate. This is the alternate-lookup half of today's
  `ShapeView::knot_range` (`:1672-1683`); the primary read stays in the pool.
- **`reset()`.** Sets `count = 0` and `table_ready = false`.

**`SweptRingStack`** is a struct of static helpers:

- **Constants:** `QUARTER`, `WINDOWS = 2`, `retains_endpoint(t)`,
  `other_endpoint(t)`, `sweep_fraction(p, middle)`. `sweep_fraction` calls
  `math::half_turns`, then wraps and clamps the result to the window.
- **`mark_arm_boundary<W>(basis, KnotStep, theta, n, displacement,
  thickness, motion, support_cells, uint8_t *out)`.** `mark_boundary`, with the
  thickness slider and the spline stencil (`SUPPORT`) as arguments.
- **`scan_two_windows<W, H>(pipeline, canvas, n_rings, pool, Alternates
  &alternates, table, shader, Reshade reshade)`.** The visitor is passed by
  value. The driver:

  1. draws window 0;
  2. calls `reshade.begin(ring, slot)` for each drawn ring;
  3. installs each masked late shift, then calls `reshade(knot, shift)` for
     every knot in the same loop;
  4. closes the shift row with `row[n] = row[0]`, then calls `reshade.finish()`;
  5. draws window 1, releases the pool and resets the alternates.

  DF's `begin` binds the hue row, hue mode and hue-table policy. Its call
  operator shades the installed shift; `finish` closes the hue row and sets
  the late alpha. The core driver owns the shift installation and closure.

**Fused late update (implemented).** Shift installation and hue evaluation
share one knot pass. The value-owned visitor retains its per-ring state
across `begin`, knot calls and `finish`; its hot methods carry
`HS_HOT_INLINE`. This is the fused fallback of the original design, with one
visitor object carrying both ring setup and knot shading. The final device
comparison in §9 qualifies the implemented driver.

**Test-hook counters.** `duplicate_knots` and `stack_flushes` move onto
`KnotAlternates` behind `HS_ENABLE_TEST_HOOKS`. `primary_knots`,
`primary_bakes` and `alternate_noise_samples` measure the bake, so they stay in
DF. The `DisplacementFieldWhiteBox` accessors read the moved counters through
the effect's alternates.

**Profile scopes.** DF keeps `df_fused_scan`, `df_chunk_cull`, `df_lut_bake`,
`df_octave_noise` and `df_hue_table_prep` at its call sites, so report rows
stay comparable with `docs/profiles/shipping/profile_displacementfield_teensy_2026-11-12.md`.

### 4.6 S6: `ArenaFootprint`

A constexpr accumulator of sequential aligned allocations (`add<T>(count)`,
`bytes`), moved from DF's `ScratchFootprint` into `core/memory/arena.h`. It
keeps the `alignof(T) ≤ alignof(std::max_align_t)` static assert. DF's
`SCRATCH_A_PEAK_BYTES` / `SCRATCH_B_PEAK_BYTES` use it.

`KnotAlternates` also exposes `static constexpr add_footprint(ArenaFootprint
&values, ArenaFootprint &table, int slots, int capacity)`, so DF's budget
asserts list core types instead of re-deriving their layouts. Its optional
callback overload mirrors `allocate`: caller footprint additions occur after
the prefix arrays and before the boundary bytes. The default overload supplies
an empty callback.

DF's persistent `FOOTPRINT_BYTES` sums keep their worst-case-pad form. They are
not converted.

## 5. Landing order

| Step | Kind | Depends on | Commit content |
|---|---|---|---|
| S1a | move | none | `rotation_vector`, `rotation_angle`; DF and `history_turn` call them |
| S2 | move | none | `ring_knots.h`; DF uses `SDF::KnotWalk`, `catmull_rom_strided` |
| S3 | move | none | `BumpPath`, `BumpCover`, `make_bump_path` (today's bound), `bump_field_at`, `cap_ring_arc`; DF uses them |
| S3b | behaviour | S3 | exact cover bound and antipodal fallback in `make_bump_path` |
| S6 | move | none | `ArenaFootprint`; DF budgets use it |
| S4a | move | none | `DistortedRingStack` and `CandidateTable` move verbatim from `core/render/scan/shapes.h` into the new `core/render/scan/ring_stack.h`; `scan.h` includes it |
| S4+S5 | move | S2, S4a | pool, shader, alternates and two-window scan added to `core/render/scan/ring_stack.h`; DF uses them; death cases move |
| S1b | behaviour | S4+S5 | `sweep_fraction` calls `math::half_turns` |

S4 and S5 land as one commit. `ShapeView::knot_range` reads the alternates, so
neither half compiles alone.

Each step that creates a file adds that file's row to the root `README.md`
tree fence in the same commit. The pre-commit fence rejects the commit
otherwise. New core files take the PolyForm NC header.

**Per-commit bar:**

- every step: the native build compiles;
- move steps: strict parity (§6.1);
- behaviour steps: their §6.3 property tests, plus the drift record and visual
  check (§6.1).

**Tip bar (once):** full suite, gates (doc coverage, format, docs-check, Teensy
size), `pio` phantasm and profile builds, WASM release and wasm-debug smoke,
then the device measurements (§6.2).

After pushing, confirm the daydream deploy.

## 6. Acceptance

### 6.1 Verifying moves and behaviour changes

**Isolated trees.** Steps change shared headers (`pose_table.h`,
`transformer.h`, `visibility.h`, `core/render/scan/shapes.h`, `arena.h`), so a header-shadowing
harness can compile both binaries against the same helper. Instead:

- **Two worktrees.** Base is the step's parent commit; candidate is the step
  commit.
- **One tree per binary.** The harness source is copied into each tree and
  compiled with that tree as its only include root, with no shared `-I` and no
  shared objects.
- **Strict flags.** Both binaries build with `-ffp-contract=off` and no
  fast-math.
- **Recorded in the harness log, for each binary:** `git rev-parse HEAD`, an
  empty `git status --porcelain`, the compiler path and `--version`, and the
  full compile command line.

**Sabotage check, against the candidate build.** In the candidate tree,
perturb the core code the step moved. Examples: add `1e-6f` in
`bump_field_at`; store `alternate + 1e-6f` in `KnotAlternates::retain`; flip
the sign branch in `rotation_vector`. Rebuild **only** the candidate and
confirm the hashes diverge, then revert. This proves the candidate binary runs
the candidate's core code. A step that touches several shared headers
sabotages each of them.

**Hash coverage** (each frame's pixels):

- 2,400 frames: one noise→balls→noise cycle at the default `cycle_frames()`.
- Full canvas and two quadrant clips.
- Default swept mode, `force_window_bakes`, and `force_column_reference`.
- A manual preset selection, which seeds preview balls.
- **Capacity flush.** Reuse the fixture from
  `test_displacement_field_alternate_groups_match_windows`: 288×144, 200
  warm-up frames, then `DisplacementFieldWhiteBox::arm_plane_pose` frames.
  Assert `stack_flushes > 1` in both binaries, so the extracted flush path is
  provably exercised. Ring count alone does not force a flush: retention also
  depends on resolution, pose, boundary coverage and whether alternate values
  differ.
- **Max-Speed run.** `ball_speed_min = ball_speed_max = 3.0` over one full
  ball phase.

**Move steps** pass when every hash matches.

**Behaviour steps (S3b, S1b)** are expected to change hashes. For each:

- **Drift record.** Report the first diverging frame, the count of changed
  pixels, and the maximum per-channel drift, for the review record.
- **Visual check.** In daydream, view the default cycle and the max-Speed ball
  phase, before and after. Bar: no clipped bump edges, no seams at the arm
  boundaries, and no new artifacts.

### 6.2 Device (tip only)

**Interleaved A/B.**

- **Images.** Base is the parent of the first step; candidate is the tip.
- **Captures.** Run 40 s pinned captures in A B A B A B order, three per
  image, for each of two phases: noise dwell and steady balls. That is 12
  captures. COM3 and COM4 are the same board and are directly comparable.
- **Recorded per capture, in ms:**
  - peak live render, with setup frame 1 excluded;
  - spilled frames in the steady window;
  - setup frame 1, reported separately and not judged.
- **Variation.** For each phase, the largest max − min peak among either
  image's own three captures.

**Bar, per phase:**

- the candidate's median peak ≤ the base's median peak + the variation; and
- candidate spills ≤ base spills.

Any cost beyond measured variation needs owner approval before landing. An
explanation does not pass the bar.

**The §4.5 hot-loop split** gets the same procedure, with base = fused loop
and candidate = split.

**Full-cycle capture.** Run `just profile DisplacementField` at the tip.
Bar: 0 spilled frames across the whole cycle, with setup reported separately.
Replace the shipping report.

**Size and ITCM.**

- The phantasm image passes its size gate.
- Report the RAM1 code delta in bytes against the base.
- No function outside DF leaves ITCM. A demotion of another function is a
  blocker, not a trade-off.
- The phantasm single TU sits at GCC's inline-growth cap, so compare the base
  and candidate symbol maps for newly outlined or newly inlined
  `Plot::rasterize` instances.

### 6.3 Tests

Every new test asserts a property, not a pinned output, and is wired into the
CI module list. The DF white-box suite keeps passing, with only the accessors
updated for moved counters.

| Location | Property |
|---|---|
| `tests/plot_scan/temporal_poses.h` | `rotate_by(previous-frame point, rotation_vector(a, b))` lands on the current-frame point within tolerance, including `r < 0` and near-identity steps; `history_turn` of a two-pose history equals `rotation_angle` |
| `tests/sdf/` | `KnotWalk` positions match the analytic `basis.v·cos θ + (u·cos a + w·sin a)·sin θ` within recurrence tolerance over `W` knots; `catmull_rom_strided` reproduces grid samples, reproduces linear data exactly, and is periodic |
| `tests/animation/` | `bump_field_at` at `t = 0` and `t = 1` equals `bump_field_profile` evaluated at the endpoint centres; `make_bump_path` rejects only paths with both reaches at or below the cut; **cover containment** (lands with S3b) over random endpoints with separations from 0 to just under π, near-antipodal pairs, and random reaches: for dense `t ∈ [0, 1]` and boundary points `p` at angle `reach(t)` from `c(t)` in random directions, `angle(p, cover.centre) ≤ cover.reach` within float tolerance. Endpoint field tests alone cannot catch an undersized cover. |
| `tests/plot_scan/` (with the chunk-cull tests) | brute force over random caps and rings: every knot with `dot(p, c) > cos_radius` lies inside `cap_ring_arc`; composed with the cover test, every knot that a swept bump's footprint touches lies inside the arc built from `cos_cover` |
| `tests/scan/ring_stack.h` (new) | pool lifecycle (`next` is stateless, `commit` counts, `release` destroys every slot); `knot_range` brackets both window values; `widen` leaves unmasked knots at the primary; `scan_two_windows` equals two independent `DistortedRingStack::draw` calls under the matching `TimeWindowScope`s with rows swapped by hand; `mark_arm_boundary` is conservative against brute force: any knot whose supported geometry reaches the widened arm plane is marked |
| `tests/death/scan_cases.h` | the moved pool death cases; `allocate` with `capacity < KNOTS`; `begin_frame` with `alternates->slots() < SLOTS` |

## 7. Decisions (owner, 2026-11-12)

1. **Q1, shader ownership.** `HueRowShader` moves to core with the pool. This
   reverses monolith rev 2's "shading stays in DF".
2. **Q2, pixel changes.** The bar is correctness, performance and beauty, not
   pixel identity. S3b uses the exact cover bound, and S1b lands.
3. **Q3, file placement.** Use the cleanest structure: the stack moves into
   `core/render/scan/ring_stack.h` alongside its pool and sweep (S4a).

## 8. Reassess later (not in this spec)

- Fold Plot's angle-addition recurrences (`core/render/plot/shapes.h` (baseline lines 1058),
  `:1292`) onto `SDF::KnotWalk`, if that is codegen-neutral on device.
- `select_hue_mode` as a `HueRotationRow` policy helper, once a second ring
  stack needs it.
- Convert persistent footprint sums to `ArenaFootprint`. That changes the
  static-assert margin, not behaviour.

## 9. Implementation and qualification record

The qualified code source is `182eb03a5c86145197a3723386ab0f0c66f12a3e`.
The final live-parent comparison uses baseline
`153be2802a81baead446ff23261e1e23eab50ad9`. The original extraction and
per-step parity evidence retain the historical `7a22dc38c` baseline.
Later test and documentation commits do not change the measured engine or
firmware code.

### 9.1 Native and render checks

- Native validation: 126 checks passed, followed by six affected/gate reruns
  and a passing scan-closure rerun. The added composed swept-bump/arc, reach-cut
  and row-closure checks also passed their focused module runs after the
  consolidated suite.
- Every pure extraction step passed strict comparison across 17 scenarios
  and 29,532 captured frames, with zero changed pixels. The corpus covers
  full and quadrant clips, all three sampling modes, manual preview balls,
  maximum Speed and the forced capacity flush. Both flush fixtures record two
  stack flushes. The inline-lookup variant, the value-owned visitor variant and the combined
  visitor/allocation-order variant also passed the same strict comparison.
- Both binaries use separate clean worktrees, their own engine objects and
  include roots, `-ffp-contract=off` and no fast-math. Per-binary compiler,
  command, source SHA and cleanliness records accompany the captures.
- Candidate-header sabotage changed the corresponding render hashes for the
  rotation, knot walk, swept field, cap arc, fused scan and alternate storage.
  The constexpr-only footprint sabotage failed the scratch-budget assertions.
  Production source was restored after every sabotage.
- Property mutations also failed as intended: omitted alternate widening
  caused 2 failures; omitted late installation caused 142; a 0.01 half-turn
  offset caused 1,514; zero arm-support reach caused 853. Omitting the final
  visitor `finish()` failed the closure-count assertion.
- Additional acceptance-property checks passed: composed swept-bump/arc
  containment (908,554 assertions), reach-cut rejection and output preservation
  (522), and shift-row closure (1,615). Their mutations failed as intended:
  narrowing the arc caused 231,673 failures; restoring the chord cover caused
  37,046; changing the reach-cut comparison caused 240; writing rejected-path
  outputs caused 27; omitting shift-row closure caused 2.

Release and debug WASM smoke passed for all 47 effects at both supported
resolutions. Whole-tree Doxygen coverage, the range commit checks and the
Markdown check with the live daydream checkout passed.

### 9.2 Behaviour drift

Counts below are changed pixel-frame samples, not distinct pixel locations.
Channel drift is measured in the framebuffer's 16-bit channel codes.

| Change and corpus | First divergent frame | Changed samples | Changed frames | Maximum channel drift |
|---|---:|---:|---:|---:|
| S3b, standard 17-scenario corpus | none | 0 | 0 | 0 |
| S3b, 288×144 default cycle | 909 | 70,078 | 1,000 | 2 |
| S3b, 288×144 maximum Speed | 2 | 253,145 | 1,489 | 7 |
| S1b, standard 17-scenario corpus | 22 | 1,796 | 518 | 3 |
| S1b, 288×144 default cycle | 14 | 5,423 | 1,560 | 4 |
| S1b, 288×144 maximum Speed | 10 | 5,494 | 1,644 | 92 |

At S3b's first production default divergence, 3 pixels change by at most 1
code; maximum Speed first changes 17 pixels by at most 1. Both S1b production
runs first change one pixel by one code. In the standard capacity-flush
scenario, S1b first diverges at frame 200: 590 pixels, maximum drift 3.
The containment properties pass; the standard S3b corpus's zero drift does
not replace the production-resolution checks.

### 9.3 Device qualification

Twelve valid 40-second captures compare the clean live baseline and final
candidate in A/B/A/B/A/B order for noise dwell and steady Balls. Setup frame 1
is excluded from live peaks and reported separately in the capture records.
Noise A1 used COM3; the other eleven pinned captures used COM4. This follows the cross-port comparison permitted by §6.2.

| Phase | Baseline median live peak (ms) | Candidate median live peak (ms) | Measured variation (ms) | Steady spills, all captures |
|---|---:|---:|---:|---:|
| Noise | 45.154 | 45.219 | 0.162 | 0 |
| Balls | 37.049 | 36.817 | 0.103 | 0 |

Noise increases by 0.065 ms, within measured variation. Balls decreases by
0.232 ms. Both phases pass the peak and spill bars.

| Shipping phantasm memory | Baseline (bytes) | Candidate (bytes) | Delta (bytes) |
|---|---:|---:|---:|
| RAM1 code | 186,984 | 187,160 | +176 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH code | 1,171,760 | 1,172,080 | +320 |
| FLASH data | 791,196 | 791,256 | +60 |

The phantasm size and hot-layout gates pass. No function outside DF leaves
ITCM or disappears from its ITCM symbol set. No exact outside-DF symbol moves
to flash. `Plot::rasterize` symbol groups have no additions, removals, byte-size
changes or placement changes.

**Full-cycle capture: PASS.** The final source-attested run covers 188
complete windows and a Noise→Balls→Noise cycle, with zero spills over 3,015
live frames. Its peak live render is 45.164 ms at frame 603; the Balls peak is
37.397 ms. Setup frame 1 takes 18.393 ms and is reported separately from the
live peak. Capture validity is **VALID**, with the recorded 0.384 ppm timing
check.

The maintained shipping report is
[the 2026-11-12 DisplacementField timing report](../profiles/shipping/profile_displacementfield_teensy_2026-11-12.md).

### 9.4 Daydream visual check

Actual Daydream UI captures at 288×144 cover representative default Noise,
manually selected Balls and maximum-Speed Balls before and after the change.
Both Speed controls are verified at 3.000 for the maximum-Speed images. The
installed candidate WASM SHA-256 is
`30f278eb53231beb5087d1d1461529e19f28f2393ee5028b237fe9a1321fc904`.
Inspection finds no new clipped bump edges, arm-boundary seams or visible
artifacts. Both browser runs record zero console, uncaught or request errors.
These are representative phase snapshots. The native corpus and the final
device capture supply the full-cycle coverage.
