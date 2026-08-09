# ShaderBall — unifying Liquid2D and Flyby

**Status: LANDED.** `effects/ShaderBall.h` ships, `Liquid2D` and `Flyby` are
gone from every roster, and the effect carries 15 presets. §5's lens blend is
superseded by the shipped one; where this spec and the code disagree, the code
is the source of truth. §11 is the executed roster checklist, kept as the record
of every oracle a roster change touches.

Merge `Liquid2D` and `Flyby` into a single stereographic shader effect,
`ShaderBall`, whose presets span both effects' looks and the mixed looks
neither can reach today. The two effects already share their entire skeleton
(stereographic projection → noise warp → trig pattern → pole normalize →
palette); they diverge only along five axes, and every one of those axes can
be expressed as a **continuous float parameter**. That is the core design
decision: no mode enums, no per-preset function pointers — every preset is a
point in one parameter space, so the existing preset-lerp machinery morphs
between *any* two presets, including cross-family ones, with no discrete pop.

## 1. Why merge

- The classes are ~85% structurally identical (both note this in their own
  doc comments and share `stereo_noise_warp` / `stereo_pattern_args` /
  `pole_normalize_pattern` in core). Two copies of the skeleton cost flash,
  ITCM, test surface, and roster slots.
- The interesting space is *between* them: a glitch-lensed grid fly-through,
  a hue-warped cross-coupled liquid, a half-lensed slow orbit. Today those
  are unreachable; in ShaderBall they are just presets.
- One effect replaces two in every roster (HS_EFFECT_LIST, Phantasm playlist,
  GOLDEN_ROSTER, README gallery, daydream favorites), shrinking the
  maintenance checklist permanently.

## 2. Current anatomy

Shared skeleton (identical in both):

- `Effect(W, H, {.strobe = true})`, no pixel persistence; full-canvas
  `Scan::Shader::draw<W, H, 1>`.
- FastNoiseLite OpenSimplex2 warp source (but the generators' frequency,
  seeding, and sharing diverge — see the table below).
- Per pixel: project to stereographic plane → `stereo_noise_warp` →
  `stereo_pattern_args` + fast trig pattern → `pole_normalize_pattern` →
  palette lookup.
- Noise-time accumulator wrapped to `STEREO_NOISE_TIME_PERIOD`; trig phases
  wrapped to 2π (white-box tested in both).
- `Presets<Params, N>` bank + timeline-scheduled `Animation::Lerp`.
- Baked generative palette in the persistent arena.

Divergences:

| Axis | Liquid2D | Flyby |
|---|---|---|
| Camera | two `RandomWalk`s (inner + global), no spin | fixed base orientation + looping Y rotation (2π / 300 frames) |
| Lens | glitch lens (hemisphere fold + degree-3 rational map) between the walks | none |
| Warp noise generator | **shared with the inner walk**: `RandomWalk`'s ctor reseeds it (random seed) and sets its frequency; init then re-sets `SetFrequency(0.02)` | dedicated, untouched: FastNoiseLite defaults — frequency **0.01**, fixed default seed |
| Warp time | `t * 0.5` | `t * 0.3` |
| Pattern | cross-coupled: `sin(x + c·sin(y + φ₁)) · cos(y + c·cos(x − φ₂))` | separable grid: `sin(x + φ₁) · cos(y − φ₂)` |
| Phase rates | φ₁ += dt, φ₂ += 0.8·dt | φ₁ += speed, φ₂ += speed·drift (drift = live slider, default 0.7) |
| Color | `StaticPalette` + `BreatheModifier` (amp 0.15, cycle-phase driven) | raw palette + `alpha *= (1 − value)` + `hue_rotate(−displacement · hue_shift)` via gamut LUT |
| Palette recipe | STRAIGHT / COMPLEMENTARY / CUP / VIBRANT / 75 | STRAIGHT / SPLIT_COMPLEMENTARY / FLAT / MID / 42 |
| Preset choreography | repeating 90–150-frame timer; the 60-frame **staggered** lerp runs concurrently inside the period (actual hold between blends: 30–90 frames) | perpetual 480-frame **parallel** lerp loop, `.then()`-chained with zero gap |
| Persistent arena | palette only | palette + gamut LUT (131,072 B) |

## 3. The unifying observation: one pattern formula

Generalize the pattern to

```
pattern(p, φ₁, φ₂) = sin(p.re + c·sin(p.im + φ₁) + s·φ₁)
                   · cos(p.im + c·cos(p.re − φ₂) − s·φ₂)
```

with `c = complexity` (cross-coupling depth) and `s = phase_direct`
(direct phase feed). This is an **exact** superset:

- `c = complexity, s = 0` → Liquid2D's `sample()` verbatim.
- `c = 0, s = 1` → Flyby's `sample()` verbatim.

Both parameters are continuous and lerpable, so a preset transition morphs
the *function itself* smoothly. `stereo_pattern_args` still bounds `p`
first, exactly as today.

Phase advance unifies the same way: `φ₁ += speed`, `φ₂ += speed ·
phase2_rate` per frame (Liquid2D: `phase2_rate = 0.8`; Flyby: `phase2_rate
= drift = 0.7`). Warp time becomes `noise_time · warp_time_scale`
(0.5 vs 0.3).

## 4. Pipeline

Per frame (all frame-constant work hoisted out of the shader):

```
timeline.step()                            // walks, drivers, preset lerp
noise_time = fmodf(noise_time + speed, STEREO_NOISE_TIME_PERIOD)
sin_phase  = fmodf(sin_phase  + speed, TWO_PI_F)
phase2     = fmodf(phase2 + speed * phase2_rate, TWO_PI_F)
spin_phase = fmodf(spin_phase + spin_rate, TWO_PI_F)
cycle_phase = fmodf(cycle_phase, TWO_PI_F)  // Driver-fed, breathe consumer
inner = R_y(spin_phase) ∘ BASE ∘ slerp(I, inner_walk, wander)   // per frame
outer = slerp(I, outer_walk, wander)                            // per frame
```

Per pixel:

```
rv = outer.unorient(v)
sv = lens_blend(rv, lens_mix)          // §5
z  = stereo(inner.unorient(sv))
r_sq = |z|²
(w, disp) = stereo_noise_warp(z, r_sq, noise, warp_scale, warp_strength,
                              pole_fade, noise_time * warp_time_scale)
pattern = generalized formula (§3) on stereo_pattern_args(w, pattern_freq)
value = pole_normalize_pattern(pattern, r_sq, pole_fade)
value = min(value, 1.0f − ULP)         // §6: Wrap palette folds exact 1.0 → 0
c = palette_get(value)                 // §6: bank blend + breathe
c.alpha *= (1 − value * value_fade)
if (hue_shift ≠ 0 or disp-dependent) c = hue_rotate(c, −disp * hue_shift)
return c
```

Every per-pixel branch (`lens_mix == 0`, `complexity == 0`, `hue_shift ==
0`, `palette_pos` integral) tests a **frame-constant** value, so the M7
branch predictor makes the skips nearly free without duplicating the inner
loop into ITCM-costly specializations (see §9).

## 5. Camera and lens

Three continuous controls replace the discrete camera modality:

- **`spin_rate`** — Y-axis rotation rate in rad/frame. Flyby: 2π/300 ≈
  0.0209; Liquid2D: 0. The rotation is an explicit wrapped `spin_phase`
  accumulator (not a fixed-duration `Animation::Rotation`), so the rate is
  lerpable. `BASE` is Flyby's `make_rotation(Vector(0,0,-1),
  Vector(0,-1,0))`; for wander-dominated presets it is an invisible fixed
  pre-rotation of a random walk.
- **`wander`** — random-walk mix in [0, 1]. Both `RandomWalk`s (inner +
  outer) run permanently on shadow orientations; the effective contribution
  is `slerp(identity, walk_q, wander)`, computed **once per frame**. At 0
  the walks contribute nothing (Flyby); at 1 they contribute fully
  (Liquid2D). Walks keep running while mixed out so a transition into
  wander picks up a live, well-mixed trajectory rather than a cold start.
  Two subtleties are load-bearing:
  - **Hemisphere continuity.** The walk quaternion accumulates unboundedly;
    when its total rotation angle crosses π, `dot(identity, walk_q)` flips
    sign and a naive `slerp(I, walk_q, wander)` jumps between the `w·α` and
    `w·(2π − α)` partial rotations — a visible camera pop for any wander ∈
    (0, 1). Keep a sign-adjusted copy of the walk quaternion (flip its
    stored sign whenever `dot(prev, curr) < 0`) and slerp against that.
  - **Inner mix sense.** Liquid2D applies the inner rotation *forward*
    (`orientation.orient`), while the unified pipeline uses `unorient`
    throughout. An inverted random walk is still a random walk, so this is
    accepted as one of the Liquid2D distribution-level deviations (§10)
    rather than special-cased.

  **Noise generators: three, strictly separated.** `RandomWalk`'s ctor
  reseeds (randomly) and sets the frequency of whatever generator it is
  handed, so each walk owns a private generator and the **warp generator is
  never given to a walk** — Liquid2D's inner walk sharing the warp source
  was an accident that made its warp field randomly seeded and
  frequency-coupled to walk construction order. ShaderBall's warp generator
  gets an explicit `SetFrequency(0.01)` (Flyby's effective default) and the
  default seed, making the warp field deterministic; Liquid2D-mapped preset
  `warp_scale` values are doubled to compensate for the 0.02 → 0.01 base
  frequency (`(z·k)·0.02 ≡ (z·2k)·0.01` bit-for-bit — a ×2 is float-exact).
- **`lens_mix`** — glitch-lens blend in [0, 1], with the shortcuts
  `lens_mix == 0 → rv`, `lens_mix == 1 → glitch(rv)`. The direction-space
  blend `sv = normalize(rv + lens_mix · (glitch(rv) − rv))` specified here is
  **superseded**: `glitch(v) = −v` has exactly two solutions, `(0, 0, ±1)`,
  where the blended direction flips by 180° at `lens_mix = 0.5`, and the
  field is badly distorted across roughly `lens_mix ∈ [0.4, 0.6]`. A
  squared-magnitude guard removes the NaN but not the flip, so the guard is
  not a fix. The shipped blend works in **sample space**: sample the pattern
  once on `rv` and once on `glitch(rv)` and lerp the two `PatternSample`s
  (`blend_lens_samples` in `effects/ShaderBall.h`, which is the source of
  truth for the blend). It is linear in `lens_mix` through the former
  midpoint and costs a second pattern evaluation across the open interval,
  which is the cost §9 budgets.

`apply_glitch_lens` is effect-local; nothing else uses it. It ships as a
closed-form rational map with no hemisphere fold and no pole guard — its
denominator `1 + 3y²` is at least 1, so the map is finite everywhere on the
sphere and needs neither.

## 6. Color

- **Palette bank.** Bake both existing recipes at init:
  slot 0 = Liquid2D's (STRAIGHT/COMPLEMENTARY/CUP/VIBRANT/75), slot 1 =
  Flyby's (STRAIGHT/SPLIT_COMPLEMENTARY/FLAT/MID/42). A float
  **`palette_pos`** ∈ [0, N−1] selects: integral → single LUT lookup;
  fractional → two lookups + `Color4` lerp. Lerping `palette_pos` between
  presets *is* the palette crossfade — no extra machinery, and the double
  lookup only costs during transitions. The bank is trivially extensible for
  future presets.
- **`breathe_depth`** ∈ [0, 0.3] — `BreatheModifier` amplitude becomes
  preset-driven (0.15 for Liquid2D looks, 0 for Flyby looks). At 0 the
  modifier's coordinate shift is exactly 0, so Flyby lookups are unchanged.
  Driven by `cycle_phase` ← `Driver(cycle_speed)` as in Liquid2D. The
  consuming palette keeps `Wrap = true` (BreatheModifier's `requires_wrap`).
  **Wrap-fold hazard:** `pole_normalize_pattern` can produce exactly 1.0,
  and a Wrap palette folds 1.0 → 0.0 — the palette's *other end* — where
  Flyby's raw `BakedPalette::get` clamps. Clamp `value` to `1.0f − ULP`
  before every lookup (§4) so pattern peaks can't flip color.
- **`hue_shift`** ∈ [0, 1] — Flyby's displacement-driven
  `hue_rotate(−disp · hue_shift)`; skipped when 0 (frame-constant test).
  Note `hue_rotate(c, 0)` is **not** an identity — it round-trips pixel →
  linear RGB → OKLab → gamut map → pixel — and Flyby applies it
  unconditionally, so the skip is a small deliberate visual delta on
  zero-hue-shift presets (F3), accepted as part of the §10 parity softening.
  The gamut LUT is always bought (§8) because any transition may pass
  through nonzero hue_shift.
- **`value_fade`** ∈ [0, 1] — generalizes Flyby's `alpha *= (1 − value)` to
  `alpha *= (1 − value · value_fade)`; 0 reproduces Liquid2D's untouched
  alpha.

## 7. Params, presets, choreography

### Params (all float, all preset-lerped, all registered animated)

| # | Param | Range | Liquid2D value | Flyby value |
|---|---|---|---|---|
| 1 | `warp_scale` | 0.1 – 100 | per preset | per preset |
| 2 | `warp_strength` | 0 – 30 | per preset | per preset |
| 3 | `warp_time_scale` | 0.05 – 1 | 0.5 | 0.3 |
| 4 | `pattern_freq` | 1 – 20 | per preset | per preset |
| 5 | `speed` | 0 – 5 | time_speed | speed |
| 6 | `complexity` | 0 – 3 | 0.5 / 3.0 | 0 |
| 7 | `phase_direct` | 0 – 1 | 0 | 1 |
| 8 | `phase2_rate` | 0 – 2 | 0.8 | 0.7 (was Drift) |
| 9 | `pole_fade` | 1 – 20 | per preset | per preset |
| 10 | `spin_rate` | 0 – 0.05 | 0 | 2π/300 |
| 11 | `wander` | 0 – 1 | 1 | 0 |
| 12 | `lens_mix` | 0 – 1 | 1 | 0 |
| 13 | `palette_pos` | 0 – 1 (N−1) | 0 | 1 |
| 14 | `breathe_depth` | 0 – 0.3 | 0.15 | 0 |
| 15 | `cycle_speed` | 0 – 1 | 0.05 | 0 |
| 16 | `hue_shift` | 0 – 1 | 0 | per preset |
| 17 | `value_fade` | 0 – 1 | 0 | 1 |

Conventions carried over: `params = presets.get()` before any
`register_animated_param` (register captures the default);
`preset_in_ranges` + `all_presets_in_ranges` static_assert; `sizeof(Params)
== 17 * sizeof(float)` static_assert tied to the lerp implementation;
public-first member layout.

Flyby's "Drift" stops being a live-only slider and becomes the
preset-driven `phase2_rate`. That trades a live control for lerpability;
"Pause Animation" + slider takeover still gives live control, matching how
every other preset param already works.

### Preset bank

Seven mapped presets: L0–L1 and F0–F4 carried into the superset via the table
above (unlisted fields take the other family's neutral value), with one mapping
adjustment: Liquid2D `warp_scale` values are **doubled** (1.5 → 3.0) to
compensate for the unified warp generator's 0.01 base frequency (§5).

**Gate-endpoint rule:** the skip branches in §4 stay latched only because
`Animation::Lerp`'s final step lands bit-exactly on the target, which holds
when the gate params' endpoints are exactly 0 or 1 (where `a + (b−a)·t` is
exact) and the easing satisfies `e(1) == 1.0f` (true for `ease_in_out_sin`).
Presets must therefore keep `complexity`, `phase_direct`, `lens_mix`,
`hue_shift`, `value_fade`, and `wander` at exactly 0 or 1 whenever they
intend the corresponding fast path — a near-0 value silently un-latches the
skip forever after the first transition. Pinned by test (§10).

Follow-up presets
should deliberately mix families — e.g. a lensed grid (`lens_mix 0.6,
complexity 0, spin on`) or a hue-warped liquid (`complexity 2, hue_shift
0.4, wander 1`) — that is the payoff of the merge, curated visually after
landing.

### Choreography

A `Choreo` array parallel to `PRESETS` (static_assert equal sizes), one
entry per preset, consumed when *entering* that preset:

```
struct Choreo { uint16_t dwell_min, dwell_max, blend_frames; bool staggered; };
```

Dwell is **sequential**: hold, then blend. The state machine on entering a
preset is:

- `dwell_max == 0`: chain the next `Animation::Lerp` **directly** from the
  previous blend's `.then()` — no timer, zero idle frames (Flyby's exact
  structure). A repeating `RandomTimer(0, 0, …)` is **forbidden** here: it
  fires every step and would spawn a fresh Lerp per frame.
- `dwell_max > 0`: arm a **one-shot** `RandomTimer(dwell_min, dwell_max)`
  whose callback schedules the Lerp; the Lerp's `.then()` re-enters the
  state machine for the next preset.

Entries:

- Flyby-style: `dwell 0/0`, `blend 480`, parallel → perpetual motion.
- Liquid2D-style: `dwell 30/90`, `blend 60`, staggered. (Liquid2D's
  repeating 90–150-frame timer ran the 60-frame blend *inside* its period,
  so the true hold between blends is 30–90 frames; sequential
  dwell-then-blend with 30/90 reproduces the original 90–150-frame cadence.)

Liquid2D's staggered `Params::lerp` (each changed field gets its own time
slice) generalizes unchanged to 17 fields; the `staggered` flag selects
between it and the parallel lerp. Cross-family transitions default to
parallel — a staggered walk through 15 changing fields would spend most of
the blend in unbalanced intermediate states. Scheduling stays on the
existing `timeline.add_pausable` + `Animation::Lerp(...).then(...)` /
`RandomTimer` primitives; both patterns already exist in the two effects.

## 8. State, wrap invariants, memory

Accumulators (all wrapped every frame, exactly the union of today's):
`noise_time` → `STEREO_NOISE_TIME_PERIOD`; `sin_phase`, `phase2`,
`spin_phase`, `cycle_phase` → 2π. The two walk orientations and the
timeline keep Liquid2D's declaration-order rule (orientations before
`timeline`, which clears RandomWalks on teardown).

Persistent-arena footprint:

```
gamut LUT: gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS)  // 131,072 B
palettes:  2 × BakedPalette::required_arena_bytes()                   // 2 × 2,052 B
```

`static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET)` as in Flyby
today; the delta vs Flyby is one extra palette bake (+2,052 B), well inside
the partition Flyby already fits.

No heap, arena only; no hardware asserts on the per-pixel path (the nlerp
guard is a branch, not a check).

## 9. Performance budget

Shipping device profiles (Teensy, 96×20 segmented, 62.5 ms window,
10,368 px quadrant):

- Liquid2D: `lq_shader_draw` 28.3 ms (1,639 cyc/px), ~34 ms headroom.
- Flyby: `fly_shader_draw` 37.3–44.0 ms across presets, worst preset 44.0
  ms (2,546 cyc/px). The report's headline says 17 ms of margin, but its
  own ISR accounting puts the render budget at ≈59.3 ms/window with the
  peak frame leaving **14.2 ms unspent — that is the binding margin**.

Cost deltas for a Flyby-look preset running in ShaderBall (all
frame-constant-branch skips active): one extra orientation transform
(`outer.unorient`: quaternion rotate is 15 mul + 15 add, no division), the
predicted-not-taken lens/complexity/hue branches, and the
`value_fade`/`palette_pos` arithmetic. Estimated well under 100 cyc/px
(< 1.7 ms) against the 14.2 ms margin.

Transient worst case is mid-transition between families: lens nlerp +
cross-coupled pattern (2 extra fast trig) + hue rotate + dual palette
lookup simultaneously. Rough estimate +300–500 cyc/px ≈ +5.2–8.6 ms over
the Flyby worst preset — an order-of-magnitude bound, not a budget, and it
consumes up to ~60% of the honest margin; note the worst Flyby preset is
itself a mid-lerp state, so the stack-up is real. The **first device
profile must capture a cross-family transition window**, not just preset
holds (use the teensy-profile flow; sweep, don't sample).

Mitigations if the transition spills (in order): shorten cross-family
`blend_frames`; schedule cross-family transitions so lens_mix and
complexity ramps don't peak together (stagger just those two fields); only
then consider a second shader variant (ITCM-costly — the ShapeShifter
campaign priced a duplicated inner path at ~7 KB, and the budget has no
room for that by default).

Keep the shader inside an `HS_O3` region exactly as both effects do today
(`stereo_noise_warp` is already `HS_O3_FN`). Re-measure the phantasm size
gate after landing — deleting two effect instantiations and adding one
superset should net out negative, but measure, don't assume.

## 10. Testing and validation

- **Port the white-box wrap test**: seed every accumulator past its ceiling,
  assert `noise_time ∈ [0, STEREO_NOISE_TIME_PERIOD)` and all four trig
  phases ∈ [0, 2π) from frame 0 (union of `test_liquid2d_phase_wrapped` and
  `test_flyby_phase_wrapped`).
- **Port `test_liquid2d_glitch_lens_unit_norm`** unchanged (lens moves
  verbatim), plus a case pinning the nlerp guard fallback.
- **New: formula-reduction test.** Pin `sample(w, φ₁, φ₂)` at
  `(c, s) = (complexity, 0)` against Liquid2D's closed form and `(0, 1)`
  against Flyby's, over a grid of arguments — this is the load-bearing
  equivalence claim of the whole merge. Extend it with the gate-latch
  assertion (§7): simulate a full `Animation::Lerp` between presets and
  assert every gate param lands **bit-exactly** on its 0/1 endpoint, so the
  skip branches provably re-latch after transitions.
- **Visual parity via the framebuffer dump harness** (render → PNG, don't
  theorize). Neither family is bit-exact; do not gate the migration on a
  bit-exact framebuffer diff — it will fail.
  - *Flyby-mapped presets*: near-identical, ULP-level drift. Known exact
    mechanisms: the spin trajectory (Flyby's `Animation::Rotation`
    accumulates incremental quaternion products; ShaderBall rebuilds
    `R_y(spin_phase) ∘ BASE` from a phase accumulator — equal in math,
    different in floats from frame ~2 on); the `hue_rotate` skip at
    `hue_shift == 0` (§6, affects F3); and preset-cycle timing. Compare
    visually / statistically with tight tolerances.
  - *Liquid2D-mapped presets*: distribution-identical only. Deviations:
    walk composition order, the BASE pre-rotation, the inner mix applying
    the walk inverted (`unorient` vs Liquid2D's `orient`, §5), and the warp
    field being deterministic instead of walk-ctor-reseeded (§5).
- Roster oracles update mechanically: smoke suite iterates HS_EFFECT_LIST;
  `GOLDEN_ROSTER` in `tests/test_param_marshal.h` is a hand-maintained
  parallel copy — remove Flyby and Liquid2D rows, insert ShaderBall between
  RingSpin and ShapeShifter.
- Wire new tests into the CI module list and update the roster row in
  `run_tests.cpp` for the deleted Flyby/Liquid2D cases and the added
  ShaderBall ones — the row pins **both** an exact case count (enforced by
  `tests/check_case_calls.cmake`) and the assertion floor, so both columns
  change.

## 11. Roster removal/addition checklist (both repos)

Holosphere:

- `effects/ShaderBall.h` new (`git mv effects/Flyby.h` to keep history on
  the closer parent); Liquid2D's header deleted.
- `core/engine/effects.h`: both `#include`s → one; both `X()` rows → one, at
  the case-insensitive-alphabetical slot (after RingSpin, before
  ShapeShifter) in **both** HS_EFFECT_LIST and HS_PHANTASM_EFFECT_LIST;
  backslashes at column 80. `HS_EFFECT_COUNT` 24→23; the `COUNT − 2`
  phantasm static_assert holds unchanged (both lists shrink by one).
- `tests/test_effects.h`: white-box blocks + runner calls (§10);
  `run_tests.cpp` roster row — case count and assertion floor (§10).
- `tests/test_param_marshal.h`: GOLDEN_ROSTER (§10).
- `README.md`: two gallery cards → one (`?effect=ShaderBall`,
  `docs/screenshots/ShaderBall.png`, alt, heading) **plus** the prose
  references: the pipeline walkthrough naming both effects, the
  `stereo_noise_warp` section, the `hue_rotate` doc, and the gamut-LUT
  passage ("`Flyby` and `MeshFeedback` arm a copy").
- `Holosphere.ino` commented show list.
- `tools/profile_one.sh`: `CYCLERS` list names both effects.
- `tools/profile_sweep.sh`: per-effect sweep rows for both; Flyby's carries
  `-D HS_PROFILE_EPOCH_REVS=2560` — ShaderBall's mixed choreography needs
  its own epoch decision (the preset cycle is longer and partly random).
- `docs/profiles/` READMEs: retire the Liquid2D/Flyby rows and reports; new
  ShaderBall report after the landing profile. Also
  `docs/profiles/memory/arena_high_water.md` (rows for both) and
  `docs/strobe_columns_audit.md`.
- wasm bootstrap needs no edit (it takes `WASM_EFFECT_NAMES[0]`, which stays
  BZReactionDiffusion).

daydream (separate repo, land via worktree + FF):

- `daydream.js`: Flyby appears only in `HiResFavorites`, Liquid2D only in
  `LoResFavorites` — remove each from its list and add ShaderBall to
  **both** (it spans both looks; the arrays are curated, not alphabetical,
  so the position is a curatorial choice). These hardcoded arrays feed both
  the sim's effect list and the `?effect=` validator; missing this breaks
  `just screenshots`.
- daydream's `README.md` is the Holosphere README mirror — the same gallery
  cards and prose references above must land there too (mirror-accumulator
  flow).
- `daydream/tests/engine_contract_wasm.test.js` needs no edit (verified: it
  names no effect).
- Screenshot order: land rename → install wasm → fix favorites →
  `just screenshots` from the main tree (Playwright lives only there).

## 12. Migration plan

1. **Land ShaderBall alongside both old effects** (roster temporarily 25;
   phantasm assert edited to match for the interim, or do it in one atomic
   land — one atomic land is simpler given the assert). Preferred: single
   commit series in one worktree that adds ShaderBall and removes both, so
   every roster oracle transitions once.
2. Parity validation (framebuffer harness) before the removal commit.
3. Device profile (holds **and** a cross-family transition window) +
   phantasm size gate re-measure.
4. daydream favorites + screenshots.
5. Curate 2–3 new mixed-family presets; retire any that read as redundant.

## 13. Curation decisions

- Preset count: 12 — five wandering-liquid, five spinning-grid, and two
  mixed (grid look on the liquid palette).
- `phase2_rate` ships preset-driven as the animated **Drift** param over
  [0, 2]. No extra live-only slider was added; pause plus slider takeover is
  the live control.
- `spin_rate`'s ceiling ships at the proposed 0.05 rad/frame (≈ 3× Flyby).
