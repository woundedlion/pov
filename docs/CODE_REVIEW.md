# Holosphere / daydream — Code Quality Review

**Scope:** the C++ rendering engine and firmware (`core/`, `effects/`, `hardware/`,
`targets/`, `tests/`, `scripts/`) and the daydream web simulator
(`*.js`, `tools/`, JS tests). Out of scope per request: `effects_legacy.h`,
`targets/Holosphere/Holosphere.ino`, `rotate.h`.

**Method:** 19 reviewers audited every in-scope file line by line against the
README architecture; every candidate defect was then re-examined by an
independent validator that re-read the cited code. Findings that the second pass
could not stand behind were dropped (see *Validated Non-Issues*). Only real,
actionable defects with a concrete minimal fix are listed — consistent with the
eligibility bar of the `code-review-fix` workflow.

---

## Overall Grade: **A−**

This is an exceptionally well-engineered codebase. The dominant signal across
every subsystem is that the obvious failure modes — divide-by-zero, NaN from
out-of-domain `acos`/`sqrt`, capacity overflow, pole/seam singularities, arena
OOM, ISR/main-loop races, detached WASM views — are already guarded *and the
guard's rationale is documented at the site*. The defects below are real but
overwhelmingly latent, cosmetic, or confined to dormant/extreme-edge paths; the
codebase carries no Critical defect and only two High-severity items, one of
which is dormant.

### Quality Dimensions

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A− | Hot paths verified sound across the board; defects are latent (unit-precondition traps, frozen baselines) or confined to dormant builds. One dormant High (DMA single-board section attr) and a handful of Mediums. |
| **Memory safety / UB-freedom** | A | Placement-new + `std::launder`, `nullptr+0` guards, nothrow `static_assert`s, alignment realign-on-restore, fixed-capacity pools with trap-on-overflow. No live UB found. |
| **Concurrency & ISR-safety** | A− | Hardware single-writer/single-observer model with correct atomics/barriers is textbook (A). Pulled to A− by the simulator's missing render-phase worker watchdog, which can wedge the worker pool. |
| **Performance** | A | Compile-time `<W,H>` specialization, split trig LUTs, arena bump-allocation, baked 256-entry palettes, zero-alloc type-erased callables, deliberate per-frame-cost discipline. Trade-offs are chosen and documented, not accidental. |
| **Numerical robustness** | A− | Domain clamps on every transcendental, sin²φ/speed floors, antipodal/pole stabilization, round-to-nearest fixed-point. Residual: slow non-orthonormal chroma drift, near-gray hue torsion, float-time ceilings on multi-hour runs. |
| **API design / interface expressiveness** | A | Owned/borrowed split (`ArenaVector`/`ArenaSpan`, `Fn`/`FunctionRef`), explicit coordinate-space ctors, `[[nodiscard]]`, compile-time-checked filter ordering, designated-initializer draw params. Misuse-resistant. |
| **Architectural elegance** | A | Clean shape/rasterizer split, variadic compile-time filter pipeline with automatic domain transitions, single-block partitioned arena, X-macro single-sources-of-truth, two-repo product seam. Coherent and self-consistent. |
| **Maintainability** | A− | Comments are dense but load-bearing. Docked for a few comments that overstate invariants, a couple of near-identical effects diverging on guard style, and dual-source-of-truth loop bounds. |
| **Portability / abstraction** | A | Arduino/WASM/Desktop split is disciplined; FastLED integer mocks are bit-exact against the `*_C` references; documented sim/device divergences with wrap-analysis. |
| **Test coverage & quality** | A− | Behavioral oracles (Euler characteristic, energy conservation, exact wire bytes, RNG draw-count pins), a verified fail-fast death harness, host-testable hardware logic, cross-language parity pins. Docked for a zero-assertion blind spot and NDEBUG-gated coverage that can silently vanish. |
| **Documentation** | A | The README is a genuinely outstanding systems document; inline rationale captures the *why* at the point it matters. |
| **Color-science rigor** | A | 16-bit linear-light blending, OKLCH shortest-arc hue, Ottosson chroma-preserving gamut projection, lightness-coupled chroma envelope. |
| **Cross-language contract integrity** | A | Resolutions, roster, and param streams enumerated from one source; embind shell kept thin over host-testable marshal/predicate cores; WASM parity tests drive the real module. |
| **Error handling / fail-fast discipline** | A | Consistent `HS_CHECK` cold-path trap philosophy, verified by a death harness. A couple of sites use `assert`/silent-return where the philosophy prefers an always-on trap. |
| **Build / CI robustness** | A− | Smoke-drives every effect at every resolution, gates arena/stack high-water, dual roster pins, both Teensy images size-gated. Minor: case-sensitivity and a couple of unpinned standalone test executables. |

---

## Prioritized Fix List

Each item is independently validated and carries a concrete fix. Numbering is
sequential across the whole list.

### High Priority

1. ✅ **`hardware/pov_single.h:215-223` — DMA single-board build has no `ledController_` definition (dormant).** For a `USE_DMA_LEDS` `POVDisplay` build there is only an explanatory comment, no out-of-line definition. `pov_segmented` lands its controller in OCRAM via an explicit specialization at `targets/Phantasm/Phantasm.ino`; the single-board path has none, so a naive `DMAMEM` definition would have its section attribute silently dropped by GCC (vague-linkage template static), placing the eDMA TX buffer in DTCM where the cache-flush is a no-op and the DMA engine reads stale memory. Dormant only because the shipping `Holosphere.ino` does not define `USE_DMA_LEDS`. *Fix:* require any DMA single-board target to provide an explicit `DMAMEM ...::ledController_` specialization, and add a `static_assert`/comment in `pov_single.h` directing it, mirroring the segmented contract.

2. ✅ **`daydream/segment_controller.js:771-782` — no per-render worker watchdog.** `init` (20 s) and `boot` (4 s) have deadlines, but the steady-state render phase — the only phase that runs every frame — has none. `pending` decrements solely on `frame` messages, so a worker that accepts `render` but never replies (a non-trapping JS/WASM hang that fires no `onerror`) leaves `pending` stuck and `renderInFlight` permanently `true`: the pipeline freezes with no fault overlay or recovery. *Fix:* arm an `unref()`'d render watchdog in `renderParallel()` (a few × `SLOW_FRAME_MS`) that calls `onWorkerFault(-1, 'render timed out')` if `pending > 0` when it fires; clear it in the `frameResolve` callback, reusing the existing latch/overlay/recovery path.

### Medium Priority

3. ✅ **`core/filter.h:898-899` — `AntiAlias::plot` computes its second AA tap from an unwrapped column.** `x1 = fast_wrap((int)x_floor + 1, W)`; for a 2D producer honoring the pipeline's documented `[-W, 2W)` x-contract, `x_floor == 2W-1` makes `x_floor+1 == 2W`, violating `fast_wrap`'s `x < 2W` precondition (out-of-bounds by one column in release, assert trip in debug). Every sibling that faces the same contract (`Feedback`, `ChromaticShift`) wraps first, then derives the neighbor. *Fix:* `int x0 = fast_wrap((int)x_floor, W); int x1 = fast_wrap(x0 + 1, W);`.

4. ❌ **`core/transformers.h:306-311` + `core/animation.h` — gnomonic Möbius warp folds at the equator.** *Rejected:* the input-sign hemisphere is load-bearing (it makes an identity Möbius an identity transform), and the gnomonic plane carries no hemisphere to derive from (antipodes coincide), so the proposed fix is ill-defined; the residual equator-band tearing is inherent to the projection — an artistic property of `GnomonicStars`, not a correctness defect. `gnomonic_mobius_transform` reconstructs the output hemisphere from the *input* vertex's `sign(v.y)`. Paired with `MobiusWarpEvolving` (which drives all coefficients every frame and can push points across `y=0`), near-equator vertices are projected back onto their origin hemisphere, producing a persistent crease at `y≈0`. *Fix:* derive the hemisphere from the transformed coordinate (or detect the equator crossing) rather than the input sign; minimally, bound `MobiusWarpEvolving`'s amplitude away from the equator and document the constraint.

5. ✅ **`core/animation.h:2104-2107` — `MobiusWarpEvolving` ignores live baseline edits.** `Transformer::prepare_frame()` refreshes a param entity from its template only when the param type exposes `refresh_from` (`if constexpr (requires …)`). `MobiusParams` lacks it, so the warp's `base` is frozen at spawn and GUI edits to `template_params` never reach the live warp — silently diverging from the live-update contract that `Ripple`/`Noise` honor. *Fix:* give `MobiusParams` a `refresh_from` (and have `step` re-read `base` from it), or document that this transformer's baseline is latched at spawn and requires a respawn.

6. ❌ **`core/mesh.h:429-432` — `MeshOps::compile` reaches a global scratch arena implicitly.** *Rejected:* the `scratch_arena_a` dependency is already documented at the function (`@details` "Borrows scratch_arena_a…" plus the inline remap comment), satisfying the finding's fallback option; the explicit-`Arena&` alternative would let a caller pass the same arena for scratch and `dst`, whose internal `ScratchScope` reclaim then drops the freshly-populated `dst` arrays — a disproportionate footgun for a maintainability nit. Every other `MeshOps`/Conway entry threads its scratch `Arena&` explicitly (the file's stated "no hidden state or implicit arena references" contract), but `compile` allocates its vertex-remap scratch from the global `scratch_arena_a` with no scratch parameter in its signature — invisible at the call site. *Fix:* add an explicit `Arena& scratch` parameter (every call site already has one in hand), or at minimum document the `scratch_arena_a` dependency so the implicit reference is visible.

7. ✅ **`daydream/segment_controller.js:474-476` — `setResolution` mid-render depends on the old generation completing.** It bumps `renderGen` and nulls `results` but never resets `pending`/`frameSeen`/`renderInFlight`, relying entirely on the still-outstanding old-generation render's `frameResolve` to release the in-flight latch. Combined with item 2, a resize during a hung frame is unrecoverable. *Fix:* primarily resolved by the render watchdog (item 2); additionally document that the re-render is deliberately deferred to `tick()` and bounded by that watchdog.

### Low Priority

8. ✅ **`core/filter.h:1077,1102` — `Blur::plot` has the same unwrapped-tap boundary gap as item 3** (`cx = round(x)` un-wrapped, `cx+1 == 2W`). *Fix:* pre-wrap the center, then offset: `int cxw = fast_wrap(cx, W); … fast_wrap(cxw + dx, W)`.

9. ✅ **`core/filter.h:933-935` — `Screen::Trails` overrides `crosses_segments = false` without explaining why.** The override is correct (screen-trail points are seeded from and re-emitted into the same band) but contradicts the `Is2DWithHistory` fail-safe default with no comment, unlike `Blur` which explains its opposite choice. *Fix:* add the one-line rationale comment.

10. ✅ **`core/filter.h:1213-1214` — `Feedback` coarse-grid `cy_hi` clamp is load-bearing on the per-pixel `HS_CHECK`.** Correct as written, but a thin single-coarse-row band relies on the trap rather than producing a valid field. *Fix:* assert `cy_hi >= cy_lo` right after the clamp so a future `ds`/clip change traps at the cheap setup site.

11. ❌ **`core/constants.h:62,68` — `ClipRegion::margin` is a public aggregate member.** `render_x_start/end`'s single-period wrap holds only while `margin < w`, enforced solely in `Effect::set_margin`; a direct `.margin =` assignment bypasses the trap. *Fix:* make `margin` private with a checked setter, or add a defensive `HS_CHECK` in `render_x_start/end` (called once per draw, not per fragment). *Rejected:* render_x_start/end are NOT per-draw — `contains_x` (filter.h:144,168, per plotted fragment) calls them, and they're documented as the per-fragment hot path; an HS_CHECK there violates the no-hot-path-cost rule. Private margin is disallowed (breaks aggregate init), and the `margin < w` invariant is already trapped at its sole production mutation site `Canvas::set_margin` (canvas.h:179).

12. ✅ **`core/animation.h:609` — `ParticleSystem::spawn` silently no-ops when the pool is unbound.** No log/trap, unlike the "pool full" branch two lines down and unlike the fail-fast discipline elsewhere. *Fix:* `HS_CHECK(pool.is_bound(), …)` so a spawn-before-init traps instead of silently dropping particles.

13. ✅ **`core/animation.h:1701-1710` — `RandomWalk` drift coordinate freezes after ~2²⁴ frames.** `static_cast<float>(t) * drift` loses single-frame resolution after ~78 h at 60 fps; the walk is perpetual with no rewind. *Fix:* accumulate a wrapped `float` drift phase instead of multiplying a monotonic frame counter, or document the bound as `Noise` does.

14. ✅ **`core/plot.h:1091-1099` — `sample_closed_ring` lacks `HS_CHECK(num_verts >= 1)`.** `Star`/`Flower`/`DistortedRing` compute `PI_F/num_sides` and then call with a zero/negative side count, silently rendering nothing where `Ring::sample` would trap. *Fix:* add the cold-path `HS_CHECK`, matching `Ring::sample`.

15. ✅ **`core/plot.h:675-677, 266` — two comments overstate invariants.** The `sim_dist >= total_dist` / `scale <= 1` claim is false on the capacity-backstop path (where `scale > 1`), and "the unit tangent regardless of magnitude" is only true within a `map_planar` linear piece, not across a knot. *Fix:* amend both comments to state the two regimes.

16. ✅ **`core/3dmath.h:1062` — `least_parallel_axis` applies a cosine threshold to a possibly-non-unit vector.** A non-unit vector parallel to +X with `|v.x|` below `COS_AXIS_PARALLEL` returns +X, collapsing the caller's `cross(v, X).normalized()` to a trap. All current callers pass unit vectors. *Fix:* make the test scale-invariant (`v.x*v.x > k²·dot(v,v)`) or assert the unit precondition.

17. ❌ **`core/geometry.h:830` (`make_basis`)** — feeds the raw, un-normalized `normal` to `least_parallel_axis`, the same latent trap as item 16. *Fix:* normalize the direction before the parallel test. *Rejected:* subsumed by finding 16. The scale-invariant `least_parallel_axis` now picks the same body axis regardless of `|normal|`, and `make_basis` consumes `normal` nowhere else — `v`, `ref`, `u`, `w` are each independently `normalized()`. A normalize here would be dead work.

18. ✅ **`core/geometry.h:724,751` — `Orientation::upsample` early-return uses the unclamped count.** `if (num_frames >= count) return;` runs before `count` is clamped to `CAPACITY`, so `CAPACITY <= num_frames < count` falls through and redundantly re-resamples on a hot path. *Fix:* clamp `count` before the early-return check.

19. ✅ **`core/color.h:2201-2202` — `BakedPalette::get` interpolates alpha without clamping.** An `AlphaFalloffShade` whose `fn` returns outside `[0,1]` bakes an out-of-range alpha that the linear interpolation extrapolates further, then over/under-blends at the canvas. *Fix:* `hs::clamp(…, 0, 1)` the interpolated alpha (matching `Color4::lerp`), or document that shade `fn`s must return `[0,1]`.

20. ✅ **`core/color.h:843-863,1271,1348` — near-gray palette stops bloom a faint hue.** A stop authored near-but-not-exactly gray keeps a small nonzero `cmax`; hue torsion then introduces a reddish tint (`h≈0`) into midtone grays. *Fix:* zero `colors_cmax[i]` for stops below the same `1e-4` gray threshold `lerp_oklch` uses.

21. ✅ **`core/styles.h:284-286,767` — `hue_rotate` chroma slowly drifts under feedback.** The cached `(ca, sa)` rotation pair is non-orthonormal (the code admits it), scaling chroma slightly; applied per-frame in a feedback loop this compounds. *Fix:* renormalize the cached pair once per frame in `sync_hue` (zero hot-path cost).

22. ✅ **`core/conway.h:892-912` — `relax(0)` returns un-normalized vertices.** The sphere-projection happens only inside the iteration loop, so a caller-settable `iterations == 0` skips the trailing `normalize` every sibling operator applies. *Fix:* `normalize(out_mesh)` before return (idempotent for unit input), or document the pass-through.

23. ✅ **`core/spatial.h:81-83,201-204` — KDTree index bounds are triplicated and one check is dead.** Three uncoupled bounds (`INT16_MAX+1` ctor, `INT16_MAX` node links, `UINT16_MAX` index) with no shared constant; the `UINT16_MAX` guard can never fire. *Fix:* hoist a named `MAX_POINTS` constant referenced by all sites.

24. ✅ **`core/platform.h:553` — `DATA_RATE_MHZ(x)` macro does not parenthesize its argument.** Latent (only ever invoked with the literal `6`) but a macro-hygiene footgun. *Fix:* `#define DATA_RATE_MHZ(x) (x)`.

25. ✅ **`core/presets.h:91-93,103-104` — `current_index`/`prev_index` mix `int` and `size_t`.** `(current_idx + 1) % Size` is `size_t` and narrows back into an `int` member. Functionally safe, but exactly what `-Wconversion` flags. *Fix:* make the index members/accessors `size_t`.

26. ✅ **`core/effect_registry.h:73-77` — `reserve(64)` is a magic number decoupled from `HS_EFFECT_COUNT` (27).** Purely advisory (no element pointers held), so it only stops being a no-realloc guarantee past 64. *Fix:* drop it or document it as a soft hint keyed off the roster size.

27. ❌ **`effects/MindSplatter.h:280` — particle-index bounds use `assert` + a clamp fallback instead of `HS_CHECK`.** *Rejected:* per-fragment hot path; `p_idx` comes from an interpolated float whose rounding overshoot the clamp legitimately absorbs (an `HS_CHECK` would false-trap and cost a hot-path branch); the real invariant is already trapped per-draw at line 297. On device the `assert` is stripped and the clamp masks the invariant violation — the bounded fallback the fail-fast philosophy explicitly avoids. *Fix:* promote to `HS_CHECK`, matching line 297.

28. ✅ **`effects/MobiusGrid.h:55` — `spawn_pinned()` return value dropped.** Sibling `GnomonicStars` `HS_CHECK`s the identical pinned-spawn handle; here it is ignored (cannot fail at `CAPACITY 1` today, but a latent null-handle and a cross-effect inconsistency). *Fix:* capture and `HS_CHECK` the handle.

29. ✅ **`effects/MobiusGrid.h:206` — redundant per-ring basis rebuild.** `draw_axis_rings` recomputes the identical `make_basis(Quaternion(), normal)` for every ring every frame (`normal` is constant). *Fix:* hoist the basis above the draw-curves call.

30. ✅ **`effects/GnomonicStars.h:91` — `points` read with only an upper-bound check.** `(int)params.points` is `HS_CHECK`ed against `MAX_POINTS` but has no lower clamp; a sub-minimum value would no-op or desync the cache. Not animated today. *Fix:* `hs::clamp((int)params.points, 1, MAX_POINTS)`.

31. ✅ **`effects/DreamBalls.h:289` — undocumented sprite/ping-pong overlap margin.** The two-slot `param_slots_`/`baked_palettes_` ping-pong is safe only while at most two sprites overlap (sprite life 320 < 2×period 576). True today with comfortable margin, but the dependency is load-bearing and unstated. *Fix:* add a `static_assert`/comment pinning `sprite_life < 2 * period`.

32. ✅ **`effects/Thrusters.h:178` — thrust geometry uses the previous fire's residual amplitude.** `amp = amplitude` is snapshotted before `warp_anim` is restarted to `0.7`, so the thrust points track the stale warp, not this fire's. *Fix:* sample after the restart (`warp_decay(0)`), or document the intentional residual.

33. ❌ **`effects/RingShower.h:122` — fade-in off-by-one.** *Rejected:* no off-by-one. The ramp is `ease_linear((age+1)/FADE_IN_FRAMES)` for `age+1 < FADE_IN_FRAMES` → 0.25, 0.5, 0.75 at ages 0-2, then 1.0 at age 3 where `age+1 == FADE_IN_FRAMES` (the short-circuit returns 1.0, which equals 4/4). The steps are a uniform 0.25 and the ramp reaches full opacity smoothly; the cited `0.75 → 1.0` step is the same size as every other step, not a discontinuity. The `/4` in the finding is already the named `FADE_IN_FRAMES` constant; the `<= FADE_IN_FRAMES` variant is identical, and the `FADE_IN_FRAMES-1` variant would reach 1.0 a frame early (a behavior change, not a fix). `opacity_at` returns `(age+1)/4` for ages 0-2 then jumps to `1.0` at age 3, producing a visible `0.75 → 1.0` step. *Fix:* divide by `FADE_IN_FRAMES` with `age+1 <= FADE_IN_FRAMES`, or by `FADE_IN_FRAMES-1`.

34. ✅ **`effects/Comets.h:80` — `update_path()` runs before `motion_` is constructed.** Its `if (motion_) reanchor()` is a guaranteed first-call no-op; correct only because `Motion` captures `path` by reference. *Fix:* add a one-line comment noting the order dependency.

35. ✅ **`hardware/dma_led.h:307-309` / `hardware/hd107s_frame.h:305-307` — per-instance color setters mutate per-`N`-global static state.** `setBrightness/Temperature/Correction` look per-controller but write `HD107SFrame<N>` statics, so two controllers of equal `N` would share them. Benign with one controller per image. *Fix:* document the shared scope, or make the correction multipliers non-static if independent per-strip correction is ever wanted.

36. **`hardware/pov_sync.h:1499-1523` — master commit anchors to a self-censored epoch boundary.** When a master `ZERO_EPOCH` symbol is censored (late), `on_epoch_symbol` has already opened the commit window, so the master's commit is anchored one crossing earlier than the boundary downstream boards decode. *Fix:* defer the content-side scheduling until the epoch symbol is confirmed onto the wire — validate against the sync spec before changing.

37. ✅ **`tests/test_harness.h:212-222` (+ standalone `*_check.cpp` mains) — a zero-assertion module reports pass.** Nothing guarantees a registered module actually executed an assertion, so an emptied runner stays green. *Fix:* treat `passed + failed == 0 && skipped == 0` as a module failure.

38. ✅ **`tests/test_memory.h:212-224,577-588` — NDEBUG-gated coverage vanishes silently.** The Arena generation-bump / stale-binding tests compile out under NDEBUG with no `++skipped`; the canonical preset is Debug so they run today, but a Release test build drops them invisibly. *Fix:* emit `++skipped` + a banner in the NDEBUG branch.

39. ✅ **`tests/CMakeLists.txt:146-158,170-180` — standalone check executables are not roster/include-pinned.** `fastmath_clamp_check` / `h_offset_renorm_check` hand-list their `test_*` calls with none of the drift protection `run_tests` gets. *Fix:* assert a minimum expected assertion count in each `main()`.

40. ✅ **`scripts/wasm_smoke.mjs:131` — stack-margin check has no explicit `capacity > 0` assert.** A degenerate zero capacity is caught only indirectly via the `hwm === 0` canary. Self-defending but opaque. *Fix:* add an explicit `capacity <= 0` failure for a clearer diagnostic.

41. ✅ **`scripts/check_screenshots.mjs:23-24` / `scripts/effect_roster.mjs:46` — screenshot freshness compare is case-sensitive.** A casing mismatch is masked on the case-insensitive Windows dev FS but can diverge on Linux CI. *Fix:* normalize-compare and flag case-only differences.

42. ✅ **`daydream/daydream.js:560-561` — WASM-load `.catch` dereferences a not-yet-assigned controller.** `testAllController` is `null` until a later synchronous block; a fast/synchronous module rejection throws `TypeError` and masks the real load error (the sibling `testAllInterval` access is already null-guarded). *Fix:* `if (testAllController) { … }`.

43. ✅ **`daydream/gui.js:35-36,333` — color-string validator accepts alpha forms lil-gui can't round-trip.** 4-/8-digit hex and `rgba()` pass validation but lil-gui is RGB-only and the serializer emits 6-digit, so an alpha deep link mis-hydrates. *Fix:* tighten the regex to `#rgb`/`#rrggbb`/`rgb()` only.

44. ✅ **`daydream/recorder.js:33-47` + `daydream.js:646-647` — explicit MP4/WebM silently falls back.** When the chosen container is unsupported, `selectMimeType` returns `''`, `mimeType` is omitted, and `MediaRecorder` records the browser default with no warning. *Fix:* `console.warn` (or explicit fallback) when an explicit format is unsupported.

45. ✅ **`daydream/state.js:236-237` — `URLSync.reset()` clears the debounce timer but doesn't null it.** Asymmetric with `schedule()`/`dispose()`; harmless but a latent footgun. *Fix:* `this.timer = null;` after the `clearTimeout`.

46. ✅ **`daydream/segment_controller.js:533,541` — compositor pre-pass doesn't reject degenerate rects.** It guards `x1 > w`/`y1 > h` but not `x1 <= x0`/`y1 <= y0`, so an inverted/empty rect yields a zero/negative `expectedLen` that masks layout corruption instead of faulting. `computeSegmentRange` doesn't emit such today. *Fix:* `if (r.x1 <= r.x0 || r.y1 <= r.y0) { onWorkerFault(s, 'degenerate segment rect'); return 0; }`.

47. ✅ **`daydream/segment_worker.js:96` — `init` applies pause via `if (msg.paused)`.** It can only ever enable pause and can't distinguish `false` from omitted, unlike the explicit-boolean `setAnimationsPaused` handler. Benign (engine defaults unpaused). *Fix:* `if (typeof msg.paused === 'boolean') engine.setAnimationsPaused(msg.paused);`.

48. ✅ **`daydream/segment_controller.js:530` vs `665` — compositor and stats loop over two different length sources.** `composite` uses `results.length`, `updateStats` uses `this.count`; equal via `create()` today but a fragile dual source. *Fix:* use one `const n = this.count` (or assert the invariant).

49. ✅ **`daydream/tools/color.js:97-117` — dead OKLCH mirror.** `srgbToOklch`/`oklchToLinearRgb` (and their `oklab↔oklch` helpers) are imported into `palettes.html` but never called — the browser bakes via WASM — yet `color.test.js` pins them with no WASM-parity justification. *Fix:* drop the unused exports/imports and their tests, or add a header note that they are not yet wired.

50. ✅ **`daydream/tests/segment_crosscheck.test.js:52-60` — the C++ reference port mismodels `>2` Y-bands per arm.** `cppSegmentMap` branches only `armSeg === 0` vs else, collapsing every non-zero slot onto one reversed strip; correct only because the sweep caps at `N=4`. *Fix:* guard/assert the `≤2`-bands precondition or compute the band offset generally.

51. **`daydream/tools/shared.js:36-48` — `showFatalError` is untested.** The user-facing WASM-load-failure banner on every tool page is pure, idempotent DOM logic with no test. *Fix:* extract it (or use a DOM stub like the existing tests) and assert one banner with the latest message.

---

## Validated Non-Issues

These candidates were raised during review and **rejected** by independent
validation — recorded so their absence is deliberate, not an oversight:

- **`core/plot.h` adaptive step-cache backstop "fires routinely on polar segments."** Rejected: `screen_step` is floored at `base_step·0.05` and the per-segment count is capped accordingly; tripping the `2*W` cache backstop requires an arc `> 4π` (a sphere-wrapping shape), not normal or polar segments.
- **`core/plot.h:1316` `PlanarPolygon`/`Star`/`Flower` divide by `num_sides` before guarding.** Rejected: the `PI_F/num_sides` value is a harmless IEEE `inf` passed as an *argument*; `Ring::sample` fail-fast traps via `HS_CHECK(num_samples >= 1)` before any use. (The distinct `sample_closed_ring` gap, item 14, is real.)
- **`core/animation.h:2055` `MeshMorph` cross-container loop bound.** Rejected: construction `HS_CHECK`s `mesh_B.vertices.size() == dest.vertices.size()` and binds `end_pos` to that size, so the loop bound is guaranteed in range.
- **`effects/MeshFeedback.h:211` per-morph persistent-arena growth.** Rejected: each morph's completion callback runs `carousel.compact()` + `persistent_arena.reset()`, reclaiming the slot storage every cycle — no unbounded growth.
- **`core/sdf.h` / `core/scan.h`.** A full audit of the rasterizer hot path surfaced no actionable defect: seam wrap, pole degeneracies, CSG span sizing, the Face LUT sign-purity gate, and Volume sphere-trace step safety are all correct-by-construction.

Validation also **downgraded** two items from their initial severity: the
DreamBalls sprite-overlap concern (the original "96-frame margin" rested on
miscounted sprite life; the true margin is large) and the `hd107s_frame.h`
`packPixel` `assert` (live in the canonical Debug test preset, inert only in a
non-canonical Release test build).

---

## Per-Component Grade Summary

| Component | Correctness | Performance | Maintainability | Notes |
|---|---|---|---|---|
| `core/` math, geometry, util | A | A | A− | Unit-precondition traps (16-18) are the only latent edges. |
| `core/` color, palettes, styles | A | A | A− | Color science is a standout; residuals are sub-perceptual. |
| `core/` memory, containers, concepts | A | A | A− | Exceptionally defended; no live defect. |
| `core/` sdf, scan | A | A | A− | No actionable defect after full validation. |
| `core/` plot | A | A− | A− | Cold-path guard (14) + comment fixes (15). |
| `core/` filter, canvas, constants | A− | A | B+ | The AA/Blur wrap-tap inconsistency (3, 8). |
| `core/` animation, transformers | A− | A | B+ | Möbius equator seam + frozen baseline (4, 5). |
| `core/` mesh, conway, hankin, solids, spatial | A | A | A | Hidden-scratch encapsulation break (6). |
| `core/` platform, registry, presets, led | A | A | A− | Macro/typing nits only. |
| `effects/` (all 28) | A− | A− | A | Per-effect guard-style inconsistencies. |
| `hardware/` drivers | A | A | A | Dormant DMA single-board section attr (1). |
| `targets/`, `scripts/`, CMake | A | A | A | Cross-language contract is clean. |
| `tests/` (C++) | A | — | A− | Coverage is broad; zero-assertion blind spot (37). |
| daydream core JS | A | A | A | WASM-bridge handling is exemplary. |
| daydream segmented/workers | A− | A | A− | Missing render watchdog (2). |
| daydream tools + JS tests | A | — | B+ | Dead OKLCH mirror (49). |

---

*Reviewers and validators operated under the constraint of not reading any prior
review document; all findings are derived directly from the source.*
