# Holosphere — Code Quality Review

**Date:** 2026-06-29
**Reviewer:** Expert C++/JS review panel (multi-agent, independently validated)
**Subject:** The Holosphere POV LED engine (C++ firmware + WASM) and the daydream web simulator (JS), reviewed as one product.

## Scope & Method

Every in-scope source file across both repositories was read in full by a dedicated
component reviewer (19 components spanning `core/`, `effects/`, `hardware/`,
`targets/`, `tests/`, and all daydream JS — application, segmented-worker, UI,
geometry tools, and tests). Each candidate finding was then handed to an
**independent validator** that re-read the cited code from scratch and returned
*confirmed*, *adjusted*, or *rejected*. Of **80 candidate findings, 31 were
rejected** by validators (intended design, misreading, already-handled, or not
minimally fixable) and **49 survived** as real, verifiable, single-commit-fixable
defects.

**Out of scope** (excluded per instruction): `core/effects_legacy.h`,
`core/rotate.h`, `targets/Holosphere/Holosphere.ino`, and all vendored/third-party
code (`FastNoiseLite.h`, `inplace_function.h`, `three.js/`, `node_modules/`).

Findings are restricted to *eligible* defects — each is independently verifiable
at a cited file/line and fixable as one minimal commit. Deliberate engineering
philosophies (fail-fast `HS_CHECK` traps, the two-buffer ISR design, always-on
branchless wrap guards, compile-time `<W,H>` specialization, explicit arena
threading) were treated as intended design, not defects.

## Headline

This is a **mature, exceptionally well-engineered codebase** that has clearly been
through multiple prior review passes. The severity profile of the surviving
findings tells the story: **0 critical, 0 high, 2 medium, 47 low**. No
memory-safety landmine, no data race, no correctness break on any shipped
configuration was found anywhere in scope. The defects that remain are
micro-optimizations, documentation drift, latent robustness guards for
not-yet-reachable inputs, and test-coverage gaps. The hard parts — a single-writer
real-time sync protocol, a 16-bit-linear/OKLCH color pipeline, arena lifetime
discipline, and bit-identical sim/device parity — are correct and heavily
defended.

## Quality Dimension Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness & Robustness** | **A−** | Pole/wrap/NaN handling, double-buffer flips, sync-protocol gating, and color math are meticulous and well-tested. The few real issues are narrow edge cases (a NaN window only reachable past the unit-vector tolerance, a 1-LSB rounding bias, an unfolded-coast assert blind spot) and a signed long-run frame-counter overflow. |
| **Memory Safety & Resource Mgmt** | **A** | Arena partition math is wrap-proof; `ArenaVector`/`ArenaSpan` carry debug use-after-free/regrow stamps; lifetimes are RAII-scoped and explicit at every call site. No leak or OOB found in scope; the JS side has a few recorder/offscreen cleanup gaps on rare abort paths. |
| **Concurrency & Real-Time Safety** | **A** | The Phantasm sync core is rigorously single-writer with clearly-reasoned relaxed/release-acquire atomics and IRQ-off brackets; the cycle-counter flywheel avoids 32-bit wrap structurally. The host-testable protocol headers are a model of how to make real-time logic verifiable. |
| **Performance & Efficiency** | **A−** | Hot paths use LUTs, fast approximations, packed `uqadd16`, and gated gamut search. The handful of pessimizations found are off the per-pixel path (redundant `sqrt`/angle recomputes, an O(n) trail eviction that ordering makes pointless). |
| **Architecture & Design (elegance)** | **A** | The compile-time-templated pipeline, the three-domain (World/Screen/Pixel) filter composition, the SDF/scan vs. plot rasterizer split, the Conway/Hankin mesh algebra, and the shared-engine sim/firmware story are genuinely elegant and cohesive. |
| **API / Interface Expressiveness** | **A−** | Factories, explicit narrowing casts, the `Fragment`/shader contract, and the arena-threading convention make intent legible. Minor expressiveness footguns remain (a `float in_frames` that truncates to int, two same-named `phi_to_y` overloads, an empty-tree path that skips a documented precondition trap). |
| **Code Style & Readability** | **A** | Consistent, idiomatic, and uniform across a very large surface; conventions (namespacing, RAII guards, naming) are applied without drift. |
| **Documentation** | **A** | Inline rationale is outstanding — epsilon choices, fail-fast-vs-soft-degrade policy, and the 1-wire sync datasheet are all explained with quantified bounds. A few comments have drifted from the code (a nonexistent `ScratchScope` factory, an overstated chroma-preservation invariant, a self-contradictory "free per frame"). |
| **Error Handling** | **A−** | Fail-fast discipline is consistent and well-placed on cold paths. Gaps are a capacity guard that checks total instead of remaining space, two `Config::valid()` omissions, and JS recorder paths that can surface a recoverable error as fatal or ship a truncated file. |
| **Maintainability** | **A−** | Strong factoring and single-source-of-truth design (param marshaling, segment math, LUT generators). The maintainability findings are latent traps that depend on undocumented invariants (single-TU effect registration, a `kMinLutN>=2` divide precondition, duplicated seam predicates that could silently desync). |
| **Testing & Verification** | **A−** | Coverage is broad and serious — death/trap harness, cross-run determinism, POV tiling proofs, WASM runtime smoke, param-marshal index alignment, parity tests. The gaps are specific untested guard branches (worker length/bounds throws, segment height guard, transformer slot recycling) rather than missing areas. |
| **Portability & Build** | **A−** | CMake-preset driven, cross-platform, with a Teensy size/layout CI gate and a clean native-Clang toolchain. Findings here are hygiene (an unpinned test include list) rather than build breakage. |
| **Overall** | **A−** | A polished, ambitious, real-time graphics system whose remaining defects are overwhelmingly cosmetic or latent. |

## Prioritized Findings

Each item is a validated, in-scope defect fixable as a single commit. Numbering is
sequential across all priorities. `repo/component` and category are tagged inline.

### Medium Priority

1. ✅ **`SmoothUnion::get_vertical_bounds` drops the culled (inverted-bounds) sentinel** — `core/sdf.h:974-983` *(Holosphere · performance)*. When both children are culled (`{y_min=1, y_max=0}`), the `±pad` clamp turns the inverted band into a non-empty `{0, pad}` phantom band, unlike `Union`/`Intersection` which propagate the sentinel. Each phantom row that overlaps the active Y-clip is scanned full-width with per-pixel distance evals that reject everything (~`pad+1`×W wasted evals/frame for an off-screen `SmoothUnion`). Output stays correct; the work is pure waste. **Fix:** before padding, compute `lo=min(...y_min)`, `hi=max(...y_max)`; if `lo>hi` return the inverted `{1,0}` sentinel unchanged.

2. ✅ **`applyResolution(true)` runs inside the WASM-load `.then`, so a post-load render throw is reported as an engine-load failure** — `daydream.js:544-564` *(daydream · error-handling)*. Any throw from `applyResolution`/`applyEffect` after the engine constructed (dot-mesh rebuild, sidebar, bad URL param) propagates into the `.catch` meant only for module-load failure, showing the fatal "Failed to load the rendering engine" overlay for a recoverable error. **Fix:** wrap `applyResolution(true)` in its own try/catch inside the `.then` so the outer `.catch` only fires for genuine module-load rejection.

### Low Priority

#### C++ engine — math, color, memory, raster

3. ✅ **`make_rotation(from,to)` recomputes the angle it already has in hand** — `core/3dmath.h:1083-1099` *(Holosphere · performance)*. `d = dot(from,to)` is already computed and inputs are asserted unit, yet `angle_between(from,to)` re-does a dot, two `sqrtf`, and a divide. **Fix:** `float angle = fast_acos(d);`.

4. **`angle_between(Vector)` uses two `sqrtf` where one suffices** — `core/3dmath.h:937-942` *(Holosphere · performance)*. `v1.length()*v2.length()` is two square roots; `sqrtf(dot(v1,v1)*dot(v2,v2))` is one. Hot in `plot.h`/`sdf.h` per-fragment paths. **Fix:** combine under a single `sqrtf` (keep the epsilon `HS_CHECK`).

5. **`vectorToLogPolar` can return NaN for `v.y` just above 1 that still passes the unit-vector assertion** — `core/geometry.h:525-536` *(Holosphere · correctness)*. `EPS_UNIT_VEC_SQ=0.02` admits `v.y≈1.005`; the magnitude-only north-pole guard misses the resulting negative `denom`, so `0.5*logf(numer/denom)` → NaN, defeating the function's stated finiteness promise. **Fix:** clamp `v.y` to `[-1,1]` (or use signed sentinels for the poles) before the `logf` branch.

6. **`hue_rotate` comments overstate exact L/chroma preservation** — `core/color.h:755-783` *(Holosphere · documentation)*. The fast-trig `(ca,sa)` (Bhaskara I, ~0.17% error) is non-orthonormal, so the rotation scales chroma by `sqrt(ca²+sa²)≠1`; the "preserves L and |(A,B)|" comments are literally inaccurate (negligible in the fade loop). **Fix:** soften the two comments to "within fast-trig accuracy" (do not change the per-pixel trig).

7. **Dead branch in `srgb_to_linear_interp` after the clamp** — `core/color.h:355-368` *(Holosphere · maintainability)*. Post-clamp `f∈[0,255]`, so `if (f <= 0.0f) return lut[0];` can only fire at `f==0`, which the `i=0`/`frac=0` path already yields. **Fix:** remove the redundant early return (keep the required `i>=255` guard).

8. **`ScratchScope` doc claims a temporary-`ArenaVector` factory that does not exist** — `core/memory.h:808-816` *(Holosphere · documentation)*. The class is offset save/restore only; the "@details … typed factory" sentence misleads. **Fix:** delete the factory sentence (or add the method).

9. **`ScratchScope` section banner still says "Factory"** — `core/memory.h:804-806` *(Holosphere · documentation)*. Same drift in the header. **Fix:** rename to "RAII Arena Offset Guard".

10. **`Plot::rasterize` planar pre-pass and draw loop recompute the seam test from different sources** — `core/plot.h:564-573, 744-748` *(Holosphere · maintainability)*. The cached arc-length metric (for `v0`) and the drawn strategy use duplicated, currently-identical seam predicates; a future tweak to one silently desyncs `v0`/`v1` from the rendered position. **Fix:** compute the per-segment seam/planar decision once and consume it in both places.

11. **`Face::build_distance_lut` divides by `(n-1)` on an unasserted `kMinLutN>=2` invariant** — `core/sdf.h:1851-1859` *(Holosphere · correctness)*. Safe today (`n∈[6,32]`), but lowering `kMinLutN` to 1 would divide by zero with no signal. **Fix:** `static_assert(kMinLutN >= 2, ...)` at the definition (do not add a silent runtime clamp).

#### C++ engine — animation, mesh, infra, pipeline

12. **Timeline global frame counter is a signed `int` with long-run overflow UB** — `core/animation.h:2432,2535,2556,2662` (`memory.cpp:126`) *(Holosphere · correctness)*. The per-animation counter was deliberately made `uint32_t` for defined wrap; the shared clock was not, so `++global_timeline_t` reaches signed-overflow UB (~276 days @90fps) and corrupts `< e.start` comparisons. **Fix:** make `global_timeline_t` `uint32_t` (extern + definition + arithmetic).

13. **`Timeline::add`/`add_get` take `in_frames` as `float` then truncate to `int`** — `core/animation.h:2491,2514,2535` *(Holosphere · api-design)*. The float advertises sub-frame precision the integer `e.start` cannot honor; a fractional delay is silently truncated and an out-of-`int` float is UB. No caller passes a float. **Fix:** change the parameter type to `int`.

14. **No test covers non-pinned `Transformer` slot recycling after timeline compaction** — `core/transformers.h:180-188` *(Holosphere · testing)*. The survive-relocation contract (the class's subtle core) is unexercised; existing tests use `spawn_pinned` and never step the timeline. **Fix:** add a test that relocates a non-pinned spawned event via `move_into`, runs it to completion, and asserts the slot is reclaimed.

15. **`classify_faces_impl` hard-traps on a degenerate (zero-length) face edge** — `core/mesh.h:606-621` *(Holosphere · error-handling)*. `normalized()`/`angle_between` trap on coincident consecutive vertices; sibling geometry (`relax`, `conway.h:967`) degrades gracefully via `EPS_LEN_SQ`. Latent (no shipping solid hits it) but inconsistent for a topology *hash* input. **Fix:** compute edges unnormalized, guard each with `dot(e,e)<EPS_LEN_SQ`, and only call `angle_between` when both are non-degenerate.

16. **`KDTree::nearest` skips the documented `k<=MAX_K` trap on an empty tree** — `core/spatial.h:102-108` *(Holosphere · api-design)*. The `root==-1 || k==0` early return precedes the `HS_CHECK(k<=MAX_K)`, so an out-of-contract `k` is caught only when data is present. **Fix:** move the `HS_CHECK` above the early return.

17. **`compile()` can silently corrupt output if a direct caller aliases `scratch_arena_a` as `geom_arena`** — `core/mesh.h:429-432` *(Holosphere · maintainability)*. The `ScratchScope` rewind would reclaim the destination buffers. `generate()` already traps this; direct callers (MeshFeedback, tests) are unguarded. **Fix:** add `HS_CHECK(&geom_arena != &scratch_arena_a && &geom_arena != &scratch_arena_b)` at the top of `compile()`.

18. **`node()` doc comment contradicts itself: "free per frame" for an init-only function** — `core/reaction_graph.h:36-38` *(Holosphere · documentation)*. **Fix:** reword to "off the per-frame render path (zero per-frame cost)".

19. **`Registrar` correctness silently depends on an undocumented single-TU-per-binary invariant** — `core/effect_registry.h:150-171` *(Holosphere · maintainability)*. The anonymous-namespace `_reg` registers once per including TU; a second TU including `effects.h` would double-register and trip the startup count check. **Fix:** add a one-line note that `effects.h` must be included by exactly one TU per binary.

20. **`Screen::Trails` capacity eviction does an O(MAX_PIXELS) shift its own `decay()` makes pointless** — `core/filter.h:976-982` *(Holosphere · performance)*. The FIFO shift (default 1024 elements) preserves an order that `decay()`'s unordered swap-remove already scrambles. The `World::Trails` sibling uses an O(1) ring. **Fix:** replace with an O(1) swap-drop of the oldest-by-index slot.

21. **`World::Trails::flush` guards on an unreachable `item.ttl==0`** — `core/filter.h:758-761` *(Holosphere · maintainability)*. Items are only pushed with `ttl>0` and removed the instant `ttl` hits 0, so the branch is dead and obscures the all-live invariant. **Fix:** drop the branch (or replace with `HS_CHECK(item.ttl > 0)`).

#### C++ effects

22. **`init_lattice()` capacity guard checks total capacity, not remaining space** — `effects/ReactionDiffusionBase.h:204-213` *(Holosphere · error-handling)*. `get_capacity()` ignores bytes already consumed by `cube_lut`/`state`/`palette`, so the documented "under-sized arena" failure is not actually caught here — it falls through to the generic OOM trap. **Fix:** compare `get_capacity() - get_offset() >= RD_N*sizeof(Vector)`.

23. **`MobiusGrid` binds `mobius_gen`'s timeline reference before `timeline` is constructed** — `effects/MobiusGrid.h:33,270,284` *(Holosphere · maintainability)*. Member init order constructs `mobius_gen` (decl 270) before `timeline` (decl 284); well-defined only because the ctor merely stores the reference, and inconsistent with `IslamicStars`'s deliberate ordering. **Fix:** declare `timeline` before `mobius_gen`.

24. **`BZ::blend_species` truncating cast can drop the top palette value** — `effects/BZReactionDiffusion.h:315-320` *(Holosphere · correctness)*. The convex-combination divide can yield `65534.999…` → truncates to `65534`, a systematic 1-LSB downward bias on bright pixels — the same bias `to_q8`/`to_q16` avoid with `+0.5f`. **Fix:** add round-to-nearest (`static_cast<uint16_t>(c + 0.5f)`) per channel.

25. **No `static_assert` binds `MAX_DEGREE` to `factorial()`'s safe float domain** — `effects/SphericalHarmonics.h:19-26,69-74,351` *(Holosphere · correctness)*. A comment invites maintainers to "widen `MAX_DEGREE` alone," but `factorial(2·MAX_DEGREE)` loses float precision at `14!` and overflows to `inf` (→NaN) at `35!`. **Fix:** `static_assert(2*MAX_DEGREE <= 12, ...)` near `MAX_DEGREE` (or rewrite normalization as a cancelling ratio-product).

26. **`DistortedRing` builds a `Mutation` member then copies it into the timeline, leaving a dead duplicate** — `effects/DistortedRing.h:32-38,66,138` *(Holosphere · maintainability)*. Only the timeline copy is stepped; the `amplitude_mut` member is dead state holding a live ref/lambda, inviting "edit the member, see no effect" confusion. **Fix:** construct the `Mutation` inline at the `add` site and delete the member (single owner).

27. **`Comets` palette rebake gate runs one redundant full-LUT rebake per cycle** — `effects/Comets.h:121-124,199-212` *(Holosphere · performance)*. The arming-frame rebake reproduces the existing LUT (the `ColorWipe` doesn't step until the next frame). *Note: the obvious "drop the `+1`" is wrong — it would skip the final-step rebake and leave the palette one step short of target.* **Fix:** start the countdown on the first frame the wipe actually steps, or rebake only when the wipe stepped this frame. Negligible cost; fix for clarity, not speed.

#### C++ hardware, targets, tests/build

28. **`Config::valid()` omits positivity checks on `acquire_quiet_cols` and `beacon_interdigit_timeout_cols`** — `hardware/pov_sync.h:223-249` *(Holosphere · error-handling)*. Both feed `col_cycles()` which casts to `uint32_t`; a zero/negative edit becomes 0 or a ~4-billion-cycle window, breaking the ACQUIRE guard and beacon stale-frame timeout — an inconsistent omission from an otherwise exhaustive boot gate. **Fix:** add `&& acquire_quiet_cols > 0 && beacon_interdigit_timeout_cols > 0` to the conjunction.

29. **`schedule_boundary` conflates wire-busy rejection with late-censor in telemetry** — `hardware/pov_sync.h:966-981` *(Holosphere · correctness)*. Both return `false`, and `master_on_crossing` counts every `false` as `emit_censored`, so a (defensive) wire-busy overlap would be misreported as a late-censor, undermining the "degradation visible in one glance" goal. **Fix:** distinguish `LATE_CENSOR` vs `WIRE_BUSY` (enum/out-param) and route to separate counters (or `HS_CHECK` the unreachable case).

30. **`position()` overflow assert has a blind spot for `elapsed >= 2^31`** — `hardware/pov_sync.h:604-615` *(Holosphere · correctness)*. The `int32` cast wraps negative and slips under the positive threshold (host-only; `NDEBUG` strips it). *Note: the naïve unsigned-cast fix is unsound — it false-traps the legitimate negative-delta snap path.* **Fix:** evaluate the guard on the 64-bit/unsigned elapsed magnitude **before** narrowing, preserving the negative-delta range.

31. **`run_tests.cpp` include list is unpinned against the X-macro roster** — `tests/run_tests.cpp:9-104` *(Holosphere · testing)*. A roster row missing its include is a compile error, but an include left behind after a module is removed compiles silently (header-only `#pragma once`, unused inline), compiling dead test source with no signal. **Fix:** emit the includes from `HS_TEST_MODULE_LIST` via the X-macro, or add a CTest that diffs the include set against the roster.

#### daydream simulator — application & segmented workers

32. **Resolution dropdown not re-synced when the engine rejects a resolution change** — `daydream.js:457,582-584` *(daydream · correctness)*. The global GUI controller is bound to an unretained object literal; after a rejected `setResolution` reverts state, nothing calls `updateDisplay()`, so the dropdown keeps showing the rejected value (latent — both shipped presets are supported). **Fix:** retain the controller and `setValue()`/`updateDisplay()` it in the revert branch.

33. **Aborted `start()` leaves the offscreen canvas latched, so a later `targetHeight` change records at stale dimensions** — `recorder.js:121-191` *(daydream · memory-safety)*. Three abort paths return without `cleanup()`, and `ensureOffscreen` only recreates `if (!this.offscreen)`, so a subsequent `start()` reuses the stale buffer at the old size (plus a small canvas/context leak). **Fix:** call `this.cleanup()` (or null the offscreen/ctx) on every post-acquire abort.

34. **`captureStream(0)` is unguarded and can throw out of `start()`, leaking the offscreen** — `recorder.js:133-140` *(daydream · error-handling)*. It runs after the offscreen is latched but before any try/catch (only the `MediaRecorder` ctor is wrapped). **Fix:** wrap the `captureStream` block; on throw, stop tracks, `cleanup()`, and return.

35. **Streaming `finish()` can download only the post-failure chunk tail, silently producing a truncated video** — `recorder.js:355-382` *(daydream · error-handling)*. On a mid-stream write rejection, `chunks` holds only the failing-and-later chunks; the fallback `download()` builds a blob from the tail alone, and the disk writable is never `close()`d. **Fix:** abort/close the partial writable and surface an explicit error (or buffer all chunks from session start); at minimum `close()` on the failed path.

#### daydream simulator — UI & geometry tools

36. **`prettify` renders small negative values as `"-0.000"`** — `label_format.js:22,39` *(daydream · correctness)*. The zero-snap only catches `|r|<=1e-5`; e.g. `-0.0004` → `"-0.000"` (signed zero) on near-zero axis labels. **Fix:** normalize the `"-0.000"` string to `"0.000"` (or widen the snap epsilon).

37. **`addColor` leaves a rejected invalid URL color in the query string** — `gui.js:308-316` *(daydream · api-design)*. Unlike `add()`'s enum eviction, an invalid color is warned-and-skipped but never removed, so it re-warns every reload. *Note: don't write raw `getValue()` (a `THREE.Color`/array serializes to garbage).* **Fix:** evict with `this.urlWriter(key, null)` (the writer deletes on null), or serialize via the existing onChange logic.

38. **`setEffects` schedules a `requestAnimationFrame` that `dispose()` does not cancel** — `sidebar.js:133,81-86` *(daydream · memory-safety)*. A dispose in the same frame leaves the queued `updateScrollArrows()` to run against a torn-down instance — benign today but contrary to dispose()'s stated contract. **Fix:** store the RAF handle and `cancelAnimationFrame` it in `dispose()`.

39. **`isValidColorString` regex admits malformed numeric components** — `gui.js:33-36` *(daydream · correctness)*. The `[\d.]+` per-component pattern accepts `rgb(1.2.3,4,5)` or a lone `.`, which then reaches the color parser the validation was meant to gate. **Fix:** use a proper number pattern `(?:\d+(?:\.\d+)?|\.\d+)` (optionally range-restricted).

40. **`hankin` op can emit an invalid C++ identifier from a negative angle** — `tools/solid_codegen.js:117-120` *(daydream · correctness)*. `_hk${round(angle)}` with a negative angle yields `base_hk-30`; every other suffix routes through `pctSuffix()` which rejects negatives precisely to keep a valid identifier. **Fix:** guard the hankin branch with the same non-negative check before appending `_hk`.

41. **`createSlider` does not clamp the initial value into `[min,max]`** — `tools/slider.js:84-90` *(daydream · correctness)*. An out-of-range `value` makes the native range input clamp `slider.value` while the readout span shows the unclamped value, so control and label disagree until first input. **Fix:** clamp `value` into `[min,max]` and seed both the slider and the readout from the clamped value.

42. **`formatFloatCpp` emits scientific notation for magnitudes `>= 1e21`, breaking the C++ literal** — `tools/cpp_format.js:25-42` *(daydream · correctness)*. `toFixed` only avoids exponent form below `1e21`, contradicting the "always plain decimal" contract (`(1e21).toFixed(6) === "1e+21"`). **Fix:** at `|n| >= 1e21` throw a clear range error (consistent with the existing non-finite fail-fast) or format via a non-exponential path.

43. **`initScene` `dispose()` does not free the scaffold's sphere geometry/material or clear the scene** — `tools/shared.js:185-192` *(daydream · performance)*. `renderer.dispose()` doesn't release per-geometry/material GPU buffers it created; minimal real-world impact (callers dispose once on `pagehide`) but an ownership-symmetry gap. **Fix:** add `sphere?.geometry.dispose(); sphere?.material.dispose(); scene.clear();`.

#### daydream simulator — test coverage

44. **`computeSegmentRange`'s height-vs-bands guard is never exercised** — `tests/segment_layout.test.js:65-75` *(daydream · testing)*. No test passes a valid even total with `h < total/2` to hit the "height must be >= N y-segments per arm" throw. **Fix:** add `computeSegmentRange(0, 8, 96, 3)` expecting the `/y-segments per arm/` throw.

45. **Worker render-time length & bounds invariant throws are untested** — `tests/segment_worker.test.js:144-170` *(daydream · testing)*. `FakeEngine.getPixels` always returns the exact-size buffer, so the only protection against a silently-short buffer (the throw — `extractSegment` clamps rather than throws) has no test. **Fix:** override `getPixels` to return a wrong-length buffer and assert the rethrow throws `/pixel buffer length/`.

46. **Init `paused` branch in the worker is never covered** — `tests/segment_worker.test.js:116-129` *(daydream · testing)*. `init` with `paused:true` calls `setAnimationsPaused(true)`, but tests never pass it, so a regression resuming animation on a paused-pool rebuild wouldn't be caught. **Fix:** dispatch `init` with `paused:true` and assert `paused === true`.

47. **`engine_contract_wasm` pins `getParamValues` existence but not its array-of-numbers shape** — `tests/engine_contract_wasm.test.js:22-31` *(daydream · testing)*. The worker/controller depend on it being array-like of numbers (real binding returns a `Float32Array` view, `FakeEngine` a plain array), but the contract test never invokes it. **Fix:** call it after `drawFrame()` and assert `.length` is numeric and every element is a number.

48. **`setParameter` contract is skipped when an effect exposes no params** — `tests/engine_contract_wasm.test.js:44-57` *(daydream · testing)*. The pin is gated behind `if (defs.length > 0)`; a zero-param default silently skips it (reads as passing). **Fix:** assert `defs.length > 0` for the bootstrap effect, or pick an effect known to expose a parameter.

49. **`cpp_format` default-precision rounding path is never pinned** — `tests/cpp_format.test.js:19-38` *(daydream · testing)*. Tests cover whole numbers, trimming, small-value escalation, and exponent avoidance, but never a value that actually rounds at the default 6 digits. **Fix:** add `assert.equal(formatFloatCpp(0.12345678), '0.123457f')` to lock default rounding.

---

*Methodology note: 80 candidate findings were produced by 19 component reviewers;
31 were rejected by independent validators (intended design, misreading,
already-handled, or not minimally fixable without a performance regression). The
49 above are the survivors. No finding rises above "medium," and the two mediums
are a bounded wasted-scan and a misleading-error-overlay — neither affects output
correctness on a shipped configuration.*
