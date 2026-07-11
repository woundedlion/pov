# Holosphere / daydream — Code Quality Review

**Date:** 2026-07-10
**Repositories reviewed:** `Holosphere` (C++ rendering engine + Teensy firmware + WASM target) and `daydream` (Three.js web simulator + geometry tools).
**Method:** Orchestrated multi-agent review. Fourteen sub-agents each took one component (engine infrastructure, math/geometry, color/shading, mesh, render rasterizers, animation, effects ×2, hardware drivers, targets/build, C++ tests, daydream core JS, daydream workers/tools, daydream JS tests). Every agent read its files in full, cross-checked callers, **validated each finding against the code before reporting it**, and graded its scope. The orchestrator reconciled the results and applied global grades.
**Out of scope (per review instructions):** `core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and all vendored code (`core/vendor/`, `three.js/`, `node_modules/`, the generated `holosphere_wasm.js`).

---

## 1. Executive Summary

This is an exceptionally well-engineered, mature codebase. Across ~90 in-scope source files and ~85 test files spanning two languages, three build targets, and a hard-real-time firmware layer, the review surfaced **no Critical and no High-severity defects** — no memory-corruption bug, no data race, no incorrect result on a live path. The findings are **6 Medium** items (latent footguns, a preview-only visual artifact, and one high-value test-coverage gap) and roughly **60 Low** items (silent-degradation-over-trap preferences, cross-file duplication, test seams, and comment verbosity).

The engineering culture on display — fail-fast invariants verified by a SIGILL death harness, compile-time capacity/footprint budgets, single-source-of-truth X-macro rosters cross-checked in both directions, host-testable pure cores split cleanly from un-mockable device layers, and C++↔JS parity pinned by golden fixtures — is at or above professional production standard. The dominant weaknesses are *maintainability-of-consistency* (a handful of hand-maintained patterns that could be derived, most notably the `full_frame`/`any_crosses_segments` coupling) and *comment discipline* (rationale essays that run against the repo's own terse-comment convention).

**Overall grade: A−** (top of its class; held just below A by API sharp-edges in the mesh layer, the bimodal effects-test coverage, and cross-cutting duplication/verbosity).

---

## 2. Grade Summary

| # | Quality Dimension | Grade |
|---|---|---|
| 1 | Correctness & Robustness | **A** |
| 2 | Architecture & Design | **A** |
| 3 | Interface Expressiveness / API Design | **A−** |
| 4 | Architectural Elegance *(subjective)* | **A** |
| 5 | Memory Safety | **A** |
| 6 | Error Handling / Fail-Fast Discipline | **A** |
| 7 | Performance & Efficiency | **A** |
| 8 | Concurrency & Real-Time Correctness | **A−** |
| 9 | Testing & CI | **A−** |
| 10 | Documentation | **A** |
| 11 | Code Style & Readability | **A−** |
| 12 | Maintainability | **A−** |
| 13 | Portability | **A** |
| | **Overall** | **A−** |

---

## 3. Dimension Assessments

**1. Correctness & Robustness — A.**
No functional defect was confirmed on any live path across fourteen independently reviewed scopes, including dense fixed-point color math, singularity-laden spherical geometry, half-edge mesh operators (all capacity bounds recomputed by hand and found exact), a scanline rasterizer stress-traced through its wrap/pole/AA seams, and a cycle-counter flywheel sync protocol. Every singularity (pole, antipode, degenerate cross-product, near-parallel slerp) has a named threshold and a documented fallback. The Medium items are latent (guarded by current call sites) or preview-only. This is the strongest dimension.

**2. Architecture & Design — A.**
The compile-time-`<W,H>` engine, the three-domain filter pipeline, the arena memory model with explicit `Arena&` threading, the CRTP animation families, and the clean split between pure host-testable cores (`pov_sync.h`, `param_marshal.h`, `segment_layout.js`) and un-mockable device/DOM layers are all principled and consistent. Single-source-of-truth X-macros (`HS_EFFECT_LIST`, `HS_RESOLUTIONS`, `MESHOP_LIST`) generate registration, dispatch, and bindings from one declaration and are cross-checked against drift in both directions.

**3. Interface Expressiveness / API Design — A−.**
Mostly excellent: owning/borrowing split (`Fn` vs `FunctionRef`), deleted-rvalue borrow contracts that turn use-after-free into compile errors, RAII scopes (`ScratchScope`, `Persist`, correction guards), the fluent `SolidBuilder`, and the zero-overhead `StaticPalette` composition. Held below A by genuine sharp edges: the `MeshState` owned/borrowed dual-mode pushes manual view-clearing onto every caller, the `SDFShape` concept under-specifies its own scan contract, a few APIs default to the lossy/dangerous path (`inv_gnomonic`, `Style` copy), and the effects layer restates a derivable pipeline trait by hand.

**4. Architectural Elegance *(subjective)* — A.**
Several designs are genuinely elegant: the one-wire, count-coded, odd-only distance-2 sync alphabet with its "fail-to-missed-never-wrong" property; the once-per-frame `orientation_id` motion-blur collapse; linear-light blending with `explicit`-only gamma egress; the quarter-turn `cos = sin[x+W/4]` LUT; and the `sdf_max_spans` compile-time span-budget traits. The least elegant seams are the multi-virtual `IAnimation` surface and the double-duty segment `results` array.

**5. Memory Safety — A.**
Arena bump-allocation with wrap-proof subtractive bounds math, `n*sizeof(T)` overflow traps on every allocation entry point, a dual-stamp (arena-generation + rebind-counter) use-after-free detector for `ArenaSpan`, a `static_assert` enforcing the trivially-destructible arena contract, compile-time per-effect footprint budgets against the *device* arena literal, and — on the JS side — best-in-class handling of the zero-copy WASM view-detach hazard plus symmetric Three.js resource disposal. The one soft spot is the unenforced borrowed-view lifetime in the mesh transform path (Fix #1).

**6. Error Handling / Fail-Fast — A.**
`HS_CHECK` is placed exactly per the stated doctrine — cold seams only (allocation, config, registration, commit deadlines), never the per-pixel loop, which uses a stripped `assert` backed by a cold trap at the bind site. The transient-vs-invariant distinction is applied deliberately (DMA overrun drops a frame; arena OOM traps). The death harness proves 41 of these invariants actually trap. Minor deductions for a few APIs that prefer silent degradation to a trap (Fixes #17–#20, #22–#25).

**7. Performance & Efficiency — A.**
Split trig LUTs, fast-approx trig confined to hot paths with exact trig at setup, no-heap render path, branchless per-segment ISR packing, zero-copy WASM→GPU color, `noinline` placement to dodge a measured register-spill cliff, warp-field caching with pure-function invalidation, and honest, documented cost accounting (Crossfade's two rasterized meshes/frame, MeshMorph's O(V²) constructor). Minor per-frame JS allocations under specific UI states.

**8. Concurrency & Real-Time — A−.**
The flywheel timebase (column derived from `DWT->CYCCNT`, not interrupt counts), the single-writer ISR ownership model with one fused synchronized seam, and the three-index double-buffer atomics are all sound and carefully annotated. Held below A because the device concurrency *shell* — the release/commit handshake that gates against a use-after-free of a deleted effect — is validated only on hardware (Fix #6), and the JS worker pool has a two-generation overrun re-blit artifact (Fix #3).

**9. Testing & CI — A−.**
Near-complete unit coverage with drift-proof rosters, a best-in-class shell-free SIGILL death harness with a minimum-case floor, self-testing determinism (globals deliberately dirtied between captures), WASM↔C++ parity pinned by absolute goldens, and three CI layers (pre-commit hook, presubmit on Linux+Windows with runtime WASM smoke, gated Pages deploy). Held below A by the bimodal effects coverage (only 4/14 in group 1 expose white-box seams), the untested device concurrency shell and embind write seam, and a pre-commit gate that runs strictly fewer frames than CI (Fixes #6–#16).

**10. Documentation — A.**
The 2,172-line README is a masterwork — architecture, data-flow diagrams, a full 1-wire signal datasheet with AC timing, and rationale for every philosophy. Per-function contracts explain *why*, and cross-target divergences are pinned to bit-exactness. Deductions only for localized drift (IslamicStars params, ID-strap N=8 wording) — Fixes #68, #69.

**11. Code Style & Readability — A−.**
Uniform naming, small focused functions, clear hot/cold separation, and the JS no-underscore convention are all honored. The single pervasive knock — flagged independently by nearly every agent — is comment *verbosity*: essay-length rationale blocks that occasionally dwarf and obscure the code they annotate, running against the repo's own terse-comment style (Fix #66).

**12. Maintainability — A−.**
Strongly resisted drift via X-macro SSOTs, CMake↔roster pins, and compile-time budgets. Dragged down by cross-file duplication (shrunk-face loops, wipe-rebake state machines, footprint boilerplate, `MAX_RINGS`), the virally-repeated `MeshState` view-clearing, and the hand-maintained `full_frame` coupling that has no cheap regression test.

**13. Portability — A.**
The `platform.h` abstraction, pinned ARM cross-toolchain and native-Clang toolchain file (transparently handling Windows `lld-link`/`rc.exe`), CMake presets, and host/device/WASM triangulation make the build host-independent. Minor cross-repo coupling in the screenshot capture script (Fix #65).

---

## 4. Component Grade Matrix

Selected dimensions per reviewed scope (— = not separately graded for that scope):

| Scope | Correct. | API | Test | Perf | Maint. | Elegance |
|---|---|---|---|---|---|---|
| core/engine infrastructure | A | A | A− | A | A | A |
| math / geometry / transformers | A | B+ | A− | A | A− | A |
| color / shading | A | A− | A− | A | A− | A |
| mesh subsystem | A | **B** | A− | A | B+ | A− |
| render rasterizers | A− | B+ | B+ | A | B+ | A− |
| animation | A | A− | A | A− | A− | A− |
| effects (group 1) | A− | B+ | **B−** | A− | B | A− |
| effects (group 2) | A | A− | A− | A | B+ | A− |
| hardware drivers | A− | A− | **B** | A | A− | A− |
| targets / build scripts | A | A− | B+ | A | A− | A− |
| C++ test suite + CI | — | — | A/A+ | — | A | — |
| daydream core JS | A− | A | A | A− | A− | A− |
| daydream workers / tools | A− | A | A | A− | A− | A− |
| daydream JS test suite | — | — | A− | — | B+ | — |

---

## 5. Prioritized Fix List

Every validated defect is listed below, numbered sequentially and grouped by priority. Severity of each finding is shown in brackets. No Critical or High findings exist; Priority 1 collects the Medium-severity items plus the one cross-cutting correctness footgun.

### Priority 1 — Correctness & robustness (address first)

1. ✅ [Medium] **Mesh — unenforced borrowed-view lifetime.** `MeshOps::transform` populates the output `MeshState`'s `*_view` spans to point into the *input* mesh's arena storage, so the output aliases memory it does not own (`conway.h:247`, `spatial.h:285`). Correct only while the input stays alive; a caller that rewinds/reuses the input arena gets silent dangling reads. Add a debug-only provenance tag or a type-level lifetime contract.
2. ✅ [Medium] **Render — solid+stroke CSG renders wrong AA with no guard.** `SDF::Union`/`SmoothUnion` compute `is_solid = A::is_solid || B::is_solid`, so a stroke winner is shaded through the solid AA branch (`sdf.h:824`, `sdf.h:940`). `Union<Face, Ring>` compiles and mis-renders. Add `static_assert(A::is_solid == B::is_solid)` so the unsupported mix fails loudly.
3. ✅ [Medium] **daydream workers — overrun re-blit composites two generations.** `renderParallel()` deliberately does not clear `results`, so during an in-flight generation the array mixes new and leftover quadrants; on a render overrun `tick()` re-blits the mix, tearing quadrants across sim frames (`segment_controller.js:943-949`). Preview-only, but violates the stated "cached results stay valid" invariant. Stage frames into scratch and publish atomically at `pending === 0`.
4. ✅ [Medium] **Effects (cross-cutting) — `.full_frame` / `any_crosses_segments` coupling is manual and inconsistently applied.** Some effects pass `.full_frame = decltype(filters)::any_crosses_segments`; others (`MobiusGrid`, `Moire`, `DistortedRing`, `GnomonicStars`, `Voronoi`, `HankinSolids`, `Raymarch`) omit it. Every omission is *currently* correct, but adding a history-bearing filter to an omitting effect silently mis-renders on the segmented Phantasm driver, with no regression test. Derive the flag from the pipeline type instead of restating it by hand.
5. ✅ [Medium] **Mesh — fragile `MeshState` owned/borrowed dual-mode.** Discrimination on `is_bound()` is subtle and every consumer must manually null stale views (`transform` ×2, `update_hankin`, `clone`, move ctor/assign) — one forgotten `clear_views()` reads a stale borrowed span (`spatial.h:355-406`). Provide a single `set_owned()`/`set_view()` mutator that always clears the other mode.

### Priority 2 — Test coverage gaps & latent footguns

6. ✅ [Medium] **Hardware — segmented device concurrency glue untested on host.** The `release_req_`/`release_ack_` teardown handshake, the acquire/release adoption of `pending_effect_`, clip alternation, and `dark_latched_` retry — the code whose commit-time `HS_CHECK` is the only guard against a use-after-free of a deleted effect — run only on hardware. Extract the handshake state machine into a pure helper the way `pov_sync.h` was.
7. ✅ [Medium] **Tests — local pre-commit gate is strictly weaker than CI.** The hook runs `DEFAULT_FRAMES=8`, but the lifecycle transitions the harness targets (RingShower slot reuse, Thrusters FIFO expiry, ShapeShifter 48-frame cut) surface only at `HS_SMOKE_FRAMES=120`, which runs only in CI; death tests also skip locally if self-spawn is unavailable. A green pre-commit can still break CI — add a nightly/opt-in long run or document that pre-commit is a fast subset, not the gate.
8. ✅ [Low] **Effects (group 1) — bimodal test coverage.** Only 4/14 effects expose white-box friend seams; Voronoi coherence classification, Hopf projection, and Dynamo's band walk have no direct test. Add seams for the non-trivial numeric paths.
9. ✅ [Low] **Color — non-segue `shade_mesh_topology` overload untested.** It is the production path for IslamicStars (gain 1.0) and HankinSolids; only the segue overload is pinned (`shading.h:102`). Add a golden for the direct `c.alpha = opacity` path.
10. ✅ [Low] **Render — per-pixel `distance()` hot path has no cheap self-test.** Correctness leans on pinned end-to-end checksums; a targeted per-shape distance unit test would localize regressions.
11. ✅ [Low] **led.h — correction-guard RAII happy path untested.** Only the double-construct *trap* is exercised; nothing asserts a normal scope exit balances `correction_guard_depth()` to 0 and restores the baseline (`led.h:70-135`). The guards are host-constructible, so a depth-balance test is feasible.
12. ✅ [Low] **Targets — `wasm_smoke.mjs` never exercises the embind write seam.** `setParameter`, `setClip`, and a rejected `setResolution` are covered only by host predicate tests, not end-to-end through embind, so a binding-signature drift on those methods ships unseen. Add a clamp-readback and an in/out-of-range `setClip` assertion.
13. ✅ [Low] **daydream — `driver.js` `LabelPool` and `coordsLabel` untested.** `LabelPool` backs the "zero allocation per frame" invariant (pool reuse/exhaustion is a real concern) and `coordsLabel` is a pure formatter; both are readily unit-testable and currently uncovered.
14. ✅ [Low] **daydream workers — no test asserts the overrun re-blit produces a non-mixed frame.** This gap is exactly why Fix #3 slipped through; add the assertion alongside the fix.
15. ✅ [Low] **Tests — determinism check compares only the final frame.** `render_capture` copies out only the last buffer, so mid-run nondeterminism that reconverges is invisible (`test_effects.h:293-297`). Fold a per-frame checksum instead.
16. ❌ [Low] **Animation — MeshMorph builds its O(V·V) nearest-vertex map synchronously in the constructor**, i.e. on the swap frame (`mesh.h:98-114`), a per-transition frame-time spike on device with no budget cap. Consider precomputing off the critical frame or capping cost. _Rejected as disproportionate: a one-time per-transition cost with small V; precompute-before-swap is infeasible (the dest mesh is unknown until the swap) and a hard cost cap would alter the morph correspondence._

### Priority 3 — Robustness hardening (prefer trap/guard over silent degradation)

17. [Low] **Geometry — `inv_gnomonic(z, original_sign = 1.0f)` defaults the lossy path.** Forward `gnomonic` collapses antipodes; a caller who forgets the sign silently gets the northern hemisphere for every southern point (`3dmath.h:793`). Make the hemisphere argument required.
18. [Low] **Geometry — `Style` copy silently drops the `noise` binding + hue cache.** `Style` is trivially assignable, but a full-struct copy (`style = Style::Churn()`, `Presets::apply`) resets `noise` to `nullptr` (degrading `noise_warp` → identity) until `sync_noise()`/`sync_hue()` re-run — an unwarped, un-hue-shifted frame with no error (`styles.h:94`).
19. [Low] **Color — `gamut_clip_preserve_chroma` hinges on an unenforced `L ∈ [0,1]`.** For L slightly out of range even the achromatic floor is out of gamut and the search returns a silently out-of-gamut color (`color.h:694`). Add a guard.
20. [Low] **Color — `GenerativePalette::get` re-derived chroma not floored at 0.** A small negative `C` from `fast_sinf` flips hue 180° in `oklch_to_oklab` (`color.h:1401`). Mirror the `std::max(0.0f, …)` already used in `lerp_oklch`.
21. [Low] **Render — `AngularRepeat` axis unit-length precondition unchecked/un-normalized.** A non-unit axis silently distorts the fold (`sdf.h:1459`). Add `HS_CHECK(|ax| ≈ 1)` or normalize.
22. [Low] **Render — dead/fragile stroke-AA `size <= 0` fall-through** leaves `alpha = 1.0` with no falloff (`scan.h:79-82`). Unreachable today, but any future zero-size stroke mis-AAs; trap or set `alpha = 0`.
23. [Low] **Engine — `shortest_distance` lacks the NDEBUG divide-by-zero floor its `wrap` sibling has** (`util.h:103-108`); under NDEBUG a non-positive `m` yields NaN. Mirror the `fmax` floor or `HS_CHECK` the precondition.
24. [Low] **Engine — `ClipRegion` defaults `w`/`h` to compile-time `MAX_W`/`MAX_H`, not the active canvas** (`constants.h:29-35`); latent if any path queries the region before the driver populates it. Make the dims non-defaulted or document the ordering.
25. [Low] **Hardware — `POVSegmented`/`phantasm_config` lack `static_assert(RPM > 0)`.** `COLUMN_US` divides by `RPM`; `RPM = 0` yields `inf`, which passes the `>= 1.0` oversample guard (`pov_segmented.h:75`,`:125`). `POVDisplay` already has the guard.
26. [Low] **Targets — `GradientShape` enum integer values are unpinned across the JS ABI.** `bakeLut` clamps to `[STRAIGHT, FALLOFF]` but nothing asserts `STRAIGHT == 0 && FALLOFF == 3` (`wasm.cpp:1173`); a reorder in `color.h` silently remaps shapes. Add the `static_assert`.
27. [Low] **Targets — `generate_reaction_graph.py` doesn't guard `RD_N-1 ≤ 32767`** for the emitted `int16_t neighbors` table (`generate_reaction_graph.py:106`); a future `RD_N` past 32767 narrows indices silently. Add the assert.
28. [Low] **Render — `SDFShape` concept under-specifies the scan contract.** It requires only `is_solid`/`thickness` but the rasterizer also needs `distance<>()`, `get_vertical_bounds<H>()`, `get_horizontal_intervals<W,H>()` (`sdf.h:307`); a shape missing those fails deep in template instantiation. Widen the concept.
29. [Low] **Render — world-aware clip cull gated by an unenforced `requires`-match** can silently degrade to raw-geometry culling if a real pipeline's signature drifts, dropping geometry an `Orient` moves into a band (`plot.h:712`, `filter.h:261`). Add a `static_assert` on concrete pipeline types.
30. [Low] **Animation — MeshMorph per-frame `slerp` can approach the antipodal degenerate** on a sparse-source nearest match; the dot-max match minimizes but does not bound the angle (`mesh.h:135`). Low likelihood; the exact-antipode gate handles only the exact case.
31. [Low] **Animation — `add_get` defaults `pin = true`.** Pinning a *finite* animation is caught only at completion, not at add, and the `true` default silently marks events non-relocatable (`timeline.h:178`). Assert `!is_finite()` at add for pinned events.
32. [Low] **daydream — session/action GUI controls serialize to the shared URL.** `Test All`, `Segmented POV.*`, and `Recording.*` value controls persist to the query string, so a copied link auto-starts the cycler / spawns a worker pool on load (`daydream.js:631-687`). Bypass `attachUrlWriter` for session controls.
33. [Low] **daydream — `gui.js` `applyOnLoad` latch re-fires on every later `onChange`** registration, diverging from lil-gui's contract that registration does not fire the handler (`gui.js:176-191`); a footgun for the tool pages that reuse the wrapper.
34. [Low] **daydream — one-frame stale blend on the heap-growth frame.** When `Daydream.pixels` detaches, `fill(0)` is correctly skipped but no JS-side clear runs that frame, so an additive/persist effect can blend one stale frame (`driver.js:485-486`). Rare and self-healing; note it explicitly at the call site.

### Priority 4 — Duplication, consistency & maintainability

35. ✅ [Low] **Mesh — "shrunk primary face" re-walk copy-pasted across `expand`/`chamfer`/`snub`** (~60 lines, `conway.h:691-736`,`:776-808`,`:1006-1033`); differs only by the vertex-position formula. Factor an `emit_shrunk_face(...)` helper.
36. ✅ [Low] **Mesh — two `transform` overloads duplicate view-setup/vertex-bind** (`conway.h:247-305`); the zero-transformer base is a prefix of the variadic body and could delegate with an empty fold.
37. ✅ [Low] **Mesh — `Candidate` duplicates `Neighbor` in `nearest()`** (`spatial.h:116-163`): collect-sort-rebuild in two passes when `Neighbor` already carries the fields. Remove the intermediate type.
38. ✅ [Low] **Render — `apply_pole_containment` duplicates the N/S point-in-polygon block** (`sdf.h:2266-2307`), two ~18-line loops differing only by `center.y` sign. Factor a helper.
39. ✅ [Low] **Effects — duplicated engine-level facts across effects:** the `full_frame` restatement (Fix #4), `MAX_RINGS = H>1?H:1` (Moire + ShapeShifter), the footprint-budget boilerplate (`GnomonicStars`/`SplineFlow`/`HopfFibration`/`Dynamo`), and `radius_at()` growth math (RingShower + Thrusters). Hoist to engine helpers. Scoped to the budget boilerplate only: hoisted `DEVICE_PERSISTENT_BUDGET` (memory.h), collapsing the repeated split expression across all 8 effects that carried it. The others left as-is: `full_frame` is already derived where cheap (Fix #4) with no forget-proof fix short of the rejected CRTP; `MAX_RINGS` and `radius_at` are 1-2 line idioms whose hoist would couple self-contained effects for negative readability.
40. [Low] **Effects — Comets and MobiusGrid duplicate the wipe-rebake state machine** (`wipe_pending_`/`wipe_frames_remaining_` gate, `Comets.h:118-123`, `MobiusGrid.h:82-87`). Factor a "rebake-while-wiping" helper.
41. ✅ [Low] **Effects — Liquid2D and DreamBalls reinvent `register_animated_param`** via `register_param` + a `mark_animated` string-loop that re-lists every param name (`Liquid2D.h:56-59`, `DreamBalls.h:74-75`). Use the single call.
42. ✅ [Low] **Effects (group 1) — `std::` vs `hs::` scalar helpers mixed within single files** (Voronoi uses both `std::clamp` and `hs::clamp`; DistortedRing mixes `std::abs`/`hs::clamp`). Standardize on the `hs::` surface. Converted Voronoi's `std::clamp` to `hs::clamp`; `std::min`/`std::max`/`std::abs` have no `hs::` equivalent (`hs::` provides only `clamp`/`lerp`, and `std::abs` is the engine-wide idiom across 6+ effects), so those remain.
43. [Low] **Effects (group 1) — `Params`/sub-struct visibility and field order inconsistent.** Public-at-top vs private-at-bottom varies; HopfFibration's struct field order differs from its registration order (`HopfFibration.h:95-101` vs `:39-43`). Cosmetic but invites a wrong "fix."
44. [Low] **Effects (group 1) — DistortedRing's amplitude wave ignores the pause gate** (`DistortedRing.h:57`); "Pause Animation" cannot hold the rings still as it can in every sibling. Gate the mutation or document the absence.
45. [Low] **Effects (group 1) — Voronoi dead `sharpness > 0` guard** (`Voronoi.h:105`): the slider floor is `1.0`, so the branch is unreachable. Lower the floor or drop the comparison.
46. [Low] **Effects (group 1) — Moire dead struct-default silently duplicated** with the `W<=96` branch value (`Moire.h:169`,`:56`); if the resolution branch is removed, the stale default becomes live with no warning.
47. [Low] **Effects (group 2) — MindSplatter vestigial `opacity` parameter** always passed `1.0` (`MindSplatter.h:247`). Per the repo's remove-dead-code-over-documenting rule, drop the parameter.
48. ✅ **Effects — Dynamo `color()` is O(boundaries) per plotted/trail point** (`Dynamo.h:171-214`); bounded (≤15) but the one non-lookup per-fragment path. Add a `size == 0` fast path or note the typical count.
49. ✅ [Low] **Color — operator asymmetry.** `Pixel16` has `operator*` but no `operator+`; `Color4` has `operator+=`/`*=` but no free `+`/`*`. Callers can only accumulate in place. Round out the set.
50. [Low] **Color — split the 2,289-line `color.h`.** The modifier/`StaticPalette` composition layer (~900 lines) is logically separable from the color-math core and would ease navigation.
51. [Low] **Animation — motion-blur virtuals leak into the universal interface.** `collapse_orientation`/`orientation_id` sit on `IAnimation` but only ~2 of ~20 families implement them (`animation.h:87-96`). Isolate via an `OrientedAnimation` mixin.
52. [Low] **Animation — single-slot `.then()` traps a second chain** (`HS_CHECK(!post)`, `animation.h:178`) rather than composing; limits sequencing expressiveness.
53. [Low] **Geometry — `angle_between` overloads have asymmetric accuracy.** The `Vector` overload uses `fast_acos`; the `Quaternion` overload uses exact `acosf` (`3dmath.h:927`,`:1018`) — same name, silently different precision. Rename one or share a policy parameter.
54. [Low] **Geometry — transformer composition order is slot-index order, not spawn order** (`transformers.h:90-110`); recycling a freed low slot inserts a temporally-later, non-commutative warp earlier in the composition, making the result depend on slot-reuse history. Document or order by spawn time.
55. [Low] **Mesh — one-shot `hankin()` reverses arena polarity** so `target` must hold `max(compile-scratch, output)` (`hankin.h:364-389`), unlike every other operator; a caller sizing `target` for the output alone under-provisions and traps deep inside `compile_hankin`. Assert the headroom at the call boundary.
56. [Low] **Mesh — composition-polarity input/output arena sharing has no static capacity guard** (`solids.h:218`, `conway.h:324`); a new recipe combining a composed op with a high-expansion successor could overflow only at runtime.
57. [Low] **Geometry — `Style` downsample lerp snaps to the denser grid at t=0.5**, spiking coarse-grid scratch demand (~16× for 8→2) mid-transition (`styles.h:115`,`:148`); an effect sized for the coarser endpoint traps. Surface the peak requirement at the provisioning site.
58. [Low] **daydream — dead per-frame axis writes.** When `labelAxes` is on, `stepSimulation` sets static axis positions and pushes six fresh literal label objects every frame (`driver.js:488-492`). Drop the redundant writes.
59. [Low] **daydream — `disposeApp` is not idempotent and never removes the `pagehide` listener** (`daydream.js:809-811`); add a `disposed` guard.
60. [Low] **daydream — `statsGroup` caches element handles (including nulls) once** for the page lifetime (`driver.js:784-793`); brittle if the stats DOM is built after first paint.
61. [Low] **daydream workers — `setEffect`/`setParameter`/`setAnimationsPaused` don't rebuild a faulted pool** the way `setResolution` does (`segment_controller.js:568-603`); the action produces no render and the fault overlay stays frozen. Route them through the same rebuild-on-fault path or document the asymmetry.
62. [Low] **daydream workers — per-frame `Uint16Array` transfer-buffer allocation** (`segment_worker.js:191`,`:236`); unavoidable without a controller→worker return channel, but a per-worker per-frame GC source worth eliminating if worker churn shows in profiling.
63. [Low] **daydream workers — stats freeze during overruns** (`segment_controller.js:938-952`); only the `pendingFrame` branch calls `updateStats()`. Cosmetic.
64. [Low] **Targets — `getArenaMetrics` stack "usage" is sampled at the metrics call, not render peak** (`wasm.cpp:637`), misleading next to the cumulative arena usage. Only `high_water_mark` is meaningful; drop or relabel the field.
65. [Low] **Targets — `capture_screenshots.mjs` brittly scrapes daydream's `<option>` text and `?effect=` URL scheme** with no cross-repo version pin (`capture_screenshots.mjs:129`,`:182`); a daydream UI change silently degrades the gallery. Fails safe with a loud summary, but an asserted contract would harden it.

### Priority 5 — Documentation & comment discipline

66. ✅ **Cross-cutting — pervasive essay-length comment blocks.** Flagged independently in nearly every scope; several rationale blocks dwarf and obscure the code they annotate, running against the repo's own terse-comment convention. Trim to fact-focused comments.
67. ✅ **Engine — stale comment cites `__LINE__` but the macro pastes `__COUNTER__`** (`platform.h:852`); the `__COUNTER__` choice is deliberate (two uses on one source line must not collide), so the comment contradicts the rationale.
68. ✅ **Doc — IslamicStars registered params drifted from README §9.** Code registers `Fade, Face Fade Lo/Hi, Burst, Ripp Amp/Decay/Dur, Debug BB` (`IslamicStars.h:47-61`); README lists a nonexistent "Duration" and omits the Face-Fade sliders.
69. ✅ **Doc — README §7.10 advertises ID pins 21–23 and an N=8 path,** but the driver hard-caps `static_assert(N <= 4)` with only `PIN_ID0`/`PIN_ID1` (`pov_segmented.h:90`,`:102-106`). Align the datasheet to the shipped N≤4 reality.
70. ✅ **Hardware — single-board `pov_single.h` reads `effect_` in the column ISR as a plain static pointer** (`pov_single.h:194`), relying implicitly on `IntervalTimer` begin/end quiescence. Add the one-line ownership-contract comment its segmented sibling has, to deter a future mid-show mutation.
71. ✅ **Tests — `print_operand` prints "?" for Quaternion/Color4-shaped aggregates** (`test_harness.h:133-135`); a failing `HS_EXPECT_EQ` on a quaternion shows "? vs ?" in a math-heavy suite. Add a `w`-component/tuple-like branch.
72. ❌ **daydream tests — `@ts-check`/`@ts-nocheck` split on the DI-heavy files** (`segment_worker`, `segment_controller`, `engine_contract_wasm`) loses type-safety exactly where the mocked contract most needs to track the real one. Rejected: the repo has no `tsconfig`/`@types/node` and no CI typecheck, so enabling `@ts-check` on these files pulls the generated `holosphere_wasm.js` blob and node globals into checking (hundreds of unfixable errors) for editor-only benefit — disproportionate; the `@ts-nocheck` pragma is deliberate.
73. ❌ **daydream tests — `gui.test.js` window-teardown flake** is largely mitigated (mock timers, `afterEach` window restore) but worth keeping as a watch item, not treated as resolved. Rejected: no code change warranted — the finding is itself a watch note on already-mitigated behavior.

### Accepted / known gaps (documented; no action required — listed for completeness)

- **`daydream.js` entry orchestrator (811 lines) untested** — inherently DOM-wiring with no exports; the right response is the existing pattern (keep extracting pure logic into thin modules), not a giant DOM harness.
- **No JS-side golden-frame regression** — `engine_contract_wasm` asserts shapes, not pixels; frame content is covered on the C++/WASM side via the framebuffer-dump harness. Reasonable division of labor.
- **`GnomonicStars` `warp_speed` default (0.035) vs pinned spawn speed (0.0)** — `draw_frame` overwrites the seed before the first `timeline.step`, so the 0.0 is a genuine don't-care. Not a bug (the two literals merely look like they should agree).
- **eDMA register/ISR layer untested** — deliberately non-mockable; covered indirectly via the transport-templated `DMALEDController` with `MockStrip`.

---

## 6. Notable Engineering Strengths

Recorded for balance and to inform the grade rationale:

- **A verified fail-fast safety net.** `HS_CHECK` traps at the violation site on cold paths only, survives `NDEBUG`, pulls in no stdio, and is proven to actually fire by a shell-free SIGILL death harness with a minimum-case floor (41 invariants).
- **Drift-proof single-sources-of-truth.** `HS_EFFECT_LIST`, `HS_RESOLUTIONS`, `MESHOP_LIST`, the test-module roster, and the include block are each cross-checked in both directions (runtime registry vs X-macro, CMake vs runner, generator vs committed artifact) so coverage cannot silently drop.
- **Provably-safe type erasure.** `TimelineEvent` inline storage uses real `static_cast` + `std::launder`, an `alignof`/`sizeof` audit, and a build-wide `LARGEST_CONCRETE_ANIM_SIZE` check; borrow contracts are enforced by deleted rvalue overloads at compile time.
- **A genuinely elegant real-time protocol.** The one-wire, count-coded, odd-only distance-2 sync alphabet degrades a lost/spurious edge to a *missed* symbol, never a *misclassified* one, and derives column position from the cycle counter so masked-IRQ windows cannot drop columns.
- **Verified C++↔JS parity.** `color_parity_wasm`, `spline_math_wasm`, `engine_contract_wasm`, and `segment_crosscheck` pin JS ports against absolute goldens and the real WASM, so even a *coordinated* drift is caught; every intentional divergence is documented at the site.
- **Best-in-class WASM/GL memory hygiene on the JS side** — the zero-copy view-detach hazard is centralized, routed, and re-verified each frame, with symmetric Three.js resource disposal and complete event-listener teardown.

---

*End of review.*
