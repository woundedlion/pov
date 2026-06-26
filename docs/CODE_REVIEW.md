# Holosphere — Code Quality Review

A persistence-of-vision LED sphere engine (templated C++17, Teensy 4.x firmware) and its
WebAssembly + Three.js simulator (`daydream`). This review covers the rendering engine
(`core/`), effects (`effects/`), hardware drivers (`hardware/`), per-target entry points
(`targets/`), the test suite (`tests/`), build/CI tooling (`scripts/`, `CMakeLists.txt`,
`justfile`), and the daydream web client.

**Out of scope (per request):** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/rotate.h`. Vendored third-party code (`core/FastNoiseLite.h`, `daydream/three.js/`,
`node_modules/`, the generated `holosphere_wasm.{js,wasm}`) was not audited.

**Method.** The README was read for architectural intent; the codebase was partitioned across
fourteen examination sub-agents (one per component group), and every concrete defect claim was
then re-checked by independent validation sub-agents. Several initially-reported issues were
**refuted** during validation and are deliberately excluded (e.g. a suspected `palette_idx`
out-of-bounds — constrained by a shuffle-permutation invariant; a suspected `Feedback` flush
off-by-one — provably unreachable; a suspected `Plot::PlanarPolygon` projection-basis bug —
mathematically equivalent to its siblings).

---

## Executive Summary

This is an exceptionally mature, internally-consistent codebase. It reads like a long-lived,
heavily-reviewed engine rather than a hobby project: invariants are enforced where they live,
the fail-fast doctrine is applied with discipline (cold-path traps, hot-path stripped asserts
backed by setup-site traps), numerical singularities are handled deliberately and uniformly,
and the host/device/WASM split is bit-parity-conscious throughout. The single most material
finding is a **supply-chain** one (CDN scripts with no integrity pinning on the public deploy);
everything else is Low-severity — latent robustness gaps, documentation/accuracy nits,
localized performance opportunities, and a handful of weak test assertions. **No Critical or
High-severity correctness, memory-safety, or concurrency defects were found and confirmed.**

### Overall Grade: **A**

| Dimension | Grade | Rationale |
|---|---|---|
| Correctness | A | Edge cases (poles, antipodes, seam-wrap, first-frame state, fixed-point parity) are handled deliberately; confirmed defects are narrow and mostly latent. |
| Memory safety | A | Containers are exact-capacity-bound with `HS_CHECK` overflow traps; arena lifetimes are tracked with debug generation/rebind stamps; no confirmed UAF/OOB path. |
| Concurrency / ISR & real-time safety | A | Single-writer model is documented and enforced; release/acquire placed exactly where cross-context publication needs it; cycle-counter-derived column position is structurally drop-tolerant. |
| Numerical robustness | A | Pervasive radicand floors, `acos` domain clamps, denom epsilons, sign-preserving stereographic clamps, and a `-ffast-math`-guarded NaN→hi clamp contract. |
| API design & interface expressiveness | A | Compile-time filter-ordering `static_assert`s, type-erased borrow-vs-store callables with deleted rvalue overloads, designated-initializer param structs, X-macro single-sources-of-truth. |
| Architectural elegance | A | Clean SDF/scan split, one rasterizer for geodesic+planar edges, variadic compile-time filter pipeline with automatic coordinate-domain transitions, three-arena partition. |
| Readability & naming | A | Precise names; load-bearing invariants documented at the exact site a refactor would break them. |
| Documentation | A+ | Doxygen + rationale comments explain *why*, not just *what*; among the best-documented embedded code in this class. |
| Error handling (fail-fast doctrine) | A | Consistent cold-seam `HS_CHECK` traps, hot-path stripped `assert`s backed by setup-site traps, soft handling reserved for genuinely transient conditions — verified by a death harness. |
| Performance & resource efficiency | A | LUT-backed trig, per-face adaptive LUT sizing, coherence-grid culling, on-demand render in the client; a few uncapped per-pixel paths noted below. |
| Test coverage & quality | A− | Broad subsystem coverage with strong oracles (Euler invariant, brute-force KNN cross-check, fixed-point round-trips, a multi-board clock-modeled sync simulator) and a best-in-class death harness; a minority of assertions are count-only or tautological. |
| Build system & reproducibility | A | CMake encodes hard-won rationale (stack-size ordering, LTO fast-math attribute flow, SHA provenance with dirty-tree guard); honest minimum-version pin. |
| Portability (host/device/WASM parity) | A | Arduino-gated device code, host-testable protocol/index cores, FastLED integer-semantics mocks preserving sim/device bit-parity. |
| Web client quality | A− | Exemplary WASM memory-view detach handling and Three.js disposal discipline; minor URL-state and recorder-labeling edges. |
| Security / supply-chain | B | XSS-safe overlays and JS-boundary input validation are strong, but the committed deploy pulls executable ESM from a CDN with no Subresource Integrity. |
| Maintainability | A− | Strong factoring overall; a few deliberately-diverged sibling families and duplicated boilerplate carry drift risk. |

### Per-Component Grades

| Component | Grade | Notes |
|---|---|---|
| `core/` math, color, geometry | A | Singularity discipline and clamp-before-cast are systematic; OKLCH envelope caching is strong perf-vs-quality engineering. |
| `core/` pipeline, filters, type erasure | A | Compile-time pipeline safety and the borrow/store callable split are textbook. |
| `core/` rasterizers (sdf/scan/plot/anim) | A− | Centralized seam handling; minor doc/perpetual-counter caveats. |
| `core/` mesh, conway, hankin, spatial | A | Exact capacity sizing with trap-on-overflow; documented scratch-polarity contract. |
| `core/` memory & effect registry | A | Wrap-proof bounds arithmetic; Meyers-singleton avoids static-init-order hazards. |
| `effects/` (27 effects) | A− | Careful lifecycle/singularity handling; localized perf and duplication items. |
| `hardware/` drivers | A | Among the best-justified embedded concurrency code reviewed; no confirmed races. |
| `targets/` + build + scripts | A | Single-source-of-truth marshaling; provenance-tracked WASM install. |
| `tests/` | A− | Broad and doctrine-aligned; tighten a handful of weak assertions. |
| `daydream/` client + tools | A− | Solid bridge/disposal/worker-fence design; supply-chain and minor UX edges. |

---

## Prioritized Items to Fix

Items are numbered sequentially and grouped by priority. Severity reflects real-world impact on
this project (a hobbyist/art LED installation plus a public read-only web demo), not abstract
severity.

### Critical

_None found._

### High

_None found._

### Medium

1. `daydream/vendor-importmap.js:48,56-57` — The committed default importmap loads `three.js` and `lil-gui` from `cdn.jsdelivr.net` with **no Subresource Integrity** (import maps cannot carry an `integrity` attribute, and no integrity-map fallback is provided). A CDN compromise or MITM would inject arbitrary executable ESM into the live GitHub Pages deploy. Version pinning (`@0.183.1`, `@0.21.0`) reduces but does not remove the mutable-tag risk. Mitigate by self-hosting the vendored copies on the deploy (the `npm run importmap:local` path already exists), or document the accepted risk explicitly. Blast radius is limited (no auth/secrets/user data — a defaced canvas), hence Medium not High.

### Low — Correctness & Robustness

2. ✅ `effects/MobiusGrid.h:243` — The longitude basis uses a bare `cross(v, w).normalized()` that yields NaN if the line normal becomes parallel to `Y_AXIS`. Unreachable today (`theta ∈ [0,π)` keeps the normal in the XZ-plane), but sibling code in the same file uses `normalized_or(...)`; add the same fallback for robustness against future parameter changes.
3. ✅ `effects/Dynamo.h:354` — `drag()` runs `while (shortest_distance(...) > gap) move(follower);` with no iteration guard. Safe under the current Gap slider max (20) vs `W ≥ 96`, but it relies on an unasserted `gap < W/2` invariant; add an iteration cap or a cold `HS_CHECK` on the gap range.
4. ✅ `effects/Dynamo.h:228,241` — `color()` reads `baked_palettes_[i + 1]` with no guard; correctness depends entirely on the push-paired invariant in `color_wipe()`/`reap_completed_wipes()`. In-bounds today, but a cold assert at the band-walk entry would harden the seam against future edits.
5. ✅ `effects/MindSplatter.h:101-102` — `draw_frame` overwrites every attractor's `strength` with the scalar `params.well_strength` each frame, making the per-attractor `strength` argument to `add_attractor` dead/misleading. Either honor per-attractor strength or remove the parameter.
6. ✅ `core/geometry.h:486` — `vector_to_pixel` projects via `fast_atan2`/`fast_acos` (approximate), while the inverse `pixel_to_vector` uses exact trig, so `vector→pixel→vector` is sub-pixel asymmetric. Document the projection tolerance, and correct the README wording (§6) that implies an *exact* `theta/2π·W` projection.
7. ✅ `targets/wasm/wasm.cpp:308-327` — `setResolution` tears the effect down without calling `stack_paint_canary()` (unlike `setEffect`), so a `setResolution` not followed by `setEffect` leaves a stale stack high-water mark. Metrics-accuracy only — no memory-safety impact.
8. ✅ `daydream/segment_controller.js:505-524` — The boundary overlay (`showBoundaries`) is stamped even when `blitted === 0` (a fully generation-fenced frame), drawing cyan seam lines on an otherwise-black sphere that is still shown on-screen during a resolution change. Gate the overlay on `blitted > 0`.
9. ✅ `daydream/recorder.js:280` — `_extension()` returns `mp4` only for `video/mp4`, else unconditionally `webm`; an exotic MediaRecorder fallback mimetype (e.g. `video/x-matroska`, or empty) produces a `.webm` file/Blob type that disagrees with the real container. Derive the extension from the actual mimetype.
10. ✅ `daydream/gui.js:236` — An out-of-range **enum** value hydrated from the URL sets `urlApplied = false` but does not rewrite the URL, whereas numeric clamping does; a shared link keeps the stale enum value in the address bar. Make the enum path symmetric.
11. ✅ `daydream/state.js:113` — `_notify` iterates the live `_listeners` array; a subscriber that calls `subscribe()` during a notification would be invoked for the current event (add-during-notify is not snapshotted). No current caller does this, but snapshot the array or document the contract.
12. ✅ `daydream/daydream.js:116` — `EngineHost.refresh()` re-points `instanceColor.array` without setting `needsUpdate`; it works only because both current callers flag it afterward. Document the "caller must flag `needsUpdate`" contract or set it inside `refresh()`.

### Low — Documentation & Consistency

13. ✅ `core/util.h:27,54` — The `wrap(T,U)` doc says `m <= 0` yields "NaN/garbage", but the `fmax(m, min())` floor actually produces a silent tiny-positive modulus (≈0) under `NDEBUG`. Fix the comment to match behavior (the debug `assert(m > 0)` already covers dev builds).
14. ✅ `effects/GSReactionDiffusion.h:210-211` — The comment claims `interpolate_b`'s weight guard prevents a 0/0 NaN from reaching `palette.get()`, but the actual NaN cull (`b < B_CULL_THRESHOLD`) lives in `render()`. The guard is correct (returns 0 before dividing); only the comment misattributes the filter location.
15. ✅ `core/effects.h:48` — `HS_EFFECT_LIST` (X-macro) and the `#include` list are two hand-maintained lists; a forgotten `X()` row only *silently* drops native smoke coverage (the anti-drift `static_assert` exists only on the WASM `HS_EFFECT_COUNT` path). Add a registry-size check to the native test harness.
16. ✅ `effects/SphericalHarmonics.h:103,232,320` — Magic literals `6` (seed mode) and `24` (`rand_int(1,24)`) implicitly couple to "modes through l=4" with no named constant or `static_assert`; widening one without the other silently changes the visual range.
17. ✅ `daydream/tools/mobius.html:343-348` — GLSL `cmult`/`cadd`/`cdiv` (including the `1e-6` divide guard) duplicate `mobius_transforms.js`; the on-screen warp runs the GLSL copy while presets use the JS copy, so an edit to one silently diverges. Add a cross-reference comment (GLSL cannot import the JS module).
18. ✅ `daydream/tools/clipboard.js:67` — `copyWithFeedback`/`copyToClipboard` return a success boolean that no caller inspects (`solids.html:1184`, `splines.html:494`, `lissajous.html:370`), so a copy failure is fully silent — the "Copied!" label never flips and nothing surfaces. Surface a failure state.

### Low — Performance

19. ✅ `effects/Moire.h:147` — `draw_layer` rasterizes `ceil(params.density)` rings per layer with no per-row budget cap (Density slider max 100, ×2 layers), unlike `ShapeShifter` which caps at `H`. On the 96×20 target many rings collapse below a row and waste rasterizer work; add a resolution-aware cap.
20. ✅ `effects/PetalFlow.h:221` — `draw_ring` recomputes `expf`/`fast_cosf`/`fast_sinf` per sample for every active ring every frame (~9,200 `expf`/frame at `W=288`, up to `MAX_RINGS=64`). Ring geometry depends only on `rho`/`twist`/static shift; cache per-rho buckets for the device budget.
21. ✅ `effects/SphericalHarmonics.h:86-90,297` — Each pixel evaluates `associatedLegendre` up to twice over the full sphere with no horizontal-interval narrowing and full-height vertical bounds; acceptable for `l ≤ 4` but it is the group's heaviest per-pixel cost and is uncapped. Cap `l` or document the cost.

### Low — Test Quality

22. ✅ `tests/test_styles.h:318` — `test_sync_noise_pushes_scalars` asserts `HS_EXPECT_TRUE(true)` on the unbound (`noise == nullptr`) branch — a tautology that cannot fail and does not prove the null guard skipped the dereference. Assert observable state instead.
23. ✅ `tests/test_animation.h:1299` — The collapsed-frame `deep_tween` test asserts `gts.size() >= 1` (always true) rather than the documented exact count `M + (N-1)*(M-1)` (which the sibling non-collapsed test at `:1269` does pin), so the boundary-skip logic is not actually verified.
24. ✅ `tests/test_mesh_raster.h:341` — The solid-fill coverage threshold is `> 99%` for a closed convex octahedron, tolerating ~1% dark pixels that could hide edge holes or clipping artifacts; a closed convex solid should fill far tighter.
25. ✅ `tests/test_sdf.h:1007` — The interior-cull check uses an arbitrary loose aggregate bar (`total_interior > 1000` summed over 84 cases); a regression dropping hundreds of interior pixels still passes. Assert per-shape coverage.
26. ✅ `tests/test_conway.h:339-378` — The `dual`/`kis` cube tests assert vertex/face **counts** (plus all-triangles for `kis`); the suite's `check_euler_genus0` manifold oracle covers much of the structural risk, but there is no per-vertex degree or face-type-histogram oracle, so a count-correct topology error could slip through.
27. `tests/test_effects.h:67,71` — The default `kDefaultFrames = 8` is below several effects' lifecycle windows (RingShower slot reuse, Thrusters FIFO, ShapeShifter 48-frame cut, DreamBalls 288-frame respawn); those paths run only when CI raises `HS_SMOKE_FRAMES`, so a local pre-commit run is weaker than it appears.

### Low / Trivial — Maintainability & Nice-to-Have

28. ✅ `effects/BZReactionDiffusion.h` / `effects/GSReactionDiffusion.h` — The odd-parity ping-pong land-back (BZ's `advance_substeps()` vs GS's inlined copy at `:230-243`) and the cubemap-LUT vertex/fragment-shader preamble are near-duplicated; hoist both into `ReactionDiffusionBase` to prevent drift.
29. ✅ `effects/FlowField.h:46` — `rand_int(0, 65536)` uses a different upper-bound convention than sibling `ShapeShifter.h:99` (`rand_int(0, 65535)`); align the half-open-vs-inclusive convention across call sites.
30. ✅ `effects/DistortedRing.h:96-98` — `shiftFn` is rebuilt as a fresh `this`-capturing lambda inside the per-ring loop every frame; hoist it out of the loop.
31. ✅ `justfile` (`screenshots` recipe) — Shells through `npm run screenshots` while every sibling recipe (e.g. `smoke`) calls `node scripts/...` directly, adding an undocumented npm layer; align for consistency.
32. ✅ `daydream/tools/solid_codegen.js:69,120-123` — Function-name parameter suffixes quantize (0.01, and 1° for hankin), so two solids whose params round equal collide on one `funcName` and the later C++ paste overwrites the earlier; `snub` additionally ignores `t`/`twist` params. Documented tradeoffs, but worth surfacing a collision warning.
33. ✅ `daydream/driver.js:166` — The aspect guard `this.canvas.width / this.canvas.height || 1` yields `Infinity` (truthy, bypassing the `|| 1`) when `height === 0`; transient and corrected by `setCanvasSize()`, but the guard doesn't cover the case it appears to intend.
34. ✅ `daydream/geometry.js:35` — `pixelToSpherical` hard-codes the `H_OFFSET == 0` simulator assumption (acknowledged in-comment); it is an unenforced cross-repo invariant that a future device-parity mode would silently violate.
35. ✅ `core/animation.h:350-355` — Perpetual animations (`duration == -1`) increment a signed `int` frame counter forever; `t++` is signed-overflow UB past ~2³¹ frames. Documented as an accepted limit ("restart before then"); consider an unsigned counter for defensiveness.
36. ✅ `hardware/pov_single.h:116` — `duration * 1000` (`unsigned long`) overflows only past ~49.7 days; unreachable in practice, but add a guard for symmetry with the segmented driver's `run_show`.
37. ✅ `hardware/pov_single.h:171,180` — The single-board DMA path discards `submitFrame()`'s result, giving a weaker overrun watchdog than the segmented driver. Acceptable given the ~1.3 ms single-board column period (an overrun is realistically impossible), but the asymmetry is undocumented here.
38. `tests/test_memory.h:209,572` — The ArenaVector stale-binding and generation-bump (use-after-reset) contract tests are `#ifndef NDEBUG`-gated. They *do* run under the canonical Debug `tests` preset, so this is not a live gap; note it only because a non-canonical Release/NDEBUG test build would skip the only automated cover for those contracts.

---

## Notable Strengths

These materially raise the grades and are worth preserving through future changes:

- **Fail-fast placed by doctrine, and *verified*.** Cold bind/alloc/config seams trap via always-on `HS_CHECK`; hot per-pixel paths use stripped `assert`s backed by a cold trap at the corresponding setup site. An out-of-process death harness asserts the traps actually fire (`SIGILL`/`STATUS_ILLEGAL_INSTRUCTION`) with a shape-probe sentinel — the safety net is tested, not assumed.
- **Compile-time pipeline safety.** `filter.h`'s `static_assert` chain (terminal-last, no-2D-before-3D, history-flush signatures) turns whole classes of filter-composition mistakes into build errors; the `FunctionRef` (borrow, accepts rvalue) vs `StoredFunctionRef` (`= delete`d rvalue) split encodes the dangling-store hazard directly in the type system.
- **Singularity discipline is systematic.** Every forward stereographic/gnomonic/Möbius projection emits a shared `STEREO_INF` sentinel that every inverse recognizes with margin; radicand floors, `acos` domain clamps, and denom epsilons are applied uniformly, and a `__FINITE_MATH_ONLY__` `#error` protects the load-bearing NaN→hi clamp contract against a `-ffast-math` regression.
- **Real-time hardware design.** Column position is derived from the CPU cycle counter (`now − epoch`) with per-half-rev rebase, making masked-IRQ windows structurally unable to drop columns or observe 32-bit wrap; the single-writer concurrency model fuses test-and-take (`EdgeMailbox::try_claim`) so a split-race cannot be reintroduced, and memory-order choices are individually justified rather than cargo-culted.
- **Single-source-of-truth everywhere.** `HS_WASM_RESOLUTIONS`, `MESHOP_LIST`, `HS_EFFECT_LIST`, and `param_marshal.h`'s def/value ordering each drive multiple consumers from one declaration; the geometry tools route real geometry through the actual WASM exports rather than re-implementing it.
- **WASM bridge & client hygiene.** The typed-memory-view detach-on-heap-growth hazard is factored into one tested module with re-fetch guards at every consumer; Three.js disposal is symmetric and complete (geometry/material/controls/observers/listeners/workers), and the WASM-aliased `instanceColor.array` is nulled before `dispose()` so the GPU cannot re-upload freed memory. The segmented-render generation fence is re-checked at four independent layers.
- **Documentation as engineering.** Non-obvious invariants (scratch-arena polarity, the double-buffer TOCTOU reasoning, find_nearest_node chaining, OKLCH envelope caching) are explained precisely where a future refactor would break them — the comments routinely capture genuine hazards rather than restating code.
