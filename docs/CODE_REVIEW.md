# Holosphere + Daydream — Code Quality Review

**Date:** 2026-06-06
**Reviewer:** Expert C++/TypeScript code audit
**Method:** 16 parallel component agents, each reading the assigned source directly and citing `file:line`, judged against the project's own stated philosophy (fail-fast traps, 16-bit linear color, compile-time `<W,H>`, single partitioned arena). Findings were cross-checked against the prior review; the highest-impact new findings (the shipped `-ffast-math` guard loss, the `Lerp` dangling-pointer, `Test`/`TestShapes` in the registry, the absent `deploy.yml`) were hand-verified by the lead reviewer.
**Scope:** ~47k lines: the hand-written C++ engine (`core/`), 27 effects, 3 hardware drivers + `platform.h`, the WASM/embind bridge, a 23-module native test suite, the CMake/CI build; plus the sibling vanilla-JS `daydream` web simulator. Generated tables (`reaction_graph.cpp`, `color_luts.h`) and third-party code (`FastNoiseLite.h`) were judged on integration/generation correctness, not line-by-line.
**Out of scope (per review rules):** `effects_legacy.h` and `targets/Holosphere/Holosphere.ino`.

---

## Overall grade: **B+**

This remains an exceptional solo engineering project — in most places it reads like production firmware from a top embedded studio. The architecture is principled and, more rarely, *consistently applied* across a very large surface: every render class is templated on `<W,H>` for zero-overhead specialization; all color is 16-bit linear light with OKLCH-correct interpolation; memory is a single partitioned arena with explicit `(Arena&, Arena&)` plumbing and no hidden state; the rendering core is a variadic `Pipeline<W,H,Filters...>` that performs automatic 2D↔3D coordinate coercion at compile time; and the fail-fast doctrine (`HS_CHECK` traps surviving `NDEBUG`) is genuinely wired into the memory and lookup layers. The performance engineering is expert-level — split trig LUTs, anisotropic per-face distance LUTs with correctness-preserving validity gates, fast-reject culling, a non-blocking DMA path with a single cache flush, zero-copy WASM readback. The 1,857-line README is outstanding documentation.

The prior review's headline fixes **hold up under independent re-audit**: the `H_OFFSET` double-apply is fixed and now regression-tested; `SDF::Subtract`/`Intersection` sort their child intervals and overflow now traps through `push_interval`; the silent-fallback *lookup* cluster (solids, spatial-hash, registry, param-list, palette modifiers) is converted to traps or `static_assert`s; the `relax` squared-length unit bug is fixed; `geometry.js` no longer aliases mutable singletons; the segment lifecycle is cleanly extracted into `SegmentController`; the `plot.h` vertex-shader loop is genuinely deduplicated; and the snake_case migration is complete and consistent.

What keeps it at **B+ rather than A−** is that this fresh, independent pass surfaced more *open* correctness items than the prior self-review implied remained — and they cluster into a coherent, fixable theme rather than being scattered noise:

1. **The fail-fast doctrine stops at the math/geometry core.** `3dmath.h` and `geometry.h` contain *zero* `HS_CHECK`. Degenerate cases log-and-substitute a plausible-wrong value (zero-vector `normalize` → `(1,0,0)`, `angle_between` → `0`), and the only invariant guards present (`unit_inverse`, `operator/=`) use `assert()`, which `NDEBUG` strips. This is the exact "mask the violation with plausible-but-wrong data" anti-pattern the project's own memory calls out — now the largest remaining instance of it.

2. **The shipped WASM artifact compiles its own safety nets out.** The live demo is built `-O3 -ffast-math` (→ `-ffinite-math-only`, which folds `std::isfinite()` to constant `true`) with the default `-DNDEBUG` (which strips every `assert`). The pipeline's documented NaN/bounds guards (`filter.h:92,106`) are therefore *both* dead in the only build users run, while the test suite (Debug, finite-math on) sees neither failure mode.

3. **A handful of genuine latent bugs at the hard seams** — pole/wrap/antipode/planar-projection edges — plus one borrow-safety hole (`Lerp`) inconsistent with the rest of the codebase's discipline.

4. **The live-art interface keeps regrowing inert sliders.** The prior cluster was fixed, but this pass found a fresh one (Flyby `Falloff` dead, SplineFlow `Drift` and Comets `Cycle Dur` read-once, ChaoticStrings' stub timer, Liquid2D's `lerp` silently dropping a field, GS `dt`'s default above its own slider max). This is a recurring hygiene gap, not a one-time cleanup.

None of these is fatal; all are fixable; most are localized. The bones are excellent and the *core* (math LUTs aside, color, pipeline, memory, mesh, WASM bridge, hardware concurrency) is genuinely shippable.

---

## Dimension grades (aggregated across all components)

| Dimension | Grade | Justification |
|---|---|---|
| **Architecture & Design** | **A−** | The variadic `Pipeline` with transparent domain coercion, the Conway-operator algebra + fluent `SolidBuilder`, the "animations mutate state / effects read state and draw" separation, the single-`CMakeLists` dual-mode build, the X-macro resolution single-source-of-truth, and the `SegmentController` DI are all elegant and *orthogonal*. The convention discipline across 60+ files is the standout. Docked by the two brute-force effects, IslamicStars leaking arena/carousel lifecycle into the effect body, and `Test`/`TestShapes` shipping in the registry. |
| **Performance & Efficiency** | **A−** | The strongest dimension. Microsecond-budget awareness everywhere: split sin/cos LUTs (~145× RAM win), anisotropic distance LUTs with exact fallback, fused linear-correction→sRGB wire packing with no intermediate buffer, branchless segmented ISR, zero-alloc steady-state loops, single-draw-call instanced web mesh, provably-correct precomputed reaction graph. Dragged only by Metaballs/Voronoi (`O(W·H·N)` per-pixel loops bypassing the rasterizer) and a per-call `initialized` branch on the trig LUTs. |
| **Correctness & Robustness** | **B−** | The genuine limiter. The math core ships silent degenerate fallbacks and `NDEBUG`-stripped asserts; the shipped build strips its NaN guards; `Lerp` can dangle; `Canvas` construction is a TOCTOU on relaxed atomics; `scan_region`'s coalescer isn't seam-wrap-aware; `Face` truncates >MAX_VERTS to wrong geometry; non-manifold edge fans and >32767 vertex indices escape the fail-fast net; `plot.h` planar mode feeds a chord length into arc-length thresholds; `setClip` forwards unvalidated JS ints into `%` arithmetic. Most are latent at current params — but they are precisely the regressions the test suite cannot catch. |
| **Readability & Maintainability** | **A− / B+** | Outstanding Doxygen and *why*-not-*what* comments at every subtle seam (relaxed-atomic rationale, buffer-detach contract, `Persist` member ordering, H_OFFSET warnings). Held back by recurring copy-paste (`DistortedRing`≈`Ring`, BZ/GS physics the CRTP base didn't capture, the slew-rate IOMUX block ×3) and profiling cruft shipped in `IslamicStars::draw_frame`. |
| **Memory Safety & Resource Mgmt** | **B+** | Disciplined arena ownership; `Arena::allocate`/`ArenaVector` are model fail-fast citizens, making raw `allocate()` derefs safe by construction; `Persist`/`ScratchScope` RAII and Three.js GPU disposal are correct. Undercut by `ArenaSpan` lacking the generation tracking its sibling has, the hardware `END_FRAME_BYTES` odd-N spec mismatch + in-ISR busy-wait on DMA overrun, and the test-only 8 MB arena widening hiding MindSplatter's known on-device over-subscription. |
| **API / Interface Expressiveness** | **B+** | `Pipeline<W,H,...>`, `SolidBuilder`, `REGISTER_EFFECT`, the compile-time `=delete` rvalue guards on non-owning animations (where present), the X-macro resolution SSOT, and `AppState`/`URLSync` are a genuine pleasure and largely misuse-resistant. Pulled down by the recurring inert/read-once effect params, the `Lerp` raw-pointer lifetime hole, and the untyped `setClip(y0,y1,x0,x1)` crossing the JS boundary with no validation. |
| **Testing** | **B+** | Excellent native harness asserting *real* algebraic/topological invariants (Hamilton algebra, slerp composition, Euler V−E+F=2, half-edge pair symmetry, Conway winding oracle, OKLab round-trips, machine-verified reaction-graph integrity). Held back because the two things the project markets most — the always-on trap surface and the 335 KB device budget — are exactly the two untested: no death test ever fires a trap, the smoke harness runs on an 8 MB arena so it *cannot* catch MindSplatter's over-subscription, the terminal smoke assertion (`acc+1>0`) is tautological, and `fake_canvas()` is UB. |
| **Documentation** | **A−** | The README is a genuinely outstanding engineering document. Downgraded from the prior A because it has drifted: it describes a `deploy.yml` gated-deploy CI layer that does not exist in the repo, `VectorTrail` as int16/uint8-quantized when it stores full-precision `Vector`, and "DMAMEM placement" for DMA buffers that aren't tagged `DMAMEM`. |

### Per-component snapshot

| Component | Arch | Correct | Perf | Mem | API | Overall | Standout |
|---|---|---|---|---|---|---|---|
| Core Math & Geometry | B+ | C+ | A− | B | B+ | **B** | `TrigLUT<W,H>` split-trig reconstruction (~145× RAM win) |
| Color & Palettes | A− | C+ | A | B+ | A | **A−** | OKLCH shortest-arc interp with achromatic guards |
| SDF + Scan Rasterizer | A− | B | A | A− | A | **A−** | Anisotropic distance LUT with sign+margin validity gate |
| Plot (Curve) Rasterizer | A− | C+ | A | A | A | **A−** | Simulate-then-replay `sin(φ)` adaptive stepping |
| Pipeline / Filters / Canvas | **A** | C+ | A− | A− | A | **A−** | EBO variadic pipeline + `if constexpr` domain dispatch |
| Animation / Transformers | **A** | C+ | A− | B− | A | **A−** | Compile-time `=delete` rvalue borrow guards (Motion/MeshMorph) |
| Memory / Arena | **A** | A− | A | A− | A− | **A−** | `Persist<T>` RAII compaction; verified-correct alignment math |
| Mesh / Conway / Hankin / Spatial / Solids | **A** | B | A− | A− | A | **A−** | `SolidBuilder` ping-pong-arena fluent DSL; DRY half-edge pairing |
| Effects A (RD / Hopf / Islamic / SH) | B− | B | B+ | — | C | **B** | BZ cubemap-LUT + 1-ring nearest-node walk; correct SH recurrence |
| Effects B (Mesh / Field / Particle) | B | A− | C+ | — | A− | **B+** | DreamBalls "declare state, let engine drive"; live re-seed gating |
| Effects C (Shader / Raymarch / Curve) | A− | B+ | A− | — | B | **A−** | Lipschitz-safe warped-SDF sphere tracing in Raymarch |
| Hardware Drivers + Platform | **A** | C+ | A | C | B+ | **B+** | Honest relaxed-atomic + explicit-cache-flush handoff |
| WASM Bridge + Reaction Graph | **A** | A− | A | B+ | A− | **A−** | X-macro resolution SSOT; provably-correct K-NN table |
| Tests | — | — | — | — | — | **B+** | Operator-agnostic Conway winding oracle |
| Build / CI / Cross-cutting | A− | C+ | — | — | — | **B** | Self-bootstrapping `toolchain-native-clang.cmake` |
| Daydream JS | A | A− | A | — | B | **A−** | Lazy-getter DI at the reassignable WASM seam |

---

## Cross-cutting themes

### 1. Fail-fast is enforced at the edges but not in the numerical core

The prior review's biggest finding — silent bounded fallbacks at *lookup/registration/config* — is genuinely closed and verified (solids `get_entry`/`get_by_name` trap, `SpatialHash::insert` traps, `registerParam` traps, `AnimatedPalette::add` traps, `get_fill_fn` is a compile error, `push_interval` traps). But the same anti-pattern survives in the place hardest to test and most numerically load-bearing:

| Location | Silent behavior | Should |
|---|---|---|
| `3dmath.h:192-212,392-411` | zero-vector `normalize`/`normalized` → `(1,0,0)`; quaternion → identity | `HS_CHECK(m > eps)` |
| `3dmath.h:718-726` | `angle_between` degenerate → `0` | trap or document |
| `3dmath.h:170,363` | `operator/=`, `unit_inverse` guarded by `assert()` (NDEBUG-stripped) | `HS_CHECK` |
| `geometry.h:186-187,252-254` | `PhiLUT`/`TrigLUT`/`Orientation` indexed by caller `y`/`x`/`i`, no bounds trap | `HS_CHECK` index |
| `sdf.h:1039-1042` | `Face` clamps `count` to MAX_VERTS, logs, builds *wrong* polygon | trap |
| `mesh.h:122-130` | >2 half-edges sharing an edge → arbitrary pairing, no trap | trap |
| `conway.h:423-428,552-554` etc. | `int16_t` vertex-index casts overflow silently past 32767 | trap or widen |

This is the single highest-leverage cluster remaining: mechanical to fix, and it converts a class of silent-wrong-geometry bugs into bench-time crashes — the entire point of the stated philosophy.

### 2. The shipped artifact disables the guards the engine relies on

`CMakeLists.txt:51` applies `-ffast-math` to the WASM **Release** build, implying `-ffinite-math-only` (so `std::isfinite()` is constant-folded to `true`); the default Release `-DNDEBUG` strips `assert`. The pipeline's NaN/coordinate guards (`filter.h:92,106`) are therefore *both* inert in the live demo, while the Debug test build sees neither. A NaN coordinate the tests would trap on becomes silent OOB UB in production. The one-line mitigation is to append `-fno-finite-math-only` after `-ffast-math` (keeps SIMD/reassoc, restores `isfinite`); the durable fix is to promote the hot finite-check to an always-on `HS_CHECK`. This theme *compounds* theme 1: the few guards that exist get compiled out exactly where it matters.

### 3. Latent correctness at the hard seams

The math core's pole/wrap/antipode handling is where the remaining real bugs live: `mod_tau` only wraps one period (`rotate.h:18-24`); `wrap()`'s `abs(x)<TOLERANCE` clamp maps small negative angles to the wrong column (`util.h:25-35`); `plot.h`'s planar strategy feeds a *projected chord* into arc-length-calibrated thresholds; `Line::sample` takes an unguarded `cross()` for antipodal points; `scan_region`'s coalescer isn't seam-wrap-aware for CSG output that crosses θ=0. The `Canvas` constructor is a check-then-act on relaxed atomics whose correctness rests on an external single-writer/single-ISR invariant rather than the memory model. None bite at current parameters — which is exactly why they need tests, not just fixes.

### 4. The live-art interface keeps accreting inert controls

`registerParam` correctly binds a live pointer; a slider is "live" iff `draw_frame` actually *reads* that member each frame. This pass found a fresh batch that don't: Flyby `Falloff` (registered + lerped, never read — `Flyby.h:23,66,121`), SplineFlow `Drift` (read once at init — `SplineFlow.h:34,46`), Comets `Cycle Dur` (read once into the animations — `Comets.h:35,42,44`), ChaoticStrings' empty `next_preset()` stub wired to a live `PeriodicTimer` (`ChaoticStrings.h:117-119`), Liquid2D's `Params::lerp` hardcoding `N=6` for a 7-field struct so `cycle_speed` is never tweened (`Liquid2D.h:146`), and GS `dt` defaulting to `2.5` above its own registered max `2.0` (`GSReactionDiffusion.h:48,178`). The fix is the same re-seed/re-read discipline already used correctly by Metaballs/Voronoi/DreamBalls; the deeper fix is a debug assertion (or a GUI "unused param" lint) so a registered-but-unread param is caught at build time.

### 5. Documentation and CI have drifted from the repo

The README §11 describes a 3-layer test gate culminating in a `deploy.yml` gated deploy to GitHub Pages; `.github/workflows/` contains only `ci.yml` — there is no `deploy.yml`, no Pages publish step, no secret usage. The `wasm-release-install` step is a local copy into `../daydream`. So the "gated deploy" layer is aspirational, and nothing in this repo gates or performs the publish. The README also still describes int16/uint8-quantized trails (the code stores full-precision `Vector`) and DMAMEM-placed DMA buffers (they're 32-byte-aligned static members, not `DMAMEM`). Documentation this good is worth keeping honest.

---

## Standout engineering

Excellence should be named, not just defects:

- **The variadic filter pipeline** (`filter.h:129-191`). `Pipeline<W,H,Head,Tail...>` derives from `Head` (EBO, zero storage) and dispatches `plot()` via `if constexpr (Head::is_2d)`, transparently lifting 2D↔3D coordinates at the mismatch boundaries. A filter author writes in whatever domain is natural and the chain stitches mismatched domains together at compile time with no runtime cost. The cleanest expression of the engine's compile-time-resolution thesis.
- **The anisotropic distance LUT with a validity gate** (`sdf.h:1460-1490`). Bilinear interpolation is trusted only when all four corners share sign *and* the smallest magnitude exceeds the per-face cell-diagonal safe distance; otherwise it falls back to exact register-packed edge evaluation. Textbook-correct LUT speed without interpolation artifacts near edges.
- **The honest relaxed-atomic rationale** (`dma_led.h:294-317`). It explicitly separates the single-observer ISR/thread flag (no barrier needed on a single core) from DMA buffer coherence (handled by `arm_dcache_flush_delete` + DSB) and warns against a bogus acquire/release "fix". Staff-level reasoning that most embedded code gets wrong — and it is correct; do not "fix" it.
- **OKLCH interpolation with achromatic guards** (`color.h:469-486`). Correctly special-cases near-zero chroma on either endpoint before the ±π hue wrap, avoiding the classic hue-spin artifact when interpolating to/from gray. Perceptual interpolation done right, and the known 8-bit regression is *actively defended* by `srgb_to_linear_interp` (`color.h:214-223,812-826`), not merely avoided.
- **Lipschitz-safe warped-SDF sphere tracing** (`sdf.h:2204-2207,2267` + `scan.h:762`). Raymarch's twisted-torus SDF supplies an analytical Lipschitz bound and the marcher applies a 0.9 safety factor — correct sphere tracing of a *warped* field, which naive marchers overshoot.
- **The X-macro resolution single-source-of-truth** (`wasm.cpp:90-98,166-172,325-328`). `setResolution`, `setEffect`, and `getEffectSizes` all dispatch through one `HS_WASM_RESOLUTIONS` list, so the supported set provably cannot drift between entry points — and the reaction-graph table is machine-verifiable (7680 rows, uniform K=6, all indices in range, no self-refs).
- **`SolidBuilder`** (`solids.h:182-268`). Each Conway op runs `MeshOps::op(mesh, *a, *b)` then `std::swap(a,b)`, hiding a two-arena ping-pong behind a fluent chain so Islamic-pattern recipes read as a declarative geometry DSL while staying single-allocation and arena-deterministic.
- **Lazy-getter DI at the reassignable WASM seam** (`segment_controller.js:25-31`). Because `ALLOW_MEMORY_GROWTH` reassigns the engine/view, passing `getWasmEngine()`/`getMemoryView()` instead of snapshots is the correct design and eliminates a whole class of stale-reference bugs; the detached-buffer guard (`daydream.js:88`) and re-fetch ordering are spec-correct.

---

## Prioritized fix list

### P0 — Shipped-artifact / real-time correctness

1. **✅ FIXED (2026-06-06) — `-ffast-math` made the live WASM build assume no NaN/Inf module-wide** (`CMakeLists.txt:51`; sink at `filter.h:84-115`). Validation refined the finding: the sink's `isfinite`/bounds checks are `assert`s, *deliberately* debug-only (documented at `filter.h:89-91,104-105`) since this is the per-pixel hot path where `HS_CHECK` is explicitly disallowed (`platform.h:25`) — and WASM Release already strips them via `NDEBUG` independent of `-ffast-math`. The real exposure was narrower: `-ffast-math` implies `-ffinite-math-only`, which folds `std::isfinite()` to `true` and lets the compiler assume finiteness **module-wide**, so a NaN an effect might emit could UB-cast at `filter.h:107` and bypass even the runtime clip/wrap guards. **Fix:** append `-fno-finite-math-only` after `-ffast-math` in the WASM Release compile options, restoring real finite-math semantics across the module while keeping the larger `-ffast-math` wins (reassoc/FMA/reciprocal/`-fno-math-errno`) for SIMD throughput. The flag is a compile-time float-semantics attribute that bakes into the LTO bitcode, so it lives on the compile line. **Perf:** a module-wide cost was accepted up front by the project owner (the `-ffinite-math-only`-specific optimizations are the smallest part of `-ffast-math`). The hot-path asserts are intentionally left debug-only per the documented design. *(Remaining follow-up, tracked under Testing: a CI check that a known NaN-producing input faults in the guarded test build.)*
2. **✅ FIXED (2026-06-06) — DMA overrun busy-waited inside the column ISR; the end-frame "spec mismatch" was not a bug.** Two coupled claims, validated separately:
   - **`END_FRAME_BYTES` odd-N mismatch — VALIDATED, NOT A BUG.** `(N/2)+1` (`dma_led.h:55`) is *exactly* `ceil((N+1)/2)` (the documented form, `:42`) for **all** N: even `N=2k` → `k+1 = ceil((2k+1)/2)`; odd `N=2k+1` → integer `N/2=k` → `k+1 = ceil((2k+2)/2)`. The comment and code agree, the buffer is correctly sized (constructor indexing verified in-bounds), and `ceil((N+1)/2)` *bytes* far exceeds the `ceil(N/2)` *bits* the protocol requires. **No code change.**
   - **In-ISR busy-wait on overrun — FIXED.** On overrun, `submitFrame()`/`show()` incremented `overrunCount_` then fell through to `transmitAsync()` (`dma_led.h:268-275`), which spins in `waitComplete()` (`:281-283`) on a flag set only by the DMA-completion ISR (`:286-291`). Called from the column `IntervalTimer` ISR, that spin can **hard-hang** (the DMA-completion IRQ can't preempt an equal/lower-priority ISR) or at best overrun the next column tick. **Fix:** drop the frame on overrun (early `return` before `transmitAsync`) in both `submitFrame()` and `show()` — exactly the transient/soft-handling the engine's own philosophy prescribes for dropped frames (`platform.h:27`). The in-flight DMA keeps the previous column; the next tick repacks and retries. **No perf cost** (drops work on the overrun path; normal path unchanged). *Hardware-only path (`#ifdef ARDUINO`), not in the native suite — validate on a Teensy build.*

### P1 — Latent bugs, fail-fast gaps, borrow safety

3. **The math/geometry core has no fail-fast and two real wrap bugs** (theme 1 table; `rotate.h:18-24`, `util.h:25-35`). Replace silent degenerate `normalize`/`angle_between` fallbacks and `NDEBUG`-stripped `assert`s with `HS_CHECK`; add bounds traps on LUT/`Orientation` indexing; fix `mod_tau` to wrap any period (`x - floorf(x/τ)*τ`) and drop `wrap()`'s `abs(x)<TOLERANCE` clamp.
4. **`Lerp` stores `const T*` to its `start`/`target` with no rvalue `=delete`** (`animation.h:702-711`, verified). `tl.add(0, Lerp(subj, T{...}, compute(), dur, ease))` compiles and dangles — the exact footgun `Motion`/`MeshMorph` are hardened against. **Fix:** add the matching `=delete` rvalue overload so a temporary start/target is a compile error.
5. **`Canvas` construction is a TOCTOU on relaxed atomics** (`canvas.h:296-304`). `while(!buffer_free()){}` then `advance_buffer()` are not in one interrupts-disabled window, and the relaxed ordering doesn't synchronize-with the ISR's `advance_display()`. Correct by external single-writer invariant, not by the memory model. **Fix:** wrap the check-plus-advance in the same `disable_interrupts()` window already used by `queue_frame`, or use acquire/release on `prev_`/`next_`.
6. **SDF seam/truncation gaps** (`scan.h:120-153`, `sdf.h:1039-1042`). The `scan_region` coalescer is a single forward sweep that isn't seam-wrap-aware (latent for CSG output crossing θ=0), and `Face` silently truncates >MAX_VERTS to wrong geometry. **Fix:** normalize emitted intervals into `[0,W)` (splitting seam-crossers) before coalescing; replace the `Face` clamp+log with `HS_CHECK`.
7. **Mesh non-manifold + index-width gaps** (`mesh.h:122-130`, `conway.h` int16 casts). >2 half-edges sharing an edge are mispaired with no trap, and vertex-index casts to `int16_t` overflow silently past 32767. **Fix:** `HS_CHECK` on the non-manifold fan and after large vertex-generation passes (or widen to `int32_t`).
8. **Plot planar/antipode edges** (`plot.h:90-108` planar `dist`; `Line::sample` antipodal `cross`). Planar mode feeds a projected chord into arc-length-calibrated thresholds; `Line::sample` normalizes a near-zero cross for antipodal points. **Fix:** compute `total_dist` as the true angular length for the planar strategy; add the antipodal branch (or route `Line` through the geodesic strategy).
9. **`setClip` forwards untyped JS ints with no validation** (`wasm.cpp:189-192` → `canvas.h` `set_clip` → `constants.h` `%` arithmetic). Negative `x0`/`x1` produce negative-modulo garbage clip logic. **Fix:** `HS_CHECK`/clamp bounds at the WASM boundary before forwarding.
10. **`Test`/`TestShapes` ship in the live registry** (`Test.h:108`, `TestShapes.h:200`; `REGISTER_EFFECT` appends under `__EMSCRIPTEN__`, verified). Dev scaffolding appears as selectable effects with sliders in the demo. **Fix:** gate both behind `#ifdef HS_DEV_EFFECTS` or drop the registration.
11. **Inert / conflicting effect params** (theme 4): Flyby `Falloff` (dead), SplineFlow `Drift` (read-once), Comets `Cycle Dur` (read-once), ChaoticStrings `next_preset()` stub, Liquid2D `lerp` dropping `cycle_speed`, GS `dt` default above its slider max. **Fix:** wire each to a live read / re-seed (mirroring Metaballs/Voronoi), or remove the control; bring GS's default into range.
12. **Color robustness/portability** (`color.h:399-406,570-573`). `constexpr srgb_to_linear_float`/`linear_to_srgb_float` call non-`constexpr` `powf` (compiles only because never used at compile time — a latent hard break), and `Gradient`'s linear fill casts to `uint16_t` without a clamp. **Fix:** drop the misleading `constexpr`; clamp the gradient channels.
13. **The test suite never exercises the philosophy it markets.** No death test fires any trap; the smoke harness runs on an 8 MB arena (`HS_TEST_BUILD`) so it *cannot* catch MindSplatter's documented 335 KB over-subscription; the terminal smoke assertion `HS_EXPECT_TRUE(acc+1 > 0)` is tautological; `fake_canvas()` is a reinterpret_cast to never-constructed storage (UB); the 28-entry smoke list is hand-maintained and can silently drift from the registry. **Fix:** add a process-per-case death test for one representative trap per surface; after each smoke render assert `get_high_water_mark() ≤ 335 KB`; replace the tautology with a real post-condition; cross-check the smoke count against the registry.

### P2 — Quality, consistency, cleanup

14. **IslamicStars hygiene** (`IslamicStars.h:54-93,178-199`): delete the `HS_OS_CYCLES`/`hs::log` profiling cruft from `draw_frame` (or gate behind `#ifdef HS_PROFILE`), and push the `persistent_arena.reset()` + `Persist` compaction + palette rebake into a `MeshCarousel` helper instead of the effect body.
15. **BZ/GS duplication** (`ReactionDiffusionBase.h`): the CRTP base hoisted only constants/orientation; the graph-Laplacian loop, kernel accumulation, and render skeleton are still copy-pasted. Lift a templated `graph_laplacian` / `render_skeleton(substeps, shader)` into the base.
16. **Metaballs/Voronoi brute force** (`Metaballs.h:70-93`, `Voronoi.h:63-119`): replace the `O(W·H·N)` per-pixel field loops with bounded-AABB rasterization (Metaballs) and a `core/spatial.h` KD-tree/hash nearest+second-nearest query (Voronoi), hoisting the per-pixel `acosf` behind the near-boundary test.
17. **MindSplatter over-subscribes the 335 KB device arena** (`MindSplatter.h:24-25`): ~316 KB pool (1024 × `Particle<23>`) against a ~324 KB persistent partition before attractors/emitters/LUT. Reduce `NUM_PARTICLES`/`TRAIL_LEN` to leave headroom; the smoke harness can't catch this until P1 #13 lands.
18. **README/code drift**: remove (or implement) the `deploy.yml` gated-deploy description; correct the `VectorTrail` "int16/uint8 quantized" claim (it stores full `Vector`); correct the "DMAMEM placement" claim for the DMA buffers.
19. **CI/build polish**: add `actions/cache` for ccache + emsdk (the repo wires ccache but CI never benefits); document that configuring the `tests` preset silently repoints `core.hooksPath`.
20. **Daydream polish** (`segment_controller.js:300-305`, `vendor-importmap.js:21-30`): the per-segment render stat reads microseconds but renders `… ms` (divide by 1000); replace the deprecated synchronous-XHR HEAD probes with a build-time importmap or async preflight; consolidate the two independent URL-sync writers.
21. **De-dup / consistency**: extract the slew-rate IOMUX block (3 sites) into one helper; factor the shared `Ring`/`DistortedRing` annular core; switch `Polygon`/`Star`/`SphericalPolygon` per-row `acosf`/`atan2f`/`cosf` to the `fast_*` family used by sibling shapes; fold the scattered `1e-5`/`0.001`/`0.0001` epsilons in `plot.h`/`3dmath.h` into the existing named `math::EPS_*` vocabulary; add `ArenaSpan` generation tracking to match `ArenaVector`.

---

## Comparative appraisal

**Versus hobbyist / maker LED-art projects (the immediate peer group):** Not comparable — several leagues above. The typical FastLED POV or addressable-LED project is a single `.ino` of gamma-corrected 8-bit blends, `delay()`-driven timing, and copy-pasted effect loops. Holosphere has a typed compile-time rendering architecture, perceptual color science, an arena memory model with fail-fast invariants, a half-edge mesh kernel with the full Conway-operator algebra, a WASM digital-twin simulator, and a real invariant-driven test suite. Top-0.1% for the category.

**Versus professional / commercial offerings:** It compares favorably with the *engine* quality of commissioned generative-LED installation software. Where a shipping product pulls ahead is the operational envelope this project hasn't needed: CI on real hardware, death-tested invariants, a build that doesn't compile its own guards out of the production artifact, a content-creation tool that isn't five hand-rolled HTML pages, and the QA pass that would have flagged the inert sliders and the dev-scaffold effects in the registry. The *core* is genuinely shippable; the *seams* (the math-core fail-fast gap, the build-flag guard loss, effect-param hygiene, test rigor) are where it reads as one expert's project rather than a team's product.

**Versus academic work:** The numerical maturity (matched projection sentinels, analytic Lipschitz bounds for safe sphere-tracing of warped SDFs, OKLab/OKLCH shortest-arc interpolation, Conway-operator + Hankin-pattern geometry, a provably-correct Fibonacci-lattice K-NN graph) is at or above a strong graphics/geometry-processing course project or a SIGGRAPH-adjacent demo. What it lacks is the academic apparatus — no novelty claim, no baseline evaluation, no ablation. It is *engineering* of known-good algorithms executed unusually well, not a research contribution. As the artifact accompanying a thesis on real-time spherical rendering it would be a standout.

**Versus the state of the art:** The individual techniques are well-chosen current practice rather than frontier (split LUTs, instanced rendering, embind zero-copy, EBO variadic templates). The *synthesis* — one unchanged C++ engine driving both microsecond-budget ISR firmware and a browser digital twin from a single source — is genuinely sophisticated systems design and is where the project is most impressive relative to anything off-the-shelf.

### Technical and artistic merit

**Technically**, this is the work of an expert who understands the whole stack — Cortex-M7 cache coherency, template metaprogramming, color appearance models, computational geometry, the WASM/JS FFI — and who maintains a *consistent design philosophy* across most of it. The recurring weakness is not ignorance but a boundary problem: the fail-fast rule that governs the memory and lookup layers stops at the numerical core and is compiled out of the shipped build, and the test suite never proves the rule holds. Close that gap — `HS_CHECK` in the math core, `-fno-finite-math-only` in the release build, one death test, one device-budget assertion — and the correctness grade moves from B− to A−, and the project with it.

**Artistically**, the effect catalog is ambitious and coherent: reaction-diffusion on a Fibonacci-lattice graph, the Hopf fibration with a 4D tumble, procedurally-generated authentic Islamic star patterns via Hankin's method, volumetric raymarched metallic tori at the vertices of a dodecahedron, Möbius-warped wireframes, real spherical-harmonic lobes with zero-crossing crossfade. These are not stock patterns, and the engine primitives (the SDF/scan/plot families, the transformer pool, the feedback styles) are expressive enough that the strongest effects (`Raymarch`, `MeshFeedback`, `HopfFibration`, the Liquid2D/Flyby shaders) achieve real visual sophistication in ~25–80 lines apiece. The medium itself — a 1,200-RPM persistence-of-vision sphere painting full-color 16-bit-linear imagery — is novel, and the engineering is in genuine service of the art.

**Bottom line:** **B+.** A remarkable, near-professional codebase. The prior review's hard-won fixes hold, but an independent pass shows the work isn't finished: the same fail-fast discipline that makes the memory layer excellent needs to reach the math core and survive into the shipped build, a short list of seam-edge bugs needs closing, and the test suite needs to finally exercise the invariants the project is built on. None of it is architectural — it is one philosophy, applied to two more boundaries, away from an unqualified A−.
