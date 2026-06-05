# Holosphere + Daydream — Code Quality Review

**Date:** 2026-06-05
**Reviewer:** Expert C++/TypeScript code audit
**Scope:** ~50k lines C++ (engine, firmware, 27 effects, WASM bridge, tests) across `Holosphere`, plus the ~2.4k-line vanilla-JS `daydream` web simulator. Conducted by **8 parallel component audits**, each reading source directly and citing `file:line`. Generated tables (`reaction_graph.cpp`, `color_luts.h`) and third-party code (`FastNoiseLite.h`) were judged on integration and correctness of generation, not line-by-line. **Out of scope (per review rules):** `effects_legacy.h` and `targets/Holosphere/Holosphere.ino`.

---

## Overall grade: **A− / B+**

This is an unusually sophisticated and well-engineered codebase for what is effectively a solo generative-art project. The architecture is principled and consistent (compile-time `<W,H>` specialization, 16-bit linear color, partitioned arena allocation, a variadic filter pipeline with automatic coordinate-domain coercion), the performance engineering is genuinely expert-level (LUT-driven trig, fast-reject heuristics, zero-copy WASM readback, non-blocking DMA, ISR double-buffering), and the documentation is exceptional. The project also has a clearly-articulated and *largely well-executed* **fail-fast philosophy**: invariant violations hard-fault via an always-on `HS_CHECK` trap that survives `NDEBUG`, while genuinely transient conditions degrade gracefully.

It is held back from an unqualified A by three things: (1) a **pixel-sink hardening gap** where projected/fractional coordinates can in principle write out of bounds on the `NDEBUG` device before any trap fires; (2) a recurring **`assert`-vs-`HS_CHECK` inconsistency** — several invariant checks still use `assert`, which is stripped on device, partially undermining the otherwise-excellent fail-fast posture; and (3) a small cluster of **concrete correctness bugs** — two 16-bit→8-bit color truncations in shipping effects, a WASM format-string defect, three effects that abort at 288×144, and a couple of dangling-reference hazards in the animation layer.

---

## Dimension grades (aggregated across all components)

| Dimension | Grade | Justification |
|---|---|---|
| **Architecture & Design** | **A−** | Compile-time resolution, clean `platform.h` abstraction, the variadic filter pipeline with transparent 2D↔3D domain coercion, explicit-`Arena&` purity, the Conway-operator algebra + `SolidBuilder`, the "animations mutate state, don't render" separation, and the JS pub/sub state layer are all elegant and orthogonal. Consistent across the whole stack. |
| **Readability & Maintainability** | **B+** | Outstanding Doxygen, named tolerances, error-bounded comments, intent-revealing rationale. Dragged down by copy-paste (edge-pairing heap-sorts ×5, MeshOps clone divergence, interval-merge duplication, duplicated SSAA/metrics lambdas), a few 200–270-line monolith functions (`sdf.h` `Face` ctor), and profiling cruft in hot paths. |
| **Correctness & Robustness** | **B / B−** | Strongly defensive at poles/wrap/NaN/div-by-zero in the math core, and most edge cases carry comments showing prior fixes. But real hazards remain: the pixel-sink single-step wrap, two effect color truncations, dangling path/draw references in `animation.h`, a WASM format-string defect, and 3 effects that trap at full resolution. |
| **Memory Safety & Resource Mgmt** | **A− / B+** | Disciplined arena ownership; `Arena::allocate` and `ArenaVector` are model fail-fast citizens (trap at the violation site, never hand back null/OOB), which makes the many raw `allocate()` derefs safe by construction. `Persist`/`ScratchScope` RAII and Three.js GPU disposal are correct. Undercut only by the residual `assert`-only guards on cold paths and a couple of effect/arena over-subscriptions. |
| **Performance & Efficiency** | **A−** | The best dimension. Microsecond-budget awareness everywhere: split TrigLUTs, fast-reject culling, anisotropic per-face distance LUTs, non-blocking DMA with single cache-flush, direct `Pixel16`→wire packing, zero-alloc steady-state render loops, single-draw-call instanced web mesh. |
| **Testing** | **B / B+** | Excellent native harness with *real algebraic invariants* (Euler characteristic, half-edge symmetry, partition-of-unity, Hamilton algebra, error-bounded epsilons vs `std::`). Held back by a tautological effect smoke-test post-condition, 3 quarantined effects, untested WASM-bridge/ISR/integration layers, and fail-fast traps that are never themselves exercised. |
| **Documentation** | **A−** | The README is genuinely outstanding — architecture diagrams, data-flow lifecycles, register conventions, rationale for every major decision. Code comments quantify error bounds. Gaps: the JS-facing zero-copy bridge contract isn't documented in a JS-readable place, and a few operator preconditions (closed-manifold `E=I/2`) are implicit. |

### Per-component snapshot

| Component | Arch | Correctness | Perf | Standout note |
|---|---|---|---|---|
| Core Math & Color | A− | B | A | `3dmath.h` is the strongest file in the repo; provably non-overflowing `lerp16` |
| Rendering Pipeline | A− | C+ | A− | Elegant variadic pipeline; the x-sink wrap is the one real hazard |
| Animation / Transformers | A− | B | A | Allocation-free type-erasure; two dangling-reference hazards |
| Mesh System | A− | B | A− | Conway algebra + `SolidBuilder` are beautiful; orbit `assert`s should be `HS_CHECK` |
| Memory / Platform / HW | A− | B | A | Branchless segmented ISR; `effect_` not `volatile`, driver ownership divergence |
| Effects (×27) | A− | B− | A− | World-class art; two color truncations + 3 capacity aborts |
| WASM Bridge / Build / Tests | A− | B | A | Correct detachment contract; format-string bug; shallow effect tests |
| Daydream Simulator | A− | B | A | Detached-view guard done correctly; worker-after-resolution blank-frame bug |

---

## Failure philosophy (the standard these findings are measured against)

This codebase runs on hardware in a **fail-fast** posture, and the recommendations below assume it:

- **Invariant violations must hard-fault, not degrade.** Arena over-allocation, container capacity overflow, OOM, or an out-of-bounds index is a logic/sizing bug with no valid recovery — a truncated or aliased allocation just relocates the corruption into the rendered output. These must crash at the violation site so they're caught on the bench, not ship silently.
- **The crash must survive the release build.** The hardware runs `-O3 -DNDEBUG`. `assert()` is therefore the *wrong* tool — it's stripped, and historically pulled in newlib's `__assert_func → fprintf`. The right tool is the project's existing `HS_CHECK(c)` (a `__builtin_trap()` on a predicted-not-taken cold-path branch, zero hot-loop cost).
- **Bounded/soft handling is correct only for genuine *transient* conditions** — a DMA overrun dropping a frame, a missed sync pulse, an unbound particle pool early-returning. Distinguish "should never happen" (trap) from "happens occasionally by design" (degrade).

The engine *mostly* gets this right: `Arena::allocate`, `ArenaVector::push_back/append_bulk`, the `vertex_orbit` anti-hang trap, and the plot `_steps_cache` guard are all proper `HS_CHECK` traps. The remaining gaps — pixel-sink bounds, in-orbit `assert`s, `Canvas::operator()` bounds, math-core silent fallbacks — are where invariant violations still slip past the device build, and they are the spine of the P0/P1 list below.

---

## Prioritized fix list

### P0 — Can corrupt memory or crash a live show

1. **✅ FIXED (2026-06-05) — debug-only precondition tripwire, zero device cost.** Traced every producer feeding the terminal sink before changing anything: the 3D path full-wraps via `vector_to_pixel` (`geometry.h:281` uses the safe `wrap()`), and `AntiAlias`/`Blur`/`ChromaticShift` keep x within `[-1, W]` before forwarding (`filter.h` AntiAlias `x0/x1` come from `modf` of an already-wrapped coord; Blur/ChromaticShift use full `wrap()`). So `fast_wrap`'s `[-W, 2W)` precondition **holds for every finite coordinate** — the OOB is unreachable in practice and the only real trigger is a non-finite (`NaN`) coord reaching the float sink, where `static_cast<int>(round(NaN))` is UB. Because `fast_wrap` is on the hottest path of all, an always-on `HS_CHECK` branch was rejected (per the perf constraint). Instead the precondition is now enforced **debug-only** (`assert`) at both terminal sinks (`filter.h` int sink: `assert(x >= -W && x < 2*W)`; float sink: `assert(std::isfinite(x) && std::isfinite(y))` + range): stripped under `NDEBUG` on device (literally zero hot-loop cost), but fires in the native test suite and WASM-debug build — exactly where a NaN-producing effect or a future out-of-range filter is first exercised. **Verified:** full native suite (every effect rendered at 288×144 with asserts on, plus all filter/plot/scan tests) passes with the asserts active, confirming the invariant holds and the tripwire is non-spurious. The genuinely correct place to stop a `NaN` is its source (an effect's 0/0), which this tripwire will surface; the per-pixel sink stays lean. Original finding below.

   ~~**Harden the pixel sink against out-of-range coordinates.** `filter.h:85-103` (`Pipeline<W,H>::plot`). The x-sink uses `fast_wrap(xi, W)`, which subtracts/adds `W` *exactly once*. Any coordinate with `|x| ≥ 2W` (a runaway procedural value, a future scaling 2D filter, or a `NaN→int` cast) leaves `xi` still out of range; `contains_x` returns `true` on a full canvas, `plot_virtual` checks **y only**, and `canvas(xi,y)` then indexes the buffer out of bounds — guarded solely by an `assert` in `Canvas::operator()` (`canvas.h:317`) that is stripped on device. **Fix:** use the full `wrap()` (not `fast_wrap`) for x in both sinks, reject non-finite coordinates, and add an always-on `HS_CHECK(xi>=0 && xi<W && y>=0 && y<H)` immediately before the store. This is the single path where projected/fractional coordinates can corrupt memory on hardware.~~

### P1 — Wrong output, latent crash, or fail-fast gap

2. **✅ FIXED (2026-06-05).** Two 16-bit→8-bit color truncations in shipping effects, both fixed by widening the cast to `uint16_t` (the channel-wise lerp/scale stays in `[0, 65535]`, so this is behavior-preserving; the separate rounding-bias item is tracked at P2 #9):
   - `Voronoi.h:90-95` — smoothed cell color was `static_cast<uint8_t>(...)` on `uint16_t` linear channels (0–65535). With ~200 sites the smoothing path runs on essentially every pixel, so colors wrapped mod-256 → mostly-dark garbage. Now `static_cast<uint16_t>`.
   - `SphericalHarmonics.h:167` — negative-lobe green channel was force-cast to `uint8_t`, truncating up to 65535. Now `static_cast<uint16_t>(pos.color.g * 0.8f)`.
   - **Verified:** native suite passes (both effects render in the 288×144 smoke harness).

3. **Three effects abort at 288×144 (quarantined: SplineFlow, TestShapes, Thrusters).** All trip *fixed-capacity* guards in the plot path at full resolution: SplineFlow overflows its `MAX_TRAILS`/trail buffer with closed-geodesic sample fragments; TestShapes overflows `_steps_cache` (`plot.h:192`); Thrusters faults in the `DistortedRing`/`Ring` plot path. These are correct fail-fast traps, but they make three shipping effects unrenderable standalone at Phantasm resolution. **Fix:** size the plot caches/trail buffers off `W` (as the Face dist-LUT already does, commit `f0ac2ea`) or clamp sample counts; then remove from the quarantine list. **Confirm whether they also reproduce on-device** — under the failure philosophy these are live-show crashes waiting to happen.

4. **WASM `setEffect` format-string defect.** `wasm.cpp:154-156` builds a string with `snprintf` then passes it as the *format* argument to `hs::log` (which is `vprintf`-based). An effect name containing `%` makes `vprintf` read nonexistent varargs → UB/crash. **Fix:** `hs::log("WASM: setEffect called with %s", name.c_str());` and drop the intermediate buffer.

5. **Dangling references in the animation layer.**
   - `animation.h:799` — `Motion` captures its `Path`/`ProceduralPath` by reference; a temporary argument dangles and `step()` reads freed memory every frame. **Fix:** own the path by value / move it in.
   - `animation.h:1253,1331` — `MeshMorph` stores non-owning `FunctionRef` draw callbacks but lives across many frames; temporary lambdas dangle after the constructor returns. **Fix:** use an owning callable type (like `Sprite::SpriteFn`) or hard-document the outlive requirement.

6. **`MindSplatter` over-subscribes the 335 KB device arena.** Its 1024-particle system + baked palette + pools only fit because `HS_TEST_BUILD` widens the arena to 8 MB; on-device it is a hard arena-overflow fault. **Fix:** reduce particle count / trail history under the device budget.

7. **Promote device-stripped `assert`s on invariant paths to `HS_CHECK`.** The fail-fast posture is undermined wherever invariant checks still use `assert` (stripped on `NDEBUG`):
   - `Canvas::operator()` bounds (`canvas.h:317,329,340`) — the buffer's most safety-critical accessor.
   - Hand-rolled mesh-orbit guards (`hankin.h:193,200`; `conway.h` `relax` loop) that bypass `vertex_orbit`'s always-on trap, so a non-manifold input can overflow scratch silently.
   - Math-core silent fallbacks in `Vector::normalize` / `Quaternion::normalize` / `angle_between` (`3dmath.h:198,398,715`) that log-and-substitute instead of trapping a zero-length input.

8. **Daydream: segments render blank/stale after a resolution change.** `workerSetResolution` clears each worker's effect, but `applyResolution` never re-sends the current effect, so segmented mode shows blank until the next effect switch. **Fix:** always `workerSetEffect(currentEffect)` after `workerSetResolution` (or have the worker re-`setEffect` itself).

### P2 — Quality, robustness hardening, maintainability

9. **Systematic rounding bias in color lerps.** `color.h:159,1190` (`frac` truncates) and `color.h:367` (`lerp8` floors and can wrap for out-of-range `t`). Add `+0.5f` and clamp `t`.
10. **`effect_` is a non-`volatile`, non-atomic pointer shared with ISRs** (`pov_single.h`, `pov_segmented.h`). Teardown ordering is currently safe, but the compiler may cache it across the ISR body. Qualify as `volatile Effect*` (or `std::atomic`).
11. **Unify effect ownership across the two POV drivers** — `pov_single::show` deletes after `run()`, `pov_segmented::run` deletes internally. Pick one owner to remove the latent double-free/leak divergence; drop the redundant `effect_ = e` (`pov_single.h:92`).
12. **`POVSegmented` ID decode assumes power-of-two `N`** (`& (N-1)`) but only `static_assert`s `N % 2 == 0`. Add `static_assert((N & (N-1))==0 && N <= 4)`. Also gate the per-frame `Serial.print` in the `pov_single` main loop behind `hs::debug`.
13. **`Star` vertical-bounds `h_virt` mismatch** (`scan.h:363` passes `H,H` without `H_OFFSET`), shifting its extent by the offset and clipping/AA-truncating at one pole. Pass `H + hs::H_OFFSET`.
14. **`Path::append_segment` unchecked ring-buffer overflow** (`animation.h:44`) silently overwrites oldest points instead of trapping. Add an `HS_CHECK` capacity guard.
15. **`compile_hankin::getMidpointIdx` boundary deref** (`hankin.h:104`) reads `halfEdges[HE_NONE]` when both `prev` and `pair` are `HE_NONE`. Unreachable on closed seeds but internally inconsistent; guard it.
16. **Missing div-guards / dead checks in effects:** `MobiusGrid.h:115` (`/(size-1)` without the `size>1` guard present at line 156); `MindSplatter.h:181` dead `size_t < 0` check. Sweep stray `<map>/<memory>` includes from MobiusGrid/Moire/RingShower.
17. **De-duplicate hot maintenance liabilities:** the edge-pairing `make_heap`+`sort_heap` block (5+ copies across `mesh.h`/`conway.h`), the two divergent `clone` paths (`spatial.h:486` vs `mesh.h:248`), the duplicated WASM `add_metrics` lambda, and the SSAA loops in `Scan::Shader`.
18. **Testing depth:** the effect smoke-test asserts only `acc+1 > 0` (`test_effects.h:85`) — add a mockable clock seam and at least one non-trivial property (golden hash / not-all-black). Add death-tests that actually fire the memory traps. Add a small harness for the WASM bridge's `fromData` parsing and detachment contract (currently validated only in the separate JS repo).
19. **Daydream hardening:** replace the synchronous-XHR importmap probe (`vendor-importmap.js:23-25`, blocks page load, deprecated) with a build-resolved map; re-validate `wasmMemoryView.buffer.byteLength` in `compositeSegments` (`daydream.js:243`) or assert the main engine is idle in segment mode.
20. **Reconcile gnomonic equator threshold** (`3dmath.h:568` 1e-9 band vs the `STEREO_INF` 1e4 sentinel) so near-equator points round-trip correctly.

---

## Overall impression — technical and artistic merit

**Technical merit: very high.** This is the work of an engineer who understands their constraints to the microsecond and has built an architecture that earns its complexity. The compile-time `<W,H>` specialization is the right call for embedded firmware; the 16-bit linear color pipeline is correct where most LED-art codebases are casually wrong; the partitioned-arena memory model with explicit `Arena&` threading is disciplined functional-style memory management rarely seen outside systems code; and the variadic filter pipeline that automatically lifts 2D↔3D coordinates at compile time is a genuinely elegant piece of template design. The fail-fast philosophy is not just stated but *mostly implemented* — `Arena::allocate` trapping at the violation site is what makes dozens of "unchecked" downstream derefs provably safe. The reaction-diffusion-on-a-Fibonacci-lattice work, the Hopf fibration, the Hankin-method Islamic star patterns generated live from Conway-operated polyhedra, and the SDF volumetric raymarcher are all serious computational-geometry achievements running in real time on a 600 MHz microcontroller.

The codebase's weaknesses are the *normal* weaknesses of an ambitious solo project: a handful of places where the discipline lapsed (`assert` where `HS_CHECK` was meant, two stray 8-bit truncations, dangling references from convenience-captured temporaries), some copy-paste that wants extracting, and a test suite whose harness quality (real algebraic invariants) outruns its coverage (effects, integration, and the bridge are thin). None of these are architectural; they are the punch-list of a strong system that needs a hardening pass.

**Artistic merit: exceptional for the medium.** A persistence-of-vision LED sphere is an unusually demanding canvas — the geometry is spherical, the time budget is brutal, and the output is physical light, not pixels behind glass. The effect library is not a collection of demos but a coherent aesthetic: mathematically-derived motion (Möbius warps, spherical harmonics, fibrations), perceptually-correct color in OKLCH, and motion blur derived for free from the orientation history. The decision to render in linear light and blend in 16-bit is an *artistic* choice as much as a technical one — it's why the gradients and multi-layer composites look clean rather than muddy. The web simulator compiling the identical C++ to WASM so the art can be previewed and shared without the hardware is the kind of end-to-end thinking that elevates the whole project.

**Bottom line:** an **A− / B+** codebase — top-decile engineering and artistry for a solo embedded-art project, with a clear, bounded, and mostly-mechanical path to an unqualified A. Land the P0 sink hardening and the P1 correctness/fail-fast cluster, and this is an A-grade system end to end.
