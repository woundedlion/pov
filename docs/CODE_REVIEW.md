# Holosphere — Code Quality Review

*Independent multi-agent audit, 2026-06-09. Ten reviewers each read one component end-to-end and cross-referenced the README to recover the intended contract; this report synthesizes their findings and forms its own overall appraisal. Out of scope by request: `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, and `core/rotate.h`; the third-party `core/FastNoiseLite.h` is treated as a vendored dependency. Letter grades (A+ … F) are measured against a high professional / graphics-middleware bar — not "does it run." The comparative appraisal in §7 is derived independently and does not reuse any prior review's wording or conclusions.*

---

## 1. Executive Summary

Holosphere is a single C++20 rendering engine (~50k LOC across `core/`, `effects/`, `hardware/`, `targets/`, plus the sibling `daydream/` JavaScript simulator) that compiles to two Teensy 4.x persistence-of-vision LED-sphere firmware targets — a 96×20 single-controller *Holosphere* and a four-microcontroller 288×144 *Phantasm* — and to a WebAssembly module that drives a Three.js browser simulator. The same effect classes, the same 335 KB arena, and the same per-frame `Pixel16[]` buffer run unchanged across all three targets; resolution is a compile-time `<W,H>` template parameter, so the device build carries no overhead from generality.

The audit's dominant finding is **uniformity of judgment**. Every component reflects the same authorial discipline:

- a fail-fast `HS_CHECK` doctrine placed deliberately on *cold* seams (allocation, registration, config, container growth) and kept out of the per-pixel hot loop, where a stripped `assert` backed by a cold bind-site trap is used instead;
- comments that explain *why* a non-obvious decision was made and which bug class it forecloses (the `lerp16` signed-MAC unsoundness, the seam-split rasterization rationale, the relaxed-atomics single-observer proof, the `-fno-finite-math-only` ordering);
- a 16-bit-linear / OKLCH perceptual color pipeline and a compile-time `Pipeline<W,H,Filters…>` filter chain that **multiple reviewers independently named the best things they saw**.

The safety story is not merely asserted — it is **mechanically verified**. A re-exec death harness asserts the *exact* illegal-instruction trap status (`SIGILL` / `STATUS_ILLEGAL_INSTRUCTION`) for 21 cold-path invariants; a determinism harness requires byte-identical re-renders under an injected mock clock and pinned RNG; the segmented-POV index math is *proven* to tile the canvas exactly once; every Conway operator is checked against the Euler characteristic V−E+F=2; and the WASM zero-copy memory-view lifetime — the single most dangerous contract in the whole product — is correctly engineered (pre-sized, never-resized backing buffers) on both the C++ and JS sides.

The weaknesses are equally uniform and unsurprising for a fast-moving solo project. Three themes recur across nearly every component:

1. **Duplication.** The Liquid2D/Flyby stereographic-shader pair (~70% shared, copy-pasted), the "preset-cycle + Lerp" idiom re-implemented across several effects, the "retain handle + live-apply slider" idiom in three or four, the `Scan::Shader` SSAA loop duplicated across its two overloads, an inline `quintic_kernel` re-implementation in DistortedRing, and the two near-identical Trails filters.
2. **Internal-consistency gaps against the project's own documented standards.** Five `int → uint8_t` face-count pushes bypass the `narrow_face_count` trap the codebase says should guard them (one site, `hankin.h:196`, is reachable in principle); several face-loop walks lack the orbit-spin guard their siblings carry; a few effects diverge on patterns two of them should share (HankinSolids' manual rotate loop vs IslamicStars' `OrientTransformer`; Liquid2D's `Driver` time source vs Flyby's plain accumulator).
3. **Testing and documentation that stop at the easy side of a boundary.** The engine harness is excellent, but a few small logic-bearing headers (`util.h`, `presets.h`, `styles.h`) have no dedicated tests, and the README's "divide by 65535 in JS" sample for the WASM→Three.js path contradicts the actual implementation (the divide is GPU-side via a normalized attribute).

Three latent correctness edges rise above cosmetic and deserve action regardless of priority order:

- **`Pixel::Feedback::flush` populates the full coarse warp grid ignoring the x-clip** (`filter.h:813-830`), so a feedback effect on segmented Phantasm hardware has every board recompute the whole sphere's warp field and use only its band — 4× wasted `space_fn` work, and a latent semantic trap if the clip were ever relied on for correctness rather than just culling.
- **The Phantasm segmented frame-sync ISR re-aligns the column counter without driving the corresponding buffer flip** (`pov_segmented.h:420-451`), which can make a board miss an `advance_display()` and show a stale buffer for up to half a revolution. This is the one finding that can visibly mis-paint real hardware.
- ~~**`Liquid2D`'s pattern time comes from an `Animation::Driver` that wraps at 1.0**~~ — *disproven on verification: the `Driver` is constructed with `wrap=false` (`Liquid2D.h:35`), so `accumulated_time` already accumulates unbounded exactly like Flyby's `phase`. No hitch, no fix.*

Neither of the first two crashes today; both are real (the third, a Liquid2D time-wrap, was investigated and disproven).

**Net: A− work.** A cohesive, philosophically consistent, exceptionally documented and verified engine whose remaining debt is ordinary — duplication, a few argued-rather-than-enforced invariants, doc drift, and a glue-testing backlog — not architectural rot.

---

## 2. Component Grades

| Component | Grade | One-line |
|---|---|---|
| Math & geometry (`3dmath`, `geometry`, `util`, `waves`, `easing`, `concepts`) | **A / A−** | Outstanding singularity/pole sentinel system and documented fast-math error bounds; unguarded `fast_wrap`, duplicated raw/templated `phi↔y` overloads, `Vector↔Spherical` trig-precision asymmetry. |
| Color system (`color`, `palettes`, `color_luts`) | **A−** | Correct 16-bit-linear / OKLCH numerics with the `lerp16` derivation documented; `BakedPalette::get` truncation bias contradicts the file's own rounding rule, `Color4` unbounded alpha, `GenerativePalette` low-`t` mis-map. |
| Memory & mesh (`memory`, `mesh`, `conway`, `hankin`, `spatial`, `solids`, `generators`, `static_circular_buffer`) | **A** | Wrap-safe arena math, fully-trapped truncation seams, exemplary scratch-arena contract docs; five un-narrowed face-count casts and a few unguarded face-loop walks deviate from the project's own discipline. |
| Rendering pipeline (`canvas`, `filter`, `sdf`, `scan`, `plot`, `transformers`, `styles`, `animation`) | **A** | The showpiece — compile-time variadic filter Pipeline with automatic 2D/3D domain bridging; feedback-flush ignores clip, `Scan::Shader` SSAA duplicated, the `Timeline` retained-handle contract is the one leaky abstraction. |
| Engine infrastructure (`platform`, `effect_registry`, `effects`, `presets`, `reaction_graph`, `led`, `constants`) | **A−** | Reference-quality `REGISTER_EFFECT` macro hygiene and X-macro SSOTs with bit-exact FastLED mocks; the 301 KB reaction-graph table's `PROGMEM` is a no-op on Cortex-M7 (RAM-placement unverified), `check_fail` breadcrumb can truncate. |
| Hardware drivers (`dma_led`, `hd107s_frame`, `pov_single`, `pov_segmented`, `pov_segment_map`) | **A− (corr. B+)** | Top-tier memory-ordering reasoning and host-testable protocol/index splits; the segmented frame-sync vs `advance_display` desync is a real on-hardware buffer-flip bug. |
| Effects, batch 1 (BZ/GS/RD base, Hopf, IslamicStars, HankinSolids, SphHarmonics, Metaballs, MobiusGrid, Moire, FlowField, Voronoi, Raymarch, Liquid2D, Flyby) | **A−** | Genuinely advanced math executed with understanding (S³ fibration, Legendre recurrence, graph reaction-diffusion, volumetric raymarch); Liquid2D/Flyby ~70% duplication, one unverified time-wrap, scattered dead declarations. |
| Effects, batch 2 (PetalFlow, DreamBalls, Comets, RingSpin, RingShower, ChaoticStrings, MeshFeedback, MindSplatter, Dynamo, Thrusters, GnomonicStars, SplineFlow, DistortedRing, ShapeShifter) | **A− / B+** | Disciplined determinism, byte-level pool accounting, explicit lifetime contracts; Comets rebakes a 256-entry palette every frame, RingShower's slot reuse is unguarded, repeated live-param-sync boilerplate. |
| Build targets & bridge (`wasm.cpp`, `param_marshal`, `Phantasm.ino`, `CMakeLists`, `CMakePresets`) | **A** | Correctly-engineered zero-copy view lifetime, SSOT roster X-macros, build-provenance SHA gate; the exceptions-abort-the-module contract is undocumented to JS callers. |
| Test suite (`run_tests`, `test_harness`, `test_death`, `test_effects`, +25 modules) | **A−** | Verified fail-fast death harness, byte-exact determinism, canvas-tiling proof, SPI wire-format pinning; no dedicated `util.h`/`presets.h`/`styles.h` tests, some misuse seams uncovered in every layer. |
| Daydream web simulator (`daydream`, `driver`, `segment_controller`, `segment_worker`, `state`, `recorder`, +tools) | **A / A−** | Disciplined frameworkless web-graphics: generation-fenced worker pool, transferable correctness, zero-alloc hot paths; fragile non-segment pixel-buffer re-aliasing, README divide-by-65535 divergence. |

---

## 3. Cross-Cutting Quality Dimensions

These grades are the whole-product synthesis across all components.

### Correctness & Robustness — **A−**
No memory-safety defects (use-after-free, OOB, alignment, arena overflow) were found in the engine core; the arena bounds math is provably wrap-safe and the comment proves it. Numerical edges are guarded with uncommon consistency — every stereographic projection, normalize, and division on a hot path has a documented epsilon or `normalized_or` fallback, and the NaN/Inf fault cluster is death-tested to *trap* rather than propagate. Held back from the A band by two real latent edges (feedback-flush clip, segmented frame-sync buffer-flip; a third candidate, a Liquid2D time-wrap, was investigated and disproven) and a scatter of low-severity issues: `BakedPalette` truncation bias, unbounded `Color4` alpha, and the unverified `PROGMEM`/RAM placement of the 301 KB reaction-graph table. *(The un-narrowed `hankin.h:196` face count has since been routed through `narrow_face_count`; see fix list #4.)*

### Architecture & Design Elegance — **A**
The compile-time `Pipeline<W,H,Filters…>` with automatic world/screen/pixel domain bridging via discarded `if constexpr` branches is a genuinely excellent abstraction — better than most production filter/scene-graph systems, and the reviewers did not grade it on a curve. The one-engine/three-targets seam is clean and leaks nothing into effect code; the owned-vs-borrowed `MeshState` duality, the purely-functional Conway operators over ping-ponged arenas, the CRTP `ReactionDiffusionBase`, and the POD `Feedback::Style` are all cohesive. The single leaky abstraction is the `Timeline` global-singleton retained-handle contract, whose "infinite events added before any finite event" invariant is enforced only by a trap deep in compaction.

### Interface Expressiveness — **A−**
`registerParam` → WASM → auto-generated GUI is a clean reflection seam; the `Fragment` shader interface with documented per-rasterizer register conventions, `SolidBuilder`'s fluent Conway chaining, the `.then()` animation sequencing, and `ArenaSpan`'s type-level borrow-vs-own distinction are all expressive. Docked for ergonomic seams a new contributor can misuse: the duplicated raw-`int`/templated `phi↔y` overloads, the `setClip(y0,y1,x0,x1)` Y-before-X argument order, the global `random()` mocks that compile fine but diverge from `hs::rand_*`, and the single-`onChange`-slot collision in the GUI wrapper.

### Performance & Resource Discipline — **A**
This is a standout. Compile-time resolution dispatch, `FunctionRef` zero-allocation type erasure, split-trig LUTs, fast-reject heuristics (ripple ~90% vertex skip, Voronoi KD-tree with a correctness proof, the `cos_eh` particle reject), branchless saturating color math, `fast_cbrt`/`fast_sinf` deployed *exactly* on the OKLab hot path and nowhere else, a single instanced draw call with GPU-side normalize and cull on the web side, and worker parallelism measured as max() not sum(). The memory story is rigorous to the byte (MindSplatter's 316 KB pool accounted against the 335 KB budget). Minor waste only: the feedback full-grid warp under clip, Comets' per-frame palette rebake, GnomonicStars' 2000 uncull'd star rasterizations, and a 2× cache flush in the segmented ISR.

### Documentation & Comments — **A+**
Consistently the highest-scoring dimension across every component. Comments encode *contracts and rationale*, not restatement — the smlad unsoundness, the sentinel-margin math, the LUT thread-safety contract, the scratch-arena asymmetry that a well-meaning refactor would otherwise "fix" and break, the relaxed-atomics proof, the `-fno-finite-math-only` ordering, and numerous "do NOT fix this" guards that document past bugs. The 1900-line README is genuinely architectural documentation, not a marketing page. The only material doc defects are *drift*: the WASM divide-by-65535 sample, a few stale effect parameter lists, and a couple of row-range wording mismatches.

### Testing & Verification — **A−**
A verified fail-fast death harness (re-exec, exact trap status), a determinism harness enforcing byte-identical re-renders, an auto-generated smoke pass that constructs and renders every effect at 288×144 with asserts on, hardware wire-format and segment-tiling proofs, and a WASM memory-view stability guard. Three continuous layers (pre-commit hook, presubmit CI, gated deploy) make a regression structurally unable to reach the live demo. Gaps are narrow and known: no dedicated `util.h`/`presets.h`/`styles.h` tests, several misuse preconditions uncovered in any layer, and no test for the `SegmentController` one-frame state machine or the eDMA handoff.

### Consistency with Stated Philosophies — **A**
The four design philosophies (16-bit linear color, compile-time resolution, fail-fast cold-only `HS_CHECK`, double-buffer-by-design) are applied with rare uniformity and *enforced at build/startup*, not just documented. The deductions are precisely the places where the code deviates from its *own* rules: the five un-narrowed face-count casts, the inline `quintic_kernel`, the `BakedPalette` truncation, and a handful of effect-to-effect pattern divergences.

---

## 4. Prioritized Fix List

**P0 — correctness on real hardware / latent semantic traps**

1. Fix the Phantasm frame-sync vs `advance_display` desync (`pov_segmented.h:420-451`): drive the buffer flip from the authoritative frame-sync pulse, or perform `advance_display()` when a snap crosses a boundary. *Only finding that mis-paints hardware.*
2. ✅ Make `Pixel::Feedback::flush` honor the x-clip when populating the coarse warp grid (`filter.h:813-830`) — correctness-of-intent on segmented hardware plus a 4× perf win. *Fixed: step 1 now marks and populates only the coarse columns the clipped sampling pass reads (`cx0`/seam-wrapping `cx1`); sampled output is byte-identical, the full-grid `space_fn` waste on out-of-band columns is removed.*
3. ✅ Confirm/unify `Liquid2D`'s `accumulated_time` Driver wrap semantics against `animation.h`; if it wraps at 1.0, switch to Flyby's unbounded accumulator (`Liquid2D.h:35,71`). *Confirmed — not a bug, no change. The `Driver` is constructed with `wrap=false` (`Liquid2D.h:35`), so `Driver::step` skips the `wrap_t` fold and `accumulated_time` grows unbounded — already equivalent to Flyby's inline `phase += speed`. The premise that it wraps at 1.0 is incorrect; nothing hitches.*
4. ✅ Route `hankin.h:196`'s `face_counts.push_back(count)` through `narrow_face_count` — the one reachable un-narrowed cast; apply to the other four sites (`conway.h:465/657/747/953`) for consistency. *Fixed: all five `int` face counts now narrow through the trapping `narrow_face_count`, matching the existing `conway.h:145/583` sites; valid solids are unaffected (tests green), an over-wide face now traps at the bench instead of silently wrapping.*

**P1 — verification & resource risks**

5. Confirm where the 301 KB reaction-graph table physically lands on Teensy (bare `PROGMEM` is a no-op on Cortex-M7); use an explicit flash-section attribute or document/`static_assert` the placement (`reaction_graph.cpp:4`).
6. ✅ Add the non-segment-mode pixel-buffer alias assertion in `daydream.js` mirroring `composite()`'s `dst === Daydream.pixels` check, or unconditionally re-point the three aliased views after a resize. *Fixed: the normal (single-engine) `drawFrame` path now asserts, right after `refreshPixelView()`, that all three aliases (`wasmMemoryView`, `Daydream.pixels`, `daydream.dotMesh.instanceColor.array`) point at the one WASM view — the same invariant `composite()` already enforces on the segment path. The lazy re-aliasing was previously unguarded here, so a future change re-pointing only some of the three would have silently cleared/displayed a buffer the engine never rendered into; it now throws loudly. One reference-equality check per frame, no hot-path cost.*
7. ✅ Document the WASM exceptions-abort-the-module contract on the JS-facing MeshOps surface (allocation failure aborts, not throws), or enable exception catching for the tooling path. *Fixed by documentation (no perf/size cost — exception catching is deliberately left disabled). The solids editor's `update()` now carries a contract comment at the point where MeshOps ops are invoked: because the module is built `-fno-exceptions`, an op that overflows the tooling arena or trips an HS_CHECK calls `abort()` and tears down the whole module rather than throwing, so the surrounding `try/catch` only catches Embind marshalling errors and must not be mistaken for overflow recovery. The note also records that the C++ wrapper pre-rejects the known recoverable cases (unknown solid name, oversized `fromData`) while an arbitrary op chain that outgrows the arena has no such guard (`tools/solids.html`).*
8. ✅ Enlarge `check_fail`'s format buffer to 256 and/or basename `__FILE__` so the on-device breadcrumb can't truncate (`platform.h:617`). *Fixed: did both. `check_fail` now strips the directory from `__FILE__` (the compiler bakes in a long absolute build path that otherwise crowds out the cond/message tail — the part that names the failure) before formatting, and the on-device `hs::log` buffer that the breadcrumb passes through grew 128→256 (the WASM `console.error` buffer 176→256 to match). A full path + long condition + 95-char message previously overran the 128-byte device buffer and dropped the message tail; basename + 256 covers the realistic worst case. Both changes sit on the cold trap path only — no hot-path cost. Tests green.*

**P2 — consistency & cosmetic correctness**

9. ✅ Make `BakedPalette::get` round via `float_to_pixel16` for parity with `Gradient::get` and the file's own anti-truncation rule (`color.h:1316`). *Fixed: the blend weight now rounds (`frac * 65535.0f + 0.5f`) instead of truncating, matching `float_to_pixel16`/`Gradient::get` and removing the ~1/65535 down-bias on every interpolated sample. `float_to_pixel16`'s clamp is intentionally not reused at this call site — `frac` is already in `[0,1)` after the bounds checks and `get()` is a per-pixel hot path — so the fix adds no per-pixel cost; tests green.*
10. Add the always-on orbit-spin `HS_CHECK` to `face_centroid`/`face_normal`/`face_side_count` (`conway.h`).
11. ✅ Add a stripped `assert` to `fast_wrap` and the non-negativity preconditions in `AntiAlias`/`ChromaticShift`. *Fixed: `fast_wrap` now traps its documented `[-W, 2*W)` single-step precondition, and `AntiAlias`/`ChromaticShift` `plot()` trap `age >= 0 && alpha >= 0`. All three are debug-only `assert`s — stripped under NDEBUG on device (zero hot-loop cost), fire in the native tests / WASM-debug; tests green.*
12. ✅ Gate Comets' per-frame palette rebake on an active color wipe. *Fixed: Comets now arms a frame countdown when the cycle timer schedules a `ColorWipe` and rebakes the 256-entry LUT only while that wipe is mutating `palette` — static between wipes, so the per-frame rebake is skipped the vast majority of frames (`Comets.h`).*
13. ✅ Make RingShower's recyclable-slot lifetime explicit (adopt Thrusters' age-driven model). *Fixed: RingShower drops the per-ring `Sprite`+`Transition` pair (which split a slot's lifetime across two captured animations and freed it implicitly via `.then`) for Thrusters' model: radius growth, fade-in and expiry are pure functions of a per-slot `age` advanced in `draw_frame()`, so a recycled slot can never be drawn by a stale animation (`RingShower.h`). Render values are preserved (radius/opacity sequences unchanged, whole ring shifted one frame earlier); tests green.*
14. ✅ Replace `innerHTML` interpolation of worker-fault messages / tool labels with `textContent` (the codebase already knows better elsewhere). *Fixed: `LabelPool.acquire` (`driver.js`) sets label text via `textContent`, and the segmented worker-fault panel (`segment_controller.js`) builds its `div`/`span` nodes and writes the arbitrary worker error message through `textContent` instead of interpolating it into `innerHTML`. Styling/layout unchanged; daydream tests green.*

**P3 — duplication & doc drift (quality, not bugs)**

15. Extract a `StereoShaderBase` CRTP for Liquid2D + Flyby (largest single dedup).
16. Factor the live-param-sync and orientation-trail-render idioms shared across several effects.
17. ✅ Unify the two `Scan::Shader` SSAA loops. *Fixed: the copy-pasted sub-pixel offset table, the per-sample theta/phi→vector projection, and the LUT-domain `HS_CHECK` are now shared `Shader` helpers (`make_sample_offsets`, `ssaa_sample_vector`, `check_lut_domain`) called by both `draw()` overloads (`scan.h`). The projection trig is kept verbatim so SSAA output is byte-identical (determinism harness green); the helpers are `constexpr`/`inline`, so codegen is unchanged.*
18. ✅ Call `quintic_kernel()` from DistortedRing instead of its inline re-implementation. *Fixed: the hand-expanded `t*t*t*(t*(t*6-15)+10)` smootherstep in the ring fragment shader now calls `quintic_kernel(1 - norm_dist)` (`DistortedRing.h`). `norm_dist` is already clamped to `[0,1]`, so `quintic_kernel`'s clamp is a no-op and the falloff is byte-identical; tests green.*
19. Reconcile README drift: the WASM divide-by-65535 sample, stale effect parameter lists (IslamicStars/Metaballs/Voronoi/MobiusGrid), and the §7.10 segmented row-range wording.
20. Add `test_util.h`, `test_presets.h`, and direct `styles.h` coverage.
21. Extend the death harness to the deferred memory/spatial/circular-buffer misuse preconditions.
22. Delete dead code: `ReactionDiffusionBase::H_VIRT`, the unused `<functional>` includes, the stale using-decls, and the unreferenced `hs::H_OFFSET`.

---

## 5. Overall Impression

Holosphere reads like the work of one strong engineer who has internalized a small number of principles and applied them everywhere, then *built the machinery to prove they hold*. The rare quality here is not any single clever file — it is that ten independent reviewers, each handed a different slice, came back with the **same** description of the author's values and the **same** shape of weakness. That coherence is itself the headline result. Most codebases of this size are an archaeology of three different past selves; this one has a single, legible voice.

The engineering is genuinely sophisticated where sophistication pays and restrained where it doesn't. The compile-time filter pipeline, the arena model with RAII compaction, the bit-exact device/simulator parity, and the verified fail-fast safety net are professional-grade. The fail-fast doctrine in particular is executed with a precision most teams never reach: not "asserts on," but *traps on cold seams, stripped-assert on hot paths backed by a cold bind-site trap, soft-handle genuine transients, and a death harness that proves the traps fire*. That is a coherent theory of failure for a headless device with no console — exactly the right theory for the problem.

The weaknesses are the weaknesses of velocity, not of judgment. Duplication accumulates because shipping the next effect beats refactoring the last two; documentation drifts because the code moves faster than the prose; a handful of invariants are argued in comments rather than enforced in types. None of it is rot. All of it is the ordinary debt of a solo project under aesthetic deadline pressure, and the audit surfaced no architectural decision that needs reversing.

---

## 6. Technical & Artistic Merit — The Achievement

The achievement is best understood as **three hard problems solved at once, sharing one codebase**:

1. **A real-time spherical rendering engine** with its own SDF rasterizer, geodesic curve plotter, volumetric ray-marcher, half-edge mesh kernel with the full Conway operator algebra, a perceptual (OKLCH) 16-bit-linear color pipeline, and a composable filter graph — all templated to zero-overhead at a fixed resolution.

2. **Hard-real-time POV firmware** that paints a sphere by spinning an LED strip at 480 RPM, firing columns on a microsecond ISR budget, double-buffered lock-free, and — in the Phantasm build — *coordinated across four microcontrollers* over two sync wires with branchless per-segment ISRs and deterministic cross-board RNG so all four paint an identical canvas.

3. **A browser digital twin** that compiles the identical C++ to WebAssembly, renders it through Three.js with zero-copy pixel readback, and reproduces the four-board partitioning in Web Workers so the hardware layout can be exercised in software before fabrication.

The artistic layer is not decoration bolted onto an engine — it is mathematically literate. The effects implement the *actual* Hopf fibration of S³ with a 4D tumble, the associated-Legendre recurrence for real spherical harmonics (with an integer-sqrt rounding proof for the mode decode), Belousov-Zhabotinsky and Gray-Scott reaction-diffusion on a precomputed 7680-node Fibonacci K-NN graph, authentic Hankin-method Islamic star patterns generated procedurally from Archimedean solids, twisted-torus volumetrics with metallic Blinn-Phong, and Möbius-warped grids with pole-stabilizing counter-rotation. These are not shader-toy copies; they are implemented by someone who understands the geometry.

---

## 7. Comparative Appraisal — Independent Assessment

*This section is formed from first principles against the named reference classes; it does not draw on any prior review.*

**Versus peer / hobbyist LED-art projects.** This is not in the same category. The overwhelming majority of POV and addressable-LED projects are single-target Arduino sketches with 8-bit sRGB blending, hand-tuned magic constants, no tests, and no separation between effect and driver. Holosphere has a perceptual linear-light pipeline, a templated rendering engine, a verified safety net, a CI-gated deploy, and a browser simulator that is byte-faithful to the device. On structure, rigor, and ambition it is one to two full tiers above the visible hobbyist ceiling (the FastLED/WLED ecosystem), and the multi-controller synchronized POV build is something that ecosystem essentially does not attempt.

**Versus professional offerings.** The closest commercial analogues are creative-coding and LED-mapping toolchains (TouchDesigner, MadMapper, Pixelblaze, the high-end architectural-LED controllers). Holosphere is narrower — it targets one geometry, the sphere — but within that niche its engineering practices are recognizably professional: zero-overhead compile-time specialization, a documented memory model, host-testable hardware protocol layers, build-provenance hashing, and a reflection-driven GUI. Where it trails professional middleware is breadth and team-hardening: no fuzzing, the glue layers (WASM marshaling, DMA handoff, the worker state machine) are under-tested, and there is duplication a code-review gate at a company would have caught. The *core* is professional quality; the *seams* show that one person wrote all of it.

**Versus academic research.** Individually, several techniques (geodesic SDF rasterization on a sphere, graph-based reaction-diffusion on a Fibonacci lattice, procedural Hankin patterns over Conway-operated polyhedra, multi-board POV synchronization) are each at the level of a solid graphics/fabrication workshop paper or a strong SIGGRAPH demo. What distinguishes this from research code is the opposite of the usual gap: research prototypes are typically clever and unmaintainable, whereas this is clever *and* engineered to ship and survive. As a research artifact it would be praised for reproducibility (the determinism harness, the simulator) more than for any single novel algorithm — none of the individual algorithms is new, but the *integration* is unusually complete.

**Versus state of the art.** It does not advance the rendering or color-science state of the art, and does not try to — OKLab, SDF CSG, sphere tracing, and Conway operators are all known art, applied correctly. Where it is plausibly at or near the state of the art is the *narrow intersection it actually occupies*: a perceptually-correct, compile-time-specialized, formally-fail-fast rendering engine that runs bit-identically across a four-microcontroller POV sphere and a WebAssembly browser twin. It is unlikely that a more rigorous instance of *that specific combination* exists in the open.

**Overall standing.** As a solo project, this is exceptional — top single-digit-percentile for the LED-art domain on every axis the audit measured, and professional-grade in its core engineering. The honest ceiling on the grade is set not by what is here but by what a hardened team practice would add: broader glue testing, mechanical de-duplication, and enforcement (in types or CI) of the invariants currently argued in comments. That ceiling is exactly why the synthesis lands at **A−** rather than **A** — not a deficiency of vision or craft, but the ordinary distance between superb individual work and battle-hardened collective work.

---

*Grades reflect the state of the tree at report time. Re-run the audit after the P0/P1 fixes land; several would move their component into the A band.*
