# Holosphere — Code Quality Review

*Independent multi-agent audit, 2026-06-09. Eleven reviewers each read one component end-to-end and cross-referenced the README to recover the intended contract; this report synthesizes their findings and forms its own overall appraisal. Out of scope by request: `core/effects_legacy.h`, `targets/*/*.ino`, `core/rotate.h`, and the third-party `core/FastNoiseLite.h`. Letter grades (A+ … F) are measured against a high professional / graphics-middleware bar — not "does it run." This review supersedes the prior one and re-derives every grade and the comparative appraisal from scratch.*

---

## 1. Executive Summary

Holosphere is a single C++20 rendering engine (~53k LOC across `core/`, `effects/`, `hardware/`, `targets/`) that compiles to two Teensy 4.x POV-sphere firmware targets — a 96×20 single-controller *Holosphere* and a four-microcontroller 288×144 *Phantasm* — and to a WebAssembly module that drives the `daydream` Three.js browser simulator. The same effect classes, the same 335 KB arena, and the same per-frame `Pixel16[]` buffer run unchanged across all three targets; resolution is a compile-time `<W,H>` template parameter, so the device build carries zero generality overhead.

The audit's dominant finding is **uniformity of judgment**. Every component reflects the same authorial discipline: a fail-fast `HS_CHECK` doctrine placed deliberately on *cold* seams (allocation, registration, config, container growth) and kept out of the per-pixel hot loop; comments that explain *why* a non-obvious decision was made and which bug class it forecloses; a 16-bit-linear / OKLCH color pipeline and a compile-time `Pipeline<W,H,Filters…>` that multiple reviewers independently named the best things they saw. The safety story is not merely asserted — it is **mechanically verified**: a re-exec death harness asserts the *exact* illegal-instruction trap status (SIGILL / `STATUS_ILLEGAL_INSTRUCTION`), a determinism harness requires byte-identical re-renders under a mock clock, the segmented-POV index math is *proven* to tile the canvas exactly once, every Conway operator is checked against the Euler characteristic, and the SDF interval-cull is verified against a brute-force exact-distance scan (a property the prior review listed as untested — it is now covered).

The weaknesses are equally uniform and unsurprising for a fast-moving solo project. Three themes recur across nearly every component:

1. **Duplication.** The Liquid2D/Flyby stereographic-shader pair (~80% shared, copy-pasted), the "preset-cycle + Lerp" idiom re-implemented in four effects, the "retain handle + live-apply slider" idiom in three, the SDF shape-family bounds/setup blocks, the Plot primitive close-loop skeletons, and the two Trails filters.
2. **README ↔ code drift.** Several effects' documented parameter lists are stale (IslamicStars, Metaballs, Voronoi, MobiusGrid expose registered params the README omits); the Dynamo description describes a different effect entirely; the RingSpin/RingShower blurb fits neither; `getParameterDefinitions()`'s documented shape is wrong; the AntiAlias "quintic on both axes" claim overstates the code (X is linear-after-sin); and the README oversells tool/WASM sharing (only `solids.html` actually loads the WASM `MeshOps`; the other four tools re-implement the math in JS).
3. **Testing that stops at the easy side of the boundary.** The engine harness is excellent, but the most failure-prone glue — the WASM `emscripten::val` marshaling, the `dma_led.h` eDMA handoff, and `SegmentController.tick()`'s one-frame-deep state machine — is untested.

Two latent correctness edges rise above cosmetic: `lerp_oklch` does not clamp interpolated L/C (extrapolated `t` can yield wrong/near-black colors), and `Pixel::Feedback::flush` ignores the clip region, so a feedback effect on segmented hardware would have every board composite the whole sphere into its own band. Neither crashes today; both are real.

**Net: A− work.** A cohesive, philosophically consistent, exceptionally documented and verified engine whose remaining debt is ordinary — duplication, a few argued-rather-than-enforced invariants, doc drift, and a glue-testing backlog — not architectural rot.

---

## 2. Component Grades

| Component | Grade | One-line |
|---|---|---|
| Math & geometry (`3dmath`, `geometry`, `util`, `waves`, `easing`, `concepts`, `constants`) | **A−** | Rigorous sentinel/degeneracy contracts and measured fast-math approximations; lossy gnomonic hemisphere contract, `square_wave` `fmod` negative-phase bug, `phi↔y` duplication. |
| Color system (`color`, `palettes`, `styles`, `color_luts`) | **A−** | Correct 16-bit-linear/OKLCH numerics with the `lerp16` derivation documented; `lerp_oklch` unclamped L/C, `blend_alpha` cast-before-clamp, residual truncation paths bypass `float_to_pixel16`. |
| Rasterization & curves (`canvas`, `filter`, `sdf`, `scan`, `plot`) | **A−** | Elegant compile-time domain auto-conversion, per-face Lipschitz LUT, arc-uniform stepping, rigorous double-buffer proof; latent scratch-arena aliasing hazard, Feedback ignores clip, heavy shape/primitive duplication. |
| Animation & transformers (`animation`, `transformers`) | **A−** | Compile-time borrow contracts and a fail-fast relocation latch; global-singleton `Timeline` with a dead `CAPACITY` param, inconsistent rvalue-delete coverage, Möbius-family duplication. |
| Memory & data structures (`memory`, `static_circular_buffer`, `spatial`, `presets`, `generators`) | **A−** | Wrap-safe arena math + dual-stamp use-after-free detection; `KDTree::build` silently drops subtrees on exhaustion, `int` count narrowing, OOM-path `printf`. |
| Mesh system (`mesh`, `conway`, `hankin`, `solids`, `reaction_graph`) | **A−** | Half-edge scaffold reuse collapses a dozen operators into a few primitives; `uint8_t` face-count narrowing unguarded, implicit closed-manifold precondition, `Solids::get(index)` is WASM-only. |
| Effects — Group A (RD/Hopf/Islamic/Hankin/Harmonics/Metaballs/Mobius/Moire/FlowField/Voronoi/Petal/DreamBalls) | **A−** | CRTP RD base and exact Voronoi KD-tree are A-grade; BZ hand-rolls a color composite that bypasses `Color4`, README params stale, init-order fragility. |
| Effects — Group B (Comets/Rings/ChaoticStrings/MeshFeedback/Liquid2D/MindSplatter/Dynamo/Thrusters/GnomonicStars/Raymarch/Flyby/SplineFlow/DistortedRing/ShapeShifter) | **B+** | Disciplined lifetime/arena hygiene; Liquid2D/Flyby and the preset-cycle idiom are heavily duplicated, `10000` triplicated in Dynamo, Dynamo/Ring docs wrong. |
| Hardware drivers (`dma_led`, `hd107s_frame`, `pov_single`, `pov_segmented`, `pov_segment_map`) | **A** | Standout DMA/cache-coherency and lock-free real-time reasoning, exact-tiling branchless ISR, host-testable pure logic; unconditional per-frame serial print, possibly-dead CRGB `show()` path. |
| Platform & WASM bridge (`platform`, `effect_registry`, `effects*`, `led`, `wasm.cpp`, `param_marshal`) | **A** | X-macro resolution SSOT, reject-don't-trap boundary, detached-view safety, bit-faithful FastLED mocks, LTO-hardened registrar; doc drift on param shape, unguarded `random(min,max)`. |
| Test suite (`tests/`) | **A** | Verified death + determinism harnesses, Euler/manifold invariants, cull-vs-brute-force equivalence, exact-tiling proof; WASM/DMA/ISR glue untested, order-coupled global fixtures, duplicated solid fixtures. |
| Web simulator (`daydream/`) | **A−** | Generation-fenced worker pipeline, fault latch, and detached-view handling are production-grade; README oversells tool/WASM sharing, `tick()` orchestration untested, Three.js-internal coupling. |

---

## 3. Per-Dimension Assessment

### Architecture & Design / Elegance — **A**
The strongest and most consistent dimension. The recursive `Pipeline<W,H,Head,Tail…>` (`filter.h:143-255`) classifies each filter by coordinate domain and re-dispatches into itself to lift (2D→3D via `pixel_to_vector`) or project (3D→2D via `vector_to_pixel`) a mismatched coordinate, with zero runtime branches and a `static_assert` that turns a non-terminal or mis-ordered pipeline into a compile error (`filter.h:213-234`). It is the cleanest possible expression of the README's "automatic world/screen/pixel conversion" claim. It is matched by the CRTP `ReactionDiffusionBase` (genuinely zero-overhead lattice iteration, `ReactionDiffusionBase.h:71-106`), the half-edge Conway scaffolds (`vertex_orbit`/`emit_vertex_orbit_faces`/`for_each_edge`, `conway.h:75-188`) that collapse a dozen operators into a handful of primitives, `WarpedVolume<Shape,Warp>` policy composition, the split-trig LUT (~145× memory reduction on the embedded target), the `Transformer<…,TransformFunc,…>` function-pointer pool, and on the JS side the `SegmentController`/worker split with generation-fenced dispatch. The lone B is Effects-Group B, where leaf-class duplication (Liquid2D/Flyby, the preset-cycle idiom) holds back otherwise-clean code.

### Correctness & Robustness — **A−**
Overwhelmingly sound, with correctness reasoning that is *better documented than most production code*: the `lerp16` round-to-nearest derivation and the removal of a broken signed-`smlad` dual-MAC (`color.h:91-119`), the arena no-wrap subtractive-bounds proof (`memory.h:66-88`), the per-face LUT 1-Lipschitz argument (`sdf.h:1290-1297`), the Canvas double-buffer TOCTOU analysis (`canvas.h:380-421`), the DMA/atomic decoupling note (`dma_led.h:145-155`), and the SDF interval-cull verified against brute force (`test_sdf.h:570-630`). Held below A by a cluster of latent edges: `lerp_oklch` not clamping L/C on extrapolated `t` (`color.h:509-526`); `Pixel::Feedback::flush` ignoring `canvas.clip()` (`filter.h:817-857`); `blend_alpha`'s float→int cast before the clamp (`color.h:268-271`); `KDTree::build` silently storing a `-1` child on pool exhaustion (`spatial.h:219-263`, currently a dead-but-silent guard); `square_wave`'s raw `fmod` misbehaving for negative phase (`waves.h:59`); and unguarded `uint8_t` face-count narrowings (`conway.h:145,583`, `hankin.h:196,238`).

### Interface Expressiveness & API Design — **A−**
`FunctionRef`/`PipelineRef` zero-allocation type erasure, the `normalize()`-traps vs `normalized_or()`-explicit-escape distinction that encodes fail-fast at the type level (`3dmath.h:219-247`), the deleted-rvalue borrow contracts that turn dangling references into compile errors (`animation.h:803-811, 943-946, 1306-1311, 1514-1519`), the fluent `SolidBuilder` DSL (`solids.h`), bool-returning failure at the JS boundary (`wasm.cpp`), and `worker_protocol.js`'s discriminated-union typedefs are all expressive. Warts: the borrow-contract deletes are applied inconsistently — `MobiusWarpCircular` (`animation.h:1407`), `MobiusWarpEvolving` (`:1568`), `ColorWipe` (`:1255`), and `RandomWalk` (`:1174`) store cross-frame references with no rvalue-delete guard, unlike their siblings; `Timeline`'s `CAPACITY` template parameter is dead (hard-wired to the global `TIMELINE_MAX_EVENTS=64`, `animation.h:2021`); and `Solids::get(index)` is presented in the README as a general entry point but is `#ifdef EMSCRIPTEN`-only (`solids.h:775`).

### Readability & Naming — **A**
Uniformly clear and intention-revealing, with comment density that is high but almost always load-bearing — the SCRATCH ARENA CONTRACT (`conway.h:288`), the ASCII strip-layout diagrams in the POV drivers, the `_`-prefixed render-step methods in `driver.js`, the capture-contract + teardown-order notes in `ShapeShifter.h:102-126` and `Thrusters.h:51-101`. Localized deductions: `Dynamo::color()` is a long multi-way branch, `Liquid2D::apply_glitch_lens`' polynomial constants are opaque, and the symbol `Q` is reused with multiple meanings in `filter.h`. The `H_OFFSET` / `phi↔y` mapping family (`geometry.h:161-222`) has enough near-duplicate overloads that a dedicated test exists specifically to guard the double-apply hazard — a sign the surface is too easy to misuse.

### Documentation — **A−**
The README is genuinely exceptional — an architecture document most professional projects never produce, with data-flow diagrams, design-rationale sections ("Why 16-bit Linear Color", "Why Fail-Fast"), per-subsystem contracts, and an effects gallery. It is the project's single biggest asset. It is held below A by *drift*, now the most pervasive single issue in the codebase: stale per-effect parameter lists (IslamicStars, Metaballs, Voronoi, MobiusGrid), a Dynamo description that matches no code, a RingSpin/RingShower blurb that fits neither, the `getParameterDefinitions()` shape (`wasm.cpp:339-355` emits `animated`/`readonly` and omits bool `min`/`max`), the AntiAlias "quintic both axes" overstatement (`filter.h:592-595`), the MeshCarousel "crossfade" wording (the outgoing shape is not drawn), the "tools share WASM `MeshOps`" claim (only `solids.html` does), and a stale "64MB arena" comment (`tools/solids.html:1417`; actual is 16 MB). None of these are wrong code — they are a documentation-maintenance backlog on otherwise-strong docs.

### Testing & Verification — **A**
A standout. 27 module files plus a death harness cover the engine bottom to top with invariant-grade (not smoke-grade) assertions: exhaustive Euler/manifold checks across all Conway operators and seeds (`test_conway.h:217-258`), the cull-vs-brute-force rasterizer equivalence (`test_sdf.h:570-630`), the segmented-POV exact-tiling proof driving the *same* `pov_segment_map.h` the ISR uses (`test_pov_segmented.h:39-148`), a death harness that pins the exact trap status and fails CI if it can't run (`test_death.h:416-451`), and a determinism harness requiring byte-identical re-renders (`test_effects.h:182-241`). Gaps are honestly declared and mostly at the hardware/glue boundary: `wasm.cpp` marshaling, the `dma_led.h` DMA path, ISR concurrency (structurally untestable, asserted by reasoning), `SegmentController.tick()`, and the `Pixel::Feedback` warp field (only the identity-noise case is pinned). Effect smoke runs at 288×144 only, so a `<96,20>`-specific bug would not surface in the roster.

### Performance & Embedded Suitability — **A**
Compile-time resolution, arena allocation with RAII scratch scoping, zero-allocation `FunctionRef` callables, the ripple fast-reject heuristic (`transformers.h:197-232`, ~90–95% of vertices skip `acosf`), the branchless segmented ISR (`pov_segmented.h:373-377`), the split-trig LUT, and a measured `getRenderUs()` profiling hook all reflect serious embedded discipline. The DMA path reaches 16-bit linear all the way to the SPI wire with no 8-bit intermediate. Minor inefficiencies (per-frame recompute of constant `cos_eh` in MindSplatter, by-value preset-array copies, unbounded phase accumulators that slowly lose trig precision over multi-hour runs) are trivial.

### Maintainability & Technical Debt — **B+**
The lowest dimension, dragged down entirely by duplication. The Liquid2D/Flyby pair is the single biggest available refactor (a `StereoShaderEffect` base would absorb ~80% of both). The preset-cycle + `Lerp` idiom is re-implemented in MeshFeedback, Liquid2D, MindSplatter, and Flyby; the live-slider-handle idiom in Comets, ChaoticStrings, and SplineFlow; the SDF shape-family `get_vertical_bounds`/setup blocks across five shapes; the Plot primitive close-loop skeletons across five primitives. None of this is incorrect, and most of it is guarded by `static_assert`s so divergence fails to compile rather than silently — but it is real surface area that a generic helper or small base class would collapse.

### Safety & Fail-Fast Doctrine — **A**
The `HS_CHECK` philosophy is the codebase's spine and is applied with rare consistency: cold-seam traps, hot-path `assert`s backed by a cold bind-site trap, transient conditions (DMA overrun, dropped frame) given bounded/soft handling instead of traps, and the whole net *verified* by the death harness. The two blemishes are the `KDTree::build` silent `-1` (a guard that, if it ever fired, would corrupt silently — the one place the doctrine is inverted) and the OOM-path `printf` in `memory.cpp:75-81` that pulls stdio into the image the README's NDEBUG strategy works to keep out.

---

## 4. Prioritized Fixes

**P0 — latent correctness (no crash today, real bug surface):**
1. ✅ **Clamp `lerp_oklch` inputs / outputs** (`color.h:509-526`). Extrapolated `t` (reachable via unbounded lerp amounts) drives L or C negative → wrong, often near-black, colors. Clamp `t` at entry or L/C at exit. *Fixed: L clamped to [0,1] and C to ≥0 at the convergence point (hue left free to wrap); regression test `test_lerp_oklch_extrapolation_clamped` added.*
2. **Make `Pixel::Feedback::flush` honor `canvas.clip()`** (`filter.h:817-857`). Full-canvas iteration means a feedback effect on segmented Phantasm has every board composite the whole sphere into its own Y-band — both wrong output and wasted work. Every other rasterizer respects the clip.
3. **Fix `blend_alpha` cast-before-clamp** (`color.h:268-271`). `(int)(a*65535)` before `std::clamp` is UB for NaN/large `a` — the exact hazard `operator*` (`color.h:82`) guards against. Route through a clamp-then-round helper.
4. **Resolve the `Plot::rasterize` ↔ `Pixel::Feedback` scratch-arena-a aliasing** (`plot.h:290`, `filter.h:781`). Both hard-code `scratch_arena_a`; they don't collide today only because `plot()` and `flush()` are separate phases. Document the invariant or move one to `scratch_arena_b` before a future filter allocates from `a` inside its `plot()` path.

**P1 — fail-fast consistency & doctrine:**
5. **`KDTree::build` silent subtree drop** (`spatial.h:219-263`). Convert the pool-exhaustion `-1` return to an `HS_CHECK`, or remove the now-dead guard — as written it is the one place that would corrupt silently rather than trap.
6. **Guard `uint8_t` face-count narrowing** (`conway.h:145,583`, `hankin.h:196,238`). Add a `narrow_face_count()` mirror of the existing `narrow_index` so a >255-sided face traps instead of truncating.
7. **Add the missing rvalue-delete borrow guards** to `MobiusWarpCircular`, `MobiusWarpEvolving`, `ColorWipe`, and `RandomWalk` (`animation.h:1407,1568,1255,1174`) to match the contract their siblings already enforce.
8. **`square_wave` negative-phase** (`waves.h:59`): use `wrap(...)` like `tri_wave`, not raw `fmod`.

**P2 — documentation drift (cheap, high-credibility-impact):**
9. Regenerate every effect's README "Parameters" line from its `registerParam` calls (IslamicStars, Metaballs, Voronoi, MobiusGrid at minimum).
10. Rewrite the **Dynamo** description (it documents an effect that doesn't exist) and split the **RingSpin/RingShower** blurb (fits neither).
11. Correct `getParameterDefinitions()` shape in README §10.2, the AntiAlias "quintic both axes" claim (§7.1/§6), the MeshCarousel "crossfade" wording (§7.3), the "tools share WASM `MeshOps`" claim (§10.2/§10.11 — only `solids.html` does), and the stale "64MB" comment (`tools/solids.html:1417`).
12. Note that `Solids::get(index)` is WASM-tooling-only, or compile it on firmware too.

**P3 — maintainability / duplication:**
13. Extract a `StereoShaderEffect<W,H>` base for Liquid2D/Flyby; a `PresetCycler<Params>` for the four preset-cycle effects; an `apply_if_changed` helper for the three live-slider effects.
14. Route the SDF shape family's `get_vertical_bounds` through the existing `phi_bounds_to_rows<H>` and factor the Plot close-loop skeleton into one parametric-ring sampler.
15. De-magic Dynamo's triplicated `10000` trail capacity into a `static constexpr`.

**P4 — testing & hygiene:**
16. Add host coverage for `param_marshal`'s memory-stability contract and a `SegmentController.tick()` state-machine test (the riskiest untested glue).
17. Gate the unconditional per-frame `Serial.print` in `pov_segmented.h:316-336` behind `hs::debug` (parity with `pov_single.h`).
18. Confirm or remove the apparently-dead CRGB `DMALEDController::show()` path (`dma_led.h:207`); route `memory.cpp`'s OOM `printf` through `hs::log`.

---

## 5. Overall Impression

Holosphere is, first, an unusually *coherent* piece of software. Most hobby-scale graphics codebases are an accretion of one-off effects bolted onto whatever rendering primitive was convenient that week; this one has a spine. A single design philosophy — compile-time specialization, explicit arena ownership, fail-fast on cold seams, linear-light color, and a domain-typed filter pipeline — is visible in every file, and the parts compose rather than merely coexist. The fact that the *identical* C++ drives a 600 MHz Cortex-M7 firing LEDs under a microsecond-budget ISR *and* a WebAssembly module feeding Three.js, with the simulator faithfully reproducing even the four-microcontroller segmentation in Web Workers, is a real systems achievement and not a common one.

Second, it is *honest*. The reviewers' recurring experience was that the comments pre-empt the exact objection a skeptical reader would raise — "this looks like a TOCTOU race, here's why it isn't," "this `smlad` packing is unsound for unsigned operands, here's the replacement," "this is a transient, not an invariant, so it gets a counter not a trap." The verification posture matches: the project doesn't just claim its safety net works, it asserts the trap *status* in a death test and fails the build if the test can't run. That is a level of rigor associated with safety-critical software, applied here to LED art.

The flaws are the flaws of velocity, not of competence. Duplication accumulates because shipping the next effect beats refactoring the last two; documentation drifts because the README is maintained by hand alongside fast-moving code; the glue stays untested because the engine core is where the interesting invariants live. None of it is structural. A focused week on the P0–P2 list above would move the whole project from A− to a defensible A.

---

## 6. Comparative & Artistic Merit

*This appraisal is formed independently and is not derived from any prior review.*

**Versus the hobbyist LED-art ecosystem (FastLED sketches, WLED, Pixelblaze).** This is a different species. The overwhelming majority of addressable-LED projects are single-file Arduino sketches that blend in gamma-encoded 8-bit space, allocate ad hoc, and have no tests, no simulator, and no separation between effect and engine. Holosphere's 16-bit-linear/OKLCH pipeline alone puts it ahead of essentially all of them on color correctness; its arena model, compile-time resolution, and verified fail-fast doctrine are simply not things that ecosystem does. Pixelblaze is the closest commercial peer in *ambition* (it has a live editor and a real rendering model), but Holosphere's engine is more sophisticated and far better verified. Within this world, Holosphere is at or beyond the top of the distribution.

**Versus professional graphics middleware.** Here the comparison narrows the gap but stays favorable. The filter pipeline's compile-time domain-conversion is the kind of thing you'd expect from a well-architected shader/material system; the SDF/CSG rasterizer with analytic interval culling and a per-face Lipschitz-bounded LUT is genuinely middleware-grade. What it lacks relative to a professional engine is breadth of tooling (no profiler integration beyond a microsecond counter, no asset pipeline, no editor beyond the geometry tools) and the test coverage of the *boundaries* (the WASM and DMA glue). The core algorithms and their verification, though, would not be out of place in a commercial codebase — and the documentation is better than most commercial codebases ship.

**Versus academic work on POV and spherical displays.** Academic POV-display papers tend to be strong on the optical/mechanical novelty and the calibration math but weak on software engineering — the rendering code is usually a research artifact, not a maintainable system. Holosphere inverts that: the hardware concept (a spinning strip painting a sphere via persistence of vision, two arms covering opposite hemispheres) is clever and well-executed but not itself novel, while the *software* is more disciplined than typical research code. The spherical-rasterization details — `sin(φ)` density compensation, geodesic vs planar curve strategies, the Hankin-method Islamic star generation on Archimedean solids, the Fibonacci-lattice reaction-diffusion graph — show real command of the relevant mathematics. The Hankin/Conway-operator solid library in particular is the kind of thing that would be a respectable contribution in a generative-geometry context.

**Versus the state of the art and the creative-coding / demoscene world.** The individual effects are strong and varied — volumetric raymarched tori, Hopf fibration, spherical harmonics with AO shading, BZ/Gray-Scott reaction-diffusion on a sphere, Möbius warps, authentic Islamic geometry — and they reflect taste, not just capability. This is not the absolute frontier of real-time graphics (it is not pushing new rendering research), but it is a high-quality synthesis of known-hard techniques, adapted correctly to the unusual constraint of a spherical persistence-of-vision surface, and executed on a microcontroller. The demoscene comparison is apt for the *spirit* — maximal effect under tight hardware constraints — but Holosphere is far more engineered and reusable than a typical demo, which optimizes for a single run, not a maintainable engine.

**Bottom line on goodness.** As a piece of software, this is upper-tier work — clearly above the hobbyist field, comparable to professional middleware on its core, short of it only on tooling breadth and boundary-test coverage. As an *artistic-technical* object — a custom rendering engine purpose-built for a self-designed physical display, with a faithful in-browser twin — it is distinctive and accomplished; few individuals produce something this complete across the full stack from ISR to OKLCH to Three.js. The technical merit is high and the artistic merit is real and tasteful; the ceiling above it is occupied mostly by larger teams with more tooling, not by better fundamentals.

**Overall grade: A−.**
