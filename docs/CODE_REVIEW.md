# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`, `targets/wasm/`, `tests/`) and the daydream web simulator (`*.js`, `tools/`, `index.html`). Out of scope per request: `effects_legacy.h`, the `*.ino` target sketches, and `rotate.h`.

**Method:** the architecture was read from the README, then thirteen independent sub-agents reviewed each component in depth against that architecture. Every substantive finding was re-checked by a separate validation agent reading the actual source; refuted claims were dropped and over-stated severities were downgraded. Verdicts (CONFIRMED / PARTIAL / REFUTED) are reflected below.

This is, in aggregate, an exceptionally well-engineered codebase. The findings are dominated by Low-severity consistency, documentation-drift, and missed-optimization items; only a small number of genuine runtime defects survived validation. Grades are calibrated against the codebase's own very high bar, not against a hobby-project baseline.

---

## Summary of Grades

| Dimension | Grade |
|---|---|
| Architecture & Design Elegance | A+ |
| Interface Expressiveness / API Design | A |
| Correctness & Reliability | A− |
| Memory Safety & Resource Management | A |
| Concurrency & Real-Time Discipline | A− |
| Performance & Efficiency | A− |
| Testing & Verification | A |
| Documentation | A− |
| Code Style & Consistency | A− |
| Maintainability | A− |
| Hardware / Embedded Discipline | A |
| **Overall** | **A− / A** |

---

## Dimension-by-Dimension Rationale

### Architecture & Design Elegance — A+
The structural decisions are coherent and mutually reinforcing: compile-time `<W,H>` resolution specialization eliminates runtime generality cost; a single partitioned arena (persistent + scratch A/B) gives deterministic memory behavior with explicit `Arena&` plumbing at every call site (no hidden state); the rendering pipeline cleanly separates Generate → Transform → Rasterize → Filter; and the variadic `Pipeline<W,H,Filters...>` expresses three coordinate domains (World/Screen/Pixel) with compile-time domain-transition routing and `static_assert`-enforced ordering. The Phantasm 1-wire sync protocol is the standout: deriving column phase from a free-running cycle counter so masked-IRQ windows cannot drop columns, folding the epoch forward to make the 32-bit counter wrap structurally unobservable, and an odd-only distance-2 symbol alphabet that degrades any single glitch to a *missed* (never *misclassified*) symbol. This is genuinely fault-tolerant systems design.

### Interface Expressiveness / API Design — A
APIs are hard to misuse and read as intentional. `FunctionRef` vs `StoredFunctionRef`, `ArenaVector` (owning, move-only) vs `ArenaSpan` (explicit borrow), the `Fragment`/shader-register model shared across all rasterizers, named factories over ambiguous constructors, and the `SolidBuilder` fluent Conway chain are all expressive and self-documenting. The Feedback `Style` POD-with-presets design removes a whole class of template/adapter boilerplate. Minor expressiveness leaks exist (the `Fragment::aux`/`v3` register is universally zero across built-in 2D shapes yet documented as a live channel; `Ring`'s per-plane normal ctor parameter is always passed `Y_AXIS`).

### Correctness & Reliability — A−
Most subsystems are defect-free under scrutiny — the math, distance/interval generation, mesh topology, arena polarity, color round-trips, and sync codec all validated clean. Two genuine runtime defects pull this off an A: an attractor `kill_radius` that does not actually kill a particle (it resurrects the next frame), and a `Rotation` animation constructor missing the `duration >= 0` guard every sibling animation carries. A handful of Medium issues in the segmented-render path of the simulator (a `composite()` that hard-throws where its sibling self-heals; animated GUI sliders that freeze) round out the list.

### Memory Safety & Resource Management — A
Arena discipline is rigorous: RAII `ScratchScope`/`Persist`, LIFO reclaim contracts honored, compile-time `static_assert`s sizing scratch against the real device budget, and always-on `HS_CHECK` traps at every allocation/capacity/registration seam. `ArenaVector` carries debug generation-tracking for use-after-free. The death harness proves the traps actually fire. No leaks, double-frees, or aliasing defects survived validation.

### Concurrency & Real-Time Discipline — A−
The double-buffer atomics, the single-writer sync-state ownership model, and every non-obvious memory-order choice are documented with specific, correct rationale. The one gap: the single-board `pov_single.h` driver shares `effect_`/`x_` between the IntervalTimer ISR and the main loop as plain statics (not `std::atomic`/`volatile`), even though the segmented driver and `Canvas` both deliberately use `std::atomic` for the identical handoff. Low probability on a single aligned-word core, but a formal data race and an inconsistency with the codebase's own published pattern.

### Performance & Efficiency — A−
Hot paths use fast-trig approximations, baked palette LUTs, branchless ISRs, and a 16-bit-linear pipeline that reaches the SPI wire with no 8-bit intermediate. A few avoidable costs remain: `DistortedRing` evaluates a `GenerativePalette` (full OKLCH lerp) per pixel instead of baking like its siblings; planar arc length is computed twice per frame in `plot.h`; and a few shape constructors use exact `acosf` where the rest of the file uses `fast_acos`.

### Testing & Verification — A
Oracle-driven assertions throughout (independent double-precision/analytic/brute-force references rather than re-reading set values), a real out-of-process death harness asserting `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` for 35 invariant traps, standalone TUs that recompile asserts under the *shipping* fast-math flags and the *device* `H_OFFSET` geometry, and a runtime WASM smoke gate over every effect × resolution. Deductions: `dma_led.h` and `pov_single_map.h` have no host-testable seam extracted; the committed `color_luts.h` is not pinned against a fresh regeneration; native perf/arena budgets are measured but not gated.

### Documentation — A−
The README is extraordinary — 2142 lines including a datasheet-grade sync-protocol spec with AC timing tables and signal waveforms — and the in-code comments earn their place by explaining *why* rather than restating code. The grade is held back only by drift: a phantom "Ripp Width" parameter, a direct self-contradiction about `sin(φ)` density compensation, an undocumented `SlowTwist` feedback preset, a "bilinear" label for what is actually a quintic-eased splat, and a stale "cubemap face" doxygen on `MobiusGrid`.

### Code Style & Consistency — A−
Comment discipline is terse and fact-focused; naming and idiom are uniform. Minor inconsistencies: `std::clamp` vs the NaN-safe `hs::clamp` in a few color paths, a handful of `HS_CHECK` calls missing their message argument, member/initializer ordering that would trip `-Wreorder`, and an inconsistent "exceeds 16-bit range" bound convention between two guards.

### Maintainability — A−
Code is well-factored with shared scaffolds (CRTP reaction-diffusion base, half-edge operator helpers, pooled label sprites). The residual risks are a few hand-synced duplications: `Plot::Mesh::draw` vs `extract_edges` repeat the face-walk/edge-dedup loop verbatim, and the Möbius complex-arithmetic core is copied between JS and GLSL with only a comment to keep them aligned.

### Hardware / Embedded Discipline — A
The fail-fast philosophy is applied exactly as prescribed (cold-path `HS_CHECK` that survives `NDEBUG`; stripped `assert` on hot paths backed by cold bind-site traps; transient conditions like DMA overrun get bounded/soft handling, not traps). DMAMEM placement, branchless ISR, and clean-no-invalidate cache flush direction for the TX-only buffer are all correct.

---

## Prioritized Fix List

Items are numbered sequentially across all priority tiers. Each is independently validated unless noted.

### Priority 1 — Correctness defects (runtime behavior)

1. ✅ **Attractor `kill_radius` never kills a particle (resurrection bug).** `core/animation.h` ~658–726 (`step_particle`): the kill branch `if (dist_sq < attr.kill_radius * attr.kill_radius) { active = false; break; }` sets the local flag but never zeroes `p.life`. Next frame `active = p.life > 0` is true again, so the particle re-runs full physics and resumes recording its trail, fully resurrecting if it drifts back out of the kill radius before `life` naturally expires. Fix: set `p.life = 0;` in the kill branch.

2. ✅ **`Rotation` constructor missing `HS_CHECK(duration >= 0)`.** `core/animation.h` ~1450–1455: every sibling animation (`Transition`, `Mutation`, `Motion`, `MobiusFlow`/`Warp`/`WarpCircular`, `Ripple`) rejects the perpetual `-1` the base permits; `Rotation` does not, yet `step()` computes `easing_fn(t / duration) * total_angle`, so a `-1` duration divides by `-1` and feeds the easing function a negative, growing argument — a garbage, never-completing rotation. Fix: add the same guard to the parameterized ctor.

3. ✅ **Segmented `composite()` hard-throws where the single-engine path self-heals.** `daydream/segment_controller.js` ~470–483 throws if `getMemoryView() !== Daydream.pixels`, but only the single-engine branch of the adapter `drawFrame` (`daydream/daydream.js` ~490–508) re-points the alias; the segmented branch never reconciles it. A resolution-change race can throw inside the animation loop and kill rendering. Fix: re-point the aliases in `composite()` (mirroring the single-engine branch) instead of throwing.

4. ✅ **Animated GUI sliders freeze in segmented mode.** `daydream/daydream.js` ~486–510: the segmented `drawFrame` path calls only `segments.tick()` and never steps `host.engine`, yet `syncGUI()` reads `host.engine.getParamValues()` — so animation-driven sliders read static values from the un-stepped main engine while the worker-rendered sphere animates. Fix: source GUI param values from segment 0's worker in segmented mode, or annotate the controls as worker-driven.

5. ✅ **`setEffect` failure path leaves stale strobe-column mode.** `daydream/daydream.js` ~234–241: on `host.engine.setEffect(...) === false` the function returns early without calling `daydream.setStrobeColumns(...)`, so the previous effect's column-fill mode persists. Fix: set the strobe/round-dot mode on the failure branch too.

### Priority 2 — Concurrency, latent safety, and documentation defects

6. ❌ **Single-board ISR handoff is a formal data race.** `hardware/pov_single.h` ~205–208: `static Effect *effect_;` and `static int x_;` are shared between the IntervalTimer ISR (`show_col`) and the main loop (`run`) with no `std::atomic`/`volatile`/barrier, unlike `pov_segmented.h` (~601, `std::atomic<Effect*>` with release/acquire) and `Canvas` (`std::atomic<int> prev_`). Fix: make `effect_`/`x_` `std::atomic` to match the codebase's own published model. — **INVALID (validated):** no concurrent access exists. `effect_`/`x_` are written only before `timer.begin()` and (for `effect_`) after `timer.end()`, so the ISR never runs during a foreground write; `x_` is ISR-exclusive after its init. `timer.begin()` is an opaque-call compiler barrier and the PIT ISR shares the single Cortex-M7 core, so publication is correctly ordered without atomics. The cited asymmetry is backwards: segmented's ISR-owned analog (`live_effect_`, ~601) is itself a plain `Effect *`; its atomics serve a live handoff to a *running* ISR, which single-board does not do.

7. ✅ **README self-contradiction on AntiAlias `sin(φ)` density compensation.** `README.md` line ~448 says `Screen::AntiAlias` distributes "with `sin(φ)` density compensation"; line ~508 says "no `sin(φ)` density compensation". The code (`core/filter.h` ~875–921) applies no latitude term, so line 448 is wrong. Fix: delete the density-compensation clause from line 448.

8. **README advertises a non-existent IslamicStars parameter.** `README.md` §9 lists IslamicStars params including "Ripp Width", but `effects/IslamicStars.h` `init()` registers no such slider — ripple thickness is the fixed constant `kRippleThickness` ("thickness is fixed (not a slider)"). Fix: remove "Ripp Width" from the README parameter list.

9. **`vector_to_pixel` returns an unclamped south-pole `y` (contract hazard).** `core/geometry.h` ~485–491: `y` can land "a hair above `H_VIRT-1`" at the south pole and is deliberately left unclamped, with a documented "caller must clamp" contract. The primary consumer (`AntiAlias::plot`) already bounds-checks each tap, so there is no active OOB today — but every other `(int)y` row-buffer consumer inherits the hazard. Fix (defensive): clamp `y` once at the producer, or keep the contract but add a guard at any bare-index consumer.

10. **Disagreeing degeneracy thresholds feed a zero tangent into `screen_step`.** `core/plot.h`: `rasterize_geodesic_strategy` (~317) treats a segment as degenerate below `EPS_GEODESIC_SEGMENT` (1e-3) and installs a zero-tangent sampler, while `process_segment` (~581) only short-circuits below `EPS_GEOMETRIC` (1e-5); a length between the two passes the draw guard yet carries a zero tangent. Benign today (one dot) but the thresholds should be reconciled.

11. **Undocumented `SlowTwist` feedback preset.** `core/styles.h` ~184 defines `Style::SlowTwist()` (the preset anchoring the field-order `static_assert`), but the README §6 preset table omits it. Fix: document it or note it is internal.

12. **`AntiAlias` weights mislabeled "bilinear".** `core/filter.h` ~875–921 uses `quintic_kernel` (smootherstep) for a 2-tap partition, not true bilinear weights; README §6 (lines ~448/476/508) and the code comment call it "bilinear". Energy is conserved (partition of unity), so this is terminology, not a bug. Fix: describe it as a "quintic-eased 2×2 splat".

### Priority 3 — Performance & efficiency

13. **`DistortedRing` evaluates a generative palette per pixel.** `effects/DistortedRing.h` ~101–110 calls `ringPalette.get(f.v0)` (full OKLCH lerp) in the fragment shader, though `ringPalette` is an immutable `GenerativePalette`. Siblings (`Moire`, `Dynamo`) bake a `BakedPalette` LUT specifically to avoid this. Fix: bake once in `init()` and sample the LUT.

14. **Planar arc length computed twice per frame.** `core/plot.h` ~542–553 / ~733–738: the cold pre-pass and the draw loop both call `planar_arc_length` per segment for Star/Flower/PlanarPolygon every frame. Fix: cache per-segment lengths in the pre-pass.

15. **Exact `acosf` on per-shape construction paths.** `core/sdf.h` (multiple shape ctors, e.g. ~514/712/2276): `acosf(clamp(ny))` is used for a bounds estimate where the rest of the file uses `fast_acos`. Minor, but inconsistent with the file's own fast-trig discipline.

### Priority 4 — Maintainability, consistency & dead code

16. **`Plot::Mesh::draw` and `extract_edges` duplicate the face-walk/edge-dedup loop verbatim.** `core/plot.h` ~1980–2053, differing only in the per-edge action. Fix: extract a shared `for_each_unique_edge(mesh, fn)`.

17. **Dead branch in `compute_inradius`.** `core/sdf.h` ~1776: `(min_edge_dist > 1e8f) ? 1.0f : min_edge_dist` can never take the true branch because the `count > 0` loop always overwrites `min_edge_dist`. Fix: remove or document the dead guard.

18. **`Ring::sample` lacks the divisor guard its siblings have.** `core/plot.h` ~1582 area: divides by `num_samples` with no `HS_CHECK(num_samples > 0)` (the `Spiral` sampler guards with `HS_CHECK(n >= 1)`); a zero yields NaN fragments. Fix: add the guard for consistency.

19. **`MobiusGrid` doxygen describes the wrong geometry.** `effects/MobiusGrid.h` ~12–13: `@tparam W/H` documented as "Cubemap face width/height" (copy-paste); the effect renders a lat-long grid. Fix: reword to canvas width/height.

20. **`Voronoi` struct field order diverges from `registerParam`/README order.** `effects/Voronoi.h` ~238–242 (`num_sites, speed, borderThickness, sharpness`) vs ~47–51 (`Num Sites, Speed, Sharpness, Border Thick`). Harmless with named-member init, but a hazard for any future aggregate-init or preset table. Fix: align the orders.

21. **`RingSpin` spawns all four rings with an identical `Y_AXIS` normal.** `effects/RingSpin.h` ~95–98: the per-plane normal ctor parameter is always `Y_AXIS`, so rings differ only by palette/noise; the unused parameter reads as dead intent. Fix: pass varied normals or drop the parameter.

22. **`Thrusters` advances `t_global` before constructing `Canvas`** and uses a mixed declaration/init member order. `effects/Thrusters.h` ~78/80 and ~34–41/292–299. Behaviorally inert; breaks the "Canvas first" shape used everywhere else and would trip `-Wreorder`. Fix: reorder.

23. **`MindSplatter` per-pixel `pool[p_idx]` index guarded only by a stripped `assert`.** `effects/MindSplatter.h` ~273–274: the hot-path index has no cold bind-site `HS_CHECK` backstop, so a corrupt `f.v2` reads OOB silently in the device build. Fix: add a cold capacity/`active_count` invariant trap at the draw-call seam.

24. **`HS_CHECK` calls missing their message argument.** e.g. `core/hankin.h` ~241/246 and several orbit guards in `core/conway.h`; on a device trap these surface only file:line, losing the self-explanatory text the codebase otherwise standardizes on. Fix: add messages.

25. **`std::clamp` used where the file relies on NaN-safe `hs::clamp`.** `core/color.h` ~1289/1692/814: `std::clamp(NaN,...)` diverges from the documented `hs::clamp` NaN→hi contract the file depends on elsewhere. Fix: use `hs::clamp`.

26. **Redundant double `std::fmod` in `shortest_distance`.** `core/util.h` ~115: two `fmod` calls (each ~20+ soft-float cycles on M7) where one reduction plus a branchless fold suffices on the common in-range input. Minor; non-hot path.

27. **Möbius complex-arithmetic core duplicated between JS and GLSL.** `daydream/tools/mobius.html` ~345–354 vs `mobius_transforms.js` ~44–51: `cmult`/`cadd`/`cdiv` (including the `1e-6` divide guard) are hand-copied and can silently diverge. Fix: a parity unit test, or generate the GLSL from the shared spec.

28. **`recorder.js` mkv path produces a mismatched blob MIME / picker filter.** `daydream/recorder.js` ~283–335: `_extension` maps `video/x-matroska → 'mkv'` but `_download` sets the blob type to `video/webm` for any non-mp4 ext and builds an internally inconsistent `accept` map. Fix: derive the blob type from an ext→MIME map.

29. **`gui.js` `_attachUrlWriter` makes `onChange` non-idempotent.** `daydream/gui.js` ~172–181: reassigns `controller.onChange` to a single-slot setter, so a second `.onChange` registration silently replaces the first and the load-replay only fires for the original handler. Fix: support an array of user handlers or document the single-registration contract.

30. **`effectReady` worker message is received and ignored.** `daydream/segment_controller.js` ~204: the worker emits `effectReady` after every `setEffect` but the controller does nothing with it, so a frame can composite mid-switch; the message reads as if synchronization exists when it does not. Fix: gate compositing on it (mirroring `ready`) or drop it from the protocol.

31. **`URLSync` applies the `effect` URL key with no validator.** `daydream/state.js` ~159–167: a garbage/empty `?effect=` is written straight into state during hydration, overriding the validated default (it self-heals later via `resolveActiveEffect`, but briefly holds an invalid value). Fix: pass an `effect` validator alongside the `resolution` one.

32. **No runtime guard at the worker postMessage trust boundary.** `daydream/segment_worker.js` ~201 / `worker_protocol.js`: inbound `e.data` is cast to the protocol union with an unchecked JSDoc assertion and the `switch` has no `default`, so a malformed/drifted message no-ops silently. Fix: add a `default` arm that logs unknown message types.

### Priority 5 — Testing & verification gaps

33. **`dma_led.h` and `pov_single_map.h` have no host tests.** The project's own pattern (extracting `pov_segment_map.h`/`hd107s_frame.h` as host-testable pure logic) shows the index/timing math can be lifted; currently any pure-arithmetic seam in these two ships unverified. Fix: extract a host-testable index unit for `pov_single_map.h`.

34. **Committed `color_luts.h` is not pinned against its generator.** `scripts/generate_luts.py` is the "generator of record," but no test asserts the committed table matches a fresh regeneration, so a hand-edit or stale divergence would pass. Fix: a CI check that regenerating yields no diff.

35. **Native perf and per-effect arena budgets are measured, not gated.** `tests/perf_bench.cpp` / `arena_measure.cpp` are informational (not CTests), so a native perf regression fails no build. Partially backstopped by device-byte-budget regressions, the WASM HWM gate, and the gated `stack_measure`. Fix (optional): add coarse perf-regression thresholds where stable.

---

## Closing Note

The validated defect density here is remarkably low for a codebase of this ambition and surface area. The two Priority-1 engine bugs are notable precisely because they are *omissions* — a missing `life = 0` and a missing duration guard — in code that is otherwise guarded with unusual rigor, and the segmented-simulator issues are confined to one render path. Addressing Priority 1 and 2 would put correctness on par with the codebase's already A+ architecture.
