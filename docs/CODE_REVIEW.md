# Holosphere / daydream — Code Quality Review

**Scope.** The C++ rendering engine (`core/`, `effects/`, `hardware/`, `targets/`), the
native test suite (`tests/`), and the daydream web simulator (`*.js`, `tools/`). Out of
scope by request: `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/rotate.h`, and all vendored/build trees (`.pio/`, `node_modules/`, `three.js/`).

**Method.** Eighteen component-scoped reviewer agents read every in-scope file in full
against the README architecture; each candidate finding was then handed to a separate
agent that re-opened the cited code and independently confirmed, rejected, or
re-graded it. 62 candidate findings were raised; **7 were rejected in validation** and
**55 survived** and appear below. No finding is included that a second reader could not
substantiate against the source.

**Headline.** This is a mature, professionally-engineered codebase of unusually high
quality. Across ~30k lines of C++ and ~6k lines of JavaScript the validators found **no
critical or high-severity defect** — no memory-safety bug, no concurrency/ISR race, no
incorrect core algorithm. The 55 findings are dominated by documentation drift,
consistency/DRY nits, micro-performance idioms, and a small number of latent
robustness gaps that are unreachable on the shipped configuration. The engineering
discipline — fail-fast `HS_CHECK` at cold seams, compile-time resolution, arena
allocation with RAII scoping, single-source X-macro rosters, and a self-verifying death
harness — is applied consistently and is the dominant characteristic of the code.

---

## Overall Letter Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness & reliability** | **A** | No substantiated correctness bug in any core algorithm. Singular/degenerate cases (poles, antipodes, stereographic singularities, degenerate edges, NaN folds) are each handled and round-trip-tested. The handful of correctness-tagged findings are latent (unreachable on the shipped roster) or documentation-precision issues. |
| **Memory safety** | **A** | Trap-on-overflow arena allocator with wrap-proof bounds math, fixed-capacity non-allocating containers, debug generation/rebind stamps for use-after-free, and `DMAMEM`-segregated DMA buffers. Index casts route through trapping `narrow_index`. No OOB or lifetime hazard found. |
| **Concurrency & ISR safety** | **A** | The single-core ISR/main-loop double-buffer is correct under relaxed atomics, with the load-bearing `disable_interrupts()` compiler barrier and its multi-core caveat documented at the site. The Phantasm flywheel (single-writer, release/acquire handoff) and per-worker WASM isolation are equally rigorous. |
| **Performance & efficiency** | **A** | Split trig LUTs, branchless hot-path clamp/lerp/wrap, `fast_*` approximations with measured error budgets reserved for per-pixel paths, baked-palette LUTs, analytic scanline culling, and zero per-frame allocation on steady-state paths. |
| **API / interface expressiveness** | **A−** | Borrow-vs-store encoded in the type system (`FunctionRef`/`StoredFunctionRef`, `ArenaVector`/`ArenaSpan`), `explicit` constructors that block coordinate-space confusion, const-only outward `ParamList` handles single-sourcing the write gate. Minor inconsistencies in encapsulation conventions across the filter family. |
| **Architectural elegance** | **A** | Compile-time `<W,H>` specialization, the variadic filter `Pipeline` with automatic domain-transition dispatch, CRTP scaffolds that eliminate duplication at zero cost, and a clean split of pure host-testable logic from device shells. A coherent, deliberate design. |
| **Readability & maintainability** | **A−** | Doc comments explain *why*, not just *what*, and are accurate against the code. The one drag — flagged in several components — is comment density so high (often exceeding the code) that straight-through reading is slow and the comment surface invites drift. |
| **Error handling** | **A** | The fail-fast doctrine is applied with discipline: always-on `HS_CHECK` traps at cold setup/allocation seams, debug-only `assert` on per-pixel hot paths, reject-and-continue at the untrusted JS boundary, bounded/soft handling for genuinely transient conditions (DMA overrun, dropped frame). |
| **Portability** | **A−** | x86/Cortex-M7/WASM clamp backends share a documented NaN contract; sim/device bit-identity is deliberately preserved; 32/64-bit deltas are reasoned about. Deductions for a few standards-conformance nits (signed-int frame counters, `size_t` pointer math) and one unguarded `RD_N`-derived constant pair. |
| **Testing & verification** | **A−** | A standout death harness that requires the *exact* illegal-instruction trap status, a determinism pass with deliberate inter-run state perturbation, white-box numeric seams, a multi-board sync simulator, and runtime WASM smoke in CI. Deductions for thin assertion density in the dedicated memory tests and two hand-maintained roster lists. |
| **DRY / reuse** | **A** | Shared scaffolds (`least_parallel_axis`, `STEREO_INF`, `shade_mesh_topology`, half-edge orbit/edge helpers, CRTP iterators) eliminate essentially all copy-paste; the few duplications found are trivial and localized. |
| **Documentation accuracy** | **B+** | The weakest sub-area, and the source of most "fixable" findings: several README/code contradictions (AntiAlias `sin(φ)`, the arena-budget figures and a worked example that would trap on device, the `relax()` 500-vs-1000 cap, a beacon span-budget arithmetic) where the *code* is right and the *prose* is stale. |
| **Build / tooling / CI** | **A** | CMake presets, a three-layer test gate (pre-commit hook → presubmit CI on Linux+Windows → gated Pages deploy), runtime WASM smoke, install-provenance verification, and an optional headless Teensy size/layout gate. |

**Overall: A.** A codebase whose defect profile after deep multi-agent review is "tighten
the prose and a few latent guards," not "fix the bugs," is operating at a level well
above typical hobbyist or even much professional embedded/graphics work.

---

## Per-Component Summary

| Component | Files | Standout grades | Findings |
|---|---|---|---|
| Math / geometry | `3dmath.h`, `geometry.h`, `util.h`, `easing.h`, `waves.h` | Correctness A, Readability A | 3 |
| Color | `color.h`, `palettes.h`, `styles.h`, `color_luts.h` | Memory-safety A+, Elegance A+ | 2 |
| Render core | `canvas.h`, `concepts.h`, `engine.h`, `platform.h`, `led.h` | Memory-safety A, Concurrency A | 3 |
| Filter pipeline | `filter.h` | Performance A, Elegance A | 5 |
| SDF / scan | `sdf.h`, `scan.h` | Memory-safety A | 4 |
| Plot / transformers | `plot.h`, `transformers.h` | Correctness A, Memory-safety A | 3 |
| Animation | `animation.h`, `presets.h`, `generators.h` | Memory-safety A, Elegance A | 3 |
| Memory | `memory.h/.cpp`, `static_circular_buffer.h`, `inplace_function.h` | Correctness A, Error-handling A | 3 |
| Mesh | `mesh.h`, `conway.h`, `solids.h` | Memory-safety A, DRY A | 4 |
| Mesh (advanced) | `hankin.h`, `spatial.h`, `reaction_graph.*` | Correctness A, Memory-safety A | 3 |
| Effects A | RD base, BZ/GS, ChaoticStrings, Comets, DistortedRing, DreamBalls, Dynamo, FlowField | Correctness A, Arena A | 3 |
| Effects B | Flyby, GnomonicStars, HankinSolids, Hopf, IslamicStars, Liquid2D, MeshFeedback, MindSplatter, MobiusGrid, Moire | Correctness A, Pipeline A | 1 |
| Effects C | PetalFlow, Raymarch, RingShower, RingSpin, ShapeShifter, SphericalHarmonics, SplineFlow, Thrusters, Voronoi | Correctness A, Memory-safety A | 2 |
| Hardware | `dma_led.h`, `hd107s_frame.h`, `pov_segment_map.h`, `pov_single*.h`, `pov_segmented.h`, `pov_sync.h` | Concurrency A, Architecture A | 2 |
| Targets | `wasm.cpp`, `param_marshal.h`, `Phantasm.ino` | Memory-safety A, Elegance A | 2 |
| C++ tests | `run_tests.cpp`, death/effects/memory/sync harnesses | Testing A, Architecture A | 5 |
| daydream core | `daydream.js`, `driver.js`, `state.js`, `gui.js`, `recorder.js`, … | Memory-safety A, Testing A | 5 |
| daydream segmented | `segment_controller.js`, `segment_worker.js`, `segment_layout.js`, tools | Correctness A, Concurrency A | 2 |

---

## Prioritized Fix List

Every validated finding appears below, numbered sequentially. Priority reflects
behavioral/robustness/clarity impact, **not** raw severity — all findings are individually
low-to-medium severity. Each item cites `file:line`, the dimension, and the recommended fix.

### High Priority — real behavioral, UB, or actively-misleading issues

1. ✅ **README contradicts the code on AntiAlias `sin(φ)` scaling** — `core/filter.h:1031-1039` vs `README.md:508`. *(documentation-accuracy, medium severity)* The README states `Screen::AntiAlias` "scales the X fractional by `sin(φ)`"; the code deliberately does **not** (an in-code comment documents that this density compensation was removed because it collapsed horizontal AA near the poles). Update `README.md:508` to describe uniform bilinear quintic-eased weights in framebuffer space with no `sin(φ)` compensation.

2. ✅ **README arena-budget figures are stale and the worked example would trap on device** — `core/memory.h:37,44-47` vs `README §7.5`. The README documents 335 KB total / 303 KB persistent and a `configure_arenas(271K, 32K, 32K)` example summing to 335 KB, but the real device block is 330 KiB / 298 KiB. A reader copying the example verbatim over-subscribes by one 4 KiB page and trips the `configure_arenas` budget trap (`HS_CHECK`) at `init()`. Update the README figures and example to the current 330 KiB / 298 KiB; state current values only (no rename history).

3. ✅ **`Thrusters::t_global` is a raw signed `int` that overflows (UB) on a long-running piece** — `effects/Thrusters.h:330` (incremented at `:80`, read as `frame % 32`). On an explicitly long-running art installation a signed `int` at ~60 fps is UB after ~414 days. Wrap it (`t_global = (t_global + 1) % 32;`) or make it `uint32_t`, matching Raymarch's deliberate switch to wrapped accumulators.

4. ✅ **`ShapeShifter::frame_count_` is a raw signed `int` that overflows (UB)** — `effects/ShapeShifter.h:297` (incremented at `:75`, read as `% 48`). Same class as #3 on the same long-running target. Make it unsigned or wrap at the cycle period.

5. ✅ **Segmented recorder captures stale/black frames after a worker fault** — `c:/work/daydream/segment_controller.js:820-823`. The faulted-tick early return never resets `frameComposited`, which can remain `true` from the last good composite; the recorder gate then keeps treating the frozen black buffer as real content. Set `this.frameComposited = false` before the faulted return, mirroring the normal path.

6. ✅ **`Rotation::step` dereferences `orientation` with no null guard** — `core/animation.h:1560-1605`. The default ctor sets `orientation = nullptr` and `has_orientation()` exists for exactly this state but has zero callers; a default-constructed `Rotation` that is stepped is a null deref under `NDEBUG`. Latent (every shipped effect uses the bound ctor). Add `HS_CHECK(orientation)` at the top of `step()` (the fail-fast trap, not a silent early return).

7. ✅ **Recorder `getContext('2d')` return is dereferenced unchecked → silently blank recording** — `c:/work/daydream/recorder.js:275,302`. A null 2D context is tolerated by the blit guard but `start()` still creates the `MediaRecorder`, producing a permanently blank file with no diagnostic — the file's only silent-failure path. On null, `console.error` and abort `start()`.

8. ✅ **`DistortedRing` silently culls geometry when `max_distortion` underestimates `shift_fn`'s true range** — `core/sdf.h:640-733`. `max_distortion` is a hard correctness precondition (it widens the row/column/per-pixel reject band) but is documented merely as "used to widen the band," and unlike every sibling shape it carries no guard. Underestimating it drops genuine arcs with no diagnostic. (All shipped callers pass an exact bound, so this is latent.) Document it as a precondition; optionally widen the early-reject by a safety epsilon.

### Medium Priority — robustness/consistency gaps and misleading documentation

9. ✅ **`set_clip`/`set_clip_x` validate nothing, unlike the eager `set_margin` trap** — `core/canvas.h:137-156`. An inverted band (`x1<x0`/`y1<y0`) is accepted and silently yields an empty/malformed render region with no breadcrumb at the configuration site. Add a cheap `HS_CHECK(x0<=x1 && y0<=y1 && x1<=clip.w && y1<=clip.h)`, or document the deferred-validation intent. *Fixed: guard the ordering + non-negative origins in both setters; the upper bound (`x1<=w`) is intentionally left to the downstream `Scan::Shader::draw` LUT-domain trap, which a death test pins.*

10. ✅ **`Conway::dual()` uses strict `normalized()` where every other operator uses `normalized_or()`** — `core/conway.h:435`. A zero-length face centroid (centrally-symmetric face) would trap the build; `snub()` already guards the analogous case with an EPS fallback. Unreachable on the fixed 52-solid roster. Use `normalized_or(c, fallback)` or document provable unreachability. *Fixed: `normalized_or(c, first_face_vertex)` — the first vertex is already on the unit sphere, so the dual degrades gracefully like every sibling operator.*

11. ✅ **Hankin rosette emission silently drops `count < 6` windings instead of trapping** — `core/hankin.h:280-285`. A degree-<3 orbit is silently skipped, diverging from the file's own fail-fast convention (cf. the unpaired-half-edge `HS_CHECK` in the same loop). Unreachable on a closed convex solid. Replace the silent skip with `HS_CHECK(count >= 6, …)`.

12. ✅ **`srgb_to_linear_interp` low guard does not catch NaN, unlike its sibling LUT helpers** — `core/color.h:356-369`. A NaN input reaches `static_cast<int>(f)` (float→int UB) — the exact hazard `Gradient::get`/`BakedPalette::get` were hardened against via `hs::clamp`. Both current callers satisfy the [0,1] contract, so latent. Route the input through `hs::clamp` at entry (NaN→hi).

13. ✅ **`Driver` fixed-speed ctor accepts a non-finite speed the bound-source path explicitly rejects** — `core/animation.h:1006-1009`. A NaN/Inf `speed` poisons `mutant` permanently via `wrap_t(NaN)`, the very hazard the bound-source ctor/`step()` guard with `std::isfinite`. Add a matching `HS_CHECK`/`isfinite` guard at construction.

14. ✅ **Feedback divisibility trap has no compile-time backstop for `downsample`** — `core/filter.h:1390-1396`. The per-flush `HS_CHECK(W % ds == 0 && H % ds == 0)` is the sole guard; a preset author setting a non-dividing `downsample` is caught only at runtime. (The live-lerp concern is refuted — `downsample` is snapped, not interpolated.) Add a `static_assert` alongside the other preset-range asserts. *Fixed: `static_assert` in `Feedback<W,H>` that the DEFAULT style downsample is `>0` and divides `W`×`H` (every shipped preset uses the default). The general check stays the runtime `HS_CHECK` because `downsample` is a runtime field — snapped/copied between presets and settable after the resolution is fixed — so divisibility against `W`/`H` can't be made compile-time in the general case.*

15. ✅ **`checkStaleTransfer` watchdog only runs on the overrun path; a wedged channel during a sustained dark window is never surfaced** — `hardware/dma_led.h:142-164` / `hardware/pov_segmented.h`. Once `dark_latched_` is set, `submitFrame()` (the only watchdog poll site) stops being called, so a wedge in a long dark window is undetected. Symptom (dark strip) equals intended output, so it is a diagnostics gap. Poll `checkStaleTransfer()` from the dark-latched idle path and/or document the dependency.

16. ✅ **`World::Trails::set_lifetime` can shrink lifetime below buffered TTLs, relying on scattered downstream clamps** — `core/filter.h:785-788,871-882` (and `Screen::Trails` at `1182-1187`). Behavior is correct via two emission-site clamps, but the out-of-range invariant is enforced ad-hoc at each use site rather than at one chokepoint. Document the clamps as the contract (a single `set_lifetime` reconcile would not remove them — the seed path produces the same condition).

17. ✅ **In-place `ArenaVector` rebind reuse does not bump `rebind_generation_`, leaving a debug use-after-free gap** — `core/memory.h:466-469`. A span snapshotted before a reuse-bind keeps stale `data_`/`size_` and `check_alive()` will not fault, asymmetric with the grow path (which does bump). Debug-only instrumentation gap, not a release hazard. Bump on the reuse path too, or document the asymmetry in the `ArenaSpan` lifetime contract.

18. ✅ **`MeshOps::compile` retains orphan vertices when degenerate faces are stripped** — `core/mesh.h:400`. Faces with <3 sides are dropped but `src.vertices` is copied wholesale, so vertices referenced only by stripped faces remain. Harmless to face-iterating renderers but inflates vertex-count-driven consumers (e.g. `MeshMorph`'s brute-force match could animate to an unrendered point). Document the wholesale copy, or compact with an index remap.

19. ✅ **`relax()` clamp comment contradicts the actual cap constant (500 vs 1000)** — `targets/wasm/wasm.cpp:1080-1081,1096`. "the editor caps the slider at 500" sits directly above `kMaxRelaxIterations = 1000`. The 1000 is intentional headroom above the editor's real 500 cap (verified in `solids.html`), but the juxtaposition reads as a contradiction. Reword to state the engine enforces 1000 as deliberate headroom.

20. ✅ **Beacon span-budget comment arithmetic is internally inconsistent with its own worst case** — `hardware/pov_sync.h:1611-1621`. The note's `5*(7+5) = 60` cols double-counts a trailing inter-burst gap after the last digit; the true worst case is ~55 cols and the conclusion (ends before HALF) holds with *more* margin. Tighten the arithmetic to match `schedule_beacon()` (4 inter-burst gaps, no terminating gap).

21. ✅ **Local death-suite skip records zero pass/fail, so absent fail-fast coverage reads green** — `tests/test_death.h:912-919,930-961`. When self-exec is blocked locally (non-CI), the entire death module — the only coverage of the marketed `HS_CHECK` layer — contributes zero assertions and the run still reports success. CI hard-fails (the real backstop). Emit a visible counted info line or a sentinel `HS_EXPECT` on the local skip path.

22. **No guard that `kModules` covers every included test header** — `tests/run_tests.cpp:12-44,60-94`. The 33 `#include`s and the `kModules[]` roster are two hand-maintained parallel lists; a header included but never registered compiles clean and silently never runs. Drive the roster from one X-macro list (as `HS_EFFECT_LIST` already does), or add a count `static_assert`.

23. **`AppState.update()` applies all mutations before notifying, so a subscriber sees a fully-advanced sibling value** — `c:/work/daydream/state.js:66`. Diverges from `set()`'s per-mutation notify; a future multi-key reactive handler reading a sibling key would observe the post-batch snapshot. Latent (no current consumer uses `update()` with sibling reads). Document the batch-then-notify ordering or notify per key.

24. **Init-time off-list effect clears per-effect param URL entries despite `preserveParams=true`** — `c:/work/daydream/daydream.js:499`. When a hydrated `?effect=` is invalid for its resolution, the correction fires `applyEffect()` with the default `preserveParams=false`, dropping accompanying `?param=` values the init path asked to preserve. Practical harm is minimal (the params target an unavailable effect). Thread the preserve intent through the hydration `set('effect')` path.

25. **`node()` divides by `(RD_N-1)` and `D_AVG` is a hand-synced literal, neither guarded against an `RD_N` change** — `core/reaction_graph.h:43,22`. A future `RD_N` edit silently leaves `D_AVG = sqrt(4π/RD_N)` wrong (it sets the RD kernel radius). The file already traps one other `RD_N`-derived value, so the silence is inconsistent. Add `static_assert(RD_N >= 2, …)` and a CI/comment pin tying `D_AVG` to `RD_N`.

### Low Priority — naming, micro-performance, comment precision, DRY

26. ✅ **Tolerance-based `Vector`/`Quaternion operator==` is scale-dependent and non-transitive, undocumented** — `core/3dmath.h:181-184,445-447`. Absolute `1e-4` per-component compare; fine for current uses (no container keys), but the doc only says "within TOLERANCE." Add a caveat or expose a named `approx_equal()`.

27. ✅ **`ease_out_cubic` uses `powf` where sibling cubics use direct multiply** — `core/easing.h:124`. Immaterial (per-frame path), idiom-inconsistent. Use `float u = 1 - t; return 1 - u*u*u;`.

28. ✅ **`random_vector` rejection loop has no documented bound or progress guard** — `core/geometry.h:788-798`. Marsaglia rejection (~1.27 expected iterations); termination relies on `rand_f()` being well-distributed. Add a one-line note on the expected-iteration bound.

29. ✅ **`ProceduralPalette::get` evaluates exact `cosf` ×3 per sample on the non-baked path** — `core/color.h:1578-1583`. Every other per-sample trig path uses `fast_cosf`; `ChaoticStrings` samples this unbaked per fragment. Switch to `fast_cosf` (matching the file's ~0.17% budget) or document a bake-before-per-pixel contract.

30. ✅ **Modulo-based `rand_int`/`random` introduce non-uniform bias on the active path** — `core/platform.h:1151-1156,567-582`. Classic modulo bias for non-power-of-two ranges; negligible for the small spawn ranges but untreated while `rand_f` is carefully unbiased. Use rejection reduction if uniformity matters, or note the cheap-but-biased intent.

31. ✅ **`ParamList::find()` and the LED correction-guard destructors carry byte-identical duplicated bodies** — `core/canvas.h:384-413`, `core/led.h:97-129`. Implement non-const `find()` via `const_cast<…>(std::as_const(*this).find(name))`; factor a shared `restore_baseline()` for the two guard dtors.

32. ✅ **`Hole::requires_unit_world_input` rationale overstates `angle_between`'s precondition** — `core/filter.h:588-590` (and `42-44,421-422`). `angle_between` normalizes internally; Hole's real sensitivity is the zero-length degenerate, not off-unit length (unlike Mobius, which genuinely needs unit input). Reword, or drop the trait from Hole.

33. ✅ **`OrientSlice` exposes mutable `enabled`/`axis` as public members while peer filters use accessors** — `core/filter.h:573-574`. Inconsistent encapsulation; a non-unit `axis` skews bucket selection (clamped, so no UB). Add `set_axis()` that renormalizes, or document the unit-length contract.

34. ✅ **`clamp_phi_band` recomputes `clamp_phi(a1)`/`clamp_phi(a2)` redundantly** — `core/sdf.h:191-207`. Both branches recompute the same two values (4 trig calls vs 2) on a cold path. Hoist `p1`/`p2` above the conditionals.

35. ✅ **`compute_inradius` documents a "dead" branch the ctor's `HS_CHECK` already makes unreachable** — `core/sdf.h:1729-1734`. The `count==0` fallback is unreachable (the ctor traps `count>0`). Tighten the comment to note the floor only guards pathological all-coincident vertices.

36. ✅ **`Star` full-width scan fallback is under-documented relative to `AngularRepeat`'s PERF-CLIFF note** — `core/sdf.h:297,2564-2573`. A pole-spanning star falls back to a full-row scan (correct geometry); the cost-model framing is less explicit than its analog. Add a one-line note. (Borderline trivial — the docstring already mentions the fallback.)

37. ✅ **Closed rings pass `close_loop=true` *and* append an overlap-close vertex, drawing a redundant degenerate segment** — `core/plot.h:1156-1164,1254`. The final segment is short-circuited by the degenerate guard (zero measurable cost) but is conceptually a double-close; `Multiline` already does it the clean way. Pick one closure mechanism per primitive.

38. ✅ **`Plot::Mesh::draw` face-walk and edge-list overloads duplicate the per-edge sample/shade/rasterize body** — `core/plot.h:1944-1972,2080-2102`. Near-identical blocks the comments explicitly say must stay in sync. Extract a private `draw_edge(...)` helper.

39. ✅ **`edge_row_span` planar comment overstates cull/renderer sample agreement** — `core/plot.h:392-399`. The unprojection *map* is shared but the cull samples projection-uniform while the renderer samples arc-uniform; correctness comes from the Lipschitz + one-row margin, not bit-identical samples. Soften the wording.

40. ✅ **`MobiusWarpEvolving::phase` can return exactly 100.0, contradicting its documented `[0,100)` range** — `core/animation.h:2177-2180`. `(h & 0xFFFF) * (100/65535)` hits 100.0 at `0xFFFF`. No functional consequence (feeds `sinf`/`cosf`). Divide by `65536.0f` for a true half-open range and uniform bucketing.

41. ✅ **`Arena::allocate` uses `reinterpret_cast<size_t>` for pointer address math instead of `uintptr_t`** — `core/memory.h:110`. Works on every target; `uintptr_t` is the strictly-correct round-trippable type and `<cstdint>` is already included. Switch for standard-conformance/self-documentation.

42. ✅ **`pair_half_edges` and the classify `HashNode` sort use `make_heap`+`sort_heap` with no rationale** — `core/mesh.h:151,667`. Heapsort vs `std::sort` (used elsewhere in the same subsystem) on a cold path; immaterial perf, conspicuous silence. Switch to `std::sort` or add a one-line justification.

43. ✅ **`narrow_index`'s `int16_t`-scratch rationale cites scratch that does not exist in `mesh.h`/`conway.h`** — `core/mesh.h:341`. The only `int16_t` `-1`-sentinel store is `hankin.h:239`; the bound is correct but the comment is not traceable from where it is enforced. Add "see `hankin.h` `face_indices`".

44. ✅ **`hankin()` one-shot's reversed arena polarity is correct but the `@param target` understates its working-set requirement** — `core/hankin.h:410-429`. `target` must hold `max(compile scratch, output)`, not just the output; a caller sizing for output alone under-provisions. Document the working-set requirement (and optionally add a tightly-sized test).

45. ✅ **`FlowField` noise-time comment understates the z-input scaling** — `effects/FlowField.h:88-90`. The parenthetical `(p.z + t)` drops the `* noise_scale` factor present in the code. Reword to `p.z*scale + t` or drop the parenthetical.

46. ✅ **`DreamBalls::update_displaced_mesh` indexed write relies on an unstated target-sizing precondition** — `effects/DreamBalls.h:334-348`. Safe today because `MeshOps::transform` pre-sizes `target.vertices`; a future reorder could silently go OOB on device (debug `assert` catches it natively). Add an always-on `HS_CHECK(target.vertices.size() == base.vertices.size())` at the cold seam.

47. ✅ **BZ `perturb_state` couples global RNG stream position to substep count** — `effects/BZReactionDiffusion.h:246-254`. Documented and intended (single-stream determinism is a deliberate project stance); informational only. Optionally pin the per-frame draw count in a determinism-test comment so a future retune is recognized as a global-stream change.

48. ✅ **Magic stagger/recurrence/pool constants in the IslamicStars ripple burst are unnamed** — `effects/IslamicStars.h:110-112`. The `16`/`96`/`8`/`144` values are coupled by a capacity invariant (overflow safely drops ripples) but live as bare literals in different methods. Hoist to named `constexpr` members with a note, mirroring `MeshFeedback`'s `static_assert`-locked period constants.

49. ✅ **`clearToolingMemory()` bumps the generation but does not free the 16 MB tooling block** — `targets/wasm/wasm.cpp:893-902`. Correct by design (block retained for module lifetime) but the name reads as if it releases memory; a JS caller will not see linear memory shrink. Rename to `resetToolingArenas()` or add a JS-facing note.

50. ✅ **Trap-shape probe re-runs the first death case twice** — `tests/test_death.h:954,963-969`. The probe spawns `cs[0]` to detect the shell's trap-relay shape, then the loop spawns it again. Subprocess spawns dominate suite runtime. Start the loop at `i=1` after asserting the probe trapped, or reuse the probe's rc for `i==0`.

51. ✅ **Windows death-child command string is fragile to special characters in the exe path** — `tests/test_death.h:805-819`. Relies on `cmd.exe` outer-quote stripping; robust for benign paths but `CreateProcess`/`_spawnv` with an argv array avoids shell parsing entirely. Document the constraint or switch to a shell-free spawn.

52. ✅ **Dead-slider lint epsilon can mask a partial per-frame revert** — `tests/test_effects.h:234-241`. A slider slowly pulled back toward its driven value (≤0.13%/frame) could still read as "persisted" within `eps` after 3 frames. No current effect exploits it (all slider drivers are absolute writes). Render more frames or assert the value is closer to the written target than to the pre-write value.

53. ✅ **`captureFrame` stretches the source into a fixed-aspect offscreen after a mid-recording resize** — `c:/work/daydream/recorder.js:233`. The 4-arg `drawImage` fills the pinned offscreen, so a window resize mid-recording distorts geometry (track size is correctly pinned; aspect is not). Compute a letterbox/pillarbox dest rect (and clear before blit), or document the limitation.

54. ✅ **`DeepLinkGUI` numeric URL value is range-clamped but not step-snapped** — `c:/work/daydream/gui.js:230`. A deep-linked off-step value (e.g. `?Segments=3` on a `step=2` control) loads a value the slider could never produce; lil-gui snaps dragged values but the URL path does not. Effect params (no `step`) are unaffected today; some global/segment/recording controls do pass steps. Snap the clamped value to the nearest step multiple when `args[2]` is finite.

55. ✅ **`findBestRationalRatio`'s `maxDenominator` parameter also bounds the numerator** — `c:/work/daydream/tools/lissajous_math.js:63-92`. Both loop terms are bounded by it (a square grid), so the name is mildly misleading (a ratio of 9 with `maxDenominator=8` snaps to 8/1). The JSDoc and tests already note "both M and N." Rename to `maxTerm`, or document the clamp-to-grid behavior.

---

## What the validators rejected (7)

Independent re-reading rejected 7 candidate findings as non-issues — chiefly cases where
a project convention (an intentional always-on guard, single-stream determinism, snapped
vs interpolated style scalars, fail-fast-by-design deferral) made the reported concern a
non-defect. This rejection rate (≈11%) is itself a signal: most of what looks suspicious
in this codebase is deliberate and documented.

## Closing assessment

The defining quality of this codebase is **discipline under constraint**. It targets a
debugger-less, console-less spinning microcontroller and a browser simultaneously from one
template-specialized engine, and it holds both targets bit-identical while never once
reaching for an unbounded allocation, a hidden global, or a silent fallback where a trap
belongs. The review surface that remains after deep scrutiny is almost entirely
*documentation* (keep the prose in step with code that has evolved past it) and
*defense-in-depth* (promote a few latent, unreachable-today guards to explicit traps,
consistent with the project's own fail-fast doctrine). There is no firefighting to do
here — only polishing.
