# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ rendering engine and its host tooling — `core/`, `effects/`,
`hardware/`, `targets/` (wasm + Phantasm), `tests/`, `scripts/`, `tools/`, and the
CMake/PlatformIO/CI build system. Out of scope by request: `core/effects_legacy.h`,
`core/rotate.h`, and `targets/Holosphere/Holosphere.ino`. Third-party vendored headers
(`FastNoiseLite.h`, `inplace_function.h`) were assessed for integration only, not audited
line-by-line.

**Method:** thirteen component-scoped reviewer passes examined every in-scope file in full,
each consulting the README to place its scope in the wider architecture. Every candidate
finding was then handed to an independent skeptic that re-read the cited source before the
finding was admitted. Of 63 candidate findings, 14 were rejected at verification (intentional
designs, miscites, or misread invariants) and are not listed here. The 49 confirmed findings
below skew heavily toward *latent* and *documentation* issues — there is no crash, leak, or
race reachable on the shipping configuration.

---

## Overall grade: **A−**

Holosphere is a mature, unusually disciplined embedded-graphics codebase. It pairs a genuinely
elegant compile-time-resolution architecture with an exacting safety culture (fail-fast
`HS_CHECK` on cold seams, explicit `Arena&` threading, host/device bit-parity treated as an
invariant) and documentation that is among the best seen in a hobby-or-professional engine of
this size. The deductions are narrow and concentrated: CI *wiring* lags the (excellent) test
**infrastructure** that already exists, a handful of README/comment figures have drifted from
the code, and a few effects carry acknowledged copy-paste debt. None of these touch correctness
on the hardware that ships.

### Quality-dimension grades

| Dimension | Grade | One-line rationale |
|---|---|---|
| Correctness & numerical robustness | **A−** | Pole/antipodal/singularity/fixed-point edges handled with rigor; remaining issues are latent edge cases (CSG seam set-difference, one-sided `acosf` clamp). |
| Memory & resource safety | **A** | Arena model is rigorously guarded — generation stamps, rebind counters, LIFO `ScratchScope` traps, no heap in steady state. |
| Concurrency / real-time / ISR | **A** | Single-core relaxed-atomic and release/acquire reasoning is correct and documented; flywheel sync is a model design. |
| Performance & efficiency | **A−** | Tight LUT-driven hot loops and baked palettes; minor unconditional per-frame rebakes and a type-erased world-filter forwarding cost. |
| Architecture & modularity | **A** | Clean layering, single-source rosters, self-registering effects, CRTP factoring of shared physics. |
| Architectural elegance (subjective) | **A** | Compile-time `<W,H>` specialization, split-trig LUT, X-macro single sources of truth, filter-ordering `static_assert`s. |
| API / interface expressiveness (subjective) | **A−** | Expressive declarative pipelines and `.then()` chaining; a few leaky abstractions effects must hand-manage. |
| Readability & code style | **A** | Consistent naming; fact-focused comments that explain *why* — occasionally verbose to the point of essay. |
| Documentation accuracy (code vs README) | **B+** | The README is exceptional but several figures/claims have drifted (arena size, azimuth range, method table, "zero-copy" example). |
| Error handling & robustness | **A−** | Fail-fast vs soft-degrade boundary applied consistently and deliberately; WASM boundary validates-and-logs untrusted input. |
| Testing & verification | **B+** | Broad suite + portable death/trap harness + provenance diffs, undercut by three CTests in zero CI shards and an unwired warning ratchet. |
| Build system & tooling / CI | **A−** | Strong reproducibility pinning and a self-testing size/layout gate; CI gate *completeness* is the weak spot. |
| Portability / platform abstraction | **A−** | First-class host/device parity; a couple of decorative `PROGMEM` qualifiers that imply a discipline the code doesn't follow. |
| Maintainability | **A−** | Highly navigable, but `animation.h` is a 2.7k-line monolith and several effect families carry hand-propagated duplication. |

### Per-component grades

| Component | Grade | Note |
|---|---|---|
| `core/` math, geometry, color | **A** | Numerically careful, well-guarded LUTs, elegant projection/sentinel design. |
| `core/` render pipeline & filters | **A−** | Correct double-buffer/Canvas RAII; one undeclared unit-input requirement, one forwarding cost. |
| `core/` rasterization (sdf/scan/plot) | **A−** | Excellent seam/pole handling; latent CSG-subtract seam edge case. |
| `core/` animation | **A−** | Rigorous lifetimes & type-erased event storage; monolithic file, perpetual-timer overflow. |
| `core/` mesh / Conway / solids / spatial | **A** | Half-edge invariants systematically enforced; minor doc drift & hash-as-id fragility. |
| `core/` platform & memory | **A** | Foundation layer; overflow-safe arena math, exact FastLED parity. **Memory safety A+.** |
| `effects/` (28 effects) | **A−** | Disciplined arena use & numerical stability; acknowledged trail/carousel duplication. |
| `hardware/` drivers & sync | **A** | Exemplary ISR/atomic discipline and host-testable protocol split. |
| `targets/wasm` bindings | **A** | Exemplary single-source param marshaling (**A+**) and memory-view detachment handling. |
| Build / test / CI / tooling | **A−** | Reproducibility & tooling code are strong; CI gate completeness is **C+**. |

---

## Detailed rationale

**Correctness & numerical robustness (A−).** The math layer clamps every `acos`/`sqrt`
argument before the call, falls back at `slerp` poles, and shares a single `STEREO_INF`
sentinel across stereographic/gnomonic projection with documented recognition thresholds.
Fixed-point color blending is provably exact at endpoints and within 1 LSB across the range,
with round-to-nearest applied uniformly so per-frame fades don't bias channels. Reaction-
diffusion stability bounds are *derived* in the parameter-registration comments. The residual
correctness items are latent: CSG `Subtract` computes its set-difference in pre-wrap angle
space (can under-carve at the θ=0 seam), and Voronoi's border `acosf` uses a one-sided clamp
that breaks the codebase-wide two-sided convention.

**Memory & resource safety (A).** This is the strongest dimension. The arena allocator uses
subtractive, overflow-safe bounds math; `ArenaVector`/`ArenaSpan` carry debug generation
stamps that catch use-after-reset *and* a separate rebind counter that catches re-grow dangles;
`ScratchScope` traps non-LIFO teardown; `Persist<T>` traps a forgotten reset via a construction
watermark. Borrow-contract `= delete` overloads on `Lerp`/`Motion`/`MobiusFlow` reject dangling
temporaries at compile time. There is no heap in steady state and every alignment/padding path
is bounds-checked before the offset advances.

**Concurrency / real-time / ISR (A).** The single-core single-observer reasoning for relaxed
atomics in `dma_led.h`, the release/acquire publish of the freshly-constructed effect in
`pov_segmented.h`, and the single-writer flywheel with a pure-publisher edge mailbox in
`pov_sync.h` are each correct and carefully justified. The double-buffer index protocol
(`prev_`/`cur_`/`next_`) documents its single-core dependency and routes a stalled-ISR through
`HS_CHECK` rather than hanging headless. The one asymmetry: the single-board driver's cross-
context `effect_` pointer is non-atomic (safe by attach/detach ordering) where the segmented
driver uses atomics for the analogous handoff.

**Performance & efficiency (A−).** Hot paths are `FASTRUN`, branchless per-pixel, LUT-driven,
and bypass virtual `get_pixel` dispatch via a display-buffer fast path. Baked-palette LUTs
replace per-pixel OKLCH lerps. The deductions are modest and mostly per-frame waste: `Moire`
rebakes two 256-entry LUTs every frame with no static-palette skip, `HopfFibration` round-trips
`Vector`→`Spherical` per fiber per frame, and World-filter stages forward through a type-erased
`PassFn3D` (an indirect call per stage per point) which sits in mild tension with the README's
"zero-overhead" framing.

**Architecture & elegance (A).** Three build targets share one engine through compile-time
`<W,H>` specialization. Rosters are single-sourced (`HS_EFFECT_LIST`, `HS_*_RESOLUTIONS`,
`MESHOP_LIST`) and cross-checked at startup/test-gen so a target can't drift. Effects self-
register via a Meyers singleton with no static-init-order hazard. Filters live strictly in
`Filter::{World,Screen,Pixel}` with the `Feedback::Style` *data* namespace kept separate, and
compile-time ordering `static_assert`s (terminal-last, 2D-not-before-3D, non-unit-emitter
ordering) push pipeline invariants into the type system. `ReactionDiffusionBase` is a model
CRTP factoring that shares lattice/kernel scaffolding while leaving divergent physics in leaves.

**API & interface expressiveness (A−).** The variadic `Pipeline` with compile-time `get<T>()`,
the `.then()` fluent chaining, the friend-gated parameter write path, and the WASM parameter-
marshaling SSOT are all expressive and hard to misuse. The friction points are leaky
abstractions a few effects must hand-manage (e.g. `ChaoticStrings` mirroring live slider values
into a transformer's internal `entities`, `Comets` hand-snapping params into a 16-byte inline
capacity) and a couple of unconstrained generic boundaries (`Tweenable::get()` return type).

**Documentation (A readability / B+ accuracy).** Comments explain *why* at exactly the points
of subtlety — relaxed-atomic rationale, clamp-before-scale, fail-fast placement, approximation
error budgets. The accuracy deduction is for real drift the reviewers caught against the
otherwise-superb README: the `≈335 KB` arena figure (actual 330 KiB), the `Spherical`-from-
`Vector` azimuth range vs the documented `[0,2π)`, two undocumented exported WASM methods, a
"zero-copy bound once" example that contradicts the documented re-fetch contract, and a stale
Pages URL in a docs workflow.

**Testing, build & CI (B+ / A−).** The test infrastructure is excellent: 33 module suites
covering essentially every subsystem, stack/arena high-water probes, a WASM runtime smoke that
drives every effect at every resolution, generator provenance diffs, and a genuinely portable
death/trap harness that probes the shell's signal-relay shape and fails loud under CI. What
drags the grade is *wiring*, not capability: three registered CTests (`concepts`,
`h_offset_renorm`, `stack_budget`) match zero CI shard regexes — so they get no Linux/sanitizer
coverage *and* the shard-coverage self-gate that should flag them appears to be bypassed — and
the fully-built, baselined, unit-tested Teensy warning ratchet is invoked by no workflow.

---

## Prioritized fix list

Items are numbered sequentially across all priority sections. Each cites `file:line` and the
verified impact. No Critical-severity defects were found.

### Critical

*(none)*

### High

1. ✅ **Three registered CTests run in zero CI shards** — `.github/workflows/ci.yml:37-45` vs `tests/CMakeLists.txt:96-100,148,163`. `unit_concepts`, `unit_h_offset_renorm`, and `unit_stack_budget` match no shard regex, so they get no Linux/sanitizer coverage, and the `shard-coverage` self-gate (which `exit 1`s on any test matching ≠1 shard) should be failing CI — meaning either master is red or the gate is bypassed. The "partition all 32 modules" comment also undercounts the 33 registered. *Fix:* add the three to shard regexes, correct the module count, and confirm `shard-coverage` passes.

### Medium

2. ✅ **Teensy warning-hygiene ratchet is implemented, baselined, and unit-tested but wired into no workflow** — `tools/teensy_warnings.py` (whole file); cf. `.github/workflows/ci.yml:467-510`, `platformio.ini:30-33`. The spec treats the ratchet as a load-bearing gate, but the `teensy-size` job runs only a *cached* `pio run`, and the ratchet's own docstring requires a cache-disabled build (a cached TU emits no warnings). Net: new first-party Teensy warnings ride green CI. *Fix:* add a CI step that does a cache-disabled Teensy build, captures the compiler log, and runs `teensy_warnings.py --build-log <log> --github`.

### Low

3. ✅ **`Spherical`-from-`Vector` azimuth contradicts the `[0,2π)` convention** — `core/3dmath.h:374-381`. `theta = fast_atan2(n.z, n.x)` yields `(-π, π]`; consumers reading `s.theta` raw (e.g. `filter.h:1165`) get negatives and must hand-wrap. *Fix:* wrap into `[0,2π)` in the constructor, or document the `(-π,π]` range on `Spherical`.
4. ✅ **`geometry.h` uses `std::pair` without `#include <utility>`** — `core/geometry.h:818`. Compiles only via a transitive include. *Fix:* add `#include <utility>`.
5. ✅ **`Tweenable` concept doesn't constrain `get()`'s element type** — `core/concepts.h:317-321`. Non-conforming containers fail deep inside `slerp`/`lerp` with a worse diagnostic; `length()` admits signed `int` with no non-negativity contract. *Fix:* constrain the deduced `get()` result (or at least document the element contract).
6. ✅ **`concepts.h` depends on `Pixel` and the global `Fn` alias only transitively** — `core/concepts.h:212-307`. Resolves only because `canvas.h` pulls in `color.h`/`platform.h`. *Fix:* include `color.h` and `platform.h` directly.
7. ✅ **LUT generator has no round-trip / monotonicity self-check** — `scripts/generate_luts.py:64-73,124-142`. Regeneration correctness rests solely on a CI token diff; a future libm change could land a one-entry shift silently. *Fix:* add a `--check` pass asserting each table is non-decreasing and the sRGB→linear→sRGB round trip stays within ±1 code.
8. ✅ **NaN `scale` silently maps a pixel to full white** — `core/color.h:130-134`. `hs::clamp` maps NaN to `hi`, so a non-finite per-frame feedback fade would flash the buffer to white and persist it, masking the upstream NaN. Currently unreachable. *Fix:* if tightened, trap a non-finite fade at the cold sync site where it is computed.
9. ✅ **`World::OrientSlice` assumes unit input but omits `requires_unit_world_input`** — `core/filter.h:435-499`. Unlike `Hole`/`Mobius`, it doesn't declare the flag, so a `World::Trails` (which re-emits non-renormalized vectors) placed upstream compiles silently and biases bucket selection. Latent — no shipping pipeline pairs them. *Fix:* add `static constexpr bool requires_unit_world_input = true;` (or normalize `v` inside `plot()`).
10. ✅ **World-filter stages forward through a type-erased `FunctionRef` (indirect call per stage per point)** — `core/filter.h:27-29,313-324,420-425,571-579,612-617`. High-fanout filters (`Replicate`, `VertexReplicate`, `Trails::flush`) pay a per-point indirect call; only `Screen::AntiAlias`/`Blur` use an inlining forwarding-reference. *Fix:* template the World-filter pass parameter as `Screen::AntiAlias` does, or soften the README's "zero-overhead" wording.
11. ✅ **`Pipeline` int-coordinate terminal vs recursive overloads wrap via different paths** — `core/filter.h:138-147` vs `299-301`. A filter-less pipeline wraps the int directly; a filtered one promotes to float and `round`s. Both correct for in-range ints, but the divergence is a maintenance trap. *Fix:* unify the paths or add a one-line note on the deliberate divergence.
12. ✅ **CSG `Subtract` set-difference runs in pre-wrap angle space, missing seam-straddling overlaps** — `core/sdf.h:1056-1092`. When A's carve region and B's removal region describe the same θ=0 band in different wrap frames, the disjointness test misses the overlap and a solid B under-carves at the seam. *Fix:* normalize both children into `[0,W)` (seam-split) before the set-difference, and add a Subtract-straddling-θ=0 test.
13. ✅ **`DistortedRing` `max_distortion` is an unguarded correctness precondition** — `core/sdf.h:624-652,740`. If a caller under-bounds the true sup of `shift_fn`, the shape silently disappears in patches; every sibling guards its analogous invariant. *Fix:* sample `shift_fn` at construction and `HS_CHECK` the supplied bound (cold path), or at least a debug-only sampled check.
14. ✅ **Ring centerline fast path depends on `scan_region` coalescing for the azimuthal seam but doesn't document it** — `core/sdf.h:518-538`. The two padded twin arcs can both straddle θ=0; correctness relies entirely on the downstream coalescer (the pole case *is* documented; this one isn't). *Fix:* add a one-line note mirroring the `pole_wrap` comment and verify a CSG parent's unwrapped-space merge can't double-count.
15. **Mesh face/edge index stored in `Fragment::v2` as `float` loses precision past 2²⁴** — `core/scan.h:800-803`, `core/plot.h:1909,1914`. An unstated capacity ceiling on an otherwise general API; no shipped mesh approaches it. *Fix:* document the 2²⁴ ceiling on the `v2` register convention (or thread large indices through an integer channel).
16. ✅ **Perpetual timers' `int` frame counter overflows with no documented bound** — `core/animation.h:336,734,743,783,803`. `RandomTimer`/`PeriodicTimer`/`Driver`/`RandomWalk` run forever on `duration=-1`; `t++`/`next=t+period` are signed-overflow UB at ~2³¹ frames, after which `t >= next` comparisons misbehave. *Fix:* compare elapsed deltas with unsigned arithmetic (defined wrap), or document the bound as `Noise` does.
17. ✅ **Negative duration silently turns a finite animation perpetual** — `core/animation.h:325-327,228-230`. Only `duration == 0` is remapped; a negative passes through and `done()` is permanently false, so a one-shot never completes or fires `.then()`. *Fix:* `HS_CHECK(duration >= 0 || duration == -1)` in the ctor (or clamp like the 0 case), reserving `-1` as the sole perpetual sentinel.
18. ✅ **`MeshMorph` doc claims transients are freed on destruction, but only `compact()` reclaims them** — `core/animation.h:1860-1862,1974-1989`. The class has no destructor; arena bytes persist until the caller compacts/resets. *Fix:* reword to state object destruction runs only member dtors, not arena reclamation.
19. ✅ **`generate()`'s `[[nodiscard]]` is ineffective for void-returning generators** — `core/generators.h:45-46`. The attribute is silently ignored for the common in-place void case and forces a `(void)` cast for intentional discards otherwise. *Fix:* drop `[[nodiscard]]` (primary effect is the in-place write) or scope it to value-returning generators.
20. ✅ **Paused-gate boilerplate duplicated across four animation types** — `core/animation.h:911-912,988-989,1094-1095,1151-1154`. `Mutation`/`Driver`/`Lerp`/`Sprite` repeat the pause gate and member; future pause-semantics changes touch four sites. *Fix:* factor into a protected `AnimationBase::is_paused()` helper.
21. ✅ **`KDTree::nearest` doc says `k` is "capped at MAX_K" but the code traps** — `core/spatial.h:98-107`. The implementation `HS_CHECK`-traps `k > MAX_K` (correct fail-fast); only the point-count limit is a soft cap, so a caller trusting the doc hard-crashes. *Fix:* reword the doc; keep the trap.
22. **`classify_faces_by_topology` uses a 32-bit hash directly as the topology class id** — `core/mesh.h:528-533,644-652`. Two distinct topologies that collide in `fmix32` merge into one palette class; whole-degree angle rounding widens the collision surface. Acceptable for the fixed roster. *Fix:* widen the angle quantization, or add a per-solid distinct-class-count test.
23. **`compile()` copies all vertices wholesale, leaving orphans from stripped degenerate faces** — `core/mesh.h:388-394`. `num_vertices()` can exceed the referenced set — a hazard for index-mapped per-vertex data (e.g. `MeshMorph`). *Fix:* document the invariant at `MeshState`/`num_vertices()`, or add an optional compaction pass that remaps and drops unreferenced vertices.
24. ✅ **`update_hankin` re-binds output vectors every call; the "no reallocation" claim is conditional** — `core/hankin.h:334-361`. Allocation is avoided only when `out_mesh` is reused against the same sufficiently-sized arena. *Fix:* tighten the README wording to the steady-state-reuse condition.
25. ✅ **`PROGMEM` on the `neighbors[]` definition but not its declaration** — `core/reaction_graph.cpp:4` vs `core/reaction_graph.h:60`. Harmless only because all live targets flatten `PROGMEM` and use a flat address space; the qualifier implies a flash-access contract the plain `[]` reads don't honor. *Fix:* drop `PROGMEM` from the generated definition, or make the generator emit a matching declaration plus a note that direct subscripting is intentional.
26. ✅ **Stale arena-size figure in `memory.cpp` comment (335 KB vs actual 330 KiB)** — `core/memory.cpp:48`. The authoritative `DEVICE_GLOBAL_ARENA_SIZE` is `330*1024`; this is the one place a reader looks for the RAM decision. *Fix:* update to 330 KiB.
27. ✅ **Header documents a negative-neighbor sentinel the generator never emits** — `core/reaction_graph.h:56-58`. For `RD_N=7680` every entry is a valid index; the doc presents the negative case as a live table state. *Fix:* reword to state the table is always fully populated and the `ni < 0` guard is defensive belt-and-suspenders.
28. ✅ **`BZReactionDiffusion` odd-substep copy-back is dead code; the documented contract is unexercised** — `effects/BZReactionDiffusion.h:151,357-370`. `STEPS_PER_FRAME=2` (even) means the parity copy-back never runs; the sibling GS path (16 substeps) likewise never exercises it. *Fix:* `static_assert` even substep count and drop the branch, or add a white-box test driving an odd count.
29. **`Comets`/`ChaoticStrings`/`RingSpin` duplicate the trail record+`deep_tween` skeleton with hand-propagated fixes** — `effects/Comets.h:26-29,131-144`. Acknowledged in-header; any trail-fix must be edited in N places. *Fix:* extract a small CRTP/helper (as done for `ReactionDiffusionBase`) so the shared mechanics live once.
30. **`ChaoticStrings` manually mirrors live slider values into the noise transformer's `entities` each frame** — `effects/ChaoticStrings.h:137-146`. Couples the effect to the transformer's internal entity layout. *Fix:* provide a `NoiseTransformer::set_live_params(...)` (or bind to a `Params&`) so effects don't iterate the internal collection.
31. ✅ **`Comets` closing-domain floor silently rewrites any future sub-π authored entry** — `effects/Comets.h:167-172,233-245`. All 12 current entries clear the threshold; the `CometsWhiteBox` test checks closure, not floor activation. *Fix:* have the test assert every entry's `m2*domain >= PI` so the floor can never silently engage.
32. ✅ **`Flyby`/`FlowField` handle the same unbounded-float-time risk inconsistently** — `effects/Flyby.h:97-99,176`. `FlowField` wraps its accumulator to a period; `Flyby` leaves `noise_time` unbounded (ULP swallows the increment after ~2²⁴, freezing the field after days of uptime). *Fix:* apply the `FlowField`-style `fmod` wrap in `Flyby`, or document why the install runtime makes it acceptable here.
33. ✅ **`Moire` rebakes two 256-entry OKLCH LUTs every frame with no static-palette skip** — `effects/Moire.h:91-92`. 512 OKLCH evaluations/frame even during quiescent gaps; inline doc calls it "cheap" but offers no guard. *Fix:* track a dirty flag / palette generation and skip rebake when unchanged (mirror `apply_if_changed`).
34. ✅ **`FLASHMEM` applied to `HankinSolids::start_morph_cycle`, which runs on the per-cycle runtime path** — `effects/HankinSolids.h:219`. `FLASHMEM` is for cold one-time setup; this recurring path pays slower flash instruction fetch. *Fix:* drop `FLASHMEM` here (reserve it for ctors/`init()`).
35. **`MindSplatter` particle-index bound is a bare `assert` with no cold trap at the bind site** — `effects/MindSplatter.h:270-273`. The README's hot-path discipline is a stripped `assert` *backed by* a cold `HS_CHECK` at setup; here there's no backstop establishing `f.v2 ∈ [0, active_count)`. *Fix:* add a once-per-frame cold `HS_CHECK(active_count <= pool.capacity())` (or document the `v2`-is-pool-index contract).
36. ✅ **`HopfFibration` re-derives `Spherical` from each fiber `Vector` every projection** — `effects/HopfFibration.h:180,212-216`. ~210 trig-heavy `Vector`→`Spherical` round-trips/frame to recover angles it already had at init. *Fix:* store base fibers as `Spherical` (or a precomputed `{azimuth, polar}` pair).
37. **Duplicated load/classify/draw machinery between `HankinSolids` and `IslamicStars`** — `effects/HankinSolids.h:115-162`. Same carousel + `palettes_slots[2]` + near-identical `draw_shape`/classify idiom, copy-paste-diverged. *Fix:* extract a shared `MeshShapeCarousel<W,H>` mixin; effects differ only in their generate step.
38. ✅ **`IslamicStars` cross-fade window `HS_CHECK` is tautological after the preceding clamp** — `effects/IslamicStars.h:201-202`. `fade` is already `std::min`-clamped to `dur/2`, so the guard can never fire and misleads about what's possible. *Fix:* drop the `HS_CHECK` (or convert to a why-comment).
39. ✅ **`Voronoi` border `acosf` lacks a lower-bound clamp, breaking the codebase convention** — `effects/Voronoi.h:137-138`. A dot can land marginally below `-1.0f` under accumulated rotation drift; `acosf(<-1)` returns NaN and silently drops the border pixel. This is the lone one-sided `acosf` in the tree. *Fix:* `acosf(hs::clamp(maxDot, -1.0f, 1.0f))` on both.
40. ✅ **`RingSpin` retains `source_palettes` in persistent RAM after baking** — `effects/RingSpin.h:170-173`. Consumed only in `init()`; rings color from the baked LUTs. Dead footprint, not a leak. *Fix:* bake from the global `Palettes::` constants and drop the member (or document intentional retention).
41. ✅ **Cross-effect duplication of the age-driven recyclable-slot pattern is hand-propagated** — `effects/RingShower.h:68-79,120-134`. `RingShower` and `Thrusters` independently implement the same `age→growth+fade, expire when age>=life` lifecycle (RingShower's comment says it "mirrors Thrusters"). *Fix:* a tiny shared `AgedSlot` helper (`age`, `life`, `radius_at`, `expired()`), or cross-reference both directions.
42. ✅ **`PROGMEM` LUTs indexed via plain `operator[]` in the per-pixel hot path** — `hardware/hd107s_frame.h:155-157,161-163,205-207`. Correct on Teensy 4.x (unified address space) but the qualifier implies an AVR-style separate space where direct `[]` would read garbage — a portability/readability landmine. *Fix:* add a one-line comment that `PROGMEM` is a placement no-op on the only supported target and direct subscripting is intentional.
43. ✅ **README §7.10 "FastLED fallback remains" doesn't reflect that the segmented driver hard-errors without `USE_DMA_LEDS`** — `hardware/pov_segmented.h:52-56`. `POVSegmented` is a compile error without DMA LEDs; the general README statement implies a fallback exists. *Fix:* note in §7.10 that the FastLED fallback exists only for the single-board `POVDisplay`.
44. ✅ **README API table omits two exported `HolosphereEngine` methods (`getParamGeneration`, `strobeColumns`)** — `README.md:1832-1847` vs `targets/wasm/wasm.cpp:1235,1242`. `getParamGeneration()` is the documented cache-invalidation handle a JS consumer must poll. *Fix:* add both rows to the §10.2 table with a note on the param re-fetch protocol.
45. ✅ **README zero-copy example binds the pixel view "once", contradicting the documented detach/re-fetch contract** — `README.md:1857-1864` vs `targets/wasm/wasm.cpp:484-500`. `ALLOW_MEMORY_GROWTH` (e.g. the lazy 16 MB MeshOps malloc) detaches the `ArrayBuffer`; copying the snippet verbatim ships a latent black-frame-after-growth bug. *Fix:* annotate the example with the per-frame re-fetch guard (mirroring `daydream.js refreshPixelView`).
46. ✅ **Shared module-global WASM tooling scratch arenas rely on an unenforced single-call reentrancy assumption** — `targets/wasm/wasm.cpp:76-83,798-799,860-861,887-888`. Safe today (per-worker module, synchronous calls) and acknowledged, but an async refactor would silently alias scratch with no trap. *Fix:* add a cheap module-global in-op flag with an `HS_CHECK` on re-entry at each MeshOps entry point (cold tooling, so free).
47. ✅ **Windows death-harness child spawn relies on cmd.exe-safe self-exe paths** — `tests/test_death.h:797-816`. `std::system` runs `cmd /c`; a checkout path containing `&`/`%`/quotes could mis-parse and spuriously skip or misclassify the death module. Self-documented. *Fix:* switch the Windows child spawn to a shell-free API (`_spawnv`/`CreateProcess`).
48. **`teensy_budgets.json` Phantasm reaction_graph/DMA invariants are documented as unverified, conflicting with the "builds green" CI claim** — `tools/teensy_budgets.json:63-69`. The budget file still says "UNVERIFIED on a real ELF" while a real Phantasm ELF fixture now exists and CI claims green; the DMA TX-buffer OCRAM invariant remains a deferred TODO. *Fix:* reconcile the comment with the captured fixture and land the deferred `DMALEDController` eDMA TX-buffer OCRAM invariant.
49. **LUT/reaction-graph provenance jobs use the runner's ambient `python3` with no `setup-python` pin** — `.github/workflows/ci.yml:327-349,351-373`. The generators reproduce bit-for-bit only in IEEE double (CPython floats are doubles, so currently safe), but this breaks the otherwise-strict pinning discipline the teensy jobs follow. *Fix:* add `actions/setup-python@v5` with a pinned `python-version` to both provenance jobs.

---

## Notable strengths (for balance)

- **Fail-fast safety culture done right.** `HS_CHECK` traps on cold seams (arena OOM, capacity,
  registration), survives `NDEBUG` in the optimized device build, and is verified by a death
  harness that asserts the traps actually fire. Hot paths use stripped `assert`s backed by cold
  traps at setup — a deliberate, documented split.
- **Memory discipline.** Generation stamps, rebind counters, LIFO `ScratchScope` enforcement,
  and compile-time `= delete` borrow contracts make the no-heap-in-steady-state guarantee
  *checkable*, not aspirational.
- **Single sources of truth.** `HS_EFFECT_LIST`, `HS_*_RESOLUTIONS`, `MESHOP_LIST`, and the WASM
  parameter-marshaling pass are each derived once and cross-checked, so a target/roster cannot
  drift unnoticed.
- **Host-testable embedded math.** Protocol/timing cores (`pov_sync.h`, `pov_segment_map.h`,
  `hd107s_frame.h`, `param_marshal.h`) are kept Arduino-free and unit-tested on the host — a
  model for making interrupt/protocol code verifiable.
- **Documentation that earns its length.** Comments justify *why* at the points of subtlety and,
  where checked, match the code; the README is a genuine architectural reference.
