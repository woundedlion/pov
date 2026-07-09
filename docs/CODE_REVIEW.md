# Holosphere / daydream — Code Quality Review

**Overall grade: A− (3.63 / 4.3)**

A two-repository review of the Holosphere C++ rendering engine + firmware and the
daydream web simulator. Every rendering, math, mesh, color, animation, hardware,
build, and test file in scope was read in full; every daydream app, worker, tool,
and test module was read in full.

> **Scope.** Out of scope per instruction: `effects_legacy*`, `targets/Holosphere/Holosphere.ino`,
> `core/math/rotate.h`. `core/vendor/FastNoiseLite.h` (third-party) was skipped;
> `color_luts.h` and `reaction_graph.cpp` were treated as generated data (their
> *generators* and invariants were reviewed, not the tables). Project-stated design
> choices — fail-fast `HS_CHECK` traps, no-heap arena allocation, the 2-physical-buffer
> ISR design, always-on wrap guards, compile-time `<W,H>` specialization — were judged
> *as intended behavior*, not flagged as defects.

## Methodology

The review ran as a deterministic multi-agent fan-out: **26 file-group reviewers**
worked in parallel, and **every finding each produced was handed to a fresh,
independent verifier** that re-read the cited code and rendered a verdict before the
finding was allowed into this report. Of **86 raw findings, 28 (33%) were rejected**
by verification as non-issues, already-handled, or deliberate design — leaving **58
confirmed or plausible findings** below. 48 are CONFIRMED (verifier reproduced the
failure by reading the code); 10 are PLAUSIBLE (real but reachability depends on
unread context). 112 agents, zero errors.

**The headline result is what is *not* here:** across ~107,000 lines, the review
surfaced **zero critical and zero high-severity defects**. The distribution is 1
medium, 38 low, 19 nit — i.e. every confirmed issue is a refinement, a latent
landmine behind an input no current caller supplies, a test-coverage gap, or a doc
drift. That profile is characteristic of a mature, disciplined codebase, not a
work-in-progress.

---

## Quality dimensions

| Dimension | Grade | Rationale |
|---|---|---|
| **Documentation** | **A** | Best-in-class. Comments name the *concrete bug each design choice avoids* rather than restating code; error bounds, epsilon magnitudes, and branch thresholds are quantified and sourced. The README is a 2,168-line architecture document that matches the code. |
| **Safety / robustness** | **A−** | Fail-fast is applied with judgment: `HS_CHECK` traps on cold setup/growth/OOM paths, stripped `assert` on hot paths backed by a cold trap, float→int casts routed through a NaN→hi clamp that is itself compile-time-guarded against `-ffast-math`. Ownership encoded in the type system (move-only `ArenaVector`, borrow-only `ArenaSpan`, rvalue-rejecting `StoredFunctionRef`). |
| **Performance** | **A−** | First-rate embedded discipline: split sin/cos LUTs, `fast_*` transcendental approximations with stated error bounds, trivially-copyable types for memcpy paths, zero-indirection inlined filter pipeline, no heap on any render path. A few local inconsistencies (exact `expf` beside `fast_acos`; redundant per-frame planar projection). |
| **Elegance** | **A−** | Compile-time polymorphism used tastefully (variadic filter `Pipeline`, `requires`-detected segue hooks, CRTP reaction-diffusion base, X-macro anti-drift). A handful of god-objects / grab-bag headers keep it from an A. |
| **Architecture** | **A−** | Clean layering and one engine driving three targets (device / segmented device / WASM) with genuinely shared code. Weak spots: `Effect` base class does too much; several headers (`platform.h`, `spatial.h`, `concepts.h`, `sdf.h`) mix concerns. |
| **Correctness** | **B+** | No shipping-path crashes found. The confirmed correctness issues are latent (guards missing for inputs no roster entry currently supplies) or minor visible artifacts (seam over-plotting). Numerical singularity handling is systematic. |
| **API design** | **B+** | Lifetime intent is largely enforced by type, not convention — but several load-bearing contracts (`get_pixel`↔`overrides_get_pixel`, `ArenaVector` no-destructor, Arena non-copy) rest on documentation only, a recurring "forget-proof-ness gap." |
| **Testing** | **B+** | Wide, layered, anti-drift-obsessed (X-macro module lists, a professional death harness, determinism perturbation, parity tests driving real WASM). Gaps: per-*module* (not per-case) assertion floors, adaptive-regime coverage holes, no coverage instrumentation. |
| **Maintainability** | **B+** | High cohesion locally, strong single-source-of-truth mechanisms. Held back by cross-file duplication (three trail effects, two shader/dark-path pairs, gen-script helpers) that the code's own comments flag as "propagate fixes by hand." |

## Per-component grades

| Component | Grade | | Component | Grade |
|---|:--:|---|---|:--:|
| `core/math` | A− | | `core/render/filter.h` | A− |
| `core/engine` memory + platform | A− | | `core/animation` base/params | B+ |
| engine callables/registry | A− | | `core/animation` motion/mesh | A− |
| transformers/styles | A− | | hardware sync + segmented POV | **A** |
| mesh core (PolyMesh/KDTree) | A− | | hardware single/DMA/HD107S | A− |
| mesh ops (Conway/Hankin/solids) | A− | | WASM target | A− |
| color system + palettes | A− | | effects A / B / C | A− |
| render canvas/scan/shading | A− | | build system + tooling | B+ |
| render SDF (`sdf.h`) | B+ | | C++ test suite | A− |
| render plot (`plot.h`) | A− | | daydream core / gui / segment / tools / tests | A− |

The lowest-graded component (`sdf.h`, B+) is a 3,600-line header that mixes two
rendering models; it still had exactly one confirmed finding, and that finding only
over-scans a few rows (never drops geometry).

---

## What the codebase does exceptionally well

- **Fail-fast with a scalpel.** The hot/cold split is applied per call site with the
  rationale written down: cold paths trap via `HS_CHECK` (survives `NDEBUG`, no stdio),
  hot paths use a stripped `assert` backed by a cold trap at the bind site. The death
  harness *verifies the traps actually fire* (`SIGILL`) rather than assuming it.
- **Sim/device bit-parity as a first-class invariant.** Shared `Pcg32(1337)` seed,
  bit-exact FastLED reimplementations, `uint32`-narrowed `millis`, an exhaustively-verified
  `beat88` argument — and a cross-run byte-identical determinism test that pins it.
- **The Phantasm 1-wire sync protocol** (`pov_sync.h`) is a standout: a count-coded,
  odd-only, distance-2 symbol alphabet that provably degrades "to missed, never to wrong,"
  a single-writer concurrency model with one fused claim point, and a cycle-counter
  flywheel that makes masked-IRQ windows structurally unable to drop columns. Fully host-tested.
- **Anti-drift by construction.** X-macro rosters, bidirectional CMake↔source module-list
  pins, param-order zip checks, `install(CODE)` WASM provenance hashing — the codebase
  systematically makes "add it in one place, forget the other" a compile error or CI failure.
- **Numerical honesty.** Every approximation states its error bound; every singularity
  (stereographic/gnomonic poles, antipodal slerp, gimbal, degenerate cross products) has
  an explicit, documented fallback.

---

## Prioritized findings

Every confirmed/plausible defect is listed below, numbered sequentially, most-impactful
first. Tags: **[severity/verdict/category]**. Each item cites `file:line`, the concrete
failure, and a fix — sized for the `code-review-fix` workflow.

### Priority 1 — Fix first (reachable today: visible artifacts, data loss, or a test gap over the production regime)

1. ✅ **Closed `Multiline` double-plots the seam vertex** — [low/CONFIRMED/correctness] `core/render/plot.h:1080`. Closed mode appends a duplicate of vertex 0 *and* passes `close_loop=false`, so the seam pixel is deposited twice; with alpha<1 `over`-blend it over-saturates into a persistently brighter dot at every closed-loop join (trail/ribbon effects). Fix: keep the wrap vertex but pass `close_loop=true` so the degenerate closing segment is omitted.
2. ✅ **Closed `SplineChain` double-plots the seam control point** — [low/CONFIRMED/correctness] `core/render/plot.h:2320`. Same root cause; `SplineFlow` (closed, alpha<1 through `World::Trails`) exercises it today, over-blending the join. Fix: emit the closing knot once (stop the final segment at `j<S`, or route closure through `close_loop=true`).
3. ✅ **`Lerp` omits the `duration>=0` guard its five siblings enforce** — [low/CONFIRMED/correctness] `core/animation/params.h:271`. A `-1` duration is accepted, `step()` clamps progress to 0 forever, `done()` never fires, the slot is never freed, and any `.then()` chain (preset cycling depends on this) stalls permanently. Fix: add `HS_CHECK(duration >= 0, ...)` matching `Transition`/`Mutation`/`ColorWipe`/`MobiusWarp`.
4. ✅ **Voronoi coherence test covers only the `B=8` low-density regime** — [medium/CONFIRMED/test-coverage] `tests/test_effects.h:1400`. The test hardcodes `B=8`/6 sites, but production defaults to 200 sites (`B≈6`) up to 400 (`B=4`); the entire adaptive `B=4..7` regime — including the default and the `COHERENCE_BLOCK_MIN=4` corner-indexing worst case — is never exercised, so a regression in adaptive block sizing passes CI. Fix: add a high-site-count case that reads `B` from the effect's own `cell_px`/`COHERENCE_BLOCK_MIN` computation and asserts fast-path equivalence.
5. ✅ **`showSaveFilePicker` truncates the user's chosen file before any data is confirmed** — [nit/CONFIRMED/correctness] `daydream/recorder.js:366`. `createWritable()` runs eagerly under the start gesture and zero-truncates the target; a recording that then produces no data still calls `writable.close()` on the success branch, leaving a 0-byte file over the user's previously-existing file. Real data loss. Fix: only `createWritable()` once the first chunk arrives, or detect the no-data streaming case in `finish()` and skip/remove the empty file. (Nit by engine severity, elevated for user-data impact.)

### Priority 2 — Latent correctness landmines and API/concurrency footguns (real defects; undefended but not reached by a current caller)

6. ✅ **`long_way` Quaternion `slerp` traps on (near-)identical endpoints** — [low/CONFIRMED/correctness] `core/math/3dmath.h:1292`. `slerp(q,q,0.5f,/*long_way=*/true)`: the sign fixup drives `d≈-1`, the degenerate-lerp branch reaches the zero quaternion, and `normalized()` traps. The `Vector` overload handles this antipodal case; the Quaternion one doesn't. Fix: synthesize a half-turn about an arbitrary orthogonal axis (mirror the Vector path) or `HS_CHECK` the precondition.
7. ✅ **`ClipRegion` margin expansion over-wraps** — [low/CONFIRMED/correctness] `core/engine/constants.h:99`. When the margin-expanded band width lands strictly between `w` and `2w`, `contains_x()` folds to a thin wrap sliver and blanks most of the frame (e.g. `x_start=5,x_end=283,margin=8,w=288` → column 100 wrongly clipped). Fix: test `(x_end-x_start)+2*margin >= w` for full-coverage, not just the raw display width.
8. ✅ **SDF arc-extremum refinement reuses a stale compacted plane on a degenerate edge** — [low/CONFIRMED/correctness] `core/render/sdf.h:2474`. A zero-length edge pushes no plane, so `refine_phi_from_arc_extremum` indexes the *previous* edge's normal against this edge's endpoints. Only ever over-scans rows (never drops geometry), but a genuine index/data mismatch. Fix: guard the refine call with the same non-degeneracy test used to push the plane.
9. ✅ **`Scan::Shader::draw` populates the vertical render margin but not the horizontal one** — [low/PLAUSIBLE/correctness] `core/render/scan.h:891`. The `y` loop is margin-expanded while the `x` loop uses the raw display band (unlike `Scan::rasterize`, which clips both axes). A future segment-clipped shader effect with a horizontal neighbor-reading filter would see a left/right seam. Fix: make both axes agree (iterate full row with `x_clip()`, or drop the `y`-margin expansion).
10. ✅ **Single-board strip mapping silently misregisters on odd effect width** — [low/CONFIRMED/correctness] `hardware/pov_single.h:108`. `run()` traps `height != S/2` but nothing constrains width parity, while `strip_opposite_col`'s `(x+w/2)%w` antipode and the `x==w/2` swap cadence both truncate `w/2`. An odd-width effect misregisters the bottom hemisphere and skews frame cadence. Fix: add a sibling `HS_CHECK(width % 2 == 0, ...)`.
11. ✅ **`get_pixel`-override roster incompatibility traps mid-show in the ISR, not at boot** — [low/CONFIRMED/correctness] `hardware/pov_segmented.h:478`. The boot path validates height/width but not `overrides_get_pixel()==false`; a mis-rostered effect boots fine and `__builtin_trap()`s in ISR context when its epoch turn arrives, crashing the spinning device in the field. Fix: validate at boot alongside the height/width checks (ideally sweep all factories at `run_show()` entry).
12. ✅ **Head-first flush lets an upstream 2D-history stage's emissions be overwritten by terminal `Feedback`** — [low/CONFIRMED/architecture] `core/render/filter.h:392`. `Pipeline<W,H, Screen::Trails, Pixel::Feedback>` passes every static_assert, yet `Feedback::flush`'s opaque store clobbers the trails the preceding stage just re-emitted. Fix: static_assert no emitting 2D-history stage may precede a terminal filter, or have `Feedback` over-blend instead of opaque-store.
13. ✅ **Beacon in flight when the master crosses HALF trips the overlap trap before its own self-abort** — [low/CONFIRMED/concurrency] `hardware/pov_sync.h:1538`. The fold loop calls `schedule_boundary` (which `HS_CHECK`s the queue drained) *before* `emitter_.tick()` runs its stale-beacon abort. Reachability needs a multi-ms ISR mask, but the mitigation is defeated purely by call ordering. Fix: run the stale-beacon drain before the fold loop schedules new symbols, or drop a provably-stale beacon instead of trapping.
14. ✅ **`classify_faces_impl` builds edge records for degenerate faces without the side-count guard** — [low/CONFIRMED/correctness] `core/mesh/mesh.h:594`. The hash loop guards `count>=3` but the record loop doesn't, so the public `classify_faces_by_topology(PolyMesh&)` overload on an uncompiled mesh self-pairs a 2-gon (silent misclassification) or trips the non-manifold trap far from the cause. Fix: skip `count<3` faces in the record loop too, or `HS_CHECK count>=3` up front.
15. ❌ **Interior-angle quantization uses float round-to-whole-degree — a cross-target fork hazard** — [low/PLAUSIBLE/portability] `core/mesh/mesh.h:561`. `topo_id` hashes `(int)round(angle*180/π)`; a true angle within FP noise of an x.5° boundary can round differently on host/WASM/device, changing which faces cluster and thus which get LUTs — contradicting the "sim and device cannot fork" claim. Trigger is unproven (needs a registry polyhedron near a boundary). Fix: quantize to a coarser tolerance-aware bucket, or document topo_id as advisory. _Rejected: the angle flows through the codebase's own bit-identical `fast_acos` (not libm) under pinned no-fast-math, so the round is identical on every target and the fork is precluded; the trigger stays unproven and a coarser bucket would perturb the perf-sensitive congruence clustering (no sign-off)._
16. ✅ **`Arena` is implicitly copyable, allowing two allocators to alias one buffer** — [low/CONFIRMED/api-design] `core/engine/memory.h:60`. `Arena a = scratch_arena_a;` yields a second allocator over the same buffer with an independent offset, handing out overlapping memory with no trap — defeating the no-overlap guarantee its siblings enforce. No caller does this today. Fix: `=delete` copy ctor/assignment.
17. ✅ **`ArenaVector`'s no-destructor contract is enforceable only by documentation** — [low/CONFIRMED/api-design] `core/engine/memory.h:297`. Storing a type owning out-of-buffer state (`std::function` with a large capture, `std::string`) compiles cleanly and leaks silently on reset. `append_bulk` already static_asserts `is_trivially_copyable`; the destructor contract has no such guard. Fix: `static_assert(is_trivially_destructible_v<T> || is_sanctioned_inplace_function)` in `bind`/`push_back`/`emplace_back`.
18. ❌ **`get_pixel` override requires a paired `overrides_get_pixel()` override with no compile-time enforcement** — [low/PLAUSIBLE/api-design] `core/render/canvas.h:237`. Forgetting the second override compiles clean, renders correctly in the sim, and produces a *wrong device image* (ISR fast paths index the raw buffer). Same forgotten-override footgun class as `needs_full_frame`. Fix: collapse the two into a single CRTP/tag-derived trait, or add a debug sample check `display_buffer()[i]==get_pixel(x,y)`. _Rejected: same forget-proof class as `needs_full_frame`, whose only robust fix is a disproportionate CRTP; the debug sample-check is false-positive-prone (an effect legitimately equal on the sampled pixel while overriding elsewhere) and PLAUSIBLE — no mis-rostered effect exists today._
19. ✅ **Rapid stop→start reuses the previous session's offscreen at stale dimensions** — [low/CONFIRMED/correctness] `daydream/recorder.js:121`. After a stop → resolution-change → start before the async `onstop` fires, `ensurePinnedOffscreen` returns the old-height buffer and the new recording encodes at the wrong size. Fix: reset `this.offscreen` at session boundaries, or validate dimensions before reuse.
20. ✅ **`parseFloat` accepts trailing garbage in numeric URL params and leaves the malformed value in the URL** — [low/CONFIRMED/correctness] `daydream/gui.js:204`. `?Effects.Speed=5px` parses to `5`, and because raw==clamped the URL is never rewritten, so `5px` persists across reloads and masks genuinely wrong links. Fix: use `Number(val)` (strict) so partial-numeric strings fall through to the existing warn path.
21. ❌ **No per-test-case assertion floor: an emptied test loop drops coverage while staying green** — [low/CONFIRMED/test-coverage] `tests/test_harness.h:218`. `end_module` only flags a module whose *entire* tally is zero; a single test function whose asserts stop executing (emptied range-for, early return) is invisible while siblings still assert. The correct guard (`MIN_ASSERTIONS`) exists only for the 2 standalone TUs, not the ~40 in-process modules. Fix: assert a per-module minimum assertion delta. _Rejected: a per-module minimum is arbitrary (each module's count differs by orders of magnitude) and still cannot detect the described failure — an emptied single test-function among asserting siblings keeps the module tally high. The real per-case floor needs a harness redesign (a declared expected count threaded through all ~40 module entry points), out of scope for this fix._

### Priority 3 — Performance, portability, and maintainability

22. ✅ **`ripple_transform` uses two exact libm `expf` per active-ring pixel, beside a deliberate `fast_acos`** — [low/CONFIRMED/performance] `core/engine/transformers.h:353`. Both `expf`s only scale a warp angle already approximated by `fast_acos`; on Cortex-M7 (no HW transcendentals) that's ~100–200 cycles per pixel inside every live ripple. Fix: add a bounded fast-exp (argument is always ≤0) matching the existing `fast_*` convention. *(Do not land without the perf A/B the skill requires.)*
23. ✅ **`Line::sample` degenerate branch leaves progress registers unset** — [low/CONFIRMED/correctness] `core/render/plot.h:883`. A coincident-endpoint line pushes fragments without setting `v0/v1/v2`, so a shader keying on line progress reads stale caller registers (the sibling `Polyline` branch zeroes them). Fix: set `v0=v1=v2=0` in the degenerate branch.
24. ❌ **Planar polylines recompute the azimuthal projection and arc length several times per segment per frame** — [low/PLAUSIBLE/performance] `core/render/plot.h:588`. The pre-pass, rasterizer, and cull each re-project the same endpoints and rebuild the same cumulative-arc table. Cold per-segment (small N), so low impact; a real fix must cache the whole arc table, not just `dist`. Fix: reuse the projected endpoints/table between cull and sampler. **Rejected:** the pre-pass already caches per-segment arc length; caching the full per-segment arc table + projected endpoints spends scarce device persistent-arena bytes (count × PLANAR_LEN_SAMPLES floats) to save cold, low-impact CPU and reworks the seam/arc code findings 1–2 just repaired — a poor, risky trade.
25. ✅ **`deep_tween`'s documented/concept-advertised `Orientation` support is a hard compile error** — [low/CONFIRMED/docs] `core/animation/trails.h:159`. `Orientation<CAP>` satisfies `Tweenable`, so `deep_tween(orientation, cb)` passes the constraint then fails deep inside on `trail.get(i).length()` (a `Quaternion` has no `length()`). Fix: tighten `Tweenable` to reject a non-`Tweenable` element type, or correct the doc to `OrientationTrail`-only (a single `Orientation` is `tween()`'s job).
26. ✅ **Namespace-scope palettes use `static constexpr` (internal linkage) and are address-taken across TUs** — [low/CONFIRMED/maintainability] `core/color/palettes.h:15`. Each TU gets its own copy; `MeshPaletteBank::sources()` (inline, multi-TU) hands out per-TU addresses — duplicated flash, no pointer identity. Fix: `inline constexpr` (C++17+), matching four other core headers.
27. ✅ **Budget-degrade / class-drop / low-quality-class branches in `build_mesh_class_bake` lack targeted tests** — [low/CONFIRMED/test-coverage] `core/mesh/mesh_classes.h:330`. The trickiest arithmetic (shrink recompute, must-not-decrement-budget discard) is unverified; a regression in that accounting passes the census tests. Fix: a synthetic mesh (or tiny `CLASS_LUT_BUDGET`) forcing the shrink/drop and low-quality paths, asserting exact budget accounting.
28. ✅ **`MobiusGrid` over-provisions scratch arenas it barely uses** — [low/CONFIRMED/architecture] `effects/MobiusGrid.h:46`. Reserves 128 KiB scratch (64+64) but the per-frame path uses ~6–7 KiB of `scratch_arena_a` and never touches `scratch_arena_b` — ~120 KiB dead reservation stolen from persistent. Fix: right-size to the per-curve staging need (`scratch_b = 0`), matching MindSplatter/IslamicStars.
29. ✅ **Voronoi render timer excludes the coherence corner-classification pre-pass** — [low/CONFIRMED/performance] `effects/Voronoi.h:188`. `ScopedRenderTimer` is constructed *after* the `~W*H/B²` KD-query pre-pass, so `getRenderUs` under-reports exactly in the dense regime where the pre-pass is costliest. Fix: construct the timer above the classify loop.
30. ✅ **`perf_bench` is `EXCLUDE_FROM_ALL` and not a ctest — zero CI compile coverage** — [low/CONFIRMED/maintainability] `tests/CMakeLists.txt:225`. It pulls the full engine barrel, so an API rename compiles green everywhere while silently breaking `perf_bench.cpp`. Fix: add a build-only CI step (`cmake --build --target perf_bench_*`); running it as a ctest is correctly avoided (noisy wall time).
31. ✅ **URL parse/validate logic duplicated between `daydream.js` bootstrap and `URLSync`** — [low/CONFIRMED/maintainability] `daydream/daydream.js:162`. Two independent encodings of "which effect/resolution values are valid" that can silently diverge. Fix: make `URLSync` the single URL reader; seed `AppState` with plain defaults.
32. ✅ **Device pixel ratio hard-capped at 1 softens rendering on HiDPI displays** — [low/PLAUSIBLE/performance] `daydream/driver.js:155`. `setPixelRatio(min(dpr,1))` renders at CSS-pixel resolution on Retina, softening dots and "Native" recordings, with no comment on whether it's a deliberate fill-rate budget. Fix: document the cap, or raise to `min(dpr,2)`.
33. ✅ **`warmModules` primes URLs against a different base than the Worker load** — [low/CONFIRMED/portability] `daydream/segment_controller.js:77`. Warm-up resolves `./segment_worker.js` against `import.meta.url` while the pool spawns it resolved against the document base; if the module and document base diverge (sub-app mount/bundler), the prewarm caches the wrong entry and becomes a no-op. Fix: resolve both via `new URL('./segment_worker.js', import.meta.url)`.
34. ✅ **`PRNG.nextInt` docstring falsely claims it reproduces device RNG values** — [low/CONFIRMED/docs] `daydream/tools/palette_math.js:141`. The JS preview is a Numerical-Recipes LCG; the device is PCG32 — unrelated streams. The comment contradicts the module's own divergence caveats. Fix: state only that the range *literals* mirror the engine's call sites; the stream is tool-local.
35. ✅ **Fractional sample count silently truncates the exported curve** — [low/PLAUSIBLE/correctness] `daydream/tools/spline_math.js:38`. `numSamples=2.5` drops the endpoint (t never reaches 1). Every other numeric input uses `requireCount` (`Number.isInteger`); this is the gap. Fix: add `|| !Number.isInteger(numSamples)` to both guards.
36. ✅ **`engine_contract_wasm`'s FakeEngine-divergence guarantee is unimplemented; the return contract has already drifted** — [low/CONFIRMED/test-coverage] `daydream/tests/engine_contract_wasm.test.js:5`. The test never compares against the `FakeEngine`; meanwhile the test pins `setEffect`/`setParameter` as boolean while `FakeEngine` returns `undefined`. Harmless only because the worker ignores those returns today. Fix: make `FakeEngine` return the pinned types, or assert the mock's return types against the real engine.
37. ✅ **Dead heavyweight STL includes in the `engine.h` umbrella** — [low/CONFIRMED/maintainability] `core/engine/engine.h:23`. `<list>` (used nowhere) and `<map>` (used only in `test_conway.h`) are pulled into all 27 effects and the no-heap device build. Fix: delete both; have `test_conway.h` include `<map>` directly.
38. ✅ **Circular-buffer iterator advertises random-access category but fails the C++20 concept** — [low/CONFIRMED/portability] `core/engine/static_circular_buffer.h:410`. No default ctor → not semiregular → `std::random_access_iterator` is false, so any ranges algorithm on the buffer fails to compile despite the advertised tag. Fix: add a defaulted default ctor (`m_buffer=nullptr, m_index=0`).
39. ✅ **BZ/GS arena static_asserts bound the resident footprint, not the CubemapLUT build transient peak** — [nit/PLAUSIBLE/maintainability] `effects/BZReactionDiffusion.h:84`. The true peak fits only because `build()` runs before `init_lattice()` *and* the transient happens to equal the node-array size; reordering or growing the scratch would exceed the asserted bound (caught only by the runtime trap). Fix: model the build-time peak explicitly, or static_assert the init ordering.
40. ✅ **`long_way` slerp degenerate (identical-endpoint) path is untested** — [low/CONFIRMED/test-coverage] `tests/test_3dmath.h:767`. Pairs with #6 — the crashing path is uncovered. Fix: add the case (defined valid result once #6 lands, or `HS_EXPECT_DEATH` if the precondition route is chosen).

### Priority 4 — Nits, doc drifts, and micro-cleanups

41. ✅ **`platform.h` pulls `<iostream>` (and its static-init) into every host/WASM TU to back `SerialMock`** — [nit/CONFIRMED/maintainability] `core/engine/platform.h:238`. The rest of the header already formats via `cstdio`. Fix: route `SerialMock` through `printf`/`fputs`; longer term, split the 1,569-line header.
42. ✅ **`construct_in_place` recovery guidance is unachievable after a throwing constructor** — [nit/CONFIRMED/correctness] `core/engine/static_circular_buffer.h:385`. The `@warning` implies post-throw recovery, but every retry re-destroys the dead slot (UB). Unreachable on device (`-fno-exceptions`, trivial T). Fix: state plainly that a throwing element ctor is unsupported.
43. ✅ **`concepts.h` uses `std::addressof`/`std::forward`/`std::nullptr_t` without including their headers** — [nit/CONFIRMED/maintainability] `core/engine/concepts.h:6`. They arrive transitively; a future prune breaks it, against the codebase's self-sufficiency discipline. Fix: add `<utility>`, `<memory>`, `<cstddef>`.
44. ✅ **`Category::Simple` assigned to solids that run heavy `relax(50)`/`snub` pipelines** — [nit/PLAUSIBLE/api-design] `core/mesh/solids.h:1073`. The docstring says Complex marks long pipelines, but several `relax(50)` Archimedeans are tagged Simple. The only consumer is a JS picker label (no runtime cost gate exists). Fix: restate the docstring as registry-membership, or derive the tag from pipeline depth. *Resolved: measured every roster `relax` — the `1e-8` convergence gate fires at 14–39 passes, so the caps (`relax(50/100/217)`) are never reached and these solids are genuinely cheap (the Simple tag is correct); rewrote the `Category` docstring as a display-only label. Separately the `1e-8` gate was ~200× tighter than sub-pixel, so loosened it to `1e-7` (residual <0.05 px, saves a few passes on every relax caller).*
45. ✅ **`AnimationBase::t` documentation misidentifies which animations increment `t` forever** — [low/CONFIRMED/docs] `core/animation/animation.h:252`. Lists `Driver` (duration=1, rewound each frame → bounded) but omits `Noise` and `MobiusWarpEvolving` (genuine `-1` perpetuals with documented 2²⁴ limits). Fix: correct the list to the never-rewound animations.
46. ✅ **`Breakdown` silently maps out-of-range face classes onto the first class** — [nit/PLAUSIBLE/correctness] `core/animation/mesh.h:445`. A mis-declared class count fades every OOB face in `rank[0]`'s window rather than trapping. `Breakdown` is currently unwired. Fix: `HS_CHECK(cls in range)` in `face_offset` (cold per-face), or document the class-count contract. *Fixed: `reorder()` now derives `num_classes` from the mesh's per-face classes (`max+1`) instead of taking a declared count, so it can't be mis-declared; a class at or past `MAX_CLASSES` still folds into `rank[0]` (intended, per the class-count cap).*
47. ✅ **`HD107SFrame` class doc presents the `load()`-only 5-step pipeline as universal** — [nit/CONFIRMED/docs] `hardware/hd107s_frame.h:53`. The shipped `packPixel()` path takes already-linear input and skips step 1 (sRGB→linear). Fix: note step 1 applies only to `load()`.
48. ✅ **`transmitAsync` comment says the completion ISR 'clears' the flag it actually sets** — [nit/CONFIRMED/docs] `hardware/dma_led.h:112`. `dmaISR` stores `true`; `transmitAsync` stores `false`. Fix: reword to "marks `transferComplete_` true."
49. ✅ **Explicitly-bound mesh operators sit outside the `MESHOP_LIST` drift guard** — [nit/CONFIRMED/maintainability] `targets/wasm/wasm.cpp:1309`. `snub`/`hankin`/`relax` are hand-defined and hand-bound; a future explicit op with a forgotten `.function()` line silently never reaches JS. Fix: a second small list macro for the irregular ops, expanded at both sites. *Fixed: added a `MESHOP_IRREGULAR_LIST` name-list macro; the three bindings now expand from it, so a new irregular op is bound the moment it joins the list.*
50. ✅ **`FlowField` particle-emitter lambda is needlessly marked `mutable`** — [nit/CONFIRMED/style] `effects/FlowField.h:68`. Captures only `[this]`; `mutable` is a no-op that mildly misleads. Fix: drop `mutable`.
51. ✅ **`blend_species` carries a concentration-sum guard already enforced by its only caller** — [nit/CONFIRMED/maintainability] `effects/BZReactionDiffusion.h:319`. The inner branch is unreachable in production. Fix: drop it and document the precondition, or note it exists for standalone/test callers.
52. ❌ **Voronoi open-codes `Scan::Shader::draw`'s clip iteration and telemetry** — [nit/CONFIRMED/maintainability] `effects/Voronoi.h:189`. Re-implements the clip band, `check_lut_domain`, and timer inline to reach integer pixel coords; diverges silently if the engine's iteration contract changes. Fix: add a coordinate-exposing `Scan::Shader::draw` variant and route Voronoi through it. *Rejected: not behavior-preserving. `Scan::Shader::draw` iterates the margin-expanded, wrap-aware `x_clip()` band (default `margin=1`); Voronoi deliberately renders the raw `[x_start, x_end)` display band, so routing it through such a variant would add `margin` columns per segment (overlapping writes across the parallel-worker bands). A variant faithful to the raw band would be single-caller-specific with a redundant `check_lut_domain` — a second engine iteration semantics for one caller, disproportionate for a nit. The open-coding is intentional.*
53. ✅ **`run_drc` error-count parsing is fragile and duplicates its REAL-fault classification** — [nit/PLAUSIBLE/maintainability] `hardware/phantasm/gen/analyze_candidates.py:118`. Double regex search + dead `or [0,0]`, and the real-fault classification is duplicated (line 247 vs 319). The headline disqualification does *not* actually trigger (gated on robust counts). Fix: search once; base counts on the `by_type` Counter; factor the classification into one helper.
54. ✅ **Ambiguous single-letter `l` in `linear_to_srgb` (flake8 E741)** — [nit/CONFIRMED/style] `scripts/generate_luts.py:49`. Fix: rename to `v` or `lin` (sibling uses `s`).
55. ✅ **All-black smoke exemption is a hardcoded effect-name string with no pin to the roster** — [nit/CONFIRMED/maintainability] `tests/test_effects.h:107`. `strcmp("RingShower")` — renaming the class silently breaks the exemption (and, with `DEFAULT_FRAMES=8<30`, fails the smoke run). Fix: key off a per-effect trait, or static_assert the string resolves to a registered effect.
56. **Global Space-key handler toggles pause without `preventDefault`, allowing page scroll** — [nit/CONFIRMED/correctness] `daydream/daydream.js:744`. On the `<=900px` scroll layout, Space both pauses and scrolls the sphere out of view. Fix: `e.preventDefault()` for the handled keys.
57. **Default `repointDisplayAliases` contradicts the divergence message and skips the Three.js color buffer** — [nit/CONFIRMED/correctness] `daydream/segment_controller.js:116`. The fallback repoints only `Daydream.pixels`, but the self-heal log claims it also repoints `instanceColor.array`. Production injects the full re-pointer, so this bites only tests/embedders. Fix: make the message accurate, or require the dep (no default).
58. **Per-dispatch `frameSeen` allocation adds avoidable GC churn** — [nit/CONFIRMED/performance] `daydream/segment_controller.js:631`. `new Array(...).fill(false)` every frame; sibling arrays are sized once and `.fill()`-reset. Fix: keep a persistent `this.frameSeen`, reset with `fill(false)`.

---

## Notes for the fix workflow

- Per the `code-review-fix` skill: **validate each finding in a worktree before fixing**,
  one finding = one commit landed FF-only on `master`, mark `N. ✅`/`N. ❌` in this doc in
  the same commit, no finding numbers in code comments, no co-author line.
- **#22 (fast-exp) trades an approximation for speed — get the user's sign-off and an A/B
  measurement before landing.** It is the only finding that touches the performance/accuracy
  contract.
- Several findings are **daydream-repo** (`daydream/...`): those fix on `C:/work/daydream`
  master, with the ✅ flip committed here on Holosphere master.
- The three trail effects (Comets/ChaoticStrings/RingSpin) and the Flyby/Liquid2D shader pair
  each carry duplicated logic the code's own comments flag; no single finding captures the
  duplication, but it is the codebase's most consistent maintainability theme and a good
  candidate for a follow-up refactor.
