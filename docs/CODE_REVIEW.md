# Holosphere — Code Quality Review

**Scope:** the full C++ rendering engine (`core/`, `effects/`, `hardware/`, `targets/wasm/`), the
build/CI/tooling surface (`CMakeLists.txt`, `platformio.ini`, `scripts/`, `.github/workflows/`), the
native test suite (`tests/`), and the sibling **daydream** web simulator (Three.js + WASM front-end).
Out of scope by request: `effects_legacy.h`, `*.ino` sketches, and `core/math/rotate.h`.

**Method:** eleven component reviewers read every in-scope file against the README's stated design
philosophy, generated candidate findings across every quality dimension, and then dispatched *fresh,
independent* validator sub-agents to confirm each finding against the actual code before it was
admitted. A finding was kept only if a second agent (or a direct re-read) confirmed it was real and
not an intentional, documented design choice. The overwhelming majority of candidate defects were
**refuted** — the code was arithmetically correct or a deliberate, documented decision consistent with
the engine's fail-fast / arena / compile-time-resolution philosophies.

**Headline:** this is exceptional, production-grade embedded C++. Across a deliberately adversarial
audit spanning correctness, memory safety, concurrency/ISR-safety, UB, and numerical robustness, **not
a single validated correctness, memory-safety, or undefined-behavior defect survived.** The residue is
fourteen low-severity maintainability, comment-accuracy, CI-hardening, and test-coverage items. The
dominant risk in this codebase is *over-documentation density*, not incorrectness.

---

## Overall Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | **A** | Zero validated correctness bugs after adversarial validation of ~90 candidate defects spanning quaternion/OKLCH/Möbius math, reaction-diffusion integration, rasterization intervals, Conway operators, and the sync protocol. Every singularity (poles, antipodes, 0/0, denormals) is deliberately guarded. |
| **Memory & Resource Safety** | **A** | Bump-arena bounds math is wrap-proof; every scratch carve is `ScratchScope`-rewound and pre-bound by a compile-time `static_assert` on the real 330 KiB budget; narrowing casts trap; the WASM zero-copy view-detachment hazard is defended in three independent places. |
| **Concurrency / ISR-Safety** | **A** | Strict single-writer model on Teensy; relaxed atomics justified per single-observer Cortex-M7 with the one true cross-thread edge (`pending_effect_`) correctly using release/acquire; the worker pool uses a generation fence + bounded watchdogs. |
| **API / Interface Design & Expressiveness** | **A−** | Borrow-vs-store enforced by type (`FunctionRef`/`StoredFunctionRef`), `Arena&`-everywhere explicitness, X-macro single-source rosters, designated-initializer aggregates. Minor: a duck-typed segue-policy hook whose arity diverges silently, and `transform`'s output silently carrying empty topology (a latent footgun). |
| **Architectural Elegance** | **A** | Compile-time `<W,H>` specialization, the variadic `Pipeline` with compile-time domain routing, the three-arena partition, CRTP animation/effect bases, and pure host-testable cores split cleanly from Arduino/DOM shells. Coherent and deliberate throughout. |
| **Performance** | **A** | Hot loops move OKLCH lerps, Q-format conversions, and trig-heavy builds off the per-pixel path; branchless ISR; noinline register-spill avoidance; LUT-domain checks hoisted per-draw not per-pixel; fast/exact-trig dichotomy applied consistently. |
| **Portability** | **A** | Meticulous host↔device parity (FastLED integer mocks, LP64 wrap analysis, DMAMEM-via-explicit-specialization for the GCC vague-linkage attr-drop), target-agnostic `platform.h`. |
| **Maintainability** | **A−** | Load-bearing invariants are named and pinned by `static_assert`; shared scaffolds prevent operator/effect drift. Costs: very dense files, hand-maintained cross-effect trail coupling, and a couple of stale comments. |
| **Documentation & Comments** | **A−** | Comments explain *why*, not *what*, and follow house style (no history, no finding-refs). Blemishes: two stale comments (a `~75%` gate figure, a wrong uniform name) and occasional essay-length rationale blocks. |
| **Test Suite & Testability** | **A−** | Assertions pin *correct* behavior against independent oracles (double-precision refs, Ottosson OKLab triples, brute-force KNN, golden noise hashes); the roster is triple-pinned; a death harness verifies traps fire. A few leaf render helpers and the `<96,20>` determinism diff are thin. |
| **Build / CI Robustness** | **A** | CI compiles both Teensy images + WASM (release *and* debug) + a sharded native suite with ASan/UBSan and a real Windows run, and verifies the deploy provenance trio in-job. |

**Composite: A / A−.** A mature, unusually self-aware codebase where nearly every place a reviewer
would expect a race, overflow, or protocol-timing hazard already carries a precise invariant comment,
a `static_assert`, or an `HS_CHECK` trap.

---

## Per-Component Grades

### core/engine — platform, arena, callables, rosters, styles, transformers, reaction graph

| Dimension | Grade |
|---|---|
| Correctness | A |
| Memory Safety | A |
| API/Interface | A |
| Architecture | A |
| Performance | A |
| Portability | A− |
| Maintainability | A |
| Comments/Docs | A+ |
| Testability | A− |

No validated defects. Every candidate (signed overflow in `beatsin16`, div-by-zero in `sync_hue`,
integer division in `node()`, unsigned wrap in `Presets::prev()`, arena padding wrap) was independently
refuted as correct or documented-intentional. Strongest comment discipline in the tree.

### core/math + core/color — 3dmath, geometry, easing, waves, color, palettes

| Dimension | Grade |
|---|---|
| Correctness | A |
| Numerical Robustness | A |
| API/Interface | A |
| Architecture | A |
| Performance | A |
| Portability | A |
| Maintainability | A |
| Comments/Docs | A |
| Testability | A |

No validated defects across **20** candidates. Notable near-misses that proved correct: `fast_cbrt`'s
Halley iteration, `fast_atan2`'s four-quadrant sign/range, `hue_rotate`'s unit-renormalization making
the OKLab rotation exactly chroma-preserving, and `Complex::operator/`'s denominator-zero sentinel.

### core/mesh — mesh, mesh_classes, spatial, conway, hankin, solids

| Dimension | Grade |
|---|---|
| Correctness | A |
| Memory Safety | A |
| API/Interface | A− |
| Architecture | A |
| Performance | A |
| Maintainability | A |
| Comments/Docs | B+ |
| Testability | A |

One LOW finding (stale `~75%` gate comment vs the real `0.4` constant). Conway operator pool sizing,
KDTree/k-NN correctness, arena polarity, and index narrowing all validated exact.

### core/render — canvas, scan, plot, filter, sdf, shading, led

| Dimension | Grade |
|---|---|
| Correctness | A |
| Memory Safety | A |
| API/Interface | A |
| Architecture | A |
| Performance | A |
| Maintainability | A |
| Comments/Docs | A− |
| Testability | A− |

No validated defects across **20** candidates covering interval seam-splitting, CSG span bounds, AA
ramps, filter-pipeline domain conversions, feedback bilinear sampling, and pole/antipode guards.

### core/animation — animation, timers, params, motion, trails, sprites, timeline, mesh

| Dimension | Grade |
|---|---|
| Correctness | A |
| Lifetime/Memory Safety | A |
| API/Interface | A− |
| Architecture | A |
| Performance | A |
| Maintainability | A− |
| Comments/Docs | A |
| Testability | B+ |

Two LOW findings (sub-frame fade-in truncation; duck-typed segue-hook arity divergence). The motion
co-driver/repeat-seam/underflow and timeline gap-fill/pin/overflow paths were all refuted, several
pinned by existing tests.

### effects (group A) — BZ/GS reaction-diffusion, ChaoticStrings, Comets, DistortedRing, DreamBalls, Dynamo, FlowField, Flyby, GnomonicStars, HankinSolids, HopfFibration, IslamicStars, aged_slot

| Dimension | Grade |
|---|---|
| Correctness | A |
| Memory Safety | A |
| API/Interface Use | A |
| Architecture | A |
| Performance | A |
| Maintainability | A− |
| Comments/Docs | A |
| Consistency | A− |

No validated defects across **14** candidates (RD physics, orbit integration, Hopf S³ lift, arena
reuse, index bounds). Only soft note: deliberately hand-maintained trail coupling between sibling
effects.

### effects (group B) — Liquid2D, MeshFeedback, MindSplatter, MobiusGrid, Moire, PetalFlow, Raymarch, RingShower, RingSpin, ShapeShifter, SphericalHarmonics, SplineFlow, Thrusters, Voronoi

| Dimension | Grade |
|---|---|
| Correctness | A− |
| Memory Safety | A |
| API/Interface Use | A− |
| Architecture | A− |
| Performance | A− |
| Maintainability | B+ |
| Comments/Docs | A |
| Consistency | A− |

Two LOW + one INFO finding (duplicated `kMaxRings`; raw `allocate` vs `allocate_n` in RingSpin;
`MeshState` slot-clear idiom inconsistency). The one investigated correctness candidate (Voronoi
corner-grid overrun) was refuted — indexing is provably in-bounds.

### hardware — dma_led(_core), hd107s_frame, pov_single(_map), pov_segment_map, pov_sync, pov_segmented

| Dimension | Grade |
|---|---|
| Correctness | A |
| Concurrency/ISR-Safety | A |
| API/Interface | A− |
| Architecture | A |
| Performance | A |
| Portability | A− |
| Maintainability | A |
| Comments/Docs | A |
| Testability | A |

No validated defects across a full concurrency/atomics/timing audit. Color-correction overflow bounds,
eDMA end-frame `ceil(N/16)` latch, flywheel `position()` int32 window, and the relaxed-vs-release
atomic choices all carry proven-correct invariants and `static_assert`s.

### targets + build + scripts — wasm bridge, CMake/PlatformIO/toolchain, CI, generators

| Dimension | Grade |
|---|---|
| Correctness | A |
| Build/CI Robustness | A |
| WASM Bridge Design | A |
| Tooling Quality | A− |
| Portability | A |
| Maintainability | A− |
| Comments/Docs | B+ |

Three LOW + one INFO finding (dead `-?` regex branch; CI-only roster scripts lacking local npm entries;
un-smoked debug WASM; uncaught `FileNotFoundError` on a bogus `CLANG_FORMAT` override). All correctness
candidates (param-view lifetime, `>32` params, finiteness single-sourcing, `getPixels` overrun) were
refuted.

### tests — native host harness

| Dimension | Grade |
|---|---|
| Coverage | A |
| Assertion Quality | A |
| Correctness | A |
| Determinism/Robustness | A |
| CI Integration | A |
| Maintainability | A− |

No defects in the tests themselves. Coverage-gap observations: `shading.h` has no dedicated module;
`ChromaticShift`/`Blur` are not isolated; the cross-run determinism diff runs only at 288×144, not the
`<96,20>` device specialization (the ISR/DMA shell in `dma_led.h` being host-untestable is correct
factoring, not a gap).

### daydream (web simulator) — Three.js + WASM front-end

| Dimension | Grade |
|---|---|
| Correctness | A |
| Resource/Memory Safety | A |
| Concurrency | A |
| API/Interface | A |
| Architecture | A |
| Maintainability | A− |
| Comments/Docs | A− |
| Testability | A |

One LOW finding (stale uniform name in a JSDoc comment). The zero-copy `typed_memory_view` detachment
hazard is defended three ways; the segmented-worker pipeline has a rigorous generation fence with
bounded watchdogs; every resource owner has a symmetric bfcache-aware `dispose()`.

---

## Prioritized Fix List

Every validated item is listed below under its priority band, numbered sequentially. Severity reflects
that **no correctness, memory-safety, or UB defect was found** — the list is maintainability,
comment-accuracy, CI-hardening, and test-coverage work.

### Priority 1 — Critical / High

None. No high-severity or critical defect was identified in scope.

### Priority 2 — Medium

1. ✅ **Unify the `Segue` policy hook signature** — `core/animation/mesh.h:392` vs `232,336,355`. `Segue::Breakdown::face_offset(const Vector&, int, int cls)` takes three params while `Base`/`TerminatorSweep`/`Shockwave` take two; because policies are duck-typed (no `virtual`/`override`), the divergence is invisible until a call site instantiates against the specific policy, and nothing structurally guarantees a uniform contract. Give all policies one signature (pass the class index everywhere, defaulting `cls` in `Base`) so the contract is compiler-checkable.

2. ✅ **Correct the stale `~75%` gate figure** — `core/mesh/mesh_classes.h:107,356`. Two comments describing the class-LUT retention gate say "~75%", but the actual constant is `kMinClassHitShare = 0.4f` (a 40% floor). Restate the comments as ~40% / reference `kMinClassHitShare` so a future tuner isn't misled into thinking the bar is nearly double what it is.

3. ✅ **Fix the stale uniform name in the strobe doc comment** — `daydream/driver.js:569`. The `setStrobeColumns` JSDoc says the mode is applied via `uColumnFillScale`, but no such uniform exists — the real uniform is `uColumnFillArc` (a repo-wide grep for `uColumnFillScale` returns zero hits). Rename in the comment.

4. ✅ **Runtime-smoke the debug WASM build in CI** — `.github/workflows/ci.yml:260-267`. `wasm-debug` is compiled (with `-sASSERTIONS=1`) but never run; `wasm_smoke.mjs` already accepts a `WASM_JS` override, so an assertion-guarded regression on a debug-only branch currently rides a green build. Add one `node scripts/wasm_smoke.mjs build/wasm-debug/holosphere_wasm.js` step.

### Priority 3 — Low / Informational

5. ❌ **De-duplicate the `kMaxRings` constant** — `effects/Moire.h:134` and `effects/ShapeShifter.h:32`. Rejected: these are different effects with independent per-effect raster budgets that happen to share a formula today; coupling them through a shared helper would falsely tie two unrelated tuning knobs together.

6. ✅ **Route the RingSpin ring pool through `Arena::allocate_n`** — `effects/RingSpin.h:73-74`. The pool is allocated as raw bytes + placement-new instead of the typed `allocate_n<Ring>` idiom introduced to centralize exactly this pattern. No correctness impact (persistent arena, never freed); consistency only.

7. ✅ **Add a one-line note on sub-frame fade-in truncation** — `core/animation/sprites.h:49-53`. When `fade_total > duration`, the integer rescale `(long long)duration * fade_in / fade_total` floors to 0 at tiny durations (e.g. `duration=1`), dropping the fade-in. `fade_out` stays non-negative, so there is no UB — this is acceptable degradation; document the truncation rather than change behavior.

8. ✅ **Drop the dead `-?` branch in the reaction-graph provenance regex** — `.github/workflows/ci.yml:382`. Fibonacci-lattice neighbor indices are always in `[0, 7679]`, so `grep -oE '-?[0-9]+'` never matches a sign; simplify to `[0-9]+` (matching the sibling LUT job) so the pattern doesn't imply signed data that doesn't exist.

9. ✅ **Give the roster-check scripts local npm entries** — `package.json` / `scripts/check_effect_roster.mjs` / `scripts/check_screenshots.mjs`. These run only from the `screenshot-gallery` CI job (so they are wired into CI, not orphaned), but a developer cannot cheaply reproduce the gate before pushing. Add `check-roster` / `check-screenshots` npm scripts for parity with `screenshots` and `wasm-smoke`.

10. ✅ **Pick one carousel slot-clear idiom** — `effects/MeshFeedback.h:103,211`. `init()` uses `carousel.slot(...).clear()` while `start_morph()` uses `carousel.slot(new_slot) = MeshState();`; the two achieve the same "empty before compile" goal by different means and would diverge if `MeshState` assignment and `clear()` ever differ (e.g. one retains capacity). Standardize on one.

11. ✅ **Handle a bogus explicit `CLANG_FORMAT` override gracefully** — `scripts/generate_luts.py:115`. `subprocess.run([cf, ...])` raises an uncaught `FileNotFoundError` if `CLANG_FORMAT` points at a non-existent path (the auto-detect path degrades with a warning). Arguably fail-loud-is-correct for a power-user override; catch and emit a clear message if a friendlier failure is wanted.

12. ❌ **Add isolated pixel-output tests for `Pixel::ChromaticShift` and `Screen::Blur`** — `tests/`. Rejected: already handled. `tests/test_filter.h` has direct functional tests asserting the emitted taps/pixel output on known inputs — `test_blur_factor_zero_is_identity`, `test_blur_full_kernel_sums_to_alpha`, `test_blur_update_changes_kernel`, `test_blur_pole_row_renormalizes`, and `test_chromatic_shift_fanout`. The premise ("only transitively through effect smoke tests") is not accurate.

13. **Add a dedicated `shading.h` test module** — `tests/`. `Fragment`/null-shader helpers are covered only incidentally through rasterizer tests; a small isolated module would pin their contract directly.

14. **Extend the cross-run determinism frame-diff to the `<96,20>` device resolution** — `tests/`. The byte-identical determinism oracle runs only at 288×144; the `<96,20>` specialization gets a smoke pass but no frame-diff, so `PhiLUT<20>`/`H_OFFSET`-specific nondeterminism would not be caught by the determinism oracle.
