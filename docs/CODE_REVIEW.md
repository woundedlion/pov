# Holosphere + Daydream — Code Quality Review

**Scope:** the Holosphere C++ engine/firmware (`core/`, `effects/`, `hardware/`, `targets/`, `tests/`, `scripts/`, `tools/`) and the daydream web simulator (`*.js`, `tools/`, `tests/`).
**Out of scope (excluded by request):** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/rotate.h`, and vendored/generated artifacts (`three.js/`, `node_modules/`, `FastNoiseLite.h` internals, `holosphere_wasm.{js,wasm}`, generated LUT/graph tables).

**Methodology.** 21 component reviewers each read the relevant README sections and every in-scope file in their slice, producing findings tagged by dimension and severity. Every candidate finding was then handed to an independent verifier that re-read the actual code (and the README design rationale) and returned *confirmed / adjusted / rejected*. Only findings that survived verification appear below. Of **98** candidate findings, **73** survived (60 confirmed, 13 severity-adjusted) and **25** were rejected as misreads or deliberate, documented design choices.

---

## 1. Overall Assessment

**Overall grade: A− (high).**

This is a mature, unusually disciplined codebase. The dominant impression across every component — engine, effects, hardware drivers, the WASM boundary, the test suite, and the JS simulator — is *intentionality*: load-bearing design choices are documented with their rationale rather than narrated, invariants are enforced at compile time where possible and trapped fail-fast where not, and the simulator is engineered to be bit-faithful to the device (shared C++ compiled to WASM, integer-exact FastLED mocks, explicit host/device divergence ledgers). The 16-bit linear color pipeline, compile-time-resolution templating, single-block partitioned arena, ISR double-buffer, and `HS_CHECK` fail-fast policy are coherent, consistently applied, and verified by tests (including a death harness that proves the traps actually fire).

The defect profile reflects that maturity: **no Critical issues, one genuine latent memory-safety bug (High), three Medium issues, and 69 Low items that are overwhelmingly documentation drift, minor style/idiom, dead code, and test-coverage gaps** rather than functional defects. The most common real-world risk is *documentation that has drifted from code* — the worst kind, because accurate-looking comments mislead. The single High finding (a use-after-relocation in `Transformer::spawn()`) is latent and untested; it is the one item that should be fixed before the next device deployment.

---

## 2. Quality Dimension Grades

Grades are aggregated across all 21 components (each scored only the dimensions relevant to its code).

| Dimension | Grade | Rationale |
|---|---|---|
| Architecture & design elegance | **A** | Clean layering with no cycles, compile-time-resolution specialization, single-source-of-truth X-macros (effect roster, resolutions), CRTP base for shared effect scaffolding, borrow-vs-own callable wrappers, and erasure pushed to boundaries. Polarity/ownership of composed mesh operators is the one place architecture leans on prose over types. |
| Memory safety & resource management | **A** | Arena bounds math is provably wrap-proof, no hidden heap on the hot path, `std::launder` used correctly, nothrow guarantees enforced for destroy-then-construct paths, WASM-aliased buffers detached before THREE dispose. One latent use-after-relocation (see P1 #1) and a couple of "realloc-only-when-null" gaps in the simulator. |
| Performance & efficiency | **A** | Branchless wrap on the SDF hot path, squared-magnitude comparisons, LUT-based trig/vector reconstruction with documented error bounds, compile-time span budgets. A handful of avoidable O(n) loops where O(1) resets or fill primitives would do. |
| Documentation & comments | **A** | Exceptional intent-focused documentation that explains *why*. Knocked down by the largest single defect class: comments that have drifted from the code they describe (atan2 origin, Gradient rounding, roster-enforcement scope, hankin index bound, pin polarity). |
| Code style & idiom | **A** | Consistent naming, namespace discipline, and formatting. Minor inconsistencies: `std::sqrt` vs `sqrtf`, hand-written copy ops defeating triviality, a mid-class member declaration, duplicated move-transfer blocks. |
| Correctness & bug-freedom | **A−** | Math and protocol logic are correct across quaternions, projections, splines, sync state machines, and segment tiling. Blemishes are contained: a non-repeating timer freezing one aesthetic dimension, a NaN-coerces-to-false path, and a few interfaces that accept out-of-contract inputs without guarding. |
| Maintainability & readability | **A−** | Self-explanatory, centralized constants, clear contracts. Pulled down by verbatim boilerplate duplication (per-effect strobe/persist blocks, duplicated test predicate lambdas) and a little dead code. |
| Error handling & robustness | **A−** | Fail-fast on cold invariant violations, soft handling of genuinely transient conditions, JS boundary treated as untrusted. Gaps: silent mod-256/grayscale/`false` coercions in a few spots that violate stated contracts instead of failing loudly, and a couple of unhandled-rejection paths in CI scripts. |
| Interface expressiveness & API design | **A−** | `explicit` constructors, named factories, `[[nodiscard]]`, type-encoded lifetime contracts. Weak spots: a public mutable `entries` array undercutting an encapsulation invariant, `set_offset` permitting forward jumps, and "any non-X means Y" format fallbacks that should reject. |
| Testing & verification | **B+** | The strongest dimension *relative to peer projects* and the weakest *relative to this codebase's own bar*: a rigorous death harness, WASM↔JS parity tests with golden pins, and host-testable protocol math. But several written tests are never run, key paths (recorder letterbox math, URLSync reset) are untested, and a couple of tests mask the very out-of-range writes they should catch. |

---

## 3. Per-Component Grades

| Component | Repo | Grade | Component | Repo | Grade |
|---|---|---|---|---|---|
| core-infra | Holosphere | A | core-effectsys | Holosphere | A− |
| core-color | Holosphere | A | effects-1 | Holosphere | A− |
| core-rasterizers | Holosphere | A | effects-2 | Holosphere | B+ |
| core-memory | Holosphere | A | core-math | Holosphere | A− |
| effects-3 | Holosphere | A | core-render-pipeline | Holosphere | A− |
| hardware | Holosphere | A | core-mesh-a | Holosphere | A− |
| targets-build | Holosphere | A | core-mesh-b | Holosphere | A− |
| dd-segment-recorder | daydream | A | cpp-tests | Holosphere | A− |
| dd-app-core | daydream | A− | build-scripts | Holosphere | A− |
| dd-ui | daydream | A− | dd-tools | daydream | A− |
| dd-infra-tests | daydream | A− | | | |

---

## 4. Notable Strengths

- **Bit-faithful simulation.** The host FastLED mocks reproduce exact integer `sin8`/`sin16`/`scale8` semantics rather than approximating, and the WASM build shares the device C++ — so the browser preview matches the sphere.
- **Compile-time invariant enforcement.** Scanline span budgets (`sdf_max_spans` + per-combinator `static_assert`), `TrigLUT` quarter-turn offset (`static_assert(W%4==0)`), `MeshFeedback` coprimality (`std::gcd`), and `Liquid2D`'s `sizeof(Params)` guard all reject misuse at compile time instead of trapping at runtime.
- **Fail-fast discipline, verified.** `HS_CHECK` guards cold paths and traps on invariant violation even under `NDEBUG`; hot paths drop to stripped asserts backed by cold traps at the bind site; the `test_death.h` harness proves the traps actually fire (`SIGILL`/`STATUS_ILLEGAL_INSTRUCTION`).
- **Provably-correct arena math.** The subtractive-form bounds checks in `Arena::allocate()` are genuinely wrap-proof, and the persistent/scratch partitioning is explicit at every call site.
- **Host-testable protocol cores.** `pov_sync.h`, `pov_segment_map.h`, `pov_single_map.h`, and `param_marshal.h` split load-bearing arithmetic from the Arduino-only shells so the entire Phantasm sync state machine is unit-tested without hardware.
- **DOM-free logic extraction in the simulator.** `sidebar_logic.js`, `pixel_view.js`, `label_format.js`, `param_sync.js` are pure, single-responsibility, JSDoc'd, and individually tested; the segment-worker generation fence is correct and thoroughly tested including subtle destroy-ordering.
- **Anti-drift tooling.** CI scripts enumerate the effect roster, resolutions, and effect sizes from the X-macro single source of truth rather than hand-maintained lists.

---

## 5. Prioritized List of Items to Fix

Items are numbered sequentially 1–73 across all priority sections, ordered by severity. Each finding cites file `@` location, the quality dimension, and a concrete fix. (`▹` marks findings whose severity was adjusted down from the reviewer's original rating by an independent verifier.)

### P1 — High (fix before next device deployment)

1. **Use-after-relocation in `Transformer::spawn()` completion callback** — `core/transformers.h` @ `spawn_impl` lines 144–174 (callback 168–173) · *Memory safety*. The non-pinned `spawn()` path leaves `TimelineEvent::handled == false`, so the event is freely relocated by `Timeline::step` compaction when any *earlier* event in the process-wide `global_timeline_events` array completes. The completion lambda captures the raw storage pointer `p` by value; after relocation it dereferences the old, destroyed slot (`p->repeats()`), taking a virtual call + member read through reclaimed inline storage. This path has no trap on the fail-fast target and can leave a warp slot stuck active forever or free garbage. The documented use case (multiple `Ripple`/`MobiusWarp` spawns) is exactly the trigger, and the existing test never calls `step()` through completion, so the bug is latent and untested. **Fix:** do not capture the storage pointer; capture only `this` + `idx` and re-query through the live interface, or deactivate the slot keyed by `idx` from the manager/destroy path.

### P2 — Medium

2. **`Moire` palette-wipe timer is non-repeating** ▹ — `effects/Moire.h` @ `init()` line 67; `color_wipe()` 115–131 · *Correctness*. `PeriodicTimer(80, …)` omits the trailing `repeat` arg (defaults `false`), so the timer fires once and cancels itself; the two generative palettes cross-fade exactly once (~frames 80–160) then freeze, contradicting the class doc and `color_wipe()`'s own "called periodically to keep the colors drifting." Sibling effects (`MobiusGrid`, `IslamicStars`, `MeshFeedback`) pass `true` explicitly, confirming this is an omission. **Fix:** pass `repeat=true`: `Animation::PeriodicTimer(80, [this](Canvas&){ color_wipe(); }, true)`.
3. **Written test `test_correct_multifactor()` is never run** — `tests/test_hd107s_frame.h` @ defined line 142; `run_hd107s_tests` 266–276 · *Testing*. This is the *only* test exercising the shipped correction+temperature gain combination together (255,176,240 × 255,147,41 → R=65535, G=26195, B=10121) and the dead-clamp invariant, but `run_hd107s_tests()` never calls it (one definition, zero callers repo-wide). A regression in multi-factor compounding would go uncaught. **Fix:** add `test_correct_multifactor();` to `run_hd107s_tests()` after `test_correct_pipeline()`.
4. **`GenerativePalette` silently produces black/grayscale on an unknown profile** — `tools/palette_math.js` @ ctor, satProfile switch 225–238, brightnessProfile switch 240–265 · *Error handling*. `s1=s2=s3=0` / `v1=v2=v3=0` are pre-initialized and neither switch has a `default`, so a typo'd profile yields a silent grayscale (sat=0) or pure-black (value=0) preview — inconsistent with the same ctor's `gradientShape` throw and the sibling `generativePaletteCpp` which validates via `reject()`. **Fix:** add `default: throw new Error(...)` to both switches (or validate up front against the exported `SATURATION_PROFILES`/`BRIGHTNESS_PROFILES` sets).

### P3 — Low (polish, docs, style, test coverage)

#### core/ — math, infra, color

5. **`fast_atan2` doc says it returns 0 at the origin but returns ~π/2** — `core/3dmath.h` @ 321–341 (comment 322) · *Documentation*. **Fix:** correct the comment to state the bias only prevents a NaN; if a defined 0 is wanted, special-case `abs_x==0 && |y|<=1e-10`.
6. **`Quaternion` copy ctor/assignment hand-written for no functional reason, defeating trivial copyability** — `core/3dmath.h` @ line 407, 429–433 · *Style*. **Fix:** replace both with `= default` to match `Vector`.
7. **Inconsistent sqrt precision: `std::sqrt` (double) vs `sqrtf`** ▹ — `core/3dmath.h` @ `quaternion_from_basis` 1116/1122/1128/1134 · *Style*. **Fix:** use `sqrtf` for consistency with the rest of the file.
8. **`slerp(Vector)` and `slerp(Quaternion)` use different trig backends without note** — `core/3dmath.h` @ 1195–1212 vs 1222–1249 · *Style*. **Fix:** add a one-line comment on each stating why fast-vs-exact trig is intentional, so the divergence is on-record.
9. **`StaticCircularBuffer::clear()` is O(n) where O(1) suffices** — `core/static_circular_buffer.h` @ 204–210 · *Performance*. **Fix:** `head = tail = count = 0;` (no destructors are being skipped).
10. **`construct_in_place` leaves a destroyed slot if `T`'s ctor throws** ▹ — `core/static_circular_buffer.h` @ 381–387 · *Memory safety*. **Fix:** add `static_assert(std::is_nothrow_constructible_v<T, Args&&...>)` (mirroring `inplace_function`) or narrow the docstring claim.
11. **`EVERY_N_SECONDS` multiplies before widening, can overflow the interval** — `core/platform.h` @ 935 (and 927–929) · *Correctness*. **Fix:** `#define EVERY_N_SECONDS(N) EVERY_N_MILLIS((N) * 1000UL)`.
12. **`rand_int` modulo bias is undocumented** — `core/platform.h` @ 1132–1137 · *Documentation*. **Fix:** note the result is modulo-biased for ranges that don't divide 2³², acceptable because callers use small setup-time ranges (matching the `Pcg32` determinism stance).
13. **`GenerativePalette` manual_seed mishandles values < −1** — `core/color.h` @ 1046–1050 · *Error handling*. **Fix:** change the predicate to `if (manual_seed >= 0)` so the auto-seed path covers every sentinel below 0.
14. **Stale `Gradient` doc: claims `static_cast<int>(pos*255)` but code rounds with `+0.5f`** — `core/color.h` @ 895–896 vs 924/936/937/965 · *Documentation*. **Fix:** update the doc to `static_cast<int>(pos * 255 + 0.5f)`.

#### core/ — render pipeline, rasterizers, memory

15. **`Motion::set_duration` accepts negative durations, walking the path backward** — `core/animation.h` @ 1299 · *Error handling*. **Fix:** clamp to `frames < 1 ? 1 : frames` (matching `PeriodicTimer::set_period`) or `HS_CHECK(frames >= 0)`.
16. **`DistortedRing` distortion-bound check is debug-only on device** ▹ — `core/sdf.h` @ ctor 724–737 · *Error handling*. **Fix:** keep the 256-sample scan always-on (it is a cold per-ring cost) or document in `draw()` that an under-bound silently drops geometry.
17. **Duplicated insertion-sort and seam-normalize logic** — `core/scan.h` @ `scan_region` 196–223 vs `sdf.h` `normalize_intervals_to_range`/`sort_intervals_by_start` · *Maintainability*. **Fix:** have `scan_region` call the `SDF::` helpers instead of inlining bit-identical copies.
18. **`Plot::Flower::draw` recomputes `get_antipode` that `sample` already computes** ▹ — `core/plot.h` @ 1811–1817 / 1767 · *Performance*. **Fix:** compute the antipode once and pass it in, or add a comment that both must use identical mapping.
19. **`ScratchScope::make_vector` has no production callers** — `core/memory.h` @ 856 · *Maintainability*. **Fix:** remove it, or note it is currently unused.
20. **`Arena::set_offset` permits a forward jump past the current offset** — `core/memory.h` @ 138–141 · *API design*. **Fix:** tighten to `HS_CHECK(new_offset <= offset)`, or document the forward-seek case and update `high_water_mark`.
21. **`TriangularBitset::data` is uninitialized until `clear()`; reading a pair before `clear()` is UB** ▹ — `core/memory.h` @ 207/212/239 · *Correctness*. **Fix:** document the "must `clear()` first" precondition on the struct, or add an in-class default member initializer.
22. **`ArenaVector` move ctor and move-assign duplicate the full transfer block** — `core/memory.h` @ 351–368 / 375–397 · *Style*. **Fix:** factor into a private `steal_from(ArenaVector&) noexcept` helper.

#### core/ — mesh, effect system

23. **`MeshOps::compile` copies vertices wholesale, leaving orphans that inflate vertex-count consumers** — `core/mesh.h` @ 395–400 · *Correctness*. **Fix:** compact vertices with a used-vertex remap, or expose `num_referenced_vertices()` and pin the behavior with a test.
24. **Comment claims a single shared flash copy of `sort_edge_records`, but internal linkage gives each TU its own** — `core/mesh.h` @ 135–147 · *Documentation*. **Fix:** correct the comment (per-TU COMDAT-folded copies) or give the function external linkage via a single anchor.
25. **Composed-operator return-arena polarity is encoded only in prose** ▹ — `core/conway.h` @ `gyro/meta/needle/zip/bevel` 1118–1184 · *API design*. **Fix:** make polarity legible at the type level (return `{PolyMesh, Arena&}`) or re-clone into `target` to restore primitive polarity.
26. **`compile_hankin` index-range comment states a looser bound than the code checks** — `core/hankin.h` @ 108–114 · *Documentation*. **Fix:** reword so the prose matches the `(I/2)+I` guard and explain the `I-1` underflow it conservatively avoids.
27. **`effects.h` doc overstates roster enforcement — a forgotten `X()` row is silently unshipped on native** — `core/effects.h` @ 43–53 · *Documentation*. **Fix:** scope the guarantee to WASM, or add a native `static_assert` tying the include list to `HS_EFFECT_COUNT`.
28. **`Presets` exposes a public mutable `entries`, undercutting its documented index-invariant encapsulation** — `core/presets.h` @ 102 vs 83–85/77/104–106 · *API design*. **Fix:** make `entries` private with a non-const accessor for the rebind use case, or drop `get_entries()` and document direct access as supported.
29. **`Presets` ctor takes the entry array by value, forcing a full-array copy** — `core/presets.h` @ 43 · *Performance*. **Fix:** `: entries(std::move(e))` (future-proofs non-trivial `Params`).
30. **Stale hardcoded roster count in `kReserveHint` comment** — `core/effect_registry.h` @ 71 · *Documentation*. **Fix:** drop the parenthetical, or `v.reserve(2 * HS_EFFECT_COUNT)` so it tracks the SSOT.
31. **`HS_REG_FILL_PTR` macro left defined in the global preprocessor namespace** — `core/effect_registry.h` @ 175 · *Style*. **Fix:** rename to a clearly-internal `HS_DETAIL_REG_FILL_PTR`.

#### effects/

32. **`FlowField` channel-offset constants named for Y/Z axes but applied to the X coordinate** — `effects/FlowField.h` @ 78–96 · *Documentation/Style*. **Fix:** rename to e.g. `kChannelDecorrelationOffset1/2` so the name stops implying a per-axis offset.
33. **Repeated `strobe_columns()`/`persist_pixels` boilerplate duplicated verbatim across effects** — `effects/{ChaoticStrings,Comets,DistortedRing,DreamBalls,Dynamo,FlowField,Flyby}.h` + `ReactionDiffusionBase.h` · *Maintainability*. **Fix:** a thin intermediate base or defaulted member for the `strobe=true` family; at minimum stop copy-pasting the rationale comment.
34. **`init_lattice()` ordering contract enforced only by prose** ▹ — `effects/ReactionDiffusionBase.h` @ 134–148 · *API design*. **Fix:** `HS_CHECK` that the persistent arena is configured, or fold derived-state allocation into a base-driven template-method init.
35. **`IslamicStars` ripple-pool `static_assert` under-models the staggered burst peak** — `effects/IslamicStars.h` @ 90–93; `ripple()` 116–122 · *Testing*. **Fix:** tighten the comment (overflow degrades to dropped ripples, not corruption) or use a ceil-based bound accounting for the stagger window.
36. **`GnomonicStars` passes a hardcoded `0.0f` time to `fib_spiral` with no explanation** ▹ — `effects/GnomonicStars.h` @ 106 · *Documentation*. **Fix:** name the argument inline (`// phase=0: base lattice is animation-independent`) or hoist a named constant.
37. **`Voronoi` has an unreachable empty-buffer guard after the per-site rotation loop** — `effects/Voronoi.h` @ 89 · *Maintainability*. **Fix:** remove the dead guard, or move it above the loop as a documented invariant.
38. **`PetalFlow` recomputes `move_dist` identically in two methods that must stay in lockstep** — `effects/PetalFlow.h` @ 145 / 187 · *Style*. **Fix:** extract `float move_dist() const { return params.speed * RHO_PER_SPEED; }`.
39. **`PetalFlow::next_hue` declared mid-class between member functions** — `effects/PetalFlow.h` @ 159 · *Style*. **Fix:** move the declaration up with the other data members.

#### hardware/, targets/, build

40. **`read_id()` docstring misstates pin polarity** — `hardware/pov_segmented.h` @ 335–352 · *Documentation*. **Fix:** "grounding a pin reads LOW; the raw reading is inverted so a grounded pin contributes a 1, and all-floating = ID 0 (master)."
41. **`bakeLut` silently mod-256-truncates out-of-range h/s/v ints** — `targets/wasm/wasm.cpp` @ 1162–1175 · *Error handling*. **Fix:** clamp to [0,255] with a log (matching `gradientShape`) or document the intentional wrap; the doxygen currently says values are "in [0,255]".
42. **Tooling `Arena` globals have external linkage, inconsistent with sibling internal-linkage state** — `targets/wasm/wasm.cpp` @ 76/81/82 · *Style*. **Fix:** mark the three `Arena` objects `static` or move module-private globals into an anonymous namespace.
43. **Dangling sentence fragment in `setEffect` inline comment** — `targets/wasm/wasm.cpp` @ 340 · *Documentation*. **Fix:** rewrite as a complete sentence explaining the `%s`-not-format-string guard.
44. **`drawFrame` no-effect clear uses a scalar loop instead of a fill primitive** — `targets/wasm/wasm.cpp` @ 429–431 · *Performance*. **Fix:** `std::fill_n(pixelBuffer.data(), pixel_width * pixel_height * kChannels, uint16_t{0});`.

#### tests (C++)

45. **POSIX death-harness spawn breaks if `argv[0]` contains a quote or shell metacharacter** — `tests/test_death.h` @ 831–834 · *Error handling*. **Fix:** use `posix_spawn`/`fork`+`execv` with `open("/dev/null")`+`dup2` (mirroring the Windows `_spawnv` path), eliminating the shell.
46. **Two-probe trap-shape agreement can falsely classify under a forked shell if `cs[1]` does not trap** — `tests/test_death.h` @ 974–985 · *Testing*. **Fix:** use a dedicated always-trapping sentinel for shape detection and still run `cs[0]`/`cs[1]` through the normal per-case loop.
47. **`check_tiling` allocates a `w*ROWS` cover grid but only touches two columns, masking out-of-range `x_col` writes** — `tests/test_pov_segmented.h` @ 51–90 · *Testing*. **Fix:** guard the write `if (x_col < 0 || x_col >= w) continue;` after the `HS_EXPECT`.
48. **Redundant `global_timeline_t` reset obscures the canonical reset path** ▹ — `tests/test_param_marshal.h` @ 47–49/122–124/165–167 (also `perf_bench.cpp:27`, `arena_measure.cpp:29`) · *Style*. **Fix:** use `reset_globals()`; drop the redundant `global_timeline_t = 0`.
49. **Large boot-join / pre-train predicate lambdas duplicated across `pov_sync` scenarios** — `tests/test_pov_sync.h` @ multiple (e.g. 992–999, 1061–1068, 1839–1846) · *Maintainability*. **Fix:** promote the `boot_join`/`to_pre_train` helpers to file scope and reuse across scenarios.

#### scripts & Python tooling

50. **`check_screenshots.mjs` throws an unhandled rejection if `docs/screenshots/` is absent** — `scripts/check_screenshots.mjs` @ 16 · *Error handling*. **Fix:** try/catch the `readdir`; treat ENOENT as "all roster effects missing" with the same `::error::` guidance and non-zero exit.
51. **JSONC stripper silently truncates the rest of the file on an unterminated `/*`** — `tools/teensy_gate.py` @ 344–347 · *Error handling*. **Fix:** raise `ValueError('unterminated block comment in JSONC')` if no closing `*/` was found.
52. **`region_totals_from_size_a` computes `region_for_address(addr)` twice per row** — `tools/teensy_gate.py` @ 162 · *Style*. **Fix:** hoist `r = region_for_address(addr)` (or use `defaultdict(int)`).
53. **`wasm_smoke.mjs` uses `process.exit()` after console output, risking truncated logs** — `scripts/wasm_smoke.mjs` @ 46–47/68–69/321–322 · *Error handling*. **Fix:** prefer `process.exitCode = 1; return;` so buffered output flushes.

#### daydream — app core

54. **Uncaught alias-invariant throw inside the animation loop permanently halts rendering** — `daydream.js` @ `adapter.drawFrame()` 491–497 (reached via the `setAnimationLoop` callback) · *Error handling*. **Fix:** degrade per the file's own policy — `console.error` and re-point the aliases, or set a one-shot fatal overlay like the context-loss path, rather than throwing every frame.
55. **`testAll` auto-advance silently no-ops at a resolution whose effect list is missing** — `daydream.js` @ 573–588 · *Correctness*. **Fix:** mirror `applyResolution`'s fallback: `const currentList = effectsByResolution[appState.get('resolution')] || HiResFavorites;`.
56. **`instanceColor` buffer reallocated only when null, masking a count/size mismatch** ▹ — `driver.js` @ `precomputeMatrices()` 722–733 · *Memory safety*. **Fix:** reallocate when `array.length !== count*3`, or assert the lengths match.
57. **`disposeApp` does not destroy the global GUI or the active-effect GUI** — `daydream.js` @ 717–728 · *Memory safety*. **Fix:** destroy `guiInstance` and `activeEffect.gui` (mirroring `applyEffect`'s teardown, including draining `activeDragEnds` and removing the window pointer listeners).
58. **`AppState.update` batch-then-notify can re-enter `set()` and interleave notifications nondeterministically** — `state.js` @ `update()` 74–84, `_notify` 105–107 · *Architecture*. **Fix:** document the re-entrancy hazard, or snapshot the notification queue and skip a queued tuple whose key has since changed.

#### daydream — UI

59. **JSDoc block for `_createSortBtn` is detached from its method by an intervening method** — `sidebar.js` @ 174–207 · *Documentation*. **Fix:** move the JSDoc (174–181) down to immediately precede `_createSortBtn` at 195.
60. **NaN engine value for a boolean controller writes a spurious `false`** — `param_sync.js` @ `resolveParamSync` 36–40 · *Correctness*. **Fix:** `if (Number.isNaN(incoming)) return { update: false, value };` before coercion.
61. **Out-of-range numeric deep link leaves a stale value in the URL after clamping** — `gui.js` @ `add()` 210–218 / 248–251 · *Error handling*. **Fix:** call `this._urlWriter(key, val)` once after attaching when the hydrated value was clamped/snapped.
62. **Mouse-click selection does not move the roving tabindex** — `sidebar.js` @ `setEffects()` 122 · *Correctness*. **Fix:** `btn.onclick = () => { this._setRovingTabbable(btn); this.onSelect(name); };`.
63. **Unused mutable local `precision` in `prettify`** — `label_format.js` @ 22/40 · *Style*. **Fix:** inline the literal `r.toFixed(3)` or hoist a module-level `const PRECISION = 3;`.

#### daydream — segment workers & recorder

64. **`init` handler commits `canvasW/canvasH` to a rejected resolution before the success check** ▹ — `segment_worker.js` @ `handleMessage 'init'` 71–81 · *Correctness*. **Fix:** defer the `canvasW/canvasH` assignment until after the `setResolution` success check, mirroring the `setResolution` handler.
65. **`captureFrame` letterbox/scaling math has no unit test** — `recorder.js` @ `captureFrame` 189–211 · *Testing*. **Fix:** add a test recording a known source aspect into a differently-shaped offscreen, spy the `drawImage` args, and assert the letterbox math for both aspect directions.
66. **Boundary-overlay seam set rebuilt with two passes and `Set` allocations every composite call** — `segment_controller.js` @ `composite()` 497–525 · *Performance*. **Fix:** cache the boundary line sets when the layout is established and invalidate on generation bump.

#### daydream — tools

67. **`splineExportCode` treats any non-`'vectors'` format as `'fragments'` instead of rejecting unknowns** — `tools/spline_math.js` @ 91–111 · *API design*. **Fix:** validate `format` against an explicit allow-set and throw on anything else (matching `solid_codegen.js`/`palette_math.js`).
68. **`GenerativePalette.getChannelValue` recomputes all three channels per call** ▹ — `tools/palette_math.js` @ 365–367 · *Performance*. **Fix:** have plotting callers call `get(t)` once and index the triple, or document the minor cost.

#### daydream — infra & tests

69. **`URLSync.reset()` and the `setParam` deletion-marker path are entirely untested** — `tests/state.test.js` @ URLSync section 80–140 · *Testing*. **Fix:** add tests for `setParam(k, null)` removal, `reset(['keep'])` preservation, and `reset()` cancelling a pending debounced flush.
70. **`composite()` out-of-bounds test asserts a guarantee the code only gives for a leading bad segment** — `tests/segment_controller.test.js` @ 421–433 · *Testing*. **Fix:** weaken the comment, or add a two-result case (good segment then overflowing segment) and assert the partial-blit-then-throw behavior.
71. **`generate-importmap` assumes `pkg.dependencies.three`/`lil-gui` exist; opaque crash if relocated** — `scripts/generate-importmap.mjs` @ 27–30 · *Error handling*. **Fix:** validate both versions are non-empty strings and `console.error` + `process.exit(1)` with a named message otherwise.
72. **Redundant `defer` on a `type=module` script** — `index.html` @ 74 · *Style*. **Fix:** drop `defer`; module scripts already defer.
73. **`vendor-importmap` `EXTRA` can silently clobber core entries like `'three'`** ▹ — `vendor-importmap.js` @ 59–64 · *API design*. **Fix:** merge base last (`Object.assign({}, EXTRA, base)`) or warn on collision with a known core key; otherwise document the override as intentional.

---

## 6. Rejected Candidate Findings

25 candidate findings were rejected by independent verification — most because they flagged a deliberate, README-documented design choice (fail-fast traps, compile-time resolution, arena ownership, relaxed-atomics double-buffer gating, determinism-over-uniformity RNG) as a defect, or misread the code. They are excluded here to keep the fix list actionable and evidence-backed.
