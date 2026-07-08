# Holosphere / daydream — Code Quality Review

**Date:** 2026-07-08
**Scope:** The full two-repo product — the Holosphere C++ engine + firmware (`core/`, `effects/`, `hardware/`, `targets/`, `tests/`, `scripts/`, build) and the daydream web simulator (`*.js`, `tools/`, `tests/`, build scripts). Out of scope by request: `effects_legacy.h`, the `*.ino` sketches, `core/math/rotate.h`, and all vendored code (`core/vendor/`, `three.js/`, `node_modules/`).

**Method:** The codebase was partitioned into 29 cohesive subsystem units. An independent reviewer read every in-scope file in each unit (cross-referencing callers where needed) and produced candidate findings graded across all software-quality dimensions. Every candidate finding was then handed to a *separate, skeptical validator* that re-read the cited code and confirmed, adjusted, or rejected it — defaulting to rejection for anything that turned out to be a deliberate design choice, unreachable, or not a concrete fixable defect. Only findings that survived independent validation appear below. 10 candidate findings were rejected in validation and are not listed.

**Headline:** This is an exceptionally engineered, production-grade codebase. Across 29 units and ~60k LOC of engine + simulator code, validation surfaced **zero critical and zero high-severity defects** — only 5 medium and 43 low-severity items, most of them dead code, documentation drift, micro-efficiency, or test-coverage gaps. The non-obvious core math (quaternions, OKLab/OKLCH color, CSG seam handling, the flywheel sync protocol, the relaxed-atomic double buffer) was cross-checked and found correct. The defect density is among the lowest this reviewer has measured in a codebase of this ambition.

---

## Letter Grades by Dimension

| Dimension | Grade | One-line rationale |
|---|:---:|---|
| Correctness & Reliability | **A** | Non-obvious math verified correct throughout; only low-severity edge-case logic gaps (pole-morph tie-break, endpoint tangent, one-shot-timer invariant). No critical/high defects. |
| Architecture & Design | **A+** | Clean layered pipeline (generate → transform → rasterize → filter), CRTP + type-erasure split, compile-time variadic filter chain, explicit-arena discipline — repeatedly rated "excellent/exemplary." |
| Architectural Elegance *(subjective)* | **A+** | Deliberate, load-bearing abstractions with no ornamental complexity; effects compose shared primitives instead of reimplementing transitions. |
| Interface & API Expressiveness *(subjective)* | **A+** | Type-level enforcement of borrow-vs-store lifetimes, `static_assert` pipeline invariants with actionable messages, named factories that close coordinate-space and zero-value footguns. |
| Memory & Resource Safety | **A+** | Wrap-proof subtractive arena arithmetic, use-after-free generation stamps, capacity `static_assert`s, and uniform handling of the WASM ArrayBuffer-detach hazard. One latent unasserted `CANVAS_W ≤ MAX_W` relation noted. |
| Concurrency & Real-Time Correctness | **A+** | 3-index relaxed-atomic ISR double buffer proven correct with bounded watchdog; single-writer flywheel sync model rigorously documented; worker message chain serialized. |
| Performance & Efficiency | **A** | Measured hot/cold split, baked LUTs, branchless folds, `noinline` to dodge register-spill cliffs. Docked for two per-frame/per-pixel hot-path costs (items 2, 4) plus micro-inefficiencies. |
| Robustness & Error Handling | **A** | Principled two-tier policy (fail-fast trap vs. explicit soft-degrade); exemplary clamping at the untrusted JS boundary. A handful of gaps (missing null guard, bfcache freeze, 0-byte download). |
| Portability (Teensy / WASM / desktop) | **A** | Target divergences are deliberate and documented (NDEBUG, integer-only printf, H_OFFSET, clamp NaN contract). Minor latent ODR / section-attribute concerns flagged for awareness. |
| Documentation | **A+** | Consistently "among the best-documented C++ reviewed" — every non-obvious numeric choice, invariant, and rejected alternative is justified inline; Doxygen used as real API reference. |
| Maintainability & Readability | **A** | Strong single-sourcing and X-macro anti-drift. Docked for several dead-code items and two misleading comments (inverted member-ordering rationale, phantom `build_from_flat`). |
| Testing & Verification | **A−** | Deep and adversarial: 38 in-process modules + 39 death cases + roster/include/stack/arena gates, and WASM parity tests that drive the *real* shipped binary against absolute goldens. Docked because it carries the most identified gaps (untested shading math, untested render clock, weak mid-tone color assertion). |
| Build, Tooling & CI | **A−** | Triple-pinned effect roster, sha256 WASM provenance gate, sanitizer host build. Docked for two code-generators whose output diverges from the committed file (breaking "regenerate to fix"), and two perf targets built on every commit. |
| **Overall** | **A / A+** | A mature, meticulously engineered artifact. Every surviving defect is low or medium severity and individually small; the review is a punch-list of polish, not repair. |

---

## Notable Strengths

- **Correctness of the hard parts.** OKLab/OKLCH matrices match Ottosson's canonical constants; quaternion-from-basis signs, stereographic/gnomonic round-trips, and Marsaglia unit-vector generation were all verified. The flywheel sync epoch inference, beacon codec/checksum, and segment-clip math were traced with no off-by-one found.
- **Real-time safety.** The relaxed-atomic double buffer brackets the flip in a compiler barrier, is idempotent once `prev_==next_`, and traps if the writer would alias the ISR's live buffer. The Phantasm sync core enforces a single-writer model with a fused `try_claim()` that cannot open a race.
- **Memory discipline.** No heap on any render path; arena bounds math is written in overflow-proof subtractive form; every carve is guarded by a `static_assert` against the target arena size plus a runtime capacity trap.
- **Anti-drift engineering.** `HS_EFFECT_LIST` / `HS_RESOLUTIONS` X-macros single-source dispatch; the effect roster is triple-pinned (compile-time assert + `--check-modules` + on-disk registration witness); WASM parity tests pin JS ports against the real binary *and* absolute golden values so a shared drift is still caught.
- **Documentation as design record.** Rationale comments explain not just what but why — and why rejected alternatives were rejected — without leaking rename history or review-finding references.

---

## Prioritized Fix List

Every surviving defect is listed exactly once, numbered sequentially. No defect reached Critical or High severity, so those sections are empty and retained only for completeness.

### Critical Priority

*(none)*

### High Priority

1. ✅ **LUT generator emits the wrong include path** — `scripts/generate_luts.py:86` writes `#include "platform.h"` but the committed `core/color/color_luts.h:2` has `#include "engine/platform.h"`; regenerating via the documented command breaks the build, and CI (numeric-token diff only) can't see it. Fix: emit `#include "engine/platform.h"`.
2. **`DistortedRing` runs a 256-sample validation loop per frame** — `core/render/sdf.h:753`. The constructor's `BOUND_SAMPLES=256` `shift_fn` sweep is labeled "once per cold-path construction," but `Scan/Plot::DistortedRing::draw` reconstructs it per ring per frame for Moire/Thrusters. Fix: shrink `BOUND_SAMPLES` to 32–64 (the trap is best-effort anyway) and correct the comment; keep the trap always-on.
3. ✅ **`edge_row_span` under-covers antipodal arcs, can cull a rendered diameter** — `core/render/plot.h:412`. The geodesic clip-cull widens by cross-product axis and skips the extremum block for near-antipodal edges, but the renderer draws a real semicircle about `stable_perpendicular_axis`; a diameter line can vanish on Y-band (Phantasm) boards. Fix: mirror the renderer's axis selection gated on the same `TOLERANCE`, and add an antipodal case to `test_edge_row_span_covers_arc_bulge` built with `stable_perpendicular_axis`.
4. **`MobiusGrid` evaluates palette OKLCH per pixel** — `effects/MobiusGrid.h:266`. `draw_longitudes()`'s shader calls `GenerativePalette::get()` (per-call OKLCH lerp) once per rasterized pixel — up to ~12·H/frame — the exact cost every sibling effect bakes into a `BakedPalette`. Fix: bake into a `BakedPalette` member, sample the LUT, and rebake while the 60-frame `ColorWipe` is in flight (as Comets does).
5. **Tool pages freeze on bfcache back-navigation** — `tools/palettes.html:1275` (and `lissajous.html:161`, `mobius.html:733`, `solids.html:762`/`:793`, `splines.html:424`). All five register a `pagehide` teardown that runs even when `event.persisted` is true, so a Back-button restore resumes with a cancelled RAF / disposed renderer and no re-init. Fix: guard each teardown on `if (!e.persisted)`, matching `daydream.js:790`.

### Medium Priority

6. **Leftover single-step credits fast-forward the sim after unpausing** — `driver.js:410`. Queued Right-arrow `stepFrames` while paused are still spent after pressing space, producing a brief ~60 fps fast-forward. Fix: reset `this.stepFrames = 0` when transitioning to running.
7. **`MeshMorph` pole-symmetry twist checks the wrong mesh** — `core/animation/mesh.h:81`. `has_poles` scans `source` but the twist is applied to `dest` query vertices, so a cube→octahedron morph never gets its tie-break and produces a degenerate vertex correspondence. Fix: detect poles on the dest (query) set (or both).
8. **`Motion::path_frame` roll degenerates at `s≥1` for a buffer Path** — `core/animation/motion.h:296`. The forward finite-difference tangent vanishes at the clamped endpoint, giving the final/seam orientation an arbitrary roll (a per-cycle "up" snap for closed repeating paths). Fix: fall back to a backward difference only when the forward tangent is degenerate.
9. **One-shot timers report `is_finite()==false`, holing the pin-safety guard** — `core/animation/timers.h:25`. A non-repeating timer *is* removed on trigger, but the timeline add-time guard treats it as perpetual, converting a clear add-time diagnostic into a confusing deferred `move_into` trap. Fix: override `is_finite() { return !repeat; }` in `RandomTimer`/`PeriodicTimer`.
10. **`Scan::Shader` two-callback draw lacks the cold null-shader guard** — `core/render/scan.h:959`. Unlike `Plot::rasterize`/`Point`/`Vertices`, it calls type-erased shaders per pixel with no `HS_CHECK`; a null shader becomes a null-thunk indirect call under NDEBUG on device. Fix: add `HS_CHECK` on both shaders once at the top of the draw body.
11. **Empty recording downloads a corrupt 0-byte file** — `recorder.js:357`. The non-File-System-Access sink's `finish()` calls `download()` unconditionally; Record-then-Stop with no captured frames emits a 0-byte unplayable file. Fix: mirror the FSA path's `if (chunks.length)` guard.
12. **Fitted camera distance is not clamped to the controls/far range** — `driver.js:372`. An extreme tall/narrow aspect can push `fittedDistance` past `maxDistance`/`CAMERA_FAR`, framing the sphere too small and grazing the far plane. Fix: `THREE.MathUtils.clamp` to `[minDistance, maxDistance]` before `setLength`.
13. **`getArenaMetrics` failure logs unthrottled per frame, per worker** — `segment_worker.js:211`. A thrown readback emits a fresh `console.warn` every frame from every worker (64–128/s at the Phantasm preset), retaining Error objects and growing devtools memory. Fix: latch the warn with a module-level boolean.
14. **Reaction-graph generator include diverges from the committed file** — `scripts/generate_reaction_graph.py:104` emits bare `#include "reaction_graph.h"` while the committed `.cpp` has `engine/reaction_graph.h`; the numeric-token provenance gate can't detect it, so "regenerate to fix" produces a spurious diff. Fix: emit the `engine/`-relative path.
15. **Tool camera aspect guards `height=0` but not `width=0`** — `tools/shared.js:144` (and `:91`). A zero-width container yields `aspect=0` → a non-finite projection matrix, the exact failure the height guard prevents. Affects all four scene tools. Fix: floor width symmetrically with `Math.max(1, w)`.
16. **Chromium process leaks if browser/page setup throws** — `scripts/capture_screenshots.mjs:79`. Only `chromium.launch()` is inside the try/catch; a throw from `newContext`/`newPage` leaves an orphaned headless process and a raw stack. Fix: wrap the run in `try { … } finally { if (browser) await browser.close(); }`.
17. **`ShapeShifter` orientation members declared after `timeline` with an inverted comment** — `effects/ShapeShifter.h:255`. The three `Orientation` members are destroyed before `~Timeline`, opposite to what their comment claims; benign today only because no animation carries a `.then` hook, but a latent use-after-free and actively misleading. Fix: move the members above `timeline` (matching sibling effects) and correct the comment.

### Low Priority

18. **Unqualified `floor()` forces double promotion in a per-pixel path** — `core/math/geometry.h:340`. `floor(x)/floor(y)` bind to `::floor(double)` in `pixel_to_vector(float,float)`, wasting double work on the render hot loop. Fix: use `std::floor` (float overload), matching line 189.
19. **`CycleCounter::log_node` truncates 64-bit µs to 32-bit** — `core/engine/platform.h:1492`. `(unsigned long)us` re-introduces the wrap the 64-bit accumulator was chosen to avoid, after ~71 min of summed profiled time. Fix: print with `%llu` / `(unsigned long long)`.
20. **`presets.h` pulls heavy `platform.h` for a nonexistent `HS_CHECK`** — `core/engine/presets.h:11`. No `HS_CHECK` is used; only `size_t` is actually needed. Fix: replace with `#include <cstddef>`.
21. **`concepts.h` uses `assert()` without `#include <cassert>`** — `core/engine/concepts.h:162`. Relies on transitive inclusion; a non-self-contained header. Fix: add `#include <cassert>` after the `platform.h` include (NDEBUG ordering).
22. **`inplace_function` cross-type copy/move assignment untested** — `core/engine/inplace_function.h:166`. The destroy-then-construct path between two populated functions holding different closures has no test. Fix: add a case to `test_fn_copy_move_empty` (already in CI).
23. **`FlatView` is dead code with a phantom `build_from_flat` doc** — `core/mesh/mesh.h:177`. Never instantiated; its comment references a symbol that doesn't exist. Fix: delete the struct and comment.
24. **EMSCRIPTEN `Solids::get(index)` is dead code the README presents as active** — `core/mesh/solids.h:1295`. WASM drives solids by name; nothing calls `get(index)`, yet `README.md:1135` describes it as the tooling path. Fix: delete the overload and correct the README to the name-driven path.
25. **`update_hankin` computes `cn` unconditionally per vertex** — `core/mesh/hankin.h:278`. `normalized_or(...)` (a sqrt + divides) is used only on the flat/degenerate fallback branches but runs for every dynamic vertex each frame. Fix: move the computation into the branches that use it.
26. **`shade_blinn_phong` has zero test coverage** — `core/render/shading.h:163`. The metallic model on MorphBlob/Raymarch's hot path (half-Lambert², `^32` specular via five squarings, cubed Fresnel) is unpinned. Fix: add `test_shade_blinn_phong` with independently recomputed golden values.
27. **`shade_mesh_topology` segue overload untested** — `core/render/shading.h:131`. Its cover/grade/opacity routing is exercised only incidentally through rasterizer tests. Fix: add a direct test with stub `PaletteBank`/`SegueT`.
28. **Correction-guard cross-type coexistence trap is not death-tested** — `core/render/led.h:49`. The death suite covers same-type double-construct but not the `NoColorCorrection` + `NoTempCorrection` case the "either type" contract exists to guarantee. Fix: add a cross-type death case.
29. **`pov_single.h` doc claims unsupported 288×144 single-board** — `hardware/pov_single.h:37`. The fixed 12 MHz SPI clock cannot sustain a 288-LED frame within the column period (continuous DMA overrun). Fix: drop the 288×144 claim from the class doc (do not add a clock path for a nonexistent build).
30. **`prev_burst_end_` is never aged against the cycle-counter wrap** — `hardware/pov_sync.h:1382`. Unlike every sibling timestamp, a ~7.16 s master silence can collapse the wrapped difference and misroute the first post-silence symbol (~0.1% window). Fix: clear `have_prev_burst_` in `tick()` once the quiet window elapses.
31. **Dead public default constructor `MeshOpsWrapper()`** — `targets/wasm/wasm.cpp:741`. Never used (embind registers no `.constructor<>()`) and contradicts the no-empty-wrapper invariant. Fix: delete it and its doc comment.
32. **`drawFrame` fast path copies channel-by-channel** — `targets/wasm/wasm.cpp:458`. Three scalar stores per pixel where source/dest byte layouts are identical (`Pixel16` is packed RGB16). Fix: `std::memcpy` the contiguous run with a `sizeof` `static_assert`.
33. **Redundant `SetNoiseType` on the orientation-walk generator** — `effects/ReactionDiffusionBase.h:179`. `RandomWalk`'s ctor already sets type/frequency/seed; the base call has no observable effect. Fix: delete the line.
34. **`Liquid2D` lacks the wrap-invariant test seam its sibling `Flyby` has** — `effects/Liquid2D.h:275`. Four unbounded phase accumulators are unpinned; smoke+determinism can't see a dropped `fmodf`. Fix: add a `Liquid2DWhiteBox` seam and `test_liquid2d_phase_wrapped`.
35. **Dead `isfinite(dt)` guard in `Liquid2D::draw_frame`** — `effects/Liquid2D.h:115`. `dt` is a clamped slider value that can never be non-finite; `Flyby` has no such guard. Fix: delete the two guard lines.
36. **`perf_bench_{Os,O3}` build on every test build with no consumer** — `tests/CMakeLists.txt:223`. Two full-barrel TUs compiled on each `just test` / pre-commit hook, never run automatically. Fix: `EXCLUDE_FROM_ALL` (or gate behind an OFF option); still buildable by target name.
37. **Duplicated segmented-param source + detach guard in `syncGUI` and `export`** — `daydream.js:110`. Two hand-copied copies of the same subtle invariant risk silent divergence. Fix: extract one `liveParamValues()` helper.
38. **Duplicate `timeAccumulator` initialization** — `driver.js:252`. Assigned twice in the constructor with no intervening read. Fix: delete the second assignment.
39. **Main (non-segmented) render loop and its frame clock have no test** — `driver.js:453`. `advanceFrameClock` (backlog clamp, no-accrue-while-paused, one-frame-per-call) is pure over stubbed state but untested — exactly where the item-6 bug lives. Fix: add `tests/driver_clock.test.js` driving `advanceFrameClock` via `prototype.call`.
40. **`quadW`/`quadH` cross the worker boundary but are never read** — `segment_controller.js:314`. Dead data on the structured-clone boundary and dead `FrameResult` fields. Fix: remove them from the protocol, worker post, typedef, store, and the two test assertions.
41. **`computeSegmentRange` validates the Y-band floor but not the X width** — `segment_layout.js:51`. `w<2` collapses both arms to an empty rect and returns silently instead of failing at source. Fix: add `if (w < NUM_ARMS) throw …`, mirroring the height guard.
42. **`generate-importmap.mjs` bakes the version with no exact-semver check** — `scripts/generate-importmap.mjs:32`. A range (`^0.183.1`) from `npm update` would silently lose reproducible pinning and break the CDN URL. Fix: assert an exact `\d+.\d+.\d+` semver after the non-empty check.
43. **`String.replace` on the generated block uses string replacement** — `scripts/generate-importmap.mjs:65`. A `$` in an interpolated value would be treated as a replacement pattern. Fix: use a function replacer `() => block`.
44. **`require-tests` throws raw ENOENT when `tests/` is missing** — `scripts/require-tests.mjs:5`. The directory-removed variant of the guarded scenario dies with a stack instead of the intended message (still exits non-zero). Fix: wrap `readdirSync` to treat ENOENT as an empty match.
45. **`spline_math` sample-count guard admits `Infinity`** — `tools/spline_math.js:38` (and `:63`). `!(n >= 1)` rejects NaN/0/negatives but `Infinity >= 1` passes → infinite sampling loop. Not reachable from the UI slider. Fix: add `|| !Number.isFinite(n)`.
46. **`linearRgbToHex` mid-tone assertion checks only format, not value** — `tests/color.test.js:42`. The lone interior case asserts a `/^#[0-9a-f]{6}$/` shape, so a wrong transfer curve passes (this JS-only fn is not pinned by the WASM parity suite). Fix: assert the exact golden `'#bc89e1'`.
47. **`shared.test.js` only re-tests a re-export** — `tests/shared.test.js:12`. It exercises `showFatalError` (already fully covered by `banner.test.js`) while `initScene()` — `shared.js`'s real export — stays untested. Fix: narrow to a module-load smoke, or cover `initScene`'s DOM-free branches.
48. **`prettify` negative-symbolic branches and `-0.000` normalization untested** — `tests/label_format.test.js:19`. Six negative symbolic constants (`-π/2`, `-φ`, `-√3⁻¹`, …) and the negative-zero guard are unexercised. Fix: add assertions for each.
