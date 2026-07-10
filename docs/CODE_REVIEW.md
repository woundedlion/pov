# Holosphere + daydream — Code Quality Review

**Scope:** Both repositories that ship as one product — the Holosphere C++ rendering
engine + Teensy firmware + Emscripten/WASM target, and the daydream web simulator
(Three.js renderer, segmented-POV workers, geometry tools). ~200 first-party
source files.

**Out of scope (excluded from review):** `core/engine/effects_legacy.h`,
`targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and all vendored code
(`core/vendor/*`, `.pio/*`, `daydream/three.js/*`, `daydream/node_modules/*`).

**Method:** The codebase was partitioned into 26 component slices, each audited by
an independent reviewer that read every in-scope file in its slice against the
architecture in `README.md`. Every candidate finding was then re-checked by a
separate validator that opened the cited code with fresh, skeptical eyes; 13
candidate findings were rejected at that gate and are not reported here. What
remains is 56 validated findings.

---

## Executive Summary

This is an exceptionally well-engineered codebase. Across ~200 first-party files
spanning hard-real-time firmware, a templated rendering engine, a multi-board
sync protocol, and a browser simulator, the independent audit surfaced **zero
critical and zero high-severity defects on any live path**. The 56 findings are
7 medium and 49 low: latent guards that are currently unreachable, dead code the
project's own conventions say to delete, documentation that has drifted a few
lines from the code, and test-coverage gaps. Nothing here produces wrong LED
output, corruption, or a crash on a path exercised today.

The engineering philosophy is coherent and rigorously executed: compile-time
`<W,H>` resolution specialization for zero-overhead generality; a single 330 KiB
arena with explicit `Arena&` plumbing and no render-path heap; a fail-fast
`HS_CHECK` doctrine that traps invariant violations at the violation site rather
than shipping garbage to the LEDs; 16-bit linear color end-to-end with OKLCH
perceptual palette interpolation; and a variadic filter pipeline that lifts
coordinates between world/screen/pixel domains at compile time and makes illegal
filter stacks a build error. The multi-board 1-wire flywheel sync protocol is
designed to *fail to "missed," never to "wrong."* Documentation is reference-grade
and, where sampled against the code, accurate.

The honest weak spot is **test breadth**, not test quality: several delicate
primitives and decision cores are exercised only indirectly, and a handful of
guards are soft where they should be hard. These are cheap to close.

**Overall grade: A−** (composite of the per-dimension grades sits at the A / A−
boundary; held one notch below a flat A to reflect the test-breadth gaps and a
recurring tendency toward comment over-documentation — not any correctness doubt).

---

## Grades by Dimension

| Dimension | Grade | One-line rationale |
|---|---|---|
| Correctness & reliability | **A** | No live-path defect found across 26 slices; degeneracies (poles, seams, antipodes, NaN) handled explicitly and reasoned in comments. |
| Architecture & design elegance | **A** | Compile-time resolution, single-arena discipline, variadic domain-lifting pipeline, fail-fast doctrine — coherent and consistently applied (unanimous A across reviewers). |
| Interface & API expressiveness | **A** | Ownership/lifetime encoded in the type system (`Arena&` everywhere, owned-vs-borrowed vectors, trap-vs-graceful `normalized()`/`normalized_or()`); a few API asymmetries and a generic global `generate()`. |
| Readability & code style | **A** | Uniform naming, terse fact-focused rationale comments; density trends high but earns its place in bit-identical-sim/hardware code. |
| Documentation | **A** | The 2168-line README is reference-grade and matches the code where sampled; in-file doxygen is exhaustive. Minor drift in a few §9/§10 entries. |
| Test coverage | **A−** | Broad host suite + death harness + WASM runtime smoke + cross-run determinism; gaps are real (untested `fast_expf`, `screen_step`, Euler oracle on Archimedean/Catalan, sidebar DOM lifecycle, `daydream.js` sequencing) and a few soft guards. |
| Performance & efficiency | **A** | Zero-overhead by construction; measured hot-path tuning (split-trig LUTs, congruence-class fast path, downsampled cached feedback); zero render-path heap (unanimous A). |
| Memory & resource safety | **A** | Wrap-proof arena math, debug-only use-after-free instrumentation, `ScratchScope` LIFO, `Persist<T>` watermark checks; the few gaps are latent and unreachable by current callers. |
| Portability & build | **A** | Clean Teensy/host/WASM tri-split, compile-time index-width `static_assert`s, double-precision LUT/lattice provenance CI-gated; a couple of host-only parity nits. |
| Error handling & robustness | **A** | Fail-fast `HS_CHECK` on cold seams, device-stripped `assert` on hot paths, out-of-process death tests that prove the traps fire; a few missing finite guards. |
| Maintainability | **A** | X-macro single-sources-of-truth with anti-drift guards, shared skeletons preventing operator/loop drift; drags are dead surface and comment volume. |

---

## Rationale

### Correctness & reliability — A
After independent validation, not one finding is a bug on a path exercised today.
The reviewers probed the hard cases — stereographic pole caps, gnomonic-equator
floors, slerp near-parallel/antipodal branches, Conway output-pool sizing derived
from the closed-manifold identity `E = I/2`, the feedback warp seam-unwrap, and
the double-buffer index ownership — and found them correct and, unusually, reasoned
in the comments. The residual correctness findings are latent: `Segue::hash01` can
return exactly `1.0f` against its `[0,1)` contract (F25), a congruence-class census
double-counts on a documented-but-unreached rebake path (F12), and a `DeepLinkGUI`
key can collide for identically-named sibling folders that the shipping app never
creates (F49). All have downstream clamps or single safe callers today.

### Architecture & design elegance — A (unanimous)
Every reviewer independently graded this A. The five stated philosophies —
compile-time resolution, arena allocation, the ISR double buffer, fail-fast
`HS_CHECK`, and 16-bit linear color — are not slogans; they are visible in the
structure of every subsystem. Standouts: the variadic `Pipeline<W,H,Filters...>`
with compile-time World/Screen/Pixel domain lifting and trait `static_assert`s that
make an illegal filter stack fail to build; clip-culling routed through the same
world transforms the renderer applies, so segmented rendering culls by *rendered*
latitude; the congruence-class LUT facility architected to degrade to the exact
path at every seam so it is never load-bearing; and the 1-wire sync protocol's
odd-only distance-2 symbol alphabet that degrades a glitch to a missed symbol, never
a misclassified one.

### Interface & API expressiveness — A
Types carry intent: `Arena&` on every allocating call site (no hidden state),
move-only `ArenaVector` vs borrowed `ArenaSpan`, `normalized()` (traps) vs
`normalized_or()` (graceful edge), `spawn`/`spawn_pinned` encoding the retained-handle
contract, and a trusted (`get_by_name`, fail-fast) vs untrusted (`has_name`) registry
boundary the WASM target actually honors. The blemishes are minor and were caught by
the review itself: a converting-constructor exclusion narrow enough to let a derived
wrapper double-erase (F10), `generate()` living in the global namespace under a very
generic name, and the composed-Conway-operator arena polarity that a direct caller
must internalize.

### Readability & code style — A
Naming is uniform (snake_case, with `noiseGenerator` the one caught exception, F45),
structure is idiomatic, and the rationale comments explain the *why* of every epsilon
and branch. The recurring critique across reviewers is comment *over-*density — several
function docblocks are longer than their bodies — but the comments are accurate and, in
code where a subtle numeric drift ships to physical LEDs, the verbosity is defensible.
No finding-reference or rename-history comments were found (consistent with project
conventions).

### Documentation — A
The README is genuinely a reference document: the frame-sync datasheet alone (§7.10)
specifies symbol waveforms, AC timing, and failure modes to a professional hardware
standard, and it matches `hardware/pov_sync.h`. In-file doxygen is exhaustive. The
deductions are a handful of drifted claims the review pinned to line numbers: README
§7.6 overstates "all palette interpolation in OKLCH" (Gradient is linear-light, F38),
the §9 Raymarch entry names the wrong solid and parameter (F40), MeshFeedback's
parameter line is missing (F41), and §10.2 omits half the WASM bridge's export surface
(F42). All are doc-only edits.

### Test coverage — A− (the lone sub-A dimension)
The suite is broad and multi-layered: a native Clang unit build, an effect smoke
harness that renders every effect at 288×144 with asserts on, a cross-run determinism
diff, an out-of-process death harness that proves each `HS_CHECK` traps, the sync-protocol
core, HD107S wire-format tests, POV tiling proofs, and a CI WASM runtime smoke that drives
every effect at every resolution. The gaps are specific and enumerated: `fast_expf` is the
one `fast_*` primitive with no accuracy test though it runs per-pixel (F35); `screen_step`'s
adaptive-density formula has no direct pin (F36); the Euler oracle is not run over the
Archimedean/Catalan registries (F37); `sidebar.js`'s DOM-lifecycle/leak-prevention contract
is untested (F52); `daydream.js`'s effect/resolution sequencing is unextracted and unexercised
(F54); and two guards are soft where they should halt (`face_topo_record`'s `count<=24`, F9;
the `MIN_DEATH_CASES` floor sits one below the roster, F49).

### Performance & efficiency — A (unanimous)
Zero-overhead is achieved by construction, not by tuning after the fact: compile-time
`<W,H>` collapses every coordinate transform and LUT index, bump allocation is pointer
advancement, and all debug bookkeeping compiles out under `NDEBUG`. The measured hot-path
work is careful — split-trig `TrigLUT` recovering `cos` from a quarter-turn tail, the
congruence-class LUT fast path with an exact-walk fallback, the feedback warp field
downsampled and cached on a pure-function key, ARM `uqadd16` saturating adds, `fast_*`
transcendentals sized to sub-pixel error. The only inefficiencies found are cold and
trivial (a doubled per-frame `prepare_thresholds`, F30; a redundant stored phase array,
F24; an `O(n²)` classifier in an offline PCB tool, F32).

### Memory & resource safety — A
This is the codebase's core competency. Arena OOM and capacity math is written in
subtractive, wrap-proof form so a colossal size can't slip past the bounds check;
debug builds stamp generation + per-vector rebind counters to turn use-after-free and
source-vector re-grow into asserts at zero release cost; `ScratchScope` enforces LIFO
teardown and `Persist<T>` verifies a post-restore watermark. The findings are latent:
`allocate_n` omits the `size_t`-overflow guard its sibling `bind` has (F3, unreachable
by current small counts), `MobiusWarp::bind_scale` stores a `const&` with no rvalue
protection against a temporary (F1), and DreamBalls strands ~5 KB of raw preset meshes
in the persistent arena where the IslamicStars idiom would reclaim them (F2).

### Portability & build — A
The Teensy / native-host / WASM split is handled deliberately and single-sourced through
`platform.h`, with compile-time index-width `static_assert`s that fail the build rather
than silently truncate on a budget bump, and double-precision provenance (the reaction-graph
lattice folded in `double` and narrowed once) that is CI-gated so the runtime table and its
Python generator cannot diverge. CMake presets drive WASM/tests uniformly. The nits are
host-only: `map()`'s `long` intermediate differs on LP64 vs the device's 32-bit `long`
(F56), and two PCB scripts hardcode the Windows KiCad path (F47).

### Error handling & robustness — A
Fail-fast is applied with judgment, exactly per the stated doctrine: always-on `HS_CHECK`
traps on cold setup/bounds seams (arena OOM, capacity, non-manifold edges, config
over-subscription), device-stripped `assert` on the per-pixel hot path, and genuine
transient conditions (DMA overrun, dropped frame) given bounded/soft handling instead.
The death harness proves the traps actually fire. The gaps are missing finite guards on
one Möbius warp variant (F4) and a browser label formatter that renders `NaN` as literal
text (F14).

### Maintainability — A
`HS_RESOLUTIONS` and `HS_EFFECT_LIST` X-macros give single sources of truth with
compile-time and startup anti-drift guards; shared skeletons (`emit_vertex_orbit_faces`,
`draw_fragments`, `fill_edge_record`, one `rasterize` core) structurally prevent parallel
loops from drifting. The drags are the dead surface the review catalogs (below) and the
comment volume.

---

## Prioritized Fix List

Every validated defect, numbered sequentially. No critical or high-severity items were
found; **Priority 1** collects the latent correctness/safety/robustness defects worth
fixing first (a change to surrounding code or an unlucky input turns several of them live),
**Priority 2** is dead code and simplification the project's own "delete over document"
rule targets, and **Priority 3** is documentation, test, tooling, and naming hygiene.

### Priority 1 — Correctness, safety & robustness

1. `core/animation/params.h:471` — `MobiusWarp::bind_scale(const float&)` stores the address of a parameter that binds to a temporary; add `void bind_scale(const float&&) = delete;` (matches the file's Lerp/MobiusFlow discipline).
2. `effects/DreamBalls.h:214` — each preset strands its raw `PolyMesh` (~5 KB total) in the persistent arena; build the raw solid into scratch via `generate(persistent, [&](target,a,b){…})` like IslamicStars, and fix the misleading "does not outlive this loop" comment.
3. `core/engine/memory.h:126` — `allocate_n<T>(n)` computes `n*sizeof(T)` with no overflow guard; add `HS_CHECK(n <= SIZE_MAX / sizeof(T), …)` mirroring `ArenaVector::bind`.
4. `core/animation/params.h:563` — `MobiusWarpEvolving` guards no finiteness; add `HS_CHECK(std::isfinite(scale) && std::isfinite(speed))` in the ctor and clamp-to-last-good in `set_speed`/`set_scale` (siblings already do this).
5. `core/animation/mesh.h:247` — `Segue::hash01` can return exactly `1.0f`, breaking its `[0,1)` contract; use the top 24 bits: `float(h >> 8) * (1.0f/16777216.0f)`.
6. `effects/ShapeShifter.h:204` — `SphericalPolygon` is routed through `default:` in both dispatch switches, suppressing `-Wswitch`; use an explicit `case` so a future `ShapeType` fails to compile until wired.
7. `effects/Thrusters.h:34` / `effects/RingShower.h:31` — hardcode `full_frame` instead of `decltype(filters)::any_crosses_segments`; derive it (as ShapeShifter/SplineFlow do) so a later segment-crossing filter can't silently tear on Phantasm.
8. `core/mesh/mesh_classes.h:282` — `build_mesh_class_bake` census counters accumulate with `+=` and never reset at entry; zero the five scalars so a documented in-place rebake doesn't double-count telemetry.
9. `tests/test_mesh.h:554` — `face_topo_record`'s `HS_EXPECT_TRUE(count<=24)` doesn't halt, so a future >24-side face overruns `int angles[24]`; make it a hard early-out (or size the array to the true max + `static_assert`).
10. `core/engine/concepts.h:125` — `FunctionRef` converting ctors exclude only the exact type, so a `StoredFunctionRef` lvalue double-erases; change the exclusion to `!std::is_base_of_v<FunctionRef, std::decay_t<Callable>>` on both ctors.
11. `daydream/segment_controller.js:564` — `setEffect()` leaves `this.animationsPaused` stale, so a pool rebuilt after an effect switch re-pauses while the main sphere animates; reset it to `false` alongside `this.paramValues = null`.
12. `daydream/driver.js:616` — the picture-in-picture camera copies the main camera's position+orientation, rendering a duplicate of the front instead of the opposite hemisphere; reflect the camera through the orbit target and `lookAt(target)` (back view valid only with `cullBackSphere` off).
13. `daydream/gui.js:145` — `DeepLinkGUI.getKey` collides deep-link keys for sibling folders sharing a non-empty name; disambiguate only genuine duplicates (or fail-fast on duplicate names) — do **not** change the key format unconditionally (breaks existing shared links).
14. `daydream/label_format.js:39` — `prettify` renders `NaN`/`Infinity` as literal text on sphere labels; add `if (!Number.isFinite(r)) return "—";`.
15. `daydream/tools/slider.js:94` — the initial readout uses the raw value while the thumb uses the rounded value, so first paint can disagree; derive the readout from the same rounded value.

### Priority 2 — Dead code, simplification & efficiency

16. ✅ `core/math/geometry.h:391` — the log-polar trio (`LogPolar`, `LOGPOLAR_RHO_SENTINEL`, `logPolarToVector`, `vectorToLogPolar`) has zero production callers; delete it and its tests (it also carries a latent over-wide pole cap).
17. ✅ `daydream/tools/palette_math.js:165` — `hsvToRgb` + `CPixel` are dead production code kept alive only by a WASM parity test; delete them and drop the orphaned assertions. (Dead JS mirror deleted; the two mirror-dependent parity tests replaced with a golden absolute-pin on the live `hsv_to_rgb` WASM export so its coverage — otherwise untested — is retained.)
18. ✅ `core/mesh/solids.h:61` — `finalize_solid` copies `temp.topology`, which no generator ever populates (it's filled later by `classify_faces_by_topology`); drop the two topology lines.
19. ✅ `core/color/color.h:943` — the `Pixel lerp_oklch(CPixel,CPixel,float)` overload has no callers; delete it.
20. ✅ `core/render/sdf.h:515` — per-shape `nx`/`nz` (and Star/Flower `scan_nx`/`scan_nz`) are stored but never read; drop the dead members and their construction-time assignments.
21. ✅ `core/render/plot.h:93` — `Plot::Point` and `Plot::Vertices` are unreferenced in both repos; delete them (and their README rows).
22. ✅ `core/render/shading.h:204` — `NullVertexShader`/`NullFragmentShader` are referenced only by their own tests; wire them in as canonical no-op defaults or delete both.
23. ❌ `core/mesh/solids.h:1213` — `Collections::get_archimedean_solids()` is a dead accessor (its four siblings are live); delete it (and its README mention). **Rejected:** kept as the parity accessor to its four live siblings; a future effect may need it.
24. ✅ `effects/HopfFibration.h:184` — `fiber_phase[i]` is always `i * PHASE_STEP`; delete the 840-byte array + init loop and compute the phase inline.
25. ✅ `hardware/pov_sync.h:1535` — `master_on_crossing` takes an unused `TickActions&`; drop the parameter and update the single call site.
26. ✅ `targets/wasm/wasm.cpp:461` — dead store `idx += count*CHANNELS;` on the memcpy fast path (and `&pixelBuffer[idx]` where `idx==0`); delete the store and use `pixelBuffer.data()`.
27. ✅ `daydream/daydream.js:8` — `SLOW_FRAME_MS` is imported but never used; drop it from the import.
28. ✅ `daydream/driver.js:206` — `this.resources = []` is never read or written; delete the field.
29. ✅ `core/engine/transformers.h:346` + `core/animation/params.h:640` — the ripple half-width floor is duplicated across two files with a load-bearing coupling; extract a single `RippleParams::half_width()` accessor used by both.
30. ✅ `core/engine/transformers.h:234` — `prepare_thresholds()` runs twice per frame per ripple (`Transformer::prepare_frame` and `Ripple::step`); drop the redundant call from the ripple path.
31. ✅ `core/math/3dmath.h:19` — an unused `#include "vendor/FastNoiseLite.h"` in the universal math header; delete the line (consumers get it from `engine.h`/`animation.h`).
32. ✅ `hardware/phantasm/gen/analyze_candidates.py:141` — arc-vs-segment classification is `O(n²)` via `if b in arcs`; tag the kind while iterating instead.
33. ✅ `effects/GnomonicStars.h:50` — add the persistent-footprint `static_assert` its sibling effects use, so a future `MAX_POINTS` bump fails at compile time rather than trapping at runtime.
34. ✅ `hardware/phantasm/gen/shorts.py:79` — the geometric short-checker ignores T-junctions and `(junction …)` records, so it can't fulfill its stated guarantee (false "clean"); union wire endpoints colinear on other spans, or narrow the docstring and defer to `check.py`.

### Priority 3 — Documentation, tests, tooling & naming

35. ✅ `core/math/3dmath.h:1209` — `fast_expf` (a per-pixel bit-hack approximation) has no accuracy test while every sibling `fast_*` does; add a dense-grid relative-error sweep vs `expf` plus boundary cases.
36. ✅ `core/render/plot.h:495` — `screen_step`'s adaptive-density / pole-clamp formula has no direct regression pin; add a unit test comparing against the analytic `SCREEN_STEP_PX/|v_screen|` in an unclamped regime (sign is squared into speed — test magnitude, not sign).
37. ✅ `tests/test_solids.h:214` — the Euler oracle (`V−E+F==2`) runs only on Platonic + Islamic slices; extend it over the Archimedean and Catalan registries.
38. `core/color/color.h:1061` / `README.md:1008` — README §7.6's "all palette interpolation in OKLCH" is false for `Gradient` (linear-light); qualify the claim.
39. `core/render/scan.h:844` — `ssaa_sample_vector`'s doc claims byte-for-byte parity with `pixel_to_vector`, but the two use different trig paths; correct the comment to state the actual (harmless) invariant.
40. `README.md:1785` — the §9 Raymarch entry names "20 vertices of a dodecahedron / Core Size"; the code builds a 26-vertex disdyakis dodecahedron with a "Fill" param. Fix both.
41. `README.md:1715` — the §9 MeshFeedback entry omits its **Parameters** line (Fade, Distort Amp/Freq/Speed, Noise Scale, Hue Shift, Feedback); add it.
42. `README.md:1901` — §10.2 omits the bridge's eight color/palette/lissajous free exports and the `PaletteOps`/`bakeLut` class; document them.
43. `core/mesh/mesh_classes.h:335` — a comment references the nonexistent `kClassLutMaxN`; use the real `CLASS_LUT_MAX_N`.
44. `effects/HankinSolids.h:36` — the comment says the H=144 arena high-water "isn't in CI"; it is (`test_hankinsolids_arena_budget_covers_every_solid`). Fix the stale comment.
45. `core/animation/motion.h:606` — rename the camelCase `noiseGenerator` member to `noise_generator` to match the file's snake_case convention.
46. ✅ `targets/wasm/wasm.cpp:1186` — `bakeLut`'s boundary clamps are inline and host-untested, breaking the extract-pure-logic pattern this scope establishes; move them into a `wasm_predicates.h`-style header with tests.
47. `hardware/phantasm/gen/pcb.py:19` / `check.py:11` — hardcode the Windows KiCad 10.0 path with no env override; reuse `fab.find_kicad_cli()` / `KICAD_CLI`.
48. ✅ `tests/CMakeLists.txt:224` — `perf_bench` is `EXCLUDE_FROM_ALL` with no test wiring, so it can bit-rot silently; keep one opt level in the default build as a compile guard.
49. ✅ `tests/test_death.h:1175` — `MIN_DEATH_CASES = 40` sits one below the 41 registered cases, defeating the anti-drop floor for a single deletion; set it to 41.
50. ✅ `tests/CMakeLists.txt:143` — `unit_color_luts_check` vanishes with no skip marker when Python is absent (unlike its sibling); register it unconditionally and route the no-Python case through `SKIP_RETURN_CODE 77`.
51. ✅ `tests/test_platform.h:260` — `test_serial_printf` leaks its temp file + `FILE*` on the `dup` failure path and uses a fixed CWD-relative name; close/remove before the early return and use a unique name.
52. `daydream/sidebar.js:83` — the DOM-lifecycle/leak-prevention methods (`dispose`/`setEffects`/`applySortOrder`/`setActive`/`updateScrollArrows`) are untested; add a fake-stub test (no jsdom — the runner is plain `node --test`).
53. `daydream/tests/geometry.test.js:21` — `pixelToSpherical` is validated only against a co-located re-derivation of the same formula; add 2–3 hardcoded golden `[x,y,z]` pins from an independent source (as the other perceptual ports do).
54. `daydream/daydream.js:417` — the 790-line app orchestrator has no test; extract the pure `applyEffect`/`applyResolution` sequencing + skew-guard predicates into a testable module (mirroring `resolveParamSync`).
55. `core/engine/static_circular_buffer.h:195` — `pop()` removes the *front* element, asymmetric next to `push_front`/`push_back`/`pop_back`; rename to `pop_front()` and update the two call sites.
56. `core/engine/platform.h:601` — the host `map()` mock computes its intermediate in `long`, which is 64-bit on LP64 hosts vs 32-bit on device/Win32, so a wide-range `map()` diverges sim-vs-device; compute in an explicit 32-bit type or document the legacy-only gap (as `beat88` does).
