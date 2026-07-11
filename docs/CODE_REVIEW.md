# Holosphere + daydream — Code Quality Review

**Scope:** The Holosphere C++ engine and firmware (`c:/work/Holosphere`) and its WebAssembly browser simulator, daydream (`c:/work/daydream`), reviewed as one product.
**Out of scope:** `effects_legacy`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and the compiled `holosphere_wasm.{js,wasm}` artifacts.
**Method:** 26 independent Opus reviewers each owned one subsystem, read every file in their scope against the README architecture, gathered candidate defects, and validated each one (re-reading the cited code and its callers) before inclusion. Findings that dissolved on inspection, restated a deliberate design decision, or were purely stylistic were discarded. Every surviving finding is concrete, reachable, and fixable with a minimal, non-performance-regressing change.

---

## Executive Summary

**Overall grade: A.**

This is an exceptional codebase. Across ~91k lines of C++ and ~19k lines of JavaScript, the validated defect list is **34 findings, every one low-severity** — no correctness bug, no undefined behavior on a live path, no race, no resource leak survived validation. For a real-time, no-heap firmware engine that also cross-compiles to WebAssembly and drives a distributed four-microcontroller display, that defect density is remarkable and is itself the headline result: the surface has already been swept thoroughly enough that what remains is dead code, missing test pins, duplication, and doc drift.

The engineering culture is coherent from the platform layer to the browser: a single deterministic sim/device parity invariant (down to bit-identical PRNG and integer math), a fail-fast trap discipline (`HS_CHECK` on cold seams, stripped asserts on hot paths, both proven to fire by a SIGILL death harness in CI), arena-only allocation on every render path, and a test suite engineered specifically against false-green tests via independent oracles and non-vacuity guards. Documentation is publication-grade and lives at the point of use.

The weak spots are narrow and consistent: **interfaces** and **testing** are the lowest dimensions (both A−), held back by ergonomic burdens that rest on comments/cold-traps rather than the type system, and by a handful of honest coverage gaps (untested cold Persist paths, a few hand-derived numeric branches with only smoke coverage, one vacuous assertion). The PCB-generation Python tooling is the least-polished corner (its `testing` grade is C+) but is peripheral to the shipped board's correctness.

---

## Dimension Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | **A** | No validated correctness bug or live-path UB across either repo. Overflow/wrap math is uniformly in subtractive wrap-proof form; every div-by-zero and singularity seam (stereographic poles, degenerate cross products, Möbius cancellation, factorial precision) is identified and guarded rather than hoped away. The only latent items are a legacy-path `map()` division edge and a per-session teardown race behind a modal dialog — both low-reachability. |
| **Architecture** | **A / A+** | First-rate structural design: the three-arena single-block allocator with `configure_arenas()` repartitioning; the compile-time-resolution templating that specializes the whole pipeline per `<W,H>` at zero runtime cost; the variadic filter `Pipeline` that auto-inserts domain conversions and turns ordering mistakes into compile errors; the type-erased `PipelineRef` as a code-size lever; the phase-driven `Segue` policy system that vaporizes unused shading paths. Elegance is engineered, not incidental. |
| **Interfaces** | **A−** | Expressive where it counts — the `FunctionRef`/`StoredFunctionRef`/`inplace_function` trio encodes lifetime intent in the type system, concepts are tight, the SDF/`Fragment` ABI is shared cleanly across five rasterizer families. Marked down because several load-bearing contracts (arena-polarity for Conway composition, borrowed-vs-owned `MeshState`, unit-position preconditions on fast paths, the `DistanceResult` overloaded-register convention) demand out-of-band knowledge that a stronger type-level encoding could enforce; the enum-param-as-float convention and hand-tuned callable byte budgets add friction. |
| **Readability** | **A** | Doxygen on nearly every member with the load-bearing *why* spelled out (NaN→hi clamp operand order, determinism/reentrancy contracts, destructor-skipping contract). Consistent naming and idiom. The one recurring cost is comment density: a fraction of the rationale essays border on excessive and slow a first read — but given a no-debugger hardware target and genuinely subtle invariants, the trade is defensible and largely earns its place. |
| **Testing** | **A−** | The suites are the strongest most reviewers had seen in an embedded/graphics codebase: independent oracles (Newell normal, Euler characteristic, azimuthal-equidistant arc reconstruction), non-vacuity guards that assert the interesting case actually fired, a self-probing death harness, a multi-board sync event simulator, device-only hot paths reconstructed in flag-specific TUs, and a cross-repo layout crosscheck. Held to A− by honest gaps: untested cold `clone_from`/`step_wipe_rebake` Persist paths, a few hand-derived effect branches (MobiusGrid poles) with only smoke coverage, one vacuous assertion, and the peripheral PCB tooling being largely ungated. |
| **Performance** | **A** | Microsecond-budget discipline throughout: branchless clamp with a documented architecture-specific NaN contract, `always_inline` on hot scalar helpers, LUT-driven trig, arena-only render paths, hot/cold code placement tuned to the Teensy ITCM granule cliff, and profiling that compiles to nothing unless enabled. The reaction-diffusion solvers fit a 7,680-node Laplacian in a 330 KiB arena at POV frame rates. The few efficiency findings are all off the per-pixel path (per-layer palette re-eval, per-segment reprojection). |
| **Error handling** | **A** | The fail-fast philosophy is applied with unusual discipline and, crucially, consistency: cold invariant seams trap (survive `NDEBUG`, no stdio), hot paths use stripped asserts backed by cold-site traps, and genuinely transient conditions (dropped frames, full pools) get bounded soft handling. The code reliably distinguishes "programmer error, crash now" from "expected saturation, degrade gracefully." `check_fail` formats into a fixed stack buffer so it is safe under OOM. |
| **Documentation** | **A / A+** | Exceptional and internally consistent: a 2,000-line README that matches the code, spec documents for the frame-sync protocol / PCB / CI gate, and inline contracts that state invariants and host/device divergences with their rationale. The only drift found is cosmetic (palette-name casing, one hard-pinned V/F/I count, an inaccurate factorial-precision note). |

### Per-component grade matrix

Grades are the reviewing agent's assessment of each subsystem on the same eight dimensions (`-` = not observable in that scope).

| Component | Corr | Arch | Intf | Read | Test | Perf | Err | Doc |
|---|---|---|---|---|---|---|---|---|
| engine-core | A | A+ | A | A | A− | A+ | A+ | A+ |
| engine-support | A | A | A | A− | A | A | A | A |
| math | A | A | A | A | A+ | A | A | A |
| mesh-topology | A | A | A− | A | A− | A | A | A |
| mesh-patterns-solids | A | A | A− | A | A+ | A | A | A− |
| color | A | A+ | A | A | A | A | A | A− |
| render-sdf | A | A+ | A− | A− | A | A+ | A | A |
| render-scan-shading | A | A | A− | A | A | A+ | A+ | A |
| render-plot | A | A | A− | A | C | A | A | A+ |
| render-filter | A | A+ | A | A− | A | A+ | A | A |
| animation | A+ | A | A− | A | A+ | A | A+ | A |
| effects-1 | A | A | A− | A+ | A− | A | A | A+ |
| effects-2 | A | A− | A | A− | A− | A | A | A |
| effects-3 | A | A | A− | A | A− | A | A− | A |
| hardware | A | A+ | A | A | A | A+ | A | A+ |
| wasm-target | A | A | A | A+ | A | A | A | A+ |
| build-ci-tools | A− | A | A− | A | A | A | A | A+ |
| pcb-gen | A− | A− | B+ | A | C+ | A | B | A |
| hs-tests (harness) | A/A+ | A/A+ | A | A/A+ | A/A+ | A−/B+ | A | A+ |
| dd-app-core | A− | A | A− | A− | B+ | A | A | A |
| dd-app-support | A | A | A | A− | A | A | A | A |
| dd-segments | A | A | A− | A− | A+ | A | A | A |
| dd-tools | A− | A | A | A− | A− | A | A | A+ |
| dd-tests | A | A+ | A | A+ | A | A | A | A+ |

> The two low outliers in the matrix are worth reading in context: render-plot's **C in testing** reflected that `plot.h`'s subtle dual-metric arc math was validated only indirectly through scan/effect tests; finding 35 adds the dedicated `plot.h` numerical suite that pins it against independent oracles. pcb-gen's **C+ / B** reflect ungated diagnostics and a stale-report bug in the *candidate-analysis tooling*, not the emitted schematic/board, which round-trips clean through KiCad's own netlist exporter.

---

## Prioritized Fix List

All 34 items are low-severity; the tiers below rank by impact, not by any critical/high classification. Numbering is sequential across all tiers. (Finding 35, added after the review, is a test-coverage follow-up rather than a defect.)

### Priority 1 — Correctness & Functional Risk

1. ✅ **daydream/tools/palettes.html:767** — Keyboard zoom (Shift+ArrowRight) on the Generative tab silently corrupts the hidden Procedural `C_R..D_B` params. The mouse drag-zoom path guards with `if (activeTab !== 'procedural') return;` but `zoomAroundPhase()` does not. *Fix:* add the same guard as the first line of `zoomAroundPhase()`.

2. ✅ **daydream/recorder.js:387** — The `showSaveFilePicker().catch(AbortError)` branch calls `this.stop()` unguarded, so a cancelled save dialog from session A can tear down a superseding session B. Diverges from the module's `if (this.mediaRecorder === recorder)` per-session discipline. *Fix:* guard the abort-path stop with the same session check.

3. ✅ **Holosphere/hardware/pov_single.h:78** — Single-board `show()` constructs the effect with unguarded `new E()`/`delete`, diverging from Phantasm's `new (std::nothrow)` + `HS_CHECK(e != nullptr, ...)` fail-fast OOM guard; an OOM here aborts with no diagnostic. *Fix:* mirror the segmented path's nothrow + named trap.

4. ✅ **Holosphere/hardware/phantasm/gen/analyze_candidates.py:117** — `run_drc` writes every candidate's report to one fixed path with `check=False`, so a `kicad-cli` early-exit makes a stale neighbor's DRC verdict be attributed to the current candidate — a broken candidate can be scored "clean" and recommended. *Fix:* per-candidate temp file (or `os.remove` first) and verify the return code / that the report was rewritten before parsing.

5. ✅ **Holosphere/hardware/phantasm/gen/shorts.py:106** — The script detects real shorts via union-find and prints them but always exits 0, so it cannot act as an automated gate despite being listed as a validation step. *Fix:* count surviving conflict groups and `sys.exit(1)` when nonzero.

6. ✅ **Holosphere/core/engine/platform.h:557** — The host `map()` mock's final signed division is UB for `product == INT32_MIN && divisor == -1` (reachable with an inverted range differing by 1), despite the comment claiming it avoids signed-overflow UB; the device's ARM SDIV returns `INT_MIN` without trapping. Latent (only legacy effects call unqualified `map()`). *Fix:* guard `divisor == -1` (or divide in `int64_t`) to match the device's non-trapping result.

7. ✅ **Holosphere/tests/test_effects.h:1250** — `test_distorted_ring_palette_mod_selection` calls `hs::set_mock_time()` but never `hs::clear_mock_time()`, leaking a frozen clock into the next test (and the smoke/determinism roster) until the first `render_capture` re-pins it. Harmless today only by luck of the downstream effect being frame-clocked. *Fix:* add `hs::clear_mock_time();` at the end, matching sibling wrapped-phase tests.

8. ✅ **daydream/segment_controller.js:976** — On the overrun re-blit tick, `updateStats()` runs after the next generation's `renderParallel()` has zeroed the per-segment arrays, so a multi-tick effect shows `Range='?'` / `0.0 ms` for not-yet-landed segments while a real previous-generation image is on screen. Debug-overlay only; compositing is unaffected. *Fix:* skip `updateStats()` on the re-blit branch, or source its columns from a publish-time snapshot.

### Priority 2 — Test Coverage & Robustness Gaps

9. ✅ **Holosphere/effects/Flyby.h:214** (and **effects/Liquid2D.h:280**) — Unlike MeshFeedback/MindSplatter, the Flyby and Liquid2D preset tables have no compile-time `preset_in_ranges` + `static_assert` proving each preset value lies inside its registered slider range; a future range tightening could silently push a param out of range with no build failure. *Fix:* add the sibling effects' `preset_in_ranges` helper and `static_assert`.

10. ✅ **Holosphere/effects/MobiusGrid.h:259** — The conformal-radius pole branch and the counter-rotation singularity guard are hand-derived numeric branches that still render a plausible frame if wrong, and have only generic smoke coverage. *Fix:* add a `MobiusGridWhiteBox` seam and a test sweeping `z→±1` (finite coord in [0,1]) and `mid≈0` (q stays identity), mirroring `test_liquid2d_glitch_lens_unit_norm`.

11. ✅ **Holosphere/core/color/composition.h:999** — `BakedPalette::clone_from` (raw LUT memcpy, reached via arena compaction) and `step_wipe_rebake` (ColorWipe arming-skip + per-frame rebake) have no direct unit test; a wrong size or missed allocation would only surface at runtime after compaction. *Fix:* add a clone-and-compare test and a small `step_wipe_rebake` skip-once/decrement test.

12. **Holosphere/core/mesh/solids.h:1250** — Nothing enforces that registry names are unique across the three registries; the WASM picker enumerates by index but builds by first-name-match, so a future duplicate name would make the two paths silently diverge. *Fix:* add a test (or constexpr check) asserting name uniqueness and `get_entry(i).name` round-trips through `find_entry`.

13. **daydream/scripts/require-tests.mjs:7** — The zero-test guard hardcodes `'tests'` and `'.test.js'` independently of the `package.json` test glob; changing the glob while stale files remain lets the guard pass while `node --test` runs nothing. *Fix:* derive the pattern from the `package.json` `test` script instead of restating literals.

14. **daydream/.github/workflows/js-tests.yml:8** — The PR test workflow has no `concurrency` group, so rapid pushes spawn overlapping redundant runs and a stale run can report after a newer one. *Fix:* add `concurrency: { group: js-tests-${{ github.ref }}, cancel-in-progress: true }`.

15. **Holosphere/.github/workflows/ci.yml:284** — The debug WASM runtime-smoke (the only `-sASSERTIONS=1` build) omits `WASM_SMOKE_FRAMES` and runs 3 frames, never reaching the frame-48 ShapeShifter cut / arena compaction that the release smoke bumps to 120 frames to hit — so no CI config exercises late-lifecycle code *with assertions on*. *Fix:* set `WASM_SMOKE_FRAMES: 120` (or ≥60) on that step.

16. **Holosphere/tests/test_color.h:597** — `test_oklch_to_pixel_bounded` asserts each `uint16_t` channel is `<= 65535`, a tautology that passes for any implementation including a broken overflowing cast. *Fix:* assert a non-tautological property (a deeply out-of-gamut color saturates a channel to exactly 65535, or round-trips in-gamut), or delete it as redundant with the two clamp tests that already have teeth.

### Priority 3 — Maintainability: Dead Code & Duplication

17. ✅ **Holosphere/core/math/3dmath.h:44** — `math::EPS_UNIT` (scalar unit-length tolerance) is declared and documented but referenced nowhere; live checks use the squared-length `EPS_UNIT_VEC_SQ`/`EPS_UNIT_QUAT_SQ`. Reads as a live knob and invites the wrong (non-squared) tolerance. *Fix:* delete it and its doc-block line.

18. ✅ **Holosphere/core/engine/transformers.h:170** — In `spawn_impl`'s non-pinned reclaim callback, the `repeats` snapshot and `if (!repeats)` guard are dead: an always-on `HS_CHECK(!p->repeats())` a few lines above proves `repeats` is always false here. *Fix:* drop the snapshot and unconditionally deactivate; keep the pinned branch's live re-query.

19. ✅ **Holosphere/core/mesh/mesh.h:196** — `HEVertex::half_edge` (and the whole `HalfEdgeMesh::vertices` array) is written on every build but never read; the `relax` operator explicitly scans `half_edges` directly instead. Wastes budget-sensitive temp arena on each Conway build. *Fix:* remove the struct, member, and its bind/populate loop.

20. ✅ **Holosphere/core/render/sdf.h:1764** — `Face::max_r2` is a persistent member written and read only inside `setup_frame_and_polygon`; it enlarges every `Face` and implies post-construction validity it doesn't have. *Fix:* make it a local; `radius`/`max_dist_sq` are the real persisted derivatives.

21. ✅ **Holosphere/core/color/color.h:397** — `blend_max`/`blend_over`/`blend_under`/`blend_mean` and the free `blend_add` have no production callers (only `blend_alpha` is on the render path); vestigial from a removed selectable-blend-mode feature. *Fix:* delete them (and their orphaned tests), or document them as an intentionally-retained public library.

22. ✅ **Holosphere/core/mesh/solids.h:28** — `MAX_VERTS`/`MAX_INDICES` are documented capacity budgets whose only references are their own self-referential `static_assert`s; they size no arena and validate no mesh (runtime enforcement is `narrow_index` against `INT16_MAX`). *Fix:* wire them into the counts they claim to bound, or remove them and keep the index-width relationship documented at `narrow_index`.

23. ✅ **Holosphere/hardware/phantasm/gen/sexp.py:119** — `_find` is defined but no module calls it (all callers use `_val`); it can drift out of sync. *Fix:* delete it, or add a caller/test if it is intended public API.

24. ✅ **Holosphere/targets/wasm/wasm.cpp:1373** — The `spline_catmull_rom_tangents` binding re-implements the finite-check + `EPS_NORMALIZE_SQ` degeneracy loop that `eval_cubic_spline` centralizes for the other spline exports, so a future guard change must be made twice or the catmull export diverges. *Fix:* extract a shared `spline_inputs_valid` predicate and call it from both.

25. **Holosphere/hardware/phantasm/gen/pcb.py:34** — The `kicad-cli` netlist-export sequence and the `F(n,k)`/`fmt`/`uid` helpers are duplicated across `pcb.py`, `check.py`, `shorts.py`, and `builder.py`; a KiCad-flag or schema change drifts the copies. *Fix:* hoist `export_netlist` and the pure helpers into a shared module and import them. ✅

26. ✅ **Holosphere/core/render/plot.h:1445** — `Ring::sample` and `DistortedRing::sample` contain byte-identical angle-addition-from-`TrigLUT` blocks plus the tangent-plane vector build; a future LUT-recovery correction must be applied in both or they diverge. *Fix:* extract a `static inline` helper (kept `inline` so codegen is unchanged).

27. ❌ **Holosphere/core/render/plot.h:521** — Each planar edge endpoint is run through `azimuthal_project` (acos + fast_atan2 + fast_cosf/sinf) up to 3× per frame across the arc pre-pass, clip cull, and draw. Per-coarse-segment (not per-pixel), so cost is small. *Fix:* cache the per-vertex projection once in the pre-pass and have all three consume it. — REJECTED: render-path projection caching not provably bit-identical / disproportionate for a per-segment win. The three sites do NOT receive identical inputs: the clip-cull site (`edge_row_span`) is routed through the pipeline's `could_intersect_clip`, where a `World::Orient` stage re-emits the edge under each rotation it applies at plot() time (transformed positions and a transformed basis, possibly multiple invocations per edge), so a pre-pass cache would diverge from the rendered latitude whenever an orientation filter meets a clip band. The pre-pass also projects only non-seam edges. A shared cache is only provable for the pre-pass↔draw pair (a partial win) and would require a new per-vertex arena buffer plus a signature change to `rasterize_planar_strategy` — invasive for a small per-segment cost on a behavior-pinned hot path.

28. ✅ **Holosphere/effects/ShapeShifter.h:139** — `draw_all()` calls `Palettes::RICH_SUNSET.get(t)` (a virtual call with 3× `fast_cosf` + 3 sRGB LUT interps) once per visible layer every frame (up to ~144), where a once-baked 256-entry `BakedPalette` — the pattern its siblings already use — is behaviorally identical. Per-layer, not per-fragment, so impact is small. *Fix:* bake `RICH_SUNSET` in `init()` and look it up.

### Priority 4 — Documentation & Cosmetic

29. ✅ **Holosphere/core/render/filter.h:1080** — The `Screen::Trails` saturation-eviction comment claims `World::Trails` "swap-removes likewise" at saturation, but `World::Trails` evicts the oldest via FIFO `pop_front()`; it only swap-removes dead slots during `flush()` (an unrelated aging path). No runtime effect. *Fix:* reword to state World=FIFO-oldest, Screen=drops slot 0, both acceptable because per-point ttl fade absorbs the transient.

30. ✅ **Holosphere/README.md:1002** — README §7.6 lists palettes as camelCase code identifiers (`darkRainbow`, …) in backticks, but the constants are `Palettes::DARK_RAINBOW` (SCREAMING_SNAKE); a copied identifier fails to compile. The camelCase spellings also leak into `DreamBalls.h`/`GnomonicStars.h` prose. *Fix:* rewrite to the real identifiers, or drop the backticks and mark them informal display names.

31. ✅ **Holosphere/core/mesh/solids.h:690** — `truncatedIcosahedron_hk58_chamfer63`'s `@brief` pins exact counts `(V=990, F=452, I=2880)` — the only recipe that does — which contradicts the project's own stance (`test_solids.h` deliberately does *not* golden per-entry counts because the generators are actively tuned). A retune silently makes the doc lie. *Fix:* drop the parenthetical to match the other 23 briefs.

32. ✅ **Holosphere/effects/SphericalHarmonics.h:344** — The `MAX_DEGREE` comment claims `factorial()` "loses precision at 14!"; float loses exact-integer representation above 2^24 (by 13!), so the stated headroom is wrong (functionally harmless at the shipped `MAX_DEGREE=4`, whose max argument is 8). *Fix:* correct the comment to describe the real bound (≈7 sig-digit relative precision; exact-integer loss above 2^24; overflow near 35!).

33. ✅ **Holosphere/core/animation/motion.h:17** — `Path<W, RESOLUTION>`'s first parameter `W` is never used; every instantiation writes `Path<32>`, so `32` binds to the dead `W` and `RESOLUTION` silently defaults to 1024 (a 12 KiB ring, not the apparent 32-point path). The `@tparam W` doc is inaccurate; the class is also unused in production (effects use `ProceduralPath`). *Fix:* drop `W`, make `RESOLUTION` the sole first parameter, update the four test instantiations — or reconsider whether the test-only `Path` earns its place.

34. ✅ **daydream/daydream.js:74** — `resolutionPresets` keys read H×W (`"Phantasm (144x288)"`) while the rest of the project uses W×H (README "288×144", `setResolution(p.w, p.h)`); the label is also the user-visible `?resolution=` deep-link token. Internal math is consistent — UI/doc only. *Fix:* relabel to W×H, optionally with a one-time alias map for old deep-links (validators already reject unknown values gracefully).

### Priority 5 — Test Coverage (post-review follow-up)

35. ✅ **Holosphere/core/render/plot.h:150** — The dual-metric planar-arc primitives (`azimuthal_project`/`azimuthal_unproject`, `planar_arc_cumul`/`planar_arc_length`) were validated only indirectly: existing tests reuse them as their own ground-truth helpers or pin the rasterizer against `planar_arc_length` itself, so a systematic error in the projection or the anisotropic arc-length integral would move oracle and code together and stay green — the gap behind render-plot's C testing grade. *Fix:* add a dedicated numerical suite to `test_plot_scan.h` pinning the projection round-trip, the projection-radius = geodesic-angle identity, unprojection against an independent libm great-circle reconstruction, `planar_arc_length` against a fine libm quadrature, and the radial-isometric vs azimuthal-bow dual-metric discriminator.

---

## Notes on Scope and Deliberate Designs

The following were examined and **not** filed, because they are validated-intentional decisions (several confirmed load-bearing in prior cycles): the 2-buffer + `buffer_free()` double-buffer; always-on branchless `wrap()`/`wrap_t()` boundary guards; slerp's antipodal raw-dot gate; the near-pole azimuth interval gap in `sdf.h`; `classify_faces_impl`'s `always_inline`; `assert`-not-`HS_CHECK` on the `ArenaVector`/`MeshPaletteBank` per-pixel `operator[]`; Conway 2-gon self-pairing and degree-≥2 hankin; ChaoticStrings' param-member placement; the Feedback warp tap anchoring on `d00`; the un-mockable eDMA register layer; the `needs_full_frame` forgotten-override footgun (known, no cheap test); daydream's deliberate absence of a typecheck and its `@ts-nocheck` DI test files. A handful of near-singularity behaviors (stereographic pole-cap magnitude discontinuity, gnomonic recognize-band, non-march-safe cubic `smin`) were validated as correct consequences of representing infinity with a finite sentinel and left unfiled.

Two out-of-scope observations worth a project-level note: (a) `World::Mobius` does not itself guard the stereographic pole singularity — a sample landing exactly at the pole could produce a non-finite vector reaching the finite-only sink; the correct fix belongs in `geometry.h`'s `stereo()`/`inv_stereo()`, and reachability depends on whether `MobiusGrid` can place a sample exactly at the pole. (b) `SHMath`'s closed-form recurrences (`associated_legendre`, `normalization`, `decode_lm`) are exercised only through the render smoke path; a direct numerical unit test would be cheap insurance.
