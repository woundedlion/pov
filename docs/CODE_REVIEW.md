# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware and the daydream TypeScript/JavaScript web simulator.
**Excluded by request:** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/rotate.h`, and third-party vendored code (`core/FastNoiseLite.h`, `.pio/`, `node_modules/`).
**Method:** the codebase was partitioned across ~19 independent examination agents (one per subsystem), each grounded in the README architecture. Every material finding was then handed to a separate validation agent that read the cited code and returned CONFIRMED / REFUTED / PARTIAL with evidence. Two findings were refuted on validation and dropped; one was narrowed. Only validated defects appear below.

---

## Overall Grade: **A−**

This is a remarkably disciplined codebase. It is a hard-real-time embedded rendering engine that *also* compiles to WebAssembly and runs bit-identically in a browser, with a fail-fast safety doctrine that is consistently applied, compile-time-enforced where possible, and verified by an out-of-process death-test harness. The defects found are overwhelmingly low-severity polish: documentation drift between a very large README and the code, a single wrong-direction boolean in CSG anti-aliasing selection, and a few footprint/robustness nits. No critical or high-severity correctness, memory-safety, or concurrency defects were found anywhere in scope — including the areas most likely to harbor them (the ISR double-buffer, the multi-board sync protocol, the arena allocator, and the WASM zero-copy boundary).

### Quality-Dimension Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A− | Numerical edge cases (poles, antipodes, 0/0, wraparound) are handled with explicit, documented policy. The only genuine code-logic bug is the CSG `is_solid` direction (cosmetic AA impact). |
| **Memory & Concurrency Safety** | A | Arena lifetime model is auditable and trap-guarded; ISR/DMA double-buffer reasoning is airtight; type-erasure is heap-free and alignment-correct; WASM detach contract is correct. |
| **Architecture & Elegance** | A | Compile-time resolution parameterization, the variadic filter pipeline with static-assert-enforced domain ordering, explicit `(target, temp)` arena passing, and X-macro single-sources-of-truth are model designs. |
| **Interface Expressiveness** | A− | The shader/Fragment model, SDF/Scan/Plot families, and Feedback::Style presets read cleanly. Minor: a few generic templates lack `requires`-clause guards, so misuse fails deep in instantiation. |
| **Readability** | A− | Comments are load-bearing and explain *why*. The one tax is over-documentation in a few hot spots (multi-paragraph proofs where a line would do). |
| **Performance** | A− | Zero-overhead abstractions, LUTs, pre-baked palettes, and pole-aware sampling. The reaction-diffusion effects are the heaviest per-frame budget and aren't runtime-gated at H=144. |
| **Documentation** | B+ | The README is exceptional in ambition and mostly accurate, but it has drifted from the code in several places (Conway templating, plot stepping, Archimedean completeness, six effects' parameter lists). |
| **Test Coverage & Rigor** | A− | Deep oracle-driven white-box tests, a rigorous death harness that checks the *specific* illegal-instruction status, and dedicated TUs that recompile under shipping math flags. Thin spots: mesh-raster and generators modules. |
| **Build / CI** | A− | Pervasive toolchain pinning, self-gating coverage jobs, both Teensy images size/layout-gated, provenance hashing. Minor CI-minute and gate-precision nits. |

### Subsystem Grades

| Subsystem | Correctness | Safety | Elegance | Docs | Notes |
|---|---|---|---|---|---|
| Math / geometry (`3dmath`, `geometry`, `util`, `waves`, `easing`) | A− | A | A | A− | Exemplary singularity discipline; split-trig LUT. |
| Color (`color`, `color_luts`, `palettes`, `styles`) | A | A | A | A | OKLCH shortest-arc verified; rigorous rounding. |
| Memory / containers (`memory`, `concepts`, `inplace_function`) | A | A− | A | A+ | Dual debug-stamp staleness detection is standout. |
| SDF + Scan (`sdf`, `scan`) | A− | A | A | A | Compile-time span-cap proof; one `is_solid` bug. |
| Plot / Filter / Canvas / Transformers | A− | A− | A | B+ | Static-assert-enforced pipeline ordering. |
| Animation / Generators / Presets | A− | A− | A | A | Compile-time-enforced borrow contracts. |
| Mesh + Conway | A− | A | A | B | Excellent arena hygiene; README API claim wrong. |
| Hankin / Spatial / Solids / Reaction graph | A− | A | A− | B+ | Reproducible 90 KB table provenance. |
| Platform / Engine / Registry | A− | A | A | A | Bit-exact FastLED host mocks; X-macro SSOTs. |
| HW drivers (`dma_led`, `hd107s_frame`, `pov_single`) | A | A− | A | A+ | Reference-grade concurrency reasoning. |
| Phantasm sync (`pov_sync`, `pov_segmented`) | A− | A− | A | A | Wrap-proof flywheel; race-free mailbox. |
| WASM target + marshaling | A | A | A | A | Correct zero-copy contract; pre-sized buffers. |
| Effects (27) | A− | A− / A | A / B+ | A− / B+ | Mature; arena budgets static-asserted. |
| daydream core JS | A− | A− | A | A | Disciplined detach handling + dispose symmetry. |
| daydream workers + tools | A | A− | A | A | No protocol drift; shared-WASM math reuse. |
| Test suite | A | A | A− | A | Rigorous death harness; perturb-between-runs. |

---

## Prioritized Fix List

Items are numbered sequentially across all priority sections. Each cites the location, the defect, and the suggested fix.

### Critical

*None found.*

### High

*None found.*

### Medium

1. ✅ **`core/sdf.h:805` & `core/sdf.h:1031` — CSG `is_solid` boolean is the wrong direction.** `Union<A,B>::is_solid = A::is_solid && B::is_solid`, but a union is `min(dA,dB)` and keeps *both* interiors, so a solid∪stroke is still solid — it should be `||`. `Subtract<A,B>::is_solid` is likewise `A && B`, but carving a stroke out of a solid leaves it solid — it should track `A` only. Consequence is not lost coverage (interiors still fill) but a wrong anti-aliasing branch in `core/scan.h:61` (`solid ? pixel_width : 0.0f`): mixed-solidity CSG composites get the stroke opacity-falloff ramp at their boundary instead of the 1-pixel quintic solid-edge AA. *Fix:* `Union` → `||`; `Subtract` → `A::is_solid`.
2. ✅ **`README.md:1074` (§7.7) — false Conway-operator API claim.** The README states "All Conway operators are templated on input mesh type and take `(const MeshT& mesh, …)`; `MeshT` can be either `PolyMesh` or `MeshState`." Every operator in `core/conway.h` (`dual`, `kis`, `ambo`, `truncate`, `expand`, `chamfer`, `relax`, `snub`, and the composed `gyro`/`meta`/`needle`/`zip`/`bevel`) takes a concrete `const PolyMesh&` — none is templated on mesh type (the code's own contract at `conway.h:311` confirms it). A reader trusting the README would write `dual(meshState, …)` and hit a compile error. *Fix:* describe the actual `const PolyMesh&` signature.
3. ✅ **`README.md:421` (§6) — stale Plot-path sampling description.** The README says the plot rasterizer uses "adaptive step size scaled by `sin(φ)` for uniform sampling." `core/plot.h:441-468` (`screen_step`) explicitly *replaced* that heuristic with full 2-D screen-velocity tracking, and the in-code comment calls the `sin(φ)`-only model deficient ("ignored the curve's latitudinal motion"). *Fix:* update the README to describe 2-D screen-velocity stepping.

### Low

4. ✅ **`README.md` §9 — six effect entries omit their `**Parameters**` line.** Comets, ChaoticStrings, GnomonicStars, MindSplatter, PetalFlow, and Moire all register live parameters (verified in their headers) but their Effects Reference entries lack the parameter list every sibling entry carries, so users can't discover the sliders. (DistortedRing and ShapeShifter, initially flagged, *are* already documented — excluded.) *Fix:* add the missing parameter lists.
5. ✅ **`README.md:1546` (§9) — HankinSolids description inaccurate.** It claims the effect "sequences through the full Archimedean solid library," but `effects/HankinSolids.h:84` uses `get_simple_solids()` = Platonic (5) + Archimedean (11). It is neither Archimedean-only nor the full Archimedean set. *Fix:* "sequences through the Platonic and Archimedean solids."
6. **`core/solids.h:1056-1072` — Archimedean registry returns 11 of 13 solids.** `truncatedIcosidodecahedron` (`solids.h:549`) and `snubDodecahedron` (`solids.h:561`) are defined and used as Catalan/Islamic bases but are absent from `simple_registry`; `get_archimedean_solids()` returns count 11, while the registry comment ("Platonic and Archimedean solids") implies the full set. *Fix:* add the two solids to the registry or note the deliberate omission.
7. **`effects/MindSplatter.h:270` — particle-index guard uses debug-only `assert`.** The per-fragment bounds check is stripped under `NDEBUG` (defined in the device build), so an out-of-range `p_idx` would read `pool[p_idx]` out of bounds on hardware, contrary to the project's always-on fail-fast doctrine. Because this is a per-fragment hot path, the project-consistent fix is *not* a per-fragment `HS_CHECK` (the doctrine forbids that) but validating the index once on the cold spawn/emit path where `p_idx` originates.
8. **`core/filter.h:1028,1067` — `Screen::Blur` is not flagged as crossing segments.** `Blur` derives from `Is2D` (`crosses_segments = has_history = false`) yet taps rows `cy±1` with only a canvas-edge clamp. In segmented rendering it samples unwritten neighbor-segment rows (reads black at the seam), producing a 1-pixel seam artifact. *Fix:* give spatially-neighborly stateless filters a `crosses_segments = true` trait, or clamp taps to the segment band.
9. **`effects/SplineFlow.h:21` — `MAX_TRAILS = 30000` is ~3.5× over-provisioned.** At 8 bytes/item this is ~234 KB (~73% of the 330 KB device arena) with no `configure_arenas` override, while the live trail set is ~144 samples/frame × 60-frame lifetime ≈ 8.6 k points. It fits (an overflow would trap) but leaves almost no headroom. *Fix:* right-size to the lifetime×rate budget or repartition.
10. **`core/transformers.h:286-288` — fragile sentinel recognition in `gnomonic_mobius_transform`.** Output stays finite (the equator singularity *is* guarded in `gnomonic`), so this is not a blow-up — but `inv_gnomonic` recognizes the pole only at *exactly* `±STEREO_INF` (no margin, `3dmath.h:788`), while a Möbius map runs between the projections and can scale the sentinel off its exact value. The stereographic path already solved this with `STEREO_INF_RECOGNIZE = STEREO_INF * 0.5f`; the gnomonic path should mirror it.
11. **`tools/mobius.html:638,642` (daydream) — touch handlers read `e.touches[0].clientX` with no empty-list guard.** A multi-touch/`touchend` race could dereference `undefined`. *Fix:* guard `touches.length` before indexing.
12. **`daydream.js` — documented `getParamGeneration()` re-fetch contract is never consumed.** README §10.2/§10.6 specify caching `getParamGeneration()` and re-fetching definitions when it changes; the symbol has zero references in JS. It is safe today only because the GUI is fully rebuilt on every effect/resolution switch. *Fix:* consume the generation token, or document the rebuild-everything reliance so a future partial-update path doesn't silently mis-map the value stream.
13. **`segment_controller.js:492` (daydream) — compositing seam is unguarded.** `composite()` validates the rect against display bounds but never asserts `r.quadW === x1-x0 && r.quadH === y1-y0`; a protocol drift would silently `subarray` a short/long compact buffer instead of trapping, unlike the fail-fast rigor elsewhere. *Fix:* assert the quad dimensions match the rect.
14. **`effects/RingSpin.h:97,158` — `spawn_ring`'s `normal` parameter is effectively dead.** `init()` passes `Y_AXIS` for all four rings, so the per-ring great-circle-plane distinction is illusory; variety comes only from random-walk seeds. *Fix:* pass distinct axes if planes were intended, or drop the parameter.
15. **`effects/ShapeShifter.h:48,127` — "Count" slider range exceeds the render cap.** The slider goes to 128 but `drawAll` caps rings at `kMaxRings = H` (20 at 96×20), so the slider's upper ~84% is inert at the hardware resolution. *Fix:* track `H` in the registered max, or document the cap.
16. **`core/led.h:80-119` — correction-guard nesting contract is mis-described.** `NoColorCorrection` and `NoTempCorrection` share one depth counter, so *any* two live guards (even of different types) trap, but the contract comment says only "do not nest." *Fix:* clarify that at most one correction guard may be live at a time.

### Informational / Polish

17. **`effects/Flyby.h:37` — stale `// TODO:` preset comment left in shipped `init()`.** Remove dead commentary.
18. ✅ **`core/effect_registry.h:71` — `kReserveHint` comment says "~30 effects"; the roster is 27.** Harmless (reserve is a hint) but stale against the X-macro SSOT.
19. ✅ **`.github/workflows/docs.yml:8` — docs deploy has no `paths:` filter.** Every push to master rebuilds and redeploys the full Doxygen+Graphviz site. *Fix:* scope to `**.h`, `**.cpp`, `README.md`, `Doxyfile`.
20. **`.github/workflows/ci.yml:329-381` — provenance gates diff *all* numeric tokens.** The LUT/reaction-graph drift check tokenizes the whole file (includes, type widths, array bounds), so it's self-consistent but coarse. *Fix:* compare only the array-body tokens.
21. ✅ **`.githooks/pre-commit:90` — local Teensy gate builds only the `phantasm` image.** A Holosphere-only budget/layout regression passes the hook and fails only in CI. *Fix:* build both images (`pio run -e holosphere -e phantasm`), matching CI's teensy-size job.
22. **Cross-effect duplication — the `strobe_columns()` override and its identical doxygen block are copy-pasted across ~9 effects** (and `needs_full_frame()` across two). A shared base default would remove ~50 lines of boilerplate.
23. ✅ **`effects/Voronoi.h:208-220` — coarse-grid coherence reuse is an exact-shading approximation.** When all four block corners agree on the nearest pair, interior pixels skip the k=2 query and can mis-shade where a third site is genuinely nearest (near poles). Documented in-code as a deliberate speed/accuracy trade; the README's "smoothly blended" wording slightly overstates exactness.
24. **`effects/GSReactionDiffusion.h:115` — 16 substeps/frame is the heaviest per-frame budget in the engine** (16 × 7680 nodes × 6-neighbor Laplacian × 2 buffers) and the H=144 path is not runtime-gated in CI. Confirm it meets the 16 FPS device budget; the same untested-at-H=144 caveat applies to HankinSolids' arena high-water.

---

## Notable Strengths

- **Fail-fast as an enforced doctrine, not a slogan.** `HS_CHECK` is always-on (survives `NDEBUG`), stdio-free, placed on cold seams only, and the per-pixel hot paths use stripped `assert` backed by a cold trap at the bind site. A death harness (`tests/test_death.h`) verifies the *specific* `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` shape across 34 trap seams, so the safety net is tested rather than assumed.
- **Correctness pushed into the type system and the build.** Compile-time span-cap proofs in the rasterizer, `static_assert`-anchored filter-pipeline ordering, deleted rvalue overloads that make dangling borrows un-compilable, X-macro single-sources-of-truth for resolutions/effects/mesh-ops, and a CI that proves every test lands in exactly one shard.
- **A genuinely hard concurrency problem solved cleanly.** The single-core ISR double-buffer (2 buffers + `buffer_free()` gate) and the Phantasm multi-board flywheel (wrap-proof 64-bit position math, count-coded sync alphabet, race-free `try_claim` mailbox) are reasoned to the metal, with relaxed atomics correctly justified and cache coherence separated from atomic ordering.
- **Device/simulator parity by construction.** FastLED host mocks are bit-exact reimplementations (including the Cortex-M7 UDIV-by-zero=0 behavior), the WASM zero-copy readback contract is correct and pre-sized so the backing buffer never moves, and the daydream tools route geometry/mesh math through the shared WASM engine rather than re-implementing it.
- **Reproducible generated artifacts.** The 90 KB reaction-diffusion K-NN table and the sRGB↔linear LUTs are regenerated bit-for-bit by committed Python generators and CI-gated against the checked-in copies.
