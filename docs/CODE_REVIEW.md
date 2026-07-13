# Holosphere + daydream Comprehensive Code Review

**Review date:** 2026-07-12  
**Overall engineering grade:** **B+**  
**Holosphere engine/effects grade:** **A-**  
**Phantasm hardware/firmware readiness:** **B-**  
**daydream runtime grade:** **A-**  
**daydream tools/presentation grade:** **B**

## Scope and method

This review covers the Holosphere C++ engine, effects, tests, firmware, hardware design, PCB generators, build system, CI, documentation and visual assets, plus the daydream JavaScript/WASM runtime, workers, tools, tests, deployment, documentation and visual assets.

The following were explicitly excluded and were not opened:

- `core/engine/effects_legacy.h`
- `targets/Holosphere/Holosphere.ino`
- `core/math/rotate.h`
- The previous contents of `docs/CODE_REVIEW.md`
- The previous contents of `docs/PROJECT_ASSESSMENT.md`

Untracked daydream `.claude/**` and `prompts/**` content was treated as user/editor WIP rather than product source and was not read. Existing modified daydream WASM artifacts were inspected without being changed.

Reviewers consulted the repository README first, then exhaustively inspected disjoint components. Candidate findings were retained only after they were checked against the implementation, tests, generated artifacts or hardware files. Non-issues and already-defended cases were rejected rather than recorded as defects.

### Validation evidence

- Holosphere native suite: **49/49 CTest tests passed** after a fresh build; the effect lifecycle also passed with `HS_SMOKE_FRAMES=120`.
- daydream: **335/335 Node tests passed**, with no failures, skips or TODOs.
- Teensy gate self-tests: **42/42 passed**.
- KiCad: **ERC clean; DRC clean with zero unconnected pads**; netlist and short checks passed.
- Holosphere roster/screenshots: **26/26** effects represented; no corrupt or blank source assets.
- Current daydream WASM hash matches `holosphere_wasm.wasm.sha256`; its recorded source commit is an ancestor of Holosphere HEAD.
- Current Phantasm ELF was independently inspected with Teensy `nm`, `size` and `objdump`.

## Executive appraisal

The project is substantially better engineered than the typical hobbyist LED/POV codebase. Its strongest qualities are the spherical rendering model, arena discipline, linear-light color pipeline, mesh/effect vocabulary, pure-core/device-shell hardware split, synchronization state machine, WASM boundary, artifact provenance and unusually deep tests. The APIs often read as a compact visual language, not merely a collection of drawing routines.

The overall grade is held below A by one release-blocking firmware resource failure, release-only mesh UB at a malformed public seam, incomplete physical/manufacturing readiness, several false/red CI paths, and browser failure paths that sit above the otherwise excellent unit-test layer. Documentation breadth is exceptional, but parts of it have accumulated as historical strata and now contradict the shipping code and board.

## Quality grades

| Quality dimension | Holosphere | daydream | Combined rationale |
|---|---:|---:|---|
| Functional correctness | B+ | A- | Normal paths are heavily validated; Phantasm OOM, malformed mesh input and browser exception paths are material exceptions. |
| Reliability / long-run stability | B | B | Strong deterministic lifecycle design, offset by the Phantasm playlist trap, a 13-day effect freeze and incomplete browser rollback/cleanup. |
| Memory / resource safety | B+ | A- | Arenas, fixed-capacity storage and detached-view healing are excellent; target heap sizing and a few invalid-input seams remain unsafe. |
| Numerical robustness | A- | B+ | Spherical seams, poles and finite math are handled carefully; lifetime NaN, biased axis sampling and two tool-domain discontinuities remain. |
| Concurrency / real-time behavior | A- | A- | ISR/DMA and worker generation fences are standout work; one volatile ISR race and partial worker-spawn leak remain. |
| Performance / efficiency | A | A- | Compile-time pipelines, LUTs, DMA, zero-copy views, pooling and staged workers are professional; uncapped tool DPR is a localized miss. |
| Architectural elegance | A- | A- | Coherent domain boundaries and pure cores; globals/singletons and oversized composition/state-machine modules constrain composition. |
| Interface expressiveness | B+ | B+ | Filters, shapes, palettes, mesh recipes and protocol types are expressive; stringly state, parallel mesh arrays and weak callable constraints leak contracts. |
| Maintainability | B+ | B | Strong naming/invariants/tests, but dense template headers, inline tool controllers, stale IDE metadata and documentation sediment increase cost. |
| Readability / comment discipline | A- | B+ | Most comments record real invariants; a few are overlong or stale, and one review number remains in build source. |
| Test quality | A | A- | Exceptional host, death, differential, lifecycle and WASM parity coverage; browser E2E and HIL are the principal blind spots. |
| Build / CI quality | B | B+ | Broad multi-target gates and good presets, but current Holosphere sharding is red and a fallback size parser can false-green. |
| Reproducibility / provenance | B+ | B+ | Generated tables and WASM hashes are strong; routed PCB regeneration and binary-to-source proof are incomplete. |
| Portability | B | B | Deliberately good across Teensy/Clang/WASM/modern browsers, but not general C++, MSVC or legacy-browser portable. |
| Security / defensive design | B+ | B+ | Fail-fast checks, text-only DOM insertion and exact npm pins are strong; two supply-chain paths and missing governance need work. |
| Accessibility | N/A | C+ | Main keyboard behavior is respectable, but contrast, live regions, motion preferences and tool controls are incomplete. |
| UX / error recovery | B+ | B | Rich controls and diagnostics in normal use; bootstrap and transactional failure states can strand the user. |
| Documentation breadth | A+ | A | The architecture, protocols and build are documented at unusual depth. |
| Documentation accuracy | C+ | B- | Several active specs, diagrams, links and copied claims no longer match current code or hardware. |
| Electrical engineering | B+ | N/A | ERC/DRC-clean design with thoughtful power and signaling; sync-node placement and decoupler placement need correction. |
| Mechanical / manufacturing readiness | C | N/A | No rotor mounting holes, weak assembly markings and non-gating fab output prevent production readiness. |
| Artistic system design | A | A- | A coherent mathematical/generative visual instrument with sophisticated color and geometry; gallery curation has minor drift. |
| Community / reuse readiness | C+ | B- | Documentation is valuable, but restrictive split licensing and missing contributor/security guidance limit adoption. |

## Architectural and interface assessment

### What is especially good

- The world → screen → pixel pipeline makes coordinate-domain changes explicit and statically rejects invalid filter ordering.
- Sixteen-bit linear-light compositing, perceptual palette work and spherical seam/pole handling are technically and artistically sophisticated.
- Arena ownership and scratch polarity are visible at call sites; device hot paths avoid hidden allocation.
- Mesh recipes, Conway/Hankin operations, SDFs, curves, transformers, filters and animation policies form a genuine visual DSL.
- Hardware logic is separated into host-testable arithmetic/state machines with a thin Teensy register/ISR shell.
- Phantasm synchronization handles modular time, missed symbols, reacquisition, beacons and exactly-once boundaries with exceptional rigor.
- daydream generation fences old worker frames, stages complete generations, heals detached WASM views and versions its worker protocol.
- Generated LUTs, reaction graph data, effect rosters, screenshots and WASM artifacts have unusually strong mechanical provenance.

### Architectural costs

- Global arenas, timeline, framebuffers, RNG state and mutable daydream singleton state assume one active engine. That is efficient on device but limits concurrent previews, reentrancy and isolated composition.
- The flat mesh interchange format exposes parallel counts/indices/views whose consistency is not fully embodied in a validated type.
- Several objects are default-constructible but unusable until `bind()`, `bake()` or `init_storage()`, spreading temporal validity rules.
- Units such as radians, frames, normalized radii, phases and pixels are plain numeric types; this keeps the embedded ABI small but weakens self-documentation.
- daydream’s 818-line composition root, 1,004-line worker controller and large inline tool controllers centralize too much implicit state.
- The copied cross-repository README and historical audit/spec documents do not have a reliable freshness model.

## Prioritized findings

### P0 — Release blocker

1. ✅ **Phantasm cannot allocate `Dynamo` and will trap during the default playlist.** `targets/Phantasm/Phantasm.ino:68-69`, `tools/teensy_budgets.json:41-57`, current `.pio/build/phantasm/firmware.elf`. The ELF exposes `_heap_start=0x2027ed80` and `_heap_end=0x20280000`, only **4,736 bytes**, while `construct_effect<Dynamo<288,144>>` loads **5,276 bytes** before nothrow `operator new`; allocator metadata and the earlier `g_pov` allocation make the real margin worse. The seventh roster entry is reached at roughly 12 minutes with the current 120-second cadence. Reclaim OCRAM or use fixed maximum-sized effect storage, then add a RAM2 free-floor gate derived from the largest target `sizeof` plus allocator/boot margin.

### P1 — High priority

2. ✅ **Flat mesh length inconsistencies can reach release-only UB.** `core/mesh/mesh.h:214-275,348-416`. `build_half_edge_mesh` and `MeshOps::compile` do not require the sum of face counts to equal the flattened index length; short sums sort uninitialized records or ignore trailing indices, and long sums can index through an `assert` that disappears on device. Add overflow-safe preflight accumulation and exact equality checks at both public seams, with short/long mismatch death tests.

3. **The 480-RPM rotor PCB has no defined mounting holes.** `hardware/phantasm/phantasm.kicad_pcb:7617`, `hardware/phantasm/gen/pcb.py:407-415`. The 58.28 × 32 mm outline has no NPTH/mounting footprint and the generator creates none. Define the rotor fastener pattern and swept-envelope clearances, add suitable holes with copper/courtyard keep-outs, and rebalance/regenerate all four boards.

4. **The sync receiver’s high-impedance RC node is spread across the board.** `hardware/phantasm/phantasm.kicad_pcb:4029,6384,7469`. R1, R2 and `C_SYNC` are widely separated; R2-to-cap is about 29.7 mm and `FRAME_SYNC` copper about 70.6 mm. In a BLDC/LED-noise environment this defeats the intent of the spike filter. Cluster the divider/cap at Teensy pin 3/U1, keep the high-impedance node minimal and give the cap a short ground return.

5. ✅ **Fabrication packaging continues after DRC failure.** `hardware/phantasm/gen/fab.py:171-177`. DRC is informational: the subprocess return, violation count and parse failure do not stop BOM/CPL/Gerber ZIP production. Require successful execution, zero error-severity violations and zero unconnected items; withhold/delete the package on failure.

6. **Assembly-critical board-ID and shield instructions are absent from visible silkscreen.** `hardware/phantasm/phantasm.kicad_pcb:1327,3766,4177,5677,7628-7836`. Jumper references are hidden and there is no ID truth table, large board-role field, “MASTER = both open” or “shield jumper: master only” marking. A duplicate master creates push-pull sync contention. Add permanent unobscured assembly markings.

7. ✅ **`DistortedRing` modulation freezes after about 13.1 days at 16 FPS.** `effects/DistortedRing.h:175,263-273`; `tests/test_effects.h:1276-1315`. The unbounded `float mod_time += 0.05f` reaches 1,048,576 after 18,073,721 frames, where the addition rounds to itself, freezing four palette modes. Use periodic/bounded time or integer-cell-plus-fraction representation without a reset seam, and test near the precision threshold.

8. ✅ **The current CI shard matrix omits `unit_poi`.** `.github/workflows/ci.yml:42-48`. A fresh build registers 49 tests, but the four anchored shard regexes match `unit_poi` zero times, so the workflow’s shard-coverage job is guaranteed to fail. Add `poi` to a shard and keep the mechanical coverage check.

9. ✅ **The published circuit diagram specifies an obsolete over-voltage divider.** `docs/phantasm_circuit.svg:126,634,642,690,702,710`; `docs/phantasm_pcb_spec.md:166-177`. The diagram/PNG says 10 kΩ/18 kΩ and DIP/through-hole, while the current design requires 10 kΩ/15 kΩ SOIC/SMD; the spec itself calculates that 18 kΩ can drive the Teensy input to about 3.38 V. Regenerate both diagram assets from the current KiCad/spec source and mark their revision.

10. ✅ **The tracked Visual Studio/Visual Micro project is structurally broken.** `targets/Holosphere/Holosphere.vcxproj:77-78,85-121` and `.filters`. Thirty-nine of 69 source/header entries resolve to nonexistent pre-reorganization paths, including both core translation units, and the project embeds contributor-specific absolute paths. Regenerate it against the current tree, remove machine state and add a project-item existence check.

### P2 — Medium priority

11. ✅ **`Effect` accepts dimensions larger than its fixed shared buffers.** `core/render/canvas.h:55,270,612,614,716`. Current factories use supported sizes, but the public template seam can drive out-of-bounds copies/fills/indexing. Fail fast on positive dimensions within `MAX_W/MAX_H` before activating the effect and add boundary death tests.

12. ✅ **Particle lifetime normalization accepts NaN and zero.** `core/animation/sprites.h:221-222`, `core/render/plot.h:2287`. NaN survives the nested min/max and is converted to `uint16_t` with UB; zero permits `0/0` during rendering. Require a finite lifetime in `[1,65535]`, or explicitly define/guard disabled zero lifetime, with NaN/zero tests.

13. ✅ **Negative timeline delays have frame-dependent unsigned-wrap semantics.** `core/animation/timeline.h:200`. Adding signed `in_frames` to a `uint32_t` timeline makes `-1` mean “near wrap” at frame zero but “already due” later. Trap negative delays and define/check positive addition overflow.

14. ✅ **Voronoi spin axes are biased and have a valid zero-vector trap state.** `effects/Voronoi.h:314-317`. Normalized cube samples make body diagonals about 5.2× as likely as coordinate axes, and the all-midpoint draw normalizes zero. Use the existing tested Marsaglia `random_vector()` helper and add delegated/distribution coverage.

15. **The external Teensy VIN decoupler misses its placement requirement.** `hardware/phantasm/phantasm.kicad_pcb:3470,7156,7458`. `C_DEC2` is about 6.1 mm from VIN versus the specified `<3 mm`, increasing input-loop inductance. Move it adjacent to VIN/GND with a short return.

16. ✅ **Documentation advertises an unsupported N=8/ID2 mode.** `hardware/phantasm/README.md:7,86,98`, `docs/phantasm_pcb_spec.md:238-260`, `hardware/pov_segmented.h:92-108,366-418`. Firmware rejects `N > 4`, reads only two straps and models top/bottom segments per arm. Remove the functional claim and label ID2 reserved, or implement/generalize/test the eight-segment architecture.

17. ✅ **`build_word_` is an intentional C++ data race hidden behind `volatile`.** `hardware/pov_sync.h:1302,1612-1615,1637-1641`. Aligned Cortex-M7 loads are tear-free, but ISR writes and foreground reads are not race-free in the C++ memory model. Use relaxed `std::atomic<uint32_t>` with an explicit reset path, or a minimal interrupt bracket.

18. **The PCB generator cannot reproduce the committed routed four-layer board.** `hardware/phantasm/gen/pcb.py:348-438`, `hardware/phantasm/gen/stackup.py:12-89`. It emits an unrouted board with random UUIDs and targets the unplaced stackup, losing manual routing, zones, silk and fabrication metadata. Separate logical generation from routed-board maintenance, stabilize identifiers and gate structural equivalence.

19. **Hardware-only DMA/sync behavior has no automated HIL validation.** `tests/test_dma_controller.h:10-14`, `tests/test_hd107s_frame.h:9-11`, `hardware/pov_segmented.h:497-590`. Host tests cannot validate LPSPI/eDMA registers, cache coherency, pulse widths, ISR latency or real-strip behavior. Add a Teensy/logic-analyzer smoke rig covering stream bytes, buffer integrity, sync pulses, overrun and long-run wrap.

20. ✅ **Manufacturing metadata and board documentation have drifted.** `hardware/phantasm/README.md:121-138`, `hardware/phantasm/phantasm.kicad_pcb:10-14,35-103,11583,12222`. Track/via counts, pour layers, finish and copper descriptions disagree with the committed board. Generate measured statistics and fabrication parameters from the manufacturing source of truth.

21. ✅ **BOM/CPL completeness is warning-only.** `hardware/phantasm/gen/fab.py:191-224`. Missing supplier mappings or centroid data can still yield orderable-looking output with blank fields. Fail packaging for any assembled SMD lacking a part number, centroid, numeric placement or expected side.

22. ✅ **Malformed `size -A` output can false-pass the Teensy gate.** `tools/teensy_gate.py:420-430`, `tools/teensy_gate_extra.py:81-104`. Unrecognized output becomes synthesized zero-use/full-free regions and can return PASS. Require recognized allocated sections/expected buckets before synthesis and add empty/malformed fallback tests.

23. ✅ **The fast-forward hook permits deletion of `master`.** `.githooks/reference-transaction:24`. An all-zero new OID is skipped despite the hook’s stated non-fast-forward protection. Reject deletion or require the same explicit one-shot override.

24. ✅ **The Teensy gate documentation is contradictory and malformed.** `docs/teensy_ci_gate_spec.md:3-30,890,944`, `tools/teensy_budgets.json:9`, `docs/RAM1_AUDIT.md:5`, `docs/strobe_columns_audit.md:109-136`. It mixes obsolete report-only/overflow/FastLED claims with current active-gate state, contains an unmatched fence and broken moved-file links, and describes already-implemented work as pending. Date historical audits and refresh active truth, links and fences.

25. ✅ **The docs workflow executes an unverified binary with broad permissions.** `.github/workflows/docs.yml:37-40,65`, `justfile:61,65`. A downloaded Doxygen archive is executed without a digest while Pages/OIDC permissions apply workflow-wide; local docs clone a mutable tag. Verify SHA-256, pin the theme commit and scope elevated permissions only to deployment.

26. **Documentation gates allow known broken output.** `Doxyfile:52`, `.github/workflows/docs.yml:21-24`. Doxygen warnings are nonfatal and there is no Markdown fence/link validation, so the unmatched fence and dead links publish green. Add lint/link checks and make warnings fatal after an explicit baseline.

27. **Whole-app bootstrap failures leave an endless loading spinner.** `daydream/index.html:62-73`, `daydream/daydream.js:7-17,153,585-606`, `daydream/driver.js:146-150`. Static import/CDN or synchronous WebGL/constructor failures happen before the only WASM catch. Use a small local dynamic-import bootstrap with a fatal/reload UI and catch WebGL initialization.

28. ✅ **Worker construction failure leaks a partial unmonitored pool.** `daydream/segment_controller.js:278-419`. A reproduced second-worker `SecurityError` leaves the first worker alive, with no fault latch or watchdog. Guard constructor/initial post operations, terminate the partial pool, latch a visible fault and test worker-N synchronous failure.

29. ✅ **Thrown effect/resolution switches are not transactional.** `daydream/daydream.js:223-239,463-487`. GUI is destroyed before engine mutation, and thrown errors do not roll back state/URL/controls. Preserve prior applied state until success or implement guarded rollback and recovery UI.

30. ✅ **Initial application removes the loader before initialization succeeds.** `daydream/daydream.js:577-584`. `applyResolution(true)` runs after overlay removal, ignores `false` and only logs throws. Remove the overlay only after successful resolution/effect initialization; route rejection to the fatal UI.

31. ✅ **`MediaRecorder.start()` failure leaks capture resources.** `daydream/recorder.js:198-224`. A reproduced thrown `NotSupportedError` escapes while retaining recorder, stream, offscreen and an unstopped track. Catch start failure, close/stop all session resources, reset state and test the path.

32. ✅ **Lissajous closed-curve mode can export an open curve.** `daydream/tools/lissajous.html:126,241-264`. Required rational periods such as 8/7 need `14π` but are silently clamped to the `8π` domain. Grow the domain or restrict rational choices to representable closing periods.

33. ✅ **Flat solids rendering mis-triangulates non-convex Hankin faces.** `daydream/tools/solids.html:922-927,1902-1914`. Non-geodesic paths use an `f[0]` fan even though the adjacent code acknowledges non-convex star faces, spilling triangles outside the face. Use centroid fans for the supported star-shaped faces or robust triangulation in both main and thumbnail paths.

34. ✅ **Palette-strip interaction is mouse-only.** `daydream/tools/palettes.html:643-660,1239-1243`. Mobile users cannot set phase or drag-zoom. Replace mouse handlers with pointer events, capture and cancel handling.

35. ✅ **Core geometry-tool controls are not keyboard/screen-reader operable.** `daydream/tools/slider.js:79-91`, `daydream/tools/solids.html:948-967` and related palette/saved-card controls. Labels are unassociated spans, clickable thumbnails are divs, reorder is mouse-only and status/fatal UI lacks live semantics. Prefer native controls, associated labels, keyboard reorder and live regions.

36. ✅ **Tool renderers allocate unbounded high-DPI backbuffers.** `daydream/tools/shared.js:98-100,141-147`. Raw DPR 3 produces roughly 9× fragment/backbuffer work. Apply an intentional DPR cap consistent with the main simulator.

37. **Deployment provenance does not prove the WASM came from the recorded source SHA.** `daydream/.github/workflows/deploy.yml:11-23,55-99`. The job proves binary/hash consistency and tests the recorded engine commit, but does not rebuild/compare or verify an attestation. Rebuild deterministically and byte-compare, or consume a signed build attestation.

38. **Mutable remote dependencies remain in deployment paths.** `daydream/.github/workflows/*.yml` and each geometry tool’s Tailwind include. Actions use major-version tags and tools fall back to unversioned `cdn.tailwindcss.com`. Pin actions to commits and publish versioned/self-hosted production CSS.

39. **There is no real-browser integration/a11y/visual gate.** `daydream/package.json:6-11`, `daydream/.github/workflows/js-tests.yml`. Node tests cover pure logic and actual WASM well but never boot `index.html`, WebGL, real workers, large inline tool controllers, mobile pointers, accessibility or failure UI. Add a small Playwright/axe matrix for boot, switches, segmented mode, tools and key failure paths.

### P3 — Low priority

40. ✅ **Horizontal clip invariants are only partially enforced.** `core/render/canvas.h:147,165`; `core/render/scan.h:905-906`; `tests/test_canvas.h:551`. Setters do not require `x1 <= clip.w` despite promising enforced invariants; later raster code catches it and the test comment is stale. Enforce at the setter, retain downstream defense and update tests/comments.

41. ✅ **`AlphaFalloffShade` accepts a null mandatory callback.** `core/color/composition.h:670,678`. The raw pointer is invoked unchecked. Add a cold-path `HS_CHECK(fn != nullptr)` and death test.

42. ✅ **Callable-wrapper constraints disagree with thunk behavior.** `core/engine/concepts.h:106-133`, `core/engine/inplace_function.h:41-46,129`. Return-incompatible callables pass declared constructibility constraints and fail inside thunk instantiation, including non-void expressions returned from a void thunk. Align invocable-return constraints and use an explicit void discard path, with compile-time tests.

43. ✅ **A review-history reference remains in build source.** `tests/CMakeLists.txt:151`. Remove the `Finding 369` prefix while preserving the useful factual explanation.

44. ✅ **Main UI small text fails normal-text contrast.** `daydream/styles/index.css:260-268,316-319,334-346,639-645`. Ten-to-eleven-pixel `#555/#666` text on black/dark gray is about 2.5–3.7:1 rather than 4.5:1. Raise muted colors or size/weight.

45. ✅ **Dynamic loading/fault state is not announced.** `daydream/index.html:62-66`, `daydream/segment_controller.js:829-854`. Add polite `role=status` for loading and assertive alert semantics/focusable recovery for failures.

46. ✅ **Infinite UI animations ignore reduced-motion preference.** `daydream/styles/index.css:94-112,182-193`. Provide static recording/loading alternatives under `prefers-reduced-motion`.

47. ✅ **JSDoc/type contracts are not enforced in CI.** `daydream/worker_protocol.js:9-14`, `daydream/package.json:6-11`. Add `allowJs`/`checkJs` TypeScript configuration and `tsc --noEmit` to CI.

48. ✅ **Möbius pole guards use inconsistent thresholds.** `daydream/tools/mobius.html:369-380`. One branch checks `|c| < 1e-6`, while `cdiv` checks `|c|² < 1e-6` (`|c| < 1e-3`), causing a discontinuous zero result between thresholds. Compare the same squared magnitude or use one shared epsilon.

49. ✅ **The copied daydream README has broken repository-relative links and one false tool claim.** `daydream/README.md:2038` plus links to Holosphere-only hook, CMake, toolchain, Teensy spec and smoke script paths. Generate repository-aware links or use absolute Holosphere URLs; remove the nonexistent Lissajous amplitude-slider claim.

50. ✅ **daydream screenshot installation leaves stale and weak gallery assets.** `daydream/docs/screenshots/`. `Metaballs.png`, `NoiseRings.png` and `SplineFlow.png` are unreferenced leftovers; several temporal captures are overly dark/sparse at gallery scale and `HankinSolids` is not representative. Delete destination orphans and use deterministic per-effect capture offsets. Hankin timing at 2, 30 and 32 seconds was visually phase-stable, so no unsupported replacement was committed.

51. ✅ **Committed import-map drift is not directly gated.** `daydream/scripts/generate-importmap.mjs`, `daydream/tests/generate_importmap.test.js`, `daydream/vendor-importmap.js`. Tests prove generation behavior but do not regenerate and diff the committed deployment map. Run the default generator in CI and require a clean tree.

52. ✅ **Expected worker-fault logging overwhelms test output.** `daydream/tests/segment_controller.test.js` and related worker tests. Capture/stub expected logging so unexpected CI diagnostics remain visible.

53. **Contributor and security governance is missing.** There is no `CONTRIBUTING`, security policy, release/version policy or issue template. Add concise guidance and make the noncommercial-core/proprietary-effects scope prominent for potential contributors.

## Rejected candidates worth noting

- WASM heap-growth view detachment is detected and healed.
- Old/duplicate worker frames are generation-fenced and cannot settle a frame twice.
- ISR/DMA buffer ownership and relaxed completion ordering are correct for the current publication design.
- Effect handoff release/ack generation matching prevents stale use-after-free.
- Spherical seam, pole, antipode and CSG interval-capacity concerns are extensively bounded and tested.
- `FunctionRef` temporary borrowing is separated from retained `StoredFunctionRef` use.
- The current KiCad netlist is electrically connected and ERC/DRC clean; the findings concern physical layout, manufacturing safety and documentation.
- Sparse Holosphere source screenshots are valid compositions when viewed at full size; daydream’s three extra PNGs are separate stale destination assets.

## Recommended repair order

1. Fix the Phantasm allocation model and enforce a RAM2 heap floor.
2. Add mesh-length fail-fast checks and correct physical rotor attachment/sync-node layout.
3. Restore green/trustworthy release gates: CI shard coverage, DRC packaging, size-parser failure, BOM/CPL completeness.
4. Correct the unsafe circuit diagram and assembly markings before any PCB fabrication.
5. Make daydream bootstrap/switch/worker/recorder failure paths transactional and leak-free.
6. Add browser E2E and Teensy HIL coverage at the two boundaries unit tests cannot model.
7. Repair reproducibility, documentation freshness and accessibility/tool interaction gaps.
