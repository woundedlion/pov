# Documentation and comment fix report

Review date: 2026-10-01. Requested scope: documentation/comment drift beyond the existing code review, including excessive comments under the code-style skill. Original audit evidence is retained below; the completed fixes are recorded here.

## Resolution — 2026-10-01

All 13 accuracy groups and 84 comment-style groups are closed. Fixes landed through isolated implement worktrees and fast-forward integration, in one commit per repository:

- Holosphere: `b26a7e3a8dbb2085655ced20ddf106518d1b5983` — 56 files.
- daydream: `7e09e7fdde98c8d513d0198b729fd28e9737caaf` — 36 files.

C001's `CMakeLists.txt:309–331` subrange was a false positive: those comments are inside emitted CMake strings, which the code-style skill protects. That subrange is withdrawn and preserved unchanged; the remaining C001 source comments were fixed. This correction does not withdraw the C001 group.

The D003 correction also updates `setDisplayCaps`'s Doxygen contract in `targets/wasm/engine_bindings.h`. Final source review confirmed the authoring tables include azimuthal kaleidoscope symmetry and scope transition-pause behavior to composed effects.

Validation passed: executable-source/token identity across 50 changed Holosphere source files and all 35 changed daydream JavaScript files; staged commit hooks; 124 tracked Markdown files with sibling links supplied; 42 image references; Doxygen 1.17.0 with warnings treated as errors; daydream lint and typecheck. No runtime behavior changes were introduced.

The local simulator was rebuilt and installed from the clean isolated Holosphere checkout. Its source pin is exactly `b26a7e3a8dbb2085655ced20ddf106518d1b5983`, and both artifact checksums pass. Focused runtime probes confirmed retired-alias rejection, schema-2 snapshot restoration, and identical complete ShaderChain snapshots before and after a display-cap rebuild.

## Results and evidence baseline

**97 additional fix groups: 13 documentation/comment accuracy groups (D001–D013), and 84 comment-style groups (C001–C084).** Groups collect related locations; they are not a count of individual defective sentences. All are P4 documentation/style work. Accuracy corrections should precede style cleanup.

- Holosphere source snapshot: [d1b45d901bcec38c58599d7921699f8bc7cd8128](https://github.com/woundedlion/pov/tree/d1b45d901bcec38c58599d7921699f8bc7cd8128).
- daydream source snapshot: [f2cae14c37afb48b60cccb7006bad085ca1f8875](https://github.com/woundedlion/daydream/tree/f2cae14c37afb48b60cccb7006bad085ca1f8875).
- Exclusion baseline: all 261 numbered findings in the historical docs/CODE_REVIEW.md, including completed findings. Exclusions were checked by subject as well as location; this report does not reopen them.
- Evidence paths and line numbers below refer to these pinned source snapshots. The other fixing agent advanced both shared checkouts during this audit; re-anchor line numbers against the landing head before applying a fix.
- The live-comment comparison found the listed style targets still present. The changed whole-file headers in `analyze_candidates.py` and `profile_one.sh` retain the independently flagged rationale paragraphs; their existing-review fixes are excluded.

## Scope and review method

Inventoried 1,087 tracked first-party text files across both repositories and extracted 21,369 comment/docstring groups. Reviewed current user/developer documentation, subsystem/effect references, specs and their implementation banners; searched comments and docstrings for historical names, version claims, absolute scope claims, obsolete restrictions, restatement, numbered narration and defensive/performance essays. Checked candidates against implementation, catalog fields, existing tests and the original finding list. Script/tool comments and specs were included even though the original review excluded those areas.

This is a broad semantic audit with source-backed findings, not a claim that a heuristic comment extractor proves every remaining sentence correct. The extractor is a search aid, not a complete C++/JavaScript parser. Vendored/generated content, licenses, directives, actionable TODO/FIXME/HACK notes, commented-out code and text inside strings/template literals are not cleanup targets. Untracked board candidates and other agents' working files were excluded. Historical measurements and explicitly unimplemented or superseded designs were interpreted according to their status, rather than compared as current promises.

Applied the local code-style and Holosphere workflow skills. In particular, “Comments state **facts**, tersely” and “A genuinely non-obvious correctness gotcha … **One terse line**.” Descriptive Doxygen/JSDoc remains API reference; length alone is not a finding.

## Accuracy findings

### D001 — README advertises three rejected effect aliases

**Locations:** Holosphere README.md:1087

**Mismatch:** The setEffect row calls `Shader`, `ShaderBall` and `ShaderWorkbench` aliases for `ShaderChain`. These names are rejected.

**Evidence:** Holosphere targets/wasm/effect_factory.h:65-75 resolves only a roster name or stable ID. Holosphere scripts/wasm_smoke.mjs:626-628 explicitly requires UNKNOWN_EFFECT for all three retired names.

**Fix:** Remove the alias sentence. Document the current ShaderChain identity and, if useful, say retired names return UNKNOWN_EFFECT while keeping the prior effect. Do not restore compatibility aliases in code.

### D002 — README gives the wrong snapshot schema version

**Locations:** Holosphere README.md:1095

**Mismatch:** The getSnapshot row calls the result a “complete version-one program.” The executable snapshot schema is version 2.

**Evidence:** Holosphere workbench/shader/chain_snapshot.h:36-37 sets SCHEMA_VERSION to 2; Holosphere docs/specs/chain_snapshot_spec.md:18 documents schemaVersion: 2.

**Fix:** Replace version-one with version-two/current-schema, and link the snapshot contract. Do not confuse snapshot schema 2 with the v1 names of individual runtime payload kinds.

### D003 — Geometry rebuild documentation loses the ShaderChain preservation exception

**Locations:** Holosphere README.md:1119; Holosphere docs/specs/display_geometry.md:64-69; daydream docs/development.md:85-87

**Mismatch:** The documents describe animation/trail or simulation history resetting during an effect rebuild without describing the complete runtime restored for ShaderChain.

**Evidence:** Holosphere targets/wasm/engine_bindings.h:1018-1035 captures adapter rebuild state and restores it after replacement. Holosphere targets/wasm/workbench_bindings.h:150-163 captures and restores chain.snapshot(). Holosphere docs/specs/chain_snapshot_spec.md:77-79 explicitly includes hidden palette clocks and noise seeds.

**Fix:** Distinguish geometry-dependent effect reconstruction from ShaderChain's full executable snapshot restoration. State that its program, parameters, typed runtime and palettes survive; newly built geometry caches do not constitute a reset of those clocks/seeds.

### D004 — README claims every design spec is outside Doxygen

**Locations:** Holosphere README.md:1460-1463

**Mismatch:** The blanket statement that design specs are outside the Doxygen reference disagrees with the configured inputs.

**Evidence:** Holosphere Doxyfile:13-17 explicitly includes chain_snapshot_spec.md and parameter_description_spec.md.

**Fix:** Describe the specs index separately while naming the two contracts included in the API reference, or say most specs remain outside it.

### D005 — Face-LUT header describes an unused facility as shipping behavior

**Locations:** Holosphere core/render/sdf/face.h:28-34

**Mismatch:** The header claims every Islamic mesh is near-canonical, MeshOps bakes one class LUT at spawn, and Scan::Mesh binds it per frame. IslamicStars does not wire this facility.

**Evidence:** Holosphere docs/specs/congruence_class_lut_spec.md:3-9 says FACILITY ONLY and no effect wiring. Holosphere effects/IslamicStars.h:313-315 specifies the exact face-distance path because ripple/segue deformation invalidates canonical LUT shading.

**Fix:** Replace the production claim with a short description of the optional static-mesh facility and its source-of-truth pointer. Do not imply all Islamic meshes satisfy the deformation restriction. This block is counted once as drift, not again as overkill.

### D006 — Landed morph spec describes a contact-angle sweep instead of baked-direction interpolation

**Locations:** Holosphere docs/specs/opchain_morph_spec.md:27

**Mismatch:** The HANKIN_SWEEP table describes movement “across the contact-angle sweep.” This suggests the leg interpolates the contact angle or re-solves the star at each angle.

**Evidence:** Holosphere core/animation/opleg.h:458-463 expressly says the drive is a slerp fraction, not a contact angle; Holosphere core/animation/opleg.h:1539-1560 interpolates corner directions toward baked arrival points.

**Fix:** Describe HANKIN_SWEEP as spherical interpolation from corner directions to baked arrival directions. Reserve contact-angle wording for constructing its endpoints, and keep the behavior distinct from a live Hankin re-solve.

### D007 — Authoring vocabulary still presents a fixed two-warp host

**Locations:** Holosphere docs/effects.md:650-662; Holosphere docs/effects.md:678-682; Holosphere docs/effects.md:694-697

**Mismatch:** The current authoring section presents universal Planar Warp 1/2, a global Surface Noise Placement and an always-registered Camera Wander. The chain workbench uses selected, labeled operator instances and does not require a camera or two warp slots.

**Evidence:** daydream src/workbench/shader/chain_document_store.js:67-116 derives parameters only from the supplied chain. The compiler admitted `sample.spherical-noise.v3 → colorize.generated-palette.v3` with no camera, projection or warp. Holosphere docs/specs/shader_workbench_chain_spec.md:54-79 defines ordered instance-bearing chains.

**Fix:** Separate the fixed ComposedEffect derivation vocabulary from the current chain editor. Describe arbitrary admitted operator instances, per-instance controls and explicit ordering; remove the universal two-slot/camera claims. Keep the fixed-effect derivation limits accurately labeled.

### D008 — Projected-noise whitelist and sphere-source projection contract are obsolete

**Locations:** Holosphere docs/effects.md:673-676

**Mismatch:** The section restricts projected noise to four projections, then says spherical noise/rings work with every projection and merely reject planar warps. In a typed chain, sphere sources consume Sphere directly; they cannot be placed after a Plane-producing projection.

**Evidence:** The compiler admitted `project.bonne.v3 → sample.projected-noise.v2 → colorize.generated-palette.v3`. It refused `project.stereographic.v2 → sample.spherical-noise.v3 → colorize.generated-palette.v3` with FAMILY_ORDER (“steps back from the plane carrier to sphere”). The carrier authority is Holosphere scripts/engine_catalog.json:1 and semantic validation is Holosphere scripts/shader_workbench.mjs:436-441.

**Fix:** Document source input carriers rather than the retired host whitelist. Sphere sources bypass projection; projected noise follows compatible Plane producers. Preserve specific edge-distance/bounds requirements without inventing a general ban on Bonne/Peirce/Airocean.

### D009 — Affine authoring docs promise snapping and cross-stage restrictions absent from the chain operator

**Locations:** Holosphere docs/effects.md:667

**Mismatch:** The paragraph derives period solely from Primitive Lattice, promises nearest-integer translation snapping, and requires that lattice source, no later warp and a particular hue mode. These are not the current warp.affine.v3 chain contract.

**Evidence:** Holosphere core/render/pullback/operators_warp.h:52-99 exposes an explicit lattice-period parameter and passes it to Warp::prepare. Holosphere core/render/pullback/warp.h:211-216 declares translation as float LERP; Holosphere core/render/pullback/warp.h:399-407 uses it without rounding. The compiler admitted a speed-0.005 affine stage with translation-x 0.5 followed by sample.rings.v2.

**Fix:** Describe explicit lattice-period units/range and fractional translations for the chain operator. If the whole-cell seamless-scroll restrictions remain useful for a fixed composed derivation, label that scope explicitly instead of promising chain-level snapping or rejection.

### D010 — Polar Chart docs promise restrictions the typed chain does not enforce

**Locations:** Holosphere docs/effects.md:704

**Mismatch:** The paragraph requires Polar Chart to be the only planar warp except for a following Wave Shear, requires Grid/Lattice, and promises whole-seam product admission checks.

**Evidence:** The current compiler admitted `project.stereographic.v2 → warp.polar-chart.v2 → warp.vortex.v2 → sample.rings.v2 → colorize.generated-palette.v3`. Holosphere core/render/pullback/operators_warp.h:264-306 defines the Plane endomorphism and its own mode/harmonic checks.

**Fix:** Describe typed composition admission separately from conditions needed for a visually seamless periodic pattern. Label any surviving fixed-derivation constraints by scope. Do not introduce new validation or silently narrow the engine to make the old prose true.

### D011 — Authoring docs say incompatible chain edits stay pending

**Locations:** Holosphere docs/effects.md:704; Holosphere docs/effects.md:708

**Mismatch:** The documentation says incompatible stage combinations remain pending until repaired. The chain document store discards a refused candidate and retains the prior committed document.

**Evidence:** daydream src/workbench/shader/chain_document_store.js:247-276 validates before changing doc or history; daydream src/workbench/shader/chain_document_store.js:603-621 routes preset writes through that commit.

**Fix:** Describe transactional refusal and retained committed state for the chain editor. Distinguish any requested/accepted value reporting on ordinary effect parameters from a pending structural chain. This is a documentation issue, not a new implementation of review finding 8 or the test gap in finding 217.

### D012 — Parameter-deactivation JSDoc wrongly gates depth on brightness

**Locations:** daydream src/workbench/shader/chain_presentation.js:55-64

**Mismatch:** The comment says brightness endpoints and depth require a brightness envelope. The current generated-color operator gates brightness-bottom/top but not phase-oscillation-depth; there is no brightness-depth parameter in its schema.

**Evidence:** daydream src/workbench/shader/chain_presentation.js:68-83 reads each catalog field's gated_by. The generated-palette entry in Holosphere scripts/engine_catalog.json:1 has no gated_by on phase-oscillation-depth.

**Fix:** Remove the depth example and preferably describe catalog-driven deactivation generically. Keep the documented distinction between stored parameter values and dimmed controls. This is different from original finding 205's stageControlLabel branch.

### D013 — Specs direct readers to removed supporting artifacts

**Locations:** Holosphere docs/specs/spherical_perspective_spec.md:8-14; Holosphere docs/specs/hyperlattice_octet_58ms_spec.md:6-9; Holosphere docs/specs/hyperlattice_octet_58ms_spec.md:65-68

**Mismatch:** The spherical spec says “See” an admission report and describes a removed preset guide as the place recording current controls. The Octet spec points to a removed optimization ledger and says the final disassembly “is preserved” in removed evidence. Those are unusable present-tense references, despite the explicit removed-artifact marker.

**Evidence:** The cited paragraphs themselves label every supporting artifact removed. Retained reports are indexed under Holosphere docs/profiles/README.md:1; the Octet spec also has working links to retained preset captures.

**Fix:** Replace actionable pointers with an existing document that actually contains the required data, summarize the missing control/admission facts in the spec, or state that the evidence is no longer retained. Remove “is preserved” if no retained disassembly exists. Keep historical numbers labeled as historical; do not fabricate a replacement capture. Excludes original finding 44's different orphan sentence in docs/effects.md:288.

## Comment-style findings

Each C item is an independent cleanup group, with exact snapshot ranges and a concrete edit boundary. Ranges locate the offending comments and may also contain surrounding code: **edit comment text only, never delete the range wholesale**. Preserve license/directive lines and every API parameter, return, lifetime or admission contract that still earns its place. Corrections below explicitly retain non-obvious constraints rather than deleting them because the surrounding explanation is long.

### C001 — Holosphere: CMakeLists.txt

**Locations:** Holosphere CMakeLists.txt:1-5, Holosphere CMakeLists.txt:64-68, Holosphere CMakeLists.txt:89-111, Holosphere CMakeLists.txt:151-156, Holosphere CMakeLists.txt:213-221, Holosphere CMakeLists.txt:239-244, Holosphere CMakeLists.txt:309-331

**Style violation:** Build comments include performance anecdotes, hypothetical regressions and repeated option explanations.

**Fix:** Keep the CMake floor, actual CSP/JIT restrictions, stack budgets, explicit-target instructions, install override and package provenance contracts in brief labels; delete the defenses and move benchmark detail to documentation. Exclude the warning-gate paragraph already covered by finding 192.

### C002 — Holosphere: cmake/toolchain-native-clang.cmake

**Locations:** Holosphere cmake/toolchain-native-clang.cmake:43-51

**Style violation:** The linker alias explanation is a walkthrough of the setup and its justification.

**Fix:** Keep one factual sentence identifying the lld multicall alias and search path.

### C003 — Holosphere: core/animation/opleg.h

**Locations:** Holosphere core/animation/opleg.h:632-639, Holosphere core/animation/opleg.h:831-841, Holosphere core/animation/opleg.h:1118-1122, Holosphere core/animation/opleg.h:1787-1792

**Style violation:** Packed-position and leg-timing comments contain implementation walkthroughs and performance explanations.

**Fix:** Keep only classifier precision, settle-before-sweep ordering, palette ownership and symmetry invariants; move measured tradeoffs to the morph spec.

### C004 — Holosphere: core/animation/recipe_build.h

**Locations:** Holosphere core/animation/recipe_build.h:240-244, Holosphere core/animation/recipe_build.h:263-269, Holosphere core/animation/recipe_build.h:817-821, Holosphere core/animation/recipe_build.h:841-845

**Style violation:** The comments narrate recipe lowering, endpoint construction, reclamation and retained palette identities.

**Fix:** Replace each with its non-obvious endpoint invalidation or prefix-identity rule; omit the step-by-step narration.

### C005 — Holosphere: core/animation/sprites.h

**Locations:** Holosphere core/animation/sprites.h:73-78

**Style violation:** The overlapping fade-duration comment walks through the normalization algorithm.

**Fix:** Keep a terse statement that overlapping fades scale proportionally and normalization uses the member duration.

### C006 — Holosphere: core/animation/timeline.h

**Locations:** Holosphere core/animation/timeline.h:399-407, Holosphere core/animation/timeline.h:473

**Style violation:** The collapse block explains each pass; the start-time label repeats the immediately following check.

**Fix:** Keep the shared-orientation collapse ordering invariant in one sentence; delete the check label.

### C007 — Holosphere: core/animation/transformer.h

**Locations:** Holosphere core/animation/transformer.h:342, Holosphere core/animation/transformer.h:524-530, Holosphere core/animation/transformer.h:829-833

**Style violation:** There is an obvious sync label and expanded explanations of pointer pinning and channel decorrelation.

**Fix:** Delete the sync label; keep the retained-pointer lifetime and channel-offset rule briefly.

### C008 — Holosphere: core/color/srgb_decode.h

**Locations:** Holosphere core/color/srgb_decode.h:52-59

**Style violation:** The LUT-placement comment mixes required placement facts with an extended performance rationale.

**Fix:** Keep only the DTCM placement and initialization constraints; move cache/performance discussion to the subsystem documentation.

### C009 — Holosphere: core/color/composition.h

**Locations:** Holosphere core/color/composition.h:57-68, Holosphere core/color/composition.h:406-414

**Style violation:** API notes repeatedly defend allowing a value outside the unit interval.

**Fix:** Preserve the API brief and return contract; state wrap behavior once and remove the justification.

### C010 — Holosphere: core/control/params.h

**Locations:** Holosphere core/control/params.h:211-216

**Style violation:** The private-overload comment explains consequences rather than stating the access constraint.

**Fix:** Keep the overload-resolution/access constraint in one sentence.

### C011 — Holosphere: core/engine/memory.cpp

**Locations:** Holosphere core/engine/memory.cpp:154-158

**Style violation:** The persistent-storage rebind comment narrates why the implementation preserves live contents.

**Fix:** Keep the content/generation preservation invariant briefly.

### C012 — Holosphere: core/engine/memory.h

**Locations:** Holosphere core/engine/memory.h:23-34, Holosphere core/engine/memory.h:55-59

**Style violation:** Budget comments contain host/device footprint rationale and a defense of the WASM budget.

**Fix:** Keep units, actual budgets and the pointer-size caveat; move detailed sizing discussion to the memory documentation.

### C013 — Holosphere: core/engine/static_storage.cpp

**Locations:** Holosphere core/engine/static_storage.cpp:58-63, Holosphere core/engine/static_storage.cpp:91-96

**Style violation:** Terminate-handler comments contain demangler/code-size rationale and a runtime-installation essay.

**Fix:** Keep the device demangler suppression and Teensy handler override contract in terse factual comments.

### C014 — Holosphere: core/mesh/conway.h

**Locations:** Holosphere core/mesh/conway.h:426-477

**Style violation:** The operator contract block mixes valuable arena/manifold rules with a long digression about why dropping degenerate faces is load-bearing.

**Fix:** Preserve the allocation, arena-distinctness, manifold and output-location contracts; shorten the dropped-face consequence to a factual rule and move the Hankin/ambo explanation to the morph spec.

### C015 — Holosphere: core/render/filter.h

**Locations:** Holosphere core/render/filter.h:7-14

**Style violation:** The file documentation defends why the wider filter roster is not dead code.

**Fix:** Keep the file brief; remove the defensive details paragraph.

### C016 — Holosphere: core/render/canvas.h

**Locations:** Holosphere core/render/canvas.h:508-516

**Style violation:** The prev_data API details justify constness already expressed by the signature and return contract.

**Fix:** Keep the pointer layout, constness and lifetime docs; remove the defense of the return type.

### C017 — Holosphere: core/render/sdf/face.h

**Locations:** Holosphere core/render/sdf/face.h:468-475

**Style violation:** The collapsed-face comment adds a long quantitative explanation to an empty-face guard.

**Fix:** Keep the area/degeneracy condition briefly; move the numerical argument to tests or the spec.

### C018 — Holosphere: core/render/sdf/common.h

**Locations:** Holosphere core/render/sdf/common.h:112-118, Holosphere core/render/sdf/common.h:203-209

**Style violation:** Emission and buffer comments expand their non-obvious contracts into caller walkthroughs.

**Fix:** Keep the no-emission-on-false rule and capacity contract without the narrative.

### C019 — Holosphere: core/render/sdf/volume.h

**Locations:** Holosphere core/render/sdf/volume.h:514-519

**Style violation:** The frame-preparation comment explains optimizations and repeats the twist-zero special case.

**Fix:** Use a short preparation label; omit the routine optimization rationale.

### C020 — Holosphere: core/render/sdf/csg.h

**Locations:** Holosphere core/render/sdf/csg.h:531-537

**Style violation:** The seam-handling comment walks through the merge procedure and its consequences.

**Fix:** Keep the common-frame seam split and doubled-buffer capacity invariant tersely.

### C021 — Holosphere: core/render/plot/raster.h

**Locations:** Holosphere core/render/plot/raster.h:472-480, Holosphere core/render/plot/raster.h:557-563, Holosphere core/render/plot/raster.h:870-874

**Style violation:** Raster comments include a measured cache bound and expanded explanations of normalization.

**Fix:** Keep capacity/backstop and unit-fragment invariants; move empirical measurements out of the implementation.

### C022 — Holosphere: core/render/scan/raster.h

**Locations:** Holosphere core/render/scan/raster.h:453, Holosphere core/render/scan/raster.h:661-666

**Style violation:** An explicit clamp label repeats the following checks; the scratch comment defends memory placement at length.

**Fix:** Delete the clamp label; retain the arena ownership and normal-buffer capacity facts briefly.

### C023 — Holosphere: core/render/scan/shapes.h

**Locations:** Holosphere core/render/scan/shapes.h:686-691

**Style violation:** The ring-walk comment defends the row-local algorithm with a performance walkthrough.

**Fix:** Keep a short label for the covering ring's row runs.

### C024 — Holosphere: effects/HankinSolids.h

**Locations:** Holosphere effects/HankinSolids.h:141-146, Holosphere effects/HankinSolids.h:209-216

**Style violation:** Duration and API details justify choreography tuning and relative palette weights.

**Fix:** Keep frame units and the factual choreography contract; remove the tuning rationale while retaining the API brief/parameters.

### C025 — Holosphere: effects/RingSpin.h

**Locations:** Holosphere effects/RingSpin.h:121-122

**Style violation:** The speed-constant comment defends keeping the value above a floor.

**Fix:** Delete the ordinary-constant rationale.

### C026 — Holosphere: effects/Thrusters.h

**Locations:** Holosphere effects/Thrusters.h:98-99

**Style violation:** The draw-loop comment repeats the advancement and drawing performed immediately below.

**Fix:** Delete the restatement.

### C027 — Holosphere: hardware/hd107s_frame.h

**Locations:** Holosphere hardware/hd107s_frame.h:78-111

**Style violation:** The frame/descriptor documentation includes a numbered color-correction procedure and a hypothetical capacity failure explanation.

**Fix:** Keep protocol layout and the eDMA 15-bit iteration bound; describe color-correction order briefly without numbered narration.

### C028 — Holosphere: hardware/pov_segmented.h

**Locations:** Holosphere hardware/pov_segmented.h:444-448, Holosphere hardware/pov_segmented.h:495-524

**Style violation:** Logging includes a code-size anecdote; strap documentation expands debounce limitations into numbered failure stories.

**Fix:** Keep log-buffer capacity and factual strap/debounce limitations in brief comments and API details.

### C029 — Holosphere: hardware/pov_sync.h

**Locations:** Holosphere hardware/pov_sync.h:733-737, Holosphere hardware/pov_sync.h:743-753

**Style violation:** The beacon/deadline comments narrate scheduling and defend skipped emission.

**Fix:** Keep only boundary precedence and the tail-before-HALF deadline invariant.

### C030 — Holosphere: hardware/pov_sync_protocol.h

**Locations:** Holosphere hardware/pov_sync_protocol.h:360-397

**Style violation:** Protocol constant comments repeatedly derive and defend the chosen margins and timeout.

**Fix:** Keep units and the few required spacing/timeout relationships; move the full numerical derivation to the sync specification.

### C031 — Holosphere: targets/Phantasm/phantasm_playlist.h

**Locations:** Holosphere targets/Phantasm/phantasm_playlist.h:168-174

**Style violation:** The roster guard comment defends the check with hypothetical selection failures.

**Fix:** Keep a short label identifying the segmented roster admission check.

### C032 — Holosphere: targets/wasm/effect_factory.h

**Locations:** Holosphere targets/wasm/effect_factory.h:78-87, Holosphere targets/wasm/effect_factory.h:95-96, Holosphere targets/wasm/effect_factory.h:114-117

**Style violation:** Registry expansion comments explain why shared generation prevents drift and narrate ordinary table construction.

**Fix:** Keep a terse authority pointer and the API brief/return contract; remove the defenses and obvious construction labels.

### C033 — Holosphere: targets/wasm/mesh_ops_bindings.h

**Locations:** Holosphere targets/wasm/mesh_ops_bindings.h:68-79, Holosphere targets/wasm/mesh_ops_bindings.h:629-643

**Style violation:** Scratch ownership and relaxation API details use a long reentrancy explanation and hypothetical failure scenarios.

**Fix:** Keep module-global non-overlap/lifetime constraints and the admitted iteration clamp; remove the narrative while retaining API parameter/return docs.

### C034 — Holosphere: hardware/phantasm/gen/kicad_common.py

**Locations:** Holosphere hardware/phantasm/gen/kicad_common.py:1-6

**Style violation:** The module docstring defends centralizing helpers so a schema change touches one place.

**Fix:** Keep the factual helper roster; delete the organizational justification.

### C035 — Holosphere: hardware/phantasm/gen/sexp.py

**Locations:** Holosphere hardware/phantasm/gen/sexp.py:48-58

**Style violation:** The discovery API docstring explains hypothetical foreign-geometry consequences and why import must not fail.

**Fix:** Keep discovery precedence, major-version filtering, override semantics and missing-path behavior in a short descriptive contract.

### C036 — Holosphere: hardware/phantasm/gen/analyze_candidates.py

**Locations:** Holosphere hardware/phantasm/gen/analyze_candidates.py:23-32

**Style violation:** The module docstring's What matters here, and why section is an electrical/layout rationale essay.

**Fix:** Move that rationale to hardware documentation; retain the tool's scoring criteria, invocation and admission rules.

### C037 — Holosphere: tools/device_lock.sh

**Locations:** Holosphere tools/device_lock.sh:1-44

**Style violation:** The usage header adds a Why lock at all essay and repeats lock-method discussion.

**Fix:** Keep commands, environment variables, lock scope and limitations; remove the defensive rationale. Do not refile the portless-fallback wording covered by finding 175.

### C038 — Holosphere: tools/profile_one.sh

**Locations:** Holosphere tools/profile_one.sh:13-25

**Style violation:** The profiling header duplicates the device-lock explanation at length.

**Fix:** Keep lock scope and the port override; replace the repeated loader/selection failure story with a source pointer. Exclude the broken wrap at line 5 already covered by finding 185.

### C039 — Holosphere: tools/teensy_size_trail.py

**Locations:** Holosphere tools/teensy_size_trail.py:2-33

**Style violation:** The module reference mixes storage/subcommand contracts with a branch-conflict and failure story.

**Fix:** Keep the shared-history versus per-worktree pending paths and two-phase admission rule; remove the narrative defense.

### C040 — Holosphere: tools/teensy_gate.py

**Locations:** Holosphere tools/teensy_gate.py:2-34

**Style violation:** The gate reference repeatedly explains how an otherwise-green build could silently lose protection.

**Fix:** Keep the ceiling/schema/failure contracts and CLI documentation; remove hypothetical regression stories.

### C041 — Holosphere: tools/build_pins.py

**Locations:** Holosphere tools/build_pins.py:350-354

**Style violation:** The API docstring narrates exception/rename handling rather than simply describing recorded errors.

**Fix:** Retain a concise description of the check result and error recording.

### C042 — Holosphere: scripts/capture_screenshots.mjs

**Locations:** Holosphere scripts/capture_screenshots.mjs:1-13, Holosphere scripts/capture_screenshots.mjs:134-143

**Style violation:** The tool header uses a numbered walkthrough, and effect resolution includes an extended wrong-image failure story.

**Fix:** Keep CLI usage and admitted-resolution semantics; remove numbered narration and the defensive scenario.

### C043 — Holosphere: scripts/wasm_smoke.mjs

**Locations:** Holosphere scripts/wasm_smoke.mjs:290-293, Holosphere scripts/wasm_smoke.mjs:543-544, Holosphere scripts/wasm_smoke.mjs:1038-1045, Holosphere scripts/wasm_smoke.mjs:1148-1152

**Style violation:** Enum checks repeatedly justify their existence; the growth probe narrates what a weaker arrangement would miss.

**Fix:** Keep enum success distinctions and the growth-before-first-mesh ordering constraint briefly; delete the test defenses.

### C044 — Holosphere: tests/check_pullback_capture.cmake

**Locations:** Holosphere tests/check_pullback_capture.cmake:1-9

**Style violation:** The gate header argues that a wrong stable run could pass the comparison.

**Fix:** Keep its capture-contract scope and required inputs in a short header.

### C045 — Holosphere: tests/check_vendor_patches.cmake

**Locations:** Holosphere tests/check_vendor_patches.cmake:1-9

**Style violation:** The patch gate header narrates hypothetical compiling-but-bloated builds.

**Fix:** Keep the vendor patch contract and the inputs it checks; remove the regression story.

### C046 — Holosphere: tests/check_linker_patches.cmake

**Locations:** Holosphere tests/check_linker_patches.cmake:1-8

**Style violation:** The linker patch gate header defends its pinning against hypothetical template drift.

**Fix:** Keep what is pinned and which framework inputs are required.

### C047 — Holosphere: tests/CMakeLists.txt

**Locations:** Holosphere tests/CMakeLists.txt:104-110

**Style violation:** The host-diagnostic comment argues for the inherited warning gate at length.

**Fix:** Keep the target/flag inheritance fact; remove the defense.

### C048 — Holosphere: tests/test_effects.h

**Locations:** Holosphere tests/test_effects.h:4699-4707

**Style violation:** The white-box section introduction explains why a smoke test alone would not prove bounds.

**Fix:** Keep the factual bounded-output scope and remove the test justification.

### C049 — daydream: src/app/app_lifecycle.js

**Locations:** daydream src/app/app_lifecycle.js:407-413

**Style violation:** The load timeout comment defends the chosen deadline with cold-start anecdotes.

**Fix:** Keep the watchdog's scope and deadline; delete the tuning rationale.

### C050 — daydream: src/app/bootstrap.js

**Locations:** daydream src/app/bootstrap.js:59-63

**Style violation:** The hot-reload timeout comment defends its relationship to the segment-load timeout.

**Fix:** Keep a short label for the reload sweep deadline.

### C051 — daydream: src/app/daydream.js

**Locations:** daydream src/app/daydream.js:339

**Style violation:** The renderer comment repeats constructor ordering and the now-available canvas.

**Fix:** Delete the obvious constructor narration.

### C052 — daydream: src/segments/module_warmer.js

**Locations:** daydream src/segments/module_warmer.js:13-17

**Style violation:** The warming deadline comment expands best-effort waiting into an ordinary-constant rationale.

**Fix:** Keep the best-effort deadline behavior briefly.

### C053 — daydream: src/segments/segment_controller.js

**Locations:** daydream src/segments/segment_controller.js:64-79, daydream src/segments/segment_controller.js:1091-1112, daydream src/segments/segment_controller.js:1147-1152, daydream src/segments/segment_controller.js:1263-1268, daydream src/segments/segment_controller.js:1288-1292

**Style violation:** Retry, reconstruction, transfer and publication comments mix useful ownership rules with speculative causes and implementation walkthroughs.

**Fix:** Keep the bounded retry/reset rule, latch ownership, buffer-transfer ownership, replay rule and publication invariant tersely; omit the stories. Exclude the watchdog and first-fault paragraphs already covered by findings 206 and 207.

### C054 — daydream: src/segments/segment_compositor.js

**Locations:** daydream src/segments/segment_compositor.js:108-112

**Style violation:** The mismatch comment defends retaining prior output and explains downstream fault routing.

**Fix:** Keep the buffer-identity refusal contract briefly.

### C055 — daydream: src/recording/recorder.js

**Locations:** daydream src/recording/recorder.js:624-628, daydream src/recording/recorder.js:639-643

**Style violation:** Cancellation and chunk serialization comments walk through race scenarios.

**Fix:** Keep the session-ownership and lazy-file-creation ordering invariants in brief comments.

### C056 — daydream: src/renderer/driver.js

**Locations:** daydream src/renderer/driver.js:497-501, daydream src/renderer/driver.js:517

**Style violation:** Orbit fitting is described step by step; the DPR label repeats the immediately following refresh.

**Fix:** Keep the user-zoom preservation rule and delete the DPR restatement.

### C057 — daydream: src/workbench/mobius/mobius_page.js

**Locations:** daydream src/workbench/mobius/mobius_page.js:291-296

**Style violation:** Cached labels are defended with mutation counts and a 60 Hz performance anecdote.

**Fix:** Delete the rationale; keep a short cache-purpose label only if needed.

### C058 — daydream: src/workbench/palettes/palettes_page.js

**Locations:** daydream src/workbench/palettes/palettes_page.js:93-97, daydream src/workbench/palettes/palettes_page.js:555-559, daydream src/workbench/palettes/palettes_page.js:662-666

**Style violation:** Drawn hit targets, accessible labels and seed edits are explained through multi-line scenarios.

**Fix:** Keep the drawn-recipe identity, group/channel accessible-name and seed-update timing contracts briefly.

### C059 — daydream: src/workbench/solids/solid_build.js

**Locations:** daydream src/workbench/solids/solid_build.js:180-186

**Style violation:** Mesh readback documentation expands the returned-buffer ownership/lifetime constraint into a walkthrough.

**Fix:** Keep the ownership and invalidation contract in a concise API description.

### C060 — daydream: src/workbench/solids/solids_page.js

**Locations:** daydream src/workbench/solids/solids_page.js:342-348, daydream src/workbench/solids/solids_page.js:1150-1158, daydream src/workbench/solids/solids_page.js:1229-1234, daydream src/workbench/solids/solids_page.js:1382-1386

**Style violation:** Thumbnail code recounts the earlier frozen first paint; input, revision and halted-module comments give lengthy defenses.

**Fix:** Describe current scheduling only; keep actual input/revision/terminal-handle constraints tersely. Exclude the saved-shape validator wording covered by finding 213.

### C061 — daydream: src/workbench/solids/solids_page.js

**Locations:** daydream src/workbench/solids/solids_page.js:120, daydream src/workbench/solids/solids_page.js:314, daydream src/workbench/solids/solids_page.js:386, daydream src/workbench/solids/solids_page.js:421, daydream src/workbench/solids/solids_page.js:939, daydream src/workbench/solids/solids_page.js:1477, daydream src/workbench/solids/solids_page.js:1540

**Style violation:** Initialize, Create, Update and Check comments repeat clearly named calls or predicates.

**Fix:** Delete these obvious operation labels.

### C062 — daydream: src/workbench/palettes/palette_canvas.js

**Locations:** daydream src/workbench/palettes/palette_canvas.js:219

**Style violation:** The dashed-line reset comment repeats setLineDash([]).

**Fix:** Delete the restatement.

### C063 — daydream: scripts/generate-importmap.mjs

**Locations:** daydream scripts/generate-importmap.mjs:1-29

**Style violation:** The reference header mixes invocation/pinning requirements with an extended network/cache rationale.

**Fix:** Preserve usage, pin inputs and output contracts; shorten the narrative paragraphs.

### C064 — daydream: scripts/lissajous-probe.mjs

**Locations:** daydream scripts/lissajous-probe.mjs:1-16

**Style violation:** The probe header explains at length why fake geometry is insufficient.

**Fix:** Keep exercised behavior and invocation; move test-strategy rationale to development documentation.

### C065 — daydream: scripts/mobius-probe.mjs

**Locations:** daydream scripts/mobius-probe.mjs:1-15

**Style violation:** The probe header explains at length why mocks cannot replace a real browser.

**Fix:** Keep exercised behavior and invocation; remove the rationale essay.

### C066 — daydream: scripts/palettes-probe.mjs

**Locations:** daydream scripts/palettes-probe.mjs:1-15

**Style violation:** The probe header expands the real-browser need into a mock/layout explanation.

**Fix:** Keep exercised behavior and invocation; remove the rationale essay.

### C067 — daydream: scripts/panel-probe.mjs

**Locations:** daydream scripts/panel-probe.mjs:1-21

**Style violation:** The probe header narrates layout failures and why unit mocks miss them.

**Fix:** Keep exercised behavior and invocation; remove the defensive walkthrough.

### C068 — daydream: scripts/solids-probe.mjs

**Locations:** daydream scripts/solids-probe.mjs:1-14

**Style violation:** The probe header explains why fake DOM interactions do not prove browser behavior.

**Fix:** Keep exercised behavior and invocation; remove the rationale essay.

### C069 — daydream: scripts/vendor-stage.mjs

**Locations:** daydream scripts/vendor-stage.mjs:1-15

**Style violation:** The staging header includes a hypothetical CDN incident/failure rationale.

**Fix:** Keep actual staging source, destination and invocation contracts; delete the defensive prose while retaining the shared-inode ownership constraint.

### C070 — daydream: scripts/record-module-loads.mjs

**Locations:** daydream scripts/record-module-loads.mjs:1-9

**Style violation:** The module-load header explains how its roster defense prevents hidden omissions.

**Fix:** Keep which modules are recorded and the NODE_TEST_CONTEXT filtering contract briefly.

### C071 — daydream: tests/segment_crosscheck.test.js

**Locations:** daydream tests/segment_crosscheck.test.js:1-32

**Style violation:** The file header is a long numbered explanation of why the parity oracle proves the contract.

**Fix:** Keep a concise description of shared coordinate conventions and the source pointer; delete the proof essay and numbered narration.

### C072 — daydream: tests/three_contract.test.js

**Locations:** daydream tests/three_contract.test.js:1-11, daydream tests/three_contract.test.js:236-237

**Style violation:** The header and constructor-name check explain hypothetical upstream changes and their failure consequences.

**Fix:** Keep the library contracts under test and remove the defense.

### C073 — daydream: tests/tsconfig_roster.test.js

**Locations:** daydream tests/tsconfig_roster.test.js:1-9

**Style violation:** The header argues how untyped modules could leave CI green.

**Fix:** Keep its roster/type coverage scope; omit hypothetical failure narration. This is separate from finding 257's parser/exemption wording.

### C074 — daydream: tests/chain_document_store.test.js

**Locations:** daydream tests/chain_document_store.test.js:1-9

**Style violation:** The test reference repeats the store's mutation model and argues why using real fixtures matters.

**Fix:** Keep the tested store contract and fixture authority briefly.

### C075 — daydream: tests/chain_strip.test.js

**Locations:** daydream tests/chain_strip.test.js:1-9

**Style violation:** The file header narrates editor/store/backend integration and its test setup.

**Fix:** Keep the tested surface and fixture source pointer; remove the walkthrough.

### C076 — daydream: tests/shared.test.js

**Locations:** daydream tests/shared.test.js:1-9

**Style violation:** The header elaborates disposal ordering into a test-design story.

**Fix:** Keep only the required ordering and fixture scope.

### C077 — daydream: tests/apply_notice.test.js

**Locations:** daydream tests/apply_notice.test.js:12-15

**Style violation:** The test comment defends a guard with a hypothetical renamed implementation.

**Fix:** Delete the test justification.

### C078 — daydream: tests/engine_contract_wasm.test.js

**Locations:** daydream tests/engine_contract_wasm.test.js:372-374, daydream tests/engine_contract_wasm.test.js:912-915

**Style violation:** Comments explain why weaker method/roster checks would pass missing behavior.

**Fix:** Keep the concrete assertions' contract labels and remove the weaker-test stories.

### C079 — daydream: tests/engine_source_parity.test.js

**Locations:** daydream tests/engine_source_parity.test.js:373-381

**Style violation:** The helper API docs expand field-path generation into a regression scenario.

**Fix:** Keep the helper's factual input/output and source-parity contract in a short API description.

### C080 — daydream: src/workbench/palettes/palette_controls.js

**Locations:** daydream src/workbench/palettes/palette_controls.js:844-859

**Style violation:** The recipe-marshalling API docs include a renamed-option failure story and defense of the DOM boundary.

**Fix:** Keep the descriptive brief, custom-curve endpoint rule, clone ownership, parameters, return and throws contracts; delete the defensive paragraph.

### C081 — daydream: tests/helpers/fake_dom.js

**Locations:** daydream tests/helpers/fake_dom.js:1-16

**Style violation:** The double's header includes a long list of defects that can leave unit tests green and repeats the browser-probe roster.

**Fix:** Keep the modeled surface and explicit limitations briefly, with one pointer to the browser probes.

### C082 — daydream: tests/driver.test.js

**Locations:** daydream tests/driver.test.js:1735-1739

**Style violation:** A helper comment defends mirroring the app's field roster with a renamed-state scenario.

**Fix:** Keep the API description and source pointer; delete the defensive preamble.

### C083 — daydream: tests/palette_controls.test.js

**Locations:** daydream tests/palette_controls.test.js:778-779

**Style violation:** The markup-parity comment narrates renamed members and refused exports instead of labeling the check.

**Fix:** Keep a short source-parity label or delete the comment.

### C084 — Holosphere: core/mesh/solids.h

**Locations:** Holosphere core/mesh/solids.h:499-504

**Style violation:** The BaseMesh-order comment adds a hypothetical compile-failure defense to a valuable ordinal mapping.

**Fix:** Keep the registry concatenation and slice-boundary contract; delete the future insertion/removal/reordering story.

## Verification and deduplication

- `python tools/docs_check.py --checkout daydream=C:/work/temp/daydream-doc-drift-audit`: PASS, 124 tracked Markdown files.
- `python tools/docs_images.py`: PASS, 42 image references.
- Five document-compiler probes: four VALID cases and one expected FAMILY_ORDER refusal. These are document-admission probes, not device rendering or timing results.
- Matched style ranges against the original review's primary path/line references: no intersections after excluding the profile header wrap. Also reviewed same-file and same-subject exclusions manually.
- Style groups locate 135 source ranges. Several ranges contain more than one related comment; this count is not a line count or severity multiplier.
- No firmware builds, device operations, application edits, staging or commits were performed for this audit.

Document-admission reproduction: construct a schema-2 scratch document using the parameter/default mapping in `scratchChainDocument` (D007's source link), substitute each ordered chain below, and pass it to `compileShaderDocument(document, { catalog })` using the pinned `scripts/engine_catalog.json`. Labels were `stage-0`, `stage-1`, etc.; all other values were catalog defaults.

| Ordered operator IDs | Override | Observed compiler result |
|---|---|---|
| `sample.spherical-noise.v3 → colorize.generated-palette.v3` | None | VALID |
| `project.bonne.v3 → sample.projected-noise.v2 → colorize.generated-palette.v3` | None | VALID |
| `project.stereographic.v2 → sample.spherical-noise.v3 → colorize.generated-palette.v3` | None | INVALID / FAMILY_ORDER at chain[1] |
| `project.stereographic.v2 → warp.polar-chart.v2 → warp.vortex.v2 → sample.rings.v2 → colorize.generated-palette.v3` | None | VALID |
| `project.stereographic.v2 → warp.affine.v3 → sample.rings.v2 → colorize.generated-palette.v3` | `stage-1.translation-x = 0.5`, `stage-1.speed = 0.005` | VALID |

Explicitly excluded collateral repetitions of existing findings: board-revision pin names/sync prose (3), removed docs/effects artifact sentence (44), reaction-diffusion capacity docs (51), parameter unit disagreements (69), gated-swap IslamicStars claims (98), smoke clip history narration (113), removed scan identifiers (139), portless-lock fallback (175), warning baseline/ratchet wording (192), mirrored manifest header (204), stageControlLabel prose (205), pool watchdog/first-fault scope (206–207), custom hue count (209), expand morph wording (212), saved-shape validator return description (213), test mojibake/helper docs (227/232/233), worker outcome logging wording (244), wrong README reference (249), drag-slot direction wording (253), tsconfig parser/exemption wording (257), Tailwind specificity claim (259), and the tools.css wrap/consecutive blocks (261). C073 flags only the independent rationale at the top of the tsconfig test, not finding 257.

Rejected or intentionally retained candidates: current accepted_parameter_value override description (no shipping override found); current shader store API names (they exist); full-frame trait distinctions between pixel feedback and world trails; normative historical capture measurements; superseded conway_morph_spec symbols and constants; explicit unimplemented promotion design; LUT/device ABI source pointers; NaN, unsigned-wrap, lifetime and arena-rewind constraints; rasterizer mathematical cull bounds; and descriptive Doxygen/JSDoc contracts. These were not turned into style findings solely because they are long.

## Fix acceptance criteria

1. Apply D001–D013 as prose/comment corrections against the landing head. Recheck the named implementation authority; do not modify runtime behavior to accommodate obsolete wording.
2. Apply C001–C084 only to the identified narrative/restatement portions. Retain brief source pointers, ownership/UB/wrap contracts, API reference data, CLI usage and required tooling annotations.
3. Confirm source-bearing tokens are unchanged for comment cleanup, with a lexer-aware comparison or an exact diff review. Preserve strings/templates and emitted source. Run repository lint/docs checks appropriate to the touched files; no new behavior tests are needed for pure prose cleanup.
4. Rerun tracked Markdown/image checks and Doxygen/JS documentation checks where the edited reference surfaces warrant them. Treat structural green results as supplementary evidence, not proof of semantic accuracy.

This report preserves the completed audit; it is not a design specification and does not change the specs index.
