# Holosphere / daydream — Code Quality Review

**Scope reviewed:** the Holosphere C++ rendering engine + firmware and the daydream web simulator (both repositories, shipped as one product).
**Out of scope (per review charter):** `core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and all vendored third-party code (`core/vendor/FastNoiseLite*`, `daydream/three.js`, `node_modules`, the installed `holosphere_wasm.{js,wasm}`). Generated files (`color_luts.h`, `reaction_graph.cpp`) were assessed via their generator scripts, not nitpicked line-by-line.

**Method.** The two repositories were partitioned into 23 components (engine, math, mesh, color, render, animation, 28 effects, hardware drivers, WASM/build/CI, C++ tests, and the daydream JS app / segmented-worker / tools / test layers, plus the Phantasm PCB-generation Python). Each component was audited by a dedicated reviewer agent working from the README architecture, and **every candidate finding was then independently re-validated by a separate skeptical agent** that re-read the cited code before the finding was admitted. Findings were held to the project's own fix-eligibility bar: a real, code-verified defect or concrete improvement with a minimal fix that does **not** regress performance or break hardware timing. Of the raw candidate set, **52 were confirmed and 12 were rejected** on validation; the rejected set is listed at the end.

**Repo of each file below is unambiguous by path/extension:** `.js` / `.html` / `.test.js` and `tools/…` live in **daydream**; `core/…`, `effects/…`, `hardware/…`, `targets/…`, `tests/*.cpp` and `tests/*.h`, `scripts/…`, `CMakeLists.txt` live in **Holosphere**.

---

## Executive summary

This is an exceptionally well-engineered project. Across performance-critical C++ and a non-trivial WebGL/WASM/Web-Worker JavaScript simulator, the review surfaced **zero critical and zero high-severity defects in the shipped engine, renderer, or simulator**. The single High finding is a version-sort bug in the *offline PCB-generation tooling* (no runtime path). The remaining 51 findings are 6 Medium and 45 Low — overwhelmingly latent edge cases, defensive-guard asymmetries, documentation/comment drift, and micro-optimizations, not live bugs.

The codebase's defining qualities are **rigorous memory discipline** (a heap-free arena model with debug generation-stamping and an always-on fail-fast trap doctrine), **provably-reasoned concurrency** (a three-index ISR double buffer on the device; a generation-fenced, serialized one-frame-deep Web-Worker pipeline in the simulator), and **documentation that is genuinely best-in-class** — nearly every non-obvious invariant is reasoned out at the point of code, and host/device parity hazards are proven rather than hand-waved. Testing is a standout: an always-on death harness that verifies the traps actually fire, analytic geometric oracles instead of golden-image diffs, C++/JS parity tests, and an event-driven multi-board simulator for the Phantasm sync protocol, all wired into CI budget gates.

The defects that do exist cluster in the **peripheral surfaces** — the JS UI (URL/GUI coercion asymmetries, a stuck-drag on one tool page), the measurement/test harnesses, and the PCB-gen scripts — precisely where the project's otherwise-uniform rigor thins out.

### Overall grade: **A-**

A mature, internally-consistent, and unusually well-documented system whose few blemishes are concentrated away from its core. It is held back from a straight A only by the accumulation of latent-guard asymmetries and the lower rigor of its tooling/JS periphery relative to the engine.

---

## Grades by dimension

| Dimension | Grade | One-line rationale |
|---|:---:|---|
| Correctness | **A-** | Core render/math/mesh essentially defect-free with analytic invariants; the 21 correctness findings are all latent/cosmetic edge cases (mostly JS UI + tooling), none a live engine bug. |
| Memory safety | **A** | Heap-free arena model, debug generation/rebind stamps, deleted temporary-binding, death-harness-verified traps; one cold-path face-index bounds gap. |
| Concurrency & thread safety | **A** | Three-index ISR double buffer with justified atomics orderings; simulator generation fence + serialized worker queue + layered watchdogs. One 1 Hz telemetry `__enable_irq()` nit. |
| Performance & efficiency | **A** | Compile-time `<W,H>` specialization, strict hot/cold split, split-trig LUTs, one-frame-deep worker pipeline; only two micro-inefficiencies found. |
| Architecture & elegance | **A** | Clean layering; three-domain compile-time filter pipeline; arena ping-pong with explicit `Arena&` params; pure-core / device-or-DOM-shell separations throughout. |
| API design & interface expressiveness | **A-** | Lifetime contracts encoded in the type system (`FunctionRef` vs `StoredFunctionRef`, deleted rvalue overloads), explicit coordinate-space ctors; minor coarse-bool returns and one residual lvalue-borrow gap. |
| Error handling & robustness | **A-** | Consistent cold-path fail-fast doctrine; residual gaps are a few missing non-finite guards, an unbounded `generate()` recursion, and a 0-byte recorder download. |
| Testing & verification | **A** | Death harness, analytic oracles, WASM parity, pov_sync board simulator, CI budget gates; a couple of vacuous-pass edges and a probabilistic stack estimator behind a hard gate. |
| Documentation | **A+** | Repeatedly best-in-class; a 2000-line reference README plus load-bearing inline rationale. Only 2 doc defects found across the entire tree (one misattached Doxygen block, one stale CMake comment). |
| Code style & readability | **A** | Uniform conventions, terse fact-focused comments, no history/rename noise, no dead code; enforced by `-Werror` and lint rules. |
| Portability & build | **A-** | `platform.h` cleanly abstracts Teensy/WASM/Desktop; CI compiles both firmware images and size-gates them. Dinged by a dead MSVC arch-guard token, a release-tuned stack gate reused for debug, and the KiCad-version-sort bug. |
| Security | **A-** | Largely N/A for an offline engine; the web surface is disciplined — `textContent` over `innerHTML`, untrusted JS ints clamped before enum/modulo use, importmap core keys un-clobberable. No injection surface found. |
| Maintainability & modularity | **A** | Single-source-of-truth X-macro rosters with bidirectional CI pins defeat whole classes of drift; residual risk is a handful of duplicated-logic / coercion-lockstep pairs. |

---

## Dimension rationale (detail)

**Correctness (A-).** The rasterizers, color math, mesh operators, and vector/quaternion code are pinned by analytic oracles (ring-latitude bands, Euler `V-E+F=2`, half-edge symmetry, tiling-once) and were found free of live defects. Every confirmed correctness finding is an edge case that is either currently unreachable (sub-3-gon Conway faces with the shipped Platonic seeds; NaN-poisoning of Mobius params that no current input path can produce; odd-width segment split against two even hardware resolutions) or cosmetic (label-snap tolerance). The genuine user-visible ones live in the JS UI and tooling: the paused WebGL context-restore blank-sphere, the spline tool's stuck off-canvas drag, and two URL/GUI coercion asymmetries.

**Memory safety (A).** The arena allocator uses overflow-safe subtractive bounds math, `ScratchScope` enforces LIFO teardown structurally, and `ArenaVector`/`ArenaSpan` carry *two* debug generation stamps (arena-reset and per-vector rebind) so both a reset and a realloc-without-bump fault in debug. There is no heap on either backend. The one real gap is `build_half_edge_mesh` reading face indices straight into `out.vertices[v]` guarded only by a stripped `assert` — a cold path that should carry an always-on `HS_CHECK` like its siblings.

**Concurrency (A).** On the device, the ISR reads a front buffer while the main loop renders a back buffer, mediated by three `std::atomic<int>` indices with the ISR never touching `cur_`; the Phantasm flywheel timebase and its release/acquire orderings are each individually justified and covered by an event-driven simulator. In the browser, the segmented pipeline's generation fence (`renderGen`/`inflightGen`) drops stale-sized frames while still settling the barrier so the pipeline can't wedge, and the worker message queue serializes a long WASM-init await against later resolution changes. The only smell is a debug-only telemetry poll that calls `__enable_irq()` unconditionally rather than save/restoring `PRIMASK`.

**Performance (A).** Templating the whole pipeline on `<W,H>` removes all generality overhead at the hardware resolution; hot paths are guard-free single loads while `HS_CHECK` is confined to cold seams; split trig LUTs cut geometry-LUT memory ~145x; the worker pipeline's wall-clock is `max`, not `sum`, of segment times. Only two inefficiencies survived validation — a redundant per-segment `screen_step` recompute and a per-frame GUI sync that runs even for effects with no animated params.

**Documentation (A+).** This is the project's most exceptional dimension. The review repeatedly found that *why*, not just *what*, is captured inline: the modular-multiply window proof, the DMAMEM/vague-linkage GCC pitfall, the H_OFFSET device/sim divergence, the `__FINITE_MATH_ONLY__` clamp contract, the arena polarity of every Conway op. Across the whole tree only two documentation defects were found — a Doxygen block attached to the adjacent function, and one stale "never exercised by CI" comment that CI now contradicts.

**Testing (A).** The death harness re-execs a shell-free child and requires the exact `SIGILL` trap *shape* (defeating `exit(128+SIGILL)` impostors); rosters are triple-pinned (CMake list / X-macro list / include block, bidirectionally) so a dropped test goes red, not silently green; dedicated TUs recompile the engine under the *actual* shipping flags (`-ffast-math -fno-finite-math-only`, the device H_OFFSET path). Weaknesses are minor: two tests can vacuously pass (a fake that doesn't force the copy-out it means to pin; an untested empty-URL-param edge), and the stack-depth estimator is a sentinel-density heuristic behind a hard CI gate.

**The periphery (why not a straight A).** The medium findings and a disproportionate share of the low ones sit in three places where rigor thins: the daydream JS UI (coercion lockstep between `state.js` and `gui.js` that already diverges on edge inputs; an accessibility roving-tabindex edge), the measurement harnesses (`stack_measure` hand-rolls the canonical global reset and will silently drift), and the Phantasm PCB-gen scripts (a lexical version sort that selects KiCad 9 over 10; a fab NOTE that always false-fires; a CPL export that blanks bottom-side coordinates). None touch a runtime render path, but collectively they are the reason the engine's A-grade rigor isn't uniform across the product.

---

## Notable strengths

- **Fail-fast as a verified contract, not a hope.** `HS_CHECK` survives `NDEBUG`, pulls in no stdio, formats into a fixed stack buffer (safe under a corrupt arena), and a death harness asserts each trap actually fires — with a probe that pins the trap's exact signal shape.
- **Lifetime safety encoded in types.** `StoredFunctionRef`/`ArenaSpan` `=delete` construction from a temporary so a dangling bind fails to *compile*; `FunctionRef` keeps rvalue binding for call-scoped borrows. The distinction is deliberate and documented.
- **Single-source-of-truth anti-drift.** `HS_RESOLUTIONS` and the X-macro test/effect rosters, pinned bidirectionally against CMake and the WASM smoke script, make "add an effect / resolution" a one-edit change and turn any divergence into a red build.
- **Host/device parity is proven.** Numeric provenance (double-precision folding to reproduce generated tables bit-for-bit), ODR-isolated recompiles of the south-pole path, and the DMAMEM specialization macro that dodges GCC's vague-linkage section-attribute drop all show parity treated as a first-class invariant.
- **Robust simulator concurrency.** A generation fence, a serialized worker queue, a composite pre-pass that validates every segment rect before blitting any, and three layered watchdogs each targeting a distinct non-throwing hang.

---

## Prioritized remediation list

All 52 confirmed findings, grouped by priority and numbered sequentially. Severities marked *"(severity revised down by validator)"* were confirmed real but re-graded during independent validation. There are **no Critical findings**.

### High priority

1. ✅ **find_kicad_cli picks KiCad 9 over 10 due to lexical version sort** — `hardware/phantasm/gen/fab.py:88-92` · _portability-build_
   Fix: Sort by a parsed numeric version key, e.g. key=lambda p: tuple(int(x) for x in re.findall(r'\\KiCad\\(\d+)\.(\d+)\\', p)[0]) with a safe fallback, or extract the version directory and compare as (major,minor) ints. Then take max() of that key rather than lexical [-1].


### Medium priority

2. ✅ **Final LCSC 'missing' NOTE always fires for every assembled part (ignores LCSC_BY_REF)** — `hardware/phantasm/gen/fab.py:230-233` · _correctness_
   Fix: Delete the redundant lines 230-233, or recompute using the override: missing = [r for r in assembled if not (LCSC_BY_REF.get(r) or comps[r]['lcsc'])] — matching the logic on line 194-195.

3. ✅ **CPL export uses --side front only; bottom-placed SMD get blank coordinates** — `hardware/phantasm/gen/fab.py:178-218` · _correctness_
   Fix: Export both sides (drop `--side front`, or run front+back) so bottom parts carry real coordinates, and pass through the exporter's Side column (already read as p.get('Side','top')). If the fab policy truly is top-only assembly, assert/skip bottom-side assembled parts loudly instead of emitting blank-coordinate CPL rows.

4. ✅ **WebGL context restore does not request a repaint, leaving a blank sphere while paused** — `driver.js:300-304` · _correctness_
   Fix: In onContextRestored, set this.needsRender = true (or call this.invalidate()) so the next animation-loop tick repaints and re-uploads instanceColor/instanceMatrix after the context comes back.

5. ✅ **Spline point-drag can get stuck when the mouse is released outside the canvas** — `tools/splines.html:413-415, onMouseUp/onMouseMove` · _correctness_
   Fix: Attach the drag-continuation listeners to document: register onMouseMove and onMouseUp on document (keep mousedown on the canvas). Alternatively, use Pointer Events with setPointerCapture on the canvas so release is always delivered.

6. ✅ **stack_measure hand-rolls reset_globals() instead of sharing it; silently drifts when the canonical reset grows** — `tests/stack_measure.cpp:55-56 (and 63-64)` · _maintainability_
   Fix: Include tests/test_fixture.h in stack_measure.cpp and replace the four inlined reset lines with hs_test::reset_globals(); (matching arena_measure.cpp / perf_bench.cpp). test_fixture.h is header-only and pulls no extra link deps beyond what stack_measure already links (memory.cpp / reaction_graph.cpp). This also deletes the dead global_timeline_t assignment.

7. ✅ **build_half_edge_mesh does not validate face indices against num_verts before indexing out.vertices[v]** — `core/mesh/mesh.h:288-299` · _memory-safety_
   Fix: In the per-index loop, add HS_CHECK(u < num_verts && v < num_verts, "half-edge mesh face index out of range"); before writing out.vertices[v].half_edge (or validate once when reading faces_arr[face_offset + i]).


### Low priority

8. ✅ **Motion holds a bare Fn path adapter over a borrowed path, but only the temporary-rejecting overload guards lifetime; a captured non-temporary that outlives the path still dangles** — `core/animation/motion.h:160, 169-172` · _api-design_
   Fix: No code change is warranted if the borrow contract is considered sufficient (consistent with the rest of the engine). If tightening is desired, document the lvalue-lifetime requirement on the surviving ctor as explicitly as MeshMorph's, or (debug-only) stash a generation token; do not add per-frame checks on the hot path.

9. ❌ **DMALEDController::setBrightness/setCorrection/setTemperature mutate HD107SFrame<N> statics, silently coupling any two controllers sharing N** — `hardware/dma_led.h:290-320` · _api-design_ _(severity revised down by validator)_
   Fix: Make the correction state per-instance (move tempR_/corrR_/brightness_ etc. out of the HD107SFrame statics into DMALEDController members passed to correct()), or add the same singleton-per-N trap the SPI driver has. The per-instance move costs nothing on the hot path (the factors are still read once per pixel) and removes the hidden aliasing.
   Rejected: the aliasing is unreachable in a shipped image — DMALEDController is one-per-image and a second controller's begin() already traps at TeensySPIDMA's single-instance guard (dma_led.h:76) before any second controller can exist to share the statics. The per-instance move would ripple the correction-factor storage through the FASTRUN correct() path and all 12 static call sites in test_hd107s_frame.h for a latent-only condition, disproportionate to a Low finding.

10. ✅ **SphericalHarmonics negative-lobe green scale truncates instead of rounding** — `effects/SphericalHarmonics.h:276-278` · _code-style_
   Fix: Add +0.5f before the cast: static_cast<uint16_t>(pos.color.g * NEG_LOBE_GREEN_SCALE + 0.5f). No perf impact (cold per-pixel shade already does float math here).

11. ✅ **wire() adds a junction dot at every L-bend corner (promised pruning never happens)** — `hardware/phantasm/gen/builder.py:141-146` · _code-style_
   Fix: Do not emit a junction at an L-bend corner — remove the self.junctions.append((x2,y1)) on line 145 (a corner where exactly two collinear-per-axis segments meet needs no junction), or if kept for safety, actually implement the promised pruning that drops junctions with degree < 3. Update or delete the misleading 'may be pruned later' comment.

12. ✅ **telemetry() health poll runs long hs::log() under __enable_irq() unconditionally, re-enabling IRQs even if the foreground path expected them masked** — `hardware/pov_segmented.h:316-350` · _concurrency_
   Fix: Save/restore PRIMASK instead of unconditional __enable_irq(): `uint32_t m = __get_PRIMASK(); __disable_irq(); tm = sync_.telemetry(); __set_PRIMASK(m);`. Zero added cost on the hot path (this is a 1 Hz debug branch).

13. ✅ **MobiusFlow floors num_rings/num_lines but doesn't guard non-finite (NaN/Inf) live scalars, which propagate into a/d and poison MobiusParams** — `core/animation/params.h:399-420 (MobiusFlow::step)` · _correctness_ _(severity revised down by validator)_
   Fix: Replace the sign clamps with finite-and-range clamps, e.g. `float rings = std::isfinite((float)num_rings) ? std::max(0.0f, (float)num_rings) : 0.0f;` and similarly `lines = std::isfinite((float)num_lines) ? std::max(1.0f, (float)num_lines) : 1.0f;`. Two `isfinite` checks per frame on a cold animation step (not a per-pixel loop) are free relative to the expf/cos/sin already there.

14. ✅ **MobiusWarp/MobiusWarpCircular do not guard a non-finite `scale`/live `scale_ref_`, unlike MobiusFlow and Driver** — `core/animation/params.h:468-476 (MobiusWarp::step), 509-516 (MobiusWarpCircular::step)` · _correctness_
   Fix: Add `HS_CHECK(std::isfinite(scale))` in both constructors (cold path, mirrors Driver), and in MobiusWarp::step keep the last good scale on a non-finite live read (`if (scale_ref_) { float s2=*scale_ref_; if (std::isfinite(s2)) s=s2; }`), mirroring Driver::step.

15. ✅ **Sprite constructor's fade re-scaling can strand a caller-requested fade-in at 0 frames with no diagnostic** — `core/animation/sprites.h:48-54` · _correctness_
   Fix: Optional: when the scaled fade_in floors to 0 while the requested value was >0 and duration>=2, bias one frame back to fade_in (e.g. round instead of floor) so both fades survive; keep the fade_total>duration branch otherwise unchanged. Purely a UX nicety, not a correctness fix.

16. ✅ **PetalFlow can spawn more rings than the MAX_RINGS pool holds at high Speed×Density, silently dropping the excess with no visible density guarantee** — `effects/PetalFlow.h:185-217 (check_spawn / spawn_ring_at_pos)` · _correctness_
   Fix: Add a one-line comment on spawn_ring_at_pos()'s no-free-slot return path stating the pool can saturate at high Speed×Density and the dropped spawn is an accepted transient (as RingShower documents), so the MAX_RINGS sizing rationale is discoverable. No code change; purely the missing rationale that the rest of the file otherwise supplies everywhere.

17. ❌ **Ricker wavelet is hard-truncated at the fast-reject band edge, producing a displacement discontinuity** — `core/engine/transformers.h:338-353` · _correctness_ _(severity revised down by validator)_
   Fix: Either widen the reject band so it clips where the wavelet has decayed to negligible (e.g. |t|~=3, d +/- 1.5*thickness, ~e^-4.5 tail) or smoothstep the wavelet to zero over the last part of the band. The band widening is the zero-perf-cost option: only prepare_thresholds() (cold, per animation step) changes; the per-pixel path is unchanged.
   Rejected: prepare_thresholds() already sets the band to phase ± thickness, which maps to t = ±4 at the reject edge (t = (d-phase)·4/thickness), where the wavelet is (1-16)·e⁻⁸ ≈ -0.5% of peak — the truncation step is already sub-perceptible. The suggested |t|≈3 widening is mis-scaled and would *narrow* the band, *enlarging* the step to ricker(3) ≈ -9% of peak; any genuine widening also enlarges the per-pixel active region on a perf-sensitive ripple path.

18. ✅ **Flywheel::position() overflow guard is a stripped assert only; a genuinely unfolded coast wraps int32 and mis-columns silently under NDEBUG** — `hardware/pov_sync.h:615-630` · _correctness_
   Fix: Replace the assert() with HS_CHECK (always-on cold trap) so an unfolded-coast invariant break traps on the device instead of silently returning a wrong column under NDEBUG. It is one branch per wake-up, not per pixel, so it fits the cold-path trap policy.

19. ✅ **advanceFrameClock only consumes one frame interval per rAF, so a recovered backlog is capped to real-time regardless of catchup budget** — `driver.js:450-461` · _correctness_ _(severity revised down by validator)_
   Fix: Either loop the interval-consume inside stepSimulation up to the catchup budget to actually replay backlog, or simplify to a plain 'step when >= interval' and update the comment to say missed frames are dropped (not caught up). Pick one and make code+comment agree.

20. ✅ **URLSync seeds 0 for an empty numeric URL param, diverging from gui.js's parseFloat coercion** — `state.js:157-160` · _correctness_
   Fix: Use parseFloat(raw) instead of Number(raw) in the URLSync numeric branch so both deserializers agree, e.g. `const num = parseFloat(raw); if (!Number.isFinite(num)) continue;`.

21. ✅ **URLSync does not type-coerce boolean tracked keys, seeding a raw string into a boolean-typed key** — `state.js:156-163` · _correctness_
   Fix: Mirror gui.js's boolean parsing in the URLSync constructor for `typeof current === 'boolean'` keys (accept true/1/yes/on and false/0/no/off, else skip), so the sync-layer coercion matches the GUI-layer coercion for all primitive types.

22. ✅ **sidebar onKeyDown Enter/Space selects the focused option without moving the roving tabindex to it** — `sidebar.js:289-294` · _correctness_
   Fix: In the Enter/Space branch, call `this.setRovingTabbable(focused)` before `this.onSelect(...)`, matching the arrow-key path so the tab stop always tracks the last-acted option.

23. ❌ **prettify's symbolic snapping uses an absolute 1e-5 tolerance that is too coarse near small constants and too fine near large ones** — `label_format.js:21-38` · _correctness_ _(severity revised down by validator)_
   Fix: If tighter consistency is wanted, compare against the same precision the fallback uses (snap only when r.toFixed(3) equals the constant's toFixed(3)), or use a relative tolerance `Math.abs(r - k) <= 1e-5 * Math.max(1, Math.abs(k))`. Otherwise document that 1e-5 is a deliberate visual-only threshold.
   Rejected: every snapped constant is O(1) magnitude, so the fixed 1e-5 absolute tolerance is already appropriate and the JSDoc documents it; a relative tolerance would alter label output for no visible benefit.

24. ✅ **renderParallel() leaves stale timings/renderUs/arenas for segments that never reported this frame** — `segment_controller.js:629-644, 810-849` · _correctness_
   Fix: In renderParallel(), when opening a fresh dispatch, also reset the per-segment stat arrays that updateStats() reads (e.g. this.timings.fill(0), this.renderUs.fill(0), this.arenas.fill(null)) so a fenced-out or non-reporting segment shows '-'/0 rather than a stale prior-generation value. Alternatively gate updateStats() cells on results[s] being non-null (which already tracks the current generation).

25. ❌ **computeSegmentRange silently drops the trailing column for odd canvas widths** — `segment_layout.js:65-67` · _correctness_
   Fix: No code change needed for current resolutions. If defensiveness is wanted, computeSegmentRange could throw on odd w (the GUI never exposes one), converting a silent dark column into a fail-fast — consistent with the module's other Number.isInteger/positive guards. Otherwise leave as documented.
   Rejected: the floor(w/2) split intentionally matches the firmware's w/2 partition and is already documented; all shipped resolutions are even, and throwing would diverge from hardware behavior on the one path that must stay bit-identical.

26. ✅ **kis emits degenerate triangles for <3-side faces instead of skipping them like every sibling operator** — `core/mesh/conway.h:428-449` · _correctness_ _(severity revised down by validator)_
   Fix: Tighten the guard to `HS_CHECK(count >= 3, "kis: degenerate face (< 3 sides)")`, matching relax's identical guard at line 864; costs nothing on the hot path (this is HS_COLD).

27. ✅ **truncate's ambo short-circuit tolerance (1e-4) silently redirects a narrow band of near-0.5 t values** — `core/mesh/conway.h:571-573` · _correctness_
   Fix: Either tighten the comparison to an exact `t == 0.5f` (safe here because callers pass literal constants) or update the truncate doc to state the +/-TOLERANCE snap band explicitly.

28. ✅ **expand/chamfer/snub emit an edge face even when an adjacent primary face was degenerate, silently producing a malformed quad/triangle** — `core/mesh/conway.h:730-738` · _correctness_
   Fix: This is a design tradeoff for the intentional self-intersecting cases and may be acceptable; at minimum document that edge/orbit faces are NOT gated on well_formed so a future maintainer does not assume the whole operator degrades cleanly on sub-3-gons. No perf cost to a comment.

29. ✅ **apply_pole_containment ray-cast can divide through inv_edge_j == 0 for a near-horizontal straddling edge, yielding a degenerate intercept** — `core/render/sdf.h:2325-2327, 2345-2347` · _correctness_
   Fix: Skip the crossing contribution when inv_edge_j[i] == 0.0f (guard the intercept branch with a nonzero check), matching build_canonical_distance_lut/plane_dist_exact which rely on a real ey; a horizontal edge contributes no crossing anyway. No hot-path cost (this is cold ctor code).

30. ✅ **Doxygen block for pair_half_edges is misattached to sort_edge_records** — `core/mesh/mesh.h:130-165` · _documentation_
   Fix: Move the doc block (lines 130-139) down to immediately precede the pair_half_edges template at line 148, and give sort_edge_records its own one-line brief describing that it only sorts records in place by (min_v,max_v).

31. ✅ **CMakeLists debug-stack rationale contradicts CI: the wasm-debug module IS smoke-driven** — `CMakeLists.txt:62-63` · _documentation_
   Fix: Update the comment to reflect that CI does runtime-smoke the debug build (ci.yml 'Runtime-smoke WASM debug'), and keep the generous 64 KB sizing justified by -O0 frame inflation rather than by 'never exercised'.

32. ✅ **Ripple constructor does not validate `speed` finiteness; a NaN speed silently freezes the wave phase** — `core/animation/params.h:651-660 (Ripple ctor), 675 (params.phase += speed)` · _error-handling_
   Fix: Add `HS_CHECK(std::isfinite(speed), "Ripple speed must be finite");` alongside the existing duration check in the constructor.

33. ✅ **generate() recursion depth has no upper bound / overflow guard** — `core/engine/generators.h:51-58 (generate_depth / ++depth)` · _error-handling_
   Fix: After ++depth add HS_CHECK(depth <= kMaxGenerateDepth, "generate: recursion too deep") with a generous compile-time ceiling (e.g. 16). Zero hot-path cost (cold generation-time path).

34. **read_id() reports a strap fault via HS_CHECK trap but the two same-ID boards it cannot detect drive bus contention silently** — `hardware/pov_segmented.h:379-396` · _error-handling_
   Fix: No code change on the hot path; expand the read_id() doc to state that a duplicated-ID master electrically shorts the sync bus driver, so the field procedure must strap unique IDs, or move PIN_MASTER_EN to gate an open-drain/tri-stateable sync driver so a second master cannot source-fight the bus. If open-drain is adopted, the contention becomes benign and detectable.

35. **Streaming recorder can silently produce a zero-byte download when close() succeeds but no data was written** — `recorder.js:399-424` · _error-handling_
   Fix: Guard the download calls on chunk presence, e.g. in finish's `if (!writable)` branch and the plain fallback: `if (chunks.length) this.download(...)` (the catch branch at line 423 already does this), and/or log a warning when a session ends with no data.

36. **getNormalizedX yields NaN when the color strip has zero width** — `tools/palettes.html:641-650 (getNormalizedX)` · _error-handling_
   Fix: Guard the divisor: `return rect.width > 0 ? Math.max(0, Math.min(1, x / rect.width)) : 0;`

37. **update_hankin is_flat path can trap on a zero-length base vertex without the degenerate guard the non-flat path has** — `core/mesh/hankin.h:279-281` · _error-handling_ _(severity revised down by validator)_
   Fix: Use the same fallback shape as the non-flat branch, e.g. `compiled.dynamic_vertices[i] = normalized_or(p_corner, p_corner);` or reuse the existing fallback vector, keeping the branch cost identical.

38. **PeriodicTimer public `period` clamp is duplicated between ctor and set_period instead of shared** — `core/animation/timers.h:81 (ctor), 98-99 (set_period)` · _maintainability_
   Fix: Have the constructor delegate the clamp to set_period-style logic, or extract a `static int clamp_period(int p){ return p < 1 ? 1 : p; }` used by both. No runtime cost.

39. **GenerativePalette::colors[] is a write-only member; could be a local in update_stops** — `core/color/color.h:1425 (decl), 1228-1265 (update_stops)` · _maintainability_
   Fix: Make colors a local std::array<CPixel,kMaxStops> inside update_stops() instead of a member; the switch fills it and the loop consumes it in the same scope. No perf change (same fill/read), and it shrinks GenerativePalette by kMaxStops*3 bytes.

40. **HopfFibration lacks a scratch-budget static_assert that its sibling arena-configuring effects carry** — `effects/HopfFibration.h:116-127, 257-289 (render_trails)` · _maintainability_
   Fix: Either extend the retune-guard comment to state that TRAIL_LEN is bounded for scratch by the 16KB default scratch_a (not by kFootprintBytes), or add a small static_assert that TRAIL_LEN * sizeof(Fragment) <= DEFAULT_SCRATCH_A_SIZE so a TRAIL_LEN retune is genuinely fail-fast like the sibling effects. No perf impact (compile-time only).

41. **syncGUI/export length-skew guard uses paramNames vs values but export uses params vs values — an off-by-source mismatch** — `daydream.js:116-124` · _maintainability_
   Fix: Have both syncGUI and export derive the expected length from one authority (e.g. the engine's live definition count or a single cached count on activeEffect), rather than one from paramNames and the other from the params closure.

42. **size sort is not tie-broken by name, so equal-size effects list in arbitrary (insertion) order** — `sidebar_logic.js:24-27` · _maintainability_
   Fix: Add a name tiebreaker: `if (key === 'size') { const d = a.size - b.size; return (d || a.name.localeCompare(b.name, 'en')) * mul; }` (or apply the name compare unsigned so ties always sort A→Z regardless of direction).

43. **Arena-metrics setTimeout loop keeps rescheduling forever after an engine trap** — `tools/solids.html:751-766 (updateArenaMetrics)` · _maintainability_
   Fix: Return without rescheduling when MeshOpsWasm is null (the engine is halted and metrics can never update again), e.g. `if (!MeshOpsWasm) return;` before the setTimeout, and re-arm from update() if the module is ever reloaded.

44. **vector_to_pixel unit-length assert uses a bare 1e-3f literal instead of the named EPS_UNIT_VEC_SQ tolerance** — `core/math/geometry.h:374` · _maintainability_
   Fix: Replace the literal with the named constant: assert(std::fabs(dot(v, v) - 1.0f) < math::EPS_UNIT_VEC_SQ); (debug-only assert, zero device cost, aligns the slack with the other unit-vector checks).

45. **World::Trails leaves head_ member effectively unused after init, but flush() swap-remove relies on tail_ moving in lockstep** — `core/render/filter.h:742-770, 814-839` · _maintainability_ _(severity revised down by validator)_
   Fix: Add a one-line comment at the flush() cull loop noting that the live element being overwritten is the logical last (index count_-1 == (tail_-1)%Capacity) so only tail_ retreats and head_ must stay put. Documentation only; no behavior change.

46. **Per-frame syncGUI does a full Map lookup + resolveParamSync over every param even for non-animated effects** — `daydream.js:101-140` · _performance_
   Fix: Store hasAnimated (or the list of animated controllers) on activeEffect and early-return from syncGUI when there are no animated params, skipping the getParamValues() call entirely.

47. ✅ **rasterize() recomputes screen_step for the first sub-step it already computed as first_step** — `core/render/plot.h:626-627, 661-662` · _performance_
   Fix: Seed the loop with the already-computed first_step for the first iteration instead of recomputing: e.g. use `float step = (steps_cache.empty()) ? first_step : screen_step<W,H>(smp.pos, smp.tan, base_step);` or restructure so first_step is pushed before the loop. Keeps behavior bit-identical (same inputs -> same output) with one fewer sqrt/divide per segment.

48. ✅ **x86 clamp backend uses GNU always_inline attribute but the arch guard also admits MSVC (_M_X64/_M_IX86)** — `core/engine/platform.h:1305-1306 arch guard; 1352 clamp` · _portability-build_
   Fix: Drop `_M_X64`/`_M_IX86` from the arch guard (the engine builds only with GCC/Clang/emscripten), or gate the attribute behind a portable HS_ALWAYS_INLINE macro if MSVC is genuinely intended. Minimal: remove the two MSVC tokens so the guard matches actual targets.

49. ✅ **FakeEngine.getParamValues returns a plain array, so the worker's Array.from() copy-out is never exercised** — `tests/segment_worker.test.js:54, 248-261` · _testing_
   Fix: Have FakeEngine.getParamValues() return a stateful object that is NOT a plain array snapshot (e.g. a Uint16Array-like view whose contents mutate after the call, or an object that throws when structured-cloned), then assert the posted paramValues is a detached plain array equal to the value at call time. This forces the Array.from() copy to be load-bearing.

50. ✅ **URLSync numeric coercion uses Number() but the empty-string/whitespace coercion edge is untested** — `tests/state.test.js:98-110` · _testing_
   Fix: Add a case asserting the intended behavior for `?count=` (empty) — either that it is rejected and keeps the default (if that is desired, which would require switching the coercion or adding a `raw === '' ` guard in state.js) or that it deliberately coerces to 0. Pinning it documents the Number()-vs-parseFloat choice.

51. ❌ **Smoke-test stack creep gate is a fixed 2048 bytes but runs against both 8 KB release and 64 KB -O0 debug builds** — `scripts/wasm_smoke.mjs:26-27,109-120` · _testing_
   Fix: Scale the creep budget by build type — e.g. accept an env override (STACK_HWM_CEILING_BYTES) that ci.yml sets higher for the debug run, or derive the ceiling from stack.capacity so the -O0 build is measured against a proportionate budget rather than the release-tuned absolute.
   Rejected: the 2048-byte gate is a deliberate creep tripwire (documented as "not a bound"), independent of stack capacity by design; scaling it to the -O0 build's larger frames would only loosen the tripwire on the build most prone to creep, so no code change adds value.

52. ✅ **stack-measure sentinel-density heuristic can under-report peak depth if a workload frame happens to write the 0xA5 byte densely** — `tests/stack_measure.cpp:79-91` · _testing_
   Fix: Reduce false-sentinel probability by painting a 4-byte rotating pattern (or a per-address value) instead of a single repeated byte, so a workload's incidental 0xA5 runs cannot masquerade as unpainted; alternatively widen the margin by requiring several consecutive 'written' windows before accepting a peak, or add an assertion that the region immediately below top is written (catches a scan that stopped at frame 0). Keep the -Os build so codegen still tracks the size config.


---

## Appendix: findings investigated and rejected on validation

These candidate findings were **dismissed** by the independent validation pass as non-issues — already handled, intended design, a misread of the code, or a fix that would regress performance. They are recorded so the ledger is complete.

- **[engine-memory] apply_if_changed rejects enum parameters despite being the canonical live-apply idiom** — The static_assert at core/engine/util.h:139-140 does restrict `apply_if_changed` to arithmetic types. However, searching the entire codebase reveals no Params struct containing an enum-typed member that would naturally be passed to this function, and no hand-rolled if-changed workarounds for enum pa...
- **[math] Vector slerp antipodal window is narrower than fast_acos error, leaving a thin non-monotone band** — The code at lines 1237-1258 of core/math/3dmath.h confirms the structure described: `fast_acos` (peak error ~1.3e-4 rad, documented at line 1183) is used to get theta, and the antipodal gate is `theta > PI_F - math::TOLERANCE` where TOLERANCE is 1e-4 (line 49). The thin non-monotone band argument is...
- **[render-raster] scan_region full-row detection misses coverage assembled from multiple abutting spans, so a legitimately full-circle row still walks the sort/coalesce path** — Lines 181-184 of scan.h contain an explicit code comment that reads: "A single span covering the full circle (len >= W) paints every column; detect it up front and skip the seam-split/sort/coalesce path. Coverage assembled from multiple abutting spans is not caught here — it falls to the slow path, ...
- **[render-raster] Canvas buffer_free() spin-wait calls global ::micros() while the rest of the engine standardizes on hs::micros()** — The cited spin-wait at canvas.h lines 621-623 calls unqualified `micros()`. On host/test builds, platform.h defines a global `inline unsigned long micros() { return hs::micros(); }` at line 1077 (after the namespace block closes at line 1065). This global forwarder is not just an alias for legacy co...
- **[render-pipeline] Feedback::flush degrades warp field to a single coarse row when the segment band is shorter than one downsample cell** — I traced the arithmetic end-to-end: 1. The "collapsed vertical gradient" scenario requires the render band (y_lo to y_hi) to span fewer than ds rows AND land entirely within the last coarse cell (cy_lo == cy_hi == hh-1). On the actual Phantasm hardware with N=4, S=288, PPS=72, ds=4 (default): the re...
- **[effects-shader] Flyby drift live-control not in markAnimated set but registered as a param — pause semantics inconsistency is intentional but the comment mislabels the coupling** — The code at lines 50-58 of Flyby.h is exactly as described: `registerParam("Drift", &drift, 0.0f, 2.0f)` registers Drift, but it is deliberately omitted from the `markAnimated` loop. The comment at lines 52-54 already explains this precisely: "Drift is a standalone live control, not preset-driven, s...
- **[hardware] HD107SFrame::load() is a device-compiled parity path that never runs on hardware, carrying dead per-frame cost surface** — The finding claims load() "consumes flash and must be kept bit-identical to packPixel()/correct() forever" because it is "not host-gated." This misreads how C++ template instantiation works. HD107SFrame is a class template (template <int N>). Template member functions are only instantiated and emitt...
- **[js-core-a] instanceColor buffer left detached (array=null) after a lost-context repaint attempt in segmented mode** — The described race does not exist in practice, and the "transient" windows cited in the finding are not real: 1. setupDots() (line 622) sets instanceColor.array = null on the OLD mesh, then immediately builds a fresh mesh. After setupDots() returns, this.dotMesh refers to the NEW mesh whose instance...
- **[js-core-a] EngineHost.refresh() dereferences this.engine with no null guard, can throw during teardown/pre-load** — The finding assumes refresh() can be called before host.engine is assigned or after teardown nulls it. Neither path exists in the actual code. Guard against pre-load: The animation loop in daydream.js (lines 743-750) gates on `if (host.adapter)` before invoking `daydream.render(host.adapter)`. The a...
- **[js-segmented] rebuildBoundaries() scans the entire results array instead of the configured segment count** — The finding's premise requires `this.results` to have more entries than `this.count`, but the code makes this impossible. `create()` (line 247-249) sets both `this.count = numSegments` and `this.results = new Array(numSegments).fill(null)` in the same call, so the two are always equal for any live p...
- **[js-segmented] Worker render handler does not validate that getPixels() returned a live (non-detached) view** — The finding rests on two faulty premises: 1. "subarray() on a detached source silently yields an empty subarray and set() is a no-op, zero-filling the segment." This is wrong. In all modern JS engines (V8, SpiderMonkey, JavaScriptCore), calling .subarray() on a detached TypedArray throws a TypeError...
- **[tests-cpp] param_marshal check_one editable-param fallback can pick an unchanged value on a symmetric range** — The finding misreads the control flow. The retry branch at line 95 (`newv = lo + 0.25*(hi-lo)`) only executes when the initial pick `newv = lo + 0.5*(hi-lo)` equals `views[target].value`. The claim "value == 0.25-point alone defeats the retry" is incorrect: if `views[target].value` is at the quarter...


---

_Review conducted by a fan-out of per-component auditor agents with per-finding independent re-validation. Grades are assessed against professional embedded-graphics and web-application engineering standards._
