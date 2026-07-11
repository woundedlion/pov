# Holosphere / daydream — Code Quality Review

**Date:** 2026-07-11
**Scope:** The Holosphere C++ rendering engine + firmware and the daydream web
simulator, reviewed as one product across both repositories. Out of scope by
directive: `effects_legacy`, `targets/Holosphere/Holosphere.ino`,
`core/math/rotate.h`, and all vendored third-party code (`core/vendor/*`,
`daydream/three.js/*`, `daydream/node_modules/*`).

**Method:** 23 component reviewers (one per subsystem) read every in-scope file
against the project README's documented design contracts, gathered candidate
findings, and validated each against the cited source. An independent adversarial
verifier then re-opened the cited code for every finding and confirmed,
downgraded, or **rejected** it; 8 candidate findings were rejected outright and
are not carried below. Only findings that survived verification appear here. Every
finding meets the code-review-fix eligibility bar: a real defect, a concrete and
minimal fix, and no performance regression on a hot path.

---

## Overall Grade: **A**

Holosphere is a mature, exceptionally engineered real-time graphics system. The
review surfaced **39 verified findings across ~100 files, none critical or
high severity** — 9 medium (latent correctness/robustness)
and 30 low (documentation, missing tests, cosmetic). The defect
density and severity profile are what you would expect from a codebase that has
been through sustained, disciplined review. The engine's structural decisions
(compile-time resolution, arena allocation, fail-fast invariants, 16-bit linear
color, a variadic filter pipeline) are not just sound but elegantly executed and
self-defending against drift. The residual weaknesses are narrow: a recurring
class of documentation comments that have fallen out of step with the code, a
handful of reachable-but-untested paths, and small API-consistency warts.

### Grades by Dimension

| Dimension | Grade |
|---|---|
| Correctness & Bug-freedom | **A** |
| Architecture & Design Elegance | **A** |
| Interface Expressiveness / API Design | **A-** |
| Readability & Code Style | **A** |
| Documentation | **A-** |
| Test Coverage & Quality | **A-** |
| Performance & Efficiency | **A** |
| Memory & Resource Management | **A** |
| Error Handling & Robustness | **A** |
| Portability & Platform Abstraction | **A** |
| Maintainability & Modularity | **A** |
| Build System & Tooling | **A-** |


---

## Rationale by Dimension

**Correctness & Bug-freedom — A**  
Across dense arithmetic (arena bounds, OKLab/OKLCH transforms, quaternion Shepperd branch, flywheel/beacon codec, Conway operators) no functional bug was confirmed; the only correctness findings are latent degenerate-input edges (near-zero centroid in `dual`, a codegen re-round) that cannot arise from the registered solids or normal parameter ranges.

**Architecture & Design Elegance — A**  
A single templated `<W,H>` engine feeds three targets from one source; clean layering (pixel core → color-space math → palette composition; generate → transform → rasterize → filter); X-macro single-sources (`HS_EFFECT_LIST`, `HS_RESOLUTIONS`) and compile-time anchors make whole classes of drift a build error. Six components earned A+ here.

**Interface Expressiveness / API Design — A-**  
Intent is encoded in the type system throughout — `FunctionRef` vs `StoredFunctionRef` (borrow vs store), `ArenaVector`/`ArenaSpan` (own vs borrow), `spawn`/`spawn_pinned` (retained-handle contract), explicit `Arena&` at every call site. Docked a step for recurring small warts: non-`const` pure getters, a mutable-reference `operator[]`, and one angle parameter typed `float_t` against an otherwise-uniform convention.

**Readability & Code Style — A**  
Uniformly terse, fact-focused comments that state rationale for the non-obvious and nothing for the obvious; consistent naming and structure across 100+ files. Only two components dipped to A-.

**Documentation — A-**  
Doxygen coverage is extensive and usually precise, and the README is a genuine architecture text. Documentation errors are nonetheless the single most common finding category — a handful of doc comments that contradict the verified behavior of the code they describe (buffer written, range emitted, fade direction, stale filename).

**Test Coverage & Quality — A-**  
The strongest testing culture is unusual for the domain: a native suite with a forked SIGILL death-harness that asserts fail-fast traps actually fire, reference-triple oracles, full-operand sweeps, and C++/JS parity tests. The weakest dimension only because several load-bearing paths (arena-evacuation `compact`, `HoleRef`, clone-across-compaction, a few coercion branches) are reachable but unpinned.

**Performance & Efficiency — A**  
Compile-time resolution erases generality overhead; hot loops are branch-free on device; transformers iterate only active slots; amplitude/band early-rejects skip per-pixel wavelet math; baked LUTs give O(1) lookups. One documented-rule violation (an always-on check on the per-fragment path) is the sole blemish.

**Memory & Resource Management — A**  
The core competency: a fixed 330 KiB partitioned arena, zero heap on the render path, RAII `ScratchScope`/`Persist` evacuation, subtractive wrap-proof bounds, and debug-generation stamping that faults stale borrows deterministically.

**Error Handling & Robustness — A**  
Fail-fast `HS_CHECK` placed exactly per the stated policy — cold allocation/config/registration seams — with transient conditions (dropped frames, DMA overrun) correctly given soft handling. Findings are missing belt-and-suspenders validations (a clip upper bound, a beacon-fits invariant), not mishandled errors.

**Portability & Platform Abstraction — A**  
`platform.h` cleanly isolates Arduino/WASM/desktop; host mocks reproduce device quirks bit-for-bit (including Cortex-M7 divide-by-zero trap semantics); ARM DSP intrinsics are gated with verified portable fallbacks.

**Maintainability & Modularity — A**  
Small focused headers with explicit contracts; single-source-of-truth macros and compile-time anchors mean adding an effect, resolution, or preset is a bounded, self-defending edit.

**Build System & Tooling — A-**  
CMake/PlatformIO/Emscripten targets, generated-artifact provenance (LUTs, reaction graph) with self-checks, a WASM runtime smoke harness, and roster-drift CI gates. Docked for tooling-script robustness gaps (an ineffective screenshot gate, swallowed subprocess stderr, a temp file written into the tree).

---

## Notable Strengths

- **Memory model.** A fixed 330 KiB arena partitioned into persistent + two
  scratch pools, zero heap on the render path, RAII scoping, wrap-proof
  subtractive bounds, and debug generation-stamping that turns a stale borrow into
  a deterministic fault. This is the safest embedded-allocation design in its class.
- **Fail-fast discipline, placed correctly.** `HS_CHECK` traps at cold seams and
  survives `NDEBUG`; hot paths use stripped asserts backed by a cold trap; a forked
  death-harness proves the traps actually fire (`SIGILL`). Policy and enforcement match.
- **Drift-proofing by construction.** X-macro single sources (`HS_EFFECT_LIST`,
  `HS_RESOLUTIONS`), `static_assert` anchors on preset field order, and
  registry-vs-list oracles convert "someone forgot to update the other place" from a
  latent bug into a compile error.
- **Sim/device parity investment.** Host mocks reproduce FastLED and even the
  Cortex-M7 divide-by-zero trap bit-for-bit; C++/JS color and spline parity are
  pinned by tests. The simulator is a faithful twin, not an approximation.
- **Perceptual color done properly.** 16-bit linear compositing with OKLCH
  shortest-arc hue interpolation and chroma-preserving gamut mapping gated off the
  hot path — a level of color correctness rare even in professional tools.
- **Testing culture.** Reference-triple oracles, full-operand-range sweeps, a
  fast-math re-run of the NaN/overflow contracts, and adversarial identity-vs-real
  transform oracles. The tests are engineered, not decorative.

---

## Prioritized Findings

All 39 verified findings, numbered sequentially and grouped by priority. Each
is a candidate for the code-review-fix workflow.

### Critical Priority
_None. No defect in scope can crash the device, corrupt the arena, produce data loss, or trigger undefined behavior in normal operation._

### High Priority
_None. No defect produces wrong rendered output, breaks a public contract, or constitutes a live hazard under supported configurations._

### Medium Priority
_Latent correctness / robustness / convention defects: no observable impact under current usage, but each is a real gap worth hardening._

1. ✅ **dual's degenerate-centroid fallback can emit an off-sphere vertex** — `core/mesh/conway.h:379` (Holosphere)  
   dual() is the only operator that does not call normalize(out_mesh) before returning; its per-vertex output relies on normalized_or(c, first_v) to land on the unit sphere. When the face centroid c is (near-)zero (a centrally-symmetric face), normalized_or returns the fallback first_v verbatim, and first_v is the raw face vertex, never re-normalized. kis/ambo use the same raw-fallback pattern but are saved by their trailing normalize(); dual is not. On the shipped registry dual is always fed unit meshes so first_v is unit and the branch never fires, making this latent — but the code documents no unit-input precondition and the sibling operators all self-correct.  
   *Fix:* Either pass a normalized fallback (normalized_or(c, first_v.normalized()) guarded, or first project) or add a trailing normalize(out_mesh) to dual for parity with the other operators.
2. ✅ **Dynamo drops sub-step motion-blur age on strand-body line segments** — `effects/Dynamo.h:291` (Holosphere)  
   draw_nodes(canvas, age) is called once per sub-step within a multi-step frame with age = i/steps so intermediate positions age out sooner (World::Trails::plot computes ttl = lifetime - round(age), a shorter TTL for older sub-steps). The head node (i==0) is plotted with this age, but the connecting line segments (i>0) are drawn via Plot::Line::draw from f_from/f_to Fragments whose .age defaults to 0, so Plot::Line forwards age 0 and every body segment from every sub-step gets the full TTL. At the default Speed=2 (and increasingly up to the +/-10 slider max) the strand body therefore piles up as full-strength ghost copies persisting the whole trail length instead of a fading intra-frame smear, defeating the age parameter for the bulk of the visible geometry.  
   *Fix:* Set f_from.age = age; f_to.age = age; before the Plot::Line::draw call so the line segments inherit the same sub-step age as the head node.
3. ✅ **fp_bbox lists fp_poly as handled but never reads its pts geometry** — `hardware/phantasm/gen/pcb.py:201` (Holosphere)  
   In fp_bbox the graphic-outline branch matches 'fp_poly' but the body only extracts the 'start'/'end'/'center' keys (lines 202-205); an fp_poly stores its outline in a 'pts' list, which is never read, so any polygon-defined graphic (courtyard/silk, used by some stock KiCad footprints) contributes nothing to the bounding box. Today the pads dominate the extent and the 1.2 mm pack gap absorbs the small courtyard margin, so parts do not overlap — but the helper silently under-sizes any footprint whose extent is defined by a poly, which is a latent packing hazard if a poly-bounded part is added.  
   *Fix:* Add an fp_poly case that iterates its ('pts' ...) xy children and folds each point into xs/ys (or drop fp_poly from the tuple if it is truly unsupported).
4. ✅ **initScene default resize never refreshes devicePixelRatio** — `tools/shared.js:146` (daydream)  
   renderer.setPixelRatio(window.devicePixelRatio) is called once at construction (line 100) and defaultResize only calls renderer.setSize(w, h), which reuses the stored pixel ratio. Dragging the tool window from a standard display to a Hi-DPI/Retina monitor (which fires a resize event with a changed window.devicePixelRatio) leaves the canvas rendering at the stale ratio, so the scene renders blurry/under-sampled until reload.  
   *Fix:* Add renderer.setPixelRatio(window.devicePixelRatio) inside defaultResize before renderer.setSize(w, h).
5. ✅ **Screenshot resolution gate warns but exits 0, contradicting its stated intent** — `scripts/capture_screenshots.mjs:252` (Holosphere)  
   When resolveResolutions() fails and returns [], the whole gallery is captured at the app-default (unpinned) resolution. The summary block's own comment says this 'must not look like a clean success', yet unlike the sibling gates (wrongRes at 275, blanks at 288, failures at 296) it only prints a warning and never sets process.exitCode. If no effect happens to go blank or throw, the run exits 0 and CI goes green on an unverified-resolution gallery, defeating the freshness intent stated in the comment.  
   *Fix:* Set process.exitCode = 1 in the RESOLUTIONS.length === 0 block (or soften the comment to state the degradation is intentionally non-fatal).
6. ✅ **cpp_format re-round can still silently collapse a nonzero value to "0.0f"** — `tools/cpp_format.js:37` (daydream)  
   formatFloatCpp exists specifically to stop a nonzero-but-tiny value from rounding to "0.0f" (see its own header). The re-round path caps precision at min(100, ...), but never re-checks the result. For a magnitude needing more than 100 fractional digits (roughly |n| < 5e-101, a normal double), toFixed(100) yields all zeros, the trim/pad produces "0.0", and the function returns "0.0f" — silently discarding the coefficient, the exact failure it was written to prevent, with no error raised. It is unreachable for the tool's real value ranges (palette/spline/lissajous coefficients ~1e-6..1e2) but is a genuine contract gap given the function throws on other precision-loss cases (non-finite, exponential).  
   *Fix:* After the re-round, if parseFloat(s) === 0 still holds, throw an Error like the non-finite/exponential guards instead of returning "0.0f".
7. ✅ **Config::valid() does not verify the worst-case beacon frame fits the [W/4, W/2) emission window** — `hardware/pov_sync.h:233` (Holosphere)  
   maybe_schedule_beacon() only schedules a beacon when position x satisfies `x >= W/4` and `x + beacon_span_cols() < W/2`. Since the earliest schedulable x is exactly W/4, a beacon can ever be emitted only if `W/4 + beacon_span_cols() < W/2`, i.e. `beacon_span_cols() < W/4`. Config::valid() checks many inter-constant relationships (join_grid divides 64, beacon_period < 32, epoch train fits the refractory window, the beacon int32 span bound) but never checks this one. With the shipped defaults it holds (span 55 cols < W/4 = 72), so the gap is latent. But a future tuning of beacon_pitch_cols / gap_timeout_cols / W that pushed beacon_span_cols() >= W/4 would make the guard trip for every x, so the master would silently emit no beacons at all — downstream boards would never learn their identity and would stay dark forever, exactly the kind of self-inconsistent configuration valid() exists to trap at boot (run_show() HS_CHECKs cfg.valid()).  
   *Fix:* Add `&& beacon_span_cols() < W / 4` to the Config::valid() conjunction so a beacon frame that cannot fit its emission window is rejected at boot instead of silently suppressing all beacons.
8. ✅ **set_clip validates lower bounds and ordering but not y1 <= h** — `core/render/canvas.h:146` (Holosphere)  
   set_clip's HS_CHECK traps y0>=0, y0<=y1, x0>=0, x0<=x1, and its comment defers the x1<=w check to a downstream LUT-domain trap in Scan::Shader::draw. The y-upper bound (y1 <= clip_.h) is neither trapped here nor deferred with a comment; a driver passing an oversized y1 silently clamps at render_y_end() instead of trapping the sizing bug at this config seam. There is no memory-safety consequence (render_y_end clamps to h, so no OOB), but it is the one segment bound at this fail-fast seam that neither traps nor documents where it is checked, breaking the symmetry the surrounding code establishes.  
   *Fix:* Add `y1 <= clip_.h` to the set_clip HS_CHECK (mirroring the x1<=w deferral note), so an oversized segment traps at the seam rather than silently clamping.
9. ✅ **Always-on HS_CHECK in MeshPaletteBank::operator[] runs on the per-fragment mesh-shading hot path** — `core/color/palettes.h:161` (Holosphere)  
   Both MeshPaletteBank::operator[] overloads run HS_CHECK(i>=0 && i<N) on every call; that call is on the per-pixel mesh path (shade_mesh_topology at render/shading.h:108,139 does palette_bank[palette_idx[slot]].get(t) per rasterized fragment). The index is structurally always in [0,N) (a shuffle permutation), so the trap can never fire, yet its predicted-not-taken branch is emitted per fragment and cannot be hoisted (slot varies with frag.v2). This contradicts the README rule that HS_CHECK guards cold paths only and is never in the per-pixel loop; peer hot accessors (BakedPalette::get) use a stripped assert.  
   *Fix:* Replace the two HS_CHECK calls in MeshPaletteBank::operator[] with a stripped assert(), relying on the always-on HS_CHECK at the cold bake/shuffle seam.

### Low Priority
_Documentation inaccuracies, missing tests for real behavior, dead code, and cosmetic/style nits._

10. ✅ **Message-less HS_CHECK on KDTree::nearest k>MAX_K precondition** — `core/mesh/spatial.h:111` (Holosphere)  
   The public-API precondition trap `HS_CHECK(k <= static_cast<size_t>(MAX_K));` carries no message, unlike every other HS_CHECK in these three files (e.g. lines 82, 194, 268 in spatial/mesh). On a caller misusing nearest() with k>5 the trap fires with no context string, making the fail-fast crash harder to diagnose on device where there is no stdio.  
   *Fix:* Add a message, e.g. `HS_CHECK(k <= static_cast<size_t>(MAX_K), "KDTree::nearest k exceeds MAX_K");`.
11. ✅ **ChaoticStrings lacks the persistent-footprint static_assert its sibling trail effects carry** — `effects/ChaoticStrings.h:73` (Holosphere)  
   The Node (Orientation<16> plus OrientationTrail<Orientation<16>, TRAIL_LENGTH=115>) is allocated from persistent_arena, but ChaoticStrings only static_asserts the scratch budget (SCRATCH_A_BYTES) and has no compile-time guard on the persistent allocation. The sibling trail effects with the identical TRAIL_LENGTH pattern (Comets, HopfFibration) both carry a FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET assert explicitly annotated as guarding a TRAIL_LENGTH retune. It fits comfortably today (persistent is GLOBAL_ARENA_SIZE - 200 KiB), but a TRAIL_LENGTH or ORIENTATION_SUBSTEPS bump could overflow the partition at init with only a runtime arena trap, not a compile error, on the device.  
   *Fix:* Add a static_assert bounding the Node's persistent size (sizeof(Node)) against the effect's persistent partition, mirroring Comets' FOOTPRINT_BYTES guard.
12. ✅ **kicad-cli subprocess failures swallow stderr, yielding opaque gate errors** — `hardware/phantasm/gen/pcb.py:35` (Holosphere)  
   export_netlist() (and identically check.py:37) run kicad-cli with check=True, capture_output=True. On a nonzero exit the raised CalledProcessError does not print the captured stderr, so a missing/failed kicad-cli surfaces only as 'Command ... returned non-zero exit status 1' with no reason — poor for a CI validation gate and for operators reproducing locally.  
   *Fix:* Wrap the run in try/except subprocess.CalledProcessError and print e.stderr (or run with check=False and report result.stderr on a nonzero returncode) before re-raising/exiting.
13. ✅ **Exported temp netlist leaks on failure and is written into the source tree** — `hardware/phantasm/gen/check.py:13` (Holosphere)  
   check.py writes _check.net (and pcb.py writes _pcb.net) next to the schematic in the source tree, then os.remove()s it after parsing. Because the kicad-cli run uses check=True, any subprocess failure raises before os.remove, leaving the temp file behind in the repo working tree.  
   *Fix:* Write the netlist to tempfile.NamedTemporaryFile / a tempdir, or remove it in a try/finally so a subprocess failure cannot strand it in the source tree.
14. ✅ **prettify passes non-finite input through to a raw 'NaN'/'Infinity' label** — `label_format.js:39` (daydream)  
   prettify has no guard for non-finite input: every symbolic Math.abs comparison is false for NaN, and r.toFixed(3) yields 'NaN' (or 'Infinity' for +/-Infinity), so a degenerate geometry value would render the literal text 'NaN' on an on-sphere axis label rather than a sane placeholder. It is unlikely to be hit in practice but is a trivially closable gap and is the one input the label_format test suite does not cover.  
   *Fix:* Add `if (!Number.isFinite(r)) return "0";` (or return an empty/placeholder string) at the top of prettify and add a test for it.
15. ✅ **Orientation historical-index accessors orient(v,i)/unorient(v,i)/at(i) are untested** — `core/math/geometry.h:412` (Holosphere)  
   test_geometry.h exercises Orientation::set/push/get/get(i)/collapse/orient/unorient/upsample thoroughly, but the by-index rotation accessors orient(const Vector&, int) (line 412), unorient(const Vector&, int) (line 435), and the mutable at(int) (line 501) — each carrying an HS_CHECK bounds guard on the frame index — have no test. A regression in their index handling (e.g. an off-by-one into the history array) would not be caught by the existing suite.  
   *Fix:* Add a geometry test that pushes a few known rotations then asserts orient(v,i)/unorient(v,i) match rotate(v, get(i))/its conjugate and that at(i) aliases the stored quaternion.
16. ✅ **MeshCarousel::compact and compact_keep_front arena-evacuation paths are untested** — `core/animation/mesh.h:595` (Holosphere)  
   test_animation.h exercises MeshCarousel::schedule_segue and every Segue policy, but never compact() or compact_keep_front(). Both are non-trivial and correctness-sensitive: compact() relies on Persist destruction order to restore both slots after persistent_arena.reset(), and compact_keep_front() relies on running after_reset(persistent_arena) while the front slot is still evacuated (i.e. before the Persist guard restores it on scope exit). A regression in either the Persist ordering or the after_reset timing would silently corrupt the persistent arena / mesh slots with no test catching it.  
   *Fix:* Add a test that populates both MeshCarousel slots, calls compact()/compact_keep_front() (asserting the front slot survives, the back slot is dropped, and after_reset runs before the front restore), mirroring the mesh-building helpers already used by the MeshMorph test.
17. ✅ **Determinism per-frame checksum truncates 16-bit pixel channels to 8 bits** — `tests/test_effects.h:309` (Holosphere)  
   render_capture's frame_fold lambda is declared fold_byte(uint8_t byte) but is called with p.r/p.g/p.b, which are uint16_t (Pixel = Pixel16, uint16_t r,g,b — see core/color/color.h:55). Each channel is implicitly narrowed to its low 8 bits before being folded into the FNV-1a checksum (test_scan.h confirms these channels hold values like 40000). The fold exists specifically to catch mid-run nondeterminism that reconverges before the final frame — but any per-frame divergence confined to the high byte of a channel (e.g. r=0x9C40 vs 0x0040, differing by 0x9C00) folds identically and is invisible. -Wconversion is not part of -Wall -Wextra, so the truncation compiles clean under -Werror. The final-frame buffer comparison (a[i].r != b[i].r) still uses full 16-bit values, so persistent divergence is caught; only the reconverging mid-run case that the fold was added to detect is weakened.  
   *Fix:* Fold both bytes of each channel, e.g. replace each call with fold_byte(p.r & 0xFF); fold_byte(p.r >> 8); (and likewise for g/b), or change fold_byte to take uint16_t and hash both bytes.
18. ✅ **AppState.update re-entrant stale-tuple skip is untested** — `state.js:75` (daydream)  
   AppState.update batches writes then guards each queued notification with `if (this.state[key] !== value) return;` so a subscriber that re-enters set()/update() mid-batch and changes a still-queued key does not get a superseded, out-of-order notification (lines 62-79). This is the subtlest correctness invariant in the module — a regression here would silently emit stale (key,value,old) tuples — yet state.test.js only covers non-re-entrant batching ('AppState.update batches and only fires for changed keys'). No test drives a subscriber that mutates a sibling batched key.  
   *Fix:* Add a state.test.js case where a subscriber to the first batched key calls set() on a second still-queued key, asserting that key is notified exactly once with the re-entrant (live) value and the stale batch tuple is skipped.
19. ✅ **URLSync boolean tracked-key coercion branch is untested and currently unexercised** — `state.js:161` (daydream)  
   URLSync's constructor has a full boolean-coercion branch (accepting true/1/yes/on and false/0/no/off, line 161-165) for tracked keys whose seeded default is a boolean. The app currently registers only string tracked keys ('effect','resolution'), so this branch is dead in production, and state.test.js exercises only the numeric and string coercion paths. A future boolean tracked key would ship on untested logic (e.g. the round-trip through setTrackedParam's String(val) serialization).  
   *Fix:* Either add a state.test.js case seeding a boolean default and asserting URL coercion + round-trip serialization, or remove the boolean branch until a boolean tracked key exists.
20. ✅ **Horizontal Y-band boundary overlay path is untested** — `c:/work/daydream/segment_controller.js:774` (daydream)  
   rebuildBoundaries() populates boundaryYs only from segments whose y0>0, and composite() draws a full-width cyan row for each — real behavior in the 6/8-segment (multi-Y-band) modes. Every composite/boundary test in segment_controller.test.js uses full-height segments (y0=0, y1=2), so boundaryYs is never non-empty and the horizontal-seam drawing loop (lines 774-778) has zero coverage; a regression there (e.g. an off-by-one in rowStart or a wrong y-guard) would ship silently.  
   *Fix:* Add a composite() test with two stacked Y-band segments (e.g. y0=0..1 and y0=1..2 within one arm) asserting a horizontal cyan seam at the band boundary row and untouched interiors.
21. ✅ **setParam numeric non-finite drop-key path only covered via the explicit null marker** — `tests/state.test.js:177` (daydream)  
   URLSync.setParam (state.js lines 219-223) routes a non-finite number through roundUrlNumber -> null and drops the key rather than serialize a misleading 0. The suite tests key-drop only via an explicit setParam('speed', null) marker (line 177); it never passes a NaN/Infinity number, so the number-branch's non-finite-to-null conversion (and roundUrlNumber's Number.isFinite guard, which is exported and otherwise only indirectly hit) is untested.  
   *Fix:* Add a case asserting sync.setParam('speed', NaN) followed by flush() drops 'speed' from the URL, mirroring the existing null-marker test.
22. ✅ **FoldModifier negative-input folding branch is untested** — `core/color/composition.h:148` (Holosphere)  
   FoldModifier::modify contains a deliberately-documented correctness branch — 'fmodf keeps the dividend's sign, so reduce into [0,2) first; negative scaled would otherwise fold above 1' (the `if (m < 0.0f) m += 2.0f;` guard). This branch only fires for a negative `scaled` (negative coordinate or a negative phase driver). test_palette_modifiers (tests/test_color.h ~line 1266) only exercises fold at t = 0.0/0.25/0.5 with folds=2 and no phase, so the entire negative-reduction path — the reason the guard exists — has zero coverage. A regression that dropped or mis-signed the guard would produce folded values > 1 (out-of-range palette coordinates) with a green suite. The same suite also never exercises PinchModifier's documented negative-t path, QuantizeModifier's s<1 floor / dynamic driver, or RippleModifier/BreatheModifier at non-zero coordinate/phase.  
   *Fix:* Extend test_palette_modifiers to assert FoldModifier with a negative coordinate (and a negative phase driver) stays within [0,1] and matches the triangle-wave reference, plus a non-zero Ripple/Breathe and a Pinch negative-t case.
23. ✅ **Reference-origin Hole (HoleRef) variant is untested** — `core/render/filter.h:658` (Holosphere)  
   test_world_hole_masks_cap exercises only the by-value Filter::World::Hole<W>; the HoleRef alias (Hole with std::reference_wrapper<const Vector> origin storage) has no direct test. HoleRef exists specifically so the hole center can track a live vector, and its only variant-specific line is the `const Vector &o = origin;` reference-wrapper conversion in plot(); a regression that broke live-origin tracking would not be caught by the current suite.  
   *Fix:* Add a test that constructs a HoleRef over a mutable Vector, plots through it, mutates the referenced origin, and asserts the mask follows the new center.
24. ✅ **ArenaVector::data() @return says nullptr when "empty" but a bound, size-0 vector returns non-null** — `core/engine/memory.h:627` (Holosphere)  
   The @return briefs on both data() overloads state "Mutable/Const pointer to the first element, or nullptr if unbound/empty." A vector that is bound with capacity>0 but has size()==0 (the common cleared/just-bound state) returns a non-null data_ pointer, so "empty ⇒ nullptr" is false for that case — only the unbound or moved-from state yields nullptr. The @details block correctly explains this, but the @return line (the text shown in generated API docs and IDE hover) contradicts it and could mislead a caller into using `if (!v.data())` as an emptiness test. The const overload at line 641 carries the identical wording.  
   *Fix:* Change both @return briefs from "nullptr if unbound/empty" to "nullptr if unbound or moved-from" to match the actual behavior and the @details.
25. ✅ **vector_to_theta/vector_to_pixel doc overstates x range as [0,W] (can equal W)** — `core/math/geometry.h:321` (Holosphere)  
   The doc for vector_to_theta<W> states the return is "in `[0, W]` (can round up to `W`; floor before indexing)", and vector_to_pixel's doc (lines 336-339) echoes "x from wrap() can round up to W". The value comes from util.h wrap<float,int>, whose contract and guard `return (r >= mm) ? R{0} : r;` (util.h line 48) fold the boundary back to 0, so the result is strictly in [0, W) and can never equal W. The advice to floor before indexing stays safe, but the stated [0,W]/"can round up to W" range is factually wrong and could mislead a reader into thinking a column index of W is reachable.  
   *Fix:* Change the doc range to `[0, W)` and drop the "can round up to W" clause, since wrap() strictly excludes W.
26. ✅ **emit_shrunk_face pushes face_counts before the HE_NONE start guard** — `core/mesh/conway.h:202` (Holosphere)  
   emit_shrunk_face() pushes narrow_face_count(count) into out_mesh.face_counts when count>=3, then immediately returns if start==HE_NONE without emitting any face indices, which would desync face_counts from faces (a face claiming `count` sides with zero indices). Today this is unreachable because count is derived by face_centroid walking the same face.half_edge as `start`, so count>=3 implies start!=HE_NONE — but the ordering couples correctness to that invariant being preserved at every call site rather than being locally robust.  
   *Fix:* Move the `if (he_idx == HE_NONE) return;` above the well_formed face_counts push (or fold the push into the loop's first iteration) so a side count is only recorded when its indices will actually follow.
27. ✅ **Redundant dead conditional in Union::get_horizontal_intervals** — `core/render/sdf.h:891` (Holosphere)  
   The guard `if (!has_a && !has_b) return false;` is fully subsumed by the immediately-following `if (!has_a || !has_b) return false;` — every input that satisfies the first (both children fell back) also satisfies the second, returning the identical value, so the first `if` is unreachable as a distinct branch and can never change behavior. SmoothUnion::get_horizontal_intervals (the parallel combinator) correctly has only the single `||` check, confirming the first Union check is vestigial. Harmless but dead, and it makes the two otherwise-mirrored combinators gratuitously inconsistent.  
   *Fix:* Delete lines 891-892 (`if (!has_a && !has_b) return false;`), leaving the single `if (!has_a || !has_b) return false;` that already covers the both-fell-back case, matching SmoothUnion.
28. ✅ **WarpedVolume Warp concept documents an apply() signature its only implementation doesn't match** — `core/render/sdf.h:3480` (Holosphere)  
   The WarpedVolume doc block specifies the Warp concept as `Vector apply(const Vector &p, const Ctx &ctx) const;` (Ctx passed by const-reference), but the sole provided Warp, Warp::Twist, declares `Vector apply(const Vector &p, Ctx /*s*/) const` (by value). Because Ctx is a `float` typedef both compile at the `warp.apply(p, ctx)` call site, so there is no functional bug, but the documented concept signature is inaccurate and would mislead an author writing a new Warp (e.g. one with a heavier Ctx) into copying the wrong convention.  
   *Fix:* Update the concept comment to `Vector apply(const Vector &p, Ctx ctx) const;` (or make Twist::apply take `const Ctx &`) so the documented contract and the reference implementation agree.
29. ✅ **Rotation::animate takes float_t angle while every other angle parameter is float** — `core/animation/motion.h:444` (Holosphere)  
   The static one-shot helper declares `float_t angle` (from <cmath>), whereas every other angle in the file — the Rotation ctor (`float angle`), Motion/path_frame parameters, make_rotation calls — uses plain float. float_t is an implementation-defined evaluation type that widens to double under FLT_EVAL_METHOD==2, so on such a target the function's signature/overload set silently differs from the rest of the API. It is harmless today (callers pass a float that widens then narrows back for make_rotation) but is a real, needless type inconsistency in a signature.  
   *Fix:* Change the parameter type from `float_t` to `float` to match every other angle parameter in the component.
30. ✅ **Rotation default constructor documented as 'identity' but stepping it traps** — `core/animation/motion.h:334` (Holosphere)  
   The default constructor's doc says it 'Creates an inactive/identity rotation.' It leaves `orientation == nullptr`, and Rotation::step begins with `HS_CHECK(orientation != nullptr)`, so stepping a default-constructed Rotation traps rather than acting as an identity (no-op) rotation. 'inactive' is accurate but 'identity' misleads a reader into thinking it is a safe steppable no-op; only has_orientation()-guarded reassignment is valid before use.  
   *Fix:* Reword the brief to 'Creates an inactive, unbound rotation (must be reassigned before stepping; stepping an unbound Rotation traps).'
31. ✅ **Thrusters comment calls an alpha-blended fade 'additive'** — `effects/Thrusters.h:252` (Holosphere)  
   draw_thruster's fragment shader comment says it applies opacity 'for the additive thruster glow', but the rings go through Plot::Ring -> the AntiAlias pipeline -> the base Pipeline<W,H> terminal, which composites with straight-alpha (src*a + dst*(1-a)) in linear light per the README, not additive blending. Dimming the color by opacity and setting alpha = opacity*params.alpha produces a normal alpha-blended fade, so the 'additive' wording misdescribes the actual compositing.  
   *Fix:* Reword the comment to describe an alpha-blended (opacity-driven) fade rather than an additive glow.
32. ✅ **MindSplatter hardcodes full_frame instead of deriving it from the filter pipeline like every sibling** — `effects/MindSplatter.h:31` (Holosphere)  
   The constructor passes only `{.strobe = true}`, leaving full_frame at its default (false), whereas every sibling effect with a screen-filter pipeline (DreamBalls, ChaoticStrings, GnomonicStars, HopfFibration, SplineFlow, PetalFlow, MobiusGrid) sets `.full_frame = decltype(filters)::any_crosses_segments`. It is currently correct because `Filter::Screen::AntiAlias::crosses_segments` is false, so any_crosses_segments evaluates to false either way. But the value is hand-pinned rather than derived: if a stateful filter (Feedback, Screen::Trails, both with crosses_segments==has_history==true) is later added to MindSplatter's pipeline, full_frame stays false and Phantasm segmented rendering silently breaks at segment seams, with no compile-time signal.  
   *Fix:* Change the ctor options to `{.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}` to match the sibling effects and auto-track the pipeline.
33. ✅ **Dead `idx` accumulator declared at function scope in drawFrame fast path** — `targets/wasm/wasm.cpp:449` (Holosphere)  
   `int idx = 0;` is declared before the fast/slow branch but is only ever used inside the `else` (overrides_get_pixel) slow path; when the memcpy fast path is taken it is initialized and never read. It reads as if both branches share the accumulator, which is misleading.  
   *Fix:* Move the `int idx = 0;` declaration inside the `else` block where it is actually used.
34. ✅ **Pure getters getBufferLength()/getRenderUs() are not const** — `targets/wasm/wasm.cpp:480` (Holosphere)  
   getRenderUs() (line 480) and getBufferLength() (line 523) only read state yet are non-const, while the analogous strobeColumns() (line 495) is correctly const; the inconsistency is a minor const-correctness gap (harmless to embind, which does not require const).  
   *Fix:* Mark getRenderUs() and getBufferLength() const to match strobeColumns() and the read-only nature of both methods.
35. ✅ **setActive assigns activeName before confirming the target button exists** — `sidebar.js:163` (daydream)  
   setActive removes the active class/aria-selected from the previous button and unconditionally sets this.activeName = name before checking whether a button for `name` exists in the current roster. If called with a name absent from the current list, it deselects the current option and leaves the listbox with no highlighted/selected option while activeName points at a nonexistent button. It self-heals on the next valid setActive and current call sites always pass an in-list name (resolveActiveEffect guarantees it), but the method's own contract is violated for an off-list argument.  
   *Fix:* Move `this.activeName = name;` into the `if (newBtn)` block (or early-return when the button is missing) so an off-list name doesn't strip the current selection.
36. **"quadrant" terminology is inaccurate for 6/8-segment layouts** — `c:/work/daydream/segment_controller.js:686` (daydream)  
   Multiple doc comments describe segments as "quadrants" and composite() as the "quadrant model" (lines 8, 686, and the FrameResult typedef), but computeSegmentRange explicitly supports total in {2,6,8}, where segments are halves/sixths/eighths, not quadrants. The wording understates the supported layouts and can mislead a maintainer into assuming a fixed 4-way split.  
   *Fix:* Replace "quadrant"/"quadrant model" with "segment rectangle"/"segment-rectangle model" in the segment_controller.js comments.
37. ✅ **Stale source-file reference in effect_registry.h doc comment** — `core/engine/effect_registry.h:7` (Holosphere)  
   The header comment states the self-registering factory eliminates "the hand-maintained list in wasm_bridge.cpp", but no wasm_bridge.cpp exists anywhere in the tree — the WASM entry/bindings live in targets/wasm/wasm.cpp (verified: a repo-wide search for `wasm_bridge` returns only this comment, and targets/wasm/ contains wasm.cpp). The dangling reference misdirects anyone tracing where the registry replaced the old manual roster.  
   *Fix:* Change `wasm_bridge.cpp` in the comment to `targets/wasm/wasm.cpp`.
38. ✅ **Feedback::flush doc says it writes the front buffer; it writes the back buffer** — `core/render/filter.h:1277` (Holosphere)  
   The @param on Feedback::flush reads 'Target canvas (reads cv.prev, writes the front buffer).' The method reads the front/displayed buffer via cv.prev(x,y) but writes the current draw target via cv(x,y), which is the BACK buffer (bufs_[cur_]). README section 5 and the filter's own class doc correctly describe sampling the front buffer and compositing into the back buffer, so this @param contradicts the rest of the documentation and could mislead a reader into thinking flush mutates the displayed frame.  
   *Fix:* Change 'writes the front buffer' to 'writes the back (current-draw) buffer'.
39. ✅ **RingShower FADE_IN_FRAMES doc says opacity fades "from 0" but the first drawn frame is 0.25** — `effects/RingShower.h:80` (Holosphere)  
   The FADE_IN_FRAMES member comment states "Frames spent fading in from 0 to full opacity", but opacity_at() returns ease_linear((age+1)/FADE_IN_FRAMES), so the first visible frame (age 0) yields ease_linear(1/4) = 0.25, and the fade runs 0.25 -> 0.5 -> 0.75 -> 1.0. The ring never renders at opacity 0 (this is intentional and mirrors radius_at()'s age+1 "one step in" convention, which is correctly documented at line 108). Only this one comment overstates the starting value.  
   *Fix:* Reword line 80 to say the fade-in starts one linear step in (0.25 at FADE_IN_FRAMES=4) rather than at 0, matching the age+1 convention already documented on radius_at().

---

## Verifier-Rejected Candidates (not counted above)

For transparency, 8 candidate findings were raised by reviewers and **rejected**
by the adversarial verifier as non-issues or intentional design — including a
claimed `dual` 2-gon self-pairing bug (load-bearing for Conway operators), a
`bake_all` "always allocates" nit (matches its documented contract), and several
"missing test" claims for behavior already covered indirectly. They are recorded
in the run journal and deliberately excluded from the fix list.
