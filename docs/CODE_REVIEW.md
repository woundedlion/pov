# Holosphere / daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`,
`targets/`, `tests/`) and the daydream web simulator (`*.js`).
**Out of scope (per review charter):** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/rotate.h`, and vendored third-party code (`FastNoiseLite*`, `inplace_function.h`).

**Method.** A multi-agent audit: 17 component reviewers read every in-scope file against the
README architecture, and each candidate finding was then re-checked by a *fresh, independent*
validator that read the cited code before the finding was admitted. 59 candidate findings were
raised; **39 were confirmed**, **20 were rejected** as false positives, already-handled, or
intentional documented design (see [Rejected findings](#rejected-findings)). Nothing rose above
**Medium** severity — there are **no correctness-critical, memory-safety, or security defects** in
the shipped paths.

---

## Overall grade: **A**

This is a mature, unusually disciplined codebase. The engineering is coherent from the physics of
the spinning LED arm all the way up to a browser simulator that shares the *identical* C++ via
WebAssembly. The design decisions are deliberate and documented (16-bit linear + OKLCH color,
compile-time `<W,H>` specialization, a single partitioned arena with explicit `Arena&` threading,
fail-fast `HS_CHECK`, a 2-buffer ISR double-buffer, and a formally-specified single-wire multi-Teensy
flywheel sync protocol). Test coverage is broad and includes a forked death-harness that proves the
fail-fast traps actually fire. The confirmed defects are, almost without exception, latent footguns,
documentation nits, minor UX/consistency gaps, and cold-path cleanups — the residue you find only
after the real bugs are already gone.

### Quality dimensions

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | **A** | Math (quaternion/Möbius/OKLab/Legendre/KDTree/flywheel timebase) is careful and singularity-guarded; validated against `std::`/brute-force oracles. All confirmed correctness items are edge-case/robustness, not mainline logic errors. |
| **Memory safety** | **A** | Arena bump-allocation with wrap-proof overflow math, always-on `HS_CHECK` traps at cold seams, debug-only generation stamps for use-after-free, exact-fit pool sizing, no per-frame heap. No OOB/UB found. |
| **Performance** | **A** | LUT-backed trig, branchless wrap, fast-math approximations pinned to error bounds, baked palettes, SSAA vertex/fragment split, relaxed atomics + explicit cache flush on the DMA path. Tuned to the Teensy budget throughout. |
| **Architecture & elegance** | **A** | Clean layering (util → math → geometry → color → raster → filter → effects), a compile-time variadic filter pipeline with `static_assert`-enforced domain ordering, pure/host-testable cores split from Arduino shells, and a self-registering effect factory. |
| **Interface expressiveness / API design** | **A-** | `FunctionRef`/`StoredFunctionRef`/`PipelineRef` erasure, owned-vs-borrowed `ArenaVector`/`ArenaSpan`, and named palette modifiers are excellent. Docked for a scattering of documented-but-unenforced contracts and mildly under-typed registers (`DistanceResult` v0–v3, a few `int`/`uint16_t` slips). |
| **Maintainability** | **A-** | Exceptional inline rationale; X-macro single-sources for the effect roster and resolutions. Docked for a little dead API surface, three parallel deep-copy/half-edge idioms that can drift, and one stale budget figure. |
| **Testing** | **A-** | Triple-pinned module roster, a best-in-class death harness, golden hashes, energy-conservation oracles, per-effect arena budgets, a stack high-water CI gate, and Teensy+WASM builds in CI. Docked for a missing aggregate-arena CI gate, an unpinned death-case count, and `driver.js`/a few effects lacking white-box coverage. |
| **Documentation** | **A+** | Among the best-documented embedded code you will encounter: every non-obvious invariant, platform divergence, and design trade-off is explained at its site, plus a thorough README and a formal frame-sync datasheet. |
| **Portability** | **A** | A tight `platform.h` abstracts Arduino/WASM/Desktop; deterministic seeding gives bit-identical cross-target output. One latent signed/unsigned index idiom is the only blemish. |
| **Build & tooling** | **A** | `just` task runner, CMake WASM build that installs artifacts into the sibling repo, CI that compiles both Teensy images (size/layout-gated) and runs a runtime WASM smoke over every effect. |

### Component grades (per reviewer)

| Component | Corr | Perf | Mem | Arch | API | Maint | Test | Docs |
|---|---|---|---|---|---|---|---|---|
| core-math | A | A | A | A | A- | A | A | A+ |
| core-color | A- | A | A- | A | A- | A | A | A |
| core-pipeline | A | A- | A | A | A- | A- | A- | A+ |
| core-raster | A- | A | A- | A | B+ | A- | B | A |
| core-anim | A | A | A | A- | B+ | B+ | A | A+ |
| core-memory | A | A | A | A | A- | B+ | A | A- |
| core-mesh-a | A- | A | A- | A- | A- | B+ | B+ | A |
| core-mesh-b | A | A | A | A | A- | A- | A | A |
| effects-1 | A- | A- | A | A | A- | A- | B+ | A |
| effects-2 | A- | A- | A | A | B+ | A- | B+ | A |
| effects-3 | A | A- | A | A | A- | A- | A | A |
| hardware | A | A | A- | A | A- | A | A | A |
| targets (wasm/Phantasm) | A | A | A | A | A- | A- | B+ | A+ |
| tests | A | A | A | A | A | A- | A- | A+ |
| daydream-core | A- | A | A | A | A- | A- | B+ | A |
| daydream-workers | A- | A | A | A | A- | A | A | A |
| daydream-ui | A- | A | A | A | A- | A- | A | A |

---

## Prioritized fixes

Every confirmed finding, numbered sequentially. Each is eligible for the `code-review-fix` flow
(validated real, actionable, in-scope, not an intentional-design item). Mark each `✅` (fixed) or
`❌` (rejected) with a one-line reason as it is worked.

### Critical

*None.*

### High

*None.*

### Medium

1. ✅ **Death-case roster has no anti-drift count pin.** `run_death_tests()` iterates `all_cases()` and asserts each child trapped but never pins the case count, so a dropped entry shrinks coverage silently while the module stays green — the one roster here lacking the floor every other list has. *Fix:* add `HS_EXPECT_GE(n, <N>)` after computing `n` (mirror the `kMinAssertions` floors), or pin the count from CMake. (`tests/test_death.h:797-844,1112-1118`)

2. ✅ **No CI gate on total global-arena high-water vs `DEVICE_GLOBAL_ARENA_SIZE`.** `stack_measure` is a gating CTest but its sibling `arena_measure` (which computes per-effect arena totals vs the 330 KiB device budget) is built as an informational executable with no `add_test()`, so an effect that grows past the device arena only surfaces as an on-hardware OOM trap — the exact asymmetry the stack gate removed. *Fix:* make `arena_measure` a gating CTest that exits non-zero when a per-effect total exceeds `DEVICE_GLOBAL_ARENA_SIZE` (with headroom). (`tests/CMakeLists.txt:195-205`, `tests/arena_measure.cpp:59-64`)

### Low

3. ✅ **`slerp(Vector)` antipodal branch tests a threshold finer than `fast_acos`'s own error.** The branch uses `theta = fast_acos(d) > PI_F - TOLERANCE` with `TOLERANCE=1e-4`, narrower than `fast_acos`'s ~1.3e-4 rad envelope; the sibling quaternion slerp deliberately uses exact `acosf`. Output stays finite/unit-length, so impact is negligible. *Fix:* gate the antipodal decision on the raw dot `d` (e.g. `d < -1+k`) so it does not depend on the approximation. (`core/3dmath.h:1225-1232`)

4. ✅ **`shortest_distance` / `fwd_distance` have divergent range-reduction contracts.** Same signature shape, but `shortest_distance` fully reduces arbitrary input while `fwd_distance` silently requires `|b-a| < m` (guarded only by a debug-stripped assert). No live misuse today (`fwd_distance` has no callers), but a footgun. *Fix:* either make `fwd_distance` robust with the same double-`fmod` reduction, or rename it `fwd_distance_one_period`. (`core/util.h:114-137`)

5. ✅ **`frac_to_q16` wraps to 0 on `frac > 1` instead of saturating.** Unlike its siblings (`blend_alpha`, `float_to_pixel16`) it casts without clamping, so `~1.0001 → 65536 → 0` (the opposite endpoint). All current callers keep `frac` safe. *Fix:* clamp inside the helper to match its siblings (hot-path-neutral). (`core/color.h:218-220`)

6. ❌ **Base pipeline sink recomputes the cylindrical x-clip (two modulos) per plotted pixel.** The terminal `Pipeline<W,H>::plot` calls `contains_x()` (two `% w`) per emitted sub-sample, amplified by multi-tap filters, while the scan/sdf hot loops already use the precomputed modulo-free `XClip`. Only bites partial-arc segment clips (multi-board Phantasm). *Fix:* build an `XClip` once per draw at the sink, mirroring `scan.h`. (`core/filter.h:142-147,166-170`)

7. ✅ **Redundant per-pixel y-bounds re-check in `plot_virtual`.** After `contains_y(y)` passes, `y` is provably in `[0,h)`, so `plot_virtual`'s `y>=0 && y<height()` re-test is always-true dead work on every composited pixel. *Fix:* drop the guard at the (already y-clipped) sink call sites, or inline the write. (`core/filter.h:94-98,142-146,158-170`)

8. ✅ **Stray trailing semicolons after member-function bodies.** `PipelineRef::plot(...) const { ... };` and `~Effect() { s_alive = false; };` are harmless empty declarations (no warning under the current flags). *Fix:* remove the trailing `;`. (`core/concepts.h:302`, `core/canvas.h:65`)

9. ✅ **`DistortedRing` `max_distortion` is a load-bearing bound with no device-build guard.** The 256-sample sanity check is wrapped in `#ifndef NDEBUG`, so on device/wasm an underestimate silently culls genuine arcs — the lone SDF shape whose precondition is NDEBUG-gated while every sibling uses always-on `HS_CHECK`. Current caller passes an exact bound. *Fix:* replace the `#ifndef NDEBUG` sample check with an always-on cold-path `HS_CHECK`. (`core/sdf.h:725-761`)

10. ✅ **Duplicated `STAR_INNER_RATIO` constant across `sdf.h` and `plot.h`.** The filled star and stroked star share `0.382f` by hand-copy ("mirrors SDF::…"), enforced by nothing; the only test pins Plot's value against itself. *Fix:* hoist one definition into `constants.h` (both headers include it) and delete the duplicate. (`core/plot.h:27`, `core/sdf.h:33`)

11. ✅ **`sort_intervals_by_start` uses a signed `int` index over a `size_t` count.** Casts to `int` to decrement to `-1` as the shift terminator in a shared helper; capacities are tiny so no overflow today, but it couples correctness to `N < INT_MAX`. *Fix:* rewrite with an unsigned-safe idiom (`j` from `i` down while `j>0`). (`core/sdf.h:135-143`)

12. ✅ **`SmoothUnion` vertical-bounds fold leaks a culled child's `{1,0}` sentinel into the band.** A culled child folded via `min/max` can annex rows neither child occupies. Purely conservative over-scan (no coverage/visual defect); the plain `Union` fold shares it. *Fix:* short-circuit `y_min>y_max` children before folding, mirroring the existing `lo>hi` case. (`core/sdf.h:993-1005`)

13. ✅ **`MobiusWarp`/`MobiusWarpCircular` `set_scale()`/`set_easing()` are dead API.** Four setters with no callers anywhere (the live path is `bind_scale()`); one even documents a contract no caller exercises. *Fix:* delete the four unused setters (keep `bind_scale`). (`core/animation_mesh.h:154-158,201-205`)

14. ✅ **Indefinite `Sprite` (`duration = -1`) never fires `.then()`, undocumented.** Peer perpetual animations (`RandomWalk`, `MobiusWarpEvolving`) warn about this exact footgun; `Sprite`'s doc mentions the indefinite case only for fade-out. *Fix:* add the one-line `.then()`-never-fires note to `Sprite`'s doc. (`core/animation_scalars.h:437-471`)

15. ✅ **`MeshMorph` nearest-vertex correspondence assumes unit-length vertices.** "Greatest dot == nearest" and `slerp` both require unit inputs; all current mesh sources are on the unit sphere, but the invariant is implicit/unchecked. *Fix:* add a one-time cold-path `HS_CHECK` for unit-length at construction (or normalize before the dot/slerp). (`core/animation_mesh.h:308-320`)

16. ✅ **`ParticleSystem` event-horizon steering can stall a friction-drained particle.** Inside the horizon it redirects existing speed without adding acceleration, and friction has already damped it, so a slow particle can fail to reach `kill_radius` before its speed collapses (it is still reclaimed by the life timer). *Fix:* add a minimum inward pull / floor the steering speed, or document the friction stall. (`core/animation_orientation.h:345-349`)

17. ✅ **Stale arena-budget figure in header comment (335 KB vs the real 330 KiB).** The load-bearing sizing comment's headline contradicts `DEVICE_GLOBAL_ARENA_SIZE`, the rest of the same block, and the README. *Fix:* change "335 KB" to "330 KiB". (`core/memory.h:17`)

18. ❌ (already de-drifted — the edge keying/pairing goes through the same shared `fill_edge_record`/`pair_half_edges` helpers the `HalfEdgeMesh` builder uses; only a lean face-neighbor scan remains, and a full `HalfEdgeMesh` build would raise cold-path scratch high-water ~4I→~10I+ for no drift benefit) **`classify_faces_impl` reimplements half-edge pairing that `HalfEdgeMesh` already computes.** ~30 lines of parallel connectivity build (`he_to_face`/`pair_array`/records) that can drift from the single tested builder. *Fix:* build a `HalfEdgeMesh` and read neighbor faces via `half_edges[half_edges[hidx].pair].face`; re-measure scratch high-water. (`core/mesh.h:637-687`)

19. ✅ **`narrow_index()` (uint16_t) stored into a local `int` in expand/chamfer/snub.** Value-preserving but inconsistent with `ambo`/`truncate` and the `uint16_t` sinks; obscures the topology-index type intent. *Fix:* declare `uint16_t idx = narrow_index(...)`. (`core/conway.h:764,853,1081`)

20. ✅ **`CompiledHankin::clone` duplicates the shared `copy_vector` deep-copy loop.** Uses a local `push` lambda instead of the shared helper that `MeshState::clone`/`MeshOps::clone` use — a third arena-deep-copy idiom that can drift. *Fix:* route all six vectors through `copy_vector`. (`core/hankin.h:66-79`)

21. ✅ **Redundant `params = presets.get()` in `Flyby::init` before the first preset lerp overwrites it.** The immediately-scheduled `Lerp` re-derives the start from `prev_get()`, so the assignment is never observed on the render path (dead). *Fix:* drop the line (per prefer-delete-dead-code). (`effects/Flyby.h:69-70`)

22. ✅ **`BZ::blend_species` unclamped `uint16` cast relies on an implicit precondition.** Safe as a true convex blend for the sole caller, but a future caller passing pre-scaled weights would overflow silently; the doc omits the precondition. *Fix:* tighten the doc to state `a+b+c` must equal the normalizing `wsum` (or clamp before the cast). (`effects/BZReactionDiffusion.h:307-325`)

23. ✅ **`HopfFibration` records trails every frame while paused, collapsing them to a point.** With flow/tumble Drivers frozen, the trail ring fills with identical points over paused frames and each polyline becomes a dot; regrows over `TRAIL_LEN` frames on resume. *Fix:* gate the record loop on `!animationsPaused()` (keep rendering). (`effects/HopfFibration.h:91-95`)

24. ✅ **`Moire` wipe-rebake gate never turns off in steady state (dead logic).** Wipes run back-to-back so `wipe_frames_remaining_` is re-armed before reaching 0; the per-frame decrement/guard only saves work before the first wipe. *Fix:* drop the gate and always rebake (the LUT rebakes are cheap, per the existing comment). (`effects/Moire.h:94-99`)

25. ✅ **`SplineFlow` `MAX_TRAILS=30000` reserves ~240 KB of the 330 KB device arena as a bare literal.** ~80% of the persistent budget with no compile-time guard, unlike the project-wide convention (HopfFibration/MindSplatter/Voronoi/GS/BZ all `static_assert` their footprint). *Fix:* add `static_assert(MAX_TRAILS*sizeof(Item) <= budget)` mirroring the sibling effects, and a cost comment. (`effects/SplineFlow.h:21`)

26. ✅ **`Voronoi::shade()` takes the nearest index as `int` while KD/CellId use `uint16_t`.** Harmless at `MAX_SITES=400`, but diverges from the deliberate `uint16_t` index typing elsewhere in the file. *Fix:* change `i0`/`i1` to `uint16_t`. (`effects/Voronoi.h:98`)

27. ✅ **`HD107SFrame::load()` is effectively test-only surface in the shipped header.** Shipped drivers pack via `packPixel()` only; `load()` (a second flush-owning CRGB path) is exercised only by a parity test and can be mistaken for a device path. It is a load-bearing parity guard, so keep it. *Fix:* add a one-line "host/test bulk path, not the device path" marker to `load()`'s doc. (`hardware/hd107s_frame.h:156-185`)

28. **`TeensySPIDMA` single-instance invariant is enforced only at runtime.** A hard global-singleton constraint (backing the shared DMA-completion ISR dispatch) is invisible at the `DMALEDController<N>` API; two controllers compile fine and trap only at the second `begin()`. Never fires on current single-controller targets. *Fix:* document the one-per-image limit in `DMALEDController`'s class doxygen. (`hardware/dma_led.h:75-100`)

29. **`MeshOpsWrapper`/`PaletteOps` expose implementation state as public data members.** `mesh`/`generation_`/`lut` are public though embind binds only the methods; broader than necessary and weakens the `generation_` invariant. *Fix:* make the members private. (`targets/wasm/wasm.cpp:740-746,1168`)

30. **`kMinAssertions` floors in the auxiliary TUs require manual bumping.** Hand-counted magic numbers (26, 108) that drift unless someone remembers the manual step; the guard only catches wholesale gutting. *Fix:* derive the floor from the shared case list where feasible, or at least keep the "bump when adding" note tight. (`tests/fastmath_clamp_check.cpp:24-31`, `tests/h_offset_renorm_check.cpp`)

31. **`led.h` correction-guard double-live trap has no death coverage.** The `correction_guard_depth()` single-live trap is analogous to covered Timeline/Effect traps but is not exercised. (The host `FastLEDMock` in `platform.h` makes a case runnable — no shim needed.) *Fix:* add a `case_correction_guard_double_construct` death case. (`core/led.h:89-91,114-116`)

32. **`setCanvasSize()` discards the user's zoom (orbit radius) on every resize.** Each `ResizeObserver` callback runs `camera.position.setLength(fit)`, so any reflow/DPR change/sidebar toggle/devtools open snaps a user zoom back to the ~85%-coverage distance. *Fix:* only re-fit distance on first sizing / when aspect actually changes, or clamp the existing radius to the new min/max. (`daydream/driver.js:373`)

33. **`setEffect` can publish stale (old-effect) `paramValues` from an in-flight render.** `setEffect` nulls `paramValues` but does not bump `renderGen`, so a serialized in-flight old-effect frame republishes old-ordered values and `syncGUI` binds the rebuilt sliders by index to them for one frame (only when the two effects share a param count; self-heals next frame). *Fix:* bump `renderGen` in `setEffect()` so the fence drops the in-flight frame. (`daydream/segment_controller.js:454-460,250-251`)

34. **`onWorkerFault` mislabels a render-watchdog fault as "pool init".** The overlay headline is chosen purely from `segId < 0`, but the render watchdog also passes `-1`, so a mid-run render timeout shows as a pool-init fault (the message body still disambiguates). *Fix:* carry an explicit `phase`/`kind` field, or use a non-negative render-timeout sentinel. (`daydream/segment_controller.js:403-423,687`)

35. **Worker "unknown message type" only warns; a protocol-drift `setX` is silently dropped.** Contradicts the file's own fail-fast posture (it throws on pixel/bounds mismatch, posts `initFailed` on a bad resolution). A cross-version-drift state-changing message keeps the worker rendering stale under the current generation, invisible to the fence. *Fix:* throw (→ `onerror` → fault) or post an `initFailed`-style fault on the default case. (`daydream/segment_worker.js:219-221`)

36. **Listbox keyboard nav omits Home/End.** The WAI-ARIA listbox pattern expects jump-to-first/last; today they fall through to page scroll with no focus move. *Fix:* add `Home→0` / `End→len-1` to `navTargetIndex` (the wiring already actions any non-`-1` return) and extend the test. (`daydream/sidebar_logic.js:40-49`)

37. **Both scroll arrows hidden when content overflows by less than the 4px deadzone.** For `0 < maxScroll ≤ 4`, neither arrow shows though the list is scrollable (cosmetic; scrolling still works). *Fix:* gate the deadzone against `maxScroll` (e.g. `scrollLeft < maxScroll - Math.min(4, maxScroll)`), add a small-overflow test. (`daydream/sidebar_logic.js:78-85`)

38. **`label_format` constant `g` is an opaque name for the inverse golden ratio.** `g = 1/PHI`, then `1/g` is used to recover `PHI` — a confusing double-inversion (and needless reciprocal round-trip) in a file whose purpose is symbolic clarity. *Fix:* use `PHI` / `1/PHI` (or `INV_PHI`) directly and drop `g`. (`daydream/label_format.js:12`)

39. **Sidebar click handlers and appended DOM nodes are not detached in `dispose()`.** `dispose()` removes the observers/listeners but leaves per-button `onclick` closures and appended nodes; harmless today (app-lifetime singleton) but a latent asymmetric-teardown footgun if the container is ever reused. *Fix:* make teardown symmetric (clear `innerHTML`, null `onclick` refs, remove appended nodes) or document that `dispose()` assumes the container is discarded. (`daydream/sidebar.js:82-88`)

---

## Rejected findings

20 candidate findings were raised and rejected by independent validation — recorded here so the
ledger is complete and they are not re-litigated. The dominant reasons were **intentional documented
design** (the `HS_CHECK` fail-fast philosophy, ambient-vs-parameter pause gating, per-`N` shared
color-correction state, `fmodf` argument reduction for large-`i` trig precision) and **false
premises** (elements the app never inserts dynamically, arm Y-bands that are in fact always aligned,
count-invariants set atomically in one place, docs that already contain the claimed-missing note).

Notable examples:

- **Pause-gating "inconsistencies"** in GnomonicStars / Raymarch / ShapeShifter — the documented
  contract (`canvas.h:476-490`) gates only *parameter-driving* animations; ambient rotation/camera/
  palette is intentionally left running. The flagged effects follow the rule.
- **`fmodf` before periodic trig in `fib_spiral`** — a legitimate argument reduction that preserves
  float32 precision for large indices (the golden-angle test pins it); not dead work.
- **`getFaces`/internal `HS_CHECK` traps** — always-on fail-fast on an internal invariant is the
  intended design, distinct from the reject-and-return pattern used for untrusted JS inputs.
- **`StaticCircularBuffer` single-arg constructor**, **`srgb_to_linear_interp` secant bias**,
  **per-`N` shared `HD107SFrame` correction state**, **`setClip`-after-`setResolution`** — each is an
  intentional, already-documented (and often already-tested) design decision.

---

## Summary

The confirmed defect profile — 0 critical, 0 high, 2 medium, 37 low, across a codebase spanning
firmware, a rendering engine, and a web simulator — is characteristic of a project that has already
had its serious bugs found and fixed. The two Medium items are both **test-infrastructure gaps** (an
unpinned death-case count and a missing aggregate-arena CI gate), not product defects; closing them
tightens the safety net the project already values. The Low items are best worked as small,
well-scoped `code-review-fix` commits. The engine's grade is earned: it is more disciplined,
better documented, and more thoroughly tested than the large majority of both hobbyist and
professional embedded-graphics code.
