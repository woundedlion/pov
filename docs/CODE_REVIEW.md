# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`,
`targets/wasm/`), the `daydream` web simulator (`c:\work\daydream`), the unit-test suite
(`tests/`), and the build/CI/tooling layer. Out of scope by instruction: `effects_legacy.h`,
`targets/Holosphere/Holosphere.ino`, and `core/rotate.h`.

**Method.** The codebase was partitioned across 18 independent review sub-agents (one per
subsystem cluster), each grounded in the README architecture sections relevant to its scope.
Every substantive finding was then handed to a second, independent set of verification
sub-agents that re-read the cited code and returned a CONFIRMED / PARTIAL / REFUTED verdict.
Six reported findings were **refuted** on re-examination and are listed at the end so the
record is honest; the prioritized list below contains only verified defects.

---

## 1. Overall Assessment

**Overall grade: A−**

This is an exceptionally well-engineered codebase — closer to a disciplined professional
product than to typical hobby/art firmware. The defining traits are (1) compile-time
correctness: invariants are pushed into `static_assert`s, template traits, and single-source
X-macros so whole classes of drift are unrepresentable; (2) a coherent fail-fast philosophy
(`HS_CHECK`) applied deliberately to cold seams; and (3) documentation that genuinely explains
*why*, frequently anticipating the exact failure mode a reader would worry about. No critical
or data-corrupting defects were found. The verified issues are edge cases, prose-only invariants
that lack a static guard, a small number of real UB/parity bugs, and test/CI coverage gaps —
none of which undermine the architecture.

The recurring *theme* across the findings is consistent: the project's strongest guarantees are
compiler- or trap-enforced, while a thin residue of equally load-bearing invariants survives only
in comments (arena-span re-grow, ISR memory ordering, get_pixel-override assumptions, palette/slider
range agreement, the segment-overlay seam). That residue is exactly where the next regression is
most likely to hide.

---

## 2. Quality-Dimension Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A− | Math, geometry, color, sync, and double-buffer logic are overwhelmingly right; the residue is narrow edge cases plus a real UB (`ArenaSpan::end()`), a release-build dangle window, a WASM-abort vector, and one tool/engine parity divergence. |
| **Robustness / error handling** | A | Pervasive, deliberate fail-fast on cold seams with bounded soft-handling for transients; the only gaps are debug-only safety nets that vanish in the device build and a few prose-only invariants. |
| **Readability / clarity** | A | Precise naming and intent-first comments throughout; the sole drag is comment density occasionally burying short code. |
| **Maintainability** | A− | Strong single-sourcing (X-macro rosters, derived budgets, shared scaffolds); offset by intentional sibling-effect duplication and invariants enforced only by prose. |
| **API / interface expressiveness** | A | Compile-time variadic pipeline, owned-vs-borrowed memory types, RAII guards, and type-erased callables make correct use the easy path; a few footguns remain (the slider factory, `Coords`/`Wrap` pairing). |
| **Architectural elegance** | A | The flywheel sync core, the variadic filter pipeline with compile-time domain conversion, the tri-arena model, and compile-time resolution parameterization are genuinely elegant and mutually coherent. |
| **Performance** | A | LUT-backed trig, branchless ISR, zero-copy readback, zero steady-state allocation; a couple of per-pixel `acosf` hot spots and unbudgeted heavy effects are the only soft spots. |
| **Documentation** | A+ | The README is reference-grade (including a full 1-wire sync datasheet) and inline rationale is best-in-class; a handful of docblocks overstate or omit an approximation. |
| **Testability** | A− | Designed for host testing (pure-math seams, white-box friends, a death harness that proves traps fire); weakened by self-consistency-only determinism checks and some weak assertions. |
| **Test-suite quality** | A− | Predominantly discriminating, oracle-driven assertions and an exemplary event-driven sync simulator; pulled back by tautology-prone spots, hand-maintained reset/budget lists, and coverage gaps. |
| **Build / CI / release engineering** | A | Fanatically pinned, self-testing gates (the size gate has its own fixtures; generator output is diffed in CI); the screenshot gallery and one PlatformIO cache are the only unpinned/unverified artifacts. |

### Per-component grade summary

| Component | Grade |
|---|---|
| Core math / geometry (`3dmath`, `geometry`, `util`, `easing`, `waves`, `concepts`) | A |
| Color system (`color`, `palettes`, `styles`) | A |
| Memory / platform (`memory`, `platform`, `static_circular_buffer`, `inplace_function`) | A− |
| Rasterizers (`sdf`, `scan`, `plot`) | A− |
| Filter / animation / transformers | A |
| Mesh system (`mesh`, `conway`, `hankin`, `spatial`, `solids`, `generators`) | A |
| Engine glue (`canvas`, `engine`, `effects`, `effect_registry`, `led`, `presets`) | A− |
| Hardware drivers (`dma_led`, `hd107s_frame`, `pov_*`, `pov_sync`) | A |
| WASM bridge (`wasm.cpp`, `param_marshal`) | A |
| Effects (28 effect headers) | A− |
| Unit-test suite | A− |
| daydream core JS | A− |
| daydream segmented-worker pipeline | A |
| daydream geometry tools | B+ |
| Build / CI / tooling | A |

---

## 3. Prioritized Items to Fix

No **Critical** defects were found. Items are numbered sequentially and grouped by priority.

### High Priority

1. ✅ WASM exported free functions (`spline_cubic_fast/slerp`, `spline_catmull_rom_tangents`, `lissajous`, `procedural_palette_linear`, the OKLab/OKLCH converters) pass raw JS floats straight into engine math with no `isfinite` gate, unlike every class method (which uses `finite_arg`). A NaN/Inf from a tool page flows into the slerp spline path and trips `Vector::normalized()`'s `HS_CHECK`, aborting the *entire* WASM module — the exact JS-boundary trap the rest of the bridge is written to avoid. (`targets/wasm/wasm.cpp` ~1366–1485)
2. ⏭️ *(Won't-fix — documented intentional tradeoff.)* `ArenaVector::bind()`'s in-place grow path reallocates `data_` without bumping any release-visible generation; the dangling-`ArenaSpan` detector (`rebind_generation_`) is `#ifndef NDEBUG`, so on the shipping device build an `ArenaSpan` taken before a grow silently dangles with no detection — the one place the memory-safety net is absent is exactly where corruption is worst. (`core/memory.h` ~466–513, ~753) *Decision: left debug-only. The only fix is a release-visible generation check inside `ArenaSpan`'s element accessors (`operator[]`/`data`/`begin`/`end`), which adds an always-on branch to span reads on render hot paths — contradicting the code's documented deliberate choice to keep this tracking debug-only for footprint/perf. The dangle is a stale-read (not corruption) until the next arena reset/compaction.*
3. ✅ The Lissajous tool diverges from the engine: `lissajous()` ends with `.normalize()` while the engine returns the raw (already-unit) vector, and the JS comment falsely claims it "mirrors the engine's `lissajous()`, which returns `v.normalized()`." Tools exist to predict device output; fix both the divergence and the inverted comment. (`daydream/tools/lissajous_math.js` 46–51 vs `core/geometry.h` 821–826)
4. ✅ `ArenaSpan::end()` returns `data_ + size_` with no null guard, so a default-constructed span (`data_=nullptr, size_=0`) iterated by range-for performs `nullptr + 0` — formal UB / UBSan-flagged. `ArenaVector::end()` already guards this identical case; mirror it. (`core/memory.h` 848–851)
5. ✅ `setClip(y0, y1, x0, x1)` orders the Y pair before the X pair, opposite to the (x,y) convention; embind binds positionally, so a transposed JS call passes the range check whenever the values fall within both axis bounds and silently clips the wrong axis. (`targets/wasm/wasm.cpp` ~453)

### Medium Priority

6. ✅ The Plot `v1` register (documented in README §7.0 as "cumulative arc length in radians") is computed as the geodesic great-circle chord even for the always-planar primitives (`Star`, `Flower`, `PlanarPolygon`), so `v1` is inconsistent with the drawn planar position; any shader keying off `v1`/`v0` as an arc-length proxy gets wrong values. (`core/plot.h` ~614–638; `Flower` self-documents the seam)
7. ✅ With `persist_pixels`, `advance_buffer()` copies the trail base from `bufs_[prev_]` (the frame the ISR is *displaying*) rather than `bufs_[next_]` (the last *completed* frame); this is correct only because the `buffer_free()` gate forces `prev_ == next_`. Copy from `next_` (or assert `prev_ == next_` at the site) to remove the cross-method coupling. (`core/canvas.h` 254–257)
8. ✅ `HD107SFrame::correct()` chains six `(v*f)>>8` stages with no output clamp; LUT-index safety rests entirely on `factor() <= 256` with no `static_assert`. Widening `factor()`'s domain (e.g. a future >1.0 gain) would silently over-read `linear_to_srgb_lut`. Add a `static_assert`/clamp. (`hardware/hd107s_frame.h` 130–145)
9. ✅ The segmented ISR `render_column` hard-assumes no Phantasm-roster effect overrides `get_pixel` (`pov_single.h` actively probes `overrides_get_pixel()`), enforced only by prose; an override-using effect added to the roster would silently read raw-buffer pixels. Add a `static_assert`/registry guard. (`hardware/pov_segmented.h` 644–647)
10. ⏭️ *(Won't-fix — inherent to gnomonic projection.)* `gnomonic_mobius_transform` selects the output hemisphere from the sign of the *input* `v.y`, so a Möbius map that moves a point across the equator re-projects it into the wrong hemisphere — an unguarded discontinuity/seam (the sibling stereographic path has no such branch). (`core/transformers.h` 284–289) *Decision: gnomonic is a single-hemisphere chart whose equator maps to infinity (handled by the `STEREO_INF` clamp → nearer pole). A Möbius map is an automorphism of the extended plane, so it maps a hemisphere chart to itself and cannot move a point into the other hemisphere within the chart; restoring the input-`v.y` hemisphere is the only well-defined choice. The equator seam is intrinsic to any gnomonic-based warp, not a fixable defect, and the behavior is pinned by the transformer tests.*
11. ✅ Voronoi hard-codes a 64 KB per-frame scratch budget that simultaneously holds `positions`, the KDTree (scales with the site-count slider up to 400), and `cells` (scales with `(W/B)*(H/B)`), with **no** `static_assert` on the combined high-water — unlike BZ/GS, which derive and `static_assert` their budgets from `sizeof`. It relies on the runtime OOM trap instead. (`effects/Voronoi.h` 46, 200–241)
12. ❌ *(Refuted — overstated.)* Comets palette rollover is perpetually starved at low Cycle Dur: the timer period (`2*cycle_duration`, min 20) is shorter than the wipe lockout (`WIPE_FRAMES+1` = 49), so `update_palette()` always early-returns while the path keeps switching — the palette freezes and the "next tick picks it up / self-heals" comment is optimistic. (`effects/Comets.h` 113–121, 241–260) *Re-examination: `wipe_frames_remaining_` decrements every frame independent of the cycle timer, so it always drains to 0 within 49 frames. At the Cycle Dur floor the timer fires at frames 20/40 (skipped, wipe in flight) and 60 (drained → proceeds), rolling the palette over every ~3rd path switch. The palette does **not** freeze; it lags the path cycle ~3:1 — exactly the bounded self-heal the existing comment documents. The "always early-returns / freezes" claim is incorrect.*
13. ✅ `SDF::SmoothUnion<Union, Union>` (nested-CSG children) can overflow the shared `MergedIntervalBuffer` and trap at runtime; the type system permits a composition the runtime forbids. The trap is intentional fail-fast, but a nesting-depth `static_assert` would convert a runtime abort into a compile error. (`core/sdf.h` 873–928)
14. ✅ The MeshFeedback preset cycle applies field values that may exceed their registered slider ranges (Fade max 0.99, amplitude 30, scale 50), so the live parameter can leave the registered `[min,max]` — a range/data-consistency gap. Confirm every preset's driven fields fall inside the registered ranges. (`effects/MeshFeedback.h` 77–98) *Verified: all eight presets already fall inside the registered ranges. Hoisted the ranges to constexpr constants shared by `registerParam` and a new `static_assert` over every preset, so the agreement is now compile-time enforced; the guard's resolution path is to widen the range to accommodate a future preset, not to clamp the preset.*
15. ✅ Flyby and Liquid2D expose exactly 10× different `Warp Scale` / `Warp Strength` slider ranges feeding the *same* `stereo_noise_warp` core, with no documented justification for the divergence. (`effects/Flyby.h` 38–39 vs `effects/Liquid2D.h` 44–45) *Fixed: the divergence is intentional — Flyby's presets sweep warp to scale ~47.8 / strength ~11.6 and need the wide range, while Liquid2D pins warp at 1.5/0.5 across its presets and keeps a tighter range for finer control. Documented the justification at both register sites rather than unifying (unifying down would push Flyby's presets out of range).*
16. ✅ CI never runs or validates `capture_screenshots.mjs`, yet `docs/screenshots/` is installed into daydream and served live; committed PNGs can silently rot relative to the current rendering, and the gallery is the one artifact in an otherwise fanatically-pinned pipeline with neither reproducibility pins nor a freshness gate. (`.github/workflows/*`, `CMakeLists.txt` ~131, `scripts/capture_screenshots.mjs`) *Fixed: a full per-build re-capture is too slow, so added the cheap freshness gate instead — a `screenshot-gallery` CI job that `node --check`s the capture script and asserts the committed PNG set matches the `HS_EFFECT_LIST` roster (every effect has a PNG, no orphans). The roster parse is factored into `scripts/effect_roster.mjs`, shared by the capture script and the new `scripts/check_screenshots.mjs` checker so they cannot drift. The gate immediately caught two real missing PNGs (`DistortedRing`, `ShapeShifter`, never captured after the Test/TestShapes rename), now captured at the gallery resolution. It cannot detect a visually-stale PNG; that still relies on a manual re-capture.*
17. ✅ The PlatformIO toolchain cache step keys on `hashFiles('platformio.ini')` with **no** `restore-keys` fallback (unlike the build-objects cache directly below it), so any `platformio.ini` edit forces a full multi-hundred-MB toolchain re-download. (`.github/workflows/ci.yml` ~462–468) *Fixed: added `restore-keys: pio-teensy-` so a pin bump seeds the restore from the newest prior toolchain cache and pays only for changed packages. It can't mask a pin change — `pio run` installs exactly the pinned versions on top of the seed and the new exact key re-saves the result.*
18. ✅ The effect determinism pass renders each effect twice from identical reset state and asserts byte-equality — self-consistency only. It is structurally blind to nondeterminism that is a pure function of identical inputs (e.g. a `static` seeded once and never reset), and the reset list is hand-maintained with nothing asserting its completeness. (`tests/test_effects.h` 305–333) *Fixed (completeness): scramble every output-affecting global render_capture() resets (RNG seed, global_timeline_t, hue cursor) BETWEEN the two captures, so run B must recover canonical output from a dirtied state instead of re-running from pristine state twice. Each reset line is now load-bearing and self-testing — verified by negative control: removing the hue-seed reset now fails the test (DistortedRing/FlowField diverge), where the old pristine-twice check passed. The static-seeded-once class remains out of reach in-process and is documented as such; not fixable without a fresh process per run.*
19. ✅ Islamic-family solids are validated only by `check_basic` (non-empty / consistent / finite) with no exact V/E/F oracle (unlike Platonic), so a wrong-but-self-consistent generator passes; the Hankin one-shot asserts only `size() > 0` and uses a looser 1e-2 unit-sphere tolerance versus 1e-3 in `compile_hankin`. (`tests/test_solids.h` 149–158, `tests/test_hankin.h` 341–348) *Fixed: probed all 23 Islamic entries — every one is a closed 2-manifold (`V-E+F==2`), so added `test_islamic_registry_solids_are_closed` enforcing the Euler invariant (the topological equivalent of the Platonic exact-count oracle; catches a generator that opens a seam, drops a face, or duplicates geometry). Exact per-entry counts are deliberately NOT pinned — the generators are actively tuned, so golden counts would invert the signal the way `test_effects.h` rejects golden hashing. The Hankin one-shot already asserts more than `size()>0` (face-count consistency, in-range indices, unit-sphere); its 1e-2 tolerance is justified and now documented — it covers the dynamic ray-intersection vertices (off-sphere by construction), unlike `compile_hankin`'s 1e-3 which covers only the exactly-normalized static edge-midpoints.*
20. ✅ The non-finite transformer passthrough test asserts only that the result *stays* non-finite, not that it equals the input — a buggy path returning all-NaN for any input would pass while the docstring claims "passes through identity." (`tests/test_transformers.h` 369–375) *Fixed: replaced the three `!finite_vec(...)` checks with `vec_bits_equal(..., v)`, a new bit-for-bit comparison (raw IEEE-754 bits, since NaN != NaN under ==) that proves each short-circuit returns its input verbatim. Verified passing.*
21. ✅ The param-marshal by-name round-trip is silently skipped (`if (target >= 0)`) for effects with no editable float param, with no log or counter, so roster drift toward such effects erodes coverage invisibly. (`tests/test_param_marshal.h` ~89) *Fixed: `check_one` now returns whether it exercised the round-trip; the runner tallies and prints the split (currently 22/27 exercised, 5 skipped) and `HS_EXPECT`s coverage stays > 0, so erosion is visible and a total drift to zero fails loudly.*
22. Canvas out-of-bounds access and inverted/out-of-canvas clip bounds are untested, and — unlike `test_static_circular_buffer.h`, which explicitly defers OOB to the death harness — there is no comment documenting where the fail-fast trap is covered. (`tests/test_canvas.h` 317–331, 447–463)
23. `child_trapped` accepts `WIFEXITED && WEXITSTATUS == 128 + SIGILL` in addition to the exact signal path, widening the accept set so a child that genuinely `exit(132)`s reads as a passing death test (bounded, but it erodes the "specific status, not merely nonzero" guarantee). (`tests/test_death.h` ~850–852)

### Low Priority

24. Dynamo's per-fragment `color()` performs an `acosf` (`angle_between`) plus a linear band scan on the flush hot path (trail capacity 10000) — the most expensive per-pixel op in the effect set despite the baked-palette LUTs. A dot-product band test (`cos(boundary±width)` vs `dot(v,normal)`) would remove the `acos`. (`effects/Dynamo.h` 211–266)
25. ✅ RingSpin's trail length-fade comes entirely from a transparent-vignette palette sampled at `1-t`, with no `quintic_kernel(t)`/`t` fade in the draw loop (unlike sibling Comets) — a silent dependency that disappears if the palette is swapped for an opaque one. Add a comment or assert at the draw site. (`effects/RingSpin.h` 123–126) *Fixed: documented the load-bearing dependency at the draw site — the per-trail-sample comment now states the fade is supplied entirely by the palette's alpha vignette (sampled at 1-t), that there is deliberately no kernel/t fade here, and that swapping to an opaque palette would render the whole trail at full alpha; keep the vignette or reinstate an explicit t-fade.*
26. ✅ Moire's `Transition(rotation, 2π).then(rotation = 0)` reset is seamless only because the specific two-axis rotation product (and the freq-4 wobble) are 2π-periodic; no comment asserts this, so an axis/angle/coefficient change could silently introduce a snap. (`effects/Moire.h` 75–77, 99–102, 159–160) *Fixed: documented the invariant at the `Transition` site — the snap from 2π back to 0 is seamless only because every consumer of `rotation` (the fixed-axis `make_rotation` products in `draw_frame()` and the freq-4 `sin_wave` wobble in `draw_layer()`) is exactly 2π-periodic; the note states the endpoint must stay an integer multiple of 2π and the wobble frequency an integer.*
27. `make_rotation(from, to)` silently returns identity for rotations below ~0.81° (`d > 1 − TOLERANCE`); documented, but accumulated tiny steps net zero — a footgun for any future caller needing exact small rotations. (`core/3dmath.h` 1078–1117)
28. Voronoi's coarse-coherence optimization (8-px block) silently drops Voronoi cells smaller than the block at high site counts ("a dropped speck"); documented and bounded, but a real fidelity edge. (`effects/Voronoi.h` 200–205)
29. `formatFloatCpp` uses `toFixed(6/3)`, which rounds small nonzero coefficients (e.g. `1e-7`) to `"0.0f"`, silently zeroing meaningful small values in generated C++. (`daydream/tools/cpp_format.js` ~23)
30. `findBestRationalRatio(0)` returns `{M:1, N:1}`, so a deliberately-zeroed Lissajous frequency snaps up to the passive frequency instead of staying zero; it should return `{M:0, N:1}`. (`daydream/tools/lissajous_math.js` ~63)
31. `createSlider`'s `onInput` never updates `valueSpan.textContent`, so the returned readout freezes while the slider moves unless every caller updates it — a footgun in a reusable factory. (`daydream/tools/slider.js` 104–106)
32. `segment_controller` `composite()`'s boundary overlay paints column `x==0` cyan as a "wrap seam," overwriting arm A's genuine first-column image data (cosmetic; only when boundaries are shown with an x-split layout). (`daydream/segment_controller.js` 631, 647–650)
33. `setResolution`'s faulted-recovery rebuild only fires when `(faulted && active)`; a faulted-but-inactive controller falls through to broadcast `setResolution` to a dead worker pool. (`daydream/segment_controller.js` 512–515)
34. `read_id`'s three-sample debounce only rejects straps caught mid-transition; a stable-but-wrong or shorted ID strap silently becomes a second master and drives the push-pull `SYNC` wire into bus contention — undetectable from one board. Documented as out-of-scope, but it is the one fail-*wrong* path in an otherwise fail-dark design. (`hardware/pov_segmented.h` 403–435)
35. `pov_sync::position()` reinterprets `(at − epoch)` as `int32`; the documented ~±3.5 s safe window has near-zero margin to `INT32_MAX` and is unenforced — safety rests on the rebase invariant rather than a guard. Tie it to an `HS_CHECK`/static bound on `cycles_per_half_rev`. (`hardware/pov_sync.h` 611–616)
36. `queue_frame` uses relaxed atomic stores across the main-loop→ISR boundary; pixel-payload visibility rests entirely on `hs::disable_interrupts()` acting as a compiler barrier, not on the atomics — correct on the single-core target but load-bearing on a `platform.h` implementation detail with no local assertion. (`core/canvas.h` 263–268, 615–628)
37. The self-registering `REGISTER_EFFECT` order in `EffectRegistry::entries()` is link/TU-order dependent. No consumer indexes it positionally today (dispatch is by name; rosters use the `HS_EFFECT_LIST` single source of truth), but a future positional/persisted/transmitted index would be build-order-fragile. (`core/effect_registry.h` 93–96)
38. `configure_arenas` performs no guard against repartitioning while live allocations exist; it is an init-only contract today, but a mid-run call would silently orphan persistent allocations. Add an `HS_CHECK(persistent_offset == 0)` precondition. (`core/memory.cpp` 78–100)
39. `Persist<T>` correctness relies on the un-encodable precondition that `T::clone` allocates only in its scratch argument, and its watermark check is an aggregate backstop that can miss a per-`Persist` mistake in stacked use. Both limitations are self-documented; consider encoding the precondition in the `Cloneable` concept. (`core/memory.h` 1002–1046)
40. `Face::distance` converts a gnomonic planar distance to angular via `fast_atan2(d, 1)` (exact only at the face center) and its docblock calls the result "angular distance" with no approximation note, unlike `PlanarPolygon`'s. (`core/sdf.h` 2156, 2081–2087)
41. `CycleModifier`/`ScaleModifier` emit out-of-`[0,1]` coordinates that depend on `StaticPalette`'s default `Wrap=true`; pairing either with a `Wrap=false` composition (a documented but narrow opt-out) samples the source out of range with no static guard against the combination. (`core/color.h` 1685, 1909)
42. `coordsLabel` allocates a `THREE.Spherical` + `Vector3` + template string per call and runs in the per-rendered-frame label refresh; latent per-frame GC churn that would defeat the pooled-label design if an effect emitted many labels. (`daydream/driver.js` 899–907)
43. The sidebar's initial sort indicator shows the neutral `'Name ⇅'` glyph even though Name-ascending is the active default; the directional arrow appears only after the first sort click. (`daydream/sidebar.js` 29, 199)
44. `URLSync` writes a trailing `?` to the address bar when all query params are cleared (cosmetic). (`daydream/state.js` 248, 288–289)
45. `ParticleSystem::step`'s swap-remove relies on `size_t` unsigned wraparound (`i--` at `i==0` → `SIZE_MAX`, then `++i` → 0); it is well-defined and correct but a fragile idiom that breaks under a signed-type refactor or inserted logic. (`core/animation.h` 612–618)
46. BZ `perturb_state` draws from the shared RNG every substep, coupling the global RNG stream position to `STEPS_PER_FRAME`/`NUM_PERTURBATIONS`, so tuning those constants shifts the deterministic RNG ordering for every subsequent effect in a session. (`effects/BZReactionDiffusion.h` 246–253)

---

## 4. Findings Refuted on Verification

For an honest record, these reported issues were **refuted** by independent re-examination and are
**not** defects:

- *`random_vector` rejection loop can hang on pathological/NaN `rand_f`* — REFUTED: `hs::rand_f()` is a 32-bit PRNG mapped to `[0,1)` and cannot be NaN; termination is geometrically certain (`core/geometry.h` ~788).
- *`ease_in_out_cubic` is missing an endpoint clamp* — REFUTED: all easings are intentionally unclamped (documented), and the cubic is a pure polynomial exact at both ends with no NaN/overflow term needing a guard (`core/easing.h` 26–34).
- *`buffer_free` watchdog can spuriously trap on a 64-bit host due to mixed-width `micros()`* — REFUTED: the host `micros()` returns a full-width `unsigned long` matching `wait_start`; no mixed-width subtraction exists (`core/canvas.h` 647, `core/platform.h` 1010).
- *`deep_tween` drops the newest head sample* — REFUTED: the skipped trailing motionless sub-frame is the shared junction the prior frame already plotted; the head sample of a moving trail is always emitted (`core/animation.h` 2975–3002).
- *`SplineFlow` passes an unclamped negative arg to `quintic_kernel` producing garbage alpha* — REFUTED: `quintic_kernel` clamps its input to `[0,1]` internally (`core/3dmath.h` 80).
- *`Star` ctor can divide by a zero edge-normal length* — REFUTED: with the fixed inner ratio and `radius>0` guard, the edge length is provably nonzero (`core/sdf.h` ~2474).

Two further items were **PARTIAL**: the `make_rotation` tiny-angle drop and the `Persist`/watermark
limitation are *real but already self-documented as intentional tradeoffs* (carried above as #27 and
#39 at Low), and the `gyro→kis` arena-polarity hazard is real but — contrary to the original claim —
**is** exercised, because `truncatedOctahedron_gyro_kis_hk17` is in `islamic_registry` and the
Islamic-recipe budget sweep covers every entry.
