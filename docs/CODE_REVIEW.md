# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ engine/firmware (`core/`, `effects/`, `hardware/`, `targets/wasm/`, `tests/`, `tools/`) and the daydream web simulator (JS modules, `tools/`, `tests/`).
**Out of scope (excluded by request):** `core/effects_legacy.h`, `core/rotate.h`, `targets/*/*.ino`, and vendored third-party (`FastNoiseLite.h`, `inplace_function.h`, `three.js`, generated `holosphere_wasm.js`).

**Method:** every in-scope file was read by component-focused review agents consulting the README for architectural intent; every higher-severity finding was then independently re-verified against source by a separate validation pass. Findings that did not survive validation were dropped (notably an initially-reported P1 "H_OFFSET slider crash," which was **refuted** — the live code calls `precomputeMatrices()`, not `setupDots()`, so no crash exists).

**Bottom line:** This is a mature, unusually disciplined codebase. No P0 (critical) or P1 (major-correctness) defects were confirmed anywhere in scope. The engine's memory model (compile-time-sized arenas, fail-fast `HS_CHECK`, host/device parity testing), its rendering pipeline (compile-time filter ordering, seam/pole/antipode handling), and its build-gate tooling are engineered to a standard well above typical hobbyist or even much professional embedded-graphics code. The confirmed issues are edge-case robustness gaps, leak-on-misuse paths, documentation/code drift, and tooling-test coverage holes — not steady-state logic faults.

---

## Letter Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | **A−** | No confirmed P0/P1. Core math, rasterization, simulation, and sync logic are carefully reasoned with explicit edge-case handling. The P2s are narrow: OOB only on a malformed mesh, leaks only on misuse/double-init, a worker-barrier gap only on a duplicate message. |
| **Memory / Resource Safety** | **A−** | Exemplary arena discipline (owned `ArenaVector` vs borrowed `ArenaSpan` at the type level, wrap-proof bounds math, generation stamps catching use-after-reset). Dragged from A by leak-on-misuse paths (ParticleSystem double-init, MeshMorph transients, repeating-transformer slot leak), one device-stripped per-fragment bounds `assert`, and a missing JS recorder teardown. |
| **Performance** | **A** | Strong embedded discipline: LUT-driven trig, compile-time resolution, zero per-frame heap traffic, premultiplied SSAA, baked palettes, on-demand JS rendering with instanced draws. Honest about the few heavy paths (SphericalHarmonics). |
| **Architecture / Design Elegance** | **A** | X-macro single-source-of-truth rosters, CRTP iterators/animation base, compile-time filter-ordering `static_assert`s, the pure/host-testable split under device shells (`pov_sync`, `wasm_predicates`, `param_marshal`, JS `sidebar_logic`/`param_sync`). Genuinely well-factored. |
| **Interface Expressiveness** | **A−** | `explicit` conversions, `[[nodiscard]]`, deleted rvalue-borrow overloads encoding lifetime in the type system, typed worker protocol under `@ts-check`. Minor footguns (greedy variadic ctor, long positional draw-param lists). |
| **Readability / Maintainability** | **A−** | Documentation explains *why*, not *what* — among the best-commented code reviewed. Pulled from A by extreme comment density (drift risk: doc blocks longer than the methods they describe) and several knowingly-duplicated blocks (trail rendering across 3 effects, manifold checks). |
| **Error Handling** | **A−** | Consistent two-tier doctrine: untrusted boundaries (JS→WASM) reject-and-log; internal invariants fail-fast trap. Pulled down by a device-stripped `assert`, a couple of silent JS skip paths (recorder without `requestFrame`), and the warning-ratchet misclassification. |
| **Concurrency / ISR Safety** | **A−** | The single-core Cortex-M7 model is correctly reasoned (relaxed atomics, one precise release/acquire edge, dcache flush as the real coherence mechanism); JS workers use a generation fence + fault latch. One fragile dangling-pointer window and one JS barrier-count gap keep it from A. |
| **Testing** | **A−** | C++ suite is excellent (A): 34 modules, real property/oracle assertions, out-of-process death tests, host/device parity TUs, X-macro rosters that can't drift. JS suite passes 235/0. Held back by untested build-gate paths (`--size-a` fallback, JSONC stripper, `teensy_gate_extra.py`), untested JS `slider.js`/importmap generator, and a handful of weak oracles. |
| **Portability** | **A−** | Host/device parity is explicitly tested (device specialization re-run under asserts, `-ffast-math` clamp TU, PROGMEM handled correctly for Teensy memory-mapped flash). Minor: a few latent int-width `static_assert`s missing for future resolution bumps. |
| **Documentation** | **A** | README is comprehensive and accurate to the code; per-symbol doxygen documents load-bearing invariants. The only knock is volume (maintenance burden), not accuracy. |
| **Overall** | **A−** | A standout embedded-graphics project. Clean of serious defects; the work remaining is polish, robustness-on-misuse hardening, doc/code reconciliation, and closing tooling-test gaps. |

---

## Prioritized Fix List

Items are numbered sequentially. Severity: **P0** critical, **P1** major correctness/safety, **P2** significant (robustness, leak-on-misuse, tooling integrity), **P3** minor / nit / documentation.

### P0 — Critical

_None found._

### P1 — Major correctness / safety

_None confirmed. (The initially-reported H_OFFSET render-loop crash was independently refuted: the live `onChange` calls `precomputeMatrices()`, so `instanceColor` is never null on the per-frame path.)_

### P2 — Significant (robustness, leak-on-misuse, tooling/test integrity)

1. ✅ **Plot mesh face-walk reads the flat index array out of bounds on a malformed mesh.** `Plot::Mesh::for_each_unique_edge` ([core/plot.h:1955](core/plot.h#L1955)) indexes `fi[offset + k]` without checking `offset + count <= mesh.get_faces_size()`, while the sibling `Scan::Mesh::draw` ([core/scan.h:727](core/scan.h#L727)) guards exactly this. Add the matching `HS_CHECK` so a face-counts/face-indices mismatch traps at the cold seam instead of reading past the array on-device.

2. ✅ **`Timeline` pin-ordering invariant is enforced only at removal time, not at registration.** `add_get` ([core/animation.h:2458](core/animation.h#L2458)) accepts a pinned infinite animation registered after a finite one; the `HS_CHECK(!handled)` trap fires later in `move_into` when the finite event completes — a confusing, data-dependent crash site. Validate at `add_get` (scan existing events for any finite one when `pin=true`) so misuse fails deterministically at the call site.

3. ✅ **`Timeline::step` event relocation rests on an unstated `write_idx <= active_cnt` invariant.** The compaction/relocation arithmetic ([core/animation.h:2574](core/animation.h#L2574)) is correct only because that invariant holds, but nothing asserts it. Add `HS_CHECK(write_idx <= active_cnt)` and a one-line comment pinning the no-gap case.

4. ✅ **`ParticleSystem::init` double-init orphans arena allocations.** `init` ([core/animation.h:553](core/animation.h#L553)) calls `pool.bind(arena, CAPACITY)` etc. with no re-init guard; a second `init()` leaks the first (large) allocation since the arena has no per-allocation free. Add `HS_CHECK(!pool.is_bound(), ...)` (the accessor already exists and is used in `spawn`).

5. **Repeating non-pinned transformer spawns leak their pool slot, and the guard's rationale comment is wrong.** The leak guard ([core/transformers.h:157](core/transformers.h#L157)) is `HS_CHECK(p->is_finite() || p->repeats())`, which *passes* for a repeating animation; but the recycling `.then()` frees only `if (!repeats())`, so a repeating (e.g. `Driver`/`Noise`, `AnimationBase(1,true)`) non-pinned spawn never frees its slot. The comment treats repeating spawns as safe — they are not. Require repeating spawns to be pinned+cancellable, or free the slot on cancel, and correct the rationale.

6. ✅ **`MindSplatter` per-fragment pool-index bounds check is a device-stripped `assert`.** [effects/MindSplatter.h:272](effects/MindSplatter.h#L272) uses bare `assert(p_idx < active_count)` then `pool[p_idx]` in the per-fragment shader; under the optimized device build (`NDEBUG`) the check vanishes while the float-derived index still feeds `pool[]`. (It is the only bare `assert` in `effects/`.) Move the guarantee to a cold path where fragments are produced, or clamp `p_idx` in the hot path so a bad index cannot read OOB on hardware.

7. ✅ **`pov_segmented` flywheel timer has no minimum-period guard (it is even compile-time checkable).** [hardware/pov_segmented.h:237](hardware/pov_segmented.h#L237) passes `kColumnUs / kOversample` straight into `IntervalTimer::begin` with no check, while [hardware/pov_single.h:133](hardware/pov_single.h#L133) guards `interval_us >= 1`. Since `kColumnUs`/`kOversample` are `constexpr`, add `static_assert(kColumnUs / kOversample >= 1.0f, ...)` next to the other geometry asserts.

8. ✅ **`pov_segmented` leaves `pending_effect_` holding a freed pointer for a window after `delete cur`.** Foreground deletes `cur` ([hardware/pov_segmented.h:262](hardware/pov_segmented.h#L262)) and only republishes `pending_effect_` at line 276; in between, the stale pointer is observable to the ISR and is currently saved *only* by generation/handshake guards, never by nulling. Null `pending_effect_` (release store) before the `delete` so the generation guards become defense-in-depth rather than the sole protection.

9. ✅ **Warning ratchet misclassifies vendored third-party paths as first-party.** `_relativize` ([tools/teensy_warnings.py:39](tools/teensy_warnings.py#L39)) does `p.find("/" + fp)`, so `…/.platformio/lib/SomeLib/effects/reverb.h` normalizes to `effects/reverb.h` and is counted as first-party — letting library warnings pollute (or be baked into) the baseline. Anchor first-party matching to the repo root and reject paths passing through `lib/`/`.platformio/`/`packages/`. (The existing `test_library_warning_excluded` fixture has no first-party segment, so it never exercises this.)

10. ✅ **Warning-ratchet path normalization aliases distinct nested files to one key.** Same function: `/…/Holosphere/targets/Phantasm/effects/Foo.h` and a top-level `effects/Foo.h` both collapse to `effects/Foo.h` (first-occurrence `find`), so a new warning in one can be masked by a baseline entry from the other. Compute the path relative to the project root (longest first-party prefix), not the first matching segment.

11. ✅ **Build-gate test suite leaves its riskiest paths uncovered.** `test_teensy_gate.py` does not exercise (a) the `--size-a` fallback integration through `evaluate` (including the `0x80000 - ram1` free-headroom arithmetic), (b) `_strip_jsonc_comments` (the bespoke parser that guards the single-source-of-truth budgets file and explicitly replaced a buggy regex), or (c) `tools/teensy_gate_extra.py` (toolchain discovery + `exit(2)` guards). Add unit tests for each — these are exactly the fail-closed-critical and fallback-arithmetic paths. _(The core gate is genuinely fail-closed on malformed size/symbol input — verified — so this is a coverage gap, not a live fail-open.)_

12. ✅ **Segment controller's frame barrier double-decrements on a duplicate `frame` message.** The handler ([segment_controller.js:245](segment_controller.js#L245)) does `this.pending--` unconditionally per message with no per-segId received flag, so a worker emitting two `frame` messages in one generation resolves the barrier early and composites a torn/mixed-generation frame. Gate the decrement behind a first-arrival flag reset each generation so the barrier counts distinct segments.

13. ✅ **`addColor` deep links can be self-rejecting on reload.** The writer ([gui.js:309](gui.js#L309)) only `#`-prefixes a bare-6-hex string and passes other string forms verbatim, but `isValidColorString` ([gui.js:33](gui.js#L33)) accepts only `#hex`/`rgb()`. A color value serialized as a named color or bare 3-hex is written to the URL and then rejected on reload, silently reverting to default. Keep writer and reader in lockstep (normalize any emittable form to one the reader accepts).

14. ✅ **Video recorder buffers the entire recording in the encoder with no bound.** `recorder.start()` ([recorder.js:159](recorder.js#L159)) is called with no timeslice, so `ondataavailable` fires once at stop — a long recording grows encoder memory unbounded (~2 MB/s at 16 Mbps) with no max-duration. Pass a timeslice for incremental delivery and/or add a soft max-duration auto-stop.

15. ✅ **`VideoRecorder` has no `dispose()`, so mid-recording teardown leaks the stream/track.** Cleanup runs only via `onstop` ([recorder.js:152](recorder.js#L152)); there is no `dispose()`, unlike `URLSync.dispose` ([state.js:187](state.js#L187)) and `EffectSidebar.dispose` ([sidebar.js:81](sidebar.js#L81)). Add a `dispose()` that stops an active recording (tracks + offscreen) and wire it into the same `pagehide`/`disposeApp` teardown.

### P3 — Minor / nits / documentation

**core math & geometry**

16. `normalize()`/`normalized()` gate on `length() >= epsilon` (relative spacing near 1.0), admitting sub-micro-length vectors that yield ~1e6-magnitude "unit" results ([core/3dmath.h:262](core/3dmath.h#L262)). Gate on squared length against a named small constant, consistent with `make_rotation`/`make_basis`.
17. `vector_to_pixel` documents the `y` south-pole overshoot but not the symmetric `x → W` rounding hazard ([core/geometry.h:489](core/geometry.h#L489)); extend the `@return` floor-don't-round contract to `x`.
18. `Orientation::upsample` does a redundant `slerp(q, q, 0)` at the `t==1` endpoint and resamples linearly in source-index (not arc-length) space ([core/geometry.h:727](core/geometry.h#L727)); special-case the endpoint and tighten the "shape-preserving" doc.
19. `fast_cbrt` silently collapses to 0 in the denormal/tiny-normal tail ([core/3dmath.h:359](core/3dmath.h#L359)); promote the "supported for x ≳ 1e-6" note into the `@param` contract.
20. `make_basis` hand-rolls the antiparallel reference-axis pick that `least_parallel_axis` was created to centralize ([core/3dmath.h:1045](core/3dmath.h#L1045) vs [core/geometry.h:813](core/geometry.h#L813)); route through the shared helper or cross-reference why it can't.
21. `wrap(T,U)` silently promotes large/`size_t` integral pairs to `double` fmod (only the exact `wrap(int,int)` overload is precise) ([core/util.h:46](core/util.h#L46)); add an integral overload or `static_assert` so a lossy wrap fails to compile.
22. `fib_spiral`/`lissajous`/`logPolarToVector` use exact trig while per-pixel paths use fast trig, with no documented rationale ([core/geometry.h:544](core/geometry.h#L544)); add a "setup-time generator; exact trig intentional" note.
23. `gnomonic` drops the hemisphere sign at the equator and `inv_gnomonic` needs it out-of-band ([core/3dmath.h:766](core/3dmath.h#L766)); document at the `gnomonic` site that the caller must track the sign.
24. `StaticCircularBuffer`'s greedy variadic forwarding ctor is a future-maintenance footgun ([core/static_circular_buffer.h:74](core/static_circular_buffer.h#L74)); consider `explicit` and note the hazard.
25. `apply_if_changed` uses `operator!=`, which for `Vector` is the non-transitive tolerance comparator ([core/util.h:149](core/util.h#L149)); document/assert exact-equality `T` (scalar/int).
26. `fast_sinf` has no domain note; accuracy degrades for large args (only `STEREO_PATTERN_ARG_LIMIT` documents it) ([core/3dmath.h:1165](core/3dmath.h#L1165)); add a `@warning` and cross-reference.

**core rendering / raster**

27. ✅ `World::Replicate::plot` rotates each copy from the previous one, compounding float drift so copies aren't evenly spaced and the loop doesn't exactly close ([core/filter.h:590](core/filter.h#L590)); precompute per-copy quaternions or renormalize each step (cf. `VertexReplicate`).
28. ✅ README "1-pixel AA border" (§7.1) doesn't match the implemented 2-pixel coverage ramp ([core/scan.h:62](core/scan.h#L62)); reconcile the wording.
29. ✅ `Subtract` non-solid-B pass-through emits un-normalized A spans, relying on `scan_region` to re-normalize ([core/sdf.h:1141](core/sdf.h#L1141)); add a comment that seam handling is deliberately deferred (the solid-B branch documents its normalization but this one is silent).
30. ✅ `Pipeline::plot` float overload enforces the `[-W, 2W)` column bound via `assert` *after* the cast, and `fast_wrap` corrects only a single ±W offset ([core/filter.h:158](core/filter.h#L158)); note that the single-step correction is the load-bearing reason the contract is `[-W,2W)`.
31. ✅ `BoundingSphere::get_intervals` near-pole fold uses a heuristic `angular_radius/PI_F` threshold that force-scans full rows for large radii ([core/scan.h:291](core/scan.h#L291)); correct and conservative, but a future tightening should re-derive the constant.

**core color**

32. ✅ `linear_to_srgb_lut` is a 64 KB direct table ([core/color_luts.h:35](core/color_luts.h#L35)); a 4096-entry table + interpolation (the idiom already used for the forward direction) would cut it ~16× at imperceptible 8-bit-output error. Note as a deliberate flash/speed tradeoff in the RAM/flash audit if kept.
33. ✅ `srgb_to_linear_interp` interpolates in linear space between sRGB-spaced codes, a small upward (secant) bias ([core/color.h:345](core/color.h#L345)); soften the "keeps sub-8-bit precision" comment or compute exactly on this cold path.
34. ✅ `GenerativePalette::get` zeroes a stop's chroma when its lightness nears 0/1 ([core/color.h:1227](core/color.h#L1227)); correct for current callers but document the envelope-model assumption so a future authoring path doesn't silently desaturate.
35. ✅ `oklch_to_cpixel` open-codes the linear→8-bit-sRGB quantization three times ([core/color.h:846](core/color.h#L846)); extract `linear_float_to_srgb8()` paralleling `float_to_pixel16`.
36. ✅ `ProceduralPalette` toggles `protected`→`public` solely to place a trivial destructor ([core/color.h:1453](core/color.h#L1453)); move the `= default` destructor up and keep one trailing `protected:`.

**core mesh**

37. `kis`/`expand`/`snub` build a centroid apex and call `normalize()`, which traps on a zero-length (centrally-symmetric face) centroid, whereas `dual` guards the same case with `normalized_or` ([core/conway.h:482](core/conway.h#L482)); use `normalized_or` consistently rather than relying on the unenforced "no origin-centered vertex" precondition.
38. `MeshOps::transform` produces borrowed-mode output but rejects borrowed-mode *input*, so `transform(transform(x))` traps ([core/conway.h:251](core/conway.h#L251)); route the topology copy through the unified `get_*_data()` accessors or document why chaining is unsupported.
39. `relax` omits the `require_closed_manifold` precondition its sibling operators all carry, silently giving a partial relaxation on a boundary mesh instead of failing fast ([core/conway.h:936](core/conway.h#L936)); add the check or document the tolerance.
40. `compile_hankin` open-codes a closed-manifold scan duplicating `require_closed_manifold` ([core/hankin.h:130](core/hankin.h#L130)); centralize (move the helper to `mesh.h`).
41. `find_nearest_node` calls the double-precision `node()` (cos/sin/fmod/sqrt) inside the cubemap-LUT build hill-climb — millions of transcendentals at boot ([core/reaction_graph.h:214](core/reaction_graph.h#L214)); precompute the lattice once into scratch to cut boot latency. (One-time cost; never hits frame time.)

**core engine / animation / platform**

42. `MeshMorph` allocates clone+position transients into the arena with no destructor reclamation, so back-to-back morphs without an explicit compaction grow the arena ([core/animation.h:1955](core/animation.h#L1955)); documented intentional, but the required compaction cadence should be flagged at the call sites.
43. `random_to_unit`'s clamp constant is correct only for `max == UINT32_MAX`, but `max` is a runtime parameter (the `static_assert` lives in `rand_f`) ([core/platform.h:1111](core/platform.h#L1111)); document the precondition or derive the boundary from `max`.
44. `Presets::get`/`prev_get` re-`HS_CHECK` an index the class invariants already prove in `[0,Size)`, while the index accessors don't — an asymmetry signaling dead checks or missing ones ([core/presets.h:94](core/presets.h#L94)); pick one.
45. `effect_registry.h` `kReserveHint = 64 // ~2x roster` is decoupled from the actual 27-effect roster ([core/effect_registry.h:73](core/effect_registry.h#L73)); drop the misleading ratio claim (the file can't see `HS_EFFECT_COUNT`).
46. `Sprite` fade-out casts `uint32_t t` to `int` ([core/animation.h:1208](core/animation.h#L1208)) — a latent UB seam reachable only for absurd finite durations; compute the comparison in the unsigned domain like the fade-in branch.
47. `led.h` correction-guard depth counter is a non-atomic `int` with no reentrancy note, unlike `hs::random()` ([core/led.h:60](core/led.h#L60)); add the same "main-loop only" one-liner.
48. `beat16` computes `bpm << 8` on a `uint16_t` (promotes to `int`, then narrows) — intentional FastLED parity but a `-Wconversion` trip ([core/platform.h:802](core/platform.h#L802)); make the truncation explicit.

**core memory / wasm**

49. `Arena::allocate` OOM trap logs `offset + padding`, which can exceed `capacity` and reads as a nonsensical offset ([core/memory.h:94](core/memory.h#L94)); log `padding`/`size` separately.
50. `drawFrame` pixel-index accumulators are `int` with no `static_assert(MAX_W*MAX_H*kChannels <= INT_MAX)` co-located ([targets/wasm/wasm.cpp:435](targets/wasm/wasm.cpp#L435)); add the guard for a future resolution bump.
51. An `ArenaSpan` cannot track its source vector across a *move* in debug staleness checks (runtime-safe, since it snapshots `data_`) ([core/memory.h:705](core/memory.h#L705)); document that a span must be re-taken after the source is moved.

**effects**

52. BZ and GS normalize their diffusion kernels at different points (raw bytes + folded `Q8_INV` vs pre-normalized) ([effects/BZReactionDiffusion.h:330](effects/BZReactionDiffusion.h#L330)); unify for the shared-base philosophy.
53. ✅ BZ's Lotka-Volterra reaction term is outside the diffusion-only stability bound the comment states and saturates (hard banding) at slider extremes ([effects/BZReactionDiffusion.h:208](effects/BZReactionDiffusion.h#L208)); narrow the slider tops or document that the reaction term intentionally saturates.
54. ✅ `FlowField` mixes scaled position with unscaled time/offset terms into the noise z-axis, coupling three visual axes onto `Scale` ([effects/FlowField.h:85](effects/FlowField.h#L85)); scale consistently or document the intent. The `f.v0` alpha source is also undocumented ([effects/FlowField.h:140](effects/FlowField.h#L140)).
55. `Voronoi` coarse-grid fast path can drop a cell smaller than `kCoherenceBlock` near the poles ([effects/Voronoi.h:160](effects/Voronoi.h#L160)); documented speed/accuracy tradeoff — make `kCoherenceBlock` latitude-adaptive if pole artifacts appear.
56. `SphericalHarmonics` evaluates the Legendre recurrence twice per sub-sample during a morph with no horizontal narrowing — the heaviest per-pixel path, violating the rasterizer-bound assumption ([effects/SphericalHarmonics.h:204](effects/SphericalHarmonics.h#L204)); already capped at `l≤4` and documented, but profile first on-device.
57. `HopfFibration` recomputes a constant per-fiber `phase` every frame ([effects/HopfFibration.h:227](effects/HopfFibration.h#L227)); cache at seed time, consistent with its own stated philosophy.
58. ✅ `Moire` sets a resolution-dependent `density` default before `registerParam`, a load-bearing ordering that is implicit ([effects/Moire.h:60](effects/Moire.h#L60)); add an ordering comment.
59. ✅ `RingSpin::spawn_ring` carries a stale `@param normal` doxygen line for a parameter it doesn't take (it always seeds `Y_AXIS`) ([effects/RingSpin.h:158](effects/RingSpin.h#L158)); drop the stale doc (or add the parameter to restore the documented behavior). _(Validation note: the `Ring(normal, palette)` ctor is **not** dead — it's the one called, just always with `Y_AXIS`.)_
60. Trail record/`deep_tween`/`quintic_kernel`-fade skeleton is duplicated across Comets, RingSpin, and ChaoticStrings with a self-documented "propagate fixes by hand" contract ([effects/Comets.h:131](effects/Comets.h#L131)); extract a shared helper parameterized by the per-fragment draw callback, or add a shared invariant test.
61. ✅ `Thrusters` palette banding samples `angle_between(X_AXIS, v)` in world space, contradicting its "rotation preserves angles" comment ([effects/Thrusters.h:271](effects/Thrusters.h#L271)); fix the code or the comment to match the intended (world- vs ring-anchored) banding.
62. `Dynamo::color` can pick a stale palette for a few frames after a live Wipe-Dur edit ([effects/Dynamo.h:208](effects/Dynamo.h#L208)); documented self-healing cosmetic — optionally gate the early-return on a monotonicity flag.

**hardware**

63. ✅ A stably mis-wired/floating ID strap can silently elect a phantom second master driving the push-pull sync line into contention ([hardware/pov_segmented.h:353](hardware/pov_segmented.h#L353)); out of firmware scope, but document that ID0 must be a positively-driven strap.
64. ✅ `SymbolEmitter::schedule_beacon`'s int32 signed "is-due" comparison is correct only while the worst-case burst span stays well under 2³¹ cycles, enforced by prose not by a check ([hardware/pov_sync.h:975](hardware/pov_sync.h#L975)); add a `Config::valid()` clause bounding it.
65. ✅ `position()`'s int32 cast safety relies on the fold-before-position ordering in `tick()`, argued in a comment rather than asserted ([hardware/pov_sync.h:598](hardware/pov_sync.h#L598)); assert the coast bound or document that the fold loop is the invariant.
66. ✅ `HD107SFrame::load()` does a partial write that does not clear the tail `[count, N)` ([hardware/hd107s_frame.h:157](hardware/hd107s_frame.h#L157)); a surprising public-API edge — zero the tail or rename/assert full-frame.
67. ✅ `factor()` maps intermediate brightness as `(f+1)/256`, differing slightly from the documented `f/255` scale at all non-endpoint values ([hardware/hd107s_frame.h:297](hardware/hd107s_frame.h#L297)); document the exact intermediate mapping, not just the exact endpoints.
68. ✅ Comment density creates real drift risk — 35-line rationale blocks over 8-line methods, with the single-core/relaxed-atomic rationale re-derived 4–5× ([hardware/dma_led.h:126](hardware/dma_led.h#L126), [hardware/pov_sync.h:828](hardware/pov_sync.h#L828)); consolidate the repeated rationale into one referenced note.

**tools (Python gates)**

69. ✅ `parse_readelf_sections` column-shifts on the empty-Name NULL section row, mislabeling section 0 ([tools/teensy_gate.py:197](tools/teensy_gate.py#L197)); affects only the human-readable failure message, but parse by fixed columns or skip ndx 0.
70. ✅ `region_totals_from_size_a` buckets a section solely by start VMA, ignoring a size spill past a region boundary (and can print a negative "free") ([tools/teensy_gate.py:146](tools/teensy_gate.py#L146)); it's the documented cross-check path, so note the spill limitation in the docstring.
71. ✅ `_find_teensy_size` selects any same-named binary whose `--help` merely launches (non-zero exit not rejected) ([tools/teensy_gate_extra.py:46](tools/teensy_gate_extra.py#L46)); the downstream empty-parse `exit(2)` backstops it, but validate the probe output.

**tests (C++ native suite)**

72. ✅ Circular-buffer wrap test asserts only `size() <= 3` (vacuous — capacity *is* 3) and `back()`, never the surviving sequence/head-tail positions ([tests/test_static_circular_buffer.h:344](tests/test_static_circular_buffer.h#L344)); assert the exact live window.
73. ✅ HD107S wire-order test asserts lit channels `> 0` rather than exact corrected bytes ([tests/test_hd107s_frame.h:199](tests/test_hd107s_frame.h#L199)); mitigated by `test_correct_multifactor`, but assert exact values or note the delegation.
74. ✅ Easing monotonicity check uses `v >= prev - 1e-4f`, letting steady cumulative drift slip ([tests/test_easing_waves.h:51](tests/test_easing_waves.h#L51)); tighten to `>= prev` (exact-arithmetic curves) and pin both endpoints.
75. ✅ Reaction-graph determinism test compares `node(1234)` to itself — always true ([tests/test_reaction_graph.h:75](tests/test_reaction_graph.h#L75)); compare to a frozen golden or drop it.
76. Platform LUT-error ceilings (`sin8 < 6`, `sin16 < 300`) sit well above measured worst-case ([tests/test_platform.h:57](tests/test_platform.h#L57)); tighten to just above the measured peaks.
77. Melt-warp style test's `out.x < 1.0f` leg is near-tautological ([tests/test_styles.h:139](tests/test_styles.h#L139)); assert a minimum drift magnitude.
78. The absent-filter negative test is behind a default-off macro neither preset nor CI defines ([tests/test_filter.h:213](tests/test_filter.h#L213)); convert to a runtime/death assertion or wire the macro into CI.
79. The in-process double-buffer spin-wait test is the suite's one latent timing flake ([tests/test_canvas.h:253](tests/test_canvas.h#L253)); prefer a deterministic cooperative-yield model or a generous watchdog.

**daydream driver / render**

80. `isViewLive` reads `view.buffer.byteLength` with no guard for a non-typed-array argument ([pixel_view.js:25](pixel_view.js#L25)); since this is the single source of truth for the buffer-alias contract, guard `view.buffer`/`ArrayBuffer.isView`.
81. The `testAll` interval keeps firing for the page lifetime if the WASM load failed (each tick harmlessly early-returns) ([daydream.js:600](daydream.js#L600)); disable the control in the load `.catch`.
82. `_advanceFrameClock` consumes at most one interval per frame, so the `MAX_FRAME_CATCHUP` clamp + "catch up by a few frames" comment describe behavior the single-step structure doesn't implement ([driver.js:440](driver.js#L440)); loop the step for true fixed-timestep catch-up or simplify the comment.
83. While paused, `_advanceFrameClock` still drains the accumulator, so the first post-unpause frame can briefly stall ([driver.js:396](driver.js#L396)); skip consumption while paused or reset on unpause.

**daydream segment workers**

84. Odd canvas width gives the two arms unequal widths, diverging from the firmware's symmetric `w/2` split ([segment_layout.js:59](segment_layout.js#L59)); both shipped presets are even, so document/assert even-`w` or seam at exactly `floor(w/2)`.
85. `composite()` validates bounds inside the blit loop, so a valid segment ahead of a bad one is already composited when the fault latches (one partial frame) ([segment_controller.js:513](segment_controller.js#L513)); validate all results in a pre-pass.
86. A faulted worker pool only recovers on a user-initiated resolution/mode change — no bounded auto-restart for transient errors ([segment_controller.js:372](segment_controller.js#L372)); deliberate fail-fast, but consider a bounded auto-recreate or document the choice.
87. The `frame` handler captures `paramValues` before the generation fence, so a stale-generation segment-0 frame can momentarily publish param values against a new descriptor list ([segment_controller.js:230](segment_controller.js#L230)); move the capture inside the current-generation block.

**daydream GUI / state / recorder**

88. `_getKey` drops a level when an intermediate folder has no `folderName`, risking a deep-link key collision ([gui.js:144](gui.js#L144)); assert non-empty folder names or include a stable fallback segment.
89. `captureFrame` silently no-ops (one `console.warn`) on browsers without `track.requestFrame`, yielding an empty video while recording "succeeds" ([recorder.js:181](recorder.js#L181)); probe at `start()` and refuse or fall back to timed `captureStream(fps)`.
90. `selectMimeType` probes only bare `video/mp4;codecs=avc1`, which some engines report unsupported, silently downgrading to WebM ([recorder.js:33](recorder.js#L33)); add fully-qualified fallbacks (`avc1.42E01E`, `video/mp4`).
91. Object-enum URL hydration can pick a display label that doesn't match the link when two options share a value ([gui.js:236](gui.js#L236)); document distinct-value enums or match on label.
92. The standalone-tool-page URL writer's fallback `setTimeout` is not cleared by `DeepLinkGUI.destroy()` ([gui.js:56](gui.js#L56)) — the exact hazard `URLSync.dispose()` exists to prevent; expose `cancel()`/`flush()` and call it from `destroy()`.
93. URL-hydration replay fires only the first registered `onChange` handler via a single `replayed` latch ([gui.js:166](gui.js#L166)); a second fan-out consumer's load-time side effect is skipped — replay per-handler.

**daydream tools / JS tests**

94. `relax` codegen defaults to 100 iterations vs the C++ `SolidBuilder` default of 8 ([tools/solid_codegen.js:135](tools/solid_codegen.js#L135)); the emitted C++ is explicit/valid, but align the fallback to 8 (or document the divergence).
95. Procedural-palette export rounds coefficients to 3 decimals, a separate (undocumented) fidelity gap from the linearization gap ([tools/palette_math.js:400](tools/palette_math.js#L400)); widen to 6 digits or document.
96. `tools/slider.js` validation branches (NaN guards, step/scale rounding) are untested ([tools/slider.js](tools/slider.js)); add a small `slider.test.js` with a DOM stub.
97. `scripts/generate-importmap.mjs` GENERATED-block rewrite has no test for its output (only a marker-drift `exit(1)`) ([scripts/generate-importmap.mjs](scripts/generate-importmap.mjs)); add a fixture test for `--local` vs default modes.
98. `copyToClipboard`'s async-API `catch (err) {}` is empty with an unused binding ([tools/clipboard.js:30](tools/clipboard.js#L30)); use `catch {}` or add a `console.debug` for the double-failure path.

---

## Notable Strengths

- **Fail-fast doctrine, applied with discipline.** `HS_CHECK` is genuinely always-on (survives `NDEBUG`, no heap, no stdio, ends in `__builtin_trap`), placed at cold seams and deliberately withheld from per-pixel hot paths — and the one documented `assert`-not-`HS_CHECK` exception is intentional. Verified fail-closed in the build gates: malformed size/symbol input produces violations, never a silent green.
- **A memory model that encodes safety in the type system.** Owned `ArenaVector` vs borrowed `ArenaSpan`, wrap-proof subtractive bounds math, generation/rebind debug stamps catching use-after-reset, compile-time `static_assert`s tying buffer capacities to worst-case geometry. Worst-case Conway-operator capacities were independently verified tight.
- **Compile-time correctness.** X-macro single-source-of-truth rosters (resolutions, effects) cross-checked at startup; filter-pipeline ordering constraints enforced as `static_assert`s; the DMAMEM-drop-on-template-statics hazard correctly defended via explicit specialization.
- **Seam / pole / antipode handling** across the rasterizer, SDF CSG, and geodesic curves is shared bit-identically and reasoned through closed-form edge cases — a class of bug most renderers get wrong.
- **Host/device parity is tested, not assumed** — device specializations re-run under asserts, a `-ffast-math` clamp TU, integer-only determinism contracts, and a pure/host-testable split that pulls the untestable peripheral shell down to a thin edge.
- **The web simulator's WASM-detach handling** (the README's central hazard) is textbook: liveness checks, alias-refresh on growth, per-frame divergence self-heal, and nulling GPU-aliased buffers before disposal. JS tests pass 235/0.

_No P0 or P1 defects were confirmed in any reviewed component._
