# Holosphere + daydream — Code Quality Review

**Date:** 2026-07-26
**Scope:** the `Holosphere` repo (C++ engine, effects, firmware, tests, tooling, build) and the
`daydream` repo (web simulator, tools, JS tests). Reviewed by 28 independent scoped audits, one per
component, each required to validate a finding against the cited code before reporting it.

**Out of scope (by instruction):** `core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/math/rotate.h`. Vendored third-party code (`FastNoiseLite.h`, `three.js`, `lil-gui`) was
reviewed only at its integration seams.

**Measured scale:** ~54.9k lines of hand-written C++ (engine + effects + hardware + targets),
~36.8k lines of generated tables, ~49.2k lines of native tests, ~6.0k lines of Python/Node tooling,
~13.6k lines of daydream JavaScript. 4,575 commits since 2015-07-27.

---

## 1. Summary

This is an exceptionally well-engineered codebase. The engine's core abstractions — the compile-time
filter pipeline with automatic domain lifting, the three-arena memory model with triple staleness
stamping, the SDF/rasterizer split, the Conway operator scaffolds, the pure/platform split that makes
firmware logic host-testable — are at or above the standard of professional graphics middleware. The
comment discipline is the best I have audited: comments record *measurements and rejected
alternatives*, not intentions.

The weaknesses are concentrated in three places, and they are systemic rather than scattered:

- **Gate integrity.** Four registered tests — including the relax-bake reproducibility gate — run in
  zero Linux/sanitizer CI jobs. Three separate gates (stack budget, warning ratchet, profile
  validate) return PASS when their input is empty or unparsed. A gate that goes green when its
  instrument breaks is worse than no gate, and this codebase has four of them.
- **README/docs drift.** The README is the project's primary artifact and is largely excellent, but
  it documents a deleted class (`MeshMorph`), the wrong effect count, two wrong parameter lists,
  a "stateless" filter that now carries a persistent cache, and a repo map missing ~16 files.
  Doxygen publishes *nothing* for the entire WASM bridge.
- **Rasterizer conservativeness.** Two scanline cull paths (`SDF::Ring`'s centerline fast path,
   `SDF::Flower`'s missing AA pad) are non-conservative in the visible direction, demonstrated
   numerically. The existing cull-conservativeness harness structurally cannot see either.

Nothing found is a data-loss or security-critical defect in the engine. The single highest-risk item
is a firmware phase relationship (finding 2) that is model-derived and must be confirmed on device
before action.

---

## 2. Letter grades

### 2.1 Correctness and robustness

| Dimension | Grade | Rationale |
|---|---|---|
| Engine correctness (render/mesh/color/animation) | **A−** | Across ~25k lines of the hardest code in the repo, the audits found no live logic error in the pipeline, the arena model, the Conway operators, or the color math. Defects are latent preconditions and two rasterizer cull gaps. |
| Firmware correctness (POV drivers, sync, DMA) | **B** | The 1-wire sync protocol, symbol codec, epoch scheduling and segment index math verify exactly against the datasheet at N=2/4/8. Offset against one High-severity render/display window phase relationship and two latent phase bugs. |
| Numerical robustness | **A** | Degeneracy is a first-class design axis: a named epsilon table where each constant states *what* it compares, `normalized()`/`normalized_or()` trap-vs-fallback pairs, antipodal slerp that sweeps a synthesized axis rather than blending through cancellation. |
| Color science | **A−** | Linear-light 16-bit blending, OKLCH shortest-arc hue, first-exit gamut bracketing (a genuinely sophisticated choice, documented with the pathology it avoids). Two generated LUTs verified bit-reproducible during this review. |
| Memory safety / resource discipline | **A** | Arena generation + per-vector rebind counter + `covers()` rewind check is a three-stamp staleness model where each stamp documents the hole it closes. Uniform subtractive overflow bounds. No heap on the render path. |
| Concurrency / ISR safety | **A−** | Single-writer model genuinely enforced; `EdgeMailbox::try_claim` fuses test-and-take so the race is unrepresentable rather than merely documented. Docked for one sticky flag across a multi-boundary tick. |
| Determinism | **A** | Bit-identical sim/device is an enforced project invariant: `hs::shuffle` only, seeded PRNG, generator/runtime double-precision lockstep, cross-run determinism sweeps for all 24 effects at both resolutions. |

### 2.2 Design and expressiveness

| Dimension | Grade | Rationale |
|---|---|---|
| Architectural elegance | **A** | The variadic `Pipeline` with compile-time domain lifting and four ordering `static_assert`s is the centerpiece and it earns that place. `Segue` composes without `MeshCarousel`; `param_marshal.h`/`pov_sync.h`/`pov_segment_map.h` make the platform boundaries host-testable. |
| API / interface expressiveness | **A−** | Lifetime contracts encoded in *types* (`StoredFunctionRef`'s deleted rvalue overloads, `ArenaSpan(const ArenaVector&&) = delete`), `DeepLinkGUI.add` vs `addSession`, `TickActions` as a pure "what the shell must do" value. Docked for positional-argument walls (`rasterize()`'s 11 params, the 10-field `PaletteHandoff`). |
| Modularity / separation of concerns | **A−** | Pure-logic extraction is consistently applied on both sides of the WASM boundary. `composition.h`'s non-self-contained include is the one structural wart. |
| Code duplication | **B+** | Mostly excellent factoring, with four real triplications: the scan seam-split walk (3×), the trail-gate prologue (2×, already drifted by one guard), the canvas-acquire IIFE (24×), and four near-parallel tool pages. |
| Complexity management | **B+** | `IslamicStars.h` (1520 lines), `animation/mesh.h` (2711), `sdf.h` (4589) are large but each is internally well-sectioned. `Feedback::flush` is a ~430-line monolith. |

### 2.3 Verification

| Dimension | Grade | Rationale |
|---|---|---|
| Test coverage breadth | **A−** | 0.90:1 test-to-source ratio. Every `core/`, `hardware/`, `targets/wasm` header maps to a test module. Structural gaps only: Teensy-register DMA and the embind surface. |
| Test quality / assertion strength | **A−** | Contract-level, not smoke: brute-force oracles for the KD-tree, a double-precision reference for gamut clipping, anti-identity assertions, bit-parity between fast and exact paths, a 4-board event-driven sync simulator with ppm skew and EMI injection. |
| Test harness architecture | **A−** | X-macro roster pinned in both directions, `end_module` counting "no assertions ran" as a failure, a real `fork`/`spawn` death harness with an anti-shrink floor. |
| **CI / gate integrity** | **C+** | Four registered tests match no shard; three gates PASS on empty input; the gamut LUT has no provenance gate; two profiling build configurations are compiled by nothing. The *design* anticipates all of this (`shard-coverage` exists precisely to catch it) — it simply is not green. |
| Determinism / flakiness | **B+** | All randomness explicitly seeded, no wall-clock in any CTest. Residual risk in absolute float tolerances and in the four modules that never run on Linux. |

### 2.4 Documentation

| Dimension | Grade | Rationale |
|---|---|---|
| In-code documentation (Doxygen/JSDoc) | **A−** | Genuinely load-bearing: bounds carry their proofs, comments cite measured ITCM byte deltas and name the libraries that diverged. Docked for ~12 concrete statements contradicted by their own code. |
| README accuracy | **C+** | The deep sections (sync datasheet, filter tables, WASM bridge) verify line-for-line. The drift is in counts, the repo map, two parameter lists, a deleted class, and one filter's state contract. |
| README structure / navigability | **A−** | All 30 TOC anchors resolve; machine-checked by `tools/docs_check.py` before every docs build. |
| `docs/` currency | **C+** | Excellent status-line discipline overall, but the largest spec (71 KB, cited by shipped code) is labelled "not implemented", and eight engineering docs are untracked. |
| Published API docs | **C** | `Doxyfile` does not define `__EMSCRIPTEN__`, so ~1500 lines of the WASM bridge — the entire JS-facing surface — publish as an empty page. |
| Onboarding | **B+** | §11 Building is genuinely runnable; a JS contributor is told nothing about daydream's 35-file test suite. |

### 2.5 Engineering practice

| Dimension | Grade | Rationale |
|---|---|---|
| Performance engineering | **A+** | The strongest axis. Selective-O3 regions justified by measurement, an ITCM byte ledger, a device-lock/profile-sweep regime, per-effect profile reports, dead levers recorded so they are not re-attempted. This is rigor rarely seen outside commercial engine teams. |
| Build system | **B+** | Near-exemplary pinning (every action at a full SHA, doxygen verified by sha256). Offset by the WASM TU compiling with zero warning flags and three uncompiled configurations. |
| Tooling quality | **B** | `tools/teensy_gate.py` and `tools/profile_tests/test_device_lock.py` are reference-quality negative-test suites. The rest of the surface has the silent-pass failures listed above. |
| Portability / platform hygiene | **B+** | Divergences catalogued in a dedicated ledger; two stdout generators emit CRLF on Windows, and one shell tool hardcodes an absolute path. |
| Security / supply chain | **B−** | Node deps exact-pinned with a 24-package lockfile and no `innerHTML` injection path anywhere; but CDN modules load with no Subresource Integrity while the WASM artifact gets sha256 provenance. |
| Repo hygiene | **B** | LICENSE/.gitignore/.gitattributes/.githooks/.clang-format all present in both repos. Docked for 8 untracked engineering docs, a stray build probe, and ~28 MB of 50×-oversized screenshots carried by both repos. |
| Web simulator quality | **B+** | The typed-view detach contract is airtight (the classic stale-`Uint16Array` bug is closed and tested via `ArrayBuffer.transfer()`); teardown is symmetric; two lifecycle escapes and one cross-language layout divergence. |
| Accessibility | **B** | Roving tabindex, listbox/option roles, `role=status`, `prefers-reduced-motion`. Sort controls unlabelled; no `<h1>` or landmarks. |

### 2.6 Overall

> ## **A−**
>
> An engine whose core design and performance engineering are at professional-middleware standard,
> carrying a documentation layer that has drifted from it and a CI gate layer with four holes that
> the project's own guards were built to catch.

---

## 3. Prioritized defect list

Every finding validated by the audits is listed. Numbering is sequential and unique across the whole
document. Findings in the `daydream` repo are prefixed **daydream:**.

### P0 — Critical: broken safety nets and wrong output

1. ✅ **Four registered tests run in zero Linux/sanitizer CI jobs** — `.github/workflows/ci.yml:42-48` vs `tests/CMakeLists.txt:106,189` — `unit_opchain_probe`, `unit_opchain_arena_survey`, `unit_spherical_field` and `unit_relax_bake_verify` match none of the four anchored shard regexes, so they execute only on the unsharded Windows job (no ASan/UBSan). The relax-bake gate — the proof that committed geometry bakes are still reproducible by their generator of record — is one of them. Fix: add `relax_bake_verify|spherical_field` to shard 2 and `opchain_probe|opchain_arena_survey` to shard 1. (Independently confirmed by two audits and by direct inspection.)

2. ✅ **Segment clip is sampled one display window early — every clipped frame may paint the wrong half** — `hardware/pov_segmented.h:320-323` — `clip_to_segment(cur, handoff_.window_left() == 0)` runs before `cur->draw_frame()`, but `Canvas`'s ctor blocks on `buffer_free()`, which the ISR only re-satisfies at the next flip; the sampled window is therefore two windows offset from the displayed one against a one-window alternation. Fixed under the zero-spill display-window contract by selecting the clip from a Canvas buffer-ready hook after the wait and before buffer advance/clearing, preserving the existing wait profiling. A composed host test over the real Effect double buffer, both Canvas constructors, and EffectHandoff reproduces 14 of 16 mismatches before the fix and 0 of 16 after it, completing the separately-owned acquisition-seam work. Render spills remain outside this parity contract and must be detected or handled separately.

3. ✅ **The stack-budget CI gate passes when its measurement fails** — `tests/stack_measure.cpp:83-95` — `peak` is initialized to 0 and assigned only inside the density loop; if no painted window is detected the row prints `peak = 0 B`, `worst > BUDGET_BYTES` is false, and the only automated guard on the device's MPU-protected stack reports PASS. Fix: fail explicitly when `peak == 0`.

4. ✅ **The Teensy warning ratchet passes on an empty or unparsed build log** — `tools/teensy_warnings.py:112-124` — `extract_warnings("")` yields an empty set, so a total capture break is indistinguishable from today's healthy green (the baseline file is currently empty). Fix: require evidence of real compiler invocations in the log before comparing.

5. ✅ **`parse_profile.py validate` certifies a capture containing no data** — `tools/parse_profile.py:572-638` — every check is skip-on-absent, so a log of three bare header lines prints four PASS lines and `VALID`, exit 0. The docs make this the mandatory pre-trust step for every device measurement. Fix: make the ppm-exactness check unconditional on a root counter plus wall being present.

6. ✅ **`parse_profile.py frames` mode crashes on every real capture** — `tools/parse_profile.py:664` — `frame_rows` entries became 4-tuples at line 188 when per-frame preset attribution landed; the `frames` branch still unpacks three, raising `ValueError`. Fix: `for n, wall, render, _owner in w.frame_rows:`. (Confirmed by direct inspection.)

### P1 — High: correctness defects with a real failure mode

7. ✅ **daydream: N=8 segment→band assignment diverges from the firmware, justified by a false doc claim** — `segment_layout.js:26-30,60-76` vs `hardware/pov_segment_map.h:78-86` — the JSDoc states the firmware rejects N > 4; `pov_segmented.h:102` static-asserts `N <= 8` and `test_pov_segmented.h:39-41` pins the 8-board `y_base` values. Firmware tiles north bands top-down then south bands from the pole inward; the JS tiler is purely top-down, so segments 2/3 and 6/7 are swapped at N=8. An 8-board bring-up reading the sim's per-segment overlay instruments the wrong board. `segment_crosscheck.test.js:58` asserts `segsPerArm <= 2`, which is exactly why this is invisible.

8. ✅ **`SphericalField::sample_bilinear`'s single-step x wrap is unguarded and the shipped caller can exceed its domain** — `core/math/spherical_field.h:176-180,218-222` with `core/render/filter.h:1706-1711,1797-1806` — the seam correction can lift a stored offset to ±W (which also overflows the `int16_t` store at 288×128) and `composite_previous_frame` applies a second ±W/2 correction, so `x0` can remain outside `[0,W)` after one adjustment and index past the framebuffer at `y0 == H-2`. Fix: route `x0` through `hs::fast_wrap` — byte-identical codegen plus the debug assert the codebase already uses for this contract.

9. ✅ **`classify_faces_impl` binds its output inside its own scratch scope** — `core/mesh/mesh.h:526` vs `:686` — `mesh.topology.bind(persistent, F)` runs at 686, above the `ScratchScope` opened at 526; when `persistent` aliases `scratch_a` (as at `effects/HankinSolids.h:741`) the destructor rewinds the topology block before `classify` returns. It survives only because a bump arena does not scrub. Fix: hoist the bind and zero-fill above the guard — behaviour-identical, zero cost.

10. ✅ **`classify_faces_impl` never validates the flat-face invariant its two siblings both check** — `core/mesh/mesh.h:529-638` — it is the one entry documented and tested to run on *uncompiled* `PolyMesh` input, allocates at `I = get_faces_size()`, and advances once per `(face,k)` i.e. `sum(counts)` times. A mesh with `sum(face_counts) > faces.size()` overruns three arena allocations silently. Fix: add `require_flat_face_length` after the accessor reads.

11. ✅ **`CompiledHankin::clone` produces a corner-less clone from a borrow-mode source** — `core/mesh/hankin.h:76-90` — in borrow mode `base_vertices` is empty, so `dst.corner_src` points at a null buffer while the dynamic instructions still carry in-range corner indices; the comment at line 88 asserts exactly what the code fails to do. Latent (the sole borrow-mode compile is not currently routed through `Persist<>`), but one arena compaction across a hankin leg would corrupt every star point silently. Fix: `HS_CHECK` the borrow-mode source.

12. ✅ **`SDF::Ring`'s centerline fast path systematically under-covers the inner arc edge** — `core/render/sdf.h:608-633` — the emitted span is a first-order symmetric linearization of `acos` about the centerline; `acos` curvature makes the true inner half-width strictly larger. Demonstrated numerically: at W=288, axis `ny=0.2`, radius 0.675, thickness 0.01, φ=0.319, **4.48 columns are dropped** — 45% of the inner half-stroke, with stroke alpha 0.099 at the cut, i.e. a hard 0→10% step along the ring's inner rim. Fix: delete the fast path and route every row through `emit_annular_band` (one extra `fast_acos` per row, amortized over tens of pixels). The existing conservativeness harness cannot see this because its interior predicate uses a 1-px margin.

13. ✅ **`SDF::Flower` omits the one-pixel AA pad every sibling solid applies** — `core/render/sdf.h:3881-3885` — `PlanarPolygon`, `SphericalPolygon` and `Star` all pad by `cosf(radius + pixel_width)` "so the outer AA fringe is scanned"; `Flower` uses `cosf(thickness)` unpadded, so the entire outer half of the silhouette AA band emits nothing (a pixel 0.25 px past a petal tip has alpha 0.275). Reachable from `ShapeShifter`. Fix: `cosf(thickness + 2.0f * PI_F / W)` — zero per-pixel cost.

14. ✅ **`Plot::Ring::draw` plots the seam vertex twice** — `core/render/plot.h:1899-1901` — `Ring::sample` appends an overlap-close vertex at vertex 0's position, but `draw` passes neither `close_loop` nor `omit_end`; every other overlap-close primitive pairs that vertex with one of the two. The pixel at phase 0 receives two alpha composites per frame — a fixed hot dot rotating with `phase`, in RingShower and Thrusters. Fix: `.omit_end = true`.

15. ✅ **`TransformerPool::prepare_frame`'s duck-typed hook silently opts structs out** — `core/animation/params.h:635,732,796,799` with `core/engine/transformers.h:258-264` — three params structs name the same refresh hook three ways (`prepare_thresholds()`, `sync()`, `prepare_threshold()`); the detector matches only `sync()`/`refresh_from()`, so `RippleParams` receives no per-frame refresh and live slider edits never reach a spawned ripple. A new field struct following `RippleParams`' spelling renders against spawn-time bounds with no diagnostic. Fix: rename to `sync()` and delete the alias so exactly one name exists.

16. ✅ **Per-face segue hooks are documented optional but are mandatory-together at the only consumer** — `core/animation/mesh.h:2196-2204` vs `:2411-2424,:2438-2480` — `IslamicStars.h:388-414` gates only on `face_offset` then calls `face_fade_frac` and the 3-arg `face_phase` unconditionally; `Shockwave` and `Breakdown` define neither and `Segue::Base` supplies no default, so selecting either is a hard compile error rather than the documented graceful degradation. Both policies read as ready and are effectively unshippable. Fix: add a defaulted `face_fade_frac` to `Base` and widen both policies' `face_phase` to 3 args.

17. ❌ **`noise_transform`'s channel offsets do not decorrelate at low frequency** — `core/engine/transformers.h:565-579` — the 100/200 translation is applied in *pre-frequency* space, so the effective separation is `100 × frequency`. Measured over a 200×200 sphere sweep: `corr(nx,ny)` = −0.58 (WavyTrails), −0.67 (LooseWormhole), +0.16 (Miasma), +0.17 (MeltingHi). The tangential slide has a fixed directional bias instead of being isotropic — a "combed" warp on 4 of the 12 shipped presets. **Rejected: the offset is not what drives the correlation.** The mechanism claim holds (separation is `100 x frequency`), but the flagged presets sample a window spanning only 0.16-1.0 noise units, so every channel is a near-linear gradient across the sphere and the correlation is set by the window, not the offset. Measured: scaling the offsets to 8 noise units makes it WORSE (WavyTrails -0.58 -> -0.65, LooseWormhole -0.67 -> -0.69); three fully independent seeds still leave Loose/TightWormhole at corr(ny,nz) = +1.000 and push Miasma/MeltingHi from +0.16/+0.17 to -0.29/-0.24. Sub-cell sampling is what gives those presets their smooth swirl, so no fix is available short of retuning scale/frequency, which changes their identity.

18. ❌ **`stereo_noise_warp` has the same defect in its shipped configuration** — `core/engine/transformers.h:655-659` — `Liquid2D`/`Flyby` leave `mFrequency` at FastNoiseLite's default 0.01, making `CHANNEL_OFFSET = 100` exactly 1 noise unit; measured `corr(dx,dy)` = −0.246/−0.247/−0.254 across the slider range. **Rejected for the same reason as 17.** This field is sampled over 0.02-0.12 noise units, deeper sub-cell still: measured corr(dx,dy) is -0.253 as shipped, -0.733 with the offsets scaled to noise units, and +0.243 with independent seeds -- the same magnitude, a different draw of the same lottery. No scheme decorrelates; all of them change how Liquid2D and Flyby look.

19. ✅ **`bake_blend` is the one quantizer that skips the clamp-before-cast invariant** — `core/color/composition.h:999` — `static_cast<uint16_t>(w * 65535.0f + 0.5f)` with no clamp; the `w <= 0` / `w >= 1` gates are both false for NaN, so a NaN weight (e.g. a degenerate leg span) reaches the cast as UB and writes NaN alpha into all 256 LUT entries. Fix: use `frac_to_q16(w)`.

20. ✅ **The sRGB decode consumer re-derives the generator's bucket geometry with no compile-time tie** — `core/color/srgb_decode.h:12-18,33-41` vs `scripts/generate_srgb_decode.cpp:14-17` — the shifts and hence array lengths are hand-copied; a retune regenerates shorter arrays while the consumer's counts still compute 256/480, and the static-init copy loop reads past the flash arrays **before `main()`**. Fix: `static_assert(std::size(...) == ..._N)`.

21. ✅ **Unguarded persistent-array write when a build leg lands via `arrival_mesh`** — `effects/IslamicStars.h:1176,:1449` — `carry_landing_to_seed()` and `finish_build()` check `landed_faces` against the leg's own arrival count, not the `MAX_BUILD_FACES` (1152) array capacity; every other producer has the explicit capacity check. Most roster recipes end in `Op::HANKIN`, so this path runs on every build; an over-count silently corrupts the adjacent persistent arena. Fix: add the missing `HS_CHECK`.

22. ✅ **Newborn-slot crossfade origin is resolved on the closed mesh** — `effects/HankinSolids.h:301-318` with `:823,:846` — the fallback branch nearest-matches a rosette centroid computed at angle 0, the bookend where every rosette has collapsed onto its host's rim — exactly the misfiling mode `resolve_host_faces`' own comment documents and works around. A newborn strap slot opens its crossfade from a neighbouring face's colour. Fix: call `resolve_host_faces()` first and read `host_face_palette_[f]`.

23. ✅ **Scratch-A budget assertions omit `Plot::rasterize`'s own allocation** — `effects/MobiusGrid.h:56`, `effects/HopfFibration.h:138-140`, `effects/ChaoticStrings.h:61-64` — `rasterize` binds `max(64, 2*W)` floats (2304 B at W=288) from scratch-A while the caller's Fragments buffer is live. MobiusGrid hard-codes an 8 KB split with no assertion at all and uses ~5.9 KB of it. A retune becomes a device-side runtime trap instead of the compile error the assertion messages promise.

24. ✅ **`Thrusters::on_fire_thruster` snapshots `warp_phase` after re-randomizing it** — `effects/Thrusters.h:190-198` — the comment claims the snapshot makes thrust points sit on the ring as currently displayed; that holds for `amp` but not `phase`. With the timer's 16-frame minimum against a 32-frame warp decay, the two opposed thrust rings bloom off the visible ring's crest. Fix: hoist the snapshot above the randomization.

25. ✅ **Device `Fn` returns zero on an empty call while host/WASM traps** — `core/engine/platform.h:1381-1389` — the vendored Teensy `inplace_function` overrides the throw macro to `return static_cast<R>(0)`, so invoking an unbound callable silently yields zero on the one target with no console. `tests/test_death.h:1045` asserts the trap fires but only exercises `hs::inplace_function`, giving false confidence for hardware. Fix: use `hs::inplace_function` on all three targets (or, minimally, add the ledger row and correct the comment). Took the minimal option: the divergence is now row 9 of the device/host divergence ledger (the first red row) with the comment and death-case docstring corrected; switching the device to `hs::inplace_function` was rejected as a device codegen/ITCM cost for a bug-path-only gain.

26. ✅ **`resplit_arenas()` never enforces its stated precondition** — `core/engine/memory.cpp:121-146` — the contract says both scratch arenas must be empty; the body unconditionally rebinds them, and a `ScratchScope` whose `saved_offset` is 0 cannot detect the aliasing. Every neighbouring precondition in the file traps; this one is doc-only. Fix: one `HS_CHECK` at the top.

27. ✅ **`pixel_to_vector(float,float)`'s near-integer snap can index one past the end** — `core/math/geometry.h:307-314` — the float overload delegates to the integer LUT path without range-checking either coordinate, so `x == W` reads `sin_theta[W + W/4]` on an array of that exact size; the only guard is a debug `assert` compiled out on device. Asymmetric with `y_to_phi<H>(float)`, which does check. Reachable via `Pipeline::plot(Canvas&, float, float, …)` forwarding an unwrapped screen x.

28. ✅ **Doxygen publishes nothing for the entire WASM bridge** — `Doxyfile:28` — `PREDEFINED` omits `__EMSCRIPTEN__` and all of `wasm.cpp` is inside `#ifdef __EMSCRIPTEN__`, so `HolosphereEngine`, `MeshOpsWrapper`, `PaletteOps` and every exported function publish as an empty File Reference (verified against the checked-in generated HTML). The README's "API documentation" link therefore omits the whole JS-facing surface. Fix: define the macro; two entities need `@param` first, since `WARN_AS_ERROR` is on.

29. ✅ **The WASM target compiles with no warning flags at all** — `CMakeLists.txt:79,86` — native tests use `-Wall -Wextra -Werror` and Teensy uses `-Wall -Wextra` plus a ratchet; the Emscripten branch sets only optimization flags, and `wasm.cpp` is compiled by no other target. It is the one first-party TU with zero warning coverage. Fix: add `-Wall -Wextra` (no `-Werror`) to the shared `if(EMSCRIPTEN)` block.

30. ✅ **`docs/opchain_morph_spec.md` is labelled "DESIGN, not implemented"** — `docs/opchain_morph_spec.md:3` — it landed in three merges and `effects/IslamicStars.h` cites its sections as the implemented design. A reader trusts the status line first. Fix: mark LANDED with the merge SHAs. Status line corrected to LANDED (810db7c6, 1138ff8e, 661f1a78) in the working tree; the spec itself is untracked (finding 65), so the correction is not carried by this commit.

31. ✅ **KiCad board generation is not reproducible — every run rewrites every UUID** — `hardware/phantasm/gen/kicad_common.py:16` — `uuid.uuid4()` per footprint/wire/label/junction/pin. Measured on the live working tree: 8,655 changed lines in the generated board, 1,374 of them UUIDs, making the substantive change (a corrected Teensy pad map) unreviewable and letting a regression ride in invisibly. Fix: `uuid.uuid5(FIXED_NAMESPACE, stable_key)`.

32. ✅ **Netlist connections vanish silently when no pad name matches** — `hardware/phantasm/gen/pcb.py:157-164` — `pad_net.get((ref,padname))` is a one-way lookup, and the Teensy footprint keys pads by pin *label*, so a single label typo drops that net from the board with no output — precisely the class of change currently in flight in the working tree. Fix: track consumed keys and fail on leftovers.

33. ✅ **A failed DRC run is displayed as "clean" and stays eligible to win** — `hardware/phantasm/gen/analyze_candidates.py:124,259-266,328-342` — `run_drc` returns `None` for both "tool unavailable" and "this candidate failed", so a failed candidate renders as `DRCerr 0 unconn 0 clean` and can be selected as best-by-composite. Fix: distinguish the two cases and exclude un-gated candidates.

34. ✅ **`tools/gen_gamut_lut.py` has no provenance gate** — the only generated artifact with no `--check`, no ctest and no CI job, while it hand-duplicates the OKLab matrices and gamut epsilon from `core/color/color.h`. Drift there silently invalidates the brackets the runtime bisects within. Fix: add a `--check` mode wired like `reaction-graph-provenance`.

35. ✅ **`device_lock.sh` stale-lock break is a rm-then-mkdir race** — `tools/device_lock.sh:168-170` — two peers can both judge a claim stale, and the second `rm -rf`s the first's fresh lock — two sessions flashing one board, the exact failure `mkdir` atomicity was chosen to prevent. Fix: break by atomic rename.

36. ✅ **daydream: a deferred segment-pool spawn survives page teardown** — `daydream.js:713-733,:837-857` — the toggle handler awaits `warmModules()` (three fetches) before `segments.create()`; `disposeApp()` calls `destroy()` but neither clears `active` nor bumps the epoch, so the continuation resumes after teardown and spawns N workers each loading ~1.9 MB of WASM into a discarded page. Fix: bump the epoch and clear `active` in `disposeApp()`.

37. ✅ **daydream: the worker-pool toggle and segment count are deep-linked, against the documented rule** — `daydream.js:713,:726` vs `gui.js:301-309` — `addSession`'s JSDoc names "worker-pool toggles" as exactly the class that must not auto-activate from a copied link, yet both use the deep-linked `add()` and `applyOnLoad` replays them at construction. Anyone opening a shared `?…segmented=true&…segments=8` link silently spawns 8 workers before the main engine loads. Note the inversion: the harmless `Show Boundaries` toggle is correctly session-only. (Reported independently by two audits.) Resolved by changing the documented rule: render-mode controls, the segmented pool included, are deliberately deep-linked.

38. ✅ **daydream: `destroy()` leaves `frameComposited` latched, baking black frames into a recording** — `segment_controller.js:509-544` — `tick()` early-returns before touching the flag while the pool respawns, so `captureReady()` keeps reporting true; changing the segment count mid-recording encodes 8–48 fully black frames. Fix: reset the flag in `destroy()`.

39. ✅ **daydream: the desktop effect panel can overflow the viewport with no scrollbar** — `styles/index.css:539-553` — `max-height: 100vh` with no `overflow`, and panels use `autoPlace:false` so lil-gui's own height cap never applies; `.layout-container` is `position:fixed; overflow:hidden`, so the page cannot scroll either. `DisplacementField`'s 14 params at 175–200% zoom put the bottom sliders permanently out of reach. Fix: `max-height: 100vh; overflow-y: auto` on the panel.

40. ✅ **daydream: CDN modules load with no Subresource Integrity** — `vendor-importmap.js:48-66` — versions are exact-pinned but the generated importmap emits only `imports`; a compromised jsdelivr response for `three` or `lil-gui` executes with full page privileges on the deployed Pages site. This is the one unverified link in a chain that otherwise binds the `.wasm` to a recorded sha256. Fix: emit an `integrity` map from the installed `node_modules` hashes.

41. ✅ **daydream: a bfcache freeze permanently disables page teardown in all four tools** — `tools/lissajous.html:162`, `mobius.html:735-743`, `palettes.html:1443-1453`, `solids.html:795-798,:827` — `{once:true}` combined with an early return on `e.persisted` removes the listener on the freeze it deliberately ignores, so a later real unload releases nothing (RAF loop, resize listener, THREE GPU resources, a pending timeout). Fix: drop `{once:true}` — the bodies are already idempotent.

42. ✅ **daydream: `applyOp` silently skips an unknown op, so preview and exported C++ diverge** — `tools/solids.html:1652-1663` with `:1854-1858` — a null result is skipped rather than reported, so an unbound op renders as a no-op, passes `chainIsValid`, and is still emitted by the recipe generator with a wrong `V=/F=/I=` comment. Fix: throw; both call sites already sit inside a try/catch.

### P2 — Medium: latent defects, contract gaps, and documentation that misleads

43. ✅ **README documents `MeshMorph`, which no longer exists** — `README.md:212,801,827,828` — removed in `0ce7e6ac`; the tree listing, fragment table and animation-types table all still name it, while `OpLeg` — ~900 lines of public API driving the entire Conway/opchain morph engine — appears nowhere. `:828` additionally claims HankinSolids "reuses the buffered pair" (it drives `OpLeg` directly) and omits DreamBalls' standalone `Segue::Crossfade` use.

44. ✅ **README's effect count is wrong in two places** — `README.md:217,303` — both say 26 effects (217 adds "27 headers"); `effects.h` defines 24 rows and `effects/` holds 25 headers. §9 itself documents exactly 24.

45. ✅ **README's IslamicStars parameter list does not match `register_param`** — `README.md:1592` — documents "Fade" and "Debug BB" (neither exists) and omits "Trans Speed". The prose also omits the effect's headline behaviour, the op-by-op recipe build chain.

46. ✅ **README's HankinSolids entry has no Parameters line and mis-describes the effect** — `README.md:1600-1603` — it registers "Intensity" and "Angle"; the one-liner describes neither the Conway edge-graph random walk, the on-screen morph legs, nor the interlace-angle sweep.

47. ✅ **README's SphericalHarmonics entry omits its parameter list** — `README.md:1612-1616` — the only in-scope §9 entry with no Parameters line; it registers "Amplitude" and "Debug BB".

48. ✅ **README's GS entry omits its parameters and the reaction lifecycle** — `README.md:1566-1568` — registers Feed/Kill/dA/dB/Speed, and its most visible behaviour (auto-dissolve and reseed when the field stalls) is undocumented.

49. ✅ **README names the wrong low-resolution-only effect pair** — `README.md:1546` — says "RingShower, Dynamo"; `effects.h:78-112` excludes **Thrusters**, not RingShower, and the drift `static_assert` says so literally. *Diagnosis did not hold: the sentence names which simulator preset each screenshot was captured at, not the Phantasm firmware playlist. RingShower.png and Dynamo.png are the 20-row captures, so the named pair is correct. Fixed the real defect on the line instead — the capture-resolution rationale, which claimed effects are registered at only one resolution.*

50. ✅ **README mis-describes MindSplatter** — `README.md:1742-1746` — "random-walk particle system" describes the view orientation, not the particles, which are emitted from cube vertices toward octahedral attractors and punched out by signed-axis holes.

51. ✅ **README's Feedback filter row says "Stateless — no internal frame storage"** — `README.md:529` vs `core/render/filter.h:1968-1972` — it now carries a persistent warp cache with `STORAGE_BYTES` and `init_storage(Arena&)`. An effect author trusting "stateless" gets a correct image that silently rebuilds the entire control field every frame, with no assert, log or test failure. Only `MeshFeedback.h:121` calls it.

52. ✅ **README describes a coarse warp grid the code no longer builds, and omits `DirectAntiAliasSink`** — `README.md:521-523,529` vs `core/render/filter.h:1085,1565-1578` — the field is a spherical latitude-ring control lattice, not `W/DS × H/DS`; and `DirectAntiAliasSink` is a shipped, device-used terminal sink with a non-obvious per-frame `prepare()` contract, absent from the docs.

53. ✅ **README's `meta` operator description is wrong on two counts** — `README.md:1126,1139` — `meta` is `kis ∘ dual ∘ ambo` (three steps), not "kis ∘ ambo", and being odd-length it returns in **`target`**, not `temp` as the load-bearing polarity paragraph claims. A dedicated test exists precisely to pin `meta ≠ kis(ambo)`. (Reported by three audits.)

54. ✅ **README's SDF register table says `v1` is unsigned** — `README.md:661` — for `SDF::Face`, `raw_dist` is the **signed** edge distance, and `fragment_edge_dist` depends on that sign. A shader author writing `fabsf(f.v1)` per the table loses the inside/outside distinction and shades every face as if at its edge.

55. ✅ **README's repository map omits ~16 files** — `README.md:160-282` — missing `spherical_field.h`, four `core/color/` headers, three `core/mesh/` headers, four `hardware/` headers, `targets/Profile/`, `wasm_predicates.h`, and 7 of 10 scripts; the daydream side omits six modules and `tests/`.

56. ✅ **README's parameter-registration API is under-documented** — `README.md:1527-1536` — `canvas.h` exposes five registration forms; §8 documents two, and §10.2's definition schema omits the `options` field that `param_marshal.h` actually marshals. MeshFeedback's and MindSplatter's documented params use APIs the README never introduces.

57. ✅ **daydream's JS test suite is invisible in the README** — `README.md:2120-2126` — §11 says "three layers run the same suite" and lists only C++ paths, while daydream ships 35 test files and a `js-tests.yml` workflow. A contributor editing `driver.js` has no documented way to know tests exist.

58. ✅ **README's palette roster is stale** — `README.md:1006` — says twenty-one; `palettes.h` defines 25 (`POPPED_PEACH`, `BLUE_LAGOON`, `ORANGE_CRUSH`, `PLUM_SUNRISE` undocumented). `POPPED_PEACH` is not incidental — it is a live `MeshPaletteBank` slot.

59. ✅ **README's rasterizer and transformer tables omit shipped primitives** — `README.md:742-756,908-916` — missing `Scan::DistortedRingStack` and `Scan::RingGroup`, and the entire displacement-field stack (`FieldTransformer`, `BallDropTransformer`, `NoiseProductTransformer`, plus the `spawn`/`spawn_pinned`/`init_storage`/`prepare_frame` ordering contract a new effect author will get wrong).

60. ✅ **README's §7.9 credits `Presets` with interpolation it does not implement** — `README.md:1205` — `presets.h` has no interpolation; callers drive the crossfade via `Animation::Lerp`, and the class doc says so correctly.

61. ✅ **README's coordinate section omits the `H_OFFSET` clipping** — `README.md:132` — the stated `y ∈ [0,H) → phi ∈ [0,π]` mapping holds only for the host build; on Holosphere hardware `H_OFFSET == 3` clips the image short of π, which a dedicated test pins as intended. `H_OFFSET` appears nowhere in the README.

62. ✅ **README says five geometry tools; §10.11 and the repo say four** — `README.md:1847` vs `:2027` — and only `solids.html` reuses `MeshOps`, not all of them.

63. ✅ **`docs/conway_morph_spec.md` names three symbols the tree contradicts** — `docs/conway_morph_spec.md:3-11` — documents an `Animation::ConwayMorph` (renamed `OpLeg`), says `MeshMorph` "remains for MeshFeedback" (deleted), and says `Segue::Dissolve` "is deleted" (it exists). Fix: mark SUPERSEDED rather than editing the body.

64. ✅ **Two specs still document the far-star guard at its old value** — `docs/conway_morph_spec.md:548`, `docs/opchain_morph_spec.md:1205` vs `core/mesh/hankin.h:294` — both state `STAR_FAR_RATIO_SQ = 36`; the constant is 4.0 with the blend starting at 2.25, and the calibration reasoning in `core/animation/mesh.h:442-448` no longer follows.

65. ✅ **Eight engineering docs are untracked** — including the 71 KB `opchain_morph_spec.md` that shipped code cites. They are unbacked, unshared, and invisible to `tools/docs_check.py`, which walks tracked entries only — so their fences and links are never validated.

66. ✅ **The flush-before-plot contract for a frame-replacing terminal filter is undocumented** — `core/render/filter.h:34-36,1478,1853` — a `terminal_replaces` stage's `flush()` overwrites the destination at `alpha >= 1`, so it must run *before* the frame's `plot()` calls. `MeshFeedback` does this by convention only; an author mirroring `Dynamo`'s correct flush-last pattern would silently erase the entire frame.

67. ✅ **Calling the wrong-domain `flush()` overload is a silent no-op** — `core/render/filter.h:443-468` — a pipeline whose only history stage is `World::Trails` compiles and does nothing; because aging lives inside `flush`, the ring buffer then fills to capacity and never decays, with no diagnostic. Fix: add `any_2d_history`/`any_3d_history` folds and `static_assert` in each overload.

68. ✅ **`Pipeline` publicly inherits `Head`, leaking the head stage's traits onto the pipeline type** — `core/render/filter.h:239` — `Pipeline<…>::crosses_segments` reports the *first* stage's value, not the fold. Every effect correctly uses `any_crosses_segments`, but the wrong spelling compiles and silently yields `full_frame = false` → cross-segment corruption on the segmented driver. Fix: shadow the folded value in the recursive case.

69. ❌ **`Feedback::flush`'s all-black early return assumes the destination is already black** — `core/render/filter.h:1484-1488` — true only because the `Canvas` ctor zeroes it; paired with a `persist = true` effect the band holds the previous frame and feedback silently freezes the image instead of fading it. REJECTED (non-issue): the `persist = true` path does not hold a *different* frame — `advance_buffer()` memcpys `bufs_[next_]`, and `HS_CHECK(next_ == prev_)` makes that byte-identical to what `any_pixel_lit` just scanned, so the destination is black exactly where the scan found `prev` black; `persist = false` gets the same guarantee from `clear_buffer()`. The scanned region is also exactly the written region (`make_render_band` reads the same `render_y_start`/`render_y_end`/`x_clip`), and both `color_fn`s (`plain_fade`, `hue_fade`) are linear, mapping black to black. Confirmed empirically: the full 56-test suite passes unchanged with the early return deleted.

70. ❌ **The scan seam-split/coalesce walk is duplicated verbatim three times** — `core/render/scan.h:216-265,918-954,1192-1226` — plus a fourth triplication of the clip-arc splitter. Two copies explicitly claim to "mirror scan_region so walked pixels are identical", enforced only by comment. Fix: extract sink-parameterized header templates (inlining preserves the O3-region property). REJECTED (risk out of proportion to the gain): the premise that "inlining preserves the O3-region property" does not hold here. Two of the three copies *are* the pixel-loop bodies of `HS_O3` regions (`RingGroup`, `rasterize_face`), and a template takes its options at its **definition** site (docs/selective_o3_spec.md §4.2 #4), so a helper defined outside the regions is an `-Os` callee holding the hot sink; §4.2 #2 only claims plain-inline callees are *generally* inlinable into a region, and §7 records GCC 11 refusing exactly that for several of them. A refusal here compiles the DisplacementField/RingSpin and mesh-face per-pixel loops at `-Os` with no host-visible symptom — the native suite sees `HS_O3` as a no-op, so only a device profile could catch it. The alternative, wrapping the helper in its own region, option-carries it into every `scan_region` instantiation, whose ITCM cost the §3 budget cannot absorb (the shared R2 kernel measured ~28 KB). The refactor is source-only — each of the three sink/buffer type pairs instantiates its own copy either way — so it buys no bytes and no behavior, at the cost of a device-only regression class. Revisit only with a device profile in the loop.

71. ✅ **`Scan::DistortedRingStack` omits its `ScopedRenderTimer` and its documented `n_slots` check** — `core/render/scan.h:550-645` — every other draw entry point opens a timer, so DisplacementField's normal render contributes 0 µs to the profiler while its debug fallback reports the real cost — the same frame timed differently by a debug toggle.

72. ✅ **`Scan::DistortedRingStack` has no test** — its doc makes the strongest claim in the file ("bit-identical to rasterizing the rings one by one") and nothing tests it, while its sibling `RingGroup` has exactly that test.

73. ✅ **`scan_region` clears its interval buffer only on the handled branch** — `core/render/scan.h:266-269` — a producer that emits and then returns false would leak spans into the next row and eventually trip an overflow trap far from the cause. No current producer does this, and the contract forbids it, but it is unasserted.

74. ✅ **`rasterize_face`'s `MinimalFragment` path never resets `frag.color`** — `core/render/scan.h:1326-1328` — the fragment is per-face and only `v1` is written per pixel, so a shader that early-returns on a cull condition smears the previous pixel's colour across the run. Safe today only because the sole shader writes unconditionally. Fix: state the contract on the setup parameter.

75. ✅ **`Scan::Shader` overwrites the canvas and bypasses the filter pipeline, undocumented** — `core/render/scan.h:1498,1607,1683,1765` — all three overloads assign rather than plot, so an effect returning `alpha < 1` gets a darkened opaque pixel and no `Filter::World::*` applies.

76. ✅ **A caller-supplied `fragment_interpolator` is silently replaced by `Fragment::lerp` on the type-erased path** — `core/render/plot.h:1116-1127,3593-3599` — the interpolation policy survives only for pipelines declaring `direct_raster_path`. Today's only caller supplies a strict subset of `lerp`, so the drop is a pure perf regression — but any nonlinear interpolator would render differently depending on pipeline type, with no diagnostic.

77. ✅ **`ParticleSystem::draw_impl` re-implements `gate_trail_edges`, and the two have already drifted** — `core/render/plot.h:3397-3463` vs `:2870-2947` — a near-verbatim copy missing the near-horizontal-arc-pole rejection whose comment explains exactly why it is needed. No wrong cull could be constructed, but two bodies documented as producing identical verdicts differ by one guard in the place a divergence is hardest to notice.

78. ✅ **`SDF::Intersection` can emit spans and still return false** — `core/render/sdf.h:1720-1724` — combined with finding 73 this would double-shade the next row. Unreachable today only because every leaf returns false before emitting — an invariant nothing states or tests.

79. ✅ **Antipodal `SDF::Line` endpoints: rendered geometry and vertical bounds disagree** — `core/render/sdf.h:3962-4004` — the degenerate-normal fallback renders the entire great circle `x = 0` while the phi bounds stay pinned to a band around the endpoints' latitude, clipping the circle to a sliver.

80. ✅ **`SDF::Flower` is the only solid shape missing from the cull-conservativeness grid** — `tests/test_sdf.h:1526-1552` — which is why finding 13 survived.

81. ❌ **The CSG layer has no production caller** — `core/render/sdf.h:1162-1931` — ~770 lines instantiated only by tests, yet it sets the `INTERVAL_SPAN_CAP = 32` sizing that costs `scan_region` ~1.6 KB of scratch per call (leaves emit ≤ 2 spans each). Worth an explicit keep/drop decision. REJECTED (needs user sign-off): test-only status confirmed (no Union/SmoothUnion/Subtract/Intersection instantiated outside tests/), but deleting a tested ~770-line subsystem is a product decision, and the INTERVAL_SPAN_CAP retune it would enable is separable — scan_region cannot size its scratch per-shape without threading a span-bound template parameter through every call site.

82. ✅ **`build_mesh_class_bake` guards `topology` but not `face_offsets`** — `core/mesh/mesh_classes.h:174-177,208` — dereferences `fo[f]` with no check that offsets are bound or that the span fits, while `scan.h:1415` guards exactly this per face.

83. ❌ **`MeshState`'s owned/borrowed mode is discriminated per array, not per mesh** — `core/mesh/spatial.h:362-405,444-467` — a mesh with only some arrays bound would silently serve owned counts against borrowed indices with no trap. The invariant holds today by convention alone. REJECTED (design-pinned): face_offsets is deliberately optional — DreamBalls binds an edge-only MeshState without it and feeds it through MeshOps::transform, so a per-mesh discriminator or an all-arrays bind-time trap would fire on legitimate meshes; the guard belongs at the consumers that dereference it (finding 82).

84. ✅ **`Recipe::seed` indexes `simple_registry` unchecked at three sites** — `core/mesh/solids.h:1138`, `core/mesh/recipe.h:120,186`, `targets/wasm/wasm.cpp:1153` — while the parallel `get_entry` path `HS_CHECK`s. Safe only because every current recipe uses a static-asserted constant. Fix: a constexpr fold over the registry, or one `HS_CHECK`.

85. ✅ **`expand`/`chamfer`/`snub` document a `[0,1]` parameter and never enforce it, while `truncate` traps** — `core/mesh/conway.h:842,913,1170` vs `:742` — these are the operators whose parameter is *animated* at runtime; a clamp regression feeding `t < 0` inverts the mesh with no trap, on hardware with no console.

86. ✅ **The load-bearing scratch-arena contract block omits `medial` and misstates `relax`** — `core/mesh/conway.h:381-398` — `medial` allocates in both arenas and requires a closed manifold, and is absent from both lists; `relax` allocates `orbit_start` and `movements` in `temp` despite the block saying it needs no extra buffers.

87. ✅ **`relax_baked`'s doc promises a check on source "raw data" that it does not perform** — `core/mesh/conway.h:1107-1116` — only dimensions and a connectivity hash are checked; source vertex positions are deliberately unchecked, which is what lets two direct call sites replay a swept recipe onto the shipped converged mesh.

88. ✅ **`relax` drops roughly half of a boundary mesh's edges from its target-length average** — `core/mesh/conway.h:1043-1052` — the `u < v` half-edge filter is correct for interior edges but arbitrary for unpaired boundary ones, biasing the spring target for the whole mesh. Latent (production relaxes only closed manifolds) and the open-mesh test only asserts finiteness.

89. ✅ **`update_hankin`'s Doxygen block is orphaned onto a constant** — `core/mesh/hankin.h:277-311` — six constant declarations intervene between the block and the function, so the one API here with a non-obvious contract ships undocumented.

90. ✅ **`resolve_host_faces`' doc states the opposite of its implementation** — `effects/HankinSolids.h:552-558` — claims the host is found by shared vertices "and is therefore angle independent"; the body rebuilds at π/2 and matches by centroid, with an inline comment explaining why shared vertices *don't* discriminate. A reader trusting the doxygen will "fix" the code back to the broken approach.

91. ✅ **`MeshFeedback` syncs noise before the preset switch** — `effects/MeshFeedback.h:141-153` — `apply_params()` runs before `advance_transition()`, so on each switch frame the warp reads the previous preset's noise scalars while fade/hue are already the new preset's. Fix: swap the two calls.

92. ✅ **A bare `assert` sits in the DisplacementField render path** — `effects/DisplacementField.h:213` — the only bare assert in `effects/`; it compiles out on device, so the guard exists only in the host build.

93. ✅ **`PetalFlow::MAX_RINGS` equals its worst case with zero margin** — `effects/PetalFlow.h:90-91` — 64 required, 64 provided, by coincidence rather than construction; widening the Density slider or nudging `SPACING` starts dropping rings every frame. Fix: a `static_assert` deriving the bound from the four constants and the slider max.

94. ✅ **Three effects lack the persistent-footprint `static_assert` their siblings carry** — `effects/RingShower.h:43-51` (~48 KB), `effects/DisplacementField.h` (~208 KB at W=288), `effects/RingSpin.h:67,87` — runtime coverage exists via `unit_arena_budget`, but a capacity retune has no compile-time guard.

95. ✅ **`DreamBalls` silently swallows a full warp pool** — `effects/DreamBalls.h:317-319` — a failed spawn leaves the "Warp" slider inert for a full 320-frame sprite with no signal; correct today only via an undocumented event-ordering chain.

96. ✅ **`ShapeShifter`'s constructor is the only non-cold effect constructor in the repo** — `effects/ShapeShifter.h:38` — every other effect ctor is `HS_COLD_MEMBER` so one-shot setup lands in flash rather than ITCM, against the documented granule cliff.

97. ✅ **Implementation methods are public in two effects** — `effects/ShapeShifter.h:100,133,166` and `effects/Dynamo.h:107-177` — the other five expose only ctor + `init()` + `draw_frame()`; Dynamo already has a white-box friend, so nothing needs the wider surface. `color_wipe()`'s palette/boundary invariant is delicate enough that `draw_frame` guards it with `HS_CHECK`.

98. ✅ **`Voronoi::CellId::has_second` is write-only and costs ~5.4 KB of scratch** — `effects/Voronoi.h:285` — never read (the shading path derives its own), and it pushes `sizeof(CellId)` from 4 to 6 across the whole corner grid.

99. ✅ **`Voronoi` exposes internals publicly instead of using the house white-box seam** — `effects/Voronoi.h:22` — GS, BZ, Liquid2D and Flyby all solve the same test-access problem with a `friend struct …WhiteBox`.

100. ✅ **GS never uses the base's substep driver that the base documents as shared** — `effects/ReactionDiffusionBase.h:216-241` vs `effects/GSReactionDiffusion.h:455-459` — GS hand-rolls the identical ping-pong; with an even step count the shared path is bit-identical.

101. ✅ **GS's `draw_frame` override duplicates the base body with a stale stated reason** — `effects/GSReactionDiffusion.h:175-186` — the comment claims the base doesn't bracket the phases (it does); the real reason is the `grd_*` scope names the shipping profile report is written against. As written, a step added to the base silently skips GS.

102. ✅ **The RD device-scratch `static_assert`s cannot fire in the native test build** — `effects/GSReactionDiffusion.h:136`, `effects/BZReactionDiffusion.h:100` — they score against the 8 MiB test arena; GS is the binding tenant with 4,096 B of real headroom. Fix: assert against `DEVICE_GLOBAL_ARENA_SIZE`.

103. ✅ **`TickActions::zero_crossing` is sticky across a multi-boundary tick** — `hardware/pov_sync.h:1404-1406` — set on a ZERO flip and never cleared on a HALF flip, so a folded tick publishes the wrong half-window and clips one frame to the wrong hemisphere. Requires ≥62.5 ms of ISR coast, so latent, but wrong by construction.

104. ✅ **ACQUIRE hard-snaps on a beacon frame's first digit; the README overstates the guarantee** — `hardware/pov_sync.h:1461-1470`, `README.md:1449-1452` — the quiet-before guard catches digits 2–5 but not digit 0, which is preceded by silence. A board falling back to ACQUIRE mid-show can lock with a 72-column phase error and render half-rotated for ~250 ms. README corrected; the code behaviour is the spec's explicitly accepted first-digit case (`docs/phantasm_frame_sync_spec.md` §5.3, bounded by the §9.1 mis-snap row), so it is unchanged.

105. ❌ **The beacon stale-frame timeout lacks the signed-wrap re-check its siblings have** — `hardware/pov_sync.h:828-831` — after >7.16 s of silence holding a partial frame, the modular difference can read small and a stale digit survives. Bounded (one corrupted frame, rejected by checksum), but asymmetric with the two analogous comparisons in `tick()`. The re-check would be a regression, not a fix: `feed()` is event-driven (its siblings are polled every wake, so their difference never reaches 2^31) and its two burst timestamps are always ordered forward, so adding `(int32_t)diff > 0` would *keep* a stale partial frame across a 3.58–7.16 s silence that today is correctly dropped, while the ≥2^32 case cited stays unreachable to any 32-bit comparison.

106. ✅ **`HS_CHECK` sits on the per-wake hot path** — `hardware/pov_sync.h:637-640` — `Flywheel::position()` runs at ~18.4 kHz; project policy is `HS_CHECK` on cold paths only, and the invariant is already trapped in the constructor.

107. ✅ **The DMA stale-transfer watchdog is sized and documented two orders of magnitude off** — `hardware/dma_led.h:151-154` — the comment says "single-digit ms"; the real transfer is ~99 µs. At 100 ms a wedged channel freezes the strip for ~230 consecutive column drops before trapping.

108. ✅ **The LPSPI RX path is neither masked nor drained, and `NOSTALL` is undocumented** — `hardware/dma_led.h:85-88` — every transmitted byte lands in an undrained RX FIFO; for a TX-only stream `RXMSK = 1` with `NOSTALL` cleared is the safer configuration (a stall is harmless for a clocked APA102 protocol; an underrun corrupts the frame).

109. ✅ **`HD107SFrame::load()`'s tail-blanking loop is never asserted** — `hardware/hd107s_frame.h:170-173` vs `tests/test_hd107s_frame.h:222-245` — deleting the blanking (leaving a stale tail on a partial load) passes the whole suite.

110. ✅ **`stackup.py` is not idempotent and fails open in three places** — `hardware/phantasm/gen/stackup.py:39,43-53,78-80` — unconditional zone appends, a netcode lookup that returns 0 for a missing net, and a literal whitespace-sensitive `str.replace` that silently writes no stackup while still printing success.

111. ✅ **`fab.py`'s parsed-DRC diagnostic is unreachable** — `hardware/phantasm/gen/fab.py:105-140` — `run()` uses `check=True` while DRC passes `--exit-code-violations`, so the carefully worded violation message can only fire on a board that is already clean.

112. ✅ **Locked Quilter placements are never validated against the generated outline** — `hardware/phantasm/gen/pcb.py:228-247` — absolute coordinates against a computed board length; a footprint swap that shortens the board puts the locked connectors outside it with no diagnostic.

113. ✅ **`fab.parse_components` is a third independent hand-rolled S-expression reader** — `hardware/phantasm/gen/fab.py:143-174` — a depth counter over raw text plus regex, while `sexp.parse` and `export_netlist` already exist and `check.py` uses them on the same file. A `)` inside a quoted property mis-terminates the block.

114. ✅ **`shorts.py`'s colinearity tolerance is not normalized by segment length** — `hardware/phantasm/gen/shorts.py:62-67` — the effective perpendicular tolerance is `0.02/len`, so a 100 mm wire admits 0.2 µm. Orthogonal wires are exact, so any diagonal wire silently loses mid-span detection.

115. ✅ **A missing footprint yields a *perfect* ergonomics score** — `hardware/phantasm/gen/analyze_candidates.py:184-192` — `max(0, nan-5)` evaluates to 0 in Python, so every missing part removes its own penalty and a candidate that dropped decoupling caps outranks one that merely placed them badly.

116. ✅ **The `profile` and `profile_o3` environments are compiled by no gate** — `.github/workflows/ci.yml:590,632`, `justfile:99` — `targets/Profile/Profile.ino` (404 lines) consumes a wide, fragile instrumentation API and is compiled only when a developer runs `just profile`, despite being the harness the entire on-device profiling workflow depends on.

117. ✅ **`HS_WASM_DEV_BINDINGS` is compiled by nothing** — `CMakeLists.txt:95`, `targets/wasm/wasm.cpp:1167-1227` — 60 lines touching `Solids`/`PolyMesh` APIs that can break silently until someone flips the option. Fix: add the flag to the existing WASM debug compile-check step.

118. ✅ **`platformio.ini`'s toolchain rationale contradicts its own pin** — `platformio.ini:15-20,37-45,113-119` — the file asserts three times that the pin is `teensy@5.0.0` / TD 1.59 / gcc 11.3.1 and justifies the FastLED pin from that; the pin is now `teensy@5.2.0` (TD 1.62 / gcc 15.2.1). These comments are the load-bearing justification for both the FastLED pin and the bench-parity claim.

119. ✅ **Teensy targets construct the effect before `configure_arenas_default()`; WASM does the reverse** — `targets/Phantasm/Phantasm.ino:73-77`, `targets/Profile/Profile.ino:347-350` vs `targets/wasm/wasm.cpp:373-381` — both `.ino` files document the hazard ("an arena allocation in a ctor would corrupt on Teensy yet pass on WASM") — a defect class neither the native suite nor the WASM smoke can catch, because both run the safe ordering.

120. ✅ **`CANVAS_W`/`MAX_W` are two unlinked sources of truth** — `core/engine/platform.h:7-14` vs `core/engine/constants.h:15-20` — both are 288/144 with no assertion tying them, and the header openly invites overriding one; an override above `MAX_W` fails at runtime rather than at compile time.

121. ✅ **`beat16`'s host mock is not FastLED-faithful for `bpm >= 256`, and its `@pre` claims the opposite** — `core/engine/platform.h:786-792` — FastLED treats `>= 256` as already-Q8.8; the mock shifts unconditionally. These mocks exist precisely to be bit-faithful.

122. ✅ **`effects.h` understates the native anti-drift guard** — `core/engine/effects.h:44-49` — claims a forgotten roster row "silently drops that effect from native coverage"; `tests/test_effects.h:3811` runs the same registry-count oracle unconditionally, above the FULL-tier gate.

123. ✅ **`resplit_arenas()` and `Arena::set_capacity()` have zero test coverage** — the hardest arena paths and the only repartition an effect uses mid-run; a regression resetting the persistent offset would be caught only by an on-device visual. Six always-on traps in `memory.h`/`generators.h` are likewise unexercised by the death harness.

124. ✅ **`StaticCircularBuffer::for_each()` and `is_linear()` are untested** — `core/engine/static_circular_buffer.h:283,326-333` — `for_each` is the only traversal that walks raw slots rather than through `operator[]`'s fold, and it is the hot traversal in `core/animation/trails.h:120`; `is_linear()` backs two always-on traps in `sdf.h`.

125. ✅ **`Gradient` silently yields an all-black palette for an empty stop list** — `core/color/color.h:1280-1281` — the constructor traps out-of-range and unsorted stops at the same cold seam but returns early on zero stops, which is precisely the silent failure the fail-fast philosophy exists to prevent.

126. ✅ **`Color4`'s arithmetic operators document a renderer path that no longer exists** — `core/color/color.h:224-231,289-336` — both SSAA paths now accumulate premultiplied into a `Pixel`; the four operators have zero engine call sites. An author following the doc gets the classic double-darkened AA fringe.

127. ✅ **`composition.h` is not self-contained yet is directly included by two files, with no guard** — `core/color/composition.h:5-9` — it declares itself non-standalone and uses eight types it never declares; it compiles only because both includers happen to pull `color.h` first. `animation/mesh.h` already demonstrates the `#ifndef …_INTERNAL / #error` idiom.

128. ✅ **`Style` has two sync entry points with split ownership** — `core/engine/styles.h:152,170` — `sync_hue()` is called by the filter every flush; `sync_noise()` is the effect's job and nothing verifies it ran, so a missed call silently leaves the previous preset's amplitude driving the warp.

129. ✅ **`Style::downsample`'s documented scratch requirement understates the real one ~2×** — `core/engine/styles.h:106-112` — the formula omits the `WarpControl` array (~8.1 KB at 288×144, DS=4), so an effect sizing its scratch arena from this comment traps mid-flush.

130. ✅ **`generate_reaction_graph.py` duplicates three runtime constants with no cross-check** — `scripts/generate_reaction_graph.py:41-45` vs `core/engine/reaction_graph.h:15-16,54` — the header claims both are guarded; raising `RD_N` in the header alone leaves the CI provenance diff green and produces a zero-padded table.

131. ✅ **`fast_sincosf_0_pi` states a domain it does not enforce** — `core/math/3dmath.h:1397-1409` — beyond π the sine branch silently returns the wrong sign, while the sibling `fast_expf` asserts its domain.

132. ✅ **`project_div` returns the antipodal pole for a tiny-numerator quotient** — `core/math/3dmath.h:836-847` — the 0/0 test is an absolute epsilon applied after the singular branch is taken, so a legitimate point at infinity maps to the south pole. Fix: make the test exact, as the comment already claims.

133. ✅ **`stereo()`'s doc contradicts its code and its own test** — `core/math/3dmath.h:891` — the pole cap deliberately preserves the (x,z) azimuth; only the exact pole falls back to the real axis.

134. ✅ **`rotation_substeps` casts an unbounded `ceil` to `int`** — `core/animation/motion.h:127-129` — `Rotation`'s ctor checks only finiteness, so a large sweep is UB rather than the graceful degradation `Orientation::upsample` already provides.

135. ❌ **`max_life == 1` spawns particles that are removed before they ever render** — `core/animation/sprites.h:263-265,385-390` — `spawn()` reports success and `active()` never reflects it, so a "Life" slider bottoming out at 1 emits nothing at all. Note `tests/test_animation.h:988` design-pins 1 as valid, so classify before landing. — rejected: the visibility floor comes from the trail-segment render gate (`core/render/plot.h:3342` skips any particle with `trail_len < 2`), not from the life decrement, so `max_life` 1 and 2 are both invisible whichever way the decrement is ordered — the real floor is 3, and reordering it would not make 1 render. `tests/test_animation.h:988` design-pins 1 as a valid `init()` input (the `[1, 65535]` bound is a `uint16_t` range guard, not a visibility promise) and `:1008` pins the N-1 trail arithmetic, so both candidate fixes break a design-pin. No shipped effect exposes a Life parameter; MindSplatter hardcodes 160.

136. ✅ **`deep_tween` is documented as accepting a bare `Orientation`, which cannot compile** — `README.md:857`, `core/engine/concepts.h:394-396` — the `Tweenable` concept admits it, so the failure surfaces as a deep template error rather than a clean rejection; `trails.h:282-286` states the restriction correctly.

137. ✅ **README says `ParticleSystem` uses `VectorTrail`; it uses `QuantizedVectorTrail`** — `README.md:785,820,848,872` — misstating the memory and precision contract (6 B/sample snorm16, not 12 B exact).

138. ✅ **`fill`/`grade` are silently dropped on the per-face segue draw path** — `core/animation/mesh.h:2186-2209` — the protocol presents them as one composable set with the per-face hooks; in practice the two sets are mutually exclusive, so a future policy combining `face_offset` with `fill` loses its coverage mask with no signal.

139. ❌ **The immutable birth-cohort subsystem is production-dead** — `core/animation/mesh.h:148-159,1636-1863,1954-1992` — ~270 lines of orbit-signature and union-find code gated on a flag no effect sets, reachable only from tests; a test comment additionally misdescribes which model the shipped shape uses. — rejected: production-dead confirmed (no effect sets `immutable`), but removal is out of scope for a review fix — `docs/islamicstars_palette_crossfade_plan.md` deliberately retained the engine-side machinery, and `immutable = true` is load-bearing in four general-purpose OpLeg smoke tests (medial, gated kis/dual, reconcile) whose palette expectations would need redesigning, not deleting. The misdescribing comments are corrected.

140. ✅ **The CONWAY_SWEEP settle path's arrival frame is not bitwise the classified endpoint** — `core/animation/mesh.h:860-865` — the three sibling paths shortcut at `k >= 1` explicitly so the closing bookend is bitwise; the settle blend does not, violating the file's own stated invariant on the one leg kind that settles.

141. ✅ **`MeshCarousel::compact_drop_all` is untested** — `core/animation/mesh.h:2691-2705` — the riskiest of the three compaction primitives (it evacuates nothing) with five live call sites, while both siblings have tests.

142. ✅ **`worst_a` is computed and printed but never asserted** — `tests/test_opchain_arena_survey.h:197,210` — the scratch-B budget is gated and scratch-A is not, while the shipping chains do get a per-chain gate elsewhere — a clearly unintentional asymmetry.

143. ✅ **`aa_audit_main.cpp` and `aa_search_main.cpp` are in no build target** — 338 lines referenced by no CMake file, calling APIs that have since moved. The project made the opposite call for their siblings, keeping `perf_bench` in the build explicitly as a bit-rot guard.

144. ✅ **631 lines of seam gates run under a borrowed module label** — `tests/test_conway_continuity.h:59` — `test_partition_seam.h` is invoked from inside another module's runner, so it has no CTest, `--check-modules` cannot see it, and its failures are misattributed.

145. ✅ **`capture_screenshots.mjs` writes wrong-effect PNGs when resolution probing degrades** — `scripts/capture_screenshots.mjs:190` — a short-circuit skips the effect-identity confirmation when the resolution list is empty, so the app's fallback effect is saved under another effect's filename — the outcome the code's own comment calls "worse than leaving the stale one". The nonzero exit fires only after the files are on disk.

146. ✅ **A truncated relax-bake dump is silently dropped** — `tools/relax_bakes.py:34-70` — a harness crash mid-block discards the partial bake and `emit` reports success with a short asset.

147. ✅ **Two generators emit through text-mode stdout** — `scripts/generate_luts.py:175`, `scripts/generate_reaction_graph.py:106` — CRLF on Windows, LF in CI, with no `* text=auto` in `.gitattributes`; both provenance gates compare numeric tokens only, so a whole-file line-ending flip rides green.

148. ✅ **A misplaced `unittest.main()` hides 4 of 5 test classes** — `tools/profile_tests/test_parse_profile.py:74-75` — direct invocation runs 8 tests and prints OK; discovery runs ~30.

149. ✅ **`docs_check.py` passes when it finds zero Markdown files** — `tools/docs_check.py:294` — the docs workflow step goes green having validated nothing.

150. ✅ **`profile_one.sh` hardcodes one machine's checkout path** — `tools/profile_one.sh:52` — the intent (always build the main tree, never the invoking worktree) is sound; the encoding is not.

151. ✅ **daydream: the controller silently drops unknown worker messages** — `segment_controller.js:295-353` — asymmetric with the worker, which throws on protocol drift and has a test pinning that. A renamed message leaves `ready` unset until a 20 s watchdog reports a misleading diagnosis.

152. ✅ **daydream: `count` is host-mutable ahead of the pool it describes** — `segment_controller.js:750`, `daydream.js:727` — during the module-warm await the live pool is still the old size while `count` is new; going 8→4 blits only arm A, so the right half of the sphere goes black for the length of the fetch.

153. ✅ **daydream: `MediaRecorder` errors are never surfaced** — `recorder.js:209-218` — no `onerror`, so an encoder fault silently finalizes a truncated file, the button still reads "■ Stop", and the next click starts a second session instead of stopping.

154. ❌ **daydream: picture-in-picture renders a duplicate of the main view, not the documented back view** — `driver.js:610-612` vs `README.md:1931` — both position and quaternion are copied every frame, and the backface cull keys off the main camera too. The documented rear view never existed; the full 41,472-instance mesh is drawn a second time per frame for it. The standard headless screenshot harness cannot see this (the PiP is suppressed under `navigator.webdriver`). **Rejected:** the picture-in-picture is intended to duplicate the main image; the README's rear-view claim was the error and has been corrected.

155. ✅ **daydream: the AppState subscriber is never unsubscribed** — `daydream.js:512,837-857` — `disposeApp()`'s JSDoc claims it releases the listeners this module owns; any post-dispose `set` re-enters `applyEffect()` against a disposed renderer.

156. ✅ **daydream: resolution presets are hard-coded while the engine publishes `getSupportedResolutions()`** — `daydream.js:74-77` vs `targets/wasm/wasm.cpp:1366` — the exported binding has zero references; an engine-side change leaves the dropdown offering a preset that cannot load.

157. ✅ **daydream: sort controls are unlabelled and expose no sort state** — `sidebar.js:210-230` — the direction glyph is folded into the accessible name and there is no `aria-pressed`, so a screen-reader user cannot determine the list's ordering.

158. ✅ **daydream: the tool-facing WASM surface has no contract pin** — `tests/engine_contract_wasm.test.js:22-31` — `HolosphereEngine`'s methods are pinned but `MeshOps` and `PaletteOps` — the two classes the tools actually run on — are not, and `KNOWN_OPS` is a hand-maintained mirror of the C++ op list.

159. ✅ **daydream: `NAMED_PROCEDURAL_PALETTES` mirrors 21 of the engine's 25 palettes** — `tools/palette_math.js:109-131` — already drifted; two of the four missing entries carry negative coefficients the gallery has never had to display.

160. ✅ **daydream: the generative-palette harmony and profile constants are an unpinned C++ mirror** — `tools/palette_math.js:226-265,295-323` — all values currently match `color.h` exactly, but the test only asserts they are integers in `[0,255]`, so `85 → 58` would pass.

161. ✅ **daydream: `mobius.html` has no export path and uses a different projection pole than the engine** — `tools/mobius.html:359-366` vs `core/math/3dmath.h:893-905` — the tool projects from `1 - p.z` with `(x,y)`; the engine uses `1 - v.y` with `(x,z)`. It designs exactly `MobiusParams` and is the only tool with no copy button, so transcription is manual *and* frame-rotated.

162. ✅ **daydream's LICENSE lacks the effects carve-out** — `daydream/LICENSE:1-3` — PolyForm-Noncommercial only, but the repo ships a `.wasm` linking `effects/` code that the README and the effect headers declare all-rights-reserved.

163. ❌ **`HS_O3` wraps ~1800 lines of cold one-shot generator code** — `core/mesh/solids.h:29,1800` — every generator is `FLASHMEM`/`HS_COLD` off the render path, yet the region applies O3 + fast-math to the whole namespace. Per the selective-O3 discipline each region should be justified by measurement. **Measure before changing** — fast-math inlining across the boundary could shift generator float results (the relax-bake verify ctest is the oracle). **Rejected: the region already carries its measurement.** `f1debacb` landed the four mesh-construction regions (`mesh.h`, `conway.h`, `hankin.h`, `solids.h`) with the numbers and the tradeoff stated: FLASH code 284,668 -> 303,828 B (+19,160) and RAM1 code 191,368 -> 191,048 B (-320) — the `HS_COLD` bodies are exactly why the spend is flash and not ITCM, and the region *freed* 320 B of the scarce resource. "Cold" is a placement decision, not "never latency-relevant": the generators run at effect init / shape spawn (`DreamBalls::setup_presets`, `Raymarch`'s vertex directions), and the stated intent was to stop shape construction running at -Os while the rasterizer consuming it was already -O3. The proposed oracle also does not apply — `HS_O3_BEGIN` is a no-op outside the device `-Os` GCC build (`core/engine/platform.h:1079-1092`), so the host relax-bake verify ctest cannot observe any change to this region, and the only build whose float results could shift is the one with no cheap oracle. Removing it is a measured cost against an unmeasured benefit in a resource that is not scarce. The residual is documentation: `docs/selective_o3_spec.md:227-229` still says regions contain no `HS_COLD` bodies, a rule `f1debacb` deliberately superseded and did not amend.

### P3 — Low: dead code, nits, and cosmetic drift

164. ✅ **Three dead type aliases in the live alias block** — `core/engine/concepts.h:217-219` — `TransformFn`, `SpaceTransformRef` (identical types) and `ColorTransformRef` have zero references repo-wide.

165. ✅ **`hs::clamp(float,…)` is `constexpr` on Cortex-M7/WASM but not on x86** — `core/engine/platform.h:1436,1460` — one interface, two contracts.

166. ✅ **`engine.h` says `effects.h` pulls in 27 effect headers; it pulls in 24** — `core/engine/engine.h:10`.

167. ✅ **Dead includes in `util.h`** — `core/engine/util.h:17-18` — `<memory>` and `<utility>` are unused.

168. ✅ **`pop_front_internal()`'s doc names callers that do not exist** — `core/engine/static_circular_buffer.h:364-371` — `pop()` has no such function and `clear()` never calls it; the sibling `pop_back_internal` doc is correct.

169. ✅ **`plain_fade` and `identity_warp` ship in no roster** — `core/engine/styles.h:60,77` — seven presets appear only in tests, yet both transforms carry live fast-path support in `filter.h`. Note the `static_assert` anchor at `styles.h:362` pins one of them and must be repointed first.

170. ✅ **`SphericalFieldLayout`'s accessors recompute from scratch** — `core/math/spherical_field.h:57-81,263-269` — `ring()` is O(k²) and `ring_index_at_or_before` calls `ring_count()` in its loop condition, both on per-frame paths. Fixing is perf-positive with no semantic change (all values are `static_assert`-pinned).

171. ✅ **`spherical_field.h` uses `assert()` without including `<cassert>`** — `core/math/spherical_field.h:292`.

172. ✅ **Dead public types `Dot`/`Dots`/`Points`** — `core/math/geometry.h:46-75` — zero users; they are the sole reason this header includes `color/color.h` and `static_circular_buffer.h`.

173. ✅ **Dead quaternion API with an ambiguous contract** — `core/math/3dmath.h:658,1195` — `unit_inverse()` and `angle_between(Quaternion,Quaternion)`; the latter returns the 4-D half-angle and gives different results for `q` and `−q`.

174. ❌ **Six of eleven easings and two of three waves are dead** — `core/math/easing.h:25,56,63,91,101,111`, `core/math/waves.h:39,61` — flagged for a deliberate keep-or-cut decision, since these read as an intentional curve library. — rejected: keep — easing.h/waves.h are an intentional curve library; unused entries are inventory, not dead code.

175. ✅ **Stale `wrap()` references** — `core/math/geometry.h:15,320,325,345` — the header now uses `fast_wrap` or inlines the conditional.

176. ✅ **`sin_wave`'s doc overclaims interchangeability with `square_wave`** — `core/math/waves.h:18` — they are a half-cycle apart at phase 0.

177. ✅ **`sphere_log` re-declares a local epsilon** — `core/math/geometry.h:412` — duplicates `math::EPS_NORMALIZE_SQ` verbatim.

178. ✅ **`degraded_classes` counts LUTs that were later discarded** — `core/mesh/mesh_classes.h:357,399-403` — contradicting the field's own doc.

179. ✅ **`MeshState::is_bound()`'s doc overstates what it inspects** — `core/mesh/spatial.h:349-352` — it reads `vertices` only, which contradicts `set_view`'s own comment.

180. ✅ **`relax`'s `converged_out`/`iterations_out` are dead across the whole tree** — `core/mesh/conway.h:989-995,1054-1058,1093-1098` — six branches of live bookkeeping feeding nothing, inside an `HS_O3` region.

181. ✅ **Dead scaffolding in `conway.h`** — `:485` (a `ScratchScope` over an arena never allocated from), `:137-167` (`vertex_orbit`'s discarded return), `:670` (a guard made unreachable-false by the preceding `HS_CHECK`).

182. ✅ **Baked RELAX steps carry three different dead `param` values for one bake** — `core/mesh/solids.h:1365,1392,1445,1464` — unread whenever a bake is set, but they read as authored intent; dropping the bake silently picks 8 vs 217 iterations for identical geometry.

183. ✅ **A capacity trap names a knob that does not exist** — `core/mesh/hankin.h:132-134` — "MAX_INDICES" has exactly one occurrence in the tree: this string.

184. ✅ **`hankin()`'s legal `temp`-aliases-input arrangement is undocumented** — `core/mesh/solids.h:416-420` — safe only because of `ScratchScope` mark ordering; someone "fixing" the apparent aliasing would break the ping-pong.

185. ✅ **Uncovered mesh traps in the death harness** — `core/mesh/mesh.h:165,275,346,361` — the non-manifold-edge trap, the zero-side-face trap, `narrow_face_count` and `require_closed_manifold`; the zero-side case is trivially reachable.

186. ✅ **The 10k-leg walk test never exercises `seed_fix_at_start`** — `tests/test_conway_morph.h:967-1017` — and the continuity test's ADOPT gate has already drifted from production by one condition.

187. ✅ **Generated bake aggregates use positional initialization** — `tools/relax_bakes.py:98-102` — four same-typed `uint16_t` fields in a row; a struct reorder still compiles and silently mis-binds.

188. ✅ **`World::Trails::plot` lacks the zero-alpha seed gate its `Screen::Trails` sibling has** — `core/render/filter.h:872-882` vs `:1264-1283` — transparent samples consume ring slots and evict live points, visibly shortening the tail. Prefer documenting over gating (some effects may rely on it).

189. ✅ **`Canvas`'s two constructors duplicate the watchdog spin verbatim** — `core/render/canvas.h:649-658,671-680` — ten identical lines including the test-build counter.

190. ✅ **`Feedback::flush` takes two callbacks it ignores** — `core/render/filter.h:1478` — forcing `MeshFeedback` to pass a lambda that is never invoked, when the trait system already knows the stage is terminal.

191. ✅ **Branches made unreachable by always-on `HS_CHECK`** — `core/render/plot.h:1729,1813,2024-2031` — two dominated guards plus a redundant `v2` store that rewrites a bit-identical value.

192. ❌ **`rasterize()`'s eleven-parameter positional wall** — `core/render/plot.h:1108-1115` — call sites pass three consecutive `nullptr`s; the file already has the right idiom one level up in `FragmentDrawParams`. — rejected: bundling the eleven args into a struct is not codegen-neutral on the rasterizer hot path — the flag/pointer params are read inside the per-segment and per-point loops, where the opaque indirect `pipeline.plot` call forces a reload of every struct field each iteration (by-value does not help: ~28 B goes in memory under AAPCS on the M7); `FragmentDrawParams` is already deliberately exploded back into positional args at `plot.h:1484`.

193. ✅ **README's Curve-Plot table over-claims `v2` for `Plot::Line`** — `README.md:678-686` — `Line::sample` writes 0 for every sample; `Plot::Mesh` overwrites it with the edge index.

194. ✅ **README's mesh size-floor bound is not the actual bound** — `README.md:664` — the overstatement on a sliver is unbounded, not 4×.

195. ✅ **README names the wrong rasterizer for RingSpin** — `README.md:1684` — it uses `Scan::RingGroup`, the fused group rasterizer that is the whole point of its shipped design, not `Scan::Ring`.

196. ✅ **README says BZ spiral waves are "seeded periodically"** — `README.md:1556` — seeded once at init and sustained by per-substep stochastic perturbation.

197. ✅ **The segue+bank `shade_mesh_topology` overload is test-only** — `core/render/shading.h:129-143` — its "single home so the two cannot drift" rationale does not apply to a variant no effect calls.

198. ✅ **Degenerate `sides` values are admitted by four fan-shape guards** — `core/render/sdf.h:3389,3522,3681,3826` — `sides == 2` publishes `size == 0`, so a normalized shader divides by zero. Not reachable from production; hardening only.

199. ✅ **`Dynamo::dir()`'s `@note` is falsified by its own caller** — `effects/Dynamo.h:371-373` — the fractional-step accumulator routinely calls it below the claimed threshold; the conclusion still holds, the stated reason does not.

200. **The canvas-acquire IIFE is duplicated verbatim in all 24 effects** — five identical lines plus an identical comment, differing only in the profile tag. Worth an engine-side helper.

201. ✅ **`ChaoticStrings` drops the `spawn_pinned` result its sibling checks** — `effects/ChaoticStrings.h:107` vs `effects/MobiusGrid.h:65-66` — a silent null makes the transform an identity with no diagnostic.

202. ✅ **IslamicStars arena-split literals duplicate their named constants** — `effects/IslamicStars.h:57` vs `:163-167` — they agree today; a retune of either silently desyncs the pre-first-spawn split from the steady state.

203. ✅ **Dead build-chain symbols in IslamicStars** — `effects/IslamicStars.h:146,550,749` — a self-documented unused constant, a function with zero references, and a write that is never read.

204. ✅ **`strap_blend_mask_` has no width backstop** — `effects/HankinSolids.h:878` — a `uint8_t` bitmask over `NUM_PALETTES`; raising the bank past 8 truncates silently.

205. ✅ **An orphaned Doxygen block in HankinSolids** — `effects/HankinSolids.h:539-548` — attached to nothing, leaving `start_hankin_cycle` undocumented.

206. ✅ **MeshFeedback re-evaluates a constant palette every frame** — `effects/MeshFeedback.h:157` — `palette.get(0.0f)` is a virtual call costing three `fast_cosf` plus three LUT lookups for a value that never changes.

207. ✅ **MindSplatter's pool-footprint comment is stale by 64 particles** — `effects/MindSplatter.h:126-133` — says 1024 × 180 B; the cap is 1088 (~191 KiB).

208. ✅ **`MindSplatter::start_warp()` is pure indirection** — `effects/MindSplatter.h:506` — a one-line forwarder called once.

209. ✅ **`PetalFlow::spawner` is dead state after `init()`** — `effects/PetalFlow.h:144-146,171` — `Timeline::add` takes its animation by value, so the member is never touched again.

210. ❌ **`DisplacementField::step_palette_wipe()` re-implements the engine's wipe counter** — `effects/DisplacementField.h:559-565` — byte-for-byte the arming-skip + decrement half of `step_wipe_rebake`, which two other effects call. Rejected: `step_wipe_rebake` takes a `BakedPalette &` and exists to rebake it; DisplacementField has no baked LUT and samples its `GenerativePalette` directly, so calling the helper would mean adding LUT storage and a per-frame rebake this effect does not want.

211. ❌ **The aged-slot pattern is duplicated across RingShower and Thrusters** — `effects/RingShower.h:127-136`, `effects/Thrusters.h:150-159` — identical `radius_at()`/`expired()` with the same documented `age+1` convention, each pinned by its own test. Rejected: there is not much to save by factoring this out — the two differ in how they get their divisor and life bound, so a shared helper would trade four lines for a parameterization.

212. ✅ **Half of `ReactionDiffusionBase`'s protected surface is single-client** — `effects/ReactionDiffusionBase.h:120,152,248,259,316` — overstating what is actually shared. Two helpers should be private.

213. ✅ **A stale doc reference to a nonexistent `refine_nearest_node`** — `effects/ReactionDiffusionBase.h:185`, `core/engine/reaction_graph.h:132` — a reader following the contract lands nowhere.

214. ✅ **Voronoi's pixel loop does not mirror `Scan::Shader::draw` as its comment claims** — `effects/Voronoi.h:228-231` — it iterates the bare display band, not the margin-expanded clip band. Harmless today; the comment invites adding a filter and getting a seam.

215. ✅ **Flyby marks animated params by a name-lookup loop** — `effects/Flyby.h:64-66` — six repeated string literals with runtime-only typo detection, where `register_animated_param` makes the pairing unmissable.

216. ✅ **`setEffect()`'s contract omits the clip reset that `setResolution()` documents** — `targets/wasm/wasm.cpp:335-345` — real and load-bearing (a daydream worker depends on it) but documented only on the consumer side.

217. ✅ **README omits the `options` field from the parameter-definition schema** — `README.md:1885` — the enum path is live in both the bridge and the GUI.

218. ✅ **Dead accessors in `pov_sync.h` / `pov_handoff.h`** — `pov_sync.h:733,738,864`, `pov_handoff.h:196` — zero callers anywhere. *Removed `Flywheel::epoch_cycles()`, `Flywheel::cycles_per_half_rev()` and `EffectHandoff::consumed_gen()`; `BeaconParser::active()` was kept — it is asserted on by the burst-count test.*

219. ✅ **A boundary-burst drop is counted as a beacon drop** — `hardware/pov_sync.h:1578-1579` — and the associated overlap `HS_CHECK` is unreachable from its only production caller. *The drop now reports which kind of burst it discarded and an undrained boundary symbol lands in its own `boundary_bursts_dropped` counter. The overlap `HS_CHECK` was kept: it is a `SymbolEmitter` class invariant that the emitter's own tests reach directly, not just a production-path guard.*

220. ✅ **Stale "saturation clamp" language in the frame tests** — `tests/test_hd107s_frame.h:139-141,167-169` — the assertions are correct; the stated purpose describes a branch that does not exist.

221. ✅ **`HD107SFrame::load()` has no production caller** — `hardware/hd107s_frame.h:147-176` — honestly labelled a parity guard, but the path it guards is reachable only from itself.

222. ✅ **`embed()` raises `StopIteration` on a footprint with no `(attr …)`** — `hardware/phantasm/gen/pcb.py:140-142` — an opaque traceback for a plausible input.

223. ✅ **`pcb.py`'s module docstring contradicts the spec constraint it cites** — `hardware/phantasm/gen/pcb.py:5` — says ≤30 mm; the spec and README say ≤35 mm.

224. ✅ **Unchecked file I/O in the sRGB decode generator** — `scripts/generate_srgb_decode.cpp:72,96` — `fopen` never NULL-checked and `fclose`'s return ignored, so a short write still prints "wrote…" and returns 0.

225. ✅ **Test isolation is module-granular, worked around per test** — `tests/test_fixture.h:39` — eleven tests in one file hand-call `hs::random().seed(1337)`, which is evidence the per-test reset is needed and is being open-coded.

226. ✅ **Twenty near-identical `Effect` stub fixtures** — including two distinct `DeathEffect` definitions in one file, one shadowing the other.

227. ✅ **daydream: `warmModules()` never drains the fetched bodies** — `segment_controller.js:84-87` — `fetch` resolves at headers, so the 1.86 MB `.wasm` is typically neither downloaded nor cached; the "prime the module cache" half of the stated purpose does not happen.

228. ✅ **daydream: `applyResolution` reports success for an unknown preset** — `daydream.js:431` — returns `undefined`, which the transaction helper treats as applied, contradicting its own documented contract.

229. ✅ **daydream: the WASM load-failure overlay offers no recovery affordance** — `daydream.js:639-651` — a module-eval failure and a WebGL context loss both get a Reload button; the most likely real failure does not. The safe `showBootstrapFailure` helper already exists.

230. ✅ **daydream: dead CSS and one unreachable branch comment** — `styles/index.css:48-51,251`, `sidebar.js:282-283` — `.tab.active` never applies and the effect list can never receive focus, so the matching keyboard fallback is unreachable.

231. ✅ **daydream: `palettes.html` re-implements `tools/slider.js`** — `tools/palettes.html:569-699` — 130 `innerHTML`-built lines reproducing the same id/class contract without the shared module's validation; only the locked-group drag logic is genuinely page-specific.

232. **daydream: stale C++ header paths in tool JSDoc** — `tools/palette_math.js:24,50,328,403`, `tools/lissajous_math.js:28,112`, `tools/color.js:7`, plus two tests — `core/color.h`, `core/geometry.h`, `core/3dmath.h` no longer exist, and `BakedPalette::get` is attributed to the wrong header. These comments are the only pointer to the C++ that must be kept in sync.

233. **daydream: `palettes.html` comment and import hygiene** — `:441-443,491-494,533-535` claim OKLab math lives in `tools/color.js` (it is entirely engine-side); `:453` and `:459` import two unused symbols; `:180-211` places a `<style>` block inside `<body>`.

234. **daydream: a stray build probe at the repo root** — `sect.cpp`, a two-line `__attribute__((section(".flashmem")))` compiler probe in a static web app; `prompts/` and `.claude/` are also untracked and, unlike Holosphere, not ignored.

235. **Screenshots are ~50× larger than their display size** — `docs/screenshots/*.png` — 24 PNGs at 1360×1146 (up to 1.42 MB each, ~28 MB total) embedded at `width="280"`, and installed into daydream so both repos carry the payload.

---

## 4. What is exemplary

Recorded so it is not lost in the defect list. These are the patterns worth propagating, not merely
"things that work":

- **Comments as evidence.** `-funswitch-loops` cited with real ITCM byte deltas; `hs::shuffle` naming
  the three libraries that diverged; the Hankin parallel gate recording the *measured* separation
  between healthy far stars and genuine yeets. This is the rare codebase where the comments are
  load-bearing data.
- **Invariants as compile errors.** The four filter-ordering `static_assert`s name the offending
  filter *and* the fix. `PipelineRef`'s `static_assert(!requires { T::any_crosses_segments; })`
  turns a signature drift that would silently downgrade culling into a build failure.
  `ConwayGraph` proves its whole edge table at compile time.
- **Gates on gates.** `--check-modules` pins CMake↔runner in both directions; `check_includes.cmake`
  catches orphaned includes; `teensy-gate-tests` proves the size gate *fails* on a broken fixture;
  the ccache job fails when zero cacheable calls are recorded. Finding 1 is a hole the design
  already anticipates.
- **Anti-lockstep-drift testing.** `segment_crosscheck.test.js` pins its own hand port against
  independent literals so port and subject cannot drift together; `color_parity_wasm.test.js` pins
  OKLab red to Ottosson's published value rather than to the engine.
- **The death harness.** Shell-free `fork`/`spawn`, a dedicated always-trapping sentinel to probe the
  relay shape, strict single-shape matching so a genuine `exit(132)` cannot pass, `MIN_DEATH_CASES`
  as an anti-shrink floor, hard-fail under CI and quiet-skip locally.
- **Flag-sensitive contracts get their own TU.** `fastmath_clamp_check` recompiles the NaN-clamp
  assertions under the shipping `-ffast-math`; `h_offset_renorm_check` exists purely to dodge an ODR
  clash so the device's `H_OFFSET > 0` path is exercised at all.
- **The pure/platform split.** `pov_sync.h`, `pov_segment_map.h`, `pov_handoff.h`, `param_marshal.h`
  and `wasm_predicates.h` push every load-bearing decision into host-testable code, and the 4-board
  event-driven sync simulator with ppm skew, EMI injection and mid-show reboot is unusually rigorous
  for embedded firmware.
- **Provenance discipline.** Five generated artifacts each carry a banner naming their generator;
  two were re-derived bit-identically during this review; the `.wasm` ships with a sha and a
  `-dirty` marker that CI asserts is absent on a clean checkout.
