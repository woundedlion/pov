# Holosphere / Daydream — Comprehensive Code Quality Review

**Date:** 2026-06-22 (fifteenth audit; supersedes the fourteenth. That pass's fix list 1–127 is essentially fully resolved in the tree; its two open carry-overs and its standing rejection are re-confirmed and re-numbered here, and §5 reconciles the rest.)
**Tree:** Holosphere `master`; daydream reviewed at its current `master`.
**Scope:** the entire two-repo product — the Holosphere C++ engine + firmware (`core/`, `effects/`, `hardware/`, `targets/wasm`, `targets/Phantasm`, `tests/`, `scripts/`, `tools/`, build) and the daydream web simulator (`*.js`, `tools/*.js`). Out of scope per instruction: `core/effects_legacy.h`, the `*.ino` target sketches **except** `targets/Phantasm/Phantasm.ino` (in scope), and `core/rotate.h`. `core/FastNoiseLite.h` (vendored MIT) and the generated data tables (`core/color_luts.h`, `core/reaction_graph.cpp`) were assessed for integration/provenance only.
**Method:** the 2,100-line README was read in full for architecture, then fifteen independent sub-agents each audited a coherent component — reading every in-scope file with the README as shared context and grading objective quality (correctness, UB, ISR/memory-ordering, bounds, overflow, performance, testability) and subjective quality (architectural elegance, interface expressiveness, idiom). The prior `CODE_REVIEW.md` was **not** consulted during examination; it was reconciled only afterward (§5).

---

## 1. Executive Summary

Holosphere is a persistence-of-vision LED sphere rendering engine and its bit-identical
WebAssembly simulator. One header-only C++17 engine — every rendering class templated on
`<int W, int H>` — compiles to Teensy 4.0 firmware (single-board Holosphere and four-board
synchronized Phantasm) and, via Emscripten, to the browser simulator wrapped in Three.js.

Read cold against the current tree by fifteen reviewers spanning both repositories, the
product again holds the high bar its audit series has established. **Across well over a
hundred source files, not a single Critical- or High-severity correctness, memory-safety, or
concurrency defect was substantiated.** As in every prior pass, the recurring observation was
the *absence* of the defects one expects: the zero-copy WASM view-detachment hazard is
defended at every consumption site (and the subtle `instanceColor.array = null`-before-dispose
ordering is done right); the multi-board sync protocol's failure mode is provably "missed,
never wrong"; arena bounds math is wrap-safe and trap-guarded; numeric float→int casts are
uniformly NaN-clamped; and `HS_CHECK` fail-fast lands exactly where the doctrine prescribes
(cold seams) and is withheld from the per-pixel hot loop by design. Every C++-core reviewer
independently reported the same signature: the residual defects are overwhelmingly **internal
inconsistencies with the codebase's own exacting standards** — the one implicit `CRGB` cast
where its sibling is `explicit`, the one `Spiral` sampler that skips the `v2` register every
other sampler sets, the one comment copied from a sibling effect that no longer matches its
file. When a codebase's flaws are best described as exceptions to its own rules, the rules
have largely won.

The findings concentrate in five bands, none load-bearing for correctness today:

1. **Two firmware quality gates are inert in the running configuration.** The Teensy
   warning-hygiene ratchet (`teensy_warnings.py` + a committed baseline) is wired into nothing
   — it appears only in a comment — and the device size/layout gate is `if:`-disabled by default
   because Phantasm overflows RAM1. Combined with the project's "no device build in CI"
   constraint, **no automated build coverage exercises the `#ifdef ARDUINO` device path on PRs**,
   and the README's "the size gate compensates" claim is not yet true as configured.
2. **Test-coverage breadth.** The C++ suite's two headline geometry transformers
   (`mobius_transform`, `gnomonic_mobius_transform`) and `OrientTransformer` are exercised
   **only at the identity map/quaternion**, so a coefficient-ignoring regression passes; and
   ~22/28 effects remain smoke + determinism only (carry-over). The JS browser-glue modules
   reviewed here are correct but the prior audit's untested-glue High items stand.
3. **New documentation drift introduced by prior fixes** — chiefly the README still describing
   `AABB` as living in `spatial.h` after the fourteenth audit *removed* it (Finding 6), and the
   Feedback-style table listing 7 presets where the code now has 8.
4. **Simulator/hardware parity gaps in the segmented-POV layer** — the JS `computeSegmentRange`
   accepts segment counts the Phantasm hardware map cannot represent and does not model the
   reversed bottom-strip direction, partly undercutting its "exercised before fabrication" goal.
5. **Comment/naming/consistency polish and a thin band of latent edges** — an unguarded
   `Warp::Twist` Lipschitz divisor, a few unenforced unit-vector preconditions on public
   functions, and a scattering of cross-effect idiom divergences.

**The grade is A−.** As in every prior audit, the distance to a flat A lives almost entirely
at the **gate/coverage frontier** — not the test *infrastructure*, which is A+ (a forked death
harness asserting every `HS_CHECK` fires with the exact illegal-instruction signal; a
fault-injecting four-board sync simulator; injected-clock cross-run determinism; analytic
oracles over Euler characteristic, KD-tree brute force, and SDF cull-conservativeness) — but
the *running configuration* of the firmware gates and the behavioral coverage of the most
delicate paths. The engine, color, memory, mesh, hardware-protocol, WASM-bridge, and
test-infrastructure subsystems each grade A or above on their own terms.

---

## 2. Quality Dimensions — Letter Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Architectural elegance** | **A+** | Compile-time `<W,H>` specialization, the World/Screen/Pixel filter pipeline with automatic compile-time domain conversion + fold-`static_assert`ed ordering invariants, the state-mutation-vs-render split in the animation system, the SDF/scan/plot split, the CRTP reaction-diffusion base, and the three-layer flywheel sync are genuine design contributions. Subsystems compose rather than entangle. |
| **Interface expressiveness / API design** | **A** | `StaticPalette` compile-time modifier chains, `SolidBuilder` fluent Conway chaining, the `generate()` arena wrapper, `FunctionRef`/`inplace_function` zero-alloc callables, deleted rvalue overloads enforcing borrow contracts, `worker_protocol.js` as a typed JS↔worker contract, and `registerParam()` reflection are expressive and hard to misuse. Minor sharp edges (`Color4`'s implicit `CRGB` cast vs `Pixel16`'s `explicit`; `OrientTransformer`'s dead `W` parameter; the composed-Conway return-arena polarity footgun for direct callers). |
| **Correctness & robustness** | **A−** | No substantiated Critical/High bug across 100+ files; the degenerate-geometry policy (`normalized()` trap vs `normalized_or()` soft-degrade) is principled and applied per call site. Held below A by a scattering of latent edges (the `Warp::Twist` `R==0` divisor, the `SmoothUnion` register seam, unenforced unit-vector preconditions on public projection functions) whose siblings guard. |
| **Memory safety** | **A** | Arena bounds math is wrap-safe (subtractive no-wrap checks + independent multiply-overflow trap) and max-aligned; dual-stamp use-after-free detection on `ArenaVector`/`ArenaSpan`; deleted rvalue ctors prevent dangling borrows; zero per-frame heap on any path. The `append_bulk`/`back()` wrap edges are latent (no caller reaches them). |
| **Concurrency (ISR / sync / workers)** | **A** | Single-writer discipline, correct relaxed-atomic reasoning for the single-core M7, position-from-time flywheel, airtight Transferable ownership and generation-fencing in the worker pool. The `read_id` strap-settle window is the one heuristic-timed seam. |
| **Numerical robustness** | **A−** | Pole/antipode/normalization edges handled deliberately; fast-math approximations carry measured error budgets with exact-trig escapes. Held below A by multiple independently-tuned "is-unit" epsilons the `math::EPS_*` block was meant to unify, and the one unguarded `Warp::Twist` divisor. |
| **Error handling / fail-fast** | **A** | A faithful application of the doctrine: cold-path invariant violations trap via stdio-free `HS_CHECK` (survives `NDEBUG`, death-harness-verified); transient conditions get bounded soft handling; the WASM/JS boundary logs-and-rejects untrusted input rather than trapping (one gap: `PaletteOps::bakeLut`'s unvalidated `gradientShape`). |
| **Performance** | **A−** | LUTs, interval culling, branchless ISRs, packed-SIMD blits, hot/cold partitioning, zero per-frame heap. A few effects carry slider-gated cost cliffs flagged honestly in-comment (`ShapeShifter` ~256 raster passes/frame is the standing example); the per-worker per-frame `Uint16Array` alloc and `nearest()`'s O(k)-per-node worst rescan are minor. |
| **Readability & naming** | **A−** | Consistently descriptive; comments are load-bearing, not decorative. A few terse names survive and comment *volume* occasionally dwarfs the code it annotates (several effects are >60% comments). |
| **Comments & documentation** | **A** | The README is a reference-class architecture document, accurate to code in almost every particular. Inline Doxygen is near-universal. Held at A by a handful of drifts — the `AABB`/`spatial.h` reference the prior fix invalidated, the 7-vs-8 Feedback presets, the `Catalan`/registry prose. |
| **Testability** | **A−** | C++: outstanding (oracles, death harness, determinism pass, tiling proofs, parity TUs). JS: pure modules well-covered, browser glue not. |
| **Test coverage** | **B+** | The grade-limiter, as in every prior audit. C++ suite is **A** (~2,600 oracle-driven assertions) but with a real hole — the two Möbius transformers and `OrientTransformer` are identity-only — plus the behaviorally-blind effect smoke roster and untested JS render glue. |
| **Build / CI / tooling** | **B+** | The WASM/native CI is genuinely excellent (sharded suite, ASan/UBSan, shard-coverage anti-drift guard, runtime WASM smoke, install-provenance trio). Pulled down from A− this pass by two inert firmware gates (warning ratchet unwired; size gate off by default) and a non-hermetic install-time `git` dependency. |
| **Consistency** | **A** | Idioms (arena params, `HS_CHECK`/`assert` split, `Filter::` namespacing, narrowing casts, X-macro rosters, Doxygen) are uniform across a large surface; the divergences found are individually minor (param-registration idiom, timer-lambda signatures). |
| **Portability (Teensy / WASM / host)** | **A** | `platform.h` cleanly abstracts the targets with a documented divergence ledger; host-testable cores split out of Arduino-only shells; the NaN-clamp contract is `#error`-guarded; deterministic seeded RNG gives bit-identical multi-target renders. |
| **OVERALL** | **A−** | Code quality is A/A+ across the engine; the composite is held at A− by the gate/coverage frontier (two inert firmware gates, the identity-only transformer tests, ~22/28 behaviorally-unpinned effects, untested JS glue) and a thin band of latent hazards. |

---

## 3. Notable Strengths

- **Provably-safe multi-board synchronization.** `pov_sync.h` derives column position from a
  free-running cycle counter (not interrupt counting), making masked-IRQ windows
  correctness-neutral by construction, paired with an odd-only distance-2 symbol codec whose
  worst case is a *missed* symbol (self-healed next half-rev), never a *misclassified* one. The
  same-tick burst/fold `j`-inference was independently re-verified correct this pass.
- **Fail-fast that is itself tested.** `HS_CHECK` traps are verified by a cross-process death
  harness that asserts the exact `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` status and hard-fails CI
  if unrunnable — the safety net is proven, not assumed.
- **The hardest WASM hazard, handled everywhere.** The zero-copy `getPixels()` view detaches on
  heap growth; every consumption path re-fetches, the three-way buffer alias is asserted on both
  the single-engine and segmented composite paths, and `instanceColor.array` is nulled *before*
  `dispose()` so Three.js never re-uploads freed WASM memory across a resolution change.
- **Numeric discipline.** 16-bit linear-light end-to-end, OKLCH perceptual palette interpolation
  with chroma-preserving gamut clipping, round-to-nearest fixed-point, and UB-safe float→int
  casts (every cast preceded by a NaN-folding clamp). Generated LUTs round-trip losslessly and
  were spot-checked bit-for-bit against a regenerated double-precision reference.
- **Memory model as a design feature.** A single 335 KB arena, repartitionable per-effect, with
  explicit `Arena&` parameters at every call site, RAII scratch scopes, `Persist<T>` evacuation,
  and dual-generation use-after-free detection.
- **Parity-by-delegation in the tools.** `splines.html`/`solids.html` and the generative palette
  tool route their math through the *same WASM engine code* rather than reimplementing it;
  `lissajous_math.js`'s figure math was verified bit-faithful to `geometry.h`.
- **A reference-class README** documenting not just the *what* but the *why* of every major
  decision, down to AC timing tables for the 1-wire sync protocol.

---

## 4. Prioritized Fix List

Every actionable defect found this pass, grouped by severity and numbered in a single global
sequence (matching the codebase's finding-number convention). Each entry: `file:line — problem
→ fix`. Items independently re-confirmed open from a prior audit are marked **(carry-over)**.
A leading ❌ marks a deliberately-rejected change preserved so it is not re-raised. Findings that
this pass re-discovered but the tree already resolves are not listed here as open — they are
catalogued in §5.

### Critical

*No Critical findings were substantiated. No correctness, memory-safety, concurrency, or
data-loss defect was found anywhere in scope.*

### High

1. `tools/teensy_warnings.py` + `tools/teensy_warning_baseline.txt` — the entire Teensy
   warning-hygiene ratchet is **dead code**: it is referenced only in a `platformio.ini:33`
   comment and the spec, appearing in no `extra_scripts` entry and no CI step, so new firmware
   warnings are gated by nothing. (Sibling tools `teensy_map.py`/`teensy_nano.py`/`teensy_pre.py`/`teensy_gate_extra.py` *are* wired; only the ratchet is orphaned.) → Wire it into the PlatformIO `extra_scripts` and/or the `teensy-size` CI job, or delete it and the baseline if the ratchet is abandoned.
2. `.github/workflows/ci.yml` (teensy-size job) — the only automated device-compile/size gate is `if: vars.TEENSY_GATE_ENABLED == 'true'`, **off by default**, because Phantasm overflows RAM1. With the project's "no device build in CI" constraint this means no automated coverage of the `#ifdef ARDUINO` device path runs on PRs; a device-only compile break rides a green build. → Track the Phantasm RAM reduction to re-enable the gate; in the interim, run at least the Holosphere image (which fits and passes) unconditionally so the device path is compiled on every PR.
3. ✅ `tests/test_transformers.h` — `mobius_transform` and `gnomonic_mobius_transform` are tested **only at the identity map** (a=1,b=0,c=0,d=1) and `OrientTransformer` only at the identity quaternion; a transform that ignored its coefficients entirely would pass every round-trip. These are the suite's two headline geometry ops, effectively unverified against non-trivial parameters. → Add non-identity coefficient/quaternion cases with an independent oracle (a hand-computed Möbius image and a known rotation).
4. `tests/test_effects.h` — **(carry-over, prior #4/#10)** ~22 of 28 effects are covered only by `smoke_one` (no-crash) + `determinism_one` (reproducible); both are behavior-blind, so a wrong-but-non-black, reproducible regression passes. Only ~6 effects have behavioral white-box oracles. → Add one pinned numeric property each for the highest-risk untested effects (Raymarch, HopfFibration, MobiusGrid, SphericalHarmonics, Voronoi).
5. ✅ `daydream/driver.js`, `daydream/daydream.js`, `daydream/recorder.js`, `daydream/sidebar.js` — **(carry-over, prior #1–3)** the core render/Three.js loop, the main orchestrator, the `MediaRecorder` session lifecycle, and the sidebar widget remain essentially untested; this pass's code review found them correct but their delicate paths (memory-view refresh, the buffer-alias contract, the stop→start stale-session race) are pinned by nothing. → Extract the DOM-free logic and unit-test it behind stubbed Three.js / a fake `MediaRecorder`.

### Medium

6. ✅ `README.md:1059` — §7.7 states `spatial.h` contains "a speculative `AABB` (test-only — no production consumer yet)", but the fourteenth audit (its Finding 56) **removed** `AABB` from `spatial.h`; the only surviving `AABB` is an unrelated generic doc-comment in `constants.h`. Stale reference introduced by a prior fix. → Delete the `AABB` clause from the §7.7 file list.
7. ✅ `core/sdf.h:939-947` (`SmoothUnion::distance`) — the smooth-min blends only `dist`; the surviving `DistanceResult` (`t`/`raw_dist`/`size`) snaps to whichever child is nearer, so a shader keying off `v0`/`size` sees a hard seam through the visually-smooth weld. Undocumented (unlike the AngularRepeat UV note). → Documented the register discontinuity at the site (chose the doc over blending: blending every register costs a second set of lerps on the SDF hot path for a feature no shipped geometry reads — a perf regression not worth taking).
8. `core/sdf.h:2944-2955` (`Warp::Twist::lipschitz`) — `gamma` divides by `std::max(s, R*0.5f)`; a `Twist` constructed with `R == 0` (degenerate torus) floors the divisor to 0 → div-by-zero `gamma`. No shipped geometry hits it (Raymarch passes the torus `R`), but unlike `Star`/`AngularRepeat` there is no `HS_CHECK(R > 0)` guarding the degenerate input. → Add the cold-path `HS_CHECK(R > 0)`.
9. `hardware/pov_segmented.h:427,430` (`read_id`) — the 10 ms/2 ms strap settle/resample delays are heuristic; a strap settling slower than 2 ms can read identically at both samples and pass while mid-transition, mis-IDing a board into a *second* master driving the push-pull sync wire into bus contention — a fault no single board can detect. → Lengthen/repeat the debounce, or add a post-boot cross-check (e.g. detect two boards asserting master).
10. ❌ `daydream/segment_layout.js:29` (`computeSegmentRange`) — validates `total` only as positive-even, but the hardware map (`pov_segment_map.h:44`) requires N an even power-of-two ≤ 4 (it masks `segment_id & (segs_per_arm-1)`); a `total` of 6 passes and renders a 3-way Y-split per arm that has **no hardware analog**, defeating the "exercised before fabrication" purpose. → **Rejected (validated against intent):** daydream's segmentation is a *deliberately general* even-N parallel preview, not a hardware-exact mirror — the `Segments` GUI offers 2/4/6/8 (`daydream.js:683`, step 2) and `segment_layout.test.js` pins 6 & 8 as valid tilings. `computeSegmentRange` uses real division/modulo (not the hardware's power-of-two mask), so it tiles any even N correctly. Restricting it to {2,4} would break that tested feature. The hardware-fidelity gap is instead addressed honestly in the README (finding 11).
11. ✅ `daydream/segment_layout.js:46-56` — the JS layout models only an axis-aligned `[y0,y1)` rectangle and never models the bottom segment's *reversed* strip direction (`y_step=-1`, LED 0 at the S pole) that `pov_segment_map.h` warns "silently mis-paints the sphere"; the highest-risk hardware index math is the part not simulated, so the README's "mirrors that math" overstates it. → Softened the README claim (both repos) to scope it to the arm/Y-band partition, calling out that the strip-direction wiring and the power-of-two segment-count constraint are *not* simulated. (Chose the doc fix over simulating strip direction: per the finding-10 decision, daydream is a deliberately general even-N preview, not a hardware-exact mirror.)
12. ❌ `daydream/geometry.js:34-39` (`pixelToSpherical`) — maps latitude with `y·π/(H−1)` while the engine's `pixel_to_vector` uses `y/H` semantics; the file's own comment concedes this aligns only because `H_OFFSET==0` and is "an unenforced assumption," so dot positions diverge subtly from the true engine mapping at the poles. → **Rejected (false premise):** the engine does **not** use `y/H` semantics — `y_to_phi` is `(y·π)/(H_VIRT−1)` (`core/geometry.h:318`, with `H_VIRT = H + H_OFFSET`), and `vector_to_pixel`'s inverse `phi_to_y` matches. With the simulator's `H_OFFSET == 0`, `H_VIRT−1 == H−1`, so the JS denominator is *exactly* the engine's; the dot positions coincide (no pole divergence). The only residual is the device `H_OFFSET == 3` build, which the JS already documents as an unenforced assumption and which the simulator never runs. No change.
13. ✅ `daydream/recorder.js` (`_saveWithAnchor`) — setting `iframe.src = blobUrl` to a **video** blob can start inline playback or, for a `Content-Disposition` download, never fire `load`, so revoke falls back to the 60 s timeout each time; the anchor `a.click()` already triggers the download, making the iframe redundant and capable of double-triggering. → Dropped the iframe; revoke the object URL on a short (1 s) timeout after the anchor click. (Landed in the daydream repo.)
14. `tests/test_mesh.h:~348,367` — the two uniform-topology classifier tests assert only that all class counts are zero, which (the file itself notes) passes even if the classifier were broken; there is no non-manifold/boundary-mesh fixture, so the half-edge `pair == HE_NONE` path is unexercised. → Add a truncated/boundary fixture asserting a non-trivial class split and a non-manifold guard case.
15. ❌ `core/conway.h:46` (`face_centroid`) — unlike its sibling `face_normal` (which guards `he.next != HE_NONE` inside the walk), `face_centroid` advances `he_idx` and reads `mesh.vertices[...]` without the in-loop `HE_NONE` guard, so a corrupt `.next` indexes `half_edges[HE_NONE]` once before the anti-hang count catches it. → **Rejected (non-issue):** `face_centroid` only single-steps — it reads `half_edges[he_idx]` (always live: the first index is outer-`if`-guarded, every later one is gated by the `do…while (he_idx != HE_NONE …)` post-condition), then sets `he_idx = …next`. It never dereferences `.next` *in-body*, so `half_edges[HE_NONE]` is never indexed. `face_normal` needs the extra guard precisely because it *does* read `half_edges[he.next]` in-body (two-step lookahead). The asymmetry is correct and already documented at `conway.h:73-77`. No change.
16. `effects/{Liquid2D,Flyby,Raymarch}.h` — the preset-cycle + wrapped-trig-phase boilerplate (`fmodf(phase+δ, 2π)` phase wrapping, the `markAnimated`-over-a-name-list loop, and the `RandomTimer/PeriodicTimer → presets.next() → Lerp(params, prev_get(), get(), …)` sequence) is near-identical across these effects with small per-effect drift. → Extract a small shared helper for the phase-wrap and the preset-cycle-timer construction (without coupling the shaders).
17. `targets/wasm/wasm.cpp:1185-1198` (`PaletteOps::bakeLut`) — the only JS-boundary entry point here that does **not** validate its input: `gradientShape` is `static_cast` into an enum with no range check (UB on an out-of-range int), and the h/s/v ints are silently truncated to `uint8_t`. → Validate `gradientShape` against the enum range and reject out-of-band like the sibling bridge methods.
18. `CMakeLists.txt:128-164` (install CODE block) — the provenance `.sha`/`.wasm.sha256` generation runs `git rev-parse`/`status` at **install time** via `execute_process`, so `just install` `FATAL_ERROR`s in an exported tarball / non-git checkout. → Make the git stamp best-effort (warn + omit the marker) when `.git`/`git` is absent, so install works outside a working tree.
19. `tools/teensy_budgets.json:57-67` & `tools/teensy_gate_tests/fixtures/` — the Phantasm symbol/layout budgets and the `good_*`/`broken_*` gate fixtures are hand-built/synthetic (only `fixtures/real/` holds genuine Holosphere captures), so the gate's correctness for Phantasm is unproven against real toolchain output. **(carry-over, prior #107.)** → Capture a real Phantasm ELF once the image fits and replace the synthetic fixtures.

### Low

20. `core/color.h:317` — `Color4::operator CRGB()` is **implicit**, contradicting the deliberately `explicit` `Pixel16::operator CRGB()` whose doc says explicitness exists to block silent 8-bit gamma round-trips; this path reintroduces that footgun and has no core-sink caller. → Make it `explicit`.
21. `core/color.h:2311` / `palettes.h:112` — `MeshPaletteBank::N` (5) is duplicated against `BakedPaletteBank::N` and the hardcoded 5-pointer `sources()` list with no `static_assert` tying `sources().size()` to `N`; a future edit to either side silently desyncs. → Add a `static_assert(sources().size() == N)`.
22. `core/3dmath.h:495,1057,1186` — three different "is-unit"/"is-degenerate" tolerances coexist (`0.01f` in `unit_inverse`, `0.02f` in `make_basis`/`make_rotation`, raw `0.9999f`/`PI_F-0.0001f` in `slerp`), defeating the centralized `math::EPS_*` block intended to unify them. → Promote to named `EPS_UNIT_SQ`/length-epsilon constants.
23. `core/geometry.h:466,788` — `vector_to_pixel` and `make_basis` document but do not enforce their unit-/non-zero-vector preconditions on public templated functions; a non-unit `v` yields a silently-wrong `fast_acos(v.y)` phi, and a zero `normal` traps in `make_basis` without the graceful path its `normalized_or` siblings offer. → Add a cold `HS_CHECK` or a `normalized_or` fallback at each.
24. `core/geometry.h:452` — `pixel_to_vector`'s analytic fallback feeds `y` to the fractional `y_to_phi<H>` with no clamp, so a sub-pixel `y` past `H_VIRT-1` extrapolates phi outside `[0,π]` while the integer overload traps; an asymmetry easy to trip. → Document the unchecked contract at this call or clamp.
25. `core/3dmath.h:110` — `explicit constexpr Vector(int)` (an `inplace_function` compat tag) compiles `Vector v(5)` to the *zero* vector, discarding the argument — a real footgun for a math type, mitigated only by `explicit`. → Replace with a named tag type.
26. `core/memory.h:534` (`append_bulk`) — the guard `size_ + count <= capacity_` can wrap for a colossal `size_t count`, letting an overflowed sum slip under `capacity_`; `Arena::allocate` already uses the wrap-proof subtractive form. → Rewrite as `count <= capacity_ - size_`.
27. `core/static_circular_buffer.h:237` (`back()`) — `(head + count - 1) % N`'s intermediate `head + count` can reach `2*(N-1)`, which wraps for an `N` near `UINT32_MAX`; latent only (no instantiation is near the bound). → Note the `N` bound or compute wrap-safely.
28. `core/memory.h:888` (`ScratchScope::make_vector`) — returns an `ArenaVector` by value whose backing is freed on scope rewind; a plain `set_offset` rewind (no generation bump) is **not** caught by the debug use-after-free stamp, so a returned-and-outlived vector dangles even in debug. → Document the in-scope-only contract at the API, or bump the generation on rewind.
29. `core/sdf.h:1421` (`Spiral::sample`) — does not set `v2` (vertex index) while every sibling ring/star/flower sampler does; a shader relying on the documented `v2` convention reads a stale Fragment register for spirals. → Set `v2 = index`.
30. `core/sdf.h:1322-1324` (`AngularRepeat::get_vertical_bounds`) — any non-Y fold axis returns the full canvas band, fully defeating vertical culling for tilted repeats (every row scanned). → Compute a tighter analytic band for the general axis, or document the perf cliff at the call.
31. `core/filter.h:1403` (`Pixel::Feedback::flush`) — `bool col_used[W] = {}` is a per-flush stack array up to `MAX_W` (288) bytes — the one unbounded-by-template stack allocation on a per-frame path. → Use a `std::bitset<W>` or an arena slice.
32. `core/transformers.h:188-214` — `transform()`/`operator()` are `const` but correctness depends on the non-const `prepare_frame()` having run this frame, with no compile-time or assert linkage; a caller that forgets it renders stale state silently. → Add a per-frame "prepared" assert, or fold preparation into the call.
33. `core/transformers.h:228` (`OrientTransformer`) — the `template <int W>` parameter is entirely unused ("carried for interface symmetry"), forcing callers to specify a meaningless `W` and bloating instantiations. → Drop the parameter or give it a real role.
34. `core/animation.h:582` (`ParticleSystem::spawn`) — silently drops the spawn on a full pool with no `hs::log`/trap, diverging from the file's convention (cf. `Timeline::add`). → Log the soft exhaustion for parity.
35. `core/spatial.h:110-179` (`nearest`) — `get_worst_dist()` does a full O(k) linear rescan of `best` on every visited node's prune test (not just on displacement); a cached worst, invalidated on displacement, would be O(1) per node. → Cache the running worst. (Minor; k≤5 today.)
36. `effects/PetalFlow.h:264-275` — the ring is closed twice: the loop appends a duplicate of `fragments[0]` *and* `Plot::rasterize(..., close_loop=true)` reconnects last→first, yielding a redundant zero-length closing segment. → Drop the manual duplicate or pass `close_loop=false`.
37. `effects/Voronoi.h:311` — `float goldenAngle = PI_F*(3.0f - sqrtf(5.0f))` is recomputed inside the per-site `seed_sites()` loop; it is loop-invariant. → Hoist or `constexpr`.
38. `effects/SphericalHarmonics.h:376-377` — `current_idx`/`next_idx` lack in-class initializers (set in `init()` to 6); every other member here carries a default. → Add `= 0` for defense-in-depth.
39. `effects/Thrusters.h:200-209` — `on_fire_thruster` snapshots `phase = warp_phase` into `r_fn` but the two `fn_point(...)` calls pass the member `warp_phase` directly; equal today, but contradicts the file's own "snapshot so geometry can't depend on member-mutation order" rationale. → Pass the snapshot.
40. `effects/DistortedRing.h:54` — `registerParam("Thickness", …, 2.0f*px, 0.75f)` mixes a resolution-scaled min with a hardcoded literal max, so the slider's meaning silently changes across resolutions while the comment only justifies the scaled min/default. → Scale the max too, or document the fixed cap.
41. **Cross-cutting (effects) — parameter-registration idiom split:** `Moire.h:68` uses `registerAnimatedParam`, while `IslamicStars`/`HankinSolids:50`/`MobiusGrid:62` use `registerParam(...)` + `markAnimated(...)` for the identical intent. → Pick one idiom.
42. **Cross-cutting (effects) — timer-callback lambda signature split:** some `PeriodicTimer` lambdas type the canvas as `Canvas&` (`IslamicStars:165`, `DreamBalls:83`) and others as `auto&` (`Moire:71`, `MobiusGrid:70`), even within one file. → Standardize.
43. `hardware/dma_led.h:104-114` (`transmitAsync`) — retains a `waitComplete()` spin unreachable from the only caller (`submitFrame` drops on overrun first); its presence invites a future caller to spin in the column ISR — the exact deadlock the drop path avoids. → Remove it or assert "not in ISR context."
44. `hardware/pov_single.h:70-87` — the non-DMA FastLED branch calls `setCorrection`/`setTemperature` but not `setBrightness(255)` while the DMA branch sets all three; behaviorally equivalent (FastLED defaults 255) but an unexplained asymmetry. → Set all three in both branches or comment the difference.
45. `hardware/pov_sync.h:1547` (`maybe_schedule_beacon`) — the beacon start has no late-censor on its *position* (only the per-pulse `late_censor` inside `emitter_.tick`), so a long masked window pushing the first post-W/4 wake toward W/2 could let the ~26 ms beacon frame encroach toward the HALF boundary; safe under the timing budget but unguarded unlike the boundary emitter. → Add a beacon-start late-bound or document the budget margin at the site.
46. `hardware/hd107s_frame.h:35` / `hardware/dma_led.h` — both advertise as independently host-testable yet use `HS_CHECK`/`FASTRUN`/`hs::log` without including `platform.h`/`led.h`, relying on the includer pulling them first. → Add the direct includes or document the include-order requirement in the header.
47. `targets/wasm/wasm.cpp:919-929` (`classifyFaces`) — the only mesh getter returning a typed array without the read-before-next-call caveat its siblings (`getVertices`/`bakeLut`) carry, even though the same tooling-arena lifetime applies; the `.new_(Int32Array)(view)` form copies today, so it is safe but the doc gap is a latent trap if changed to a raw view. → Add the contract note.
48. `targets/Phantasm/Phantasm.ino:78` — `Serial.begin(9600); delay(1000)` blocks setup a full second waiting on a console that won't be attached on the installed art piece, contrary to the "no console on hardware" philosophy. → Drop or shorten the delay.
49. `daydream/segment_controller.js:467` — the composite final-guard `if (r.x0<0 || … || r.x1>w || r.y1>h) continue;` silently drops a whole segment, inconsistent with the loud `throw` the same function uses for the alias-break a few lines down (and compounds Finding 10). → Log/throw on an out-of-range rect.
50. `daydream/segment_worker.js:165` — each `render` allocates a fresh `new Uint16Array(qw*qh*3)` then transfers it, the per-frame GC churn the controller side works to avoid; a double-buffered reusable pair per worker would remove it. → Reuse two transferable buffers per worker.
51. `daydream/tools/palette_math.js:308` — `calcHues`' `switch` uses `case "ANALOGOUS": default:`, so an unknown `harmonyType` silently renders as analogous, while the export path (`generativePaletteCpp:419`) *rejects* unknown tokens — preview and export disagree on validation. → Throw on an unrecognized harmony in both.
52. `daydream/tools/mobius_transforms.js` — the `hyperbolic`/`loxodromic` presets compute `flowParam = (t*speed) % logPeriod`, which makes the warp scale discontinuously snap from `exp(logPeriod)` back to `exp(0)` each period — a visible pop for presets the docstrings call "seamless"/"continuous." → Triangle-wrap the flow param or correct the "seamless" wording (verify against the shader's intent first).
53. `daydream/driver.js:336-369` (`setCanvasSize`) — called from the constructor before the resize observer exists; a 0×0 initial container yields `aspect = 0/0 = NaN` → a NaN projection matrix until the next resize, with no finite-guard (unlike the recorder's aspect clamp). → Guard against a zero container size.
54. `daydream/driver.js:148,168` — `this.canvas.parentElement` is dereferenced unconditionally (`.appendChild`, `.observe`); a parentless `#canvas` throws and white-screens, contrary to the graceful-degrade doctrine applied to `#gui-container`. → Null-check the parent.
55. `daydream/state.js:159-162` — `URLSync` subscribes to `AppState` but discards the unsubscribe handle and exposes no `dispose()`; `disposeApp()` never tears it down, so a `pagehide` discard can leave a 200 ms timer firing `history.replaceState` into a dead page. → Return/store the unsubscribe and clear the pending timer on dispose.
56. `daydream/recorder.js:233-237` — the scaled-capture path calls `_ensureOffscreen()` per frame and resizes the offscreen on a mid-recording window resize, changing the captured track width — contradicting the class doc's "track frame size fixed for the whole session" invariant (which encoders require). → Pin the offscreen size at record start.
57. `daydream/sidebar.js:164` — `scrollIntoView({behavior:'smooth'})` runs on every `setActive`, including programmatic selection during init/resolution-switch and every "Test All" tick, animating an unwanted scroll. → Use `behavior:'auto'` for programmatic selection.
58. `daydream/daydream.js:407,621` — the mobile-collapse heuristic reads `window.innerWidth < 900` while the renderer uses the driver's container-width `isMobile`; on a narrow container in a wide window the two disagree. → Drive both from one source.
59. `daydream/vendor-importmap.js:29` — depends on `document.currentScript` (null inside ES modules) with no guard, so a future `type="module"` load throws an opaque "Cannot read 'src' of null" rather than a named error (unlike `shared.js`/`slider.js`). → Add a friendly guard.
60. `daydream/gui.js:227-229` — boolean URL hydration `val = (val === 'true')` silently treats `'1'`/`'TRUE'` as `false` with no warn path (unlike the numeric/enum branches). → Accept the common truthy spellings or warn.
61. `effects/ShapeShifter.h:48-51` — **(carry-over, prior #89/#9)** `Count` up to 128 dispatches ~256 SDF rasterizations/frame (a Plot *and* a Scan pass per layer) with no per-frame budget cap (author-flagged). → Gate `Count` by resolution or clamp the effective loop against a pixel budget before it ships to a device.
62. **Cross-cutting (effects) — comment volume:** several files carry multi-paragraph rationale blocks larger than the code they annotate (`RingSpin:144-159`, `Dynamo:415-466`, `Thrusters:74-80`), hurting scannability. → Tighten the longest essays.

### Rejected (preserved from prior reviews — do not re-raise)

63. ❌ `core/memory.h:237` (`TriangularBitset::index`) — re-discovered this pass as a device-stripped `assert` on a memory-writing primitive. **Standing rejection (prior #62):** `index()` backs the per-pair dedup primitives called in O(pairs) inner loops, squarely inside the hot-path-doctrine carve-out; the cheapest correct guard is a debug `assert` (zero device cost) with the ordering contract fixed at the source. Distinct from the cold-site always-on `HS_CHECK` policy. No change.
64. ❌ `effects/Dynamo.h` — re-discovered framing of `color()` cost. **Standing rejection (prior #26):** `angle_between` is computed once per point (not per boundary); the boundary loop is ≤15 cheap float compares, 0 in steady state, and the band cannot be hoisted (it depends on the per-point direction `v` and age `t`). No free win.
65. ❌ `tests/test_solids.h` — Islamic/Catalan solids not asserted off-sphere. **Standing rejection (prior #44):** all 23 Islamic registry solids build vertices *on* the unit sphere (measured max |‖v‖−1| = 1.2e-7); an off-sphere assertion would be empirically false. No change.
66. ❌ `core/color.h` — VIBRANT `GenerativePalette` exact-cusp authoring and `ProceduralPalette`-cosines-in-OKLab. **Standing rejection (prior #128/#28):** exact-cusp keys banded the gradient for near-zero vividness gain; the 21 coefficient sets are sRGB-authored, so re-interpreting them in OKLab drives every channel out of range. No change.

---

## 5. Reconciliation With the Prior Review

The fourteenth audit (graded A−, ~18 reviewers) reported the same overall profile this pass
reaches independently: no Critical/High defect, an A/A+ engine, and a gap-to-A that lives at
the coverage frontier. **The tree confirms its fix list 1–127 is essentially fully resolved.**

**Re-discovered this pass but already resolved in the tree** (catalogued here so they are not
re-listed as open in §4). Each was independently re-flagged by a sub-agent reading cold; each
already carries a verified-resolved disposition in the current source:

- The `fib_spiral` unclamped `acosf` (prior #13) — clamped; not re-flagged this pass.
- `Pixel16::operator*`/`blend_alpha` truncation bias (prior #14, #15) — now round-to-nearest.
- `color.h` `g_hue_seed` global cursor (prior #59) — re-flagged by the color reviewer (its #5); the constructor already threads an opt-out `manual_seed`, documented. Resolved.
- `relax` convergence proxy comment (prior #21) — re-flagged by the mesh reviewer (its #5); already reworded to state it tests the raw pre-normalize force. Resolved.
- `classify_faces_by_topology` hash-as-id collision (prior #68) — re-flagged by the mesh reviewer (its #7); already documented with the untrusted-caller bound. Resolved.
- `packPixel` stripped-`assert` hot-path guard (prior #80) — re-flagged by the hardware reviewer (its #3) as a Medium UB surface; already cross-referenced "do not harmonize onto HS_CHECK." The residual is the documented hot-path doctrine; no further change.
- `CMakeLists.txt` debug `STACK_SIZE` link-order fragility (prior #102) — re-flagged by the build reviewer (its #9); already set once per build-type. Resolved.
- `capture_screenshots.mjs` resolution-degradation exit 0 (prior #105) — re-flagged by the build reviewer (its #7); already emits a loud end-of-run banner warning. Resolved-as-documented (exit 0 is the deliberate degrade-don't-abort choice).
- `GnomonicStars` per-frame `fib_spiral` recompute (prior #88) — re-flagged by the effects reviewer as a strength now (cached in `init()`). Resolved.
- The Phantasm `read_id` heuristic (this pass Finding 9) and the beacon-resync precondition the prior audit closed as its #6 (`Config::valid` asserting `beacon_period_revs < 32`) are adjacent but distinct; the resync bound is resolved, the strap-settle window remains open as Finding 9.

**Carry-over open items** (independently re-confirmed open this pass, preserved with their prior
identity):

- Prior #4/#10 → **this pass Finding 4** — ~22/28 effects have no per-effect output invariant; smoke + determinism cannot catch a wrong-but-non-black regression.
- Prior #1–3 → **this pass Finding 5** — the JS browser-glue modules (`driver.js`, `daydream.js`, `recorder.js` lifecycle, `sidebar.js`) remain untested.
- Prior #9/#89 → **this pass Finding 61** — `ShapeShifter`'s unbudgeted per-frame rasterization cost, acknowledged in-code but still shipped.
- Prior #107 → **this pass Finding 19** — the Phantasm size-gate budgets/fixtures are unverified against a real fitting ELF.

**Standing rejections carried forward** (preserved so they are not re-litigated): prior #26, #44,
#62, and #128 — re-stated as this pass's Findings 64, 65, 63, 66 respectively, with their full
rationale.

**Where this pass adds information beyond the fourteenth audit:** it elevates the two **inert
firmware gates** to High (the warning ratchet wired into nothing, Finding 1; the size gate off
by default, Finding 2) — previously folded into Low notes (#108) and the README's
"compensates" claim — and substantiates that no automated device-compile coverage runs today.
It adds the **identity-only transformer test hole** (Finding 3) as the suite's one real
correctness gap, the **segmented-POV simulator/hardware parity gaps** (Findings 10–11), the
**`SmoothUnion` register seam** and **`Warp::Twist` `R==0` divisor** (Findings 7–8), the
**`geometry.js` pole-mapping divergence** (Finding 12), and the **`AABB`/`spatial.h` README drift
introduced by prior Finding 56** (Finding 6). None contradict the prior audit; they refine the
same diagnosis and surface two drifts that prior *fixes* created.

This pass numbers its own fix list from 1 (the established per-pass convention); readers tracking
history should treat §4 as the current-tree superset and prior audits' `✅` entries as closed.

---

## 6. Overall Appraisal — Magnitude of the Achievement

Holosphere is an **A− project whose distance to an A is almost entirely a gate/coverage
frontier, not a quality deficit.** One C++ engine drives a 600 MHz bare-metal Teensy painting a sphere with a strip of
LEDs at 480 RPM, a four-board distributed sibling phase-locked over a single wire, and a
WebAssembly browser twin — sharing 27 original generative effects rendered in 16-bit linear
light with OKLCH perceptual color, on a device with ~335 KB of heap-free RAM. The rare, hard
thing is the engineering coherence that lets a firmware target, a distributed real-time system,
and a browser app share one codebase without the abstractions leaking or collapsing — held to a
self-imposed standard the project then *verifies* (a death harness that proves its own traps
fire) rather than assumes. For the niche of synchronized multi-board perceptually-correct POV
spheres with a bit-identical browser twin, it appears to be the state of the art.

---

*Review generated by independent multi-agent audit of the full codebase (README + fifteen
component sub-reviews), then reconciled against the prior review (§5). No Critical or
High-severity correctness, memory-safety, or concurrency defect was substantiated; the High
items are coverage/gate-configuration gaps. §4 is the complete set of actionable items found.*
