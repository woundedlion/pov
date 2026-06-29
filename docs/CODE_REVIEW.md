# Holosphere + daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine and firmware (`core/`, `effects/`,
`hardware/`, `targets/`, `tests/`, `scripts/`, build) and the daydream web
simulator (`*.js`, `tools/`, `tests/`). Out of scope by request:
`core/effects_legacy.h`, `core/rotate.h`, `targets/Holosphere/Holosphere.ino`,
and vendored third-party code (`core/FastNoiseLite*.h`, `core/inplace_function.h`,
`daydream/three.js/**`, `node_modules`, `.pio`).

**Method:** the codebase was partitioned into 20 components and reviewed by
independent agents reading every in-scope file in full against thirteen quality
dimensions. Every raw finding was then re-validated by a separate agent that
re-read the cited code from scratch; a finding entered this report only if that
second pass confirmed it was real, accurately described, and in scope. Of 99 raw
findings, **55 were confirmed** (1 medium, 54 low) and **44 were rejected or
ruled out of scope** — a high rejection rate that itself reflects how few claims
survived adversarial scrutiny.

**Headline:** this is an exceptionally well-engineered codebase. No critical or
high-severity defect survived validation. Every confirmed correctness issue is
either latent (unreachable on the shipping configuration) or cosmetic/visual,
and the remainder are documentation drift, dead code, missing test coverage of
already-guarded invariants, and small interface asymmetries. The engineering
discipline — fail-fast invariants, compile-time budget audits, sim/device parity
contracts, and load-bearing rationale comments — is well above the norm for both
embedded C++ and hobby-scale creative-coding projects.

---

## Letter Grades by Dimension

| Dimension | Grade | Rationale |
|---|:---:|---|
| **Architecture & elegance** | **A** | Clean, consistently-applied layering everywhere: SDF-shape vs. rasterizer split, arena mechanism vs. RAII policy, owned/borrowed mesh handles, recursive variadic filter pipeline with compile-time ordering enforcement, host-testable pure-math cores split from device shells. The compile-time-resolution `<W,H>` parameterization is carried through the whole stack without leaking. Subjectively, the abstractions are expressive and earn their complexity. |
| **Correctness** | **A−** | Math is carefully reasoned and densely edge-guarded (singular cases, antipodes, poles, NaN/Inf, modular wraps all handled and largely tested). The one genuine logic bug — `World::Trails` front-only cull leaking non-monotonic-TTL points (#7) — is real but bounded. Every other confirmed correctness item is latent on the shipping config or a visual/UX divergence. |
| **Memory safety** | **A** | Bump-allocator overflow math is wrap-proof by construction; index narrowing is centralized and trapping; arena lifetime is tracked with dual debug stamps that turn dangles into faults; placement-new uses `std::launder` + proper upcast. The only exposure is the deliberate NDEBUG-stripped `packPixel()` index on the device hot path. |
| **Concurrency / ISR safety** | **A** | The single-core ISR/main-loop model is exhaustively justified; relaxed atomics + `disable_interrupts()` barrier reasoning is explicit. The Phantasm sync engine's single-writer design (fused `try_claim`, release/acquire only on the handoff edge) is a genuinely robust hardware concurrency model, not just a documented convention. |
| **Performance** | **A** | Hot paths are lean and the cold/hot split is disciplined (`HS_CHECK` confined to cold seams, fast vs. exact trig chosen per sensitivity, LUT baking, analytic row/column culling before any distance eval). A few small per-frame rebake/rebuild inefficiencies (#24, #33) are the only blemishes. |
| **Interface design & expressiveness** | **A−** | `FunctionRef`/`StoredFunctionRef`/`PipelineRef` encode lifetime intent in the type system; deleted rvalue overloads turn dangling-reference UAFs into compile errors; fluent builders and `registerParam`/timeline vocabulary read well. Minor friction: leaky `DistanceResult` register overloading, a couple of silent-conflict APIs (#14), and small asymmetries (#6, #53). |
| **Maintainability** | **A−** | Comments explain *why*, not *what*, and capture institutional knowledge future maintainers would otherwise destroy. Detractors are concrete and fixable: dead source (#4, #13, #45), parallel idioms inviting drift (#5), and one real cross-file sync hazard between three hand-maintained test-roster lists (#1). |
| **Documentation** | **A** | Among the best-documented codebases of its kind — rationale-rich, accurate, cross-referenced to the README, and honest about approximations. Doxygen/JSDoc serve as genuine API reference. Deductions are a handful of doc/code drifts (the PiP claim #44, `flush` cull #8, `CRGB(CHSV)` parity #27). |
| **Testing** | **A−** | The suites are a standout — independent analytic oracles rather than re-deriving from code-under-test, death-harness verification of trap paths, determinism tests that perturb global state between runs, device-budget pins at the real shipping config, and anti-tautology parity-plus-golden discipline on the JS side. Held below A by real coverage gaps (some effect math, `driver.js`/`sidebar.js` DOM) and the roster-list sync risk. |
| **Portability (sim/device parity)** | **A−** | Host/device divergences are deliberate, enumerated, and defended (LP64 wrap, H_OFFSET, `-fno-finite-math-only` enforced by `#error`, fixed-point HSV ported byte-for-byte). Small seams remain: float-bit golden pins computed on one toolchain, and the arena-config ordering divergence between targets (#36). |
| **Build system** | **A** | CMake presets, per-build-type stack sizing, math-flag ordering, generator provenance pinning, and X-macro single-source rosters are carefully reasoned. The one real gap is the un-cross-checked CTest/roster list pair (#1). |
| **Security** | **A** | Limited attack surface, but the untrusted JS→embind boundary is consistently hardened (finite-checks, clamps, `textContent` over `innerHTML`, validated URL params, importmap key precedence). No injection or unsafe-format path found. |
| **Code style** | **A** | Consistent naming and member conventions, terse fact-focused comments adhering to the project's own discipline (no history/finding references, Doxygen kept as API reference), no leading-underscore privates in JS. |
| **Overall** | **A−** | A mature, internally-consistent, rigorously-defended codebase whose confirmed defects are overwhelmingly latent or cosmetic. |

---

## Notable Strengths

- **Fail-fast doctrine, applied with judgment.** `HS_CHECK` traps guard cold
  structural seams (allocation, registration, config) and survive `NDEBUG`,
  while the per-pixel hot path deliberately drops to stripped asserts — the
  trade-off is stated and consistently applied, and a death harness verifies the
  traps actually fire.
- **Compile-time budget audits.** Persistent/scratch arena footprints, inline
  animation storage, span-buffer capacities, and registry layouts are pinned by
  `static_assert` against the *real device* constants, so a retune that would
  overflow on hardware fails the build rather than the sphere.
- **Sim/device parity as a first-class concern.** 16-bit linear color, branchless
  wrap/clamp, fixed-point HSV ported byte-for-byte, and a `-fno-finite-math-only`
  `#error` gate keep WASM and Teensy bit-aligned where it matters; divergences
  are documented at the seam.
- **Documentation that captures reasoning.** Singular-case handling, memory-order
  choices, hazard rationale, and "do not 'simplify' this" warnings are recorded
  inline — institutional knowledge embedded in the code.
- **Test oracles, not tautologies.** Tests assert against independently-derived
  ground truth (Euler characteristic + exact V/E/F, brute-force k-NN, analytic
  rotation images), perturb global state between determinism runs, and pair
  wasm-vs-js parity checks with external goldens to catch coordinated drift.
- **Hardware concurrency done right.** The Phantasm flywheel timebase (position
  from `now − epoch`, epoch folded forward per crossing) makes dropped columns
  structurally impossible, and the single-writer sync engine fuses the one
  cross-context handoff so a race cannot be reintroduced.

---

## Prioritized Items to Fix

Items are numbered sequentially. Severity reflects post-validation impact;
nearly all confirmed correctness items are latent on the shipping configuration
and are flagged as such.

### Medium Priority

1. ✅ **[build-system] `tests/CMakeLists.txt:96-104` vs `tests/run_tests.cpp:66-104` — CTest module list and run_tests roster are independently maintained with no cross-check.** CMake's `_hs_test_modules` (one CTest per `run_tests <name>`) duplicates `HS_TEST_MODULE_LIST`; `run_tests.cpp` only static_asserts internal consistency, which cannot catch a CMake/roster divergence. A module added to the roster but omitted from `_hs_test_modules` is *never run under ctest* and fails silently in CI sharding. *Fix:* add a CTest that runs `run_tests --list` and diffs the emitted names against `_hs_test_modules` (or generate the CMake list from the roster); at minimum add a `run_tests` mode that exits non-zero if any module name is absent from a list passed on argv.

### Low Priority

2. ✅ **[documentation] `core/3dmath.h:1198-1203` — `fast_cosf` omits the large-argument range-reduction warning `fast_sinf` carries.** It is `fast_sinf(x + π/2)` and inherits the precision collapse for large `|x|`, but its doc lacks the `@warning` and `STEREO_PATTERN_ARG_LIMIT` pointer. *Fix:* mirror `fast_sinf`'s `@warning` or add `@see fast_sinf`.

3. ✅ **[documentation] `core/3dmath.h:385-393` — `Spherical(const Vector&)` traps on a zero/degenerate vector but does not document the precondition.** It calls `normalize()`, which `HS_CHECK`-traps below `~1e-6`; every other trap contract in the file is documented at the call site. *Fix:* document the non-degenerate-input precondition (or route legitimate zeros through a `normalized_or`-style fallback).

4. ❌ **[maintainability] `core/palettes.h:15-99` — seven named `ProceduralPalette` instances are dead source.** `vintageSunset`, `lateSunset`, `lemonLime`, `algae`, `darkPrimary`, `desertRose`, `bruisedBanana` are referenced nowhere; they inflate the file and the README inventory (compiler elides them, so no flash/RAM cost). *Fix:* delete them (and trim the README), or document them as an intentional reserve bank. *Rejected:* intentional reserve palette bank — unused instances are elided by the compiler (zero flash/RAM) and remain available to effect authors; no deletion.

5. ✅ **[maintainability] `core/color.h:992-1001` vs `2153-2167` — two divergent idioms for the same `[0,1)`-frac → `uint16_t` lerp weight.** `Gradient::get` uses `float_to_pixel16(frac)` (clamps+rounds); `BakedPalette::get`/`Color4::lerp` inline `(frac*65535+0.5)` (no clamp). Both correct today; invites silent divergence. *Fix:* route all three through one helper (`frac_to_q16`).

6. ✅ **[interface-design] `core/color.h:617-627,577-581` — OKLab helpers mix struct returns with `float&` out-parameters.** The forward direction returns `LMS`/`OKLab` structs while `oklab_to_linear_rgb*` threads three out-params, forcing callers to declare loose uninitialized locals. *Fix:* introduce a `LinRGB` POD and return it, keeping an out-param overload only if a measured hot path needs it.

7. ✅ **[correctness] `core/filter.h:756-785` — `World::Trails::flush()` cull reclaims only front items, leaking non-monotonic-TTL points.** The `while (at(0).ttl==0) pop_front()` cull is correct only if TTL is monotonic with insertion order, but `plot()` seeds `ttl = lifetime − round(age)` with per-point varying age (Orient/OrientSlice/replication), so a later point can die in the *middle* of the ring; its slot is never freed and accumulating zombies cause premature eviction of live points. `Screen::Trails::decay()` already handles this correctly with swap-remove. *Fix:* cull all dead slots (compact/swap-remove), matching `Screen::Trails`; or document and enforce a monotonic-age precondition.

8. ✅ **[documentation] `core/filter.h:745-755` — `flush()` doc claims it "culls the dead" but only the oldest contiguous run is culled.** Doc overstates the guarantee given the non-monotonic-TTL case in #7. *Fix:* once #7 is fixed the comment is accurate; if front-only is kept deliberately, say "culls from the oldest end" and document the monotonic-TTL precondition.

9. ✅ **[correctness] `core/filter.h:963-975` — `Screen::Trails` silently drops the *newest* seed point at `MAX_PIXELS`, unlike `World::Trails`.** At capacity the current sample is forwarded but not recorded, so the trail tip freezes while stale points keep their slots — the opposite of `World::Trails`' evict-oldest behavior. *Fix:* pick one overflow policy for both filters (evict-oldest, or document drop-newest; or `HS_CHECK` if overflow is an authoring error).

10. ✅ **[testing] `core/sdf.h:707-743` — `DistortedRing`'s load-bearing `max_distortion` bound is validated only under NDEBUG-off, by 256-point sampling.** An underestimate silently culls genuine arcs with no device diagnostic; the sampling can miss a sharp peak. *Fix:* add a property test with adversarial high-frequency `shift_fn`s asserting no arc columns drop vs a full-scan reference.

11. ✅ **[correctness] `core/sdf.h:86-98` — top-level `Subtract`/`Intersection` conservative span bound (66) can exceed `scan_region`'s 64-interval buffer.** With both children at `kIntervalSpanCap=32` the static bound is `|A|+|B|+2=66 > 64`; safety currently rests on a runtime `push_interval` trap plus a "no SDF shape emits 32 disjoint arcs" argument, not a static guarantee. *Fix:* size the buffer to the actual top-level max, or `static_assert` a per-shape emit-count trait so `|A|+|B|+2 ≤ 64` holds at compile time.

12. ✅ **[correctness] `core/scan.h:314-315` — `BoundingSphere` half-width cap at `W/2` leaves a hemisphere-or-larger cap one column short for odd `W`.** `2*x_half ≤ W−1` when `W` is odd, so a full-coverage row paints `W−1` columns. Latent: Phantasm's `W=288` is even. *Fix:* let the clamp reach full width for odd `W` (e.g. `(W+1)/2` on one side) or rely on the `cos_dtheta≤−1` full-row branch.

13. ✅ **[maintainability] `core/animation.h:2143` — `RippleParams::frequency` is a dead field never read by any transform.** `ripple_transform()` derives the wavelet from amplitude/phase/thickness/decay/center only; the field misleads readers into thinking spatial frequency is tunable. *Fix:* remove it (and its doc), or wire it into `ripple_transform`/`prepare_thresholds`.

14. ✅ **[interface-design] `core/animation.h:1043-1076` — `Driver::set_speed` is silently overwritten every frame when a live `speed_src_` is bound.** `step()` re-reads `*speed_src_ * scale_`, clobbering a manual `set_speed` with no diagnostic; `get_speed()` reports the source value. *Fix:* `HS_CHECK(!speed_src_)` inside `set_speed`, or have it clear `speed_src_` to switch to manual mode; document the chosen semantics.

15. ✅ **[correctness] `core/animation.h:1227-1235` — `Sprite` with `fade_out_duration ≥ duration` never reaches full opacity and starts already fading.** The fade-out condition is true from `t==0`, capping peak opacity below 1.0 — easy to hit with independent GUI sliders. *Fix:* clamp the combined fade window (or scale proportionally) so the sprite hits opacity 1.0, or document the capped behavior.

16. ✅ **[build-system] `core/animation.h:2385-2402` — inline-storage budget audit relies on manually appending each new non-templated animation type to `kLargestConcreteAnimSize`.** The `static_assert` only protects listed types; a new non-templated type can escape the audit until exercised. *Fix:* drive the audit from a single type list shared with a test that round-trips every animation type through `Timeline::add`.

17. ✅ **[interface-design] `core/generators.h:46-67` — `generate(target, …)` silently resets `target` if a scratch arena is passed as `target`.** At depth 0 it unconditionally resets both scratch arenas, then `ScratchScope` rewinds them on exit — destroying the output. Latent (all callers pass `persistent_arena`) but produces silent garbage, contrary to the fail-fast philosophy. *Fix:* add `HS_CHECK(&target != &scratch_arena_a && &target != &scratch_arena_b, …)` and tighten the doc to a hard precondition.

18. **[testing] `core/memory.h:961-968` — `~Persist()` forgot-to-reset watermark trap has no death-test coverage.** It is exactly the cold-path trap class the death harness pins elsewhere. *Fix:* add a death case that omits `persistent_arena.reset()` over a trivial `Cloneable` and asserts the trap fires.

19. **[testing] `core/memory.h:222-233` — `TriangularBitset::index()` ordered-pair/range `HS_CHECK` is untested.** Guards against swapped/out-of-range pairs that "write adjacent memory"; memory-safety-relevant and unpinned. *Fix:* add a death case with a swapped pair or `large ≥ MAX_V` via `opaque()`.

20. **[correctness] `core/conway.h:528-559,374-378` — `ambo` edge-midpoint + strict `normalize()` traps if an edge's endpoints are antipodal.** `(v1+v2)*0.5` near the origin makes `normalize()` trap; `dual`/`kis` use `normalized_or` with a fallback for exactly this case. Latent on the closed-manifold roster. *Fix:* use `normalized_or(mid, mesh.vertices[v1])` to match `dual`/`kis`.

21. **[correctness] `core/spatial.h:81-82` — `KDTree` constructor bound (`≤ UINT16_MAX+1`) is misleading; the real ceiling is the `int16_t` child links (`INT16_MAX+1`).** An oversized point set trips the per-node child-link trap mid-build rather than failing clearly at construction. *Fix:* tighten the constructor bound to `count ≤ INT16_MAX + 1`.

22. **[correctness] `core/mesh.h:306` — `vertices[v].half_edge` records only the last incoming half-edge, an unreliable orbit start for boundary meshes.** Only `relax()` documents and works around this; a future operator trusting the entry point on a boundary mesh would silently lose that vertex's ring. *Fix:* document the boundary-edge caveat on `HEVertex.half_edge`/`build_half_edge_mesh` (or prefer recording a paired half-edge when one exists).

23. **[correctness] `core/hankin.h:323-330` — `update_hankin` degenerate-intersection fallback uses strict `normalized()` that can itself trap.** `(m1+m2).normalized()` traps when the shared-corner midpoints are near-opposite; the function's other fallbacks correctly use a known-unit corner. Low probability (self-intersecting `truncate>0.5` recipes), hard abort. *Fix:* `intersect = normalized_or(m1+m2, p_corner.normalized())`.

24. **[performance] `core/hankin.h:333-360` — per-frame `update_hankin` rebuilds the output mesh with element loops instead of `append_bulk`.** Called every frame during angle sweeps; all four element types are trivially copyable and `clone()`/`finalize_solid` already use `append_bulk` (single memcpy vs N push_backs). *Fix:* replace the static/dynamic vertex and faces loops with `append_bulk`, keeping the `face_counts` loop where it interleaves the prefix sum.

25. **[documentation] `core/reaction_graph.h:84-115` — `CubemapLUT::build()` doc omits the ~90 KB transient lattice scratch it requires.** The "48 KB" figure is only the persistent table; `build()` also allocates `RD_N` `Vector`s (~90 KB) from the same arena. A caller sizing from 48 KB under-provisions and traps. *Fix:* document the transient `RD_N*sizeof(Vector)` scratch.

26. **[performance] `core/hankin.h:325-330` — redundant re-normalization of an already-unit fallback intersection.** The fallback branch normalizes, then the code unconditionally normalizes again. Trivial cost; folding the branch also fixes #23. *Fix:* give the degenerate branch its own assignment+continue using `normalized_or`, leaving the trailing `normalized()` for the primary cross-product path.

27. **[documentation] `core/platform.h:298-380` — mock `CRGB(const CHSV&)` does NOT "behave identically to the device".** The host uses the 6-sector integer spectrum; FastLED's device path uses `hsv2rgb_rainbow`, giving visibly different RGB. Impact is currently nil (no modern effect uses `CHSV`), but the blanket parity claim is wrong and could mislead. *Fix:* scope the claim — note `CRGB(CHSV)` is a legacy-compat spectrum conversion that must not be used on parity-sensitive paths (or port `hsv2rgb_rainbow` if exact parity is ever needed).

28. **[performance] `core/platform.h:1121-1139` — runtime `HS_CHECK(max == UINT32_MAX)` in `random_to_unit` is effectively dead.** Its only callers pass a value already pinned by a `static_assert`; the per-call branch can never fire. *Fix:* drop the runtime check (the `static_assert` in `rand_f` is the real guard) or make the divisor a `constexpr`/template parameter.

29. **[testing] `core/effects.h:40-83` — a fully-forgotten effect silently drops from native smoke coverage.** The WASM registry-vs-`HS_EFFECT_COUNT` check only catches a mismatch; if an author forgets the `#include`, `REGISTER_EFFECT`, *and* `X()` row, counts still agree and nothing fails. *Fix:* add a CI cross-check that scans `effects/*.h` for `REGISTER_EFFECT(` and asserts the set equals `HS_EFFECT_LIST`.

30. **[interface-design] `core/canvas.h:179-182` — `set_margin` guards `m < w` but not a negative margin.** A negative margin propagates into `ClipRegion` and inverts/empties the render band (silently skipping rows). Latent (no negative caller today), but the lower bound is the easier mistake. *Fix:* `HS_CHECK(m >= 0 && m < clip.w, …)`.

31. **[maintainability] `effects/Flyby.h:195-229` — `Flyby::Params` lacks the field-drift `static_assert` that its twin `Liquid2D::Params` has.** `lerp()` hand-lists 6 fields and presets are 6 bare floats; a field add/reorder silently mis-maps both. *Fix:* add `static_assert(sizeof(Params) == N*sizeof(float), …)` (or convert `lerp` to the indexed-array form) so a field change is a compile error.

32. **[documentation] `effects/GSReactionDiffusion.h:22-31` — GS class doc omits the memory-budget table BZ carries, despite GS being the heavier arena tenant.** GS bakes an extra palette and runs 16 substeps/frame; its budget exists only inside an `init()` `static_assert`. *Fix:* add a parallel budget table matching BZ's style.

33. **[performance] `effects/Moire.h:65-98` — palette LUTs rebake nearly every frame because the wipe period equals the wipe length.** With `period == WIPE_FRAMES == 80`, a new wipe starts the frame the previous ends, so the "skip rebake when idle" guard never engages after frame 80. *Fix:* make the timer period strictly greater than `WIPE_FRAMES`, or drop the now-misleading "skip otherwise" comment.

34. **[correctness] `effects/DreamBalls.h:285-301` — single-slot Möbius warp is shared across overlapping fade-out/fade-in sprites, violating its own documented invariant.** The spawn period (288) equals the warp duration, so an incoming sprite re-arms the warp at the exact boundary while the outgoing sprite is still fading (~32 frames), warping it with the *incoming* (non-identity) warp — opposite the "relaxed to identity" assumption. Visually subtle. *Fix:* give the transformer two slots, or correct the comment to state the shared-warp overlap and argue it is acceptable.

35. ✅ **[correctness] `hardware/pov_single.h:131-136` — `cols_per_min` division is unguarded against zero (div-by-zero before the `≥1` check).** The `HS_CHECK` on `interval_us ≥ 1` fires *after* the divide; if `RPM` or `width` were 0 it is an integer divide-by-zero. Not reachable with `Holosphere<40,480>`, and the segmented driver has the symmetric guards. *Fix:* add `static_assert(RPM > 0)` and `HS_CHECK(cols_per_min > 0)` before the division.

36. ✅ **[portability] `targets/Phantasm/Phantasm.ino:68-76` vs `targets/wasm/wasm.cpp:362-379` — `configure_arenas_default()` runs after the constructor on Teensy but before it in WASM.** Benign today only because effect constructors never allocate from the engine arenas; a future constructor allocation would work in WASM and corrupt on Teensy with no compile-time signal. *Fix:* pick one ordering for both targets (preferably configure-before-construct, matching WASM), or document the "constructors must not allocate" convention at both seams.

37. ✅ **[maintainability] `scripts/capture_screenshots.mjs:65-71` — the launch-failure path uses `process.exit(1)`, risking truncation of the actionable warning it just printed.** The sibling `wasm_smoke.mjs` documents the opposite (`process.exitCode = 1; return`) so buffered output flushes. *Fix:* use `process.exitCode = 1` and fall through.

38. ✅ **[correctness] `scripts/capture_screenshots.mjs:22-25` — `numEnv` treats an empty-string env var as a finite 0, silently disabling the timing it guards.** `Number('') === 0` passes the finite check, collapsing settle/retry waits to zero — the exact failure the comment claims to defend against. *Fix:* treat empty/blank as unset and reject negatives.

39. ✅ **[correctness] `targets/wasm/wasm.cpp:860-874` — `getFaces()` bounds-checks each index but never asserts all face indices were consumed.** Over-reads trap (`flat_idx < size`) but under-reads (`sum(face_counts) < faces.size()`) silently drop the tail, reading back as valid-but-truncated geometry in the JS editor. *Fix:* add `HS_CHECK(flat_idx == mesh.faces.size(), …)` after the loop.

40. ✅ **[interface-design] `targets/wasm/wasm.cpp:660-676` — stack metrics are emitted as `unsigned` while arena metrics are `size_t`, giving one JS object mixed numeric widths.** The `static_cast<unsigned>` would truncate on a wasm64/4GB-stack build and is the only region produced by a different integral type for no stated reason. *Fix:* emit the stack fields as `size_t`/`uintptr_t`.

41. ✅ **[testing] `tests/test_animation.h:453-503` — easing endpoint/finiteness coverage is duplicated from (and weaker than) `test_easing_waves.h`.** Two places to update per easing change, with assertions buried in a module whose scope is the animation system. *Fix:* drop the duplicated cases and rely on the easing_waves module; keep only animation-specific easing usage.

42. ✅ **[testing] `tests/test_color.h:1026-1053` — `BakedPalette` test arenas size `LUT_SIZE*sizeof(Color4)+64`, encoding a single-allocation assumption with a magic `+64`.** A future second allocation in `bake()` overflows and traps inside the test instead of failing cleanly. *Fix:* expose `BakedPalette::required_arena_bytes()` and size buffers from it.

43. **[maintainability] `tests/test_plot_scan.h:44-56` (and other `test_*.h` statics) — module-static arenas/buffers are reused across invocations with no reset contract, coupling test order to correctness.** Correctness relies on every test wrapping allocations in a `ScratchScope`; a future test that forgets it fails intermittently by run order. *Fix:* document the per-test reset contract at each shared-arena accessor, or have a fixture re-seat the arena offset on each test entry.

44. **[correctness] `daydream/driver.js:589-607` — the PiP camera duplicates the main view instead of showing the opposite hemisphere it documents.** `renderPip()` copies the main camera's position *and* quaternion, so README §10.3's "front and back visible simultaneously" no longer holds; the PiP adds no new information. *Fix:* render the PiP from the antipodal vantage (negate the main position, `lookAt` origin), or update README §10.3 to describe a square-cropped duplicate.

45. ✅ **[maintainability] `daydream/driver.js:871-891` — exported `coordsLabel` is dead code that violates the unit-vector `getLabels` contract (latent double-scale + facing-test bug).** No consumer exists; it returns a `SPHERE_RADIUS`-magnitude position while `LabelPool.acquire` multiplies by `SPHERE_RADIUS` again and the facing test assumes unit positions. A latent trap for the next maintainer who wires it up. *Fix:* remove it, or return a unit direction (`setFromSphericalCoords(1, …)`).

46. ✅ **[correctness] `daydream/recorder.js:202-241` — `elapsedSeconds` diverges from real video duration on the timed-fallback capture path.** It always advances by `frameInterval`, but the timed path (`captureStream(fps)`) encodes on wall-clock, so the displayed M:SS disagrees with the file length. *Fix:* track manual-vs-timed and derive elapsed from `performance.now()` on the timed path (or document it as manual-mode-only).

47. ✅ **[correctness] `daydream/segment_layout.js:87-100` (via `segment_worker.js:160`) — `extractSegment` silently short-copies rows if the source view is shorter than `canvasW*canvasH*3`.** `subarray` clamps rather than throwing, so a stride/length mismatch corrupts the segment with zeros instead of faulting. *Fix:* assert `allPixels.length === canvasW*canvasH*3` (or fault) before extracting.

48. ✅ **[maintainability] `daydream/segment_controller.js:108` — `aliasDivergenceLogged` is never reset, suppressing the composite alias-divergence warning for the whole page lifetime after the first occurrence.** Not cleared in `destroy()`/`create()` unlike every other session flag, so a recurring divergence in a rebuilt pool is hidden. *Fix:* reset it in `destroy()`/`create()` to keep the throttle per-session.

49. ✅ **[correctness] `daydream/recorder.js:133-203` — `captureStream(0)` fallback can store an `undefined` track, silently no-op'ing capture for the whole session.** If both capture attempts yield no video track, `captureFrame` guards on `!this.track` and does nothing, with no error surfaced. *Fix:* if the track is still undefined after the fallback, stop the stream, log a clear error, and return without constructing the `MediaRecorder`.

50. ✅ **[maintainability] `daydream/segment_worker.js:162-187` — arena-metrics conversion hand-copies fields and silently drops the engine's `stack` metric.** The hand-rolled copy must stay in lockstep with the C++ binding and the `SegArenaMetrics` typedef; a new C++ field is dropped with no type error. *Fix:* iterate the metric object's keys generically, or document that `stack` is intentionally omitted and pin the typedef as authority.

51. ✅ **[interface-design] `daydream/gui.js:178-193` — `attachUrlWriter` wraps `controller.onChange` but not `onFinishChange`, so an `onFinishChange` caller bypasses URL persistence and the load-time replay.** Reveals an asymmetry in the deep-link contract. *Fix:* document that only `onChange` participates, or also wrap `onFinishChange`.

52. ✅ **[correctness] `daydream/tools/clipboard.js:84-100` — `copyWithFeedback` leaks the prior outcome's CSS classes on a rapid success-after-failure flip.** The revert closure removes only the current call's `flashClasses`, so a failure-then-success leaves `failedClasses` on the element. Invisible with the default `COPY_FEEDBACK`, but the documented API supports distinct copied/failed classes. *Fix:* on re-entry remove both `copiedClasses` and `failedClasses` before applying the new flash.

53. ✅ **[interface-design] `daydream/tools/solid_codegen.js:26-139` — `relax` requires a params object despite fully defaulting, asymmetric with `snub`.** A bare `'relax'` string is rejected though `iter` defaults to 8, while the equivalently-defaulting `snub` accepts its bare string. *Fix:* remove `relax` from `PARAMETERIZED_OPS` (guarding the param read like `snub`), or document the intentional asymmetry.

54. ✅ **[testing] `daydream/tests/color.test.js:90-99` — the non-wrapping hue midpoint is computed but never asserted.** The test pins `mid.L` and `mid.C` but not `mid.h`, so a regression in ordinary (non-seam) hue lerp would pass. *Fix:* add `near(mid.h, 0.6)`.

55. ✅ **[maintainability] `daydream/tests/solid_codegen.test.js:9-18` (+ `cpp_format.test.js`, `spline_math.test.js`) — the same shared `formatFloatCpp` is pinned redundantly under three names.** Three suites re-pin the same formatter's whole-number/`.0f`/trailing-zero behavior; the downstream copies are weaker and risk drifting from the authoritative `cpp_format.test.js`. *Fix:* keep authoritative coverage in `cpp_format.test.js` and reduce the others to a one-line identity/wiring smoke check.

---

*Generated by a 20-component parallel review with independent per-finding
validation (119 agents; 99 raw findings → 55 confirmed, 44 rejected/out-of-scope).*
