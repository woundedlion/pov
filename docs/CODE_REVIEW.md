# Holosphere / Daydream — Code Quality Review

**Reviewer:** Independent senior C++/TypeScript audit
**Date:** 2026-06-23
**Scope:** `core/`, `effects/`, `hardware/`, `targets/wasm/`, `tests/`, and the `daydream/` web simulator.
**Out of scope (excluded by request):** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/rotate.h`, vendored `three.js/` and `node_modules/`, and generated tables (`color_luts.h`, `reaction_graph.cpp`) beyond a generation-correctness spot-check.

**Method:** The README (≈2,100 lines) was read in full to establish the architectural contract. The codebase (~65 KLOC of hand-written C++/JS across ~60 files) was partitioned across 14 review sub-agents, each cross-referencing the README. Every substantive finding was then re-checked by an independent validation pass against the cited source; two candidate findings were refuted and dropped, and several severities were corrected. All findings below quote verified `file:line` evidence.

---

## 1. Executive Summary

This is an exceptionally mature, deeply-engineered codebase — among the strongest the reviewer has assessed in the hobbyist-to-professional embedded-graphics space. A single C++ rendering engine is compiled three ways (Teensy 4.x firmware, Emscripten/WASM, desktop host) and is held **bit-identical** across all of them, which is the load-bearing invariant for both simulator fidelity and 4-board hardware coherence. The fail-fast philosophy (`HS_CHECK` traps that survive `NDEBUG`) is implemented correctly *and verified by a death harness that asserts the exact trap signal*. Documentation is extraordinary and almost always load-bearing.

**No Critical or High-severity correctness, memory-safety, or concurrency defects were found** anywhere in the ~65 KLOC reviewed. The defect surface is overwhelmingly Low/Nit: a handful of silent-input-masking inconsistencies that run against the project's own fail-fast doctrine, some documentation drift (most notably an RNG-algorithm name in the determinism spec), and several stale comments. The most valuable follow-ups are two one-line guards (non-unit rotation axis; zero-smoothness CSG divisor) and reconciling the `mt19937`→`Pcg32` references in the README determinism contract.

### Overall Grade: **A**

| Quality Dimension | Grade | One-line rationale |
|---|---|---|
| Correctness & Reliability | A | No correctness bug survived validation across 65 KLOC; singularities are branched and tested. |
| Memory Safety & Resource Mgmt | A | Dual-stamp use-after-free detection; overflow-proof arena arithmetic; one `nullptr+0` UBSan nit. |
| Concurrency & ISR / Real-Time Safety | A | Single-writer-by-construction ISR model; time-derived column position immune to masked-IRQ drift. |
| Robustness & Error Handling (fail-fast) | A | Disciplined cold-trap/hot-assert split, verified by an exact-signal death harness. |
| API Design & Interface Expressiveness | A | Compile-time pipeline/palette composition; borrow-vs-store encoded in the type system. |
| Architectural Elegance & Modularity | A | One engine, three targets, bit-identical; clean domain layering; pure cores split for host testing. |
| Readability & Documentation | A− | Best-in-class rationale comments; docked for occasional staleness and over-commenting. |
| Performance & Efficiency | A | Rigorous hot/cold discipline; LUTs, fast-rejects, baked palettes, branchless ISR. |
| Testing & Verification | A− | Death harness + determinism/parity tests are excellent; per-effect render correctness is smoke-only. |
| Numerical / Domain Correctness | A | OKLCH color, geodesic rasterization, flywheel cycle math all verified; one CSG domain gap. |
| Portability & Build | A | Three-target parity via `platform.h`; faithful FastLED integer-semantics mocks. |
| Security (JS↔WASM boundary) | A | Every untrusted JS entry validates-and-rejects; zero-copy view lifetime pinned by `static_assert`. |
| Maintainability & Consistency | A− | A few divisor-guard / pause-gate / comment inconsistencies break otherwise-uniform patterns. |

---

## 2. Dimension-by-Dimension Rationale

### Correctness & Reliability — A
Every projection and interpolation singularity (poles, antipodes, Möbius/stereographic infinity, 0/0) has an explicit, tested branch — rare in graphics-math code. Conway operator composition polarity (the documented "output lands in `temp`, the opposite arena" rule) is implemented exactly as specified and protected by a high-water budget regression test. Reaction-diffusion ping-pong handles both parities; `memcpy` sizes are byte-correct (the classic element-vs-byte count bug is absent). The only validated correctness gaps are domain/guard issues (items 2, 4) that produce wrong geometry or NaN only on degenerate authoring inputs, not in shipped scenes.

### Memory Safety & Resource Management — A
The arena allocator uses subtractive no-wrap bounds arithmetic and a multiply-overflow pre-check, and a **dual-stamp** use-after-free detector (arena generation + per-vector rebind counter) that catches both reset-out-from-under and grow-out-from-under — the second of which most arena designs miss. `Persist<T>` RAII evacuation has correct destruction ordering. The sole defect is a formal `nullptr+0` in `ArenaVector::begin/end` (item 6), UBSan-only and harmless on real toolchains.

### Concurrency & ISR / Real-Time Safety — A
The firmware's defining decision — deriving the POV column from the free-running CPU cycle counter rather than counting timer interrupts — makes masked-IRQ windows (`FastLED.show()`) structurally unable to drop columns. The sync subsystem is single-writer-by-construction (edge ISR publishes into a mailbox; flywheel ISR is the sole consumer), with the one genuinely-needed release/acquire barrier placed correctly and relaxed atomics everywhere else justified for single-core Cortex-M. The `try_claim` mailbox test+reset fusion closes a real race by construction. The double-buffer's relaxed-atomic TOCTOU analysis is rigorous and correct.

### Robustness & Error Handling — A
`HS_CHECK` is a faithful multi-backend fail-fast: predicted-not-taken branch to `__builtin_trap()`, surviving `NDEBUG`, pulling in no stdio of its own, flushing a breadcrumb first. It is disciplined about placement (cold seams trap; hot paths use stripped `assert` backed by a cold trap; genuinely transient conditions get soft handling). The death harness re-execs each trap case in a child and asserts the **exact** `SIGILL`/`STATUS_ILLEGAL_INSTRUCTION` status, so the safety net is verified, not assumed. The few findings here are spots where invalid input is silently *masked* rather than trapped (items 1, 2), violating the project's own doctrine.

### API Design & Interface Expressiveness — A
Standouts: the variadic `Pipeline<W,H,Filters...>` with compile-time domain conversion and three composing ordering `static_assert`s; the `StaticPalette<Source, Coords<…>, Colors<…>, Wrap>` zero-overhead fold-expression composition; `FunctionRef` vs `StoredFunctionRef` encoding borrow-vs-store in the type system with deleted rvalue overloads that turn dangling temporaries into compile errors; and `Feedback` taking a `Style&` directly rather than transform-type template parameters. The `explicit operator CRGB()` guardrails against silent gamma round-trips are thoughtful.

### Architectural Elegance & Modularity — A
The headline achievement — one engine, three build targets, bit-identical output — is realized cleanly through `platform.h`. Compile-time `<W,H>` resolution parameterization yields zero-overhead specialization. Pure, host-testable cores (`pov_sync.h`, `pov_segment_map.h`, `hd107s_frame.h`, `param_marshal.h`, `segment_layout.js`) are deliberately carved out of the device/DOM shells. The four-stage Generate→Transform→Rasterize→Filter pipeline is a coherent mental model that the effects genuinely follow.

### Readability & Documentation — A−
Comments explain *why*, not *what*, and pre-empt plausible "cleanup" refactors that would reintroduce bugs (the LOAD-BEARING / COMPOSITION POLARITY / SCRATCH ARENA CONTRACT annotations are exemplary). Docked half a grade for (a) documentation drift — the `mt19937` references, a stale "clamp below", a "std::min" that is now `hs::clamp`, a "sample twice" that samples three times — and (b) occasional over-commenting where the rationale outweighs the risk.

### Performance & Efficiency — A
Hot/cold discipline is consistent and visible: split sin/cos trig LUTs with quarter-turn cosine recovery, dot-product fast-rejects before `acosf`/`expf`, baked-palette LUTs rebaked only while a wipe is in flight, a branchless segmented ISR with a single display-buffer fetch, exact `lerp16` shift-rounding (0 ULP) that correctly rejects the unsafe `smlad` path, and chroma-preserving gamut mapping gated behind an in-gamut test so the common pixel pays only one matrix multiply.

### Testing & Verification — A−
The death harness and the determinism/parity suite are genuinely rigorous: cross-run byte-identity at both resolutions, full-vs-banded render identity in *both* directions (justifying the segment-stateful full-frame gate), multi-board lockstep through clock wraps, and standalone TUs that recompile the *real* engine under release math flags / device `H_OFFSET` / a 12 KB stack budget. The gap: the roster-wide per-effect test asserts only stability + non-black, not which pixels are lit (a documented, defensible trade-off), and a few subsystems lack focused unit tests (items 37–39).

### Numerical / Domain Correctness — A
Color math (sRGB↔linear LUT generation, OKLab/OKLCH matrices, packed ARM saturating add) was re-derived and matches. The flywheel position math's 64-bit intermediate + half-rev rebase provably makes the 32-bit cycle-counter wrap unobservable. The sync alphabet's odd-only/distance-2 property ("fail to *missed*, never to *wrong*") is correct. One domain gap: `clamp_phi` is only valid on `[-π, 2π]` despite a docstring claiming full real-line folding (item 4).

### Portability & Build — A
`platform.h` reproduces FastLED's *integer* LUT semantics (not smooth-float approximations) and even matches Cortex-M7 divide-by-zero-returns-0 behavior, so the host can't pass on inputs that would misbehave on device. The swap from `std::mt19937` to a concrete `Pcg32` is the correct portability fix (std distributions aren't cross-stdlib-portable). A `__FINITE_MATH_ONLY__` compile-time `#error` guards the engine-wide NaN→hi clamp contract.

### Security (JS↔WASM boundary) — A
Every untrusted JS entry point (`setClip`, `setEffect`, `setParameter`, `fromSolidName`, `bakeLut`, `relax`) validates and rejects-and-logs before acting, reserving traps for true internal invariants. The zero-copy `typed_memory_view` lifetime hazard (views detach on heap growth) is handled comprehensively: all view-backed buffers are pre-sized once, pinned by a `static_assert`, the no-realloc property is unit-tested, and the JS side centralizes detachment detection and asserts a three-alias invariant at every consumer.

### Maintainability & Consistency — A−
The dominant patterns are uniform and well-enforced. Docked for localized inconsistencies that the surrounding code's own standards make conspicuous: one CSG divisor lacks the guard every sibling shape has; one rotation axis lacks the normalize a sibling animation applies; one effect ignores the pause gate its peers wire; one mesh clone omits the view-clearing its sibling writers perform.

---

## 3. Prioritized Findings

Each item is numbered sequentially and lists `file:line`. Severities reflect the independent validation pass.

### High Priority — silent masking of invalid input & spec integrity

1. ✅ **Non-unit rotation axis silently changes the rotation *angle*** — `core/animation.h:1495` (ctor stores `axis` verbatim) and `:1564` (`make_rotation(axis, angle)`). With `|axis| = L`, `make_rotation`'s trailing `.normalized()` makes the effective half-angle satisfy `tan(φ/2) = L·tan(θ/2)`, so a non-unit axis rotates by the *wrong angle* with no trap and no obvious visual tell. The ctor neither normalizes nor `HS_CHECK`s, even though a sibling at `animation.h:711` already calls `make_rotation(axis.normalized(), …)`. Fix: normalize in the ctor or `HS_CHECK` unit length. (`Rotation`, and by extension `RandomWalk`'s public axis.)

2. ✅ **`SDF::SmoothUnion` lacks the divisor guard every sibling shape has** — `core/sdf.h:854` (ctor) and `:955-956` (`h = max(k-|dA-dB|,0)/k`). A `smoothness == 0` yields `0/0 = NaN` that poisons the blend and ships NaN to the LEDs. `Star` (`sdf.h:2453`), `AngularRepeat` (`:1304`), and `Twist` (`:2964`) all `HS_CHECK` their divisors with a NaN-rationale comment; `SmoothUnion` is the lone omission. Fix: add `HS_CHECK(k > 0.0f)` to the ctor.

3. ✅ **Determinism spec names the wrong RNG algorithm** — `README.md:327` (platform table), `:1254` (§7.10), `:1935` (§10.7 worker note), and the in-header comment `core/platform.h:672`. The engine is now `Pcg32(1337)` (`platform.h:169`, `:495`), but these four authoritative references still say `std::mt19937(1337)`. The actual code parity is intact, but the README *is* the contract a re-implementer trusts; anyone re-deriving the RNG to "match the spec" would build an mt19937 and silently break bit-for-bit parity across all 4 boards and the simulator. Fix: replace all four with `Pcg32`.

4. ✅ **`clamp_phi` is incorrect outside `[-π, 2π]` and its docstring overclaims** — `core/sdf.h:137-143`. The docstring says it is "equivalent to `acosf(cosf(x))`" (fold any real into `[0,π]`), but the body does only one reflection, returning out-of-range values for `x < -π` or `x > 2π`. Reachable only when a `Ring`/`DistortedRing` `radius > 2` feeds `center_phi ± target_angle` into `clamp_phi_band` (the ctor stores `radius` unclamped at `sdf.h:407`); sampled effects stay below 2, so currently latent. Result is a wrong bounds row (clipped/over-drawn geometry), not a crash. Fix: correct the docstring at minimum; ideally add a `radius` domain guard or a true full-range fold.

### Medium Priority — latent correctness, robustness, and user-visible inconsistencies

5. **README/comment claims 288×144 on a single Teensy, which the hardware cannot do** — `hardware/pov_single.h:40`. Holosphere is `POVDisplay<40,480>` → `S/2 = 20` rows, and the ISR enforces `HS_CHECK(effect_->height() == S/2)` at `:134`; a 144-row effect would trap, and 144 physical LEDs/half-arm is the Phantasm rig, not a 40-LED board. Fix: drop the "288×144 with one Teensy" phrase.

6. **`ArenaVector::begin()/end()` compute `nullptr + 0` on an unbound/empty vector** — `core/memory.h:691-694`, `:707-710`. Formal UB ([expr.add]/4), UBSan-flagged, reachable via a range-for over a default-constructed vector. The sibling `append_bulk` (`:553-563`) was already hardened against this exact hazard. Fix: guard with `return data_ ? data_ + size_ : nullptr;`.

7. **`Motion` hard-traps on an empty or origin-crossing path** — `core/animation.h:1369` (`angle_between`) and `:1430` (`path_fn(s).normalized()`); `Path::get_point` returns `Vector(0,0,0)` for an empty path (`:82`). Both callees trap on zero length. Likely *by-design* fail-fast (an empty/origin path is not a legitimate animation state), but the trap is far from the cause. Fix: add an earlier, clearer `HS_CHECK(!path.empty())` in the `Motion` ctor.

8. **`MeshState::clone` does not clear destination borrowed views** — `core/spatial.h:460-472`. Unlike `MeshOps::transform` (`conway.h:278-284`) and `update_hankin` (`hankin.h:310-314`), which deliberately unbind stale views on a reused destination, `clone` omits it. Latent only (owned-first accessors shadow two of the three views; no current callsite reuses a borrowed destination, and `Persist`/`Morph` clone into fresh buffers). Fix: clear the three views up front for defensive symmetry.

9. **Raymarch ignores the "Pause Animation" gate** — `effects/Raymarch.h:60,62,65`. Its phase `Driver`s and `Sprite` omit `&anims_paused_` (default `nullptr` = always runs), unlike HopfFibration (`:74-79`), Liquid2D, and Flyby, and with no justifying comment. "Pause Animation" therefore does not freeze Raymarch's spin/palette-scroll/sprite. (`RandomWalk` has no pause parameter in *any* effect — that part is uniform, not Raymarch-specific.) Fix: wire the gate or document the intent. SphericalHarmonics' ungated spin/morph (`SphericalHarmonics.h:245,349`) is the same class of gap.

10. **Liquid2D differences an unbounded time accumulator** — `effects/Liquid2D.h:62` (driven with `wrap=false`) and `:121` (`dt = t - prev_time`). After multi-day continuous uptime the float ULP swallows the per-frame increment and `dt` quantizes/jitters; the field would freeze (graceful, not a crash). FlowField solved the identical problem by wrapping its accumulator to `TIME_PERIOD = 2048` (`FlowField.h:148`). Fix: derive `dt` from `time_speed`×frame-delta, or wrap like FlowField.

11. **Segmented-ISR oversampling comment over-promises the ⅛-column bound** — `hardware/pov_segmented.h:130-133`. With `kOversample=8` the wake grid is ~54 µs but a column-render wake can take ~96 µs (README:1265), and a PIT ISR cannot re-enter, so 1–2 wakes coalesce on exactly the render/emit wakes. The design is *safe* (position is time-derived, `tick()` is idempotent/skip-tolerant), but the "⅛ column" guarantee does not hold on busy wakes. Fix: correct the comment (and ideally measure real emission jitter on-device, since device code is not built in CI).

12. **`MeshMorph` correspondence is O(dest × source) with no cap** — `core/animation.h:2033-2046`. A brute-force nearest-vertex double loop at construction; high-poly morphs can hitch on Teensy at effect-spawn. A KDTree exists in `spatial.h`. Fix: use the KDTree, or at minimum document the bound.

13. **`World::Trails` mid-buffer-expired items erode effective capacity** — `core/filter.h:842-854`. Only front-contiguous dead items are popped; with heterogeneous TTLs (e.g. after `set_lifetime` shrink) dead items linger behind younger live ones, forcing premature eviction of the oldest live point. Documented as intentional, but the capacity-erosion consequence is uncalled-out and untested. Fix: note it, and add a test.

14. **Daydream module-worker load failure recovery relies on a 20 s watchdog** — `daydream/segment_controller.js:202` (module worker), `:303-311` (`INIT_WATCHDOG_MS = 20000`). A failed module-worker top-level fetch (e.g. a renamed `holosphere_wasm.js`) may fire `onerror` with an empty message or not at all in some browsers, so the only backstop for the most common deploy breakage is a 20 s delay before the fault overlay. Handled, but slow. Fix: add a faster readiness probe.

15. **`beat16(bpm)` truncates `bpm ≥ 256` with no documented ceiling** — `core/platform.h:806-808`. `static_cast<uint16_t>(bpm << 8)` discards the high byte, so `beat16(300)` behaves as `beat16(44)`. This is **parity-faithful to FastLED's own `beat16`** (not a divergence), and >255 BPM is musically nonsensical, so this is doc-only. Fix: add a `@pre bpm <= 255` note, matching the other FastLED mocks.

### Low Priority — polish, stale comments, consistency

16. ✅ **`read_id()` brief says "sample twice" but the body samples three times** — `hardware/pov_segmented.h:389` (brief) vs `:417-422` (body) and `:400-414` (correct detailed comment). Behavior is the more-robust three-sample debounce; only the one-line brief is stale.

17. ✅ **Source comment cites a review-finding number** — `tests/fastmath_clamp_check.cpp:1` ("Finding 369"). Violates the project's own "no finding refs in comments" rule; finding numbers belong in commit subjects + this doc only.

18. ✅ **README §7.0 `v2` register row vs `sdf.h` "always 0"** — `README.md:637` says `v2` = "Face index for `Scan::Mesh` (0 otherwise)" (the rasterizer injects it via a wrapper in `scan.h:805-808`), while `sdf.h:196` says `v2` is "reserved and always 0" (true at the `process_pixel` layer, `scan.h:110`). Both are correct for their layer; add a cross-reference so a reader wiring a new SDF shape isn't confused about where `v2` originates.

19. ✅ **Reaction-diffusion seeders duplicate the neighbor sentinel skip** — `effects/BZReactionDiffusion.h:188-190` and `effects/GSReactionDiffusion.h:186-188` hand-roll the `nb >= 0` skip that `ReactionDiffusionBase::for_each_neighbor` (`ReactionDiffusionBase.h:171`) exists to centralize. Add a base seed helper.

20. **`MobiusWarpEvolving` is perpetual; its `.then()` never fires, undocumented** — `core/animation.h:2117-2120` delegates to the default `AnimationBase()` (duration −1), so it never reaches `done()`. The identical RandomWalk hazard is documented at `:1606-1614`; mirror that note here (and for any code relying on its Transformer slot-recycle callback).

21. **`Noise` with a finite duration sawtooths the time axis** — `core/animation.h:2302-2324` is constructed `repeat = true`; a finite duration rewinds `t` to 0 each cycle and snaps `params.time` backward, the opposite of the smooth-flow intent. Fix: force duration −1 or document finite durations as unsupported.

22. **Stale comment: `hue_rotate` references a "per-channel clamp below" that no longer exists** — `core/color.h:751-755`. The clamp is now inside `float_to_pixel16`; the spatial reference misleads.

23. **Stale comment: `Particle::init` cites `std::min`, code uses `hs::clamp`** — `core/animation.h:466-467` and `:543`. Rename-history drift.

24. **`BZ::to_q8` open-codes `255.0f` instead of `Q8_SCALE`** — `effects/BZReactionDiffusion.h:147` vs the constant at `:117`; `GS::to_q16` correctly uses `Q16_SCALE`. Cosmetic symmetry.

25. **Voronoi dead zero-sites guard** — `effects/Voronoi.h:106-107`. `active_site_count()` clamps to `[1, MAX_SITES]` and the slider min is 1, so `sites_buffer` can never be empty; the guard and its comment describe an unreachable state.

26. **IslamicStars redundant scratch reset diverges from the HankinSolids idiom** — `effects/IslamicStars.h:182-185` does a bare `reset()` (already done by `generate()`) before `classify_faces_by_topology`, whereas HankinSolids uses `ScratchScope` for the same operation (`HankinSolids.h:102-113`). Pick one convention.

27. **MeshFeedback `apply_params()` runs before the noise type is configured** — `effects/MeshFeedback.h:71` vs `:73-74`. The first `sync_noise()` propagates a not-yet-configured noise type; benign (re-applied per frame) but backwards from intent.

28. **WASM empty-state getters return `val::array()` instead of a typed view** — `targets/wasm/wasm.cpp:704-716` (`getParamValues`) and the definitions getter. The return type discontinuity (JS Array vs Float32Array) is a latent footgun for a consumer calling typed-array methods on the empty result. Consider returning an empty `typed_memory_view`.

29. **`Persist` post-restore watermark comment oversells the guarantee** — `core/memory.h:1018-1036`. The single-inequality watermark catches the common single-`Persist` forgot-to-reset but not all stacked cases; soften the comment to match.

30. **Vestigial `<functional>` include** — `core/concepts.h:8`. `Fn` is now `inplace_function`, not `std::function`; verify and drop if unused.

31. **GnomonicStars is the lone `show_bg() == true` effect without a comment** — `effects/GnomonicStars.h:84`. Legitimate (sparse star field), but every other effect's `false` is commented; note the divergence.

32. **Thrusters draw-color premultiply asymmetry is undocumented** — `effects/Thrusters.h:284` (RGB×opacity *and* alpha) vs `:304` (alpha only). Presumably an intentional additive-glow double-darkening; add a one-line comment.

33. **Moire density default `10.0` is dead** — `effects/Moire.h:62,168` overwrites it in `init()` before registration; a future editor of the `= 10.0f` member default will see no effect at 288×144 (honestly noted in-comment already).

34. **DistortedRing `amplitude_mut` is a dead template-copy member** — `effects/DistortedRing.h:32-38,69`. Constructed only to be copied into the timeline; correct (the copy retains the `reference_wrapper` binding) but an unusual, slightly wasteful idiom.

35. **`catmull_rom_tangents` "tension" parameter is inverted vs convention** — `core/geometry.h:1373-1380`. `tension=0` → geodesic, `tension=1` → full Catmull-Rom, opposite the standard cardinal-spline convention; documented but easy to misuse. Consider renaming to `smoothing`.

36. **`least_parallel_axis` helper is not adopted by its intended call sites** — `core/geometry.h:855-869` (`make_basis`) and `core/3dmath.h:1093` (`make_rotation` antiparallel branch) still inline the reference-axis fallback the helper was created to centralize. Also: `easing.h:26-30` `ease_in_out_cubic` lacks the endpoint guard its `circ`/`expo`/`elastic` siblings carry — confirm the asymmetry is deliberate.

### Test Coverage Improvements

37. **No focused test for `concepts.h`** — `FunctionRef`/`StoredFunctionRef` overload resolution (non-const lvalue → non-const ctor; rvalue → const ctor; `StoredFunctionRef` rejects rvalues) and `inplace_function` copy/move/empty-trap semantics are exercised only incidentally and would silently regress under a refactor. Add `tests/test_concepts.h`.

38. **`BZReactionDiffusion` lacks a white-box test seam** — `GSReactionDiffusion` has `GSWhiteBox` reaching its fixed-point round-trip and step; BZ has identical fixed-point/stability concerns (`to_q8`/`advance_species`/`perturb_state`) with no isolation seam. Add a matching friend.

39. **Conway composition polarity is not directly unit-tested** — the most error-prone, explicitly-load-bearing invariant in `conway.h` (composed ops return output in `temp`) is exercised only transitively via `SolidBuilder` recipes. Add a focused test asserting `gyro`/`meta`/`needle`/`zip`/`bevel` land output in the expected arena.

40. **Per-effect rendering correctness is smoke-only** — `tests/test_effects.h` asserts `get_pixel` stability + non-black + cross-run determinism, but not *which* pixels are lit, so a regression that changes the image while staying non-black and deterministic would pass. This is a documented, defensible trade-off (golden frames at two resolutions are costly); noted as the one real coverage soft-spot, partially mitigated by per-effect white-box numeric invariants.

---

## 4. Notable Strengths

- **Bit-identical tri-target engine.** One C++ codebase drives Teensy firmware, WASM, and desktop with verified bit-for-bit parity — the foundation for both simulator fidelity and 4-board hardware coherence. The `Pcg32` swap, integer-faithful FastLED mocks, and the `__FINITE_MATH_ONLY__` guard are what make this real rather than aspirational.
- **Time-derived POV column position.** Deriving the column from the CPU cycle counter (not interrupt counts) makes masked-IRQ windows structurally unable to drop columns — the single best idea in the firmware.
- **Verified fail-fast.** `HS_CHECK` is a correct, NDEBUG-surviving, stdio-free trap, and the death harness asserts the *exact* trap signal in a child process, so the marketed safety net is tested rather than assumed.
- **Single-writer ISR concurrency by construction**, with the one genuinely-needed memory barrier placed correctly and every relaxed atomic justified.
- **Dual-stamp use-after-free detection** catching both reset- and grow-out-from-under, the latter of which most arena designs miss.
- **The sync protocol's odd-only/distance-2 alphabet** with a position-weighted beacon checksum — "fail to *missed*, never to *wrong*."
- **Comments that encode invariants** (LOAD-BEARING / COMPOSITION POLARITY / SCRATCH ARENA CONTRACT) and pre-empt the exact "cleanup" refactors that would reintroduce bugs.
- **Comprehensive WASM-view-detachment handling** on the JS side, with a three-alias invariant asserted at every consumer — the single most common Emscripten footgun, closed loudly.

---

## 5. Verdict

A professional-grade, research-adjacent codebase that substantially exceeds the quality norm for its domain. The fact that an exhaustive, independently-validated review of ~65 KLOC surfaced **no Critical or High defect** — only a short list of one-line guards, doc reconciliations, and stale comments — is itself the strongest possible quality signal. Addressing the four High-priority items (items 1–4) would bring the fail-fast doctrine and the determinism spec to full internal consistency; the remainder is polish on an already-excellent foundation.
