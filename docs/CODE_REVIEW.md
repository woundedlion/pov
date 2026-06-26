# Holosphere — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`, `targets/`, `tests/`, build/CI tooling) and the sibling **daydream** web simulator (`c:/work/daydream`). Out of scope by instruction: `core/effects_legacy.h`, `core/rotate.h`, `targets/Holosphere/Holosphere.ino`.

**Method:** the codebase was partitioned across 19 independent reviewer sub-agents (one or more per component), each briefed against the README architecture. Every material finding was then re-examined by a separate validation sub-agent that read the cited code cold and returned CONFIRMED / PARTIALLY-CONFIRMED / REJECTED. Findings that were rejected or shown to describe intended behavior were dropped; several were down-graded after validation (noted inline). Line numbers reflect the tree at the time of review.

**Headline:** this is an exceptionally well-engineered project. Across ~47 KLOC of C++ and ~4.4 KLOC of JavaScript, the review surfaced **no confirmed Critical or High-severity correctness defect in shipping code** — the single High item is a CI-coverage gap, and the most serious code issue (a `DMAMEM` placement bug) is latent dead code behind a disabled build flag. The overwhelming majority of findings are Low-severity documentation drift, micro-optimizations, and consistency nits. The correctness, memory-safety, and fail-fast discipline are at or above professional-product standard.

---

## Overall grade: **A−**

| Dimension | Grade | Rationale |
|---|---|---|
| Correctness | **A−** | Numerical edge cases (poles, wrap, antipodes, NaN/Inf, fixed-point bias) are handled deliberately and consistently. The few real bugs are narrow (an off-by-one that hides one harmonic mode; a capacity-1 transformer that defeats a documented per-sprite intent). |
| Memory safety | **A** | Arena allocator has two-layer use-after-free detection (arena generation + per-vector rebind stamp), subtractive wrap-proof bounds math, RAII `ScratchScope`/`Persist` ordering, and deleted rvalue overloads that turn dangling references into compile errors. |
| Concurrency / ISR correctness | **A−** | Single-writer double-buffer is rigorously reasoned; memory ordering (release/acquire on handoff, relaxed on monotonic counters) is correct and justified; the sync flywheel codec closes subtle wrap/race bugs. Docked only for the latent `DMAMEM` placement issue on the single-board DMA path. |
| Numerical robustness | **A** | `hs::clamp` before every `acos`/cast (with a documented NaN→hi contract and a `#error` guard against `-ffast-math` folding it away), double-precision lattice/LUT generation matched bit-for-bit to the device, Lipschitz bounds for sphere-tracing. |
| Performance / efficiency | **A−** | Compile-time resolution specialization, LUT-backed trig/color hot paths, inlined forwarding-ref filter pipeline (no `std::function` indirection), baked palettes, coarse-grid coherence. A handful of small redundant per-frame computations. |
| API design / interface expressiveness | **A** | The variadic `Pipeline<W,H,Filters...>` with compile-time domain conversion and ordering `static_assert`s, the `Fragment`/shader model, fluent `SolidBuilder`, and the `Filter::{World,Screen,Pixel}` taxonomy are genuinely elegant and uniform. |
| Architectural elegance | **A** | Clean shape-vs-rasterizer split, X-macro single-sources-of-truth (resolutions, effect roster, mesh ops), three targets over one engine, CRTP factoring of shared reaction-diffusion scaffolding. |
| Maintainability | **A−** | Anti-drift guards everywhere (registry-vs-roster counts, golden roster order pin, provenance gates). Docked for a few very large single headers (`sdf.h` 3.2k, `plot.h` 2.3k, `animation.h` 2.8k) and a few hand-copied code blocks flagged in-source as "keep in sync". |
| Readability | **A−** | Comments are dense but load-bearing and accurate — they explain *why* (stability bounds, pole branches, ordering invariants), not *what*. |
| Documentation (in-code + README) | **A−** | Doxygen and the 2,100-line README are unusually thorough and largely match the code; several small README/comment drifts are catalogued below. |
| Test coverage & quality | **A−** | 34 in-process modules + standalone CTest binaries, oracle-based (not smoke), with a hardened death-test harness, determinism passes that scramble globals between runs, and genuine sim-vs-device parity TUs. Light spots: concurrency stress, fuzz/property testing, wave-generator oracles. |
| Build / CI / reproducibility | **A−** | Hard toolchain pinning across all three targets, matched float flags, provenance gates, size/layout gates, self-guarding meta-jobs. Docked for one test module omitted from the shard matrix (High finding #1). |
| Portability | **A** | Pure-math headers cleanly separated from Arduino/Emscripten shells and host-unit-testable; Windows-host dev flow exercised in CI; explicit 32/64-bit reasoning. |
| Error handling / fail-fast | **A** | `HS_CHECK` always-on trap discipline (survives `NDEBUG`, no stdio), cold-path-only placement, transient-vs-invariant distinction, JS-boundary input validation, graceful web degradation. |

### Per-component snapshot

| Component | Grade | One-line |
|---|---|---|
| `core/` math, color, memory, geometry | A | Reference-quality numerics and allocator design. |
| `core/` SDF / scan / plot rasterizers | A− | Rigorous interval/seam algebra; a few latent buffer-bound comments weaker than their guards. |
| `core/` filter pipeline / canvas / animation | A | The showcase abstractions; correct and expressive. |
| `core/` mesh / conway / solids / hankin / spatial | A | Exact pool sizing with compile-time guards; correct topology. |
| `core/` platform / engine infra | A− | Faithful FastLED mocks; minor host/device timing-parity seams. |
| `effects/` (27 effects) | A− | Consistent authoring idioms; mostly micro-perf and two small logic nits. |
| `hardware/` drivers | A− | Excellent ISR/sync reasoning; one latent `DMAMEM` placement bug. |
| `targets/wasm` + Phantasm | A | Solid Embind/zero-copy/param-marshal contracts. |
| `tests/` | A− | Broad, deep, oracle-based; a few coverage gaps. |
| daydream core JS | A− | Correct WASM-heap-detach handling; minor init/teardown nits. |
| daydream segmented / tools JS | A | Correct generation-fence concurrency; two hardening nits. |
| build / CI tooling | A− | Strong gates; one shard-matrix omission. |

---

## Prioritized fix list

Items are numbered sequentially and grouped by priority. Each cites `file:line`, the quality dimension, and the impact.

### High

1. ✅ **`.github/workflows/ci.yml:41-47` & `tests/CMakeLists.txt:96-104` — CI coverage.** The registered CTest `unit_wasm_predicates` is matched by none of the four anchored shard regexes (`^unit_(...)$`), so it runs in zero Linux `tests` shards and zero `sanitizers` (ASan/UBSan) shards — only in the unfiltered `windows-tests` job. The `shard-coverage` meta-job (ci.yml ~405-456) asserts every test maps to exactly one shard, so it is currently RED (`test 'unit_wasm_predicates' is in 0 shards`). Fix: add `wasm_predicates` to one shard alternation (e.g. shard4 alongside `param_marshal`).

### Medium

2. ✅ **`hardware/pov_single.h:228-230` — memory placement (latent).** `ledController_` is defined as a *generic* out-of-line template static carrying `DMAMEM` (`template<int S,int RPM> DMAMEM ... POVDisplay<S,RPM>::ledController_{...}`). GCC silently drops the section attribute on vague-linkage template statics, so the HD107S DMA TX buffer would land in DTCM — which the i.MX RT1062 eDMA cannot read, and where the `arm_dcache_flush()` the file's own comment (224-227) calls "required" is a dead no-op. The segmented path already mitigates this with an explicit specialization in `Phantasm.ino:51-54`; the single-board path does not. **Latent**, not live: the definition is gated behind `#if defined(USE_DMA_LEDS)`, and the only single-board target (`Holosphere.ino`) never defines it (it compiles the FastLED fallback). Fix: mirror Phantasm — move the definition to an explicit specialization at the instantiation site so `DMAMEM` sticks the moment the DMA path is enabled.

### Low

#### Effects logic

3. **`effects/SphericalHarmonics.h:327` — correctness.** `hs::rand_int(1, MAX_MODE_IDX)` is max-*exclusive* (`platform.h:1128-1140`) with `MAX_MODE_IDX == 24`, so the top harmonic (l=4, m=+4) is never selected, contradicting the `[1,24]`-inclusive intent in the nearby comment and `decode_lm`. Fix: `rand_int(1, MAX_MODE_IDX + 1)`.
4. **`effects/DreamBalls.h:152` (+ `:285-288`, `:391`) — correctness / API.** The Möbius warp generator is `MobiusWarpTransformer<1>`, and `draw_scene` applies the single shared `mobius_gen.transform()` to every sprite with no per-sprite slot selection — so the documented intent that "the outgoing sprite warps to its own frozen magnitude" is not achievable at capacity 1. (Validation note: the visible impact is small because the outgoing warp relaxes to identity at the 288-frame boundary just as the next sprite spawns; the comment is more wrong than the pixels.) Fix: size the transformer to 2 for true per-sprite isolation, or correct the comment to state the warp is shared during the crossfade.
5. **`effects/MindSplatter.h:273` — fail-fast consistency.** A particle-index bounds check uses `assert(...)` (stripped under `NDEBUG`, i.e. on device) immediately before indexing `particle_system.pool[p_idx]`, whereas the project standard for always-on invariant traps is `HS_CHECK` (used elsewhere in this same effect family). Fix: replace with `HS_CHECK`.
6. **`effects/PetalFlow.h:178-181` — correctness (cosmetic).** Freshly spawned rings are positioned in `check_spawn()` and then advanced again by `move` in the same `draw_frame()`, so a ring's true first-frame position is `START_RHO + gap + move` rather than the intended `START_RHO + gap`, slightly desyncing spawn density from render motion. Fix: skip advancing rings spawned this frame, or spawn at `START_RHO + gap - move`.
7. **`effects/DistortedRing.h:100` — performance.** `Basis basis = make_basis(orientation.get(), normal);` is recomputed inside the per-ring loop though it is loop-invariant. Fix: hoist above the `for`.
8. **`effects/FlowField.h:64` — performance.** The emitter force loop iterates the full `active_count` prefix including dead-but-draining particles, spending ~3 noise evaluations per dead slot whose velocity `step_particle` then ignores. Fix: `if (p.life == 0) continue;`.
9. **`effects/ChaoticStrings.h:70` — clarity.** `static_assert(SCRATCH_A_BYTES >= MAX_FRAGMENTS * sizeof(Fragment))` sizes only one fragment buffer, but two coexist in `scratch_arena_a` (the effect's `vertices` stays live while `Plot::Multiline::draw` allocates its own). (Validation note: not a defect — the 200 KiB budget is ~10× the one-buffer need and the comment calls out the Multiline headroom; any overflow would fail-fast via the arena trap.) Fix: make the assert reflect both live buffers so it reads as a real bound rather than a loose floor.
10. **`effects/Flyby.h:93-119` — performance.** `Scan::Shader::draw` runs the full per-pixel project→warp→sample→`hue_rotate`→palette chain over every pixel at `<288,144>` with no early-out where pole attenuation drives alpha→0. Fix: skip the color work below an alpha epsilon.
11. **`effects/IslamicStars.h:50,93-96` — robustness (by design).** The ripple-pool `static_assert` is a soft floor: at `Burst=4, Ripp Dur=144` the simultaneous peak can exceed `kRipplePoolSize=8` and `spawn()` silently drops ripples. Documented and intentional; flagged because the assert reads as a hard bound.

#### Core engine

12. **`core/filter.h:1402` — consistency.** `Pixel::ChromaticShift::plot` does `static_cast<int>(x)` (truncation) inside `fast_wrap`, while the base sink (`:162`) and `Screen::Blur` (`:1067`) use `std::round`, so the aberration lands one column off for sub-pixel/tiny-negative `x`. Fix: round to match the family.
13. **`core/animation.h:1747,1811,1863,2153` — robustness.** `MobiusFlow`, `MobiusWarp`, `MobiusWarpCircular`, and `Ripple` divide `t/duration` but omit the `HS_CHECK(duration >= 0)` guard that `Transition` (`:860`), `Mutation` (`:927`), and `Motion` (`:1267`) carry; since `AnimationBase` accepts the `-1` perpetual sentinel, `duration=-1` yields negative progress clamped to a silent freeze. Fix: add the guard (these effects are finite/looped, not perpetual).
14. **`core/animation.h:2409-2434` vs `:2453-2477` — maintainability.** `add()` and `add_get(pin=false)` duplicate the type-erased manager lambda and `static_assert` pair. Fix: implement `add` in terms of `add_get`.
15. **`core/color.h:231` — API footgun.** `Color4()` defaults to opaque black (`alpha = 1.0f`); used as an additive/SSAA accumulator it starts at the wrong identity and `operator+=` (`:287`) saturates alpha at 1.0. No in-scope caller relies on this. Fix: document the dual role, or provide a transparent-zero accumulator constant.
16. **`core/geometry.h:489` — consistency.** `vector_to_pixel` calls the runtime `phi_to_y(phi, H_VIRT)` while the inverse `pixel_to_vector` uses the template `y_to_phi<H>`; both are correct but the asymmetry invites an off-by-`H_OFFSET` bug. Fix: use `phi_to_y<H>(phi)` for symmetry.
17. **`core/geometry.h:521-526` — robustness.** `vectorToLogPolar` can return NaN for a slightly non-unit input (`v.y > 1` → `log(negative)`); the precondition is documented in prose but, unlike its siblings, carries no assert. Fix: add an `@pre`/assert.
18. **`core/sdf.h:1169-1172,1337-1340` — robustness (comment vs guard).** The `kSeamSplitCap` `static_assert` in `Subtract`/`Intersection` bounds the *unsplit* span count, but `normalize_intervals_to_range` can double it, so the assert is weaker than its comment implies (the runtime `push_interval` trap is the real backstop). Fix: assert the post-split worst case, or reword to point at the runtime trap.
19. **`core/scan.h:160-161` — documentation.** The `static_assert` message says the buffer "must hold a single Union's full emission," but the binding case is `Subtract`/`Intersection` (`|A|+|B|+2`). Fix: reword to the max over all top-level CSG ops.
20. **`core/sdf.h:594-595` — rasterization (cosmetic).** Ring horizontal half-width uses `eff_th = 0.95*thickness` while stroke AA in `process_pixel` uses full `size`, so a grazing-row boundary column in the outer 5% of a thin tube can be clipped. Fix: pad the half-width by the column quantum or use full thickness for the interval.
21. **`core/sdf.h:1776-1777` — documentation.** `Face` floors reported `size` to `0.25×circumradius` on sliver faces, so size-normalized shaders see up to a 4× overstated inradius. Fix: note the floor where the README calls `size` "apothem for normalization".
22. **`core/plot.h:1291-1297` — dead code.** `PlanarPolygon::sample` computes a `v1` that omits the close edge, then `rasterize` overrides every `v0`/`v1` from the rendered arc, leaving only a subtly-wrong seed for an optional vertex shader. Fix: drop the loop and rely on the override.
23. **`core/plot.h:647-648,678` — dead code.** Two unreachable guards (`omitLast && steps_cache.is_empty()` early-return and a `total_dist > 0` ternary) on paths where the conditions are provably already satisfied. Fix: remove or convert to invariant asserts.
24. **`core/plot.h:184-283,352-427` — maintainability.** The azimuthal-equidistant unprojection-and-arc-sum loop is hand-copied three times (`rasterize_planar_strategy`, `planar_arc_length`, `edge_row_span`), flagged in-source with "keep in sync"; a future edit to one is easy to miss in the others. Fix: extract one shared helper.
25. **`core/concepts.h:117-124` — type-erasure precision.** The non-const `FunctionRef` ctor constrains on `std::invocable<Callable, ...>` (rvalue) but invokes the callable as a non-const lvalue; a `&&`-ref-qualified `operator()` would satisfy the constraint yet mis-invoke. Theoretical; no current callsite. Fix: constrain on `std::invocable<Callable&, ...>` to mirror the thunk.
26. **`core/memory.cpp:80` (+ `memory.h:34-39`) — coverage gap.** `configure_arenas`'s range guard is checked against the build-dependent `GLOBAL_ARENA_SIZE` (8 MiB under `HS_TEST_BUILD`), so an effect that over-subscribes only the *device* budget (>330 KiB, ≤8 MiB) passes the host suite and traps only on hardware. Fix: gate per-effect footprint static_asserts on `DEVICE_GLOBAL_ARENA_SIZE`.

#### Platform / hardware / targets

27. **`core/platform.h:927-929,967` — sim/device parity.** `EVERY_N_MILLIS` starts `last_ = 0`; on host `millis()` is a huge `steady_clock` value so the first evaluation always fires, but on device `millis()` starts near 0 so it does not fire for the first period — breaking strict parity for setup/spawn logic gated behind it. Fix: seed `last_ = millis()` at construction on both targets, or document the divergence.
28. **`core/platform.h:996-1001,952-957` — portability.** `micros()`/`millis()` cast the chrono count to `unsigned long`, which is 32-bit on a Win32 host (wrap at ~71 min / ~49 days) but 64-bit on LP64 — so wrap behavior is host-bitness-dependent at an arbitrary phase, undocumented in a file that documents every other 32/64-bit seam. Fix: cast through `uint32_t` explicitly if device-parity wrap is intended, and document.
29. **`core/platform.h:675-689` — mock fidelity / doc.** Host `random8(top)` uses `% top` while the adjacent doc block describes the device `(r*lim)>>8`; there is also no `random16(lim)` overload. The header flags these as legacy-only/non-matching, but the doc above `random8(top)` is misleading. Fix: implement the scaled form or correct the doc.
30. **`core/platform.h:928` — macro hygiene.** `EVERY_N_MILLIS` derives its throttle object name from `__LINE__`, so two expansions on one source line collide. Inherited from FastLED. Fix: note the one-per-line constraint or use `__COUNTER__`.
31. **`targets/Phantasm/Phantasm.ino:74-85` — maintainability.** The roster is instantiated at hard-coded `<288,144>` with no `static_assert(TOTAL_PIXELS == CANVAS_W)` tying the constants to the literal, so the two can silently drift on edit. Fix: add the assert or derive the template args from the constants.

#### Documentation accuracy

32. **`README.md:116` vs `:892-896` — memory size.** §2 says the arena block is "335 KB total" while §7.5 and the code say 330 KiB. Fix: change §2 to 330 KB.
33. **`README.md:1074` — mesh API overstated.** "All Conway operators take `(const PolyMesh&, Arena& target, Arena& temp)`" — but `transform`/`normalize`/`relax` (listed in the same table) have different signatures. Fix: scope the claim to the geometry operators.
34. **`README.md:1870` — wasm doc.** `getParamGeneration()` is documented as `→ int` but returns `uint32_t`. Fix: document as unsigned/`number`.
35. **`README.md:1864` — wasm doc.** `setParameter(...)` is documented as returning `false` only "on an unknown name," but the implementation returns `false` for four cases (no effect, unknown name, readonly, non-finite). Fix: broaden to "if the write was not applied."
36. **`reaction_graph.h:62` / `README.md:199` — units.** "90 KB" is correctly 90 *KiB* (92,160 bytes); decimal-KB shorthand. Cosmetic; align units if desired.
37. **`conway.h:248,367` — framing.** `transform` and `normalize` live in `conway.h` and the Conway table but are not Conway operators. Fix: a one-line note or relocation for the next reader.
38. **`.github/workflows/docs.yml:2` — stale URL.** Comment references `woundedlion.github.io/pov/` (old repo name) rather than the Holosphere Pages URL. Cosmetic.

#### Tests & build

39. **`tests/test_canvas.h:253` — flakiness.** `test_ctor_spin_waits_for_buffer_free` uses a fixed 10 ms `sleep_for` to let the main thread enter the ctor spin; under heavy CI load this could race. Fix: gate the release on an observable spin-iteration counter/flag.
40. **`tests/test_death.h` (comments) — staleness.** Comments on `set_case_env`/`spawn_child` still reference `std::system()` though both platforms now spawn shell-free (`_spawnv`/`execv`). Fix: update the comments.
41. **`tests/test_easing_waves.h:44` — depth.** `check_curve` asserts only finiteness/monotonicity; wave generators get no amplitude/period/symmetry oracle, so an amplitude or DC-offset drift could pass. Fix: add bounded-range and period/symmetry oracles for `sin_wave`/`tri_wave`/`square_wave`.
42. **`.github/workflows/ci.yml` (setup-python) — reproducibility.** `python-version: '3.x'` is a floating range, not a pin, despite a comment claiming it pins "like the teensy jobs." Fix: pin a concrete minor (e.g. `'3.12'`) or soften the comment.

#### daydream simulator

43. **`daydream.js:474-481` vs `:547` — redundancy.** Init calls `setResolution` directly and then `applyResolution(true)` calls it again with the same dimensions; the first is redundant and double-logs on the unsupported-resolution path. Fix: drop the pre-set and rely on `applyResolution(true)`.
44. **`segment_controller.js:489-502` — robustness.** `composite()` validates each segment's rect against buffer bounds but never checks `r.pixels.length === (x1-x0)*(y1-y0)*3` before `blitSegmentRect`'s `subarray`/`set`, so a rect/buffer mismatch would silently write a truncated row. Fix: add a length assertion (fail-fast).
45. **`daydream.js:740` — teardown.** `disposeApp()` calls `recorder.stop()`, whose async `onstop` download may not complete before page teardown on `pagehide`, silently losing an in-progress recording. Fix: note the tradeoff or flush synchronously on discard.
46. **`segment_controller.js:534-548` — cosmetic.** `_rebuildBoundaries()` drops the seam of any `null` segment; harmless in normal fenced flow (boundary overlay only). Fix: none required; noted for completeness.
47. **`gui.js:87` — consistency.** Standalone-tool URL writer passes a number to `URLSearchParams.set` (coerced via `String()`), inconsistent with the explicit `String(...)` in `state.js`. Cosmetic.

---

## Coverage gaps (tests)

- **Concurrency / ISR hand-off:** only one threaded test (canvas double-buffer); the relaxed-atomic `queue_frame`/`advance_display` protocol is otherwise validated single-threaded by design. No contention/stress test of the real producer/ISR race.
- **Fuzz / property-based testing:** coverage is example- and sweep-based. No randomized mesh/topology fuzzing for Conway/half-edge, no adversarial SDF/CSG input generation, no random-quaternion round-trip stress (the geometry edge-row-span test does run a 6000-iteration randomized sweep, so the pattern exists but is not applied broadly).
- **Cull tightness:** SDF interval-cull *conservativeness* is checked (cull must cover the interior) but not *tightness*, so a degenerate "always full row" cull would still pass.
- **Thinnest modules:** `test_wasm_predicates.h`, `test_presets.h`, `test_concepts.h` validate compile-time/structural contracts rather than runtime behavior — adequate for scope but the lightest in the suite.

---

## Notable strengths

- **Fail-fast as a first-class, *tested* discipline.** `HS_CHECK` is an always-on, `NDEBUG`-surviving, stdio-free single predicted-not-taken branch to `__builtin_trap()`, placed on cold paths only — and a death-test harness spawns child processes to assert the traps actually fire (`SIGILL`/`STATUS_ILLEGAL_INSTRUCTION`), defeating constant-folding and requiring the *exact* trap shape. The safety net is verified, not assumed.
- **Two-layer use-after-free detection in the arena.** The arena `generation_` catches `reset()`/`rebind()`, and a per-`ArenaVector` `rebind_generation_` catches a re-grow that swaps the data pointer without touching the arena generation — a subtle hazard most arena allocators miss.
- **Sim/device bit-parity is engineered and tested.** `lerp16` MAC-lane layout, the `-ffast-math` NaN-cast contract (with a `#error` if the flag would fold the guard away), south-pole renormalization, and `beat88` mod-2³² intermediates are each pinned by dedicated translation units compiled under the shipping device math flags.
- **Compile-time invariants over runtime hope.** Span budgets (`sdf_max_spans`), arena footprints, mesh pool sizes, palette source counts, coprime cycle periods, and resolution/roster consistency are all `static_assert`-gated, turning whole bug classes into build errors.
- **The filter pipeline.** `Pipeline<W,H,Filters...>` inlines an entire World→Screen→Pixel chain through forwarding-reference callbacks (zero `std::function` indirection), auto-converts between coordinate domains at compile time, and enforces stage ordering with instructive `static_assert`s — a genuinely expressive, zero-overhead abstraction.
- **Honest, scope-bounded tooling.** CI meta-jobs prove the shard union covers every test, the effect roster/resolutions/solids registry are enumerated from the X-macros rather than hand-mirrored, provenance gates regenerate-and-diff the LUT and reaction-graph tables, and a `.sha`/dirty-tree guard prevents a stale binary from riding a green gate to deploy.
- **The web simulator handles the hardest correctness concern correctly:** WASM linear-memory growth detaches typed-array views, and every consumer (`instanceColor`, GUI export, segmented compositor) re-points or guards its aliases before use, with a documented contract that matches the README.

---

*Generated by a multi-agent review with independent per-finding validation. No Critical or High-severity correctness defect was confirmed in shipping code; the catalogued items are predominantly Low-severity polish.*
