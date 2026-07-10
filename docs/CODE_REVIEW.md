# Holosphere + daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware (`core/`, `effects/`,
`hardware/`, `targets/`, `tests/`, `tools/`, build system) and the daydream web
simulator (`*.js`, `tools/`, `tests/`). Out of scope per instruction:
`effects_legacy.h`, `Holosphere.ino`, `rotate.h`.

**Method:** the codebase was partitioned into 13 component scopes, each reviewed
in depth by a dedicated agent. Every candidate finding was adversarially
validated against the cited code before inclusion — a finding survives only if
its failure mode is reproducible from the source. The great majority of drafted
candidates were *disproved* during validation (the surrounding code already
guarded them) and dropped. What remains below is the residue that survived that
filter: there are **no Critical or High-severity defects**, and no reproducible
crash / out-of-bounds / use-after-free / data-race anywhere in scope. The
findings are latent guards, silent-divergence hazards, CI coverage gaps, and
documentation nits.

---

## Overall Grade: **A**

This is, by defect density, one of the cleanest codebases of its size and
ambition one is likely to review: ~55k lines of hard real-time, arena-allocated,
compile-time-specialized C++ (plus a WASM bridge and a threaded JS simulator)
in which a deliberately adversarial multi-agent pass could not manufacture a
single qualifying correctness bug. The engineering discipline — fail-fast
invariants verified by a death-harness, single-source-of-truth rosters pinned by
CI, host-testable extraction of every load-bearing algorithm, and a filter
pipeline whose legality is enforced at compile time — is at or above
professional / state-of-the-art standard for embedded creative-coding engines.

The only systemic weaknesses are matters of *consistency* and *comment
volume*, not correctness: a convention (a compile-time budget guard, a
trait-derived config flag, a graceful-normalize fallback) that is applied in most
places but omitted in a few, and rationale comments dense enough to risk drift
from the code they annotate.

### Grades by dimension

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | **A** | No reproducible defect across the entire surface; every numerical / topological / concurrency edge is either provably bounded or fail-fast guarded. Only latent, currently-unreachable gaps remain. |
| **Architecture & Elegance** | **A** | Principled and cohesive: World/Screen/Pixel domain-typed filter pipeline, CRTP animation/pipeline chains that inline to zero overhead, pure host-testable cores behind thin platform shells, X-macro single-sourcing. Hardware driver split and compile-time-legal filter chain are **A+** pockets. |
| **Interface Expressiveness** | **A−** | Intent is largely encoded in the type system (owned-vs-borrowed vectors, `Tweenable`/`SDFShape` concepts, deleted rvalue overloads, `explicit` coordinate-space ctors). Docked one notch for a few conventions that live in prose rather than types (Fragment `v0–v3` / `DistanceResult` register meanings, `Path::get_point`'s non-unit return, composition polarity, `Style` copy-vs-lerp asymmetry). |
| **Error Handling & Robustness** | **A** | Disciplined two-tier policy: always-on `HS_CHECK` traps at cold invariant/capacity seams, stripped `assert` on hot paths backed by cold setup traps, and reject-and-log at untrusted JS/WASM boundaries. The distinction is reasoned at each site. |
| **Performance** | **A** | Measured, not guessed: interpolated LUTs, ARM SIMD color packing, warp-field caching keyed on purity, active-slot compaction, deliberate `noinline` to dodge a measured register-spill cliff. No render-path heap. Color path is **A+**. |
| **Testing & Verification** | **A** | Death-harness asserting traps actually fire, ASan/UBSan, byte-identical determinism replay, size/layout gates recompiled under shipping flags, a mock-transport DMA suite and an event-driven multi-board sync simulator, host-test cross-checks of the JS/C++ boundary. Rosters triple-pinned against silent drift. |
| **Build / CI / Tooling** | **A** | Sharded matrix with single-source anchors + a shard-coverage guard, pinned toolchains, full-SHA-pinned third-party actions, ccache-engagement assertions, gates that gate *themselves* on broken fixtures. A couple of narrow coverage gaps (below). |
| **Maintainability** | **A−** | Excellent structure and drift-prevention pins, but comment density is extreme — several headers are majority prose, occasionally heavier than the code they annotate — which raises edit cost and drift risk. A handful of convention divergences are the kind of thing that erodes over time. |
| **Documentation** | **A** | Reference-quality Doxygen and a genuinely exceptional (if enormous) README that matches the code, including honest statements of numerical limits. Docked only for a few doc claims that overreach their code (below). |
| **Security / Safety** | **A** | Limited attack surface (no network, no untrusted deserialization); the one boundary that matters — the WASM/JS ABI — validates and clamps every untyped input, and DOM overlays use XSS-safe text nodes. |

---

## Prioritized Fix List

Every defect found is listed below, numbered sequentially. Severity is in
brackets. Nothing here is release-blocking; Priority 1 items are those with a
real (if low-probability) behavioral or safety-net consequence, Priority 2 are
latent-guard / robustness / convention items, Priority 3 are documentation and
cosmetic nits.

### Priority 1 — Behavioral & Safety-Net Integrity

1. ✅ [Low] **`daydream/segment_worker.js:126–137`** — a `setResolution(...) === false` rejection `break`s silently instead of faulting, asymmetric with the `'init'` path's `initFailed`. The worker keeps stale `canvasW/canvasH/segRange`; on an *enlarging* resize its old-geometry frame passes the controller's generation fence and composites a wrong-size quadrant with no fault. The unchecked `engine.setEffect(msg.name)` in the `'setEffect'` handler is the same class. Fix: post a fault message (mirroring the init path) on rejection.
2. ✅ [Low] **`effects/ChaoticStrings.h:64`, `effects/FlowField.h:33–37`, `effects/HopfFibration.h:28–30`, `effects/DreamBalls.h:49–51`** — `EffectConfig` hardcodes `{.strobe = true}` without wiring `.full_frame = decltype(filters)::any_crosses_segments` as the sibling convention does. Correct today (their filters are history-free), but swapping in a history-bearing stage would silently clip a stateful effect to one segment band with no compile error and no test coverage. Fix: use the trait-derived flag.
3. ✅ [Low] **`core/engine/styles.h:90–122, 135–145`** — `Style` mixes lerpable preset data with per-instance *bound* state (the `noise` pointer + `hue_*` cache). `lerp()` deliberately preserves `noise`, but a full-struct copy (`Presets::apply`, `style = Style::Churn()`) resets it to `nullptr`, degrading `noise_warp` to identity. The two ways to mutate a `Style` diverge silently, and the doc ("POD-copyable — safe to store in `Presets<>` and lerp") invites the unsafe path. Fix: split bound state out of the lerpable payload, or document the bind-before-`apply` invariant on the struct.
4. ✅ [Low] **`core/engine/styles.h:77, 289`** — the stock `plain_fade` / `hue_fade` function *bodies* are dead on the render path: `Filter::Pixel::Feedback::flush()` (`filter.h ~1468–1488`) special-cases their addresses and runs its own inlined copy of the math. Editing `styles.h::hue_fade` silently no-ops for the built-in styles, with no test catching the divergence. Fix: have the fast path and the free function share one `static` helper so they cannot drift (the address is load-bearing and cannot be deleted).
5. ✅ [Low] **`core/render/filter.h:1148`** — `Screen::Blur::crosses_segments = true` is an unnecessary full-frame pessimization, inconsistent with the identical-situation `Screen::AntiAlias` (`= false`). Any pipeline containing a `Blur` stage disables segmented rendering, forcing every board to render the full 288×144 instead of its quadrant — a throughput loss on a real-time target. Fix: set `= false` and drop the stale comment. *Caveat:* validate against the segment driver's "render full geometry, clip output per band" model before landing.
6. ✅ [Low] **`.github/workflows/ci.yml:600`** (teensy-warnings) — the ratchet builds `-e holosphere -e phantasm` but not `holosphere_dma`, leaving `pov_single.h`'s `USE_DMA_LEDS` single-board branch un-ratcheted (the `teensy-size` job builds it but is object-cached, so it emits no warnings). A new `-Wall/-Wextra` warning in that slice compiles green everywhere. Fix: add `-e holosphere_dma` to the teensy-warnings `pio run` line.
7. [Low] **`scripts/wasm_smoke.mjs:19`** — `FRAMES_PER_EFFECT = 3` versus the native suite's `HS_SMOKE_FRAMES = 120`. The WASM smoke is the *only* runtime exercise of the shipped `.wasm` (embind seams, `ALLOW_MEMORY_GROWTH` detachment, 8 KB stack, per-frame arena high-water gate), but at 3 frames it never reaches late-lifecycle events (e.g. the frame-48 ShapeShifter cut, arena compaction), so a WASM-only regression tied to one can ride a green deploy. Fix: honor a `WASM_SMOKE_FRAMES` env override and set it to ~120 in the wasm job.

### Priority 2 — Latent Guards, Robustness & Convention

8. [Low] **`core/mesh/hankin.h:142–144`** — `compile_hankin`'s edge-midpoint `mid.normalize()` is unguarded, unlike the identical computation in `ambo` (`conway.h:498–499`, which uses `normalized_or`). A near-antipodal edge collapses the midpoint to the origin and traps instead of degrading gracefully. Unreachable on the current roster (short edges). Fix: `mid = normalized_or(mid, p_a);`.
9. [Low] **`core/mesh/mesh_classes.h:197–200`** — `build_mesh_class_bake`'s face-centroid `center.normalize()` is unguarded, unlike every Conway sibling that normalizes a centroid. A centrally-symmetric face (vertex sum ≈ 0) traps. Facility is currently unwired to a live effect. Fix: `center = normalized_or(center, mesh.vertices[idx[0]]);`.
10. [Low] **`effects/Comets.h:68–99`** — omits the `static_assert(FOOTPRINT_BYTES <= PERSISTENT_BUDGET)` compile-time budget guard its siblings (FlowField/GnomonicStars/HopfFibration) carry, despite using the default arena split and a compile-time-computable footprint (`OrientationTrail<…, TRAIL_LENGTH>` + baked LUT). Bumping `TRAIL_LENGTH` would overflow the device partition as a runtime trap rather than a build failure. Fix: add the sibling static_assert idiom.
11. [Low] **`effects/Dynamo.h`** (near the `TRAIL_CAPACITY` declaration, ~line 390) — same missing `FOOTPRINT_BYTES <= PERSISTENT_BUDGET` guard; current footprint ~126 KB / 298 KB so no live overflow, but `TRAIL_CAPACITY` is explicitly documented as a tunable and host tests inflate the arena to 8 MB, so a device overflow would surface only as a hardware trap. Fix: add the static_assert.
12. [Low] **`effects/MeshFeedback.h`** — lacks the persistent-budget `static_assert` its siblings (MindSplatter/SplineFlow/Voronoi) carry. Not a real overflow (the warp cache is ~40 KB against the ~324 KB partition and fail-fasts at `init()`), so this is a consistency/maintainability item completing the theme of items 10–11. Fix: add the guard for uniformity.
13. [Low] **`core/engine/platform.h:601–611`** — `map()` guards only `in_max == in_min`, then divides by `static_cast<int32_t>(in_max - in_min)`. On the LP64 native/test host a range difference that is a nonzero multiple of 2³² truncates to 0 → SIGFPE. Not reachable on the 32-bit device. Fix: compute the divisor once and guard it before dividing.
14. [Low] **`core/engine/util.h:84–88`** — `wrap(int, int)` lacks the `m == 0` device-parity guard its siblings `map()` and `addmod8()` carry; on the host `wrap(x, 0)` is a SIGFPE while the Cortex-M7 `SDIV` returns `x`. `m > 0` is a documented precondition, so this is a sim/device divergence rather than a live bug. Fix: `if (m == 0) return x;` (or treat the precondition as authoritative and close it).
15. [Low] **`core/animation/motion.h:54–68`** — `Path::get_point` interpolates adjacent unit-sphere samples linearly, returning a sub-sphere chord point (magnitude < 1) with non-uniform angular speed, inconsistent with the codebase's otherwise-spherical interpolation discipline. Impact is bounded (callers re-normalize), but the contract is unstated. Fix: slerp adjacent points, or document the non-unit return so callers normalize.
16. [Low] **`daydream/daydream.js:711–714`** — the `.rec-duration` div is appended to `#canvas-container` but not removed in `disposeApp()`, unlike every other owned resource. Harmless on a real `pagehide`, but breaks the file's meticulous teardown symmetry and would accumulate under SPA-style re-mounts. Fix: `durationEl.remove()` in `disposeApp()`.
17. [Low] **`tests/test_death.h:1175`** — `MIN_DEATH_CASES = 41` pins the roster *count*, not its *identity*; a remove-one/add-one edit keeps `n == 41` and silently drops the removed case's trap coverage. Every other roster in the repo is cross-pinned to an independent witness. Genuinely harder here (no independent source of truth), so this is a noted design tradeoff. Fix (optional): assert the exact set of case names.

### Priority 3 — Documentation & Cosmetic

18. ✅ [Low] **`core/color/color.h:820–832, 355`** — `hue_rotate_rgb`'s doc claims it "preserves perceived lightness exactly," but the forward transform uses `fast_cbrt` while the inverse uses exact cubes, so a null rotation is not an identity round-trip (~6 LSB at 16-bit) and can drift lightness if fed back per-frame. Fix: soften the comment to "…to within `fast_cbrt` accuracy" (an exact round-trip would need `cbrtf` forward, rejected on perf grounds).
19. ✅ [Info] **`core/color/color.h:1823–1830`** — `QuantizeModifier::modify` returns `roundf(t*s)/s`, yielding `s+1` distinct levels over `[0,1]`, not the `s` implied by "N steps." Purely a doc clarification.
20. ✅ [Nit] **`core/render/scan.h:308–316`** — `BoundingSphere::get_intervals` justifies its asymmetric `x_lo/x_hi` caps as reaching a full row "for odd W," but `TrigLUT` `static_assert`s `W % 4 == 0`, so the odd-W branch of the reasoning is unreachable dead rationale (no behavioral bug). Fix: drop the odd-W clause (or simplify both caps to `W/2`).
21. ✅ [Nit] **`core/render/shading.h:175, 179`** — bare `TOLERANCE` where the surrounding code uses `math::TOLERANCE`; compiles via a global alias but is inconsistent qualification. Fix: qualify for uniformity.
22. ❌ [Info] **`core/animation/sprites.h:75–104`** — increment-first stepping means the first *drawn* frame is `t == 1`, so a fade-in of duration `F` begins at `easing(1/F)` rather than ~0, and `fade_in_duration == 1` shows no fade at all. Consistent with the shared increment-first convention (not a bug); reported for completeness. No fix needed for correctness. **Rejected:** intended increment-first convention, no correctness impact.
23. ❌ [Info] **`daydream/driver.js:895–903`** — `coordsLabel()` allocates a fresh `THREE.Vector3` + template string per call, so an effect returning many labelled points allocates per-point-per-frame. Not fixable without changing the `getLabels` return contract (the position must outlive the call); awareness note only. **Rejected:** not fixable without changing the `getLabels` return contract.
24. ✅ [Low] **`daydream/tests/solid_codegen.test.js:53`** — a JSDoc comment claims a `relax` default of 100; the actual default is 8 in both the JS generator (`solid_codegen.js:145`) and the C++ it mirrors (`solids.h:319`), and the very next test asserts 8. A maintainer trusting the comment could "fix" the generator and diverge every emitted recipe. Fix: change "the default 100" to "the default 8."

---

## Notable Strengths (validated, not merely asserted)

- **The safety net is verified, not assumed.** A death harness asserts that
  `HS_CHECK` traps actually fire (`SIGILL` / `STATUS_ILLEGAL_INSTRUCTION`), so
  the fail-fast philosophy is tested rather than trusted.
- **Silent drift is structurally prevented.** Effect / resolution / death /
  module rosters are single-sourced via X-macros and cross-pinned against
  independent witnesses (`check_includes.cmake`, `check_effect_roster.mjs`, the
  WASM `HS_EFFECT_COUNT` startup trap, `--check-modules`).
- **The DMAMEM-on-template-static hazard is closed.** Both POV drivers emit an
  explicit specialization with strong linkage so DMA buffers land in OCRAM, and
  both targets actually invoke it — the exact footgun that silently misplaces
  memory on GCC was pre-empted.
- **Compile-time-legal pipelines.** The `Pipeline<W,H,Filters...>` chain
  enforces stage ordering (2D-after-3D, terminal-last, history-vs-replacing)
  with `static_assert`s, turning an entire class of composition errors into
  build failures.
- **The WASM/JS boundary is hardened.** Every untyped JS entry validates and
  clamps before reaching engine math; memory-view detachment across
  `ALLOW_MEMORY_GROWTH` is respected by pre-sizing view-backed vectors once.

---

## Cross-Cutting Themes for the Maintainer

1. **Convention completion.** Several findings (2 EffectConfig flags, 3 missing
   budget `static_assert`s, 2 unguarded `normalize()`s) are the same underlying
   issue: an established, correct convention applied *almost* everywhere. A
   short sweep to make each convention uniform (or a lint that enforces it)
   removes a whole latent-defect class.
2. **Comment density vs. drift.** The near-universal maintainability ding is
   rationale volume. The already-in-progress trimming of over-commented blocks
   is the right direction; the silent-divergence hazards (items 3–4) are
   precisely where a comment at the edit site would *not* help — the fix is
   structural (shared helpers), not more prose.
3. **Two CI coverage seams** (items 6–7) are worth closing because they are the
   difference between "the safety net caught it" and "it shipped green."
