# Holosphere — Code Quality Review

**Scope:** the full Holosphere C++ rendering engine (`core/`, `effects/`, `hardware/`, `targets/`), the
`daydream` web simulator (JS/WASM), the native test suite, and the build/CI/tooling. Reviewed by a
fan-out of subsystem auditors, with every concrete defect independently re-validated by a second,
adversarial pass before inclusion here.

**Out of scope** (by request): `effects_legacy.h`, the `*.ino` target shells, and `core/rotate.h`.

**Method:** ~16 independent subsystem reviews → consolidation → a second validation pass that opened
each cited line and returned CONFIRMED / REJECTED / REVISED. Roughly a third of the raw findings were
rejected or downgraded on validation (e.g. a claimed "silent arena rewind" actually fail-fast traps;
`pov_single` already caches the virtual it was accused of re-reading; a `wrap()` NaN comment was
accurate as written). Only verified items appear below.

---

## 1. Overall Verdict

This is an **exceptionally well-engineered project — top decile for its class** (embedded creative-coding /
real-time graphics). The dominant story of the audit is the *absence* of the defects one expects in a
60-KLOC C++/JS codebase: no memory-corruption paths, no use-after-free, no race conditions in shipped
logic, no correctness bugs in the core math. Edge cases that are normally landmines — the azimuth seam,
the sphere poles, antipodal vectors, arena exhaustion, 16-bit fixed-point overflow, the ISR/main-loop
hand-off — are explicitly guarded and, in most cases, tested. The defects that remain are overwhelmingly
**Low severity**: polish, documentation precision, perf-attention notes, and test-coverage gaps. There
are **no Critical or High correctness defects**. The two High-priority items are both *test-coverage
gaps* on the device-only surface (DMA transmit, the double-buffer atomics), not active bugs.

**Composite grade: A−**

---

## 2. Quality-Dimension Grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness** | A | Large surface, rigorously reasoned; remaining issues are edge-case/precision, never core-algorithm. Seam/pole/antipode/wrap handling is pervasive and mostly tested. |
| **Numerical robustness** | A− | Radicand/denominator floors, clamp-before-`acos`, NaN-blocking, antipodal guards everywhere; a few residual discontinuities (antipodal `slerp`) and unguarded `radius=0` divides. |
| **Memory safety** | A | Arena allocator with explicit `Arena&` threading, RAII `ScratchScope`, dual-stamp span-staleness detection, always-on OOB traps at cold seams. No unsafe path found. |
| **Concurrency / ISR-safety** | A− | Single-writer double-buffer model is correct and exhaustively documented; relaxed atomics justified by interrupt-disable barriers. Knocked down only because the atomic hand-off has no contention/TSan test (see #2). |
| **Error handling / fail-fast** | A+ | `HS_CHECK` (always-on under `NDEBUG`, cold-path-only, single predicted branch to `__builtin_trap`) plus a death-harness that asserts each trap actually fires (`SIGILL`). A genuinely exemplary safety doctrine, consistently applied. |
| **Performance** | A− | LUT-driven trig, compile-time `<W,H>` specialization, scanline interval culling, adaptive per-face LUTs, zero-copy WASM views. A handful of effects (SphericalHarmonics, GnomonicStars, Voronoi) carry uncapped or un-narrowed per-pixel cost. |
| **API / interface expressiveness** | A | The `Pipeline<W,H,Filters...>` with compile-time-`static_assert`-enforced domain ordering (3D→2D→terminal) is a centerpiece worth an A+ on its own. The explicit-`Arena&` operator threading is powerful but a documented foot-gun (composition polarity). |
| **Architectural elegance** | A | Clean three-domain filter split (World/Screen/Pixel), CRTP animation bases, X-macro single-source rosters, pure host-testable cores severed from device shells (`pov_*_map`, `hd107s_frame`, `param_marshal`, `segment_layout`). |
| **Readability** | A | Consistent naming and idiom; comments explain *why*, not *what*. Comment density occasionally tips into verbosity, but it earns its keep given the no-debugger target. |
| **Documentation** | A+ | Among the best-documented codebases of its size — a 150 KB architecture README plus load-bearing doxygen with rationale at nearly every non-obvious seam. Generators-of-record documented for every generated file. |
| **Testing** | A− | Oracle-driven (Euler χ, brute-force k-NN, energy conservation, hand-computed wire bytes, device-asm parity), a standout death harness, 4-way sharding with a coverage meta-gate, ASan/UBSan, determinism-under-perturbation. Gaps: device DMA path, double-buffer race, several effects' white-box pins, a noise oracle. |
| **Build / CI / tooling** | A | Pinned toolchains throughout, generators-of-record with CI provenance diffs, a self-testing Teensy size/layout gate, FF-only ref hook, runtime WASM smoke across every effect at both resolutions, both Teensy images size-gated in CI. |
| **Portability** | A− | Three platforms (Teensy/WASM/desktop) behind a disciplined `platform.h`; `uint32_t` index parity, finite-math `#error` guards, fast-math asymmetry documented and matched. Minor signed/unsigned mixing. |

---

## 3. Per-Component Grades

| Component | Grade | Notes |
|---|---|---|
| `core/` math & geometry (`3dmath`, `geometry`, `util`, `waves`, `easing`) | A | Defensively written; pole/antipode/wrap care; antipodal `slerp` discontinuity is the only real flaw. |
| `core/` color (`color`, `palettes`, `color_luts`) | A | 16-bit linear blending + OKLCH interpolation, cached stops, packed SIMD-ish adds. Generator-of-record verified. |
| `core/` memory & platform (`memory`, `static_circular_buffer`, `platform`, `presets`) | A | Arena partition math, ArenaSpan staleness stamps, fail-fast placement all correct. |
| `core/` rasterizers (`sdf`, `scan`, `plot`) | A | Distances-in-radians SDFs, CSG, geodesic/planar curves; seam/pole guards thorough. `radius=0` divide in `Flower`. |
| `core/` mesh (`mesh`, `conway`, `spatial`, `solids`, `hankin`, `generators`) | A− | Half-edge invariants, Conway operator math, KD-tree all sound; explicit-arena polarity is a foot-gun; a couple of test gaps. |
| `core/` pipeline (`filter`, `animation`, `canvas`, `transformers`, `styles`) | A | The `Pipeline`/`Filter` design is a highlight; double-buffer Canvas correct by design. |
| `effects/` (28 effects) | A− | Idiomatic engine use, good arena discipline; some perf-attention and duplication (trail skeletons, `strobe_columns` boilerplate, `Params::lerp` bodies). |
| `hardware/` drivers | A | Flywheel timebase, symbol codec, epoch/beacon sync, fail-dark latch — rigorous embedded work; pure cores host-tested, DMA shell is not. |
| `daydream` JS simulator | A− | WASM-heap detach contract, Three.js disposal, recorder/listener teardown all correct; a few SoC/coupling smells. |
| Segmented-POV worker pool | A− | Generation fence, fault latch, transferable-buffer discipline all sound; one param-publish ordering inconsistency. |
| Test suite | A− | See Testing dimension; the two High gaps live here. |
| Build / CI / tooling | A | See Build dimension; only polish-level papercuts. |

---

## 4. Prioritized Fix List

Items are numbered sequentially. Severity reflects post-validation consensus. There are no P0 (Critical)
items — no shipped correctness or safety bug was found.

### Priority 1 — High (test-coverage gaps on the no-debugger surface)

1. ✅ **`hardware/dma_led.h` — DMA transmit path has zero host-test coverage.** The entire file is
`#ifdef ARDUINO`, so `TeensySPIDMA`/`DMALEDController` (chunking, completion, SPI clock-out, overrun
counting) compile and run only on-device, where there is no debugger or console. Only the protocol
*buffer* (`hd107s_frame.h`) is host-tested. Extract the pure framing/chunk/completion math into a
host-testable header (the same split already applied to `pov_segment_map`/`pov_single_map`/`hd107s_frame`)
and pin it.

2. ✅ **`tests/test_canvas.h` — the double-buffer atomics are untested under contention.** The three
`std::atomic<int>` indices (`cur_`/`next_`/`prev_`) are the heart of the "Why the ISR Double Buffer"
invariant, but the existing tests are single-threaded or a deterministic spin-release handshake
(`test_double_buffer_handoff_no_aliasing` is explicitly single-threaded). Add a ThreadSanitizer job or a
producer/consumer hammer loop asserting no torn or doubly-displayed frame.

### Priority 2 — Medium (correctness-adjacent, perf, and consistency)

3. ✅ **`effects/Voronoi.h:160-217` — coarse-coherence silently drops sub-block cells.** The 8 px coherence
block classifies only at corners; a Voronoi cell smaller than the block missed by all four corners
vanishes. At high site counts (up to ~400) small cells near the dense pole region can disappear. The
in-code "bit-identical to the full query" equivalence claim (~L208-217) is **untested**. Add the
equivalence unit test and, ideally, a tighter block or adaptive fallback near the poles.

4. **`effects/SphericalHarmonics.h:192-220` — heaviest per-pixel eval with no horizontal narrowing.**
`HarmonicField::distance` evaluates `associatedLegendre` twice per pixel (current + blend target),
returns `false` from `get_horizontal_intervals`, and uses full-height vertical bounds. Bounded only by
`MAX_DEGREE=4`. At 288×144 this is ~41k pixels × 2 full Legendre recurrences per frame. Add interval
narrowing or document the cost ceiling.

5. **`effects/GnomonicStars.h:111-120` — uncapped per-frame per-point work.** Every frame transforms up
to 2000 points through `transformer.transform()` and issues a separate `Scan::Star::draw` per point, with
no per-frame work cap analogous to the `kMaxRings` clamp used elsewhere. Add a screen-space cull / draw
cap, or document the worst-case cost.

6. ✅ **`core/led.h:62` — `correction_guard_depth()` lacks a re-entrancy assertion.** It is a bare
`static int` documented as main-loop-only, but unlike the Canvas/Timeline/Effect singletons it carries no
`HS_CHECK` against the wrong-context use the rest of the codebase guards. Add the assertion (or a note on
why it cannot exist) to match the project's fail-fast rigor. **Resolved (already guarded):** re-entrancy —
a second live guard of either type — already traps via `HS_CHECK(correction_guard_depth() == 0)` in both
guard constructors; the main-loop-only constraint follows the codebase's note-only convention for that same
class of precondition (e.g. `hs::random()` in `platform.h`), as no ISR-context primitive exists to assert
against.

7. ✅ **`core/3dmath.h:1216` — `slerp(Vector)` is non-monotone near the antipode.** The antipodal branch
returns `normalized_or(v1 + (v2-v1)*t, v1)`; at `t≈0.5` with `v2≈-v1` the blend collapses to ~0 and falls
back to `v1`, so `slerp(p, antipode, 0.5) ≈ p` regardless of `t`. Pick a perpendicular-bisector axis for
the degenerate case, or document that antipodal callers must supply an explicit rotation axis.

8. ✅ **`daydream/segment_controller.js:236` — `paramValues` is published outside the dedup guard.** Seg-0's
`paramValues` is mirrored above the `_frameSeen[segId]` check while `pending` is mirrored inside it; a
worker emitting two `frame` messages in one generation re-publishes params twice. Benign (idempotent
overwrite) but inconsistent — move the param mirror below the `_frameSeen` guard.

### Priority 3 — Low (robustness, polish, documentation, duplication)

9. **`core/sdf.h:2685` — `Flower` (and `SphericalPolygon`) lack `HS_CHECK(radius > 0)`.** `Star` has it;
`Flower` does not, and `radius=0` drives `t = scan_dist/thickness` to infinity. Add the guard for parity.

10. ✅ **`core/color.h:826` — `lerp_oklch` mixes `__builtin_fmaxf` with `hs::clamp`.** Use one idiom
(`hs::clamp` / `std::fmax`) for portability and consistency with the rest of the file.

11. ✅ **`core/color.h:959` — `Gradient` segments redundantly double-write boundary indices.** Adjacent
segments both write the shared join index (`for i = start; i <= end`). Correctness-neutral (identical
color) but a transposed edit could silently double-fill; make the interior loop `i < end` and let each
segment own its left endpoint.

12. ✅ **`core/memory.h`/`platform.h` `random_to_unit` — no assert on the `max == UINT32_MAX` precondition.**
The clamp constant is derived for that exact divisor, but only `rand_f`'s `static_assert` enforces it; a
direct caller with a different `max` gets a silently-wrong top band. Add an assert to the pure function
(it is advertised as independently unit-testable).

13. ✅ **`core/presets.h:103-104` — `current_idx`/`prev_idx` are `int` while `Size` is `size_t`.** Gratuitous
signed/unsigned mixing in the modular index math; make them consistently one type.

14. ✅ **`core/animation.h:1208-1218` — `Sprite::step` has unguarded fade durations.** `fade_in_duration`/
`fade_out_duration` are cast to `uint32_t`; a negative value becomes a huge bound and the fade never
completes. The float `duration` overloads `HS_CHECK(duration >= 0)`; add the same guard for the fades.

15. **`core/animation.h:297-315` — `.then()` silently no-ops on perpetual animations.** A `.then()`
attached to an indefinite, non-repeating animation (e.g. `RandomWalk`, `MobiusWarpEvolving`) compiles and
never fires; only prose protects this. Add a compile-time or runtime trap when `.then()` is attached to a
type that can never complete.

16. ✅ **`core/conway.h:947-951` — `relax()` undercounts the boundary 1-ring.** It keys on `hev.half_edge`
(only the *last* incoming half-edge written during build); if that edge is a boundary
(`pair == HE_NONE`), the entire 1-ring spring force is skipped even when other incoming edges have twins.
The closed-manifold roster never triggers it, but scan for any paired edge or document the limitation.

17. ✅ **`core/mesh.h:619` (`classify_faces_impl`) — 32-bit hash used directly as a topology-class id.** An
`fmix32` collision merges two distinct face topologies into one palette class. Safe for the fixed roster,
but add a roster-wide collision-free test (or a debug assert) so it stays safe as the roster grows.

18. **`core/filter.h:1310-1319` — `Feedback::any_pixel_lit` is an O(W·H) pre-scan.** A lit feedback frame
pays a full-canvas scan *plus* the full warp pass (two passes). Fold the lit-test into the warp loop's
first row, or accept and note the cost.

19. ✅ **`core/transformers.h:376-385` — noise decorrelation comment overstates the mechanism.** The
`ny`/`nz` channels differ from `nx` by a *constant 3D translation*, not "distinct per-axis offsets";
decorrelation relies on the translation magnitude (100/200) exceeding the noise correlation length.
Tighten the comment.

20. ✅ **`effects/IslamicStars.h:100-104` — ripple pool can silently drop spawns.** The `static_assert` is a
heuristic floor; the true simultaneous-ripple peak can exceed `kRipplePoolSize`, making `spawn()` return
null and the Burst / Ripp-Dur sliders partially non-functional at their tops. Size the pool to the real
peak or clamp the slider tops.

21. ✅ **`effects/PetalFlow.h:245` — the opacity cull is alpha-dependent.** Gating on
`opacity * params.alpha <= 0.01` means lowering the Alpha slider shortens the flow's pole reach. Gate on
geometric `opacity` before multiplying by `params.alpha`.

22. ✅ **`effects/Thrusters.h:241` — redundant `frame % 32`.** `t_global` is already wrapped to `[0,32)` at
line 80; the inner modulo in `ring_fn` is dead. Remove it (or keep it only if the function is meant to be
called with unwrapped frames, and document that).

23. ✅ **`effects/RingSpin.h:50` — dead `X_AXIS` default normal.** `spawn_ring` always placement-news with
`Y_AXIS`, so the constructor's `X_AXIS` default is never used. Benign; drop it for clarity.

24. ✅ **`hardware/pov_segmented.h:155-164` — master "longest coast" counter is meaningless.**
`max_coast_halves` increments every fold, but the reset (`halves_since_snap_ = 0`) lives in
`handle_burst`, which the master (the reference board, receiving no wire bursts) never reaches — so the
master reports an ever-growing coast. Reset it for the master, or exclude the master from the metric.

25. **`daydream/gui.js:181-188` — `_attachUrlWriter` couples to lil-gui internals.** It consumes lil-gui's
single `onChange` slot via a fan-out registrar; `reset()` reads `_onChange` directly. Document the
coupling or wrap it so a future lil-gui upgrade can't silently break the URL writer.

26. **`daydream/daydream.js:564-569` — WASM-load `.catch` forward-references `testAll*`.** The catch
references `testAllInterval`/`testAllController` declared ~40 lines later; it works because the catch is
async, but it is a TDZ fragility. Move the Test-All setup above the WASM-init block.

27. ✅ **`justfile:45-48` — `just docs` is broken on a fresh clone.** It runs `doxygen Doxyfile.local`, but
`Doxyfile.local` is gitignored and untracked and no recipe generates it. Either track it or have the
recipe synthesize it (as `docs.yml` does).

28. ✅ **`scripts/generate_luts.py` — the `--check` self-test is never run in CI.** The `lut-provenance` job
only diffs numeric tokens. Add `python3 scripts/generate_luts.py --check` for a clearer monotonicity /
round-trip signal.

29. ✅ **`platformio.ini:149-153` — `[env:phantasm]` does not inherit `[env].extra_scripts`.** PlatformIO
does not merge list options across sections, so a future script added to `[env]` would silently apply only
to holosphere. Use `${env.extra_scripts}` interpolation or add a warning comment.

### Priority 3 (continued) — Test-coverage gaps (Low)

30. ✅ **`tests/test_geometry.h` — south-pole overflow boundary not pinned.** `phi_to_y`/`vector_to_pixel`
can "land a hair above `H-1`" (documented in `geometry.h:226,477`); the existing test uses
`HS_EXPECT_NEAR` (proximity), which passes even at 144.0001 and does not pin the `≤ H_VIRT-1` upper bound.
Add a bound assertion.

31. ✅ **`tests/test_sdf.h` — `Union`/`SmoothUnion` interval merge path untested.** Subtract/Intersection
seam-straddle and AngularRepeat culls are covered, but the `Union`/`SmoothUnion`
`get_horizontal_intervals` overlap/merge path is exercised only via `distance()`, never with
seam-straddling intervals — exactly where double-paint regressions hide.

32. ✅ **`tests/test_conway.h` — `snub(cube)` has no topology-count assertions.** The most complex operator
(chiral, twist pass, Newell normals) is checked only for basic winding/Euler invariants; a winding-correct
but wrong-count snub would pass. Assert expected vertex/face/triangle counts.

33. ✅ **`tests` — `palettes.h` generative/OKLCH layer is untested.** `lerp_oklch` math is pinned, but the
named `ProceduralPalette` instances, shortest-arc hue direction, and `MeshPaletteBank` have no dedicated
coverage. Pin a few named-palette endpoints and hue-arc directions.

34. ✅ **`tests` — noise paths lack an oracle.** `noise_transform` / `Feedback::noise_warp` / the
FastNoiseLite wrapper are only checked for divergence/NaN, never against a golden reference; a silent
noise regression renders deterministically and passes. Golden-hash a fixed noise sample grid.

35. ✅ **`tests/test_effects.h` — several effects have smoke/determinism coverage only.** FlowField, Flyby,
PetalFlow, ChaoticStrings, MindSplatter, RingSpin (and the white-box-light Hopf/Islamic/Raymarch/Liquid2D)
lack invariant pins for exactly the silent-drift behaviors smoke tests miss (spawn-gap accumulators,
emit-phase wrap, pattern-arg clamps). Add targeted pins.

### Priority 4 — Nits / elegance (informational)

36. ✅ **Cross-effect duplication.** `strobe_columns()`/`needs_full_frame()` doc-comment blocks are
byte-identical across ~9 effects; the trail skeleton (record → `deep_tween` → fragment shader) is
copy-pasted across Comets/RingSpin/ChaoticStrings with a "fix must be hand-propagated" note; two
field-by-field `Params::lerp` bodies (Flyby, MindSplatter) drift-risk on field reorder. Consider a
`StrobedEffect` mixin, a `TrailRenderer` helper, and a designated-aggregate lerp.

37. **Explicit-`Arena&` composition polarity (`conway.h:312-359`).** The 50-line "scratch contract /
composition polarity" comment encodes correctness the type system doesn't enforce (a composed op returns
in the *opposite* arena). A thin `OpResult{ PolyMesh, Arena& where }` wrapper would make polarity
self-documenting at call sites. Subjective; the length of the comment is itself the evidence.

38. ✅ **`make_basis` / `least_parallel_axis` duplication (`3dmath.h:1054`, `geometry.h:824`,
`MobiusGrid.h:239-246`).** The "pick a cross-safe seed axis" idiom is re-implemented inline at several
sites instead of routing through the helper; consolidate or rename the helper to `cross_safe_axis`.

---

## 5. What This Codebase Does Better Than Most

Called out explicitly because a defect list under-represents the engineering quality:

- **The fail-fast doctrine is real, not aspirational.** `HS_CHECK` is always-on under `NDEBUG`,
cold-path-only, and *verified* by a death harness that re-execs a child and asserts the exact `SIGILL`
shape. Few embedded projects test that their assertions fire.
- **The `Pipeline<W,H,Filters...>` abstraction** enforces filter ordering (3D-before-2D, terminal-last,
unit-input contracts) at compile time via `static_assert`. This is the correct place to catch composition
errors, and it is done cleanly.
- **Compile-time `<W,H>` specialization** eliminates runtime generality cost on the MCU while letting the
simulator instantiate 288×144 — one codebase, zero-overhead device builds.
- **Provenance discipline.** Every generated artifact (`color_luts.h`, `reaction_graph.cpp`, the WASM
hashes) has a generator-of-record and a CI diff gate; the Teensy size/layout gate is itself a
self-testing, stdlib-only Python tool with good+broken fixtures.
- **Host-testable cores severed from device shells** (`pov_*_map`, `hd107s_frame`, `pov_sync`,
`param_marshal`, `segment_layout`) — the pattern that makes the otherwise-untestable hardware logic
verifiable, and the basis for the two High-priority recommendations above (extend it to the DMA path).
