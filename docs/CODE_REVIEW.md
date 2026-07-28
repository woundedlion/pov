# Holosphere + daydream — Code Quality Review

**Date:** 2026-07-27
**Scope:** the `Holosphere` C++ engine/firmware repository and the `daydream` web simulator repository, reviewed as one product.
**Excluded by request:** `core/engine/effects_legacy.h`, `core/math/rotate.h`, `targets/Holosphere/Holosphere.ino`.
**Method:** sixteen independent reviewers, each scoped to one subsystem, each required to validate every candidate finding against the cited code and against the existing test suite before reporting it. Findings that could not be substantiated were dropped and are recorded in §7 so they are not re-litigated.

---

## 1. Executive summary

This is an exceptionally well-engineered codebase. Across roughly 100,000 lines of reviewed C++ and JavaScript, **no Critical defect was found** — no memory-safety hole, no data race, no correctness bug capable of corrupting output on the shipping path. That is an unusual result at this size, and it is not an artifact of shallow review: reviewers traced ISR field ownership by hand, computed reaction-diffusion stability bounds, verified OKLab matrices against published references, and audited arena lifetime through every Conway operator.

The discipline producing that result is visible everywhere and is the project's defining characteristic:

- **Invariants are executable, not aspirational.** Compile-time budget algebra rejects an over-nested SDF (`sdf.h:83`); a `static_assert` pins preset field order by resolved name so a struct reorder fails the build (`styles.h:286`); the composition-polarity contract is asserted by arena address range in a host test (`test_conway.h:800`); an empty test module is counted as a *failure* (`test_harness.h:244`).
- **Comments record rejected alternatives and measured negative results** — the hard and rare part. `hankin.h:290` documents why the chord ratio was the wrong discriminator. `color.h:746` documents why a plain bisection converges on the wrong side of a tolerance gap. `platform.h:1168` records that three standard libraries produced three different `std::shuffle` draw counts.
- **The generated-artifact story is complete.** All five generated files are CI-gated; the strongest gate re-derives the asset rather than diffing it.
- **The pure/shell split in the firmware is real.** The 1-wire sync protocol, the DMA double-buffer maths, the segment index maps and the embind predicates are all Arduino-free and host-tested, leaving the device ISR as ~90 lines of marshalling.

Weaknesses are correspondingly narrow and cluster into four themes:

1. **Unenforced side contracts.** Several subsystems require an initialisation call the type system does not demand — `Feedback::init_storage`, `Trails::init_storage`, `DirectAntiAliasSink::prepare`, `TransformerPool::prepare_frame`. Omission is silent: a performance cliff, an invisible effect, or a write through a null framebuffer base. The project already has the right pattern (`transformers.h:187` traps); it is simply not applied uniformly.
2. **A build-hygiene gap at the device boundary.** `NDEBUG` is defined in a header rather than in build flags, making device assert-stripping include-order dependent — and one header violates the ordering convention.
3. **Test-harness reach narrower than it appears.** Sanitizer and warning flags reach 1 of 9 test executables; the deploy gate runs a weaker tier than presubmit; and the excellent host simulator for the sync protocol structurally cannot model two of the hardware findings.
4. **Accessibility and second-repo CI conventions in `daydream`** lag the standard the rest of the project sets.

---

## 2. Quality dimensions — overall grades

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness & Robustness** | **A−** | No Critical defect across 16 scopes. Degenerate cases (poles, antipodes, seams, zero-length, singular Möbius) are systematically enumerated and handled at the site. Deductions: the unenforced side contracts (items 10–12, 23) and one direction-blind half-edge pairing that accepts malformed input silently. |
| **Architecture & Modularity** | **A−** | Compile-time resolution parameterisation, a variadic domain-lifting filter pipeline, a CRTP animation base, and a policy-based Segue family are strong designs that carry invariants in the type system. Held back by three oversized units (`OpLeg` at 78% of a 2750-line file, `IslamicStars.h` at 1511, `platform.h` at 1781 carrying six concerns) and one ~770-line subsystem with no production caller. |
| **Interface Expressiveness** | **B+** | Excellent where contracts are typed: `StoredFunctionRef` encodes a lifetime rule; deleted rvalue overloads make borrow violations compile errors; `spawn` vs `spawn_pinned` names the retention rule. Dragged down by two `Arena&` parameters whose roles invert between operator classes, `DistanceResult`'s overloaded register file (no shader can be written generically across shapes), and `thickness` meaning two different things. |
| **Readability & Naming** | **A** | Consistently intent-revealing names, near-zero comment noise, comments that explain *why not* and cite measurements. |
| **Documentation** | **A−** | The 2235-line README is genuinely architectural, and §7.10's 1-wire section is datasheet-grade (waveforms, AC timing, failure-mode matrix). Doxygen coverage is near-universal. Deductions: ~14 small doc/code divergences, several in tables that omit shipped types. |
| **Testing** | **A−** | ~50k lines across 67 C++ test files plus 35 JS suites, with real oracles rather than smoke: a fork/exec death harness, a 65536-input exhaustive sRGB sweep, a brute-force k-NN oracle, an energy-conservation check, a multi-board sync simulator with crystal skew and EMI injection. Deductions: sanitizer/warning flags reach one target; no coverage instrumentation; the Three.js lifecycle — the highest-risk JS — is untested. |
| **Performance Engineering** | **A** | Optimisation decisions are measured and documented. `HS_O3` regions derived from image sizes, DTCM split-decode tables replacing a 64 KB LUT, q16 blends, fused rasterizers pinned bit-exact against unfused twins. Deductions: two per-frame paths marked `cold`, a full-resolution ring sampler, an unused subsystem sizing hot-path scratch. |
| **Error Handling** | **A−** | Fail-fast is coherent and, unusually, *verified* — the death harness asserts the traps fire. Deductions: the device is the lenient side for an unbound `Fn` (inverting the intent), and six recorder failure paths never reach the user. |
| **Concurrency & Real-Time Safety** | **A−** | The single-writer ISR claim was traced field-by-field and **holds**; no 64-bit value crosses an ISR boundary; the 32-bit counter wrap is structurally unobservable and proven by test. The worker-pool generation fence is airtight. Deduction: no NVIC priority is set anywhere, silently degrading the edge ISR's timestamp. |
| **Portability** | **B+** | Divergences are catalogued in a ledger and mostly tested. Deductions: `NDEBUG` in a header, `static constexpr` at namespace scope ODR-used from templates, one hardcoded drive letter. |
| **Consistency** | **B+** | One arena idiom, one tolerance vocabulary, one HS_CHECK policy, applied broadly. Deductions: sibling effects diverge in optimisation annotation, budget-assert discipline, and the `full_frame` fold; degenerate-face policy splits between trap and silent-skip across operators. |
| **Accessibility** (daydream) | **C+** | The sidebar is genuinely well done — roving tabindex, correct ARIA roles, keyboard tests. But the stylesheet defines **no focus indicator at all** across 677 lines, the interactive 3-D viewport is announced as a static image with no keyboard path, and the page has no landmarks or `h1`. |
| **Maintainability** | **A−** | Drift guards are pervasive: roster cross-checks, include checks, count oracles, registry index pins. Held back by the oversized units and by tuning constants whose validity is prose rather than assertion. |

**Overall: A−.** Professional-grade work that would pass review at a strong systems shop, with specific narrow gaps rather than systemic weakness.

---

## 3. Per-subsystem grades

| Subsystem | Corr. | Arch. | Iface | Doc | Test | Perf | Overall |
|---|---|---|---|---|---|---|---|
| `core/engine` foundation | A− | A | A | A | A− | A | **A−** |
| `core/math` + transformers | A− | A | A− | A | B+ | A | **A−** |
| `core/render` sdf + scan | A− | B+ | B | A− | A | A | **A−** |
| `core/render` filter + plot | A− | A | B+ | A | A− | A | **A** |
| `core/color` | A | A− | A | A− | A | A+ | **A** |
| `core/mesh` core + Conway | A− | B+ | C+ | A | A− | A− | **B+** |
| `core/mesh` solids + Hankin + spatial | A− | A | A− | A | A | B+ | **A−** |
| `core/animation` | A− | B | A− | A− | A− | A− | **A−** |
| `effects/` mesh-heavy (6) | A− | B | A− | A− | A− | A | **A−** |
| `effects/` simulation (6) | A− | B+ | A− | A− | A− | B | **A−** |
| `effects/` shader + curve (13) | A | A− | A− | A | B+ | A | **A−** |
| `hardware/` drivers + sync | A− | A | A− | A | A− | A− | **A−** |
| Build / test / CI / tooling | — | — | — | A | A | — | **A−** |
| daydream worker pool + bridge | A | A | A− | A | A | A− | **A** |
| daydream renderer + UI | A− | B+ | B | A− | B | A− | **B+** |
| daydream tools + JS tests | A− | A− | A− | A | B+ | — | **A−** |
| daydream state / URL sync | B | — | B+ | — | B+ | — | **B** |

---

## 4. Prioritized findings

No Critical findings. Numbering runs sequentially across all priority sections.

### HIGH

1. ✅ **Half-edge pairing ignores edge direction, so an inconsistently-wound mesh is accepted silently** — `core/mesh/mesh.h:124-129,154-170,308-311`. `fill_edge_record` stores only the canonical `(min_v, max_v)` key, discarding direction, and `pair_half_edges` links a run of two with no orientation test. Two adjacent faces wound the *same* way therefore pair successfully; `require_closed_manifold` passes; every downstream `vertex_orbit` walk traverses a broken orbit and emits geometrically wrong faces with no trap. This contradicts the fail-fast doctrine in a file that already traps on three other malformations. Fix: one `uint16_t` comparison in `build_half_edge_mesh`'s `set_pair` lambda — not in `pair_half_edges`, which `classify_faces_impl` calls with synthetic self-unique keys. Cold path.

2. ✅ **`NDEBUG` is defined in a header, making device assert-stripping include-order dependent — and `memory.h` violates the convention** — `core/engine/platform.h:103-106`, `core/engine/memory.h:12,15`, `core/engine/memory.cpp:1`. `NDEBUG` is not in `platformio.ini`'s build flags. `<assert.h>` is deliberately re-includable and redefines `assert` per inclusion, so liveness depends on include order. Four headers document and follow the rule; `memory.h` includes `<cassert>` *before* `platform.h`, and `memory.cpp` includes `memory.h` first — so that TU compiles six unguarded asserts, including `ArenaVector::operator[]`, with a live bounds branch and a call to newlib `__assert_func`, pulling in exactly the stdio the define exists to exclude. It is also a formal ODR violation (IFNDR). Fix: add `-D NDEBUG` to `platformio.ini` and reorder `memory.h`'s includes. Confirm via a link-map check or a `-D NDEBUG` A/B on image size before treating the severity as final.

3. ✅ **`HS_SANITIZE` reaches 1 of 9 test executables** — `tests/CMakeLists.txt:84-100`. The flags gate `run_tests` only. The `sanitizers` CI job configures with `-DHS_SANITIZE=ON` then runs `ctest` over a shard regex including `unit_relax_bake_verify`, `unit_h_offset_renorm`, `unit_fastmath_clamp`, `unit_stack_budget`, `unit_arena_budget` — five **unsanitized** binaries inside a job whose entire purpose is ASan/UBSan coverage. Fix: hoist into a `_hs_sanitize(target)` function; apply at minimum to `relax_bake_verify` and `h_offset_renorm_check`. Deliberately exclude `stack_measure`/`arena_measure` and say why — ASan redzones would invalidate their measurements.

4. ✅ **The 25 named palette coefficient sets are a hand-copy with no engine pin; the "parity" test compares a mirror to a mirror** — `daydream/tools/palette_math.js:109-135`. All 25 palettes' `a/b/c/d` vec3s are transcribed from `core/color/palettes.h:15-113`. They currently match exactly, but nothing enforces it: the only test asserting engine correspondence (`tests/palette_math.test.js:71-79`) compares *names* against a second hard-coded list inside the test file. A coefficient edit in `palettes.h` silently desyncs both the gallery preview and the C++ the tool tells you to paste back — the exact failure a tool of this kind exists to prevent. Fix: export a WASM enumerator over `Palettes::` (name + 12 floats) and `deepEqual` it in `color_parity_wasm.test.js`.

5. ✅ **The op-dispatch table is duplicated between the live preview and the C++ generator with nothing pinning them** â `daydream/tools/solids.html:1653-1664` vs `tools/solid_codegen.js:125-172`. `applyOp` (what you see) and `generateFuncAndRecipe` (what you paste) are two independent per-op parameter tables. `engine_contract_wasm.test.js:138` pins op *names* against the engine, but nothing pins the two tables' arity or argument order to each other. Adding a parameter to an op and forgetting one side yields a generated solid that differs from the preview, silently, and the generated C++ still compiles. Fix: move `applyOp`/`OP_DEFS` into `solid_codegen.js` and assert both paths consume the same param key set for every `KNOWN_OPS` entry.

### MEDIUM

6. ❌ **The entire CSG combinator layer has no production caller, yet its worst case sizes hot-path scratch** — `core/render/sdf.h:1131-1899`, consumed at `core/render/scan.h:164-178`. `Union`, `SmoothUnion`, `Subtract`, `Intersection`, `AngularRepeat` (~770 lines) are instantiated only in tests and the README. `scan_region`'s per-call scratch is dimensioned from the CSG worst case (66 + 132 interval pairs ≈ 1.6 KiB from `scratch_arena_b`), a 33× headroom no shipping shape can reach — every draw of every shape pays it. **Decide this first**: if CSG is a planned public surface, fix items 7, 8 and 57; if not, deleting the layer resolves all four at once and shrinks the scratch. **Rejected:** CSG is an intended public surface for future effects; items 7, 8 and 57 are fixed instead.

7. ✅ **`Intersection` is missing the mixed-solidity `static_assert` its siblings carry** — `core/render/sdf.h:1624-1634`. `Union` and `SmoothUnion` both reject `A::is_solid != B::is_solid` with the reason spelled out. `Intersection` does not, so `Intersection<PlanarPolygon, Ring>` compiles, reports `is_solid == false`, and takes the stroke AA branch; when the solid child wins the `max`, `result.size` is its apothem and the entire interior renders as a gradient. Latent. Fix: one `static_assert`.

8. ✅ **`Subtract` inherits the subtrahend's `size`, so the carve edge anti-aliases over the cutting shape's radius** — `core/render/sdf.h:1600-1609`. When `-d_B` wins, the whole `res_b` is copied including `size`. For the README's own advertised `Subtract<Ring, PlanarPolygon>`, the stroke AA width becomes the polygon's apothem rather than the ring's thickness — a fade band as wide as the polygon instead of a one-pixel carve. Fix: preserve the minuend's metric across the swap.

9. ✅ **`DistortedRingStack`'s per-pixel ring window is one ring too wide at each end** — `core/render/scan.h:605-606`. The exact window is `ceil(...)−1 … floor(...)−1`; the code uses `−2` and `−0`. Still correct (each candidate is rejected by its own exact cosine band) but roughly doubles the inner-loop trip count in an `HS_O3` region. The asymmetry reads as an off-by-one rather than a margin, since `b_win` already carries the `fast_acos` slop. Fix: `−1` on both lines, gated by `tests/test_scan.h:429`.

10. ❌ **`Pixel::Feedback::init_storage()` omission is an undetectable performance cliff** — `core/render/filter.h:1747,2093`. Without it `cached_warp_x` stays null, `cacheable` is false, and every `flush()` re-populates the entire spherical control lattice forever. Output is pixel-identical, so no test, assert or log can distinguish it; only a device profile can. Fix: `HS_CHECK(cached_warp_x, ...)` at the top of `flush()` — once per frame, on a path that already scans the full canvas. **Rejected:** uncached flush is a supported mode, not an omission — `tests/test_filter.h:1986` deliberately runs a storage-free pipeline as the cached-vs-uncached parity oracle, and ~15 other flush tests plus death case `feedback_downsample_indivisible` would abort on the trap.

11. ✅ **`World::Trails` / `Screen::Trails` without `init_storage()` silently render no trails** — `core/render/filter.h:988,1383`. Both seed paths are gated on the arena pointer, so the effect looks alive but never grows a tail, and `size()` reads 0 — indistinguishable from "trails just expired". Fix: `HS_CHECK` at the top of each `flush()`.

12. ✅ **`DirectAntiAliasSink::plot` without `prepare()` writes through a null framebuffer base** — `core/render/filter.h:1264,1292,1315-1328`. `base_` defaults to `nullptr` and the only protection is an `assert` compiled out in release. On the Teensy 4.0, address 0 is ITCM — *writable* FlexRAM — so a forgotten `prepare()` corrupts code memory instead of hard-faulting, the opposite of the project's fail-fast policy. The same assert is also the only thing catching a stale `base_` prepared against a previous frame's buffer. Fix: add `bool prepared_for(const Canvas&)` and `HS_CHECK` it at the two `pipeline_direct_raster_path` dispatch sites — one compare per `draw()`.

13. ✅ **Feedback's warp offset differences an exact lattice coordinate against a fast-trig projection** — `core/render/filter.h:1797-1799`. `point` is the lattice's exact rational coordinate; `projected` comes from `fast_atan2`/`fast_acos`. The subtraction assumes the `vector → pixel → vector` round-trip that README:491 explicitly says does not invert: an identity `space_fn` yields a non-zero warp field of ~0.17 px at W=288, and because feedback is recursive this fixed bias is re-applied every frame. Fix: difference in one chart (`project(position)` as the origin) so the approximation cancels to first order. Costs one extra `project()` per lattice cell on the populate path only; zero on the cached path.

14. ✅ **`Plot::Ring::draw` always emits W+1 control points regardless of ring radius** — `core/render/plot.h:1852-1917`. The LUT-optimised sampler is hard-wired to `W` samples to match `TrigLUT<W,H>`, and `draw()` only calls that overload. A ring of angular radius `r` covers ≈ `W·sin(r)` pixels of circumference, so a ring at r = 0.1 rad emits 289 fragments — each running the fragment shader and a full 4-tap splat — onto ~29 distinct pixels. Fix: dispatch to the existing runtime `sample(...)` overload when `n < W`. Strictly fewer samples.

15. ✅ **`GenerativePalette` animated transitions round-trip every key through 8-bit sRGB** — `core/color/color.h:1606-1670`. The three colour keys are stored as `CPixel`, so `lerp()` converts each interpolated OKLCH key back down to 8 bits per channel every frame and `update_stops()` re-decodes it. A slow fade advances in ~1/255 sRGB steps and the stop set is quantised before the OKLCH interpolation runs — directly against the project's 16-bit-linear thesis. Fix: hold the keys as `OKLCH` or `Pixel` (36 bytes instead of 9 per snapshot). Slightly *cheaper* at runtime; `lerp` is already `HS_COLD_MEMBER`.

16. ✅ **`MAX_DEGREE` is unpinned; a table edit overflows a stack array** — `core/mesh/conway_graph.h:284,326-332,450-451`. `edges_from(node, out)` writes `out[n++]` with no capacity parameter; its only caller declares `uint8_t cand[MAX_DEGREE]` with `MAX_DEGREE = 5`. The `node_degree` static_asserts cover only 13 of 18 nodes, and nothing asserts no node exceeds the cap. Adding one edge to any degree-5 node silently smashes the stack. Fix: a `constexpr max_node_degree()` fold plus a `static_assert`.

17. ✅ **`build_mesh_class_bake` spends 8 KB of scratch unconditionally** — `core/mesh/mesh_classes.h:329-330`. `staging` is allocated before the eligibility loop and regardless of `n_elig`, which is frequently zero (a convex-only mesh). The default scratch partition is 16 KB, so a caller that has not repartitioned loses half its scratch to an unused buffer — and `Arena::allocate` traps on overflow, converting a benign no-op call into a hard trap. Fix: hoist the allocation inside `if (n_elig > 0)`.

18. **Composition polarity is enforced only by prose and one host test** — `core/mesh/conway.h:392-395`. Geometry operators take two structurally identical `Arena&` parameters whose roles invert between primitive and even-composed operators. There is no `static_assert`, no returned tag, no `HS_CHECK`. The polarity itself was verified **true for all 15 operators**, and `tests/test_conway.h:800-846` pins the operators by arena address range — but nothing pins *call sites*, and `recipe.h:78,159` applies composed operators outside `SolidBuilder`, an unshielded second consumer neither reviewer covered. Fix: make the polarity test table-driven so a new composed operator requires a row; audit `recipe.h`'s replay path.

19. ✅ **The palette-determining topology hash quantises interior angles to whole degrees with no margin guard** — `core/mesh/mesh.h:584-600`. A face whose true angle sits within an ULP of an `X.5°` boundary rounds differently under any float difference between host and device — and this region compiles `-O3 -ffast-math` on device versus a no-op on host. The hash *is* the topology id, which drives palette assignment, and the Phantasm rig requires every board to produce a bit-identical mesh. The codebase already treats host/device divergence as real enough to warrant the `RelaxBake` mechanism. Fix: a host test sweeping the registry for `min |ang − (n+0.5)|` with a margin floor.

20. ✅ **`medial` duplicates `ambo`'s body and the "bit-identical" claim is only partially pinned** — `core/mesh/conway.h:609-695`. `medial` re-implements `ambo`'s midpoint loop, shrunk-face reconstruction and orbit emission line for line. The doc states `out_a` is bit-identical to `ambo(mesh)` and the dual-morph bridge's seam depends on it, but the test pins only connectivity exactly; vertices are compared at 1e-4 *after* snorm16 quantisation. Fix: a templated shared emitter, or a host test asserting exact float equality before quantisation.

21. ✅ **No NVIC priority is set anywhere in `hardware/`, so the sync-edge ISR cannot preempt the flywheel ISR** — `hardware/pov_segmented.h:265,484`. The edge ISR exists to capture `ARM_DWT_CYCCNT` *at the edge*; at equal priority an edge arriving during a column render is latched and serviced up to the full ~96 µs ISR body later (≈0.22 column). `snap()` then rebases on that service-time stamp, injecting variable sub-column phase noise the AC timing table does not budget. The `__disable_irq()` bracket guarding the mailbox is currently decorative — raising the priority is what makes it load-bearing. Fix: one `NVIC_SET_PRIORITY` after `attachInterrupt`. Adds <1 µs to the column ISR worst case, on the ≤5 edges per revolution.

22. **A column dropped on DMA overrun is never retried** — `hardware/pov_sync.h:1306-1314`, `hardware/pov_segmented.h:626-629`. `last_rendered_x_` advances when `tick()` *decides* to render, and the shell discards `submitFrame`'s `[[nodiscard]]` result. One overrun therefore costs a whole column (a dark radial line), not a 54 µs retry, because the remaining ~7 wakes see `x == last_rendered_x_`. `render_black()` gates its latch on the return value for exactly this reason — the two paths handle the same signal inconsistently. Fix: a `resubmit_pending_` flag; the packed data is still in `backFrame()`, so no repack is needed.

23. **A stale `prepare_frame()` renders every BallDrop bump invisible** — `core/engine/transformers.h:244-266`, `core/animation/params.h:790-858`. `BallDrop::step()` writes `p.envelope` every frame but never calls `sync()`, so `cos_radius` is refreshed only by `prepare_frame()`. The constructor leaves it at `cosf(0) = 1.0`, and `bump_field`'s reject is `cos_d <= cos_radius`, true for every sample on the sphere. A dropped or mis-ordered call makes every falling ball evaluate to exactly zero — balls vanish, the frame still renders, nothing traps. All three production call sites are currently correct. Fix: a debug-only staleness stamp, scoped to pools whose params expose `sync()`.

24. **A pinned one-shot timer is destroyed silently — the pin guard's one reachable hole** — `core/animation/timers.h:57,118`, guard at `core/animation/timeline.h:305`. `RandomTimer`/`PeriodicTimer` use `cancel()` as their *natural* one-shot termination, but `cancel()` is exactly the case the pin-completion guard exempts. A pinned one-shot timer reaches the destroy branch, passes the guard, and dangles the caller's pointer with no diagnostic — while the equivalent `Transition` case traps and is death-tested. `add_get(pin=true)` checks only *predecessors'* finiteness, never that the pinned animation is itself infinite, despite the doc line above saying it must be. Fix: end via `duration` instead of `cancel()` in the one-shot branch.

25. **`Timeline::clear()` destroys pinned events and rewinds the global clock, with neither guard the other teardown paths have** — `core/animation/timeline.h:141-147`. `move_into` traps on `handled`; `step()`'s destroy branch traps; `clear()` does neither. It also resets `global_timeline_t`, desynchronising any effect deriving phase from `timeline.frame()` (a live consumer exists). It is a public runtime API — `ShapeShifter.h:98` calls it. Fix: split into a private `reset_storage()` for ctor/dtor and a public `clear()` that `HS_CHECK`s `!handled`.

26. **`core/animation/mesh.h` is 78% one class** — `core/animation/mesh.h:31-2182`. `OpLeg` spans 2152 of 2750 lines; the remaining 568 are two well-factored units. `animation.h:250` already establishes the fragment convention, and `OpLeg` cleanly satisfies it. This is a grab-bag by size, not by intent — the intent line is already drawn. Fix: move to `core/animation/opleg.h` with the same internal guard. Mechanical.

27. **GS hot paths carry no `HS_O3` annotation while its sibling's do, despite 8× the substep count** — `effects/GSReactionDiffusion.h:360,394`. GS runs 16 substeps over ~737k neighbour gathers per frame plus a per-SSAA-subsample interpolation, with no optimisation marker anywhere in its 534 lines; BZ, running 2 substeps, marks both equivalents. **Do not land on inspection** — `HS_O3` regions consume ITCM and the device is granule-cliff sensitive. If GS was excluded to pay an ITCM bill, that decision is undocumented; record it at the call site instead.

28. **`update_hankin` is marked `cold` while running every frame inside an `HS_O3` region** — `core/mesh/hankin.h:316`, region opened at `:99`. `HS_COLD_MEMBER` expands to `__attribute__((cold, noinline, noclone))`, which optimises **for size** — defeating the enclosing `HS_O3` pragma — and `tools/phantasm.ld` routes `.text.unlikely*` to FLASH. But `HankinSolids.h:613` calls it every frame inside a profile scope. `platform.h:1043` states the annotation's own contract: cold paths only. Fix: drop the attribute here (keep it on `compile_hankin`), then re-measure. Size-gate the device build.

29. **`KDTree`'s constructor and `build` are marked `cold` but run per frame** — `core/mesh/spatial.h:78,170`. `Voronoi.h:89-95` builds the tree once per frame inside `HS_PROFILE(vo_kdtree)` with a comment saying so. Same optimise-for-size + FLASH-routing consequence as item 28. Fix: drop both attributes and measure.

30. **Duplicated OpLeg-shading mesh draw across the two Conway effects** — `effects/HankinSolids.h:485-515` and `effects/IslamicStars.h:448-482`. Same operation, same inputs (`MeshState&`, `OpLeg::Shading&`): guard, profile, camera-transform, rasterize with a `shading.ramps[shading.face_ramp[fi]]` per-face lookup. The shading half is the consumer contract of `OpLeg::Shading` and belongs next to it. Fix: an inline `FacePaletteShader` adaptor in the engine. **Must stay header-inline and `HS_O3_FN`-inlinable**, preserving IslamicStars' per-face `select_face` hoisting; do not unify the two rasterizer entry points, which deliberately scan from different arenas.

31. **Class → shuffled-palette mapping reimplemented per effect** — `effects/HankinSolids.h:198-203,815-820`, `effects/IslamicStars.h:674-683`. Both independently implement `shuffle_indices` → `wrap(topology[f], NUM_PALETTES)` → per-face id. The `wrap`-into-`NUM_PALETTES` fold appears in at least seven places across the two files. This modular aliasing is precisely the class of mapping bug already recorded for Hankin ring colouring; it should live in one place with one comment. Fix: `MeshPaletteBank::assign_by_class(...)`. All sites are once-per-shape.

32. **Four effects hardcode `full_frame` false instead of folding it from their pipeline** — `effects/Comets.h:60`, `effects/SphericalHarmonics.h:218`, `effects/RingSpin.h:58`, `effects/Raymarch.h:23`. Seven of thirteen siblings derive the flag as `.full_frame = decltype(filters)::any_crosses_segments`; these four pass `{.strobe = true}` alone despite each owning a Pipeline member. The value is correct today (all four pipelines are empty), but it is not self-maintaining: adding any `crosses_segments` filter leaves `full_frame` false and the 4-board Phantasm rig renders the effect wrong. Fix: add the fold to the four constructors — a compile-time constant fold, value-preserving. Leave `Liquid2D`/`Flyby` alone; they have no pipeline.

33. ✅ **`-Wall -Wextra -Werror` reaches 1 of 9 test executables** — `tests/CMakeLists.txt:71`. The eight others get **no `-Wall` at all**. The in-file rationale (hygiene is covered by `run_tests` compiling the same headers) is true for headers and false for the eight target-specific `.cpp` files, which no `-Wall` build ever sees. Fix: `target_compile_options(<tgt> PRIVATE -Wall -Wextra)` without `-Werror`.

34. ✅ **`scripts/generate_srgb_decode.cpp` is compiled by no build and gated by no CI job** — a repo-wide grep returns zero references. It is the generator of record for `srgb_decode_lut.h`; an API drift in `color.h` silently breaks the regeneration path, surfacing only when someone next retunes the table. (The *table* is safe — `tests/test_color.h:989` sweeps all 65536 inputs.) Fix: add as an `EXCLUDE_FROM_ALL` target built on shard 1.

35. ✅ **daydream `deploy.yml` grants `pages: write` + `id-token: write` to every job** — `daydream/.github/workflows/deploy.yml:35-38`. The block is workflow-level, so the `gate` job — which checks out a *different* repository at a SHA read from a tracked file and executes its CMake/CTest — runs holding the Pages OIDC token. Holosphere's own `docs.yml` gets this right per-job. Fix: default to `contents: read`; grant per-job.

36. ✅ **daydream `deploy.yml` uses mutable action tags while Holosphere pins every action to a SHA** — `deploy.yml:51,94,119,121,142,145,162,168`. A retagged `deploy-pages@v4` executes inside a job holding `id-token: write`. Fix: pin all five to commit SHAs.

37. ✅ **The deploy gate runs a weaker test tier than presubmit CI** — `deploy.yml:112`. `ctest` runs with neither `HS_SMOKE_FRAMES: 120` nor `HS_EFFECTS_FULL: 1`, so it executes the QUICK tier: 8-frame smoke at device resolution only. That skips the 288×144 production-resolution roster passes and never reaches the effect-lifecycle transitions — the exact defect class the deploy is shipping. Fix: add the two env vars, matching `ci.yml:133`.

38. ✅ **Enabling segmented mode blanks the sphere with no indicator for up to 20 s** — `daydream/daydream.js:747-759`. `segments.active = v` is set *before* `await warmModules()` and `create()`, so `drawFrame()` takes the segmented branch immediately while `tick()` returns at the `!ready` guard and the driver keeps clearing the pixel buffer. No spinner, no message, no stats overlay for the whole warm + N×WASM-instantiate window — which the code itself budgets at `INIT_WATCHDOG_MS = 20000`. Fix: gate on `segments.active && segments.ready`; show a "spawning N workers…" row.

39. ✅ **`setParameter`/`setAnimationsPaused` auto-rebuild a faulted pool, contradicting the documented latch policy** — `daydream/segment_controller.js:648-667`, policy at `:571-573`. `onWorkerFault` documents "No auto-restart by design". But `setParameter` is wired to lil-gui `onChange`, which fires continuously during a slider drag; with a deterministic fault each drag respawns N workers that each load a multi-MB WASM module and fault again. A test pins this as intentional, so it is a design conflict rather than an oversight. Fix: restrict auto-rebuild to `setEffect`/`setResolution`.

40. ✅ **Every recorder start failure is console-only** — `daydream/recorder.js:118,140,160,169,197,254`. Six abort paths each `console.error` and return; `this.onError` — which exists so the UI can drop its recording state — fires only from `recorder.onerror`. A user on a browser without MP4/captureStream sees the button flick back with no explanation, while the codebase already has an overlay pattern for exactly this. Fix: call `this.onError?.()` alongside each `console.error`.

41. ✅ **No focus indicator anywhere in the stylesheet** — `daydream/styles/index.css` (zero matches for `focus` across 677 lines). The sidebar implements a roving tabindex and moves focus programmatically, but nothing styles `:focus`/`:focus-visible`. Keyboard users get only the UA default ring over a fully custom dark theme, and `.effect-button` sets `border: 1px solid transparent` with `transition: all`, so focused-but-not-selected is easy to lose against `.active`. Fix: one `:focus-visible` rule covering `.effect-button`, `.sort-btn`, `.tab`, `.context-lost-reload`.

42. ✅ **`#canvas` is `role="img"` on an interactive viewport with no keyboard path** — `daydream/index.html:67`. The canvas hosts OrbitControls yet is announced as a static image and is not focusable. There is no keyboard route to orbit or zoom at all. Separately the page has no landmarks and no `h1` — the only heading is the sidebar's `h3`. Fix: `<nav>`/`<main>`, a visually-hidden `h1`, and either `tabindex="0"` with arrow-key nudges or explicit "Front / Top / Side" view buttons.

43. ✅ **`Daydream`'s per-instance API is backed by mutable class statics** — `daydream/driver.js:111-121,824-832`. `W`, `H`, `DOT_SIZE`, `pixels` are `static` but written from instance methods, so a class that presents as instantiable and disposable is a hard singleton; a second instance would silently corrupt the first. `dispose()` also leaves `Daydream.pixels` dangling. Fix: move to instance fields — `pixelToSpherical`'s `dims` parameter is already duck-typed, making this near-mechanical.

44. ✅ **The Three.js lifecycle — the highest-risk JS in the repo — is untested** — `daydream/driver.js:621-721,824-877`. The only driver coverage is `advanceFrameClock` over a stub, plus `LabelPool`. Nothing exercises the rebuild/dispose path, so the ordering invariant that makes the rebuild leak-free (null `instanceColor.array` *before* `InstancedMesh.dispose()`) is enforced only by a comment. Fix: extract `dotDetailFor()` and `fitDistance()` as DOM-free functions with tests, and add a fake-mesh test asserting the dispose ordering.

45. ✅ **`update()` re-notifies a key that a re-entrant write restored, carrying a stale `old`** — `daydream/state.js:75-78`. The supersede check compares the *value*, not whether the key was already dispatched. Verified by execution: a subscriber doing `set('b',9); set('b',2)` during `a`'s notification produces `b` announced twice at value 2, the last with `old:0` — a transition that never happened. Fix: track dispatch by key.

46. ✅ **One unsubscribe removes every registration of the same callback** — `daydream/state.js:88-90`. `filter(l => l !== callback)` removes all matching entries, so two modules sharing a module-level handler both get unhooked when either disposes. Verified by execution. Fix: `indexOf`/`splice` with a double-invoke guard.

47. ✅ **A listener unsubscribed mid-dispatch still receives the in-flight event, contradicting its own comment** — `daydream/state.js:101-103`. The `slice()` snapshot makes *addition* safe, not removal. Verified by execution: a torn-down consumer gets one more callback after unhooking — the exact hazard `dispose()` exists to prevent. Fix: re-check membership per call; correct the comment.

48. ✅ **Clearing a tracked key leaves its stale value in the URL, which re-seeds state on reload** — `daydream/state.js:293-298,250-253`. `flush()` skips tracked keys whose value is nullish but never *deletes* the param. The ad-hoc writer `setParam(k, null)` correctly deletes — the two writers disagree on what "no value" means. Verified by execution. Fix: `params.delete(key)` on the tracked path.

49. ✅ **`scrollArrowState`'s deadzone goes negative for sub-1px overflow** — `daydream/sidebar_logic.js:89`. `Math.min(4, (maxScroll-1)/2)` is negative when `maxScroll < 1`, so the left arrow appears while the list is pinned to the start. Fractional `scrollWidth` is routine under browser zoom and fractional DPI; every existing test uses `maxScroll ≥ 4`. Verified by execution. Fix: `Math.max(0, ...)`.

50. **`mobius_transforms.js` has no parity pin to the C++ engine at all** — `daydream/tools/mobius_transforms.js:221-227`. The WASM parity surface pins sRGB, OKLab, HSV, procedural palette and lissajous; there is no mobius export. Both the 8-float `MobiusParams` ordering and the stereographic convention are unpinned. They currently agree (verified by hand). Fix: export `mobius_transform(...)` and add a parity test.

51. ✅ **GenerativePalette's profile/harmony constants are pinned only to literals restated in the test** — `daydream/tools/palette_math.js:229-269,299-327`. The test titled "match core/color/color.h" proves only that the JS matches numbers copied into the test file. The WASM bridge receives already-resolved h/s/v, so the engine never validates this derivation. Fix: retitle honestly, or export a resolver and diff it under a seeded RNG.

52. ✅ **Zero test coverage of the four tool pages' inline app code** — `daydream/tools/*.html` (~3470 lines). No test imports any HTML page. Substantial pure logic is stranded inline: the validator/queue machinery, `uniqueEdges`, the geodesic tessellation budget, `zoomPalette`'s frequency/phase algebra. Fix: extract the pure functions into the existing `solid_codegen.js` / `palette_math.js` modules; leave the DOM/WASM wiring inline.

53. ❌ **`daydream.js` — the app orchestrator — has no test file** — 35 KB owning the WASM lifecycle, GUI build-out, and the resolution-preset reset path. Every other root module has one. Its collaborators are all tested individually, so the gap is specifically the wiring. Fix: a test covering GUI-rebuild-on-`setEffect` and the preset-switch reset against the existing `FakeEngine`. — rejected: daydream.js has no exports and builds the DOM/WebGL/GUI stack at import time, so a headless test would need a full fake-DOM harness; the orchestrated logic is covered by its collaborators' tests.

### LOW

54. A pipeline carrying both 2-D and 3-D history silently leaves one domain unflushed — `core/render/filter.h:478-502`. Latent; no shipped effect mixes them. Fix: `static_assert(!(any_2d_history && any_3d_history))`.
55. `edge_fits_one_dot` is computed then discarded on every planar segment — `core/render/plot.h:1418-1422`. Dead by construction whenever a planar basis is in force.
56. `aux` / `Fragment::v3` is a dead register on the entire SDF path — all 17 `DistanceResult` constructions pass literal `0.0f`, yet it is copied per shaded pixel at `core/render/scan.h:105,639,1342`.
57. ✅ CSG per-row scratch is on the stack, contradicting `scan_region`'s own stack-pressure rationale one frame up — `core/render/sdf.h:1460-1461,1531-1532,1713-1716`.
58. `Scan::Circle` and `Scan::Point` render quintic radial gradients, not the "filled"/"solid" shapes the README and their doc comments describe — `core/render/scan.h:984-988,1031-1034`.
59. `DistortedRingStack`'s even-spacing precondition is unchecked — `core/render/scan.h:534-535`. A differently-spaced caller silently loses geometry.
60. Four leaf shapes duplicate the cap-bounds constructor block; twelve duplicate the non-template `distance()` wrapper — `core/render/sdf.h`.
61. `gamut_clip_preserve_chroma`'s early-out trusts the LUT cell minimum that `gamut_bracket_refine` explicitly refuses to trust — `core/color/color.h:874-876` vs `:813-830`. Narrow reachability; no triggering input constructed.
62. README describes the gamut clip as a "binary search"; the code is a walk-then-bisect and its comment says a plain bisection is **wrong** here — `README.md:1078` vs `core/color/color.h:746-750`.
63. Two independent 256-entry LUT-with-lerp implementations — `core/color/color.h:1340-1348` and `core/color/composition.h:1049-1058`. Extraction must be `always_inline`; verify with a device size/timing diff, not host timing.
64. No test pins drift under repeated `hue_rotate` application, though it feeds a recursive feedback loop — `tests/test_color.h:933-942`.
65. `require_closed_manifold` checks edge-manifoldness only; the orbit scaffolding assumes vertex-manifoldness, so a bowtie vertex yields a partial fan with no trap — `core/mesh/mesh.h:358-364`, `core/mesh/conway.h:191-217`.
66. The classification-hash collision test covers the pre-fold hash only and skips the islamic registry — `tests/test_mesh.h:603-645`.
67. README §7.7 omits `medial`, `relax_baked` and `transform_in_place`, and its blanket signature claim is false for `medial` — `README.md:1177-1195`.
68. `relax_baked`'s doc says "topology is copied through"; `PolyMesh::topology` is a per-face palette-class array and is **not** copied — `core/mesh/conway.h:1101-1145`. No operator propagates it, and that is nowhere stated.
69. `build_recipe` and `build_steps` duplicate a nine-case dispatch switch verbatim — `core/mesh/recipe.h:123-172,191-228`. A miss produces a silent replay divergence.
70. `is_jitterbug_edge` and the "no graph edge uses CHAMFER" claim are unpinned structural invariants in a file that pins four others by `static_assert` — `core/mesh/conway_graph.h:278-280,95-97`.
71. `expand` is the only parameterised operator that does not validate its parameter; at `t == 1` a centrally-symmetric face collapses to the origin and traps inside `normalize` with a message naming neither operator nor parameter — `core/mesh/conway.h:840-877`.
72. Three consecutive no-op `HS_O3` region boundaries — `core/mesh/conway.h:309-311,369-371`. Confirm with a device size diff; if non-zero, they were load-bearing.
73. `int order[MAX_CONGRUENCE_CLASSES]` is on the stack two lines from the comment explaining why the neighbouring buffers are not — `core/mesh/mesh_classes.h:307` vs `:194-196`.
74. `SphericalFieldLayout::longitude()` reads one sample past the ring at the negative-x seam — `core/math/spherical_field.h:278-286`. `util.h` guards this exact boundary in both `wrap` and `wrap_t`; this does not. Production reaches the bounded overload.
75. `NoiseProductParams::field_bound()` rests on an unverified third-party output-range claim (`|n1·n2| ≤ 1`) — `core/animation/params.h:907-911`. `field_bound` is contractually conservative and sizes culls. Fix: a host sweep asserting the bound.
76. `Orientation::set()`/`push()` accept a non-unit quaternion unchecked, unlike every peer seam (`make_basis`, `make_rotation`, `vector_to_pixel` all check) — `core/math/geometry.h:553-568`.
77. `SphericalFieldLayout::project()` returns x in `[−W/2, W/2)` while `sample_coordinates()` returns `[0, W)`; the asymmetry is undocumented and neither `project` nor its siblings carry doc comments — `core/math/spherical_field.h:101-105`.
78. Many tests draw multiple `hs::rand_f()` values inside one argument list, so the seeded sequence is evaluation-order dependent and explores a different point set per toolchain — `tests/test_plot_scan.h` (16 sites), `test_transformers.h`, `test_effects.h`, `test_animation.h`. **Production code is clean**; the bit-identical-across-boards contract is intact.
79. `StaticCircularBuffer::operator[]` carries an always-on `HS_CHECK` plus a runtime `% N` and is used inside render loops, so README:126's "never in the per-pixel hot loop" overstates the actual placement — `core/engine/static_circular_buffer.h:250,261`. The code itself routes around it in the hottest place (`sdf.h:143`).
80. An unbound `Fn` traps on host but silently returns zero on device, making the simulator the strict side and hardware the lenient one — `core/engine/platform.h:1381-1383`. Fix: use `hs::inplace_function` on both branches; also drops a vendored dependency.
81. `hs::H_OFFSET` is `static constexpr` at namespace scope but ODR-used from inline/template entities, including as a default template argument — `core/engine/platform.h:187,1019,1021`. Fix: `inline constexpr`.
82. `BlendFn` names two unrelated types; the nested one shadows a global `FunctionRef` alias whose lifetime contract would have made a stored-past-the-call use a live bug — `core/engine/concepts.h:224` vs `core/animation/mesh.h:115`.
83. `Flywheel::position()`'s documented fold-before-position precondition is violated on the burst path — `hardware/pov_sync.h:633-640,1249-1250,1448`. The maths is fold-invariant so behaviour is unaffected; the assert's stated invariant is false.
84. `POVSegmented::overrun_count()` is dead — `hardware/pov_segmented.h:217`.
85. Mailbox members are plain non-`volatile` scalars; ISR safety rests on `__disable_irq()`'s implicit memory clobber — correct in practice, but a toolchain property rather than a language guarantee — `hardware/pov_sync.h:518-523`. **Do not make them `volatile`.** Document the dependency.
86. `compile_hankin`'s `HE_NONE` fallbacks `continue` where the analogous rosette degeneracy traps, and would desynchronise `face_counts` from the face list if they fired — `core/mesh/hankin.h:193-194,228`.
87. `MeshState::set_view` accepts three parallel spans with no length-consistency check — `core/mesh/spatial.h:459-468`.
88. `copy_vector` loops `push_back` where `append_bulk` exists and is used for exactly this elsewhere; every instantiation is trivially copyable — `core/mesh/spatial.h:267-273`.
89. The relax-bake verification harness iterates two of three registries, so a future `relax_baked` in a Catalan generator would ship unverified — `tools/relax_bake_harness.cpp:26-31`. Fix: iterate `all_registries()`.
90. Registry counts are documented as README literals nothing pins (13 and 24 appear in no assert or test), and the builder's polarity note omits that the seed is built into arena `a` — `README.md:1211-1217`, `core/mesh/solids.h:199-203`. **Note:** the seed-in-`a` and composed-even-swap behaviours were independently investigated as a suspected dangling-mesh defect and found *memory-safe* — every operator binds output into `target` before opening its `ScratchScope`, and a bump arena never rewinds below an existing allocation. The real cost is peak arena, already gated at the real budget by `tests/test_solids.h:469-476`.
91. Timeline overflow soft-drops, silently and permanently terminating a self-rearming `.then()` chain — `core/animation/timeline.h:190-193,156-160`. The reject policy is deliberate and tested; only the invisibility of the consequence for chained callers is the gap.
92. `step()`'s orientation-collapse pre-pass is O(n²) with a 64-entry stack array — `core/animation/timeline.h:239-262`. Correct, `HS_COLD_MEMBER`, far below budget at real event counts. **No change recommended**; flagged so the quadratic is a known quantity.
93. A three-line sprite-scheduling body is copied across four segue schedules; `Dissolve::schedule` is byte-for-byte `Crossfade::schedule` — `core/animation/mesh.h:2240,2338,2410,2583`. The *only* meaningful duplication in an otherwise exemplary policy family.
94. Four dead public methods with full doc blocks, one carrying a six-line semantic caveat nothing can hit — `core/animation/params.h:65,119,204`, `core/animation/sprites.h:69`.
95. README §7.3 omits `BallDrop` and `NoiseProduct` from all three animation tables, though both are in the device size audit — `README.md:827,838-859,896-906`.
96. GS's stated stability bound covers diffusion only; the registered slider box exceeds it — `effects/GSReactionDiffusion.h:142-147`. At the joint feed/k corner the Jacobian's A row reaches 5.1 and B row 2.4 against a bound of 2. Consequence is bounded (the per-substep clamp saturates rather than corrupts) and the extreme is tested — but only at *default* feed/k.
97. `mean_db` is a per-frame net delta, not the per-substep mean its doc claims — `effects/GSReactionDiffusion.h:319,441-447`. A period-2 substep oscillation reads as *stalled*, arming the stabilisation detector on a field that is oscillating hard. Behaviour is arguably desirable; fix the doc.
98. DreamBalls' persistent-footprint `static_assert` excludes its runtime-sized tenant (the four presets' baked geometry) — `effects/DreamBalls.h:163-169`. The guard passes regardless of preset choice.
99. `ReactionDiffusionBase` exposes two subclass-private sampling entry points and the Wendland weight appears three times — `effects/ReactionDiffusionBase.h:116-174,316-326`. **Do not unify**: BZ's copy is the fused 4×SSAA stencil-reuse loop and routing it through the base reintroduces the per-subsample walk.
100. BZ pins its palette seed; GS draws from the live global RNG stream, so GS's colours depend on session effect-switch history — `effects/GSReactionDiffusion.h:158-161` vs `effects/BZReactionDiffusion.h:494-496`. Does not break board parity. Fix: a comment or a seed.
101. `solid_reach` and `solid_shift` are filled with the identical expression and read in the same loop — `effects/DisplacementField.h:221-222`. 224 bytes of persistent arena and a redundant call per ball per frame.
102. MeshFeedback's persistent-budget `static_assert` omits the compiled mesh, so its message over-claims what it gates — `effects/MeshFeedback.h:229-233`.
103. `device_persistent_budget_`'s member initialiser duplicates `init()`'s split and is unreachable — `spawn_shape()` overwrites it before any read — `effects/IslamicStars.h:241-246`.
104. Twin four-way dispatch switches in ShapeShifter differing only in namespace and one trailing argument — `effects/ShapeShifter.h:198-220,234-257`. Signatures are not uniform, so the cheap fix is a shared X-macro type list, not a merge.
105. Voronoi's column band is narrower than its row band, unenforced — `effects/Voronoi.h:220-232`. The code documents the consequence itself. Latent: Voronoi declares no `Pipeline`.
106. Palette-wipe arming is duplicated and only one copy has the re-entrancy guard — `effects/MobiusGrid.h:176-184` vs `effects/Comets.h:198-211`. MobiusGrid is safe only by an unstated arithmetic relationship between its timer period (120) and `WIPE_FRAMES` (60); halving the period silently reproduces the bug Comets already fixed. The *consumption* half is already an engine helper (`composition.h:1128`); the arming half is not.
107. `radius_at`/`expired` copy-paste across RingShower and Thrusters, with near-identical explanatory prose — `effects/RingShower.h:123-132`, `effects/Thrusters.h:146-155`. Only RingShower's endpoint convention is pinned by a test.
108. The same debug toggle is registered under two names — `"Show Bounding"` in `effects/RingSpin.h:71` vs `"Debug BB"` in four siblings. The GUI is auto-generated from these strings, so one control appears under two labels. Renaming touches `RingSpin.h:71,151,191-192`, `README.md:1754`, and possibly daydream favorites.
109. `params` member placement drifts from the "params last" convention in five effects — `SphericalHarmonics.h:337-340`, `MobiusGrid.h:341-345`, `Liquid2D.h:284`, `Flyby.h:232`, `Raymarch.h:177-186`. The public/private split is clean in all thirteen; only the tighter placement drifts, with no construction-order reason. Cosmetic; defer.
110. Per-sample quaternion conjugate in the harmonic field sampler — `effects/SphericalHarmonics.h:198`. `distance()` recomputes `orientation.conjugate()` for every pixel (~41k/frame at 288×144). The compiler may hoist it; that was not verified. Fix is a two-line constructor hoist, but treat the *gain* as unproven — measure, do not score from the diff.
111. No `clang-format` check in either workflow — formatting is enforced only by an opt-out local hook that every CI job explicitly disables.
112. Six of the ten heaviest CI jobs have no `timeout-minutes`; the discipline exists elsewhere in the same files and was applied inconsistently.
113. `SKIP_RETURN_CODE 77` turns a missing `python3`/`numpy` into a green skip for three generated-artifact gates — `tests/CMakeLists.txt:142,152,163`. Mitigated by standalone provenance jobs, but the ctest-side signal vanishes invisibly.
114. `lut-provenance` compares numeric tokens only, so array names, `PROGMEM`, includes and line-ending style are invisible to the diff — `.github/workflows/ci.yml:390-391`.
115. Windows SDK discovery hardcodes the `C:` drive as its fallback — `cmake/toolchain-native-clang.cmake:63-65`.
116. Five PlatformIO build hooks have no tests — notably `teensy_isystem.py`, whose silent no-op would still let `teensy-warnings` pass because vendored warnings would land in the baseline rather than the ratchet.
117. No coverage instrumentation anywhere in CI — the single largest blind spot behind the coverage map.
118. A latched worker fault never terminates the pool, leaving N idle WASM heaps resident indefinitely — `daydream/segment_controller.js:566-587`. Both other teardown paths do terminate.
119. The embind engine handle is never `.delete()`ed in `disposeApp()`, whose docblock claims to leave nothing behind — `daydream/daydream.js:582,880-907`. One page-lifetime handle, not a per-frame leak.
120. Duplicated alias-divergence self-heal, one copy bypassing the shared `repointDisplayAliases()` helper that exists for exactly this — `daydream/daydream.js:590-618` vs `:91-95`.
121. `ReadyMsg.segId` / `InitFailedMsg.segId` are declared and populated but never read; `FrameMsg.segId` *is* read and range-checked, so the unread copies invite trusting worker-supplied identity where the controller deliberately does not — `daydream/worker_protocol.js:93,99`.
122. Unreachable re-allocation guard in `precomputeMatrices`, whose comment documents an impossible state — `daydream/driver.js:751-766`.
123. `updateStats()` latches its element lookups on first call, including nulls, permanently deadening a row absent at first tick — `daydream/driver.js:780-789`.
124. `renderMainView()` clears the viewport twice per frame — `autoClear` is never set false — `daydream/driver.js:584-586`.
125. `ensureOffscreen` and `ensurePinnedOffscreen` are near-identical, differing only in dimension derivation — `daydream/recorder.js:344-364,376-392`.
126. Three small doc/code divergences in README §10: six axis labels not three, sidebar Home/End undocumented, and the recorder records the offscreen canvas's stream rather than the source canvas's — `README.md:2000,2026,2080`.
127. A recording in progress silently freezes on WebGL context loss — the button still reads Stop while no frames enter the stream — `daydream/driver.js:408`, `recorder.js:287-293`.
128. `engine_contract_wasm.test.js` asserts `getRecipe` is called by `solids.html`; a repo-wide grep finds it only in the test file. The binding is dead from daydream's side, so the contract test pins an unconsumed surface and its failure message misdirects — `daydream/tests/engine_contract_wasm.test.js:123-128`.
129. `shared.test.js` exercises almost none of `initScene`, including `dispose()`'s teardown ordering — the leak-prevention contract all four tool pages depend on — `daydream/tests/shared.test.js:24-40`.
130. `parseFloat` accepts trailing garbage in URL coercion (`?count=42abc` → 42; `?count=0x10` → 0) — `daydream/state.js:157-160`.
131. An unseeded tracked key bypasses coercion entirely, so `?flag=false` becomes the truthy string `"false"` — contradicting the doc two lines above — `daydream/state.js:154-168`.
132. `roundUrlNumber` collapses any magnitude below 5e-5 to `0` rather than dropping it, unlike the deliberate NaN handling in the same file — `daydream/state.js:16`.
133. `prettify` emits a 24-character label for large magnitudes and has no 2π case despite advertising "multiples of π" — `daydream/label_format.js:40`.
134. Constructing a second `URLSync` orphans the first, which keeps its subscription and can still arm a `replaceState` timer while being unreachable via the only handle teardown uses — `daydream/state.js:179`.
135. Integer-valued sliders show a fractional initial readout that self-corrects only after the first drag — `daydream/tools/lissajous.html:292`, `tools/slider.js:100`.
136. `drawWaveGraph` computes each palette sample three times and discards two channels each time — ~18k redundant `pow` calls per redraw, on a path triggered by every hue-slider tick — `daydream/tools/palettes.html:1093-1111`.
137. `generate-importmap.mjs` walks the gitignored `vendor/` directory — `daydream/scripts/generate-importmap.mjs:97`. Local-only; would surface as a `git diff --exit-code` failure in CI.
138. `CATALAN_BASES` and the Platonic list are engine mirrors hard-coded inside HTML, unpinnable by any test — `daydream/tools/solids.html:712-718,754`. Failure mode is a compile error, not silent.

---

## 5. Recommended sequencing

**Decide first (blocks four findings):** whether the CSG layer (item 6) is a planned public surface. Deleting it resolves 6, 7, 8 and 57 at once and shrinks the hot-path scratch.

**Land on inspection — small, safe, high value:** 2 (NDEBUG), 1 (winding check), 16 (`MAX_DEGREE`), 17 (8 KB scratch), 10/11/12 (the three missing `HS_CHECK` guards), 24/25 (timeline pin holes), 32 (`full_frame` fold), 45–49 (the `state.js` cluster, all verified by execution), 35/36/37 (deploy workflow), 3/33 (test-target flags).

**Measure on device before landing:** 27 (GS `HS_O3`), 28 (`update_hankin` cold), 29 (`KDTree` cold), 9 (ring window), 63 (LUT helper extraction), 110 (conjugate hoist). All are perf changes on a granule-cliff-sensitive target; the project's own rule is that these are validated by `pio run -e phantasm` size gates and profiles, never by host timing.

**Verify on hardware:** 21 (NVIC priority) and 22 (DMA overrun). The host simulator structurally cannot reach either — it delivers edges at their exact instant, i.e. models the *fixed* version of 21, and has no DMA transport at all. Suggested order: scope the sync pin against a column-boundary toggle to measure real edge→ISR latency; watch `overrun` in the 1 Hz health line; read `max_coast_halves` to bound the real mask window.

**Accessibility pass:** 41 and 42 together — one stylesheet rule and a small HTML restructure.

---

## 6. Coverage gaps in this review

Stated plainly, because a review's blind spots matter as much as its findings:

- **`effects/IslamicStars.h:809-1425`** — the dual-bridge / macro / reconcile scheduling block, ~40% of the file and the densest arena-lifetime region reviewed. Entry and exit were verified; the middle was not. `build_landing_`/`build_from_pal_` are raw pointers into leg-arena storage surviving resets only by a fixed-allocation-order trick — the highest-risk unexamined invariant in the codebase.
- **`core/animation/mesh.h:260-2180`** — the `OpLeg` body. Item 26 judges it structurally, not semantically. Its arena contract ("bulk state lives in an arena-allocated `Transients` that no destructor reclaims; the caller compacts between legs") is the highest-value unexamined surface after IslamicStars.
- **`core/render/sdf.h`** — `SDF::Face`'s distance path and the congruence-LUT machinery (~1200 lines), `Warp::Twist`/`WarpedVolume`, and `Scan::Volume`/`Scan::Shader`/`Scan::Mesh::draw`. The raymarching half has a different correctness regime (Lipschitz safety).
- **`core/render/filter.h:1861-2050`** — the `composite_previous_frame` inner loops, the highest arithmetic density in the file.
- **`core/color`** — the 12 coordinate modifiers and 7 of the 9 colour shades were censused but not line-audited; `GenerativePalette`'s harmony/profile construction was read only in part.
- **`core/mesh/conway_graph.h:330-631`** — walk/tour helpers and seed reconciliation.
- **`hardware/`** — `pov_segment_map.h`, `pov_single_map.h`, and `pov_segmented.h:380-440` were not read.
- **No build or test was executed.** Every finding is static reading plus call-site tracing. Items 2, 72 and the item-90 investigation explicitly recommend a device size/link-map A/B before their severity is treated as final.
- **`docs/CODE_REVIEW.md` and `docs/PROJECT_ASSESSMENT.md` were not read** at any point, per instruction. One grep incidentally surfaced a single line of the former; nothing from it informed any finding.

---

## 7. Validated non-findings

Recorded so they are not re-investigated. Each was raised as a candidate and dropped after tracing.

- **`REGISTER_EFFECT` static-initialisation order is safe.** `EffectRegistry::entries()` is a Meyers singleton; the `used, retain` anchor closes the LTO/`--gc-sections` hole; count oracles catch double or missing registration.
- **No dangling `FunctionRef` exists.** Every use outside `concepts.h` is a by-value/by-const-ref parameter or a same-call local; the one genuine storing site uses `StoredFunctionRef`, pinned by `static_assert`.
- **Composition polarity is true for all 15 Conway operators**, verified by tracing arena arguments and pinned by arena-address-range assertions.
- **The `ConwayGraph` edge table is verified against the real operators** — all 23 rows replayed and compared against the registry generator.
- **`classify_faces_by_topology` is deterministic** despite an unstable `std::sort`: ids depend only on hash equality, and `neighbor_acc` is an order-independent sum.
- **`SolidBuilder`'s arena handling is memory-safe in both polarities** (item 90). The suspected dangling-mesh defect does not reproduce.
- **The arena block is 298 KiB on both Teensy and WASM.** The 8 MiB figure is `HS_TEST_BUILD`-only; the README is correct.
- **The Canvas clear covering only the display band is not a defect** — any filter that reads outside its band sets `full_frame`, forcing the whole-buffer clear.
- **`SmoothUnion`'s bounds padding is conservative** — surface displacement ≤ k/6 against a full-k pad on both axes.
- **The worker-pool generation fence is airtight**; no stale-generation composite path exists, and there is no double-transfer or read-after-neuter.
- **Both `*_wasm.test.js` files genuinely run in CI** — the WASM trio is committed and there is no skip path.
- **All five generated artifacts are CI-gated**; `relax_bakes_generated.h` uses the strongest form (re-derivation, not diff).
- **No test target is built but never run** — every `add_executable` is either registered with `add_test` or documented as a compile-only bit-rot guard.
- **The Three.js resolution rebuild is leak-free** — correct dispose ordering, with the material deliberately built once and reused.
- **Geometry LUT initialisation ordering is clean** — every reader is behind a lazy guard or reached only from a guarded caller.
- **No retained `spawn()` result exists in production**; every retained handle uses `spawn_pinned` and `HS_CHECK`s it.
- **`HS_CHECK` genuinely survives `NDEBUG`** (distinct from item 2, which concerns `assert`).
- **`Presets`/`Animation::Lerp` cannot overlap on a shared subject** — all four call sites checked, with timing margins shown.
- **`hs::rand_int` is half-open**, so `hs::rand_int(0, RD_N)` cannot index out of bounds.
- **MindSplatter velocity cannot grow unbounded** at `friction = 1.0` — impulses are distance-bounded and surface advance is capped.
- **`Dynamo::drag()`'s unbounded `while` terminates** — a stationary leader cannot open slack.
- **AntiAlias margin spill across segment bands is not a defect** — each board renders the whole scene and clips only its output.
- **`pending_landing_` survives `std::move`** — `landing()` returns arena-backed storage.
- **`Scan::Shader::draw` does honour the clip band**; what it deliberately skips is alpha compositing and the plot-time filter stages, documented at `scan.h:1508`. **Raymarch does not bypass the pipeline** — it passes an empty pipeline to `Scan::Volume::draw` and gets normal clip handling and compositing. Only Liquid2D and Flyby bypass it, both full-screen-first with nothing beneath.
- **README §9 matches the registered parameters** in name and order for all 13 shader/curve effects.
- **`.persist` is correctly absent** from all 13 shader/curve effects — none reads back canvas pixels; every trail effect redraws from its own retained state.
- **No missing `std::launder` beyond the one logged** — `RingSpin.h:121` is the only occurrence; the other placement-news are into arena bytes via `allocate_n`, which explicitly does not construct.
