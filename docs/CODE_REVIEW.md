# Holosphere + daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine, effects, hardware drivers, WASM
bridge, native test suite, and build tooling, plus the daydream web simulator
(driver, workers, tools, JS tests). Out of scope by request:
`core/engine/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`,
`core/math/rotate.h`, all vendored third-party code (`core/vendor/FastNoiseLite*`,
`daydream/three.js/`, `node_modules/`), and the KiCad/PCB artifacts under
`hardware/phantasm/`.

**Method.** The codebase was partitioned into 18 components, each reviewed
independently against the full quality-dimension set. **Every** candidate finding
was then handed to a *fresh, independent verifier* that read the cited code
directly and defaulted to rejection when it could not substantiate the defect
from source. **81 raw findings were reduced to 35** that survived validation and
meet the fix-eligibility bar (a real, actionable defect with a correct,
proportionate fix). The rejection rate (~57%) is itself a signal: most
speculative issues did not hold up against the actual code, and two whole
components (`core/render`, `core/animation`) had *every* candidate finding
refuted on inspection.

---

## Overall grade: **A−**

This is exceptional engineering. It reads like a production embedded-graphics
engine that happens to be a personal art project: a compile-time-parameterized
rendering pipeline shared bit-identically between a 600 MHz Teensy and a
WebAssembly simulator, an arena memory model with debug use-after-free stamps, a
heap-free type-erased callable stack, a fail-fast invariant doctrine backed by a
death-test harness, and a one-wire multi-board frame-sync protocol documented to
datasheet precision. The single most striking trait is that **almost every
non-obvious decision is justified in a comment that names the exact failure mode
it prevents** — the documentation is genuinely load-bearing rather than
decorative.

It lands at A− rather than A/A+ for three reasons, all fixable: (1) a
CI-breaking stale header path that an otherwise-complete refactor missed, leaving
the screenshot-gallery job red on current master; (2) a small set of latent and
edge-case defects (a segmented-tiling config advertised but unimplemented, two
effects that skip palette-baking on the MCU hot path, a shader threshold that
drops the dimmest lit pixels); and (3) maintainability drag from a few very large,
dense files, some dead defensive code, and cross-effect duplication.

---

## Quality-dimension grades

| Dimension | Grade | One-line rationale |
|---|---|---|
| Documentation | **A+** | Dense, honest, load-bearing rationale comments + a 2,171-line README that is effectively a design spec; reviewers uniformly called it exceptional. |
| Memory safety | **A** | Arena allocator with generation/rebind UAF stamps, subtractive wrap-proof bounds math, exact-fit pool sizing with overflow `static_assert`s; no memory defect found. |
| Architectural elegance | **A** | Compile-time `<W,H>` specialization, variadic filter pipeline with compile-time domain conversion, compile-time segue policies, pure/shell hardware split, X-macro single-sources. |
| Concurrency & ISR safety | **A** | Meticulous single-writer/relaxed-atomic reasoning, cycle-counter flywheel that makes masked-IRQ windows and 32-bit wrap structurally unobservable; only issue is a JS watchdog *tuning*, not a race. |
| Correctness & reliability | **A−** | Engine core is effectively defect-free (traced hot paths held up); the few real defects are latent (N=8 tiling), edge-case (near-black cull), or in tooling (stale CI path). |
| API design & interface expressiveness | **A−** | Compile-time borrow/own contracts (`FunctionRef` vs `Fn`, `ArenaSpan` vs `ArenaVector`), fluent `SolidBuilder`; minor asymmetries and a few prose-only invariants. |
| Performance & efficiency | **A−** | LUT baking, branchless ISR, `noinline` register-spill control, measured tradeoffs; two effects evaluate a virtual palette on the per-fragment MCU path. |
| Error handling & robustness | **A−** | Consistent fail-fast `HS_CHECK` on cold seams + reject-and-return at the WASM boundary; gaps: tooling-arena OOM, `-O`-strippable generator assert, one under-tuned watchdog. |
| Portability | **A−** | Fixed-width types, host/device parity contracts guarded by compile-time `#error`; one LP64 `EVERY_N` edge and a couple of stale-path nits. |
| Testing & verification | **A−** | Death harness verifying trap shape, cross-run determinism pass, host-testable protocol cores, mock-drift/parity double-anchoring; gaps: N=8 untested, no PR trigger for the JS suite, two gates missing an empty-roster floor, screenshot CI currently red. |
| Maintainability & readability | **A−** | Mostly high; dragged by a few ~3,600-line dense files, dead defensive branches, cross-effect trail duplication, and minor doc/code drift. |
| Build system & tooling | **B+** | Excellent anti-drift and provenance design, but an incompletely-applied header-move refactor leaves a runtime path constant stale and a CI job broken. |

---

## Component scorecard

| Component | Repo | Grade | Confirmed / raw findings |
|---|---|---|---|
| core/engine | Holosphere | A | 3 / 5 |
| core/math | Holosphere | A | 3 / 4 |
| core/mesh | Holosphere | A− | 2 / 5 |
| core/color | Holosphere | A | 1 / 4 |
| core/render | Holosphere | A− | 0 / 9 |
| core/animation | Holosphere | A | 0 / 4 |
| effects (A) | Holosphere | A− | 3 / 4 |
| effects (B) | Holosphere | A− | 1 / 4 |
| effects (C) | Holosphere | A− | 0 / 5 |
| hardware | Holosphere | A− | 4 / 4 |
| targets/wasm + Phantasm | Holosphere | A− | 3 / 4 |
| tests (C++) | Holosphere | A | 1 / 3 |
| build + scripts | Holosphere | **B** | 4 / 4 |
| daydream-core | daydream | A− | 2 / 4 |
| daydream-ui | daydream | A− | 1 / 5 |
| daydream-segmented | daydream | A− | 1 / 5 |
| daydream-tools | daydream | A− | 2 / 2 |
| daydream-tests | daydream | A− | 4 / 6 |

---

## What is genuinely excellent

- **Host/device bit-parity as a first-class invariant.** FastLED math is
  reproduced bit-exactly, the NaN→hi clamp contract is protected by a compile-time
  `#error` against `-ffinite-math-only`, index widths are pinned to `uint32_t` so
  pooled-struct layouts match across 32/64-bit builds, and every board reseeds the
  same PRNG so all four Teensys render identical canvases.
- **The arena memory model.** Persistent + two scratch pools in a single fixed
  330 KiB block, RAII `ScratchScope` rewind, `Persist<T>` compaction, all-explicit
  `Arena&` parameters (no hidden state), and debug generation stamps that catch both
  reset-dangles and regrow-dangles.
- **Fail-fast discipline, applied *deliberately*.** `HS_CHECK` guards cold seams
  (arena OOM, capacity, config) and survives `NDEBUG`; the per-pixel hot path uses
  a stripped `assert` backed by a cold trap at the bind site — and a death harness
  actually spawns subprocesses to confirm the traps fire with the exact
  illegal-instruction relay shape.
- **The Phantasm one-wire sync protocol.** A cycle-counter flywheel with epoch
  folding, an odd-only distance-2 symbol alphabet that degrades to *missed* (never
  *misclassified*), and a fail-dark (never fail-wrong) content layer — documented
  to AC-timing-table precision and host-unit-tested off-device.
- **Anti-drift tooling.** X-macro rosters as single sources of truth, a
  self-registration count check at startup, resolution-driven fill-fn generation,
  positional-preset field-order `static_assert`s, and a triple-pinned test-module
  roster.

---

## Prioritized fixes

Findings are numbered sequentially. Priority reflects severity and blast radius:
**P0** = breaks CI / shipped now; **P1** = real defect worth fixing soon
(latent correctness, MCU perf, robustness/coverage holes); **P2** = polish
(readability, dead code, documentation, test hardening).

### Priority 0 — Critical (fix immediately)

1. ✅ **Stale header path breaks the roster scripts and the screenshot-gallery CI job.** `scripts/effect_roster.mjs:14` reads `core/effects.h`, but commit `ec36e8eb` moved the file to `core/engine/effects.h`; every consumer (`capture_screenshots.mjs`, `check_screenshots.mjs`, `check_effect_roster.mjs`) crashes with ENOENT, so the CI screenshot-gallery job is red on current master. **Fix:** change the path constant to `join(REPO_ROOT, 'core', 'engine', 'effects.h')` (the surrounding comments are already correct).

### Priority 1 — Medium (fix soon)

2. ✅ **FlowField evaluates a virtual palette per fragment on the MCU hot path.** `effects/FlowField.h:147` calls `palette.get()` (virtual + OKLCH lerp + `fast_sinf` + colorspace convert) once per trail fragment (~600 particles × 14 trail, per AA splat). **Fix:** add a `BakedPalette` member, bake once in `init()`, sample the LUT in the shader — matching every sibling effect (~4 KB persistent arena).

3. ✅ **GnomonicStars evaluates a virtual palette per covered pixel.** `effects/GnomonicStars.h:85-88` calls `Palettes::mangoPeel.get(t)` per fragment per star (up to 2,000 stars). `mangoPeel` is `constexpr` and never mutated. **Fix:** add a `BakedPalette` member, bake `mangoPeel` in `init()`, sample the LUT in the shader.

4. ✅ **Segmented tiling advertises N≤8 but only implements N≤4.** `hardware/pov_segment_map.h:51-66` maps only `arm_seg==0` to the top strip; for N=8 (`SEGS_PER_ARM=4`) segments 1–3 all collapse onto the same reversed bottom band, triple-painting rows [108,143] and leaving [0,107] dark — yet `pov_segmented.h:90` `static_assert`s N≤8 and the README advertises an 8-segment build. Latent (only `POVSegmented<288,4,480>` ships). **Fix:** tighten the `static_assert`, README, and strap rationale to N≤4 to match the math (per the remove-dead-over-document norm); only generalize `segment_map()` if an 8-Teensy build is actually planned.

5. **N=8 tiling is compile-permitted and README-advertised but untested.** `tests/test_pov_segmented.h:198-219` exercises `check_tiling` only for N=2 and N=4, so the broken N=8 config passes CI. **Fix:** either tighten the bound to N≤4 (companion to #4) or fix `segment_map()` and add N=8 cases — the existing harness already asserts covered-once / no-double-paint / interior-ordering and would go red on today's behavior.

6. **Segment render watchdog can misclassify a slow frame as a permanent fault.** `daydream/segment_controller.js:39,541-550`: `RENDER_WATCHDOG_MS = 8 × SLOW_FRAME_MS = 496 ms` is a single deadline armed at dispatch; a heavy effect on a throttled/mobile GPU can legitimately exceed it, latching `this.faulted` permanently with a message ("never replied") that is factually wrong for a slow-but-progressing worker. **Fix:** make the deadline adaptive — re-arm whenever any worker reports a frame while `pending > 0` — or raise the bound substantially, and correct the fault message.

7. **JS unit suite has no `pull_request` CI trigger.** `daydream/.github/workflows/deploy.yml` runs `npm test` only as a *deploy gate* on push-to-master and `workflow_dispatch`; no PR trigger exists anywhere, so a PR can regress ~30 test files undetected until after merge. **Fix:** add a small standalone `pull_request` workflow (checkout, setup-node 22, `npm ci`, `npm test`) mirroring the existing job; keep `deploy.yml` push-gated.

### Priority 2 — Low (polish / hardening)

8. **`SpriteFn`/`TimerFn` doc drift.** `core/engine/concepts.h:305-307` define these as `Fn<…,16>` with a literal, but `platform.h` and `inplace_function.h` cite them as the canonical `sizeof(void*)`-sized callsite — an idiom used nowhere. **Fix:** make the sizing explicit (`Fn<…, 2 * sizeof(void*)>`) so the documented rule actually holds and 64-bit headroom is intentional, or correct the comments.

9. **`Transformer` relies on a prose-only non-relocation invariant.** `core/engine/transformers.h:134-190`: one-shot `then()` callbacks capture `this`/`idx`, justified by "never moves," but the class remains move/copy-*constructible* (only assignment is deleted by the reference member). **Fix:** `= delete` the copy and move constructors to make the guarantee structural; no current owner moves it.

10. **`EveryNMillis::ready()` diverges host vs device for >~49.7-day intervals.** `core/engine/platform.h:998-1006` compares a `uint32_t` elapsed against a 64-bit `period_` on LP64, so huge intervals never fire on host while the 32-bit device wraps and fires. Negligible in practice (largest real interval is 10 s). **Fix:** cast both operands to `uint32_t` (or narrow `period_`), documenting the ~49.7-day ceiling via the type.

11. **Ad-hoc epsilon literal in `stereo()`.** `core/math/3dmath.h:736` uses a bare `1e-12f` pole-azimuth guard despite the file's rule that all tolerances route through named `math::EPS_*`/`STEREO_*` constants (and the identically-valued `STEREO_DIV_NUM_EPS_SQ` is a *squared* threshold, not reusable). **Fix:** add a named `STEREO_AZIMUTH_EPS` (documented as a length) and use it.

12. **Ad-hoc epsilon literal in `gnomonic()`.** `core/math/3dmath.h:787` hardcodes `1e-9f` twice for the equator divisor floor. **Fix:** add a named `STEREO_EQUATOR_EPS`; do not reuse the coincidentally-equal `math::EPS_NORMAL_SQ` (unrelated meaning).

13. **Header comment references a nonexistent "back" easing variant.** `core/math/easing.h:16` mentions `cubic/back/elastic` overshoot; no `ease_*_back` exists. **Fix:** drop "back/" from the comment.

14. **Dead `face_offsets` branch in generic `clone`.** `core/mesh/mesh.h:500-519`: the `if constexpr (requires { dst.face_offsets; })` block can never instantiate — `MeshState` returns early, `PolyMesh` lacks the member — and implies a nonexistent accessor contract. **Fix:** delete only that block; keep the live topology branch.

15. **Untyped `arena.allocate` + `static_cast` idiom repeated ~18× per file.** `core/mesh/conway.h:423-428` (and across `mesh.h`, `filter.h`, `render/*`, `reaction_graph.h`) spells element type, count, `sizeof`, and `alignof` independently, so a copy-paste can silently mis-size a buffer. **Fix:** add `template<class T> T* Arena::allocate_n(size_t n)` and route the sites through it.

16. **Envelope trig inconsistency in `GenerativePalette`.** `core/color/color.h`: `key_oklch` (1219) authors with exact `sinf`, while recovery (1278) and `get()` (1358) use `fast_sinf`. Masked today by an 8-bit sRGB round-trip. **Fix:** switch `key_oklch` (cold path) to `fast_sinf` so all three envelope sites agree, or add a note.

17. **FlowField comment drift.** `effects/FlowField.h:90-108`: comment says `p.z*scale + t` but the code samples `p.position.z * params.noise_scale + t` (no `p.z` member). **Fix:** update the comment to the real expression.

18. **`MindSplatter::Params::lerp` lacks the `sizeof` guard its sibling has.** `effects/MindSplatter.h:158-166` hand-enumerates fields with no `static_assert(sizeof(Params) == N*sizeof(float))`, unlike `Liquid2D`; adding/reordering a preset float would silently leave it un-interpolated. **Fix:** add the matching `sizeof` `static_assert` (do not factor a shared helper — the divergent per-effect lerp strategy is intentional).

19. **`segment_map()` doc bound contradicts the driver/README.** `hardware/pov_segment_map.h:44` documents N≤4 while the `static_assert`/README say N≤8. **Fix:** reconcile toward the true limit (N≤4) as part of #4/#5.

20. **Duplicated safety trap across commit/join adoption paths.** `hardware/pov_segmented.h:544-557`: both branches carry a character-identical `HS_CHECK(!p->overrides_get_pixel(), …)`; a future change must be mirrored by hand. **Fix:** extract just that invariant into a tiny `assert_render_column_safe(Effect*)` helper; keep the divergent policy inline.

21. **Unreachable divergence branch in `drawFrame()`.** `targets/wasm/wasm.cpp:450-455`: the width/height mismatch guard can never fire (resolution is bound at construction; `setResolution` nulls the effect). **Fix:** delete the branch and its comment (do *not* replace with a trap — this path deliberately reject-and-returns).

22. **`getRenderUs()` returns stale timing on the (dead) divergence path.** `targets/wasm/wasm.cpp:498-500,437-455`: `render_us` is reset only on the normal path, so the divergence-blank frame reports the prior frame's time. **Fix:** hoist `render_us = 0.0` above the divergence check (moot if #21 lands).

23. **`getFaces()` marshals one embind `push()` per index.** `targets/wasm/wasm.cpp:835-851` does O(total indices) C++→JS crossings; `getVertices`/`classifyFaces` bulk-copy via `typed_memory_view`. Tooling-only. **Fix:** emit a flattened faces buffer + face-counts array as two typed-memory-view copies and unflatten in JS.

24. **Arena-budget gate passes green on an empty effect roster.** `tests/arena_measure.cpp:60-79` gates only on `g_worst_total > budget`; an emptied `HS_EFFECT_LIST` (or dropped `measure()` calls) leaves it 0 and PASS. `stack_measure.cpp` has the same gap. **Fix:** count measured effects and assert `g_measured == HS_EFFECT_COUNT` (the constexpr is already available).

25. **Stale `--assume-filename` path in the LUT generator.** `scripts/generate_luts.py:114` points clang-format at `core/color_luts.h` (moved to `core/color/`); works only incidentally via the single root `.clang-format`, and would silently diverge if a directory-scoped config were ever added. **Fix:** point it at `core/color/color_luts.h`.

26. **`linear_to_srgb` LUT (65,536 evals) computed twice per emit.** `scripts/generate_luts.py:155-160`: `check()` and `render()` each rebuild the tables; the CI provenance job pays it every push. **Fix:** build `fwd`/`rev` once in `main()` and pass them to both.

27. **Reaction-graph completeness guard is `-O`-strippable.** `scripts/generate_reaction_graph.py:91-94` enforces its "fail loudly, never ship a biased table" invariant with a bare `assert`, which `python -O` removes. **Fix:** replace with an explicit `if not (...): raise RuntimeError(...)`.

28. **Shader black-pixel cull also drops the dimmest lit pixels.** `daydream/driver.js:668-672`: the `dot < 1e-8` cull on normalized `uint16` color discards channels 1–6 (≈0.009% brightness) that the engine legitimately emits on fades. **Fix:** test for exact zero (buffer is `fill(0)`'d and normalized 0 → 0.0), e.g. `== 0.0` or a `< 1 LSB²` threshold like `1e-10`.

29. **`setPixelRatio` capped at 1 makes recordings soft with no escape hatch.** `daydream/driver.js:160,388`: unconditional `min(devicePixelRatio, 1)`; `recorder.js` captures the same 1× backing store, so exports can't recover detail. **Fix:** allow a higher ratio while recording (e.g. `min(dpr, 2)`) or a config toggle; at minimum document the intentional fill-rate tradeoff at both sites.

30. **`sortItems` name ordering uses the ambient locale.** `daydream/sidebar_logic.js:22-28`: `localeCompare` with no locale arg; a pin test exists but collation is environment-dependent. **Fix:** `localeCompare(b.name, 'en')` (avoid `sensitivity`/`numeric` options, which change semantics).

31. **`copyWithFeedback` permanently drops the idle color class when `revertText` is empty.** `daydream/tools/clipboard.js:100`: idle classes are removed unconditionally but restored only when `original` is truthy; `lissajous.html` and `palettes.html` pass `revertText: ''`, so after the first copy the span loses `text-gray-500` for good. **Fix:** drop the `&& original` condition from the restore, fix the JSDoc, and add a `revertText: ''` test (with a `classList`-tracking stub).

32. **Unreachable second clamp in `linearRgbToHex`.** `daydream/tools/color.js:45-46`: the input is already clamped to [0,1] and `linearToSrgbFloat` maps into [0,1], so the `[0,255]` clamp on the byte can never fire. **Fix:** drop the redundant clamp.

33. **URLSync's auto-flush-on-change path has no end-to-end test.** `daydream/tests/state.test.js:187-205` only tests manual `flush()` and the cancel case, never `set(trackedKey) → tick(200) → exactly one replaceState`. **Fix:** add a `mock.timers` test asserting the debounced write fires once with the new value.

34. **Boolean threshold at exactly 0.5 is untested.** `daydream/tests/param_sync.test.js:25-30` probes 0/1/0.6/0.4 but not 0.5 — the one input where `>` vs `>=` flips. **Fix:** add a case pinning `resolveParamSync(true, 0.5, …)` → `false`.

35. **Recorder toggle test's post-`onstop` assertion is trivially satisfied.** `daydream/tests/recorder.test.js:213-220`: `isRecording` is already false after synchronous `stop()`, so the assertion validates nothing about `onstop` teardown. **Fix:** assert the observable teardown effects (`mediaRecorder === null`, track stopped), mirroring the stronger test at line 243.

---

## Notes on what was *not* flagged

Two components — `core/render` and `core/animation` — had **every** candidate
finding refuted by independent verification: the suspected defects (a soft-coverage
denominator in `Volume::probe_occluder`, `MobiusFlow` non-finite reads, a
`MeshMorph` correspondence edge, several divide-by-zero worries) turned out to be
either already guarded, intentional-by-design, or provably in-range on tracing.
Several effect-level suspicions (duration-0 `Mutation`, unbounded morph
accumulation, KDTree scratch overlap, coherence-grid boundary indexing) likewise
held up as correct and, in most cases, explicitly justified in comments. This is
worth recording: the engine's hot paths are unusually robust, and the residual
findings cluster in tooling, latent/unshipped configs, and JS-side polish rather
than in the rendering core.
