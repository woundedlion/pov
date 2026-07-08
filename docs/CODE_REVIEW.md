# Holosphere / daydream — Code Quality Review

**Scope:** the Holosphere C++ engine + firmware repository and the daydream web-simulator repository, reviewed as one product.
**Method:** 19 component reviewers (covering every in-scope source file) each graded its subsystem and gathered findings; every candidate finding was then re-checked by a *fresh, independent validator* that read the cited code before the finding was admitted. Only findings that survived independent validation appear below.
**Out of scope (per instruction):** `effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, and third-party vendored code (`FastNoiseLite`, `three.js`, `lil-gui`, `node_modules`) except its integration.

---

## 1. Executive Summary

This is a **top-tier codebase — overall grade A.** Across ~90 in-scope C++ headers/units and ~40 JavaScript modules, the review confirmed **19 defects, every one low-severity.** Zero critical, high, or medium issues survived independent validation. Nine additional candidate findings were raised by reviewers and then *rejected* on validation as non-issues, already-handled cases, or intentional documented designs — a signal of how little genuine defect surface remains.

The engineering is unusually disciplined for a solo art project: a single C++ rendering core compiles unchanged to three targets (two Teensy firmwares and a WebAssembly simulator); memory is a deterministic partitioned arena rather than a fragmenting heap; invariants fail *fast and loud* (`HS_CHECK` traps that survive `NDEBUG` and are themselves verified by a death-test harness); and correctness-critical arithmetic (POV sync, segment index math, param marshaling, color LUTs) is factored into Arduino-free headers with host unit tests wired into CI. Documentation is the standout: non-obvious invariants are explained *at the point of use* with the failure mode they prevent.

The residual findings are almost entirely **latent robustness gaps, cross-file consistency nits, test-coverage additions, and documentation corrections** — the kind of backlog a codebase accumulates only after the substantive bugs are already gone.

---

## 2. Letter Grades by Quality Dimension

| Dimension | Grade | Rationale |
|---|:---:|---|
| **Correctness & Reliability** | **A** | Degenerate cases (poles, antipodes, zero-length vectors, div-by-zero, unsigned wrap) are handled deliberately with documented sentinels and fail-fast traps; fast-math approximations stay within stated error bounds. The confirmed defects are latent edge cases, not live miscomputations. |
| **Architecture & Design (elegance)** | **A** | One templated `<W,H>` core fans out to three targets with zero runtime generality cost. Clean seams: SDF shapes vs. scan driver, owned vs. borrowed containers, pure math vs. DOM/WASM/worker glue, CRTP families over a thin virtual interface. |
| **Interface / API Expressiveness** | **A** | Concept-constrained callables, a uniform `distance()/get_*_bounds()` shape contract, type-erased `FunctionRef`/`PipelineRef` to bound instantiation count, and `@ts-check` typed worker-message unions give precise, misuse-resistant boundaries. |
| **Readability & Code Style** | **A** | Consistent naming; dense but high-signal comments that explain *why*, not *what*. The only marks against are a few cross-file convention inconsistencies (below). |
| **Documentation** | **A** | The strongest dimension. Doxygen/JSDoc is complete and accurate, the 2,000-line README is a genuine architecture reference, and load-bearing invariants are documented where they live. |
| **Testing** | **A−** | Broad native + JS suites with independent oracles, determinism/perturbation passes, death tests, WASM parity/golden pins, and CI budget gates. Held just below A by a handful of coverage gaps in stateful/DOM-binding paths and one mis-set guard threshold. |
| **Performance & Efficiency** | **A** | Hot paths are LUT/baked lookups, branchless clamps/wraps, arena scratch, and measurement-driven `noinline` splits. `HS_CHECK` is confined to cold paths by explicit policy. |
| **Memory Safety & Resource Mgmt** | **A** | No heap on hot paths; arena bounds trap on overflow; dual debug generation-stamps turn dangling arena references into faults; JS workers/streams/canvases are disposed on every exit path. |
| **Error Handling** | **A** | Coherent fail-fast doctrine on the device (trap invariant violations, softly handle genuinely transient conditions) and defensive recovery in the browser (worker faults, module-load retries, context-loss). |
| **Portability / Platform Abstraction** | **A** | `platform.h` isolates every Arduino/WASM/Desktop divergence; device/host layout parity is enforced with `static_assert`s and documented; JS guards Node-vs-browser and feature-detects APIs. |
| **Maintainability** | **A** | X-macro single-source-of-truth rosters, compile-time drift traps, and factored pure helpers keep extension low-risk. The inherent tax is the density of coupled invariants a maintainer must respect. |

**Overall: A.** Weighted across all dimensions the codebase sits at the A / A− boundary and rounds up on the strength of its documentation, architectural coherence, and the near-absence of substantive defects.

---

## 3. Component Scorecard

| Component | Grade | Confirmed findings |
|---|:---:|:---:|
| `core/engine` — platform, callables, registry, buffers | A | 0 |
| `core/engine` — memory/arena, generators, presets, reaction graph | A | 0 |
| `core/math` — 3dmath, geometry, easing, waves | A | 2 |
| `core/color` — color, LUTs, palettes | **A+** | 0 |
| `core/mesh` — mesh, classes, spatial/KDTree | A | 0 |
| `core/mesh` — conway, hankin, solids | A | 1 |
| `core/render` — canvas, scan, sdf, shading | A | 0 |
| `core/render` — plot, filter, led | A | 1 |
| `core/animation` — timeline, params, motion, trails, sprites | A | 2 |
| `effects/` — group A (14 effects incl. reaction-diffusion) | A | 0 |
| `effects/` — group B (14 effects) | A | 0 |
| `hardware/` — DMA, POV single/segmented, sync core | A | 1 |
| `targets/wasm` + Phantasm.ino + scripts | A | 2 |
| `tests/` — native C++ suite + harness | A | 1 |
| `hardware/phantasm/gen` — KiCad PCB generators | A− | 2 |
| daydream — app core (driver, host, state, geometry) | A | 1 |
| daydream — workers, sidebar, recorder, importmap | A− | 1 |
| daydream — geometry/palette tools | A | 3 |
| daydream — JS test suite | A | 2 |

The `color` subsystem is called out at **A+**: a full 16-bit-linear / OKLab / OKLCH perceptual pipeline with hue-preserving gamut mapping, every float→int cast clamped with a documented NaN contract, and comprehensive CI-wired coverage — reviewed defect-free.

---

## 4. Subjective Assessment

**Architectural elegance — A.** The defining move is the compile-time-resolution templated core that is *the same code* on a 600 MHz Cortex-M7 firing microsecond LED columns and in a browser tab at 288×144. The filter pipeline as a variadic template that automatically inserts world→screen→pixel domain conversions at compile time is genuinely elegant, and the arena model — persistent + two scratch pools, RAII scope guards, explicit `Arena&` at every call site — is a principled answer to embedded heap fragmentation rather than an ad-hoc one.

**Interface expressiveness — A.** APIs consistently make illegal states hard to reach: `explicit` constructors block silent coordinate-space coercions, strict (`normalized()`, traps) vs. graceful (`normalized_or()`) variants are named and separated, the shader `Fragment` register conventions are documented per rasterizer, and the worker protocol is a single JSDoc-typed source of truth. The expressiveness cost the project accepts knowingly — heavy template instantiation and a dense web of invariants — is the correct trade for a zero-overhead, hardware-safe core.

---

## 5. Prioritized Fix List

All 19 items are low-severity. Priority reflects *type of impact*, not urgency. Numbering is sequential across the whole list.

### Priority 1 — Correctness & Robustness (latent behavior)

1. ✅ **`ColorWipe` is missing the `duration >= 0` guard its sibling tweens enforce** — `core/animation/params.h:335-339`. A perpetual (`-1`) duration is stored verbatim; `step()` then clamps `t/-1` to 0 forever, silently freezing the palette and never firing `.then()` or being removed. `Transition` and `Mutation` both `HS_CHECK(duration >= 0, ...)`. Add the same check to `ColorWipe`.

2. ✅ **`Screen::Trails::plot` drops the current-frame sample downstream for `alpha <= 0.001`, unlike every other filter** — `core/render/filter.h:952-971`. The early-return happens *before* `pass()`, so a dim live pixel is removed from the pipeline entirely rather than merely excluded from the trail buffer (the sibling `World::Trails` forwards it). Call `pass(...)` unconditionally, then gate only the history-seeding block.

3. **`getMaxBounds()` dev binding measures mesh sizes via unowned `ArenaVector` members instead of the view-aware getters** — `targets/wasm/wasm.cpp:1102-1104`. Reads `temp.faces.size()` / `temp.face_counts.size()` directly; the rest of the file uses `get_faces_size()` / `get_face_counts_size()` because a `PolyMesh` may be view-backed. Correct today (solids return owned meshes) but would silently under-size the constants this dev tool exists to compute. Switch to the view-aware accessors.

4. **Palette-tool wave graph plots linear-light values for `GenerativePalette` but sRGB cosine for `ProceduralPalette` on the same axis** — `tools/palettes.html:1017` (daydream). `getChannelValue` returns raw sRGB for one palette type and linear-light for the other, so the value band and `>1 / <0` clamp overlay no longer correspond to what the device clamps on the generative tab. Return a consistent domain (convert the generative sample back with `linearToSrgbFloat`, or convert per-tab in `drawWaveGraph`).

5. **`computeSegmentRange` leaves an unrendered trailing column on odd-width canvases** — `segment_layout.js:65` (daydream). `armW = floor(w/2)` leaves column `w-1` uncovered (permanent black) for odd `w`, with no diagnostic. Unreachable via the two shipped even-width presets (96, 288). Either reject odd `w` in validation (symmetric with the existing height check) or extend the last arm's `x1` to `w`.

6. **`prettify()` renders `NaN` / `Infinity` as literal strings on labels** — `label_format.js:39-40` (daydream). Non-finite inputs fall through every symbolic branch to `r.toFixed(3)`, surfacing `"NaN"` in the coordinate readout for a degenerate (e.g. zero-length-direction) label point. Add an early `if (!Number.isFinite(r)) return '—';`.

### Priority 2 — Portability, Consistency & Maintainability

7. ✅ **`check.py` / `shorts.py` default to a hardcoded machine-specific schematic path** — `hardware/phantasm/gen/check.py:9`, `shorts.py:6`. Both default `SCH` to `c:\work\Holosphere\...\phantasm.kicad_sch`, which only resolves on the author's checkout, unlike `pcb.py`/`fab.py` which derive from `__file__`. Compute the default from the script location.

8. ✅ **`check.py` extra-net filter references a removed component (`Q_GATE`) and a renamed node (`CLF`)** — `hardware/phantasm/gen/check.py:64`. The board now uses a Schottky `Q_REV` and an `LF_DAMP` node, so both suppression prefixes are dead and could mask a genuinely unexpected net. Drop `"CLF"` and `"Q_GATE"` from the `startswith()` list.

9. ✅ **`capture_screenshots.mjs` and `check_screenshots.mjs` resolve the screenshots directory against different roots** — `scripts/capture_screenshots.mjs:30`. Capture writes to a `process.cwd()`-relative `docs/screenshots`; the freshness gate reads a script-anchored `REPO_ROOT/docs/screenshots`. If ever run from another cwd they diverge silently. Anchor capture on `REPO_ROOT` too.

10. **`Quaternion` degeneracy guards use raw `float` epsilon instead of the module's named `EPS_*` constants** — `core/math/3dmath.h:511,549,562`. Three different un-named thresholds for the same "too small to normalize" decision (one even compares a squared magnitude against an unsquared epsilon), while the `Vector` path uses `math::EPS_NORMALIZE_SQ`. The header's own doc block asks that tolerances be reviewed in one place. Route all three through a named constant.

11. **Hankin degenerate-fallback eagerly calls the strict `p_corner.normalized()` inside `normalized_or`, defeating its graceful-degradation intent** — `core/mesh/hankin.h:316`. C++ evaluates the argument first, so a zero corner would trap in the very fallback meant to survive it; the three `p_corner`-handling styles in one function also disagree. Latent (corners come from unit-sphere vertices). Compute one non-trapping `cn = normalized_or(p_corner, p_corner)` and reuse it.

### Priority 3 — Test Coverage Gaps

12. **Death-case floor (`MIN_DEATH_CASES`) is one below the actual count, so a single dropped case passes silently** — `tests/test_death.h:1118-1119`. `MIN_DEATH_CASES = 38` but the table has 39 cases; deleting one still satisfies `>= 38`, defeating the anti-shrink guard for the common single-deletion mistake. Set it to the exact count (39) or derive it from a compile-time `ALL_CASES_COUNT`.

13. **`EngineHost.refresh()` has no direct test for the detached-view re-fetch — the module's raison d'être** — `tests/engine_host.test.js:28-41` (daydream). Null-view and live-view paths are tested; the detached (heap-grown) path is only pinned one layer down in `pixel_view.test.js`. Add a test seeding a detached view and asserting `refresh()` re-fetches and fires `onViewRefreshed`.

14. **`banner.js` has no dedicated test while every other `tools/` module does** — `tools/banner.js:24` (daydream). `require-tests.mjs` only checks the glob is non-empty, so the gap is silent. `showFatalError` has real logic (idempotent reuse, `textContent`-not-`innerHTML` injection safety, `document.body` guard). Add `tests/banner.test.js`.

15. **`sidebar_dom.test.js` exercises only Enter/Space, not the arrow/Home/End DOM focus wiring of `onKeyDown`** — `tests/sidebar_dom.test.js:38-50` (daydream). The navigation *math* is covered in `sidebar_logic.test.js`, but the DOM binding (focus the target button, `preventDefault` on arrows) is not. A right-index/wrong-element regression would pass both suites. Drive the arrow keys against a multi-button fake and assert focus + `preventDefault`.

### Priority 4 — Documentation & Style

16. **`segment_id_` comment overstates strap-bit count ("up to 3") — design supports at most 2** — `hardware/pov_segmented.h:676`. `ID_STRAPS` is `(N<=2)?1:2` with `N<=4`; the surrounding docs already say 2. Correct the isolated comment.

17. **Antipodal-detection basis differs between `Vector` slerp and `make_rotation(from,to)`** — `core/math/3dmath.h:1108,1244`. One gates on the raw dot `d < -1 + TOLERANCE`, the other on `fast_acos(d) > PI - TOLERANCE`; `fast_acos`'s ~1.3e-4 error exceeds `TOLERANCE`, so the two gates are non-equivalent. Harmless (the fall-through stays finite) but a reader expects them phrased identically. Key slerp's branch off the raw dot.

18. **Color-strip marker uses `t*width` while the gradient is sampled over `width-1`, causing a sub-pixel marker/color offset** — `tools/palettes.html:905,933` (daydream). Cosmetic drift, worst at the right edge. Use `Math.round(t * (width - 1))` to match the gradient denominator.

19. **Trailing-underscore member convention is applied inconsistently across sibling animation classes** — `core/animation/params.h:69-76`. Adjacent classes doing the same job mix `params_`/`scale_` with bare `params`/`from`/`to`. Purely cosmetic; pick one convention and apply it uniformly within the file.

---

## 6. What Was Deliberately *Not* Flagged

Nine reviewer-proposed candidates were rejected on independent validation, and are recorded here so the ledger is complete — several are worth stating because they read like bugs but are correct by construction:

- **`Dynamo` color-wipe boundary-slot pointer aliasing** (`effects/Dynamo.h`) — bounded by the `palettes.is_full()` gate and the paired push/pop `size == boundaries.size()+1` invariant; no physical ring slot is ever reused while a `Transition` reference is live. Safe by construction.
- **`Voronoi` zero-site / KD-query edge cases** — `active_site_count()` floors at 1, `updateParameter` clamps to `[1,400]` and rejects non-finite, and every `knn[1]` consumer guards on `size() > 1`.
- **`FlowField` emitter skipping `life==0` slots** — those are deterministically the drained slots `step_particle` already ignores; fresh particles get `life = max_life`. No RNG-stream divergence.
- **`HopfFibration` trail-staging `static_assert`** — `ArenaVector::bind` allocates exactly `capacity*sizeof(T)`; bookkeeping lives in the stack handle, and an edge overflow would trap rather than corrupt.
- **`Recorder` timed-fallback stream handling** (daydream) — every early-return path stops its stream's tracks; the start-time blit is a deliberate primer.
- **Sidebar `setActive` stale `activeName`, segment `x=0` seam line, worker `handleMessage` rethrow, importmap injection timing** (daydream) — each is either a documented intentional design or a graceful-degradation path with no reachable failure.

The review also explicitly did **not** treat the project's deliberate designs as defects: fail-fast `HS_CHECK` traps, the no-bounds-checking hot paths, the two-physical-buffer double-buffer, compile-time resolution, and DTCM/DMAMEM placement are all intentional and correct for the target.

---

*Every in-scope finding above is confirmed real, minimal to fix, and does not require sacrificing performance. The list is a polish backlog on an already-shipping, production-quality codebase — not a defect triage.*
