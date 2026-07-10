# Holosphere & daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine + firmware and the daydream web
simulator, reviewed together as one product. Out of scope by request:
`effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`,
vendored third-party code (`core/vendor/*`, `three.js/`, `node_modules/`), and
the generated `holosphere_wasm.{js,wasm}` artifacts.

**Method:** twenty-four component scopes were reviewed independently against the
README architecture, then **every** candidate finding was handed to a fresh,
adversarial verifier that re-derived the behavior from source before the finding
was admitted. Of 22 candidate findings, 18 survived verification and 4 were
rejected (a `TOLERANCE`-vs-`math::TOLERANCE` naming nit, a linker-script
comment-count quibble, a fail-fast-clean Hankin guard placement, and one more).
Findings are held to the project's own eligibility bar: concrete, cited,
one-commit-actionable, performance-aware, and not one of the design decisions the
project has already adjudicated (double-buffer sizing, always-on `wrap`/`HS_CHECK`
traps, the slerp antipodal gate, arena-only allocation, etc.).

---

## Overall Grade: **A**

This is, across both repositories, professional-grade work that would pass review
at a top-tier shop — and in several subsystems (the color pipeline, the numeric
math core, the filter pipeline, the segmented-POV worker orchestration, the
hardware sync layer, and the test suite) it is materially better than that. The
engineering is unusually disciplined for a spinning-LED art object: compile-time
resolution specialization, a single-block arena allocator with fail-fast traps,
16-bit linear color with an OKLCH palette path, SDF spherical rasterizers,
Conway/Hankin polyhedral mesh operators, and a WebAssembly build of the *actual*
production engine driving a Three.js twin of the physical sphere. Recurring
strengths are single-source-of-truth macros that make whole classes of drift into
compile errors, rationale comments that document *why* rather than *what*, and a
test culture that is actively anti-vacuous (non-vacuity floors, independent
oracles, a subprocess death-harness that verifies the traps actually fire).

The defect surface is shallow. There is exactly **one functional bug** on a live
user path (the Conway `meta` operator), one **P2** edge case (a dev-tool export
path), and the remaining sixteen findings are P3 — chiefly documentation drift,
narrow UX edge cases in the browser tools, and a handful of test-coverage gaps.
The clearest *systemic* weakness the review surfaced is that compositional mesh
operators and several of the heavier effects are validated structurally or via
smoke rendering only, never against a mathematical reference — which is precisely
the gap that let the `meta` bug reach CI green.

### Grades by dimension

| Dimension | Grade | One-line rationale |
|---|---|---|
| Correctness & Robustness | **A** | Degenerate inputs are trapped or given a documented fallback everywhere; the only functional bug is `meta` (F1). Held just below A+ by that bug and by robustness resting on many subtle, interlocking invariants. |
| Architecture & Design Elegance | **A** (A+ in places) | Compile-time `sdf_max_spans` CSG budgets, the variadic filter `Pipeline`, X-macro rosters, the generation-fenced worker pipeline, and the pure/peripheral split in the hardware layer are all first-rate. |
| Interface Expressiveness | **A−** | Intent- and constraint-revealing APIs (`normalized()` vs `normalized_or()`, `spawn()` vs `spawn_pinned()`, deleted rvalue borrow overloads). Docked for heavy `draw()` overload proliferation and wide template-parameter surfaces pushing complexity to call sites. |
| Readability & Style | **A−** | Dense but consistently legible; comments state load-bearing rationale, not history. Docked for raw density (`sdf.h` is ~3,600 lines) and a scattering of comment-accuracy drifts (F5, F8, F14, F16). |
| Testing | **A−** | The test tier itself is A/A+ (independent oracles, non-vacuity floors, a verified death harness, framebuffer A/B). Pulled down by compositional operators tested structurally-only and several effects being smoke-only (F1's root cause; F17, F18). |
| Documentation | **A** | Among the best-documented embedded C++ reviewed — measured error bounds, precondition/postcondition, sentinel-threshold relationships. Blemished only by isolated doc-drift P3s. |
| Performance | **A** | `fast_*` transcendental approximations, packed `uqadd16`, baked palette LUTs, warp-field caching, drop-on-overrun DMA, zero-copy WASM readback. Tradeoffs are stated, not hidden. |
| Portability & Platform Discipline | **A** | `platform.h`-first NDEBUG ordering, DMAMEM vague-linkage hazard handled by explicit specialization, sim/device bit-identity via integer wrap guards, insecure-context/bfcache/detached-buffer handling on the web side. |
| Build, CI & Tooling | **A** | Pinned action SHAs and toolchain versions, self-testing size/layout gates, provenance-hashed deploy gate, no-orphan-test guards, roster anti-drift static_asserts. |

---

## Prioritized fixes

Every confirmed defect is listed below under its severity, numbered sequentially.
Each item is one-commit-actionable and cites `file:line`. Severity scale: **P1** =
real bug with user-visible impact; **P2** = latent bug / significant edge case;
**P3** = minor quality, doc, or test-coverage defect. Findings in `daydream/`
files land in the daydream repo; all others land in Holosphere.

### P1 — functional bug on a live path

1. ✅ **Conway `meta` operator computes `kis(ambo)` instead of `kis(dual(ambo))`** —
   `core/mesh/conway.h:1105`. `MeshOps::meta` is `kis(ambo(...))` with a comment
   asserting the false identity "join = ambo". Conway's join is `j = da` (dual of
   ambo), so `meta = kj = kda`, not `ka`. For a cube these differ (kis of the
   cuboctahedron vs kis of the rhombic dodecahedron — coincidentally the same
   V/E/F counts, which is why the structural-only `test_meta_cube` passes). `meta`
   is embind-bound and user-selectable in the WASM mesh editor
   (`targets/wasm/wasm.cpp`), so a user gets the wrong solid with no trap. Sibling
   ops (`needle = kd`, `zip = dk`) and `SolidBuilder::meta`'s own "kis composed
   with join" doc confirm the intent. Fix: `meta = kis(dual(ambo(...)))`, correct
   the comment, and — since three primitive ops flip the output-arena polarity —
   update the COMPOSITION POLARITY note and the `test_composition_polarity` `meta`
   case in the same commit. Add a definitional (not merely structural) assertion
   so the regression can't recur.

### P2 — latent edge case

2. ✅ **Non-integer relax iteration count is accepted, then breaks C++ export with a
   silent copy failure** — `daydream/tools/solids.html:1447`. `updateOpParam()`
   clamps a param to `[min,max]` but never enforces integrality, even for relax's
   count-typed `iter`. Typing `2.5` stores `iter=2.5` in `state.ops`; on export
   `requireCount()` (`tools/solid_codegen.js`) throws, and `copyCode()` has no
   try/catch, so both saved-card "C++" buttons silently do nothing (console error
   only), while `saveSolid()`'s collision-check try/catch swallows the same throw
   and persists the bad value into the live preview. Fix: in `updateOpParam`, round
   count params to an integer (e.g. `val = Math.round(val)` when `def.step` is
   integral) before storing.

### P3 — minor quality, documentation, and test-coverage defects

3. ✅ **A duration-1 Sprite with any fade-out renders fully invisible on its only
   frame** — `core/animation/sprites.h:54`. `step()` increments `t` before drawing,
   so a duration-1 sprite's sole frame is `t==1==duration`; the constructor routes
   all of `duration` into `fade_out` (fade-in floored to 0, and the floor-to-1
   rescue is gated on `duration>=2`), yielding `progress==1`, `fade_out==0`, so it
   draws at opacity 0. Unreachable via the segue schedulers (they clamp
   `fade=min(window, duration/2)=0` at duration 1), so impact is confined to direct
   callers. Fix: clamp `fade_out_duration` to at most `duration-1`, or special-case
   `duration==1` to skip fade-out.

4. ✅ **Perpetual-type enumeration in the `AnimationBase::t` comment omits
   ParticleSystem** — `core/animation/animation.h:250`. The comment lists the
   never-rewound (`duration==-1`) animations whose `t` increments unbounded, but
   omits `ParticleSystem` (constructed `(-1, false)`, `sprites.h:190-193`) and an
   indefinite `Sprite`. A reader auditing the 2^32 wrap could wrongly conclude they
   are bounded. Fix: add both to the enumeration.

5. ✅ **`ripple_transform` comment understates the wavelet's normalized range** —
   `core/engine/transformers.h:344`. Comment says "−2 to 2 range covers the whole
   wavelet," but `RippleParams::prepare_thresholds` sets the accept band to
   `phase ± half_width()*2`, and `t = (dist_from_peak/half_width)*2`, so within the
   band `t` reaches `±4`; the Ricker wavelet is still substantial at `t=2`. Behavior
   is correct (the wider band is a deliberate margin); the comment misleads. Fix:
   state the band reaches `|t|=4` (or that −2..2 is only the primary lobe).

6. ✅ **`static_assert` message prints the hardware-ID decode without its active-low
   inversion** — `hardware/pov_segmented.h:93`. The message documents the decode as
   `(raw & (N-1))`, but `read_id()` computes `(~raw0) & (N-1)` (line 417) because
   the straps are `INPUT_PULLUP`/active-low — contradicting the adjacent
   "all-floating pull-ups => ID 0" note and the README. This is the same wrong
   formula a prior fix corrected in README/CODE_REVIEW but left stale in the header.
   Fix: change the message to `(~raw) & (N-1)`.

7. ✅ **Dead write-only global in the stack-measurement harness** —
   `tests/stack_measure.cpp:35`. `volatile uint8_t *g_hi` is written in `paint()`
   (line 56) and reset in `measure()` (line 78) but never read; the deepest-reach
   scan (line 88) bounds itself with `top`, not `g_hi`. Its partner `g_lo` is
   genuinely read. Per the project's remove-dead-code convention, delete `g_hi`
   (lines 35, 55-56, 78, and the `[g_lo, g_hi)` comment reference).

8. ✅ **Misplaced doc block leaves `oklab_to_linear_rgb` undocumented** —
   `core/color/color.h:610`. The doxygen block at 610-616 documents the
   3-out-param `oklab_to_linear_rgb`, but is immediately followed by the block for
   `oklab_to_lms_cbrt`, so Doxygen attaches it to the wrong function and the actual
   `oklab_to_linear_rgb` at line 648 is undocumented. Doc drift from extracting a
   helper between the comment and its function. Fix: move the 610-616 block to
   immediately above line 648.

9. ✅ **`check_spawn` comment misstates the ring's first drawn position** —
   `effects/PetalFlow.h:190`. The comment claims the first drawn position is
   `START_RHO + gap`; the ring actually spawns at
   `START_RHO + gap_accumulator - move` and `update_and_draw_rings` adds one move
   the same frame, landing at `START_RHO + gap_accumulator` (the residual, in
   `[0, gap)`). Steady-state spacing is correct; the comment misleads. Fix: reword
   to reference `gap_accumulator`.

10. ✅ **`test_reaction_graph.h` is not self-contained** —
    `tests/test_reaction_graph.h:20`. Uses `std::mt19937`,
    `std::uniform_real_distribution`, `std::sqrt`, and `std::printf` without
    including `<random>`, `<cmath>`, or `<cstdio>`; compiles only because
    `run_tests.cpp` pulls them in transitively first. Sibling `test_spatial.h`
    includes them explicitly for identical usage, so this violates the project's
    own self-contained-header convention. Fix: add the three includes.

11. ✅ **`docs.yml` paths filter omits `docs/doxygen-theme.cfg`, a real build input** —
    `.github/workflows/docs.yml:13`. The push and pull_request `paths:` filters
    gate `doxygen-custom.css` but not `doxygen-theme.cfg`, which is `cat`'d onto the
    Doxyfile (line 76) and sets `HTML_COLORSTYLE`, gamma, and the stylesheet list. A
    commit touching only that file skips build+deploy, so the published API site
    keeps a stale theme until an unrelated push rebuilds. Fix: add
    `- 'docs/doxygen-theme.cfg'` to both `paths:` lists.

12. ✅ **Recording settings changed before WASM load are silently dropped** —
    `daydream/daydream.js:555`. The Recording GUI controls are created
    synchronously and are interactive during the async WASM-load window, but their
    setters early-return while `host.recorder` is null, and construction at line 555
    never reads back `recSettings`. A quality/resolution/format change made in the
    first ~second is lost. Fix: after constructing `host.recorder`, push the current
    `recSettings` into it.

13. ✅ **Invalid numeric/boolean URL params are neither applied nor cleaned, so they
    persist and re-warn every load** — `daydream/gui.js:285`. In `DeepLinkGUI.add()`,
    a URL param that fails to parse as a number/boolean sets `urlApplied=false` but
    leaves `valClamped=false`, so the URL-rewrite gate never fires — unlike the enum
    path, which self-heals. A shared link like `?Effects.Speed=fast` logs a warning
    on every load and never cleans itself. Fix: set `valClamped=true` in both reject
    branches so invalid params are stripped like the enum path.

14. ✅ **Cancelling the native Save dialog still auto-downloads the recording** —
    `daydream/recorder.js:382`. `openSink()` detects the `showSaveFilePicker`
    `AbortError` only to mute a warning, not for the save decision: on cancel
    `handle` stays null, chunks buffer, and `finish()` falls through to
    `this.download(...)`, writing to the browser's default Downloads folder —
    contrary to Cancel semantics. Fix: set an `aborted` flag on `AbortError` and, in
    `finish()`, discard buffered chunks without downloading when aborted (still fall
    back to download for genuinely-unavailable/`createWritable`-failed cases).

15. ✅ **Drag-to-zoom on the Generative tab silently corrupts the hidden Procedural
    params** — `daydream/tools/palettes.html:734`. `handleDragEnd()`'s zoom branch
    calls `zoomPalette()` with no `activeTab` check; on the Generative tab it has no
    visible effect but rewrites the Procedural tab's `C_R..D_B` parameters. Switching
    back shows unexpectedly changed values. Fix: guard the zoom branch with
    `if (activeTab !== 'procedural') return;`.

16. ✅ **Lock R,G,B is bypassed by keyboard adjustment of a single channel slider** —
    `daydream/tools/palettes.html:533`. The locked relative-movement branch only
    activates when `lockedDragStartValues` is non-empty, and that map is seeded only
    by `mousedown`/`touchstart`. A user who Tabs to a locked slider and arrows it
    falls into the unlocked branch and moves only that channel despite the Lock
    checkbox. Fix: also seed `lockedDragStartValues` on `keydown`/`focus`, or apply
    the locked-group delta on the empty-map keyboard case.

17. ✅ **`export_params` readonly-skip test only covers a trailing readonly param** —
    `daydream/tests/export_params.test.js:16`. `formatExportParams` indexes
    `values[]` by the `params` index, so correctness under readonly params hinges on
    the arrays staying aligned. The sole test places the readonly param last, where a
    filtered-index-zip bug produces byte-identical output and passes. Fix: add a
    middle-readonly case (e.g. params `[A,B(readonly),C]`, values `[0.1,0.2,0.3]`)
    asserting `{ 0.1f, 0.3f }`.

18. ✅ **Paused-drain invariant is untested; the fixture clock is not actually
    "one-shot"** — `daydream/tests/driver_clock.test.js:26`. `advanceFrameClock`
    calls `this.clock.getDelta()` *before* the paused-return specifically to drain
    the clock each frame; the fixture's `getDelta` is stateless and the paused test
    only asserts `timeAccumulator===0`, so a regression moving `getDelta()` below the
    paused-return would still pass yet reintroduce the paused-span-replayed-as-backlog
    bug. Fix: count `getDelta` invocations and assert it is called exactly once while
    paused.

---

## Dimension notes

**Correctness & Robustness (A).** The prevailing pattern is a two-tier degeneracy
policy — strict `normalized()`/`inverse()`/`make_basis` trap on logic-bug inputs,
while `normalized_or()`, `cubic_fast`'s p1 fallback, and the slerp antipodal
branches soft-degrade exactly at legitimate geometric edges. Division safety,
NaN-from-sqrt, pole/antipodal singularities, unsigned wrap, and narrowing casts are
each addressed explicitly and, where load-bearing, pinned by tests. The verifiers
re-derived the riskiest paths (Complex division, `fast_expf` exponent assembly,
the congruence-class LUT bilinear fetch, the flywheel fold/snap math, the
color-correction LUT bounds) and each held. The single functional defect is F1.

**Architecture & Design Elegance (A, A+ in places).** Standouts: the compile-time
`sdf_max_spans` recurrence that sizes CSG interval buffers and turns a would-be
runtime overflow into a compile error; the variadic `Pipeline` that folds
stage-ordering rules into readable `static_assert`s; the `Feedback::Style` POD with
function-pointer transforms swapped without template parameters; the
generation-fenced one-frame worker pipeline with a shared bidirectional
`blitSegmentRect`; and the hardware layer's clean split of pure index/protocol math
from the Arduino peripheral shells.

**Interface Expressiveness (A−).** APIs make contracts part of the type: deleted
rvalue overloads and `StoredFunctionRef` turn dangling-lambda bugs into compile
errors; `spawn()`/`spawn_pinned()` encode the retained-handle contract; the const/
non-const `ParamList` friend split routes every untrusted-JS write through one
clamp/validate gate. Held at A− by `draw()` overload proliferation across the
rasterizers and wide template-parameter surfaces.

**Testing (A−).** The suite is the review's quiet hero — `perturb_determinism_globals`
makes each reset line load-bearing, roster names are compile-time pinned, the death
harness requires the exact probed illegal-instruction shape, banded-vs-full
byte-identity proves the one non-fail-safe filter trait override, and arena
high-water budgets catch device-only OOM on the host. The weakness is coverage
*shape*, not test quality: compositional Conway operators (`meta`/`needle`/`zip`)
are checked only for structural invariants (F1's root cause), and
`HopfFibration`/`MobiusGrid`/`Raymarch`/`ShapeShifter`/`SplineFlow`/
`SphericalHarmonics` rest on the frame-sum smoke roster without white-box math pins.

**Documentation (A) & the rest (A).** Documentation is a genuine differentiator —
measured `fast_*` error figures, sentinel-threshold relationships, arena-polarity
contracts, and WASM-boundary detach semantics are all spelled out at the point of
use. The doc-drift P3s (F4, F5, F6, F8, F9) are isolated slips against an otherwise
exemplary baseline. Performance, Portability, and Build/CI are uniformly strong: the
`fast_*` approximations and baked LUTs keep the POV frame budget; sim/device
bit-identity is enforced by integer wrap guards; and CI pins action SHAs, gates
Teensy image size/layout, and refuses to deploy a WASM artifact it cannot reproduce
from a clean ref.
