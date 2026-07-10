# Holosphere + daydream — Code Quality Review

**Scope:** the Holosphere C++ rendering engine / firmware and the daydream web
simulator, reviewed as one product across both repositories.
**Date:** 2026-07-09.
**Out of scope (per review charter):** `core/engine/effects_legacy.h`,
`targets/Holosphere/Holosphere.ino`, `core/math/rotate.h`, the vendored
`daydream/three.js/` tree, `node_modules/`, and the generated
`holosphere_wasm.{js,wasm}` artifacts.

**Method.** Every in-scope file was read by a component-specialist reviewer
(21 components spanning both repos). Each candidate finding was then handed to a
*fresh, independent validator* that re-read the cited code and either confirmed
or refuted it, with the codebase's stated design philosophy (fail-fast
`HS_CHECK` traps, arena-only render path, two-buffer ISR, compile-time `<W,H>`
templating) treated as intentional rather than as defects. The 19 findings below
are the survivors of that adversarial pass; a large number of raw candidates were
rejected as intentional-by-design or misreadings and are not listed.

---

## 1. Executive Summary

This is a top-tier codebase. Across roughly 150 in-scope files spanning a
hard-real-time embedded renderer, a novel multi-board POV sync protocol, a
perceptual color pipeline, a polyhedral-mesh geometry kernel, and a browser
simulator that shares the exact same C++ compiled to WASM, the review surfaced
**19 real defects — 2 high, 5 medium, 12 low** — and *none* of them are memory-safety
bugs, UB, data races, or leaks. The high-severity items are a single unguarded
allocation on the device effect-swap path and a build-tooling script that crashes
in its own summary block; both are contained and mechanically fixable.

The dominant signal from the review is not the defect list but its *shape*: the
findings are overwhelmingly documentation drift, dead code, and narrow
never-yet-reached edge cases. That is the finding profile of a codebase that has
already been through many rigorous audit cycles. The engineering discipline —
explicit arena threading, cold/hot fail-fast placement, compile-time enforcement
of borrow contracts, and comments that encode *proofs* of non-obvious arithmetic
identities — is consistently excellent and, in several subsystems, exemplary.

**Overall grade: A+** (with two A-level caveats: a couple of doc-accuracy drifts,
and a CI gate that is syntax-only where it claims runtime protection).

---

## 2. Letter Grades by Dimension

| # | Dimension | Grade | One-line rationale |
|---|---|:---:|---|
| 1 | Correctness & Reliability | **A** | 19 defects across ~150 files, none memory-unsafe; the two highs are contained and mechanical. |
| 2 | Memory Safety & Resource Management | **A+** | Three-arena bump allocator, no heap on the render path, RAII `ScratchScope`/`Persist<T>`, generation-tracked handles; zero leak/UAF/overflow findings. |
| 3 | Architecture & Modularity | **A+** | Pure-logic/hardware-shell splits, `<W,H>` compile-time specialization, CRTP bases, X-macro single-source lists, DOM-free JS cores. |
| 4 | Architectural Elegance *(subjective)* | **A+** | Repeatedly and independently judged "elegant"/"exemplary": Segue policy pattern, SDF/Scan/Fragment split, congruence-class LUT that degrades to the exact path. |
| 5 | Interface Expressiveness *(subjective)* | **A+** | Ownership/lifetime encoded in types (`ArenaVector` vs `ArenaSpan`, `Arena&` everywhere); borrow contracts enforced by *deleted overloads*, not comments. |
| 6 | Concurrency & Real-Time / ISR Safety | **A** | Single-writer ISR ownership, minimal IRQ-off windows, and a genuinely novel 1-wire flywheel sync protocol; two-buffer design correct-by-design. |
| 7 | Error Handling & Fail-Fast Robustness | **A** | Cold-path `HS_CHECK` traps applied precisely per the stated philosophy; the one asymmetric gap (finding 1) is the exception that proves the rule. |
| 8 | Documentation | **A+** | Load-bearing, not decorative — comments carry correctness proofs; the README is a ~2,200-line engineering document. Minor doc-accuracy drift keeps it human. |
| 9 | Test Coverage & Test Quality | **A** | Death harness, cross-run determinism pass, wire-format tests, host-testable pure cores, self-testing test infra (anti-drift module guards). |
| 10 | Performance Engineering | **A+** | Profiling-driven cold/hot split, LUT baking, branchless ISR loops, SIMD, fast-reject heuristics — with the tradeoffs *measured and commented*, not guessed. |
| 11 | Portability & Cross-Platform Parity | **A** | One engine → Teensy / WASM / desktop with bit-identical determinism; CI on Linux + Windows toolchains. |
| 12 | Build System, CI & Tooling | **A−** | CMake presets, justfile, PlatformIO size gates, 3-layer CI; docked for the broken screenshot script and its syntax-only CI gate (finding 2). |
| 13 | Security & Input Validation | **A** | Minimal attack surface (offline firmware + static web app, no untrusted input); WASM buffer-detach hazards handled defensively. |
| 14 | Naming & Style Consistency | **A** | Conventions are documented and adhered to; a single leading-underscore violation (finding 19). |

---

## 3. Subsystem Notes

**C++ engine — core (`engine`, `memory`, `math`, `color`, `render`, `mesh`,
`animation`).** The strongest part of the codebase. `platform.h` cleanly isolates
target divergence behind a uniform `hs::` namespace and even bit-matches
FastLED's integer quirks deliberately. The memory layer's explicit `Arena&`
threading means misuse is a compile error or a debug-build trap rather than
silent corruption. The color subsystem's OKLCH chain is a real single-source-of-
truth design (forward via `fast_cbrt` on the hot path, exact inverse), and the
mesh kernel's arena-polarity contracts (`COMPOSITION POLARITY`, `SCRATCH ARENA
CONTRACT`) are the kind of comment that prevents future bugs. Findings here are
all low/medium (one AA-pad mismatch in `sdf.h`, two narrow animation edge cases,
a few doc-accuracy drifts).

**Effects (28 headers).** Each effect is a thin, effect-specific composition over
a small shared vocabulary (`Timeline`, `Pipeline<Filters...>`, `Presets<>`,
`BakedPalette`/`StaticPalette`) rather than reinvented plumbing — genuinely
elegant for embedded graphics. White-box test hooks target exactly the numeric
invariants the smoke harness can't observe. One finding (a Liquid2D preset below
its registered slider minimum).

**Hardware drivers & the Phantasm sync protocol.** Exemplary pure-math/hardware-
shell split that directly enables >2,500 lines of host tests for this component
alone. The 1-wire flywheel sync protocol (count-coded symbols, epoch-counted
playlist, "fail to missed, never to wrong" alphabet) is the technical high point
of the project. Zero findings in the driver code; the one hardware-adjacent high
is in the Phantasm *sketch's* effect-swap allocation.

**WASM target & simulator (daydream).** The `EngineHost`/`pixel_view`/`param_sync`
and `segment_layout`/`worker_protocol` DOM-free modules are cleanly separated from
the necessarily-imperative Three.js/lil-gui wiring, and the tricky invariants
(zero-copy buffer liveness under `ALLOW_MEMORY_GROWTH`, generation fencing,
sync-vs-drag) are the parts that carry unit tests. Findings here are dead code,
one arg-order footgun, and doc nits.

**Tests, scripts, tooling.** The test harness is minimal and self-verifying (it
fails loudly if a module is silently dropped or a test gutted). The one real
process gap: `capture_screenshots.mjs` crashes in its own summary block on every
run, and its CI guard is `node --check` (syntax only), which cannot catch it.

---

## 4. Prioritized Fix List

All 19 findings, grouped by priority and numbered sequentially. Each cites
`file:line`, the defect, and the recommended fix.

### High priority

1. ✅ **Phantasm effect-swap uses throwing `new` with no null/exception handling on a `-fno-exceptions` build** — `targets/Phantasm/Phantasm.ino:68`.
   `construct_effect<E>()` runs on every effect switch (~every 120 s for the life of the device) and does `E *e = new E();` with no `std::nothrow` and no null check, then immediately calls `e->init()`; `hardware/pov_segmented.h:283` derefs it again (`cur->height()`). `platformio.ini` builds with `-fno-exceptions`, so on OOM there is either an abort with no diagnostic or a silent `nullptr` → hard-fault null-deref. The file's own `setup()` guards `g_pov` against exactly this two lines earlier. Fix: mirror that pattern — `E *e = new (std::nothrow) E(); HS_CHECK(e != nullptr, "effect allocation failed (OOM)");` before `init()`.

2. ✅ **`capture_screenshots.mjs` crashes with `ReferenceError` in its own summary block on every run** — `scripts/capture_screenshots.mjs:79`.
   `RESOLUTIONS`, `targets`, `failures`, `blanks`, `wrongRes` are `const`/`let`-declared *inside* the `try` block (79–238) but read at module top level (240–292), so the entire diagnostic/gating section is dead and every run — success or failure — throws after `browser.close()` and exits non-zero with a raw stack trace. The CI guard (`node --check`) is syntax-only and cannot catch this. Fix: hoist those declarations above the `try {` and convert the inner `const` to assignments; re-run `npm run screenshots` to confirm the summary executes and only fails on a real failure/blank/wrong-resolution.

### Medium priority

3. ✅ **`Face::clip_rejects` uses a narrower AA pad than `get_horizontal_intervals`, contradicting its own comment** — `core/render/sdf.h:1908`.
   `clip_rejects` widens each azimuth interval by `2π/W` (one pixel) while `get_horizontal_intervals` widens by `1.25·(2π/W)` (one pixel + `fast_atan2` slop). A face whose reach falls in that 0.25px gap can be whole-face-culled on segmented/clipped renders, silently dropping an edge column. Fix: make `clip_rejects` use `pw = 1.25f * (2.0f * PI_F / Wd)` so the emit/cull windows agree exactly.

4. ✅ **`RandomWalk` silently freezes when the per-frame angle is below `TOLERANCE`** — `core/animation/motion.h:594`.
   `RandomWalk::step()` calls `Rotation::animate()`, which constructs a fresh one-shot `Rotation` (with `last_angle=0`) each frame. `Rotation::step()`'s sub-`TOLERANCE` accumulation relies on `last_angle` persisting across frames; with a new object per frame, any `options.speed` below `1e-4` rad/frame is dropped forever — the exact "freeze very slow rotations" case the mechanism exists to prevent. Reachable via the public `set_speed()`. Fix: drive the `Orientation` directly (unconditional `make_rotation`) or give `RandomWalk` a persistent `Rotation<W,CAP>` member.

5. ✅ **`Ripple` with duration 0 or 1 never becomes visible** — `core/animation/params.h:707`.
   `AnimationBase` normalizes duration 0→1, and `step()` increments `t` to 1 on the first call, so the envelope branch `if (t < duration)` never executes for a one-frame `Ripple`; `params.amplitude` stays 0 for its whole life with no error signal. The only shipping caller clamps to ≥30 frames, so it is latent, but `Ripple` is a general-purpose class with a public constructor. Fix: trap on a degenerate duration in the ctor, e.g. `HS_CHECK(duration == 0 || duration >= 2, ...)` — note the value is normalized, so both 0 and 1 must be rejected.

6. ✅ **`MeshMorph` silently mishandles an indefinite (`duration == -1`) construction** — `core/animation/mesh.h:58`.
   `AnimationBase` advertises `-1` as the legal "perpetual" sentinel, but `MeshMorph::step()` unconditionally computes `t/duration`; for `-1` that clamps to 0 forever, freezing the morph at its start pose with `done()` never true. Every sibling progress animation carries `HS_CHECK(duration >= 0, ...)`; `MeshMorph` is the lone omission. Fix: add `HS_CHECK(duration >= 1, "MeshMorph duration must be a positive frame count");` in the ctor.

7. ✅ **Liquid2D preset drives "Time Speed" below its registered slider minimum** — `effects/Liquid2D.h:293`.
   `register_param("Time Speed", &params.time_speed, 0.1f, 5.0f)` declares `[0.1, 5.0]`, but preset row 2 sets `time_speed=0.05f`, and `Presets::apply()`/`Lerp` write it without clamping. While that preset is active the live value sits outside the slider range, and the first GUI edit snaps it up to 0.1 — the exact drift `register_param`'s range-assert exists to prevent, and which `MeshFeedback` guards against with a compile-time `preset_in_ranges()` `static_assert`. Fix: raise the preset to `0.1f` or widen the registered minimum to `0.05f`, and consider adding a `preset_in_ranges()` guard.

### Low priority

8. ✅ **Dead unreachable guard with a misleading comment in `Orientation::upsample`** — `core/math/geometry.h:551`.
   The `if (count < 2) return;` block is unreachable: `num_frames >= 1` is a class invariant and the earlier `if (num_frames >= count) return;` already returns whenever `count == 1`. The comment describes a `0/0` path that cannot occur. Fix: delete the dead block (551–553).

9. ✅ **`HalfEdgeMesh(Arena&, const MeshState&)` doc claims topology "may be borrowed," but topology is always owned and this ctor never reads it** — `core/mesh/mesh.h:211`.
   `MeshState::topology` has no borrowed-view counterpart (only `face_counts`/`faces`/`face_offsets` do), and the ctor body never touches `mesh.topology`. Fix: reword to describe face-connectivity borrowing and drop the word "topology."

10. ✅ **`SolidBuilder::needle()` docstring describes the wrong operator** — `core/mesh/solids.h:332`.
    Says "(dual of truncate)"; `MeshOps::needle` is `kis(dual(...))` — Hart's `n = kd`, no truncate involved. Fix: change to "kis of the dual; n = kd" to match `conway.h`.

11. ✅ **README §7.7 miscounts the `islamic_registry` (23 vs actual 24) and the total (54 vs 55)** — `README.md:1139`.
    `islamic_registry` has 24 entries; total is 18+13+24 = 55. `NUM_ENTRIES` is computed at compile time so runtime is unaffected. Fix: update lines 1139 ("24 entries") and 1141 ("55 registered solids").

12. ✅ **`Filter::Pixel::Feedback` class doc claims `flush()` covers the full canvas, but it honors the clip band** — `core/render/filter.h:1232`.
    `flush()` restricts every pass to `cv.clip()`'s render band, and `test_feedback_flush_respects_clip` pins exactly that; the class-level `@details` contradicts both. Fix: reword to "iterates the full pixel grid *within the active clip band*."

13. ✅ **`Transition::from(mutant)` ctor initializer is dead, overwritten before any read** — `core/animation/params.h:29`.
    `from` is private, read only in `step()`, which unconditionally re-snapshots it on the first call (`captured` starts false). The ctor-time capture has no observable effect and is inconsistent with `ColorWipe`'s lazy `from_snap{}` idiom in the same file. Fix: default-initialize `from` and rely on the lazy capture.

14. ✅ **`builder.py`'s `transform()` docstring cites a nonexistent `calib.py`** — `hardware/phantasm/gen/builder.py:18`.
    "Verified empirically … (see calib.py)"; no `calib.py` exists in the repo (never committed), leaving a maintainer with no artifact to re-verify the load-bearing coordinate transform. Fix: commit the calibration script or replace the parenthetical with a concrete re-verification recipe.

15. ✅ **Dead static constants on `Daydream` (`PIXEL_WIDTH`, `UP`, `DOT_COLOR`)** — `daydream/driver.js:116`.
    None are read anywhere; `PIXEL_WIDTH` is even kept in sync on every `updateResolution()` with no observer, and `DOT_COLOR` misleadingly reads as if it configures dot color (the material has no `color` option; color comes from the WASM `instanceColor` buffer). Fix: delete the three fields and the `PIXEL_WIDTH` reassignment.

16. **`setResolution(w,h)` vs `updateResolution(h,w,…)` use opposite argument order in the same function** — `daydream/daydream.js:436`.
    Both calls are individually correct today, but the two "set the sphere resolution" entry points in `applyResolution()` disagree on parameter order — one transposition away from silently swapping W/H in a future refactor. Fix: align `updateResolution`'s signature to `(w, h, dotSize)` (matching the WASM bridge and `geometry.js` dims), or comment the intentional mismatch at both call sites.

17. **`composite()` return-value doc doesn't cover the fault case** — `daydream/segment_controller.js:657`.
    Doc says `0` means "every result was null/empty (a fully-fenced frame)," but `composite()` also returns 0 when the pre-pass rejects a segment (out-of-bounds/empty/inverted rect, pixel-length mismatch) and latches a fault — a distinct condition a maintainer could misdiagnose. Fix: extend the doc to note the fault path.

18. **Dead CSS for a removed slide-up code panel** — `daydream/tools/solids.html:139`.
    `#code-panel`, `#code-panel.open`, and bare `textarea` rules (139–169) style elements that no longer exist (the tool exports via per-item Copy buttons now). Fix: delete the three rule blocks.

19. **Leading-underscore variable name breaks the codebase's no-underscore-prefix JS convention** — `daydream/tools/solids.html:2082`.
    `_labelCam` is a plain module-scoped cache object, not a WASM export or numeric/string-data token, so it falls outside the documented underscore exception; it is the sole offender in `tools/`. Fix: rename to `labelCam` (and its references).

---

## 5. What NOT to "Fix" (Intentional by Design)

Recorded so future reviews don't re-raise these as defects — each was proposed and
rejected during validation:

- **Fail-fast `HS_CHECK` traps on cold paths** are the project's error-handling
  strategy, not "missing error handling." They deliberately survive `NDEBUG` and
  trap at the violation site; hot paths use stripped `assert` backed by a cold
  trap at the bind/setup site.
- **No heap on the render path / arena-only allocation** is a hard invariant;
  persistent-arena caches are expected.
- **The two-physical-buffer ISR double buffer** is correct by design (no RAM for a
  third buffer); the `buffer_free()` gate is not a tearing bug.
- **Compile-time `<W,H>` templating** trades binary size for zero runtime
  generality overhead, intentionally.
- **Bit-exact platform mocks** (`map`/`scale8`/`sin8`) match FastLED/Arduino's real
  integer semantics, including its unguarded div-by-zero/underflow behavior, on
  purpose — divergence would break sim/device parity.
- **Deterministic seed `Pcg32(1337)` across all boards/workers** is what makes the
  distributed render bit-identical; it is not a weak-RNG issue.

---

## 6. Bottom Line

The code is production-grade and, in its engine core and hardware layer,
reference-quality. The defect density is extraordinarily low for a system of this
complexity, and the defects that remain are shallow. Fixing the two high-priority
items (device-swap allocation guard, screenshot-script scoping) closes the only
issues with runtime consequence; the medium items harden general-purpose library
classes against inputs no shipping effect currently supplies; the low items are
housekeeping. This review recommends landing all 19 fixes, but none of them
threaten the integrity of the system as it ships today.
