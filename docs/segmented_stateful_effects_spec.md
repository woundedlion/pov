# Stateful Effects in Segmented Mode — Design Spec

*Status: IMPLEMENTED. This document describes the design — now in the engine —
for making history-reading pixel effects (today only `MeshFeedback`, via the
`Pixel::Feedback` filter) render correctly under the simulator's segmented
mode, without dropping pixels, while non-stateful effects keep the full
clipping win. The seam is the WASM driver boundary (`targets/wasm/wasm.cpp`
`setClip`) plus two compile-time filter traits; the device path is already
correct and does not change.*

---

## 1. Goal

`MeshFeedback` — and any future effect whose per-frame state depends on
*other segments'* pixels — must produce identical output whether rendered as
one full-canvas instance or as N segment workers. Non-stateful effects must
retain segmented rendering's clipping win (out-of-band pixels never shaded).

---

## 2. Where segmentation actually lives

Segmentation-by-clipping is a **simulator-only** mechanism. The two render
paths diverge:

- **Device (`hardware/pov_segmented.h`).** Every Teensy renders the *full*
  `CANVAS_W × CANVAS_H` canvas (`cur->draw_frame()`, pov_segmented.h:321).
  The flywheel ISR's `render_column()` (pov_segmented.h:631) packs only this
  board's LEDs out of the full buffer via `pov::segment_x_col()`. There is no
  `set_clip` on the device. `MeshFeedback` is therefore **already correct in
  hardware** — each board independently computes the whole frame and samples
  its slice.

- **Simulator (`segment_worker.js` → `targets/wasm/wasm.cpp`).** Daydream
  reproduces the *partitioning* in software: one isolated WASM module instance
  per segment, each calling `setClip(x0, x1, y0, y1)` (wasm.cpp:465) so the
  rasterizer's scanline culling skips out-of-clip rows/columns. The readback
  copies the full canvas; `segment_worker.js` then extracts just the quadrant
  rectangle (`blitSegmentRect`) before transfer (README §10.7).

So "compute the full frame per segment" is **what the device already does**.
The simulator's clipping is an optimization that is sound for non-stateful
effects and *wrong* for unbounded-reach feedback.

---

## 3. Why clipping drops pixels for feedback

The feedback flush iterates only the clip band (filter.h:1487) but reads
`cv.prev` at warp offsets up to **±W/2 horizontally and ±H vertically**
(filter.h:1436 — a melt/swirl warp can pull a row most of the way across the
canvas). That reach is **unbounded relative to a segment band**: a worker
clipped to its band has stale or zero `prev` outside the band, so cross-band
trails read as black → dropped pixels and visible seams at the band edges.

Two consequences, both load-bearing for the design:

1. **The flush must run full-frame.** Correct output in a band requires a
   correct full prior frame, which requires the flush to have advanced the
   whole canvas last frame — recursively, every frame.

2. **The mesh draw must also run full-frame.** The mesh is the *source* that
   seeds the trails. If a worker skips drawing the mesh outside its band, the
   warp pulls that missing content into the band over successive frames, and
   the band's own display loses trail content that originated elsewhere.

There is therefore **no per-worker rasterization win for an unbounded-reach
feedback effect** — both the flush and its source draw are inherently
full-frame. The clipping win survives only for (a) non-stateful effects
(clip everything) and (b) the per-segment output slice, which is already done
JS-side and shades nothing. This is not a limitation of the design; it is the
same cost the device pays for distributed memory, and segmented WASM currently
*violates* the sim-mirrors-device invariant by pretending otherwise.

A finer reach distinction still matters for *other* history filters
(§5): `Screen::Trails` decays pixels in place (reach 0 — band clipping is
already correct), whereas `Pixel::Feedback` (warp) and `World::Trails`
(reprojected under rotation) move content across bands.

---

## 4. Design

### 4.1 Derive cross-segment reach as a compile-time trait

Statefulness is already a compile-time filter trait (`has_history`, with the
`Pipeline` constexpr-folding the history/terminal static-asserts at
filter.h:350). Add a sibling trait that captures the property that actually
matters here — whether a filter's state *moves across segment boundaries*:

- `static constexpr bool crosses_segments`, defaulting to `has_history`
  (fail-safe: a new history filter is treated as cross-segment until proven
  bounded). Add it to the four trait bases (`Is2D`/`Is3D`/`Is2DWithHistory`/
  `Is3DWithHistory`, filter.h:49) as `= has_history`.
- `true` on `Pixel::Feedback` — the load-bearing case, and the only one that
  *must* be `true`. It reads `cv.prev` (other segments' pixels).
- `true` on `World::Trails` — already `true` by default (`has_history`), so
  this is documentation, not a required override. Its store happens at
  `plot()` time (filter.h:778), upstream of projection; whether band clipping
  would actually corrupt it depends on whether the rasterizer culls
  out-of-band fragments *before* that store. Left `true` as the fail-safe
  default rather than relying on that analysis.
- `false` on `Screen::Trails` — the **only non-fail-safe override**: it turns
  full-frame *off* for a history filter, justified solely by reach 0 (decays
  in place, redraws at the same screen coordinate). If that reasoning is
  wrong, it drops pixels — so it must be tested directly (§8).

`Pipeline` does **not** today expose any aggregate trait constant — the
existing folding is per-node recursive `static_assert`s on `Head::has_history`
(filter.h:350), nothing more. So this trait needs a new recursive OR-fold
written from scratch: `static constexpr bool any_crosses_segments =
Head::crosses_segments || NextPipeline::any_crosses_segments;` in the recursive
node, **plus a `false` base case in the terminal `Pipeline<W,H>`**
(filter.h:111). There is no existing `any_*` member to mirror.

### 4.2 Expose it as one runtime query on `Effect`

```
virtual bool needs_full_frame() const { return false; }
```

There is no shared base hook — each effect that wants the gate adds the
override itself:

```
bool needs_full_frame() const override {
  return decltype(filters)::any_crosses_segments;
}
```

The override's *value* is fully trait-derived (no per-effect judgment): it reads
the pipeline's compile-time `any_crosses_segments` fold, so adding or removing a
filter updates the answer automatically. Only the one-line bridge is manual,
because `Effect` is type-erased — the driver holds an `Effect*` and the compile-
time fold lives only in the derived effect's `filters` member, so a base virtual
cannot read it without the derived type surfacing it. Going fully implicit would
mean reparenting every filtered effect onto a CRTP base that reads
`Derived::filters`, which also has to reach each effect's (private) `filters` —
strictly more intrusive than the one-liner, so it is not done.

Every effect whose pipeline crosses segments carries the identical override:
`MeshFeedback` (`Pixel::Feedback`) and the two `World::Trails` effects `Dynamo`
and `SplineFlow` (→ `true`). Everything else stays `false` by default. The
roster test (§8) pins exactly this set so a new cross-segment effect that forgets
the override is caught.

### 4.3 Honor it at the driver boundary (the only behavioral change)

In `targets/wasm/wasm.cpp` `setClip` (wasm.cpp:449): if
`currentEffect->needs_full_frame()`, leave the clip at the full canvas
(optionally record the requested band for telemetry) and return; otherwise
apply the band as today.

"Leave at full" is safe because the clip is already full when this fires: the
`Effect` constructor resets `clip` to the whole canvas (canvas.h:58), and the
worker re-applies the band *only* after `setEffect` rebuilds the effect
(`segment_worker.js` `applyClip`). So a full-frame effect's clip is never
narrowed in the first place — the early return preserves the constructor's
full clip rather than relying on resetting a stale band. (If that lifecycle
ever changes, harden this to an explicit `set_clip(0, H, 0, W)` instead.)

The elegance is that **the hot-path filter code needs no change**. With a
full clip, `XClip::active` is false and the coarse-row band spans every row,
so the flush's existing band-pruning (filter.h:1448, filter.h:1597) already
degrades to full-frame — the comments there note "a full canvas... does the
same work either way." `blitSegmentRect` still slices each worker's quadrant
from the full readback, unchanged. Every worker computes the bit-identical
full frame; only the slice differs — matching the device exactly.

**Precondition (existing invariant this rests on).** "Bit-identical full
frame across workers" holds only because per-worker frame inputs are already
deterministic: animations are *frame-stepped*, not wall-clock-stepped
(`AnimationBase::step` counts frames; `drawFrame()` advances exactly one,
wasm.cpp:501 — the `elapsed`/`renderUs` timings are telemetry and never feed
animation), the RNG is fixed-seed (`std::mt19937(1337)`), and params are
broadcast to every worker. The same invariant non-stateful segmented effects
already depend on; the design adds nothing new to it but is wholly dependent
on it.

### 4.4 Reuse dormant `set_margin` for the bounded tier

`set_margin` (canvas.h:127) exists but is exercised only in a test
(test_canvas.h:461). It is the home for the genuine "clip rasterization to
save time" case that *does* apply to neighborhood-reading effects: a bounded
spatial reach (AntiAlias ±1, small blurs) renders band + a declared margin
instead of full-frame. `MeshFeedback` is unbounded and skips this tier.

---

## 5. Policy summary

| Effect class | Reach | Simulator render bound |
|---|---|---|
| Non-stateful | none | segment band (+1 AA margin) — full clipping win |
| In-place history (`Screen::Trails`) | 0 (redraws at same coord) | segment band — no margin needed; `crosses_segments = false` |
| Bounded spatial neighborhood (AntiAlias ±1, small blurs) | finite | band + declared `set_margin` (§4.4, future) |
| Cross-segment history (`Pixel::Feedback`, `World::Trails`) | unbounded | full canvas; output sliced JS-side |

---

## 6. What changes, what does not

| Layer | Change |
|---|---|
| `core/filter.h` traits | add `crosses_segments = has_history` to the trait bases; `false` on `Screen::Trails`; add a **new** recursive `any_crosses_segments` OR-fold to `Pipeline` + a `false` base case in the terminal `Pipeline<W,H>` (no existing `any_*` to mirror) |
| `core/canvas.h` `Effect` | add `needs_full_frame()` virtual (default `false`) |
| `targets/wasm/wasm.cpp` `setClip` | branch on `needs_full_frame()` → full canvas vs band |
| flush / `scan.h` / `plot.h` hot paths | **none** — a full clip already degrades correctly |
| `hardware/pov_segmented.h` (device) | **none** — already full-frame |
| `segment_worker.js` slicing | **none** |

---

## 7. Open decision: redundant workers vs single-instance

The design above runs **N redundant, identical full-frame workers** for a
stateful effect. It is the simplest path, mirrors the device's independent
per-board compute, and is bit-identical to it.

The alternative: because all N workers compute the bit-identical frame for a
cross-segment effect, render it **once** and let all quadrants slice that
single readback, recovering the (N−1)× browser compute. This is faster in the
simulator but diverges the sim's topology from the device (hardware genuinely
has N independent computes) and complicates worker orchestration (a mixed
segmented / single-instance mode keyed on `needs_full_frame()`).

**Recommendation:** take the redundant-workers path, consistent with the
project's sim-mirrors-device doctrine. The single-instance variant is the
lever to pull only if browser cost during stateful effects outweighs
topological fidelity.

---

## 8. Test plan

Implemented in `tests/test_filter.h` and `tests/test_effects.h`:

- **Trait fold** (`test_crosses_segments_trait_and_fold`, test_filter.h): pins
  the per-filter `crosses_segments` values (incl. the `Screen::Trails == false`
  override) and the `Pipeline::any_crosses_segments` OR-fold — `true` for the
  MeshFeedback stack, `false` for a non-stateful stack and a `Screen::Trails`
  stack — so a future filter addition can't silently regress the gate.
- **Roster gate** (`test_needs_full_frame_gate`, test_effects.h): constructs the
  real effects and asserts `needs_full_frame()` is `true` for exactly the
  cross-segment set (`MeshFeedback`, `Dynamo`, `SplineFlow`) and `false` for
  representative non-stateful effects. This is the end-to-end equivalent of a
  `setClip` test — the WASM driver reads exactly this query, and `setClip` itself
  lives in the Emscripten-only TU, which the native suite cannot link. A new
  cross-segment effect that forgets the override is caught here.
- **`Screen::Trails` banded-vs-full bit-identity**
  (`test_screen_trails_banded_matches_full`, test_filter.h): the load-bearing
  proof of the one non-fail-safe override. The same multi-frame seed sequence is
  driven through one full-canvas instance and two band-clipped instances; the
  stitched banded output must equal the full output byte-for-byte. This is the
  override that drops pixels if the reach-0 reasoning is wrong.
- **Feedback band-clip divergence**
  (`test_feedback_banded_diverges_from_full`, test_filter.h): the inverse for an
  unbounded-reach filter — a melt warp drips content across the segment boundary,
  so a band-clipped `Pixel::Feedback` worker's bottom band differs from the
  full-frame render. Demonstrates *why* full-frame is required: band clipping
  here is not equivalent, so the gate is load-bearing, not cosmetic.
- `test_effect_needs_full_frame_default_false` (test_filter.h) pins the base
  `Effect::needs_full_frame()` default.
