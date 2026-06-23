# `strobe_columns()` audit

`strobe_columns()` (formerly `show_bg()`) is a **POV-display control**, read only
on the hardware column-draw path ([pov_single.h](../hardware/pov_single.h),
[pov_segmented.h](../hardware/pov_segmented.h)), once per swept column:

- **`true` — strobe:** after a column is lit, the spinning strip is immediately
  re-shown black (`FastLED.showColor(black)` on the FastLED path, or a trailing
  all-black DMA frame on the HD107S/DMA path). Each column reads as a **sharp
  bright slice with dark inter-column gaps**, and integrated brightness is lower
  (reduced duty cycle).
- **`false` — persist:** the lit column is **held on the strip until the next
  column overwrites it**, filling its full angular cell — a gap-free,
  full-brightness image.

It is **not** framebuffer clearing. Whether the canvas is cleared or carried
between frames is [`persist_pixels`](../core/canvas.h) (`advance_buffer` /
`clear_buffer`) — an orthogonal mechanism the old comments conflated this flag
with. An effect can independently be persist-buffer or clear-buffer *and*
strobe or persist on the strip.

## Audit — active effects

| Effect | `strobe_columns()` | Nature / plausibility |
|---|---|---|
| ChaoticStrings | `false` | trail strokes on black — persist for brightness |
| Comets | `false` | trailed sprites — persist |
| DistortedRing | `false` | ring over prior frame — persist |
| DreamBalls | `false` | fading sprites — persist |
| Dynamo | `false` | trail-based — persist |
| FlowField | `false` | trails fill the frame — persist |
| **Flyby** | **`true`** | bright stereographic field; documented — crisp slices vs smearing across revolutions |
| **GnomonicStars** | **`true`** | sparse star field — strobe keeps star points crisp against black |
| HankinSolids | `false` | manages its own backdrop — persist |
| HopfFibration | `false` | trails over cleared frame — persist |
| IslamicStars | `false` | mesh supplies all color — persist |
| **Liquid2D** | **`true`** | full-frame shader — strobe (crisp, accepts dimming) |
| MeshFeedback | `false` | feedback flush owns the frame — persist |
| MindSplatter | `false` | manages own canvas — persist |
| Moire | `false` | rings over transparent — persist |
| MobiusGrid | `false` | no background — persist |
| PetalFlow | `false` | paints over cleared frame — persist |
| Raymarch | `false` | ray march fills the frame — persist |
| **ReactionDiffusionBase** | **`true`** | full reaction field — strobe (base class) |
| &nbsp;&nbsp;├ BZReactionDiffusion | **`true`** | inherited from base |
| &nbsp;&nbsp;└ GSReactionDiffusion | **`true`** | inherited from base |
| RingShower | `false` | rings accumulate — persist |
| RingSpin | `false` | clear background — persist |
| ShapeShifter | `false` | rings over cleared frame — persist |
| SphericalHarmonics | `false` | sphere painted directly — persist |
| SplineFlow | `false` | trails accumulate — persist |
| Thrusters | `false` | draws over persistent ring — persist |
| Voronoi | `false` | opaque, fills every pixel — persist |

## Audit — legacy effects ([core/effects_legacy.h](../core/effects_legacy.h))

| Effect | `strobe_columns()` |
|---|---|
| ChainWiggle | `false` |
| RingTwist | `false` |
| **TheMatrix** | **`true`** |
| Curves | `false` |
| **StarsFade** | **`true`** |
| **Spiral** | **`true`** |
| WaveTrails | `false` |
| RingTrails | `false` |
| Kaleidoscope | `false` |
| RingRotate | `false` |
| **Burnout** | **`true`** |
| Fire | `false` |
| DotTrails | `false` |
| Spinner | `false` |

**Pattern:** trail / accumulation effects choose `false` (persist → full
brightness, smooth sweep); sparse or full-frame "crisp" effects (star fields,
the Matrix rain, reaction fields, Flyby's stereographic projection) choose
`true` (strobe → sharp slices against black, accepting the duty-cycle dimming).
The values look internally consistent; the spinning hardware is the final
arbiter for any borderline case (e.g. whether the full-frame shaders
Liquid2D / ReactionDiffusion truly read better strobed than persisted).

## The simulator does not respect `strobe_columns()`

### Why

The WASM engine reads back the **entire display buffer** every frame
([wasm.cpp `drawFrame`](../targets/wasm/wasm.cpp)) and hands it to JS, which maps
each `(x, y)` pixel straight onto an LED dot on the sphere mesh
([daydream.js](../../daydream/daydream.js)). The simulator therefore shows the
**complete, ideal 2-D framebuffer at once** — every column lit, at full
brightness, simultaneously. It models persistence-of-vision as *perfect*.

`strobe_columns()` only changes what the **physical strip does between column
draws** (the trailing black frame). The simulator never reproduces the per-column
sweep or the inter-column blanking, so `true` and `false` render **identically**.

Concretely, relative to hardware integrated over one revolution:

- For `false` (persist) effects — the majority — the simulator is **already
  correct**: persist fills every column's cell at full brightness, which is
  exactly the gap-free full image the sim draws.
- For `true` (strobe) effects, the simulator is **wrong in the optimistic
  direction**: it shows them gap-free and full-brightness, whereas hardware
  strobes each column to a thinner, dimmer slice.

### How to make it accurate

The difference is fundamentally about **per-column sweep timing**, which the
current whole-frame readback discards. Options, cheapest to most faithful:

**Option A — static duty-cycle dimming (low value).** For strobe effects, scale
the whole frame's brightness down by an assumed duty cycle. Trivial, but the dot
mesh has one dot per column, so there is no sub-column space to render the dark
inter-column *gap* — only the dimming is representable, and the true duty cycle
isn't modeled. Rejected as cosmetic.

**Option B — temporal decay on the dot buffer (middle ground).** Stop replacing
the dot colors wholesale each frame; when the active effect strobes, fade the
previous dot colors toward black before compositing the new frame. Approximates
the strobe's darker, trailing look without a sweep model. Cheap, but it fades the
whole frame uniformly rather than per-column, so it is only a rough match.

**Option C — model the sweep with a persistent displayed-LED buffer
(recommended — the only genuinely faithful option).** Mirror the hardware ISR in
the simulator:

1. Keep a persistent *displayed* buffer (what is currently lit on the strip),
   distinct from the effect's freshly rendered frame buffer.
2. Advance a sweep column index by *N* columns per rendered frame (N chosen to
   match a target RPM / column rate).
3. For each swept column `x`: copy the effect's column `x` into the displayed
   buffer; if `strobe_columns()` is true, blank that column back to black in the
   displayed buffer (immediately, or after a short on-time fraction to model the
   duty cycle).
4. Render the persistent displayed buffer to the dots.

This reproduces **both** behaviors from one mechanism: persist effects show the
smooth full image (a column stays lit until re-swept), strobe effects show a
moving bright column trailed by black. It also naturally surfaces sweep-timing
artifacts (tearing, `advance_display` boundaries) the current sim hides.

**Where to implement it.** Best contained **inside the WASM engine**:
`drawFrame()` already holds `currentEffect`, so `currentEffect->strobe_columns()`
is in scope right at the readback, and the persistent displayed buffer + sweep
loop can live next to `pixelBuffer`. JS needs **no change** — it keeps rendering
`getPixels()`. (Alternatively, add a `.function("strobeColumns", …)` embind
binding and run the sweep JS-side, but keeping it in the engine avoids a new
per-frame round trip and keeps the sweep logic next to the hardware-mirroring
code.) Expose N / RPM as a tunable, and surface a UI note that strobe effects are
legitimately dimmer — otherwise a faithful sim will look like a regression.
