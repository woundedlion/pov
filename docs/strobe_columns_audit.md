# `strobe_columns()` audit

`strobe_columns()` (formerly `show_bg()`) is a **POV-display control**, read only
on the hardware column-draw path ([pov_single.h](../hardware/pov_single.h),
[pov_segmented.h](../hardware/pov_segmented.h)), once per swept column:

- **`true` — strobe:** after a column is lit, the spinning strip is immediately
  re-shown black (`FastLED.showColor(black)` on the FastLED path, or a trailing
  all-black DMA frame on the HD107S/DMA path). Each column is shown only for its
  own slice width and then blanked, so columns read as **sharp slices with dark
  gaps between horizontal neighbors**.
- **`false` — persist:** the lit column is **held on the strip until the next
  column overwrites it**, so it **smears horizontally** across the gap into its
  neighbor — adjacent columns bleed together into a continuous, gap-free band.

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

**Pattern:** trail / accumulation effects choose `false` (persist → columns
smear horizontally into a continuous, gap-free sweep); sparse or full-frame
"crisp" effects (star fields, the Matrix rain, reaction fields, Flyby's
stereographic projection) choose `true` (strobe → sharp slices with dark gaps
between columns). The values look internally consistent; the spinning hardware
is the final arbiter for any borderline case (e.g. whether the full-frame
shaders Liquid2D / ReactionDiffusion truly read better as gapped slices than as
a smeared band).

## The simulator reproduces `strobe_columns() == true`, not `false`

### Why

The WASM engine reads back the **entire display buffer** every frame
([wasm.cpp `drawFrame`](../targets/wasm/wasm.cpp)) and hands it to JS, which maps
each `(x, y)` pixel straight onto a **discrete LED dot** on the sphere mesh
([daydream.js](../../daydream/daydream.js)). Each dot is a separate point with
**black space between it and its horizontal neighbors** — the simulator never
fills the angular gap between adjacent columns.

That discrete, gapped rendering **is exactly what `strobe == true` looks like on
hardware**: each column is lit only for its own narrow slice of the sweep, then
blanked to black, so adjacent columns are separated by dark gaps. For strobe
effects the simulator is **already faithful** — no change needed.

What the simulator does **not** reproduce is `strobe == false`. There, each lit
column **persists on the strip until the next column is drawn**, so it smears
horizontally across the inter-column gap — neighboring columns bleed together
into a continuous, gap-free band. The simulator instead keeps those gaps dark,
so persist effects look more "pixelated"/sharper in the sim than on hardware.

So the error is the opposite of a brightness problem: it is purely **horizontal
extent**. Strobe = isolated slices with gaps (what the dots already show);
persist = each slice held/smeared forward into the next (what the dots omit).

### How to make it accurate

Gate the **horizontal inter-dot gap** on `strobe_columns()`:

- **`strobe == true`:** render as today — discrete dots, dark gaps between
  horizontal neighbors. No change.
- **`strobe == false`:** fill each horizontal gap with a **sample-and-hold of the
  trailing (earlier-swept) column** — extend each column's color *forward in the
  sweep direction* until it meets the next column, so horizontally adjacent dots
  merge into a continuous band. (The hold direction follows the sweep: a column
  persists from the moment it is drawn until the next one overwrites it, i.e.
  forward.) A short linear blend between the two columns, instead of a hard hold,
  mimics the finite LED on-time as it physically travels the gap.

Concrete renderings of the gap-fill, simplest first:

1. **Widen each dot horizontally** (in the sweep direction) to reach its next
   neighbor when `strobe == false`; leave dots point-sized when `true`.
2. **Interpolated fill geometry / a textured band** between horizontally adjacent
   dots, shaded by hold (or the short blend above) for persist, left black for
   strobe.

This is a **renderer-side** concern: the gap is an artifact of drawing discrete
dots, not of the framebuffer, so it belongs in daydream's dot geometry/shading,
not in the effect or the readback. The renderer needs the flag per active
effect — expose it with a small `.function("strobeColumns", …)` embind binding
(or carry it alongside the pixel readback) and branch the gap-fill on it. The
pixel data itself is unchanged; only how the space *between* dots is shaded.
