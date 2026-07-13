# `strobe_columns()` audit (historical)

> **Historical audit (June 2026; renderer follow-up implemented).** Holosphere exposes the active
> flag through [`HolosphereEngine::strobeColumns`](../targets/wasm/wasm.cpp),
> [`daydream.js`](../../daydream/daydream.js) forwards it after effect changes, and
> [`driver.js`](../../daydream/driver.js) fills inter-column gaps only for persistent columns.

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
between frames is [`persist_pixels`](../core/render/canvas.h) (`advance_buffer` /
`clear_buffer`) — an orthogonal mechanism the old comments conflated this flag
with. An effect can independently be persist-buffer or clear-buffer *and*
strobe or persist on the strip.

## Audit snapshot — active effects

All non-legacy effects now **strobe** (`strobe_columns() == true`): each column
blanks to black after it is shown rather than persisting into the next. The four
that already strobed (Flyby, GnomonicStars, Liquid2D, ReactionDiffusion) are
unchanged; the remaining 22 were flipped from `false` to `true`.

| Effect | `strobe_columns()` |
|---|---|
| ChaoticStrings | `true` |
| Comets | `true` |
| DistortedRing | `true` |
| DreamBalls | `true` |
| Dynamo | `true` |
| FlowField | `true` |
| Flyby | `true` |
| GnomonicStars | `true` |
| HankinSolids | `true` |
| HopfFibration | `true` |
| IslamicStars | `true` |
| Liquid2D | `true` |
| MeshFeedback | `true` |
| MindSplatter | `true` |
| MobiusGrid | `true` |
| PetalFlow | `true` |
| Raymarch | `true` |
| ReactionDiffusionBase | `true` |
| &nbsp;&nbsp;├ BZReactionDiffusion | `true` (inherited) |
| &nbsp;&nbsp;└ GSReactionDiffusion | `true` (inherited) |
| RingShower | `true` |
| RingSpin | `true` |
| ShapeShifter | `true` |
| SphericalHarmonics | `true` |
| Thrusters | `true` |
| Voronoi | `true` |

## Audit snapshot — legacy effects ([core/engine/effects_legacy.h](../core/engine/effects_legacy.h))

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

**Pattern:** the active roster is now uniformly strobe — every non-legacy effect
renders as sharp slices with dark gaps between columns rather than smearing into
a continuous band. The legacy effects above keep their original mix: trail /
accumulation effects persist (`false`) while sparse or full-frame "crisp" effects
(the Matrix rain, StarsFade, Spiral, Burnout) strobe (`true`).

## Simulator behavior (implemented)

The WASM engine exposes the active effect's flag through
[`HolosphereEngine::strobeColumns`](../targets/wasm/wasm.cpp). Daydream reads it in
[`daydream.js`](../../daydream/daydream.js) and passes it to the renderer in
[`driver.js`](../../daydream/driver.js).

- **`strobe == true`:** the renderer keeps discrete round dots and dark horizontal gaps.
- **`strobe == false`:** the renderer extends each dot in the sweep direction with the
  `uColumnFillArc` shader uniform, filling the gap as a sample-and-hold band.

The framebuffer data is unchanged; only the geometry between adjacent columns differs.
