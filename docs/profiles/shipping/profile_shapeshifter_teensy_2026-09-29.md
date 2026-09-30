# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-09-29 21:00 on COM3.
Replaces the 2026-09-29 19:01 capture of `b00994c53`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; dense stars stroke through `Plot::PlanarChords` (`core/render/plot/chords.h`), whose walk is `HS_HOT_FLASH_MEMBER` (flash, not `cold`); Flowers split petal edges against the clip band through `Plot::PlanarBandSplit` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, master `d58a55049` + the `Plot::PlanarBandSplit` Flower change in this report |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:81536, data:153896, headers:8280` / `RAM1: variables:315040, code:42312, padding:23224, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2273–2288 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.1 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` worst window 33.62 ms/f (frames 2273–2288). Peak frame render is **40.57 ms** (frames 2289–2304), and **0/2448** frames spilled. Setup frame 1 rendered 35.18 ms.

The capture of `b00994c53` (2026-09-29 19:01) recorded peak 🟢 40.45 (9); master `45a733fe6` before this work (16:05) recorded peak 🟢 58.88 (9). Both spilled 0.

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                     62.19 ms   37.32 Mcyc   100%
  pov_preserve_half       138.2 us    82.9 kcyc     0%
  ss_draw_all             32.79 ms   19.67 Mcyc    53%
    ss_plot_dispatch      32.52 ms   19.51 Mcyc    52%  x161  121396 cyc/c
  ss_timeline_step         30.2 us    18.2 kcyc     0%
  ss_buffer_wait          29.23 ms   17.54 Mcyc    47%
    canvas_clear           84.4 us    50.7 kcyc     0%
    canvas_buffer_wait    29.15 ms   17.49 Mcyc    47%
```

Wall min/avg/max = 48.74/62.19/75.90 ms; render avg/max = 33.05/40.57 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `ss_draw_all` is the costliest modal-call-count window of each entry. The cycle wrapped to entry 1 (10 `Preset:` markers).

| Entry | Shape (count) | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | Planar Star (208, screen balanced) | 🟢 40.57 | 0/470 | 33.62 | 28/30 |
| 8 | Spherical Polygon (144) | 🟢 25.40 | 0/239 | 24.00 | 14/14 |
| 9 | Flower (72) | 🟢 24.30 | 0/239 | 23.85 | 14/15 |
| 4 | Flower (70) | 🟢 23.56 | 0/239 | 22.87 | 14/15 |
| 7 | Spherical Polygon (144) | 🟢 21.87 | 0/239 | 20.23 | 14/15 |
| 6 | Spherical Polygon (128) | 🟢 18.24 | 0/239 | 17.56 | 14/15 |
| 2 | Spherical Polygon (74) | 🟢 15.06 | 0/305 | 13.19 | 18/19 |
| 5 | Planar Star (72) | 🟢 8.42 | 0/239 | 6.43 | 14/15 |
| 3 | Planar Star (43) | 🟢 6.92 | 0/239 | 5.33 | 14/15 |

Against master `45a733fe6`, the dense planar-star entries moved (1: 58.88 → 40.57, 5: ~11.4 → 8.42, 3: ~9.5 → 6.92) and the Flower entries moved (9: 37.19 → 24.30, 4: 35.64 → 23.56). The Spherical Polygon entries render through the unchanged `Plot::rasterize` path and match within run-to-run noise. Entry 1 is still the effect's peak; the next entries sit near 25 ms.

### Pinned worst case

Entry 1 holds for 240 frames per visit, so the cycle samples a narrow band of orientations. Pinned with `-D HS_PROFILE_PRESET=0` over 45 s (688 frames of the random-walk orientation), the same entry peaked at **63.28 ms with 4/688 frames spilled on master `45a733fe6`** and **42.62 ms with 0/688 spilled** with this change. The peak orientation puts the star's center near a pole, where the innermost contours are small rings the rasterizer has to walk at its pole step floor.

### Entry 1 at 288 contours

`-D HS_PROFILE_PRESET=0 -D HS_PROFILE_SHAPESHIFTER_COUNT=288` pins entry 1 at the Count slider's maximum (two contours per row) for 45 s:

| Build | Peak render ms | Mean render ms | Spilled/frames |
|---|--:|--:|--:|
| master `45a733fe6`, before this work | 🔴 87.40 | 64.36 | 263/435 (60%) |
| this report | 🟢 58.82 | 38.96 | 0/696 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## What changed and what each lever measured

All figures are the pinned entry-1 peak (45 s, `HS_PROFILE_PRESET=0`) and the `profile`-env RAM1 code (ITCM) size, each against the step before it, measured during development on master `a2c6d2be0` (63.49 ms, ITCM 44,232 B). Rebased onto `45a733fe6` the set measured 63.28 → 42.97 ms; moving the walk into `Plot::PlanarChords` (bit-identical output under IEEE float) measured 42.62 ms, ITCM −240 B.

| Lever | Peak ms | Δ ms | ITCM Δ | Kept |
|---|--:|--:|--:|---|
| Edge cull: a 1-Lipschitz bound on each chart-straight edge skips edges that cannot splat into the clip band; the chord walk plots only the sample-index range that can reach it; star vertices charted from the angle recurrence instead of per-vertex `azimuthal_project` (`cold` placement) | 52.26 | −11.23 | +544 | yes |
| `HS_HOT_FLASH_MEMBER` instead of `HS_FLASH_MEMBER` (drops `cold`), frame setup in flash | 51.24 | −1.02 | −512 | yes |
| Analytic chart alone (A/B against `azimuthal_project` vertices) | — | −1.20 | 0 | yes |
| Pole edges: only anchor intervals within 11 rows of a pole go to `Plot::rasterize`; the rest of the edge keeps the chord walk | 49.9 | −1.3 | +16 | yes |
| Pole thresholds measured from the pole rows instead of row 0 / H−1 (identical on the ideal profile; on the device's physical profile the poles sit ~3 rows outside the display) | 43.20 | −6.7 | 0 | yes (see Caveats) |
| Same pole split with the out-of-line lambdas left in ITCM | 47.80 | −0.28 vs flash | +3,648 | no |
| First pole split at 3 rows with sag bisection | 48.08 | — | 0 | no: oracle energy drift over budget |
| Rasterizer-faithful chord walk inside pole pieces (reproduces `screen_step_components` on the raw sampled position) | — | −1.17 | −288 | no: pole-centered oracle high-error pixels 365/270 |
| Batching a contour's pole runs into one `rasterize` call | — | −0.37 | 0 | no: within noise |
| Culling whole pole runs against the band before the call | — | +0.07 | 0 | no: the rasterizer's own cull already does it |

Master `45a733fe6` to this report, full cycle: peak **58.88 → 40.57 ms (−18.31 ms, −31.1%)**; shipping `phantasm` image: FLASH code +6,384 B, ITCM −96 B, DTCM 0, stack headroom unchanged at 12,896 B; `[teensy-gate] phantasm: PASS`. `Plot::PlanarChords` and `Plot::PlanarBandSplit` bind their tables (668 B and 50 B) in the persistent arena through `init_storage()`, never on the stack.

### Flower band split

Flower petal edges go through `Plot::rasterize`, whose clip cull decides per whole edge, so a long petal edge that touches the band anywhere was walked end to end: 63% of the Flower samples splatted nothing into the band. `Plot::PlanarBandSplit` cuts each chart-straight edge into 24/`sides` chart pieces, tests each piece with `ClipBand::may_reach`, and joins the runs; runs that cannot reach the band go to `Plot::rasterize` as edges flagged invisible through `RasterProjection::planar(basis, flags)`, which it skips unsampled. Masked Flower samples fall to 33%.

| Measure | master `d58a55049` | this report |
|---|--:|--:|
| Entry 9 pinned 45 s, peak / mean render | 37.88 / 34.55 ms | 24.52 / 23.22 ms |
| Entry 9 in the cycle, peak render | 37.19 ms | 24.30 ms |
| Entry 4 in the cycle, peak render | 35.64 ms | 23.56 ms |

A visible run starts at a cut that depends on the clip, which moves the walk's sample phase inside it, so a board's Flower strokes match the unclipped render to a fraction of a pixel, not bit for bit: PSNR 50–52 dB, energy within 0.2%, and no stroke dropped outside the pole rows, where a sub-pixel shift moves a dot by whole columns. `test_planar_band_split_matches_whole_polyline` pins the energy (3.5% per quadrant) and coverage; `unit_shapeshifter_tiles` exempts the candidate Flower from exact parity. The same cut applied to entry 1's pole runs (pinned 43.20 → 41.05 ms) is not included.

## Column-ISR / DMA marshaling cost

```
isr_wake         1147/frame  min/avg/max 0.5/1.6/11.2 us  cpu 3.02%
isr_pack          143/frame  min/avg/max 6.2/6.9/9.4 us  cpu 1.58%
isr_dma_submit    143/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 53% of the peak window, 32.79 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.05 ms/f.

At the pinned entry-1 peak the remaining cost splits roughly into `Plot::rasterize` samples inside pole pieces (~850 cycles each), chord-walk splats (~370 cycles each, ~270 of them the sink's `plot`), and anchor projection. The rasterizer's per-sample cost is shared core code: it pays roughly five divides and three square roots per sample plus `vmrs` syncs from float `std::min/max` in `screen_step_components`.

README cells: peak 🟢 40.57 (9), spilled 🟢 0/2448 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: this effect's own `HS_O3` regions are unchanged; the dense-star walk runs from flash (`.text.hot`), not ITCM.
- Pole thresholds: the physical display drops ~3.6° at each pole, so measuring the pole bands from the true pole rows moves them ~3 rows outward on device. The ideal-profile oracle is bit-identical; a physical-profile build of the oracle stays inside every dense-star budget but sits closer to it (pole-centered 144-contour cases: MAE 183.1/207.8, high-error pixels 206/270; master 156.0 and 168).
- Tile parity (`unit_shapeshifter_tiles`) is exact for every star lever and covers the shipping 208-contour screen-balanced star at a pole-centered and a general orientation; the Flower band split is exempt (see above).
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock. Add `-D HS_PROFILE_PRESET=<i>` to pin one entry, and `HS_PROFILE_DEEP=1` for the `plot_chord_walk` / `plot_chord_pole` per-edge scopes.
