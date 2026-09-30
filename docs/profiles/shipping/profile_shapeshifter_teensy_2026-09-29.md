# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-09-29 19:01 on COM3.
Replaces the 2026-09-29 16:05 capture of `45a733fe6`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; dense stars stroke through `Plot::PlanarChords` (`core/render/plot/chords.h`), whose walk is `HS_HOT_FLASH_MEMBER` (flash, not `cold`) |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, master `bf84131e3` + the `Plot::PlanarChords` change in this report |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:78128, data:153712, headers:8800` / `RAM1: variables:315040, code:42248, padding:23288, free:143712` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 1921–1936 root counter cyc ÷ 600 MHz matches the measured wall sum within **3.6 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `ss_draw_all` worst window 36.33 ms/f. Peak frame render is **40.45 ms** (frames 2289–2304), and **0/2448** frames spilled. Setup frame 1 rendered 35.08 ms.

The previous capture of master `45a733fe6` (2026-09-29 16:05) recorded peak 🟢 58.88 (9) and spilled 🟢 0/2447 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 9 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2289–2304)

```
frame                     62.20 ms   37.32 Mcyc   100%
  pov_preserve_half       141.3 us    84.8 kcyc     0%
  ss_draw_all             32.69 ms   19.61 Mcyc    53%
    ss_plot_dispatch      32.49 ms   19.49 Mcyc    52%  x161  121256 cyc/c
  ss_timeline_step         31.6 us    19.0 kcyc     0%
  ss_buffer_wait          29.33 ms   17.60 Mcyc    47%
    canvas_clear           84.3 us    50.6 kcyc     0%
    canvas_buffer_wait    29.25 ms   17.55 Mcyc    47%
```

Wall min/avg/max = 48.84/62.20/75.86 ms; render avg/max = 32.95/40.45 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `ss_draw_all` is the costliest modal-call-count window of each entry. The cycle wrapped to entry 1 (10 `Preset:` markers).

| Entry | Shape (count) | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | Planar Star (208, screen balanced) | 🟢 40.45 | 0/470 | 33.59 | 28/30 |
| 9 | Flower (72) | 🟢 37.19 | 0/239 | 36.33 | 14/15 |
| 4 | Flower (70) | 🟢 35.64 | 0/239 | 34.62 | 14/15 |
| 8 | Spherical Polygon (144) | 🟢 25.37 | 0/239 | 23.99 | 14/14 |
| 7 | Spherical Polygon (144) | 🟢 21.87 | 0/239 | 20.23 | 14/15 |
| 6 | Spherical Polygon (128) | 🟢 18.23 | 0/239 | 17.55 | 14/15 |
| 2 | Spherical Polygon (74) | 🟢 15.06 | 0/305 | 13.18 | 18/19 |
| 5 | Planar Star (72) | 🟢 8.35 | 0/239 | 6.36 | 14/15 |
| 3 | Planar Star (43) | 🟢 6.87 | 0/239 | 5.28 | 14/15 |

The three dense planar-star entries moved (1: 58.88 → 40.45, 5: ~11.4 → 8.35, 3: ~9.5 → 6.87); the Flower and Spherical Polygon entries render through the unchanged `Plot::rasterize` path and match the previous report within run-to-run noise. Entry 1 is still the effect's peak, with the two Flower entries 3–4 ms behind it.

### Pinned worst case

Entry 1 holds for 240 frames per visit, so the cycle samples a narrow band of orientations. Pinned with `-D HS_PROFILE_PRESET=0` over 45 s (688 frames of the random-walk orientation), the same entry peaked at **63.28 ms with 4/688 frames spilled on master `45a733fe6`** and **42.62 ms with 0/688 spilled** with this change. The peak orientation puts the star's center near a pole, where the innermost contours are small rings the rasterizer has to walk at its pole step floor.

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

Master to this report, full cycle: peak **58.88 → 40.45 ms (−18.43 ms, −31.3%)**; shipping `phantasm` image: FLASH code +3,920 B, ITCM +16 B, DTCM 0, stack headroom unchanged at 12,896 B; `[teensy-gate] phantasm: PASS`. `Plot::PlanarChords` binds its chart and scratch tables (668 B for a 16-point star) in the persistent arena through `init_storage()`, never on the stack.

### Not landed: clip-dependent edge cuts (owner decision)

Two further levers cut an edge at clip-dependent points and rasterize only the runs that can reach the band:

| Lever | Measured | Image vs master (segment clips) |
|---|---|---|
| Pole runs split into 8 chart pieces, band-culled | entry 1 pinned 43.20 → 41.05 ms | PSNR 48–51 dB, energy ratio 1.000 |
| Flower edges split into 24/`sides` chart pieces, band-culled through per-edge `RasterProjection::planar` flags (a one-line core API addition); masked Flower samples 63% → 33% | Flower entry 9 pinned 38.83 → 26.33 ms | PSNR 49–51 dB, energy +0.2% |

Both restart the adaptive walk at a clip-dependent cut, which moves the sample phase inside the visible run, so a clipped render no longer reproduces the unclipped frame pixel for pixel: `unit_shapeshifter_tiles` fails by design. Splitting at fixed points instead keeps that parity but brightens Flowers by 3–7% (every restart adds the rasterizer's start sample and end half-step), so it is not an option. They land only if the tile-parity invariant is relaxed to a budget for these two paths; together they would take the pinned entry-1 worst case to ~41 ms and the Flower entries from ~38 to ~26 ms.

## Column-ISR / DMA marshaling cost

```
isr_wake         1148/frame  min/avg/max 0.5/1.6/11.1 us  cpu 3.03%
isr_pack          143/frame  min/avg/max 6.2/6.9/9.4 us  cpu 1.58%
isr_dma_submit    143/frame  min/avg/max 0.6/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 53% of the peak window, 32.69 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `ss_timeline_step` — 0% of the peak window, 0.05 ms/f.

At the pinned entry-1 peak the remaining cost splits roughly into `Plot::rasterize` samples inside pole pieces (~850 cycles each), chord-walk splats (~370 cycles each, ~270 of them the sink's `plot`), and anchor projection. The rasterizer's per-sample cost is shared core code: it pays roughly five divides and three square roots per sample plus `vmrs` syncs from float `std::min/max` in `screen_step_components`.

README cells: peak 🟢 40.45 (9), spilled 🟢 0/2448 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: this effect's own `HS_O3` regions are unchanged; the dense-star walk runs from flash (`.text.hot`), not ITCM.
- Pole thresholds: the physical display drops ~3.6° at each pole, so measuring the pole bands from the true pole rows moves them ~3 rows outward on device. The ideal-profile oracle is bit-identical; a physical-profile build of the oracle stays inside every dense-star budget but sits closer to it (pole-centered 144-contour cases: MAE 183.1/207.8, high-error pixels 206/270; master 156.0 and 168).
- Tile parity (`unit_shapeshifter_tiles`) is exact for every kept lever and now covers the shipping 208-contour screen-balanced star at a pole-centered and a general orientation.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock. Add `-D HS_PROFILE_PRESET=<i>` to pin one entry, and `HS_PROFILE_DEEP=1` for the `plot_chord_walk` / `plot_chord_pole` per-edge scopes.
