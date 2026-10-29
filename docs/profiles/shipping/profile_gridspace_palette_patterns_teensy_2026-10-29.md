# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-29, **selective -O3**)

Palette-pattern recipe qualification at upper spec-admitted bounds (noise scale 8 baked; analytic frequency 64, chain interpreter only). These overrides are separate from the [normal GridSpace shipping profile](profile_gridspace_teensy_2026-10-23.md). They do not replace its ranked measurement.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz; live flywheel and DMA interrupts |
| Image | `profile`, shipping `-Os` base and selective `HS_O3` regions, including composed shade and column walking |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, 72 physical pixels |
| Effect | GridSpace at 288×144; single-entry playlist; runtime recipe override |
| Method | 70 s requested capture, 32-frame counter windows; display window 62.5 ms |
| Recipe | Enabled, amount 0.8; noise seed 6047, scale 8, speed 0.001; analytic frequency 64 on both axes, speed 0.02 |
| Noise variants | Single: SIMPLEX, contrast 0. Octave: RIDGED3, OCTAVE_CONTOUR, contrast 8 |
| Animation | `setAnimationsPaused(false)`; temporal runs advance both endpoint fields |
| Exception | Interference uses a temporary non-temporal Spec override; other runs are temporal |

The modulation selectors bypass effect admission and check only `Modulation::valid`. GridSpace admits baked noise only; direct noise and analytic captures predate that admission policy and do not demonstrate shipping support.

Captures use the locked `tools/profile_one.sh GridSpace profile 70 32` wrapper with `HS_PROFILE_TREE` and a pinned board. Source selectors: noise 0, sine 1, interference 2, grid 3. Mechanics: hue 0, lift 1, opacity 2. Evaluation: direct 0, baked cube 1. Supply `-D HS_PROFILE_MODULATION_SOURCE=N -D HS_PROFILE_MODULATION_MECHANIC=N`; noise defaults to baked and analytic fields to direct. Add `-D HS_PROFILE_MODULATION_EVALUATION=N` to override, and `-D HS_PROFILE_MODULATION_NOISE_OCTAVES` for the octave recipe. The non-temporal interference run additionally changes GridSpace’s Spec before building.

Compiler: Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 (20251203). These are working-tree captures: their base SHAs alone do not reconstruct the images. Per-run ELF hashes below identify the actual artifacts. The normal Phantasm image and the profile image were built separately; their RAM figures must not be interchanged.

### Image sizes

| Recipe | Profile FLASH code/data B | Profile RAM1 code/variables B | Phantasm FLASH code/data B | Phantasm RAM1 code/variables B |
|---|--:|--:|--:|--:|
| Single noise / hue / baked | 92792/157088 | 23736/313056 | 1008112/860664 | 163512/312736 |
| Directional sine / opacity / direct | 92792/157088 | 23736/313056 | 1008112/860664 | 163512/312736 |
| Smooth grid / lightness lift / direct | 92792/157088 | 23736/313056 | 1008112/860664 | 163512/312736 |
| Interference / hue / direct | 84592/156248 | 22968/313024 | 1004472/860664 | 163512/312736 |
| Octave noise / lightness lift / direct | 92936/157108 | 23736/313056 | 1008112/860664 | 163512/312736 |
| Single noise / lightness lift / direct | 92936/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |
| Octave noise / lightness lift / baked | 92952/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |
| Single noise / hue / direct | 92936/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |
| Octave noise / hue / baked | 92952/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |
| Octave noise / opacity / baked | 92952/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |
| Single noise / lightness lift / baked | 92936/157108 | 23752/313056 | 1009824/860684 | 163528/312736 |

RAM2 is 520,064 bytes of variables with 4,224 bytes free in these builds. These are separate per-recipe profile and full-Phantasm build sizes; code-placement and admission changed during the campaign.

## Frame cadence

All aggregates below use complete counter windows, including setup frame 1. Peaks are actual per-frame **render** telemetry with `canvas_buffer_wait` subtracted. A spill means render exceeds 62.5 ms. Trailing raw rows are excluded from aggregate denominators. Green means zero spills; any spill is red.

| Recipe | Peak render ms | Spilled | Frames | Startup ms | Post-startup peak ms | Post-startup spilled | Captured local |
|---|--:|--:|---|--:|--:|--:|---|
| Single noise / hue / baked | 🟢 58.221 | 🟢 0/1088 (0.00%) | 1–1088 | 58.221 | 46.488 | 0/1087 | 2026-10-29 11:42 |
| Directional sine / opacity / direct | 🟢 46.303 | 🟢 0/1088 (0.00%) | 1–1088 | 46.303 | 34.148 | 0/1087 | 2026-10-29 11:45 |
| Smooth grid / lightness lift / direct | 🟢 59.895 | 🟢 0/1088 (0.00%) | 1–1088 | 59.895 | 49.339 | 0/1087 | 2026-10-29 11:44 |
| Interference / hue / direct | 🟢 55.702 | 🟢 0/1088 (0.00%) | 1–1088 | 55.702 | 32.588 | 0/1087 | 2026-10-29 11:49 |
| Octave noise / lightness lift / direct | 🔴 109.018 | 🔴 352/352 (100.00%) | 1–352 | 109.018 | 100.422 | 351/351 | 2026-10-29 11:47 |
| Single noise / lightness lift / direct | 🔴 70.284 | 🔴 1/1088 (0.09%) | 1–1088 | 70.284 | 59.847 | 0/1087 | 2026-10-29 11:51 |
| Octave noise / lightness lift / baked | 🔴 76.126 | 🔴 271/544 (49.82%) | 1–544 | 76.126 | 65.470 | 270/543 | 2026-10-29 11:55 |
| Single noise / hue / direct | 🔴 66.274 | 🔴 1/1088 (0.09%) | 1–1088 | 66.274 | 54.302 | 0/1087 | 2026-10-29 11:58 |
| Octave noise / hue / baked | 🔴 64.235 | 🔴 1/1088 (0.09%) | 1–1088 | 64.235 | 52.066 | 0/1087 | 2026-10-29 12:01 |
| Octave noise / opacity / baked | 🟢 58.521 | 🟢 0/1088 (0.00%) | 1–1088 | 58.521 | 46.035 | 0/1087 | 2026-10-29 12:03 |
| Single noise / lightness lift / baked | 🔴 63.194 | 🔴 1/1088 (0.09%) | 1–1088 | 63.194 | 51.947 | 0/1087 | 2026-10-29 11:59 |

Times are raw-log filesystem mtimes in America/Los_Angeles (PDT), captured when this report was generated. The source-noise speed and analytic speed keep adjacent temporal fields distinct; these measurements include two completed endpoint evaluations rather than a static-field reuse case. Interference measures one endpoint only and does not establish the temporal interference cost.

There are two existing driver budgets. Steady frames have a 62.5 ms display window. Effect handoff has `commit_revs = 2` full revolutions, or 250 ms at 480 RPM, for release, teardown, construction and the first primed frame while the display is black (`hardware/pov_segmented.h`, `hardware/pov_sync_protocol.h`). Neither budget changed for this work. The table deliberately keeps setup-inclusive spill counts red when setup exceeds 62.5 ms; that does not by itself reject a recipe under the separate startup budget. Construction is not separately timed in these frame logs, and first-frame telemetry alone cannot prove the complete 250 ms handoff cost. The captures complete without the driver handoff-budget guard firing.

Direct configured single-noise hue and lift, baked configured single-noise lift, and baked octave hue have startup-only display-window spills; their measured later frames fit 62.5 ms. The direct octave/lift recipe spills throughout its capture, and baked octave/lift also has steady spills. The device policy excludes direct OCTAVE_CONTOUR noise for all mechanics and OCTAVE_CONTOUR plus lightness lift for either strategy. Other direct octave mechanics were not timed here; that exclusion is conservative. Passing rows qualify only the recorded combinations and do not validate every baked-noise combination or color mechanic.

### Maximum window means

These are inclusive complete-window means, not frame peaks; maxima can occur in different windows. Column values are total column preparation per frame (144 columns for temporal runs), not cost per column.

| Recipe | Prepare frame ms (frames) | Prepare column ms (frames) | Shader draw ms (frames) | Pass shader mean ms | Richest-window cycle/wall error ppm | Maximum error over all windows ppm |
|---|--:|--:|--:|--:|--:|--:|
| Single noise / hue / baked | 4.411 (353–384) | 2.234 (609–640) | 38.638 (1–32) | 38.242 | 0.420 | 1.778 |
| Directional sine / opacity / direct | 0.482 (353–384) | 2.377 (129–160) | 30.467 (1–32) | 30.074 | 0.167 | 7.366 |
| Smooth grid / lightness lift / direct | 0.485 (385–416) | 2.438 (1025–1056) | 45.260 (865–896) | 44.478 | 1.116 | 2.893 |
| Interference / hue / direct | 0.480 (353–384) | not instrumented | 26.259 (1–32) | 25.504 | 0.462 | 2.407 |
| Octave noise / lightness lift / direct | 0.461 (321–352) | 2.151 (321–352) | 96.553 (321–352) | 95.278 | 0.168 | 2.526 |
| Single noise / lightness lift / direct | 0.488 (385–416) | 2.167 (833–864) | 55.167 (865–896) | 54.472 | 0.109 | 2.013 |
| Octave noise / lightness lift / baked | 9.890 (385–416) | 2.131 (513–544) | 51.688 (385–416) | 51.045 | 0.643 | 1.066 |
| Single noise / hue / direct | 0.488 (353–384) | 2.134 (641–672) | 50.571 (1–32) | 50.170 | 0.711 | 3.267 |
| Octave noise / hue / baked | 9.894 (641–672) | 2.035 (257–288) | 39.066 (1–32) | 38.596 | 0.231 | 8.581 |
| Octave noise / opacity / baked | 9.894 (641–672) | 2.120 (353–384) | 33.135 (1–32) | 32.624 | 2.032 | 3.032 |
| Single noise / lightness lift / baked | 4.405 (353–384) | 2.130 (417–448) | 43.075 (1–32) | 42.698 | 1.288 | 3.034 |

`parse_profile.py ... validate --scope fx_shader_draw` passes for every capture: one GridSpace geometry, monotonic windows, complete per-frame records, and root cycles ÷ 600 MHz agreeing with the wall sum within 5 ppm in the richest shader window. The sine capture has one other window at 7.366 ppm, above that 5 ppm diagnostic threshold; the validator checks the richest shader window only. This is not a claim of sub-5-ppm agreement for every window. Raw capture coverage is:

| Raw log suffix | Complete windows | Trailing rows excluded | Raw peak ms | Raw spilled |
|---|--:|--:|--:|--:|
| `noise-hue-baked.log` | 34 | 8 | 58.221 | 0/1096 |
| `sine-opacity-direct.log` | 34 | 8 | 46.303 | 0/1096 |
| `grid-lift-direct.log` | 34 | 8 | 59.895 | 0/1096 |
| `interference-hue-direct.log` | 34 | 9 | 55.702 | 0/1097 |
| `noise-lift-direct.log` | 11 | 14 | 109.018 | 366/366 |
| `noise-single-lift-direct.log` | 34 | 8 | 70.284 | 1/1096 |
| `noise-octave-lift-baked.log` | 17 | 6 | 76.126 | 274/550 |
| `noise-hue-direct.log` | 34 | 8 | 66.274 | 1/1096 |
| `noise-octave-hue-baked.log` | 34 | 8 | 64.235 | 1/1096 |
| `noise-octave-opacity-baked.log` | 34 | 8 | 58.521 | 0/1096 |
| `noise-single-lift-baked.log` | 34 | 8 | 63.194 | 1/1096 |

All raw filenames have the prefix `palette-patterns-`; originals and provenance are local capture artifacts under `C:/work/temp`. No preset advance markers occur in these single-recipe runs.

## Phase-by-phase readout

Each recipe is a distinct capture regime. The tree below uses that recipe’s highest `fx_shader_draw` window, with inclusive per-frame averages and the nesting recorded by the harness. A tree can contain setup; the cadence table reports setup separately. Cycle counts are divided by the window frame count; percentages are shares of the inclusive root.

### Single noise / hue / baked

Window frames 1–32; wall min/mean/max 42.705/61.691/65.437 ms.

```text
frame                          61.691 ms   37.015 Mcyc  100.00%
  pov_preserve_half             0.132 ms    0.079 Mcyc    0.21%
  fx_shader_draw               38.638 ms   23.183 Mcyc   62.63%
    fx_prepare_column           2.199 ms    1.319 Mcyc    3.56%
  fx_prepare_frame              4.250 ms    2.550 Mcyc    6.89%
  fx_advance                    1.127 ms    0.676 Mcyc    1.83%
    palette_cycle_step          1.052 ms    0.631 Mcyc    1.71%
  fx_timeline_step              0.075 ms    0.045 Mcyc    0.12%
  canvas_clear                  0.092 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait           17.156 ms   10.293 Mcyc   27.81%
```

### Directional sine / opacity / direct

Window frames 1–32; wall min/mean/max 30.541/60.938/65.538 ms.

```text
frame                          60.938 ms   36.563 Mcyc  100.00%
  pov_preserve_half             0.138 ms    0.083 Mcyc    0.23%
  fx_shader_draw               30.467 ms   18.280 Mcyc   50.00%
    fx_prepare_column           2.307 ms    1.384 Mcyc    3.79%
  fx_prepare_frame              0.343 ms    0.206 Mcyc    0.56%
  fx_advance                    1.135 ms    0.681 Mcyc    1.86%
    palette_cycle_step          1.058 ms    0.635 Mcyc    1.74%
  fx_timeline_step              0.074 ms    0.044 Mcyc    0.12%
  canvas_clear                  0.092 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait           28.492 ms   17.095 Mcyc   46.76%
```

### Smooth grid / lightness lift / direct

Window frames 865–896; wall min/mean/max 58.403/62.420/66.457 ms.

```text
frame                          62.421 ms   37.452 Mcyc  100.00%
  pov_preserve_half             0.137 ms    0.082 Mcyc    0.22%
  fx_shader_draw               45.260 ms   27.156 Mcyc   72.51%
    fx_prepare_column           2.422 ms    1.453 Mcyc    3.88%
  fx_prepare_frame              0.399 ms    0.239 Mcyc    0.64%
  fx_advance                    1.101 ms    0.660 Mcyc    1.76%
    palette_cycle_step          1.066 ms    0.640 Mcyc    1.71%
  fx_timeline_step              0.070 ms    0.042 Mcyc    0.11%
  canvas_clear                  0.088 ms    0.053 Mcyc    0.14%
  canvas_buffer_wait           15.178 ms    9.107 Mcyc   24.32%
```

### Interference / hue / direct

Window frames 1–32; wall min/mean/max 25.890/61.079/68.776 ms.

```text
frame                          61.080 ms   36.648 Mcyc  100.00%
  pov_preserve_half             0.136 ms    0.082 Mcyc    0.22%
  fx_shader_draw               26.259 ms   15.755 Mcyc   42.99%
  fx_prepare_frame              0.338 ms    0.203 Mcyc    0.55%
  fx_advance                    2.788 ms    1.673 Mcyc    4.56%
    palette_cycle_step          2.594 ms    1.556 Mcyc    4.25%
  fx_timeline_step              0.057 ms    0.034 Mcyc    0.09%
  canvas_clear                  0.091 ms    0.054 Mcyc    0.15%
  canvas_buffer_wait           31.362 ms   18.817 Mcyc   51.35%
```

### Octave noise / lightness lift / direct

Window frames 321–352; wall min/mean/max 183.956/187.438/191.032 ms.

```text
frame                         187.438 ms  112.463 Mcyc  100.00%
  pov_preserve_half             0.141 ms    0.085 Mcyc    0.08%
  fx_shader_draw               96.553 ms   57.932 Mcyc   51.51%
    fx_prepare_column           2.151 ms    1.291 Mcyc    1.15%
  fx_prepare_frame              0.461 ms    0.277 Mcyc    0.25%
  fx_advance                    1.060 ms    0.636 Mcyc    0.57%
    palette_cycle_step          1.045 ms    0.627 Mcyc    0.56%
  fx_timeline_step              0.083 ms    0.050 Mcyc    0.04%
  canvas_clear                  0.086 ms    0.052 Mcyc    0.05%
  canvas_buffer_wait           88.862 ms   53.317 Mcyc   47.41%
```

### Single noise / lightness lift / direct

Window frames 865–896; wall min/mean/max 58.334/62.398/66.556 ms.

```text
frame                          62.398 ms   37.439 Mcyc  100.00%
  pov_preserve_half             0.135 ms    0.081 Mcyc    0.22%
  fx_shader_draw               55.167 ms   33.100 Mcyc   88.41%
    fx_prepare_column           2.163 ms    1.298 Mcyc    3.47%
  fx_prepare_frame              0.400 ms    0.240 Mcyc    0.64%
  fx_advance                    1.101 ms    0.660 Mcyc    1.76%
    palette_cycle_step          1.069 ms    0.641 Mcyc    1.71%
  fx_timeline_step              0.077 ms    0.047 Mcyc    0.12%
  canvas_clear                  0.088 ms    0.053 Mcyc    0.14%
  canvas_buffer_wait            5.238 ms    3.143 Mcyc    8.39%
```

### Octave noise / lightness lift / baked

Window frames 385–416; wall min/mean/max 65.048/124.908/184.551 ms.

```text
frame                         124.908 ms   74.945 Mcyc  100.00%
  pov_preserve_half             0.136 ms    0.081 Mcyc    0.11%
  fx_shader_draw               51.688 ms   31.013 Mcyc   41.38%
    fx_prepare_column           2.114 ms    1.268 Mcyc    1.69%
  fx_prepare_frame              9.890 ms    5.934 Mcyc    7.92%
  fx_advance                    1.059 ms    0.636 Mcyc    0.85%
    palette_cycle_step          1.041 ms    0.625 Mcyc    0.83%
  fx_timeline_step              0.087 ms    0.052 Mcyc    0.07%
  canvas_clear                  0.088 ms    0.053 Mcyc    0.07%
  canvas_buffer_wait           61.752 ms   37.051 Mcyc   49.44%
```

### Single noise / hue / direct

Window frames 1–32; wall min/mean/max 50.725/62.192/66.277 ms.

```text
frame                          62.192 ms   37.315 Mcyc  100.00%
  pov_preserve_half             0.136 ms    0.081 Mcyc    0.22%
  fx_shader_draw               50.571 ms   30.343 Mcyc   81.31%
    fx_prepare_column           2.102 ms    1.261 Mcyc    3.38%
  fx_prepare_frame              0.344 ms    0.206 Mcyc    0.55%
  fx_advance                    1.139 ms    0.684 Mcyc    1.83%
    palette_cycle_step          1.064 ms    0.638 Mcyc    1.71%
  fx_timeline_step              0.080 ms    0.048 Mcyc    0.13%
  canvas_clear                  0.091 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait            9.629 ms    5.777 Mcyc   15.48%
```

### Octave noise / hue / baked

Window frames 1–32; wall min/mean/max 48.665/62.055/65.024 ms.

```text
frame                          62.055 ms   37.233 Mcyc  100.00%
  pov_preserve_half             0.134 ms    0.081 Mcyc    0.22%
  fx_shader_draw               39.066 ms   23.439 Mcyc   62.95%
    fx_prepare_column           2.012 ms    1.207 Mcyc    3.24%
  fx_prepare_frame              9.737 ms    5.842 Mcyc   15.69%
  fx_advance                    1.071 ms    0.643 Mcyc    1.73%
    palette_cycle_step          0.996 ms    0.598 Mcyc    1.61%
  fx_timeline_step              0.082 ms    0.049 Mcyc    0.13%
  canvas_clear                  0.092 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait           11.652 ms    6.991 Mcyc   18.78%
```

### Octave noise / opacity / baked

Window frames 1–32; wall min/mean/max 42.722/61.690/65.284 ms.

```text
frame                          61.690 ms   37.014 Mcyc  100.00%
  pov_preserve_half             0.135 ms    0.081 Mcyc    0.22%
  fx_shader_draw               33.135 ms   19.881 Mcyc   53.71%
    fx_prepare_column           2.098 ms    1.259 Mcyc    3.40%
  fx_prepare_frame              9.739 ms    5.843 Mcyc   15.79%
  fx_advance                    1.050 ms    0.630 Mcyc    1.70%
    palette_cycle_step          0.981 ms    0.589 Mcyc    1.59%
  fx_timeline_step              0.082 ms    0.049 Mcyc    0.13%
  canvas_clear                  0.091 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait           17.240 ms   10.344 Mcyc   27.95%
```

### Single noise / lightness lift / baked

Window frames 1–32; wall min/mean/max 47.381/61.997/66.611 ms.

```text
frame                          61.998 ms   37.199 Mcyc  100.00%
  pov_preserve_half             0.135 ms    0.081 Mcyc    0.22%
  fx_shader_draw               43.075 ms   25.845 Mcyc   69.48%
    fx_prepare_column           2.110 ms    1.266 Mcyc    3.40%
  fx_prepare_frame              4.255 ms    2.553 Mcyc    6.86%
  fx_advance                    1.144 ms    0.686 Mcyc    1.84%
    palette_cycle_step          1.072 ms    0.643 Mcyc    1.73%
  fx_timeline_step              0.077 ms    0.046 Mcyc    0.12%
  canvas_clear                  0.092 ms    0.055 Mcyc    0.15%
  canvas_buffer_wait           12.999 ms    7.800 Mcyc   20.97%
```

### Per-pixel figures

Blend-count instrumentation is absent. A nominal quadrant contains 10,368 logical pixels; that geometry is not a measured blend count. No cycles-per-blended-pixel estimate is inferred from it.

## Column-ISR / DMA marshaling cost

The following ISR rows correspond to the selected recipe trees above. The ISR interval has its own wall timestamp; counts divided by 32 are calls per rendered frame over that interval, not a fixed hardware rate.

| Recipe | ISR scope | Calls/frame | Min/mean/max µs | CPU share |
|---|---|--:|--:|--:|
| noise-hue-baked | isr_wake | 1392.09 | 0.488/1.693/33.662 | 3.12% |
| noise-hue-baked | isr_pack | 138.09 | 7.308/9.142/11.622 | 1.67% |
| noise-hue-baked | isr_dma_submit | 138.09 | 0.650/0.930/5.063 | 0.17% |
| sine-opacity-direct | isr_wake | 1385.19 | 0.488/1.678/26.333 | 3.09% |
| sine-opacity-direct | isr_pack | 137.25 | 6.825/9.027/11.980 | 1.64% |
| sine-opacity-direct | isr_dma_submit | 137.25 | 0.605/0.928/1.758 | 0.16% |
| grid-lift-direct | isr_wake | 1151.94 | 0.523/2.008/25.888 | 3.70% |
| grid-lift-direct | isr_pack | 143.97 | 8.455/9.323/11.997 | 2.14% |
| grid-lift-direct | isr_dma_submit | 143.97 | 0.627/0.930/5.148 | 0.21% |
| interference-hue-direct | isr_wake | 1382.41 | 0.470/1.663/21.992 | 3.06% |
| interference-hue-direct | isr_pack | 136.91 | 7.180/9.083/11.967 | 1.65% |
| interference-hue-direct | isr_dma_submit | 136.91 | 0.665/0.932/1.248 | 0.17% |
| noise-lift-direct | isr_wake | 3456.53 | 0.475/1.978/23.395 | 3.64% |
| noise-lift-direct | isr_pack | 432.03 | 8.515/9.473/13.200 | 2.18% |
| noise-lift-direct | isr_dma_submit | 432.03 | 0.602/0.927/2.858 | 0.21% |
| noise-single-lift-direct | isr_wake | 1151.66 | 0.543/2.008/32.570 | 3.70% |
| noise-single-lift-direct | isr_pack | 143.94 | 8.415/9.370/11.907 | 2.15% |
| noise-single-lift-direct | isr_dma_submit | 143.94 | 0.608/0.940/1.222 | 0.21% |
| noise-octave-lift-baked | isr_wake | 2304.00 | 0.547/1.998/21.583 | 3.68% |
| noise-octave-lift-baked | isr_pack | 287.97 | 8.692/9.522/13.467 | 2.19% |
| noise-octave-lift-baked | isr_dma_submit | 287.97 | 0.662/0.940/1.992 | 0.21% |
| noise-hue-direct | isr_wake | 1396.75 | 0.487/1.690/26.817 | 3.11% |
| noise-hue-direct | isr_pack | 138.69 | 7.238/9.170/12.152 | 1.67% |
| noise-hue-direct | isr_dma_submit | 138.69 | 0.727/0.942/7.798 | 0.17% |
| noise-octave-hue-baked | isr_wake | 1395.53 | 0.488/1.687/29.597 | 3.10% |
| noise-octave-hue-baked | isr_pack | 138.53 | 7.097/9.152/11.635 | 1.67% |
| noise-octave-hue-baked | isr_dma_submit | 138.53 | 0.670/0.942/6.420 | 0.17% |
| noise-octave-opacity-baked | isr_wake | 1392.16 | 0.488/1.680/28.697 | 3.09% |
| noise-octave-opacity-baked | isr_pack | 138.12 | 7.230/9.125/12.282 | 1.66% |
| noise-octave-opacity-baked | isr_dma_submit | 138.12 | 0.597/0.940/2.888 | 0.17% |
| noise-single-lift-baked | isr_wake | 1394.97 | 0.487/1.700/32.077 | 3.13% |
| noise-single-lift-baked | isr_pack | 138.47 | 7.400/9.257/11.835 | 1.69% |
| noise-single-lift-baked | isr_dma_submit | 138.47 | 0.653/0.942/5.530 | 0.17% |

`isr_wake` includes `isr_pack` and `isr_dma_submit`; their shares overlap and must not be added. Submit measures CPU execution, not asynchronous LED wire-transfer duration. DMA completion and unrelated interrupts are not separately measured. Cycle scopes absorb interrupt time because CYCCNT free-runs. ISR CPU share cannot be subtracted again from these measured render costs.

## Summary ranking

`fx_shader_draw` dominates every recorded regime. Direct octave/lift reaches roughly 97 ms per frame in its worst shader window, while the passing analytic recipes have substantially lower measured shader costs. Baked octave/lift adds roughly 10 ms per-frame preparation at its worst window and still spills. These source/mechanic costs differ enough that a source-only timing claim would hide admission failures.

The setup-inclusive red counts remain visible even when later frames meet the steady 62.5 ms budget. Startup admission uses the existing separate 250 ms construction-and-first-frame budget; it is not inferred from the color of the aggregate row. These are GridSpace qualification observations, not a roster-wide performance guarantee and not a replacement for the normal shipping ranking.

## Caveats

- Firmware admission evolved during the campaign. Rejected recipes were measured before or with experimental admission; their captures are evidence for rejection, not proof that the final runtime admits them.
- Different recipes were built at different source snapshots and on two boards. They are not controlled same-image comparisons.
- Amount 0.8, upper spatial bounds and moving clocks stress these particular combinations; inactive recipes, zero amounts and static shared fields can have different costs.
- No full preset cycle is involved: GridSpace has one authored preset and the harness applies one recipe override per capture.
- Raw logs, envdumps, source-copy manifests, build logs and ELF artifacts are local outputs. Base commit IDs identify ancestry but do not make dirty capture inputs reproducible from a fresh clone.

This direct-lightness recapture also bypassed effect admission and predates the baked-noise-only GridSpace policy; it does not qualify direct noise for shipping.

The final temporal-column storage implementation was recaptured on COM4 using
host `0448cff3c039f64c8dba5253266e17df065c2083`: configured-single noise,
DIRECT lightness, unchanged controls and 70 seconds. All 34 complete windows
validate: startup is 70.385 ms, subsequent render peaks at 59.831 ms, and
only frame 1 spills (1/1,088 total). Maximum window means are 479.165 us/frame
for preparation, 2,173.241 us/frame for column preparation (15.092 us/column),
and 55.278 ms/frame for shader drawing. Pattern and color kernel files are
byte-identical to the campaign sources.

The final copy manifest is `palette-final-resource-profile-copy.json`, SHA-256
`e22ba2b2fd5be7a5e32db84fe178317d9e6a6ebc96bdaf269789452e326bd926`.
The raw `palette-patterns-final-noise-lift-direct.log` SHA-256 is
`db8082e7f5c1a94bae7a2fb034e7feb6ada0a6a364c98d6634f3767b8dd44602`;
its profile ELF SHA-256 is
`5425963364ee6e32d7a1090f052a047c3aec9c77993721ec775ccf96a9abfd4f`.
This measurement includes the extra retained column context that reduces
strict-FP WebAssembly stack use; it does not raise either driver deadline.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`, shipping selective-O3 profile environment. The capture wrapper enforces the shared-board lock and verifies image identity. No hardware was accessed while generating this report.

### Capture provenance

| Recipe | Base SHA recorded by wrapper | Profile ELF SHA-256 |
|---|---|---|
| noise-hue-baked | `3a9de6a7cd1b26e302e33b388e38b02cca3df997` | `dc2f382f315e957abb7a2fc48ea1596bcd7674d0888a1ca8a1c5870011bb8ae7` |
| sine-opacity-direct | `3a9de6a7cd1b26e302e33b388e38b02cca3df997` | `42c1a159939b0bdbc72e01f3e911c464c43fa2da75a90987a7433ad845ff576d` |
| grid-lift-direct | `1a56116631a2955909d76b464c9f8bbe2ebe27ef` | `15aade042e4b9db1dfcebd819c4a7fc326e08f6512e60cdb579d2a719fbbcdbd` |
| interference-hue-direct | `1a56116631a2955909d76b464c9f8bbe2ebe27ef` | `4ee02acc117a25611493435750d3495558a4d5553b1bd40dc4da1445c7e6b98b` |
| noise-lift-direct | `6dc6b73c12fb5f105f5f90baf2166c17e77d0e0d` | `a8042f93d8ba3f7da486a8eaec13b1335b0507f8c0f734bfb495e3dec502887f` |
| noise-single-lift-direct | `8c5ec46349b041af0da44fbf9fbf4e3ae8be9a41` | `b77979103f000976b7f6db3f8f5a19a9800bcef51365e2a8464023cb73d8be45` |
| noise-octave-lift-baked | `eca907daf87f6517d1aa23f268849f7bb586a711` | `77b3a2e014a02604dee3549f464a0143e9f1ea5de6d926c62f84eadd594ad4bc` |
| noise-hue-direct | `1a56116631a2955909d76b464c9f8bbe2ebe27ef` | `445f17163db00d4ae5cc47e55e8732f30f8526f7a92c0146ee632a131921373f` |
| noise-octave-hue-baked | `1a56116631a2955909d76b464c9f8bbe2ebe27ef` | `610c79c49c61a15a4fcdc945621e109fed76b8fbd24e556e9ec8f55cc6d86840` |
| noise-octave-opacity-baked | `1a56116631a2955909d76b464c9f8bbe2ebe27ef` | `b63ae505481889f069049e966c78f13b14e6b72431e5171092d7749a36dd98b9` |
| noise-single-lift-baked | `b0d9e69bd22364dd835447eb85b0947ea416cec9` | `a8090e96ae1c508aae21bf6e6c02193c20c2ed9347187784f157a93e942a0b5d` |

The noise-hue and sine captures used the host tree at `3a9de6a7c` with uncommitted integration inputs on COM3. Grid and interference used the sources worker at `1a5611663` on COM4, populated from the host’s `3a9de6a7c` dirty inputs. `palette-resource-profile-copy.json` records the copied paths and hashes. `palette-profile-temporal-override.json` records interference’s GridSpace change from temporal to non-temporal; its before/after SHA-256 values are `24249c0c5031ea8f5419bacbb831575822488eb7fa86331ed7bc5d9bc7877f65` and `71635b5fd81d62d5b7cfcd8950a996ae6e2dbb7c85c53423540bdd3ac5cb4b36`. Later noise captures used the host tree on COM3 at the per-row bases above, with ongoing integration inputs.

Additional noise/hue and noise/opacity captures use the sources worker on COM4, populated from host `eca907daf87f6517d1aa23f268849f7bb586a711` plus its dirty spec. `palette-noise-final-profile-copy.json` (SHA-256 `0081e5547d4f6bc8d89b727833bf326da5c68963137839b5ec40fba51721a389`) records those exact copies; no temporal override is applied.

### Raw-log identity

| Recipe | Raw-log SHA-256 |
|---|---|
| noise-hue-baked | `dfe39ac9336cecb8278a8f40a4aa21a9917661a7da6cadcba5201e1f64a9bfa6` |
| sine-opacity-direct | `41338d50aabc9d15ab39e562c1b72a1782c00789e0202f6632814ba485f3ac16` |
| grid-lift-direct | `5c870ff5000916b4d775018203452ebdcb3f7bfbc0ae7e918463b97e6a074d36` |
| interference-hue-direct | `a621151dfb523ee8ad6362bb1bf5d29b50c5455aa8322f9967cc4a77126d71dd` |
| noise-lift-direct | `1298f5c51dcbef2ea6d15922f5316e5741107c3e4889abec86d4db14b2121e23` |
| noise-single-lift-direct | `451bf08c956cf7de44a4fa72952dc2e2ff5c40c11ae974601ecc8968e9243ee4` |
| noise-octave-lift-baked | `8a00177fa1e468a6d65fea5627d76968d7408abedffef97a119959c5cbfc9391` |
| noise-hue-direct | `1b0472805e7a958260eb1d0d7de270ae6c58ffa1ad169f14fb756020f7e2c7f4` |
| noise-octave-hue-baked | `b454005a280814b975d83dca0bbd5bb23b567c5682b2c13216d9a20b49dab1e9` |
| noise-octave-opacity-baked | `0ea1eef64c5d310970f1a21413479d10ceb5bae19d8f13ef0704ad470e56aaea` |
| noise-single-lift-baked | `61dfae2d87c3c5cda87f8baefef7d74ce742d84f1dd19508b32171781427ef95` |
