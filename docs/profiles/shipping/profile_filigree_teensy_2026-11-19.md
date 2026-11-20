# Filigree on-device profile — Teensy 4.0, segmented mode (2026-11-19, **selective -O3**)

Raw capture: `C:/work/temp/holosphere-20261119-filigree-profile/optimized/filigree_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` scan, spherical source, and pullback color paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Filigree, 288×144, source `47bceb75c37852a8ea5278fc27c72d65ff453d20` |
| Method | 108 s capture, 16 frames per readout, default one-hour profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Filigree profile 108 16` |

Profile image size: `FLASH: code:113016, data:349292, headers:8732; free for files:1560576 / RAM1: variables:314144, code:34680, padding:30856; free for local variables:144608 / RAM2: variables:520064; free for malloc/new:4224`.

Shipping Phantasm image size: `FLASH: code:1110896, data:617824, headers:8984; free for files:293912 / RAM1: variables:313760, code:195400, padding:1208; free for local variables:13920 / RAM2: variables:520064; free for malloc/new:4224`.

Compared with the separate-effect Phantasm build, FLASH code changed by -5,608 B and FLASH data by -188 B. RAM1 code changed by +656 B and variables by +0 B.

Captured 2026-11-19 22:42, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `293b2f98921a3f2bc4cdcb6d496d897416dfadab2936245ce87864d963177db4`. Profile ELF SHA-256 `76632e24f185e71a93c2c49a0199d8a33aef061ecbc1839336f9a9ee260d602f`. Phantasm ELF SHA-256 `846952e49df0a339c1664a21097114b01ee51f2ae0ad870bd0486fa79456164b`. The capture retains attested ELF, map and toolchain configuration paths.

`parse_profile.py ... validate` passed. In frames 1073–1088, root cycles divided by 600 MHz match the foreground wall sum to 0.85 ppm.

## Frame cadence

The Phantasm playlist gives Filigree 105 s. The measurement includes 1,672 frames confirmed complete within the slot, including setup. The last included frame completed by 104.980 s. Analysis excludes 32 frame records beyond the cutoff, including up to 2 frames with an uncertain completion time at the boundary.

Setup frame 1 took **110.198 ms**. Frames 2–1672 averaged **56.458 ms** render and peaked at **71.527 ms** on frame 1080 (Ember Lattice). **1/1671 (0.06%)** exceeded the 62.5 ms display budget.

README cells: peak 🔴 71.527 (1)<br>🟢 57.614 (1), spilled 🔴 1/1073 (0.09%)<br>🟢 0/598 (0.00%).

README setup frame 1: 110.198 ms.

## Phase-by-phase readout

The authored show holds each preset for 600 frames and uses 480-frame sinusoidally eased transitions. Preset markers attribute a transition to the incoming preset. Counter trees give inclusive averages per frame. Whole in-show readouts provide the window statistics.

### Startup (frames 1–16)

```
frame                         58.943 ms   35.366 Mcyc 100.0% x1.00
  pov_preserve_half            0.124 ms    0.075 Mcyc   0.2% x0.94
  fx_shader_draw              57.611 ms   34.567 Mcyc  97.7% x1.00
    fx_prepare_column          1.561 ms    0.937 Mcyc   2.6% x144.00
  fx_palette_base_bake         0.729 ms    0.437 Mcyc   1.2% x1.00
    fx_palette_rgb_bake        0.647 ms    0.388 Mcyc   1.1% x1.00
  fx_prepare_frame             0.032 ms    0.019 Mcyc   0.1% x1.00
  fx_advance                   0.189 ms    0.113 Mcyc   0.3% x1.00
  fx_timeline_step             0.004 ms    0.003 Mcyc   0.0% x1.00
  canvas_clear                 0.098 ms    0.059 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 52.285/58.942/110.200 ms. Mean render was 58.941 ms and the shader draw averaged 57.611 ms.

### Transition into Ember Lattice (frames 737–752)

```
frame                         58.530 ms   35.118 Mcyc 100.0% x1.00
  pov_preserve_half            0.133 ms    0.080 Mcyc   0.2% x1.00
  fx_shader_draw              57.835 ms   34.701 Mcyc  98.8% x1.00
    fx_prepare_column          2.122 ms    1.273 Mcyc   3.6% x144.00
  fx_palette_base_bake         0.082 ms    0.049 Mcyc   0.1% x1.00
  fx_prepare_frame             0.022 ms    0.013 Mcyc   0.0% x1.00
  fx_advance                   0.054 ms    0.033 Mcyc   0.1% x1.00
  fx_timeline_step             0.134 ms    0.081 Mcyc   0.2% x1.00
  canvas_clear                 0.089 ms    0.054 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 58.062/58.529/58.897 ms. Mean render was 58.528 ms and the shader draw averaged 57.835 ms.

### Peak live render (frames 1073–1088)

```
frame                         57.164 ms   34.298 Mcyc 100.0% x1.00
  pov_preserve_half            0.133 ms    0.080 Mcyc   0.2% x1.00
  fx_shader_draw              56.359 ms   33.815 Mcyc  98.6% x1.00
    fx_prepare_column          1.724 ms    1.035 Mcyc   3.0% x144.00
  fx_palette_base_bake         0.075 ms    0.045 Mcyc   0.1% x0.94
  fx_prepare_frame             0.212 ms    0.127 Mcyc   0.4% x1.00
  fx_advance                   0.057 ms    0.034 Mcyc   0.1% x1.00
  fx_timeline_step             0.069 ms    0.041 Mcyc   0.1% x1.00
  canvas_clear                 0.091 ms    0.055 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 55.359/57.163/71.528 ms. Mean render was 57.162 ms and the shader draw averaged 56.359 ms.

### Show end (frames 1649–1664)

```
frame                         56.578 ms   33.947 Mcyc 100.0% x1.00
  pov_preserve_half            0.133 ms    0.080 Mcyc   0.2% x1.00
  fx_shader_draw              56.034 ms   33.620 Mcyc  99.0% x1.00
    fx_prepare_column          1.322 ms    0.793 Mcyc   2.3% x144.00
  fx_palette_base_bake         0.076 ms    0.045 Mcyc   0.1% x1.00
  fx_prepare_frame             0.023 ms    0.014 Mcyc   0.0% x1.00
  fx_advance                   0.056 ms    0.034 Mcyc   0.1% x1.00
  fx_timeline_step             0.003 ms    0.002 Mcyc   0.0% x1.00
  canvas_clear                 0.089 ms    0.053 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 55.537/56.577/58.412 ms. Mean render was 56.576 ms and the shader draw averaged 56.034 ms.

### Per-preset table

Hold figures use the costliest whole readout by shader draw within the authored dwell. Live peaks and spill counts include that preset’s incoming transition.

| Preset | Hold windows | Representative frames | Draw ms/f | Render ms/f | Live peak ms | Live spilled |
|---|--:|---|--:|--:|--:|--:|
| 1 · Jade Weave | 36 | 97–112 | 56.049 | 57.189 | 57.614 | 0/598 (0.00%) |
| 2 · Ember Lattice | 36 | 1153–1168 | 56.833 | 57.373 | 71.527 | 1/1073 (0.09%) |

### Per-pixel figures

The clipped quadrant contains 10,368 output pixels. This capture has no blended-pixel counter.

## Column-ISR / DMA marshaling cost

Peak window frames 1073–1088. CPU percentages use the separate 1,001,299 μs ISR snapshot interval. Times are min/avg/max.

```
isr_wake          1153.6/f 0.500/1.625/20.851 us 2.99%
isr_pack           144.2/f 5.060/6.388/10.638 us 1.47%
isr_dma_submit     144.2/f 0.656/0.920/6.463 us 0.21%
```

CPU packing averages 6.388 μs per call and CPU DMA submission averages 0.920 μs per call. At 24 MHz SPI the 600-byte composite has an asynchronous transfer budget of 230 μs including framing gaps. Pack and submit run inside the inclusive wake scope, which uses 2.99% of the interval.

## Summary ranking

1. `fx_shader_draw` — 98.6% of the foreground frame, 56.359 ms/f.
2. `fx_prepare_column` — 3.0% of the foreground frame, 1.724 ms/f.
3. `fx_prepare_frame` — 0.4% of the foreground frame, 0.212 ms/f.

## Edge sampler comparison

At source `480f3cdc95879b44edc7ff2cefee62b56d56b315`, the same Filigree preset
bank peaked at 75.007 ms live with 32 of 1,671 live frames over budget. The
measured build compares candidate distances by squared sine before converting
the selected distance to radians. Its compiled edge sampler runs in ITCM.

A native comparison covered 987,400 samples across all three baked solids.
Palette coordinates matched on visible samples, and maximum absolute coverage
difference was 0.000000358.

## Caveats

- Foreground cycle counts include interrupt time. Render telemetry subtracts the `canvas_buffer_wait` scope.
- The cutoff uses the first readout’s elapsed time, subsequent ISR snapshot intervals, and individual draw durations in the boundary readout. Per-frame render measurements come directly from the log. The full-playlist handoff was outside this measurement.
- Source and build provenance identify the clean measured checkout. The shipping preset order, dwell and easing were active. No timing overrides or deep instrumentation were used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Filigree`, `HS_PROFILE_WINDOW=16`. Acquire the board through `tools/profile_one.sh`.
