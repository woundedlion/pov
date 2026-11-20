# EmberTracery on-device profile — Teensy 4.0, segmented mode (2026-11-19, **selective -O3**)

Raw capture: `C:/work/temp/holosphere-20261119-tiling-profiles/column-time/embertracery_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` scan, spherical source, and pullback color paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Ember Tracery, 288×144, source `a0ae093c47a3852317f802d4e59366c9bc4885f1` |
| Method | 176 s capture, 16 frames per readout, default one-hour profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh EmberTracery profile 176 16` |

Profile image size: `FLASH: code:113424, data:349296, headers:8320; free for files:1560576 / RAM1: variables:314144, code:33832, padding:31704; free for local variables:144608 / RAM2: variables:520064; free for malloc/new:4224`.

Shipping Phantasm image size: `FLASH: code:1119248, data:618036, headers:8628; free for files:285704 / RAM1: variables:313760, code:194744, padding:1864; free for local variables:13920 / RAM2: variables:520064; free for malloc/new:4224`.

Captured 2026-11-19 19:59, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `8cc2a15fe9880693152af45c912c0d752309d98b54767ef39b354ce883352f2f`. Profile ELF SHA-256 `16e27b1af7fd1536bad121462beb73158b0477de4df7abd2ce06a1b50e37f1e1`. Phantasm ELF SHA-256 `30aeea9589b8dda81022339b91911802f3f8bd3fe8e2a79c056cab545d3b655e`. The capture retains attested ELF, map and toolchain configuration paths.

The trace begins with Golden Mesh at frame 1, advances to Jade Weave at frame 600 and reaches Ember Lattice at frame 1679. The raw validator expects a wrap marker, so preset coverage was checked against the initial state and both advance markers. Its effect, frame-continuity, reset, and cycle-counter checks passed. In frames 1777–1792, root cycles divided by 600 MHz match the foreground wall sum to 0.00 ppm.

## Frame cadence

The Phantasm playlist gives Ember Tracery 173 s. The measurement includes 2,659 frames confirmed complete within the slot, including setup. The last included frame completed by 173.000 s. Analysis excludes 32 frame records beyond the cutoff, including up to 0 frames with an uncertain completion time at the boundary.

Setup frame 1 took **46.864 ms**. Frames 2–2659 averaged **51.992 ms** render and peaked at **82.277 ms** on frame 1792 (Ember Lattice). **485/2658 (18.25%)** exceeded the 62.5 ms display budget.

README cells: peak 🔴 82.277 (1)<br>🟢 54.891 (2), spilled 🔴 485/981 (49.44%)<br>🟢 0/1677 (0.00%).

README setup frame 1: 46.864 ms.

## Phase-by-phase readout

The authored show holds each preset for 600 frames and uses 480-frame sinusoidally eased transitions. Preset markers attribute a transition to the incoming preset. Counter trees give inclusive averages per frame. Whole in-show readouts provide the window statistics.

### Startup (frames 1–16)

```
frame                         24.654 ms   14.792 Mcyc 100.0% x1.00
  pov_preserve_half            0.130 ms    0.078 Mcyc   0.5% x0.94
  fx_shader_draw              23.063 ms   13.838 Mcyc  93.5% x1.00
    fx_prepare_column          1.554 ms    0.932 Mcyc   6.3% x144.00
  fx_palette_base_bake         0.694 ms    0.416 Mcyc   2.8% x1.00
    fx_palette_rgb_bake        0.604 ms    0.362 Mcyc   2.4% x1.00
  fx_prepare_frame             0.333 ms    0.200 Mcyc   1.4% x1.00
  fx_advance                   0.070 ms    0.042 Mcyc   0.3% x1.00
  fx_timeline_step             0.004 ms    0.002 Mcyc   0.0% x1.00
  canvas_clear                 0.095 ms    0.057 Mcyc   0.4% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 23.055/24.654/46.866 ms. Mean render was 24.652 ms and the shader draw averaged 23.063 ms.

### Transition into Jade Weave (frames 1057–1072)

```
frame                         51.336 ms   30.802 Mcyc 100.0% x1.00
  pov_preserve_half            0.134 ms    0.080 Mcyc   0.3% x1.00
  fx_shader_draw              48.194 ms   28.917 Mcyc  93.9% x1.00
    fx_prepare_column          3.499 ms    2.100 Mcyc   6.8% x144.00
  fx_prepare_frame             2.534 ms    1.520 Mcyc   4.9% x1.00
  fx_advance                   0.048 ms    0.029 Mcyc   0.1% x1.00
  fx_timeline_step             0.132 ms    0.079 Mcyc   0.3% x1.00
  canvas_clear                 0.090 ms    0.054 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 51.131/51.336/51.554 ms. Mean render was 51.335 ms and the shader draw averaged 48.194 ms.

### Transition into Ember Lattice (frames 1777–1792)

```
frame                         80.142 ms   48.085 Mcyc 100.0% x1.00
  pov_preserve_half            0.137 ms    0.082 Mcyc   0.2% x1.00
  fx_shader_draw              76.348 ms   45.809 Mcyc  95.3% x1.00
    fx_prepare_column          3.466 ms    2.080 Mcyc   4.3% x144.00
  fx_palette_base_bake         0.680 ms    0.408 Mcyc   0.8% x1.00
    fx_palette_rgb_bake        0.595 ms    0.357 Mcyc   0.7% x1.00
  fx_prepare_frame             2.447 ms    1.468 Mcyc   3.1% x1.00
  fx_advance                   0.068 ms    0.041 Mcyc   0.1% x1.00
  fx_timeline_step             0.131 ms    0.078 Mcyc   0.2% x1.00
  canvas_clear                 0.097 ms    0.058 Mcyc   0.1% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 78.533/80.141/82.279 ms. Mean render was 80.139 ms and the shader draw averaged 76.348 ms.

### Show end (frames 2641–2656)

```
frame                         60.242 ms   36.145 Mcyc 100.0% x1.00
  pov_preserve_half            0.134 ms    0.080 Mcyc   0.2% x1.00
  fx_shader_draw              59.605 ms   35.763 Mcyc  98.9% x1.00
    fx_prepare_column          2.406 ms    1.444 Mcyc   4.0% x144.00
  fx_palette_base_bake         0.080 ms    0.048 Mcyc   0.1% x1.00
  fx_prepare_frame             0.026 ms    0.015 Mcyc   0.0% x1.00
  fx_advance                   0.078 ms    0.047 Mcyc   0.1% x1.00
  fx_timeline_step             0.002 ms    0.001 Mcyc   0.0% x1.00
  canvas_clear                 0.104 ms    0.063 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.002 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 59.418/60.241/61.067 ms. Mean render was 60.240 ms and the shader draw averaged 59.605 ms.

### Per-preset table

Hold figures use the costliest whole readout by shader draw within the authored dwell. Live peaks and spill counts include that preset’s incoming transition.

| Preset | Hold windows | Representative frames | Draw ms/f | Render ms/f | Live peak ms | Live spilled |
|---|--:|---|--:|--:|--:|--:|
| 1 · Golden Mesh | 36 | 321–336 | 21.934 | 23.177 | 23.351 | 0/598 (0.00%) |
| 2 · Jade Weave | 36 | 1121–1136 | 53.174 | 54.452 | 54.891 | 0/1079 (0.00%) |
| 3 · Ember Lattice | 31 | 2321–2336 | 60.420 | 61.048 | 82.277 | 485/981 (49.44%) |

### Per-pixel figures

The clipped quadrant contains 10,368 output pixels. This capture has no blended-pixel counter.

## Column-ISR / DMA marshaling cost

Peak window frames 1777–1792. CPU percentages use the separate 1,286,701 μs ISR snapshot interval. Times are min/avg/max.

```
isr_wake          1482.4/f 0.495/1.621/24.586 us 2.98%
isr_pack           185.2/f 5.361/6.343/9.736 us 1.46%
isr_dma_submit     185.2/f 0.748/0.923/8.605 us 0.21%
```

CPU packing averages 6.343 μs per call and CPU DMA submission averages 0.923 μs per call. At 24 MHz SPI the 600-byte composite has an asynchronous transfer budget of 230 μs including framing gaps. Pack and submit run inside the inclusive wake scope, which uses 2.98% of the interval.

## Summary ranking

1. `fx_shader_draw` — 95.3% of the foreground frame, 76.348 ms/f.
2. `fx_prepare_column` — 4.3% of the foreground frame, 3.466 ms/f.
3. `fx_prepare_frame` — 3.1% of the foreground frame, 2.447 ms/f.

## Caveats

- Foreground cycle counts include interrupt time. Render telemetry subtracts the `canvas_buffer_wait` scope.
- The cutoff uses the first readout’s elapsed time, subsequent ISR snapshot intervals, and individual draw durations in the boundary readout. Per-frame render measurements come directly from the log. The full-playlist handoff was outside this measurement.
- Source and build provenance identify the clean measured checkout. The shipping preset order, dwell and easing were active. No timing overrides or deep instrumentation were used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=EmberTracery`, `HS_PROFILE_WINDOW=16`. Acquire the board through `tools/profile_one.sh`.
