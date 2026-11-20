# EmberTracery on-device profile — Teensy 4.0, segmented mode (2026-11-19, **selective -O3**)

Raw capture: `C:/work/temp/holosphere-20261119-tiling-profiles/single-preset/embertracery_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` scan, spherical source, and pullback color paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Ember Tracery, 288×144, source `eb05f411a2d479d566c270fa3ced8168d41e30c4` |
| Method | 42 s capture, 32 frames per readout, default one-hour profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh EmberTracery profile 42 32` |

Profile image size: `FLASH: code:111056, data:348932, headers:9004; free for files:1562624 / RAM1: variables:314144, code:33688, padding:31848; free for local variables:144608 / RAM2: variables:520064; free for malloc/new:4224`.

Shipping Phantasm image size: `FLASH: code:1117568, data:618024, headers:8256; free for files:287768 / RAM1: variables:313760, code:194744, padding:1864; free for local variables:13920 / RAM2: variables:520064; free for malloc/new:4224`.

Captured 2026-11-19 20:07, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `21a00f5d048ef7589a47a010a2241cf0de834ff5535d7bc682f5c6a0ec582f23`. Profile ELF SHA-256 `75857e543cddcb88dd76d4c79bc65dd65527f92f6f4bdff5e0ef9bf0b0d93652`. Phantasm ELF SHA-256 `909eeb65a6b7eff58dd79e2af1dbec9c5215f63ebc80744e1b30261f7e5ea7e4`. The capture retains attested ELF, map and toolchain configuration paths.

`parse_profile.py ... validate` passed. In frames 513–544, root cycles divided by 600 MHz match the foreground wall sum to 0.77 ppm.

## Frame cadence

The Phantasm playlist gives Ember Tracery 38 s. The measurement includes 600 frames confirmed complete within the slot, including setup. The last included frame completed by 37.964 s. Analysis excludes 48 frame records beyond the cutoff, including up to 2 frames with an uncertain completion time at the boundary.

Setup frame 1 took **115.671 ms**. Frames 2–600 averaged **58.462 ms** render and peaked at **60.930 ms** on frame 517. **0/599 (0.00%)** exceeded the 62.5 ms display budget.

README cells: peak 🟢 60.930, spilled 🟢 0/599 (0.00%).

README setup frame 1: 115.671 ms.

## Phase-by-phase readout

Ember Lattice runs throughout the slot under camera wander, sine lightness and displacement-driven hue. Counter trees give inclusive averages per frame. Whole in-show readouts provide the window statistics.

### Startup (frames 1–32)

```
frame                         60.523 ms   36.314 Mcyc 100.0% x1.00
  pov_preserve_half            0.127 ms    0.076 Mcyc   0.2% x0.97
  fx_shader_draw              59.899 ms   35.939 Mcyc  99.0% x1.00
    fx_prepare_column          1.481 ms    0.888 Mcyc   2.4% x144.00
  fx_palette_base_bake         0.060 ms    0.036 Mcyc   0.1% x1.00
  fx_prepare_frame             0.022 ms    0.013 Mcyc   0.0% x1.00
  fx_advance                   0.131 ms    0.078 Mcyc   0.2% x1.00
  fx_timeline_step             0.010 ms    0.006 Mcyc   0.0% x1.00
  canvas_clear                 0.091 ms    0.055 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.003 ms    0.002 Mcyc   0.0% x1.00
```

Wall min/avg/max = 57.999/60.522/115.674 ms. Mean render was 60.520 ms and the shader draw averaged 59.899 ms.

### Peak live render (frames 513–544)

```
frame                         58.786 ms   35.272 Mcyc 100.0% x1.00
  pov_preserve_half            0.133 ms    0.080 Mcyc   0.2% x1.00
  fx_shader_draw              58.230 ms   34.938 Mcyc  99.1% x1.00
    fx_prepare_column          1.492 ms    0.895 Mcyc   2.5% x144.00
  fx_palette_base_bake         0.058 ms    0.035 Mcyc   0.1% x1.00
  fx_prepare_frame             0.021 ms    0.012 Mcyc   0.0% x1.00
  fx_advance                   0.070 ms    0.042 Mcyc   0.1% x1.00
  fx_timeline_step             0.010 ms    0.006 Mcyc   0.0% x1.00
  canvas_clear                 0.087 ms    0.052 Mcyc   0.1% x1.00
  canvas_buffer_wait           0.003 ms    0.002 Mcyc   0.0% x1.00
```

Wall min/avg/max = 57.983/58.786/60.932 ms. Mean render was 58.784 ms and the shader draw averaged 58.230 ms.

### Show end (frames 545–576)

```
frame                         58.898 ms   35.339 Mcyc 100.0% x1.00
  pov_preserve_half            0.132 ms    0.079 Mcyc   0.2% x1.00
  fx_shader_draw              58.335 ms   35.001 Mcyc  99.0% x1.00
    fx_prepare_column          1.447 ms    0.868 Mcyc   2.5% x144.00
  fx_palette_base_bake         0.059 ms    0.035 Mcyc   0.1% x1.00
  fx_prepare_frame             0.020 ms    0.012 Mcyc   0.0% x1.00
  fx_advance                   0.070 ms    0.042 Mcyc   0.1% x1.00
  fx_timeline_step             0.009 ms    0.005 Mcyc   0.0% x1.00
  canvas_clear                 0.087 ms    0.052 Mcyc   0.1% x1.00
  canvas_buffer_wait           0.003 ms    0.002 Mcyc   0.0% x1.00
```

Wall min/avg/max = 58.038/58.898/59.431 ms. Mean render was 58.895 ms and the shader draw averaged 58.335 ms.

### Per-pixel figures

The clipped quadrant contains 10,368 output pixels. This capture has no blended-pixel counter.

## Column-ISR / DMA marshaling cost

Peak window frames 513–544. CPU percentages use the separate 1,999,749 μs ISR snapshot interval. Times are min/avg/max.

```
isr_wake          1151.9/f 0.500/1.620/23.071 us 2.98%
isr_pack           144.0/f 5.033/6.171/8.918 us 1.42%
isr_dma_submit     144.0/f 0.620/0.920/9.196 us 0.21%
```

CPU packing averages 6.171 μs per call and CPU DMA submission averages 0.920 μs per call. At 24 MHz SPI the 600-byte composite has an asynchronous transfer budget of 230 μs including framing gaps. Pack and submit run inside the inclusive wake scope, which uses 2.98% of the interval.

## Summary ranking

1. `fx_shader_draw` — 99.1% of the foreground frame, 58.230 ms/f.
2. `fx_prepare_column` — 2.5% of the foreground frame, 1.492 ms/f.
3. `pov_preserve_half` — 0.2% of the foreground frame, 0.133 ms/f.

## Caveats

- Foreground cycle counts include interrupt time. Render telemetry subtracts the `canvas_buffer_wait` scope.
- The cutoff uses the first readout’s elapsed time, subsequent ISR snapshot intervals, and individual draw durations in the boundary readout. Per-frame render measurements come directly from the log. The full-playlist handoff was outside this measurement.
- Source and build provenance identify the clean measured checkout. The shipping slot uses only the Ember Lattice preset. No timing overrides or deep instrumentation were used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=EmberTracery`, `HS_PROFILE_WINDOW=32`. Acquire the board through `tools/profile_one.sh`.
