# JadeFiligree on-device profile — Teensy 4.0, segmented mode (2026-11-19, **selective -O3**)

Raw capture: `C:/work/temp/holosphere-20261119-jade-weave-profile/jadefiligree_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` scan, spherical source, and pullback color paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Jade Filigree, 288×144, source `af3a1eb54d0d027029d84cac0db62752c1a44b9e` |
| Method | 42 s capture, 32 frames per readout, default one-hour profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh JadeFiligree profile 42 32` |

Profile image size: `FLASH: code:91512, data:345324, headers:8604; free for files:1586176 / RAM1: variables:314144, code:27576, padding:5192; free for local variables:177376 / RAM2: variables:520064; free for malloc/new:4224`.

Shipping Phantasm image size: `FLASH: code:1116504, data:618012, headers:8332; free for files:288768 / RAM1: variables:313760, code:194744, padding:1864; free for local variables:13920 / RAM2: variables:520064; free for malloc/new:4224`.

Compared with the preceding Phantasm build, FLASH code changed by -1,064 B and FLASH data by -12 B. RAM1 code changed by +0 B and variables by +0 B.

Captured 2026-11-19 22:09, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `2451b1b61af0352281a13f16a4903f823a49a1f69574c31047aa5480873625b0`. Profile ELF SHA-256 `32fa1b690bd152522649477e2816b8c349544ec28c91a940d40b0b949c41e884`. Phantasm ELF SHA-256 `4b6202bd3814027339bf39c6c3105ca6995fcf22e322b608ad5288b1dd37809c`. The capture retains attested ELF, map and toolchain configuration paths.

`parse_profile.py ... validate` passed. In frames 417–448, root cycles divided by 600 MHz match the foreground wall sum to 3.26 ppm.

## Frame cadence

The Phantasm playlist gives Jade Filigree 38 s. The measurement includes 599 frames confirmed complete within the slot, including setup. The last included frame completed by 37.989 s. Analysis excludes 49 frame records beyond the cutoff, including up to 9 frames with an uncertain completion time at the boundary.

Setup frame 1 took **95.949 ms**. Frames 2–599 averaged **47.873 ms** render and peaked at **48.421 ms** on frame 433. **0/598 (0.00%)** exceeded the 62.5 ms display budget.

README cells: peak 🟢 48.421, spilled 🟢 0/598 (0.00%).

README setup frame 1: 95.949 ms.

## Phase-by-phase readout

Jade Weave runs throughout the slot with an animated Direct Lens and the Algae palette. Counter trees give inclusive averages per frame. Whole in-show readouts provide the window statistics.

### Startup (frames 1–32)

```
frame                         49.373 ms   29.624 Mcyc 100.0% x1.00
  pov_preserve_half            0.129 ms    0.077 Mcyc   0.3% x0.97
  fx_shader_draw              48.326 ms   28.995 Mcyc  97.9% x1.00
    fx_prepare_column          0.734 ms    0.441 Mcyc   1.5% x144.00
  fx_palette_base_bake         0.664 ms    0.399 Mcyc   1.3% x1.00
    fx_palette_rgb_bake        0.611 ms    0.367 Mcyc   1.2% x1.00
  fx_prepare_frame             0.004 ms    0.003 Mcyc   0.0% x1.00
  fx_advance                   0.071 ms    0.043 Mcyc   0.1% x1.00
  fx_timeline_step             0.004 ms    0.002 Mcyc   0.0% x1.00
  canvas_clear                 0.093 ms    0.056 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.001 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 47.444/49.373/95.953 ms. Mean render was 49.372 ms and the shader draw averaged 48.326 ms.

### Peak live render (frames 417–448)

```
frame                         47.848 ms   28.709 Mcyc 100.0% x1.00
  pov_preserve_half            0.135 ms    0.081 Mcyc   0.3% x1.00
  fx_shader_draw              46.895 ms   28.137 Mcyc  98.0% x1.00
    fx_prepare_column          0.724 ms    0.434 Mcyc   1.5% x144.00
  fx_palette_base_bake         0.643 ms    0.386 Mcyc   1.3% x1.00
    fx_palette_rgb_bake        0.592 ms    0.355 Mcyc   1.2% x1.00
  fx_prepare_frame             0.004 ms    0.002 Mcyc   0.0% x1.00
  fx_advance                   0.007 ms    0.004 Mcyc   0.0% x1.00
  fx_timeline_step             0.003 ms    0.002 Mcyc   0.0% x1.00
  canvas_clear                 0.089 ms    0.053 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.001 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 47.314/47.847/48.422 ms. Mean render was 47.847 ms and the shader draw averaged 46.895 ms.

### Show end (frames 545–576)

```
frame                         47.755 ms   28.653 Mcyc 100.0% x1.00
  pov_preserve_half            0.133 ms    0.080 Mcyc   0.3% x1.00
  fx_shader_draw              46.802 ms   28.081 Mcyc  98.0% x1.00
    fx_prepare_column          0.715 ms    0.429 Mcyc   1.5% x144.00
  fx_palette_base_bake         0.643 ms    0.386 Mcyc   1.3% x1.00
    fx_palette_rgb_bake        0.592 ms    0.355 Mcyc   1.2% x1.00
  fx_prepare_frame             0.004 ms    0.002 Mcyc   0.0% x1.00
  fx_advance                   0.006 ms    0.004 Mcyc   0.0% x1.00
  fx_timeline_step             0.003 ms    0.002 Mcyc   0.0% x1.00
  canvas_clear                 0.089 ms    0.053 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.001 ms    0.001 Mcyc   0.0% x1.00
```

Wall min/avg/max = 47.427/47.754/48.030 ms. Mean render was 47.753 ms and the shader draw averaged 46.802 ms.

### Per-pixel figures

The clipped quadrant contains 10,368 output pixels. This capture has no blended-pixel counter.

## Column-ISR / DMA marshaling cost

Peak window frames 417–448. CPU percentages use the separate 2,000,243 μs ISR snapshot interval. Times are min/avg/max.

```
isr_wake          1152.2/f 0.490/1.613/22.451 us 2.97%
isr_pack           144.0/f 5.363/6.336/8.988 us 1.46%
isr_dma_submit     144.0/f 0.635/0.936/9.978 us 0.21%
```

CPU packing averages 6.336 μs per call and CPU DMA submission averages 0.936 μs per call. At 24 MHz SPI the 600-byte composite has an asynchronous transfer budget of 230 μs including framing gaps. Pack and submit run inside the inclusive wake scope, which uses 2.97% of the interval.

## Summary ranking

1. `fx_shader_draw` — 98.0% of the foreground frame, 46.895 ms/f.
2. `fx_prepare_column` — 1.5% of the foreground frame, 0.724 ms/f.
3. `fx_palette_base_bake` — 1.3% of the foreground frame, 0.643 ms/f.

## Caveats

- Foreground cycle counts include interrupt time. Render telemetry subtracts the `canvas_buffer_wait` scope.
- The cutoff uses the first readout’s elapsed time, subsequent ISR snapshot intervals, and individual draw durations in the boundary readout. Per-frame render measurements come directly from the log. The full-playlist handoff was outside this measurement.
- Source and build provenance identify the clean measured checkout. The shipping slot uses only the Jade Weave preset. No timing overrides or deep instrumentation were used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=JadeFiligree`, `HS_PROFILE_WINDOW=32`. Acquire the board through `tools/profile_one.sh`.
