# JadeFiligree on-device profile — Teensy 4.0, segmented mode (2026-11-19, **selective -O3**)

Raw capture: `C:/work/temp/holosphere-20261119-tiling-profiles/column-time/jadefiligree_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` scan, spherical source, and pullback color paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Jade Filigree, 288×144, source `a0ae093c47a3852317f802d4e59366c9bc4885f1` |
| Method | 108 s capture, 16 frames per readout, default one-hour profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh JadeFiligree profile 108 16` |

Profile image size: `FLASH: code:93344, data:345688, headers:8456; free for files:1584128 / RAM1: variables:314144, code:27720, padding:5048; free for local variables:177376 / RAM2: variables:520064; free for malloc/new:4224`.

Shipping Phantasm image size: `FLASH: code:1119248, data:618036, headers:8628; free for files:285704 / RAM1: variables:313760, code:194744, padding:1864; free for local variables:13920 / RAM2: variables:520064; free for malloc/new:4224`.

Captured 2026-11-19 19:58, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `710097746c4fb3b8234f70ea8ad76cae3c0f867f002251b23a647e8002c61e62`. Profile ELF SHA-256 `6f8bf8c1b21da216423bf547c1f22ce721acd5fe0c151fe9054e71769b0bacb6`. Phantasm ELF SHA-256 `30aeea9589b8dda81022339b91911802f3f8bd3fe8e2a79c056cab545d3b655e`. The capture retains attested ELF, map and toolchain configuration paths.

`parse_profile.py ... validate` passed. In frames 1073–1088, root cycles divided by 600 MHz match the foreground wall sum to 0.41 ppm.

## Frame cadence

The Phantasm playlist gives Jade Filigree 105 s. The measurement includes 1,671 frames confirmed complete within the slot, including setup. The last included frame completed by 104.988 s. Analysis excludes 33 frame records beyond the cutoff, including up to 4 frames with an uncertain completion time at the boundary.

Setup frame 1 took **46.193 ms**. Frames 2–1671 averaged **40.164 ms** render and peaked at **52.181 ms** on frame 1080 (Jade Weave). **0/1670 (0.00%)** exceeded the 62.5 ms display budget.

README cells: peak 🟢 52.181 (2), spilled 🟢 0/1670 (0.00%).

README setup frame 1: 46.193 ms.

## Phase-by-phase readout

The authored show holds each preset for 600 frames and uses 480-frame sinusoidally eased transitions. Preset markers attribute a transition to the incoming preset. Counter trees give inclusive averages per frame. Whole in-show readouts provide the window statistics.

### Startup (frames 1–16)

```
frame                         25.423 ms   15.254 Mcyc 100.0% x1.00
  pov_preserve_half            0.127 ms    0.076 Mcyc   0.5% x0.94
  fx_shader_draw              24.417 ms   14.650 Mcyc  96.0% x1.00
    fx_prepare_column          0.405 ms    0.243 Mcyc   1.6% x144.00
  fx_palette_base_bake         0.666 ms    0.399 Mcyc   2.6% x1.00
    fx_palette_rgb_bake        0.591 ms    0.355 Mcyc   2.3% x1.00
  fx_prepare_frame             0.008 ms    0.005 Mcyc   0.0% x1.00
  fx_advance                   0.008 ms    0.005 Mcyc   0.0% x1.00
  fx_timeline_step             0.010 ms    0.006 Mcyc   0.0% x1.00
  canvas_clear                 0.096 ms    0.058 Mcyc   0.4% x1.00
  canvas_buffer_wait           0.000 ms    0.000 Mcyc   0.0% x1.00
```

Wall min/avg/max = 23.990/25.422/46.196 ms. Mean render was 25.422 ms and the shader draw averaged 24.417 ms.

### Transition into Jade Weave (frames 929–944)

```
frame                         51.019 ms   30.612 Mcyc 100.0% x1.00
  pov_preserve_half            0.132 ms    0.079 Mcyc   0.3% x1.00
  fx_shader_draw              49.470 ms   29.682 Mcyc  97.0% x1.00
    fx_prepare_column          2.353 ms    1.412 Mcyc   4.6% x144.00
  fx_prepare_frame             1.157 ms    0.694 Mcyc   2.3% x1.00
  fx_advance                   0.009 ms    0.005 Mcyc   0.0% x1.00
  fx_timeline_step             0.091 ms    0.054 Mcyc   0.2% x1.00
  canvas_clear                 0.087 ms    0.052 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.000 ms    0.000 Mcyc   0.0% x1.00
```

Wall min/avg/max = 50.835/51.019/51.179 ms. Mean render was 51.019 ms and the shader draw averaged 49.470 ms.

### Peak live render (frames 1073–1088)

```
frame                         48.948 ms   29.369 Mcyc 100.0% x1.00
  pov_preserve_half            0.134 ms    0.080 Mcyc   0.3% x1.00
  fx_shader_draw              47.511 ms   28.507 Mcyc  97.1% x1.00
    fx_prepare_column          1.468 ms    0.881 Mcyc   3.0% x144.00
  fx_palette_base_bake         0.366 ms    0.219 Mcyc   0.7% x0.50
    fx_palette_rgb_bake        0.331 ms    0.198 Mcyc   0.7% x0.50
  fx_prepare_frame             0.590 ms    0.354 Mcyc   1.2% x1.00
  fx_advance                   0.124 ms    0.074 Mcyc   0.3% x1.00
  fx_timeline_step             0.051 ms    0.031 Mcyc   0.1% x1.00
  canvas_clear                 0.088 ms    0.053 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.000 ms    0.000 Mcyc   0.0% x1.00
```

Wall min/avg/max = 40.659/48.947/52.181 ms. Mean render was 48.947 ms and the shader draw averaged 47.511 ms.

### Show end (frames 1649–1664)

```
frame                         48.015 ms   28.809 Mcyc 100.0% x1.00
  pov_preserve_half            0.134 ms    0.080 Mcyc   0.3% x1.00
  fx_shader_draw              47.027 ms   28.216 Mcyc  97.9% x1.00
    fx_prepare_column          0.587 ms    0.352 Mcyc   1.2% x144.00
  fx_palette_base_bake         0.652 ms    0.391 Mcyc   1.4% x1.00
    fx_palette_rgb_bake        0.589 ms    0.353 Mcyc   1.2% x1.00
  fx_prepare_frame             0.009 ms    0.005 Mcyc   0.0% x1.00
  fx_advance                   0.008 ms    0.005 Mcyc   0.0% x1.00
  fx_timeline_step             0.012 ms    0.007 Mcyc   0.0% x1.00
  canvas_clear                 0.087 ms    0.052 Mcyc   0.2% x1.00
  canvas_buffer_wait           0.000 ms    0.000 Mcyc   0.0% x1.00
```

Wall min/avg/max = 47.686/48.015/48.413 ms. Mean render was 48.015 ms and the shader draw averaged 47.027 ms.

### Per-preset table

Hold figures use the costliest whole readout by shader draw within the authored dwell. Live peaks and spill counts include that preset’s incoming transition.

| Preset | Hold windows | Representative frames | Draw ms/f | Render ms/f | Live peak ms | Live spilled |
|---|--:|---|--:|--:|--:|--:|
| 1 · Golden Mesh | 36 | 417–432 | 23.099 | 24.049 | 24.150 | 0/598 (0.00%) |
| 2 · Jade Weave | 36 | 1441–1456 | 47.195 | 48.190 | 52.181 | 0/1072 (0.00%) |

### Per-pixel figures

The clipped quadrant contains 10,368 output pixels. This capture has no blended-pixel counter.

## Column-ISR / DMA marshaling cost

Peak window frames 1073–1088. CPU percentages use the separate 996,936 μs ISR snapshot interval. Times are min/avg/max.

```
isr_wake          1148.5/f 0.483/1.556/19.280 us 2.86%
isr_pack           143.6/f 4.708/5.836/9.803 us 1.34%
isr_dma_submit     143.6/f 0.663/0.915/1.780 us 0.21%
```

CPU packing averages 5.836 μs per call and CPU DMA submission averages 0.915 μs per call. At 24 MHz SPI the 600-byte composite has an asynchronous transfer budget of 230 μs including framing gaps. Pack and submit run inside the inclusive wake scope, which uses 2.86% of the interval.

## Summary ranking

1. `fx_shader_draw` — 97.1% of the foreground frame, 47.511 ms/f.
2. `fx_prepare_column` — 3.0% of the foreground frame, 1.468 ms/f.
3. `fx_prepare_frame` — 1.2% of the foreground frame, 0.590 ms/f.

## Caveats

- Foreground cycle counts include interrupt time. Render telemetry subtracts the `canvas_buffer_wait` scope.
- The cutoff uses the first readout’s elapsed time, subsequent ISR snapshot intervals, and individual draw durations in the boundary readout. Per-frame render measurements come directly from the log. The full-playlist handoff was outside this measurement.
- Source and build provenance identify the clean measured checkout. The shipping preset order, dwell and easing were active. No timing overrides or deep instrumentation were used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=JadeFiligree`, `HS_PROFILE_WINDOW=16`. Acquire the board through `tools/profile_one.sh`.
