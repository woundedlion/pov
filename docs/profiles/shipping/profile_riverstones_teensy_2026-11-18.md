# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
The raw capture is at `C:\work\temp\riverstones-profile-20261118\riverstones_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz with live flywheel and DMA ISRs |
| Image | `profile` with shipping `-Os`, newlib-nano, and selective `HS_O3` in `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `3af89d143576e1952d569f799d1d8f4a65e3fb42` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1–305. Counter windows cover 32 frames. The raw capture ran for 34 s in a single 3600 s profiling epoch. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 HS_UPLOAD_BACKEND=vendor bash tools/profile_one.sh RiverStones profile 34 32` from a clean checkout of the pinned source. Analyze the first 30 s. |

The profile image uses `FLASH: code:119328, data:156232, headers:9112; free:1746944 / RAM1: variables:314144, code:30712, padding:2056; free:177376 / RAM2: variables:520064; free:4224`.

The shipping Phantasm image uses `FLASH: code:1198776, data:430764, headers:8860; free:393216 / RAM1: variables:313760, code:175960, padding:20648; free:13920 / RAM2: variables:520064; free:4224`.

For frames 193–224, root cycles / 600 MHz gives 3,233,992.932 us. The measured wall sum was 3,233,993 us, a 0.02 ppm difference. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 08:46, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw capture SHA-256 `63dc805ba4e1f7fd218f2b69bee5c41ccd0d6397b553ca890c3b6dee71d6e01b`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

RiverStones opens the playlist for **30 s**. The setup draw took **120.350 ms**. Live render times ranged from **88.831 ms** to **102.670 ms**, with the peak on frame 211. All **304 live frames** spilled past the 62.5 ms target.

At 480 RPM each half revolution provides a 62.5 ms display window, or 16 windows per second. RiverStones renders one quadrant of 10,368 pixels. Render telemetry measures each draw with `canvas_buffer_wait` subtracted.

The ranked index covers frames 2–305. The raw capture contains 334 complete frame records. The counter trees use complete windows inside the slot.

README cells: peak 🔴 102.670, spilled 🔴 304/304 (100.00%).

README setup frame 1: 120.350 ms.

## Phase-by-phase readout

The single preset runs throughout the slot. The counter trees show per-frame averages for startup, the window containing the live peak, and the final complete window in scope.

### Startup (frames 1–32)

```
frame                    93.281 ms  55.968 Mcyc 100.0%
  pov_preserve_half      133.53 us   80.12 kcyc   0.1% x0.96875 138 us/call
  fx_shader_draw         85.410 ms  51.246 Mcyc  91.6%
    fx_prepare_column     3.421 ms   2.052 Mcyc   3.7% x144 24 us/call
  fx_palette_base_bake    75.75 us   45.45 kcyc   0.1% x1 76 us/call
  fx_prepare_frame        6.037 ms   3.622 Mcyc   6.5% x1 6037 us/call
  fx_advance              1.091 ms  654.73 kcyc   1.2%
    palette_cycle_step    1.029 ms  617.61 kcyc   1.1% x1 1029 us/call
  fx_timeline_step       112.13 us   67.28 kcyc   0.1% x1 112 us/call
  canvas_clear            98.89 us   59.33 kcyc   0.1% x1 99 us/call
  canvas_buffer_wait     127.02 us   76.21 kcyc   0.1% x1 127 us/call
```

Wall min/avg/max = 88.833/93.280/120.354 ms. `fx_shader_draw` averages 85.410 ms/f. Total render averages 93.154 ms/f.

### Peak live show render (frames 193–224)

```
frame                   101.062 ms  60.637 Mcyc 100.0%
  pov_preserve_half      135.66 us   81.40 kcyc   0.1% x1 136 us/call
  fx_shader_draw         93.285 ms  55.971 Mcyc  92.3%
    fx_prepare_column     3.674 ms   2.204 Mcyc   3.6% x144 26 us/call
  fx_palette_base_bake    74.10 us   44.46 kcyc   0.1% x1 74 us/call
  fx_prepare_frame        5.972 ms   3.583 Mcyc   5.9% x1 5972 us/call
  fx_advance              1.092 ms  655.39 kcyc   1.1%
    palette_cycle_step    1.030 ms  618.02 kcyc   1.0% x1 1030 us/call
  fx_timeline_step       115.04 us   69.02 kcyc   0.1% x1 115 us/call
  canvas_clear            95.35 us   57.21 kcyc   0.1% x1 95 us/call
  canvas_buffer_wait     100.27 us   60.16 kcyc   0.1% x1 100 us/call
```

Wall min/avg/max = 98.533/101.062/102.879 ms. `fx_shader_draw` averages 93.285 ms/f. Total render averages 100.962 ms/f.

### Last complete window (frames 257–288)

```
frame                    96.026 ms  57.616 Mcyc 100.0%
  pov_preserve_half      137.15 us   82.29 kcyc   0.1% x1 137 us/call
  fx_shader_draw         88.235 ms  52.941 Mcyc  91.9%
    fx_prepare_column     3.573 ms   2.144 Mcyc   3.7% x144 25 us/call
  fx_palette_base_bake    73.00 us   43.80 kcyc   0.1% x1 73 us/call
  fx_prepare_frame        5.975 ms   3.585 Mcyc   6.2% x1 5975 us/call
  fx_advance              1.087 ms  652.38 kcyc   1.1%
    palette_cycle_step    1.023 ms  614.04 kcyc   1.1% x1 1023 us/call
  fx_timeline_step       118.97 us   71.38 kcyc   0.1% x1 119 us/call
  canvas_clear            96.11 us   57.66 kcyc   0.1% x1 96 us/call
  canvas_buffer_wait     113.48 us   68.09 kcyc   0.1% x1 113 us/call
```

Wall min/avg/max = 93.452/96.026/98.209 ms. `fx_shader_draw` averages 88.235 ms/f. Total render averages 95.913 ms/f.

### Per-pixel figures

The clipped quadrant has 10,368 output pixels. A per-blend cost would require a blended-pixel counter for the direct render path.

## Column-ISR / DMA marshaling cost

The ISR snapshot for frames 193–224 spans 3,238,789 us. CPU shares use that interval.

```
isr_wake       1865.7/f 0.491/1.773/18.853 us 3.27%
  isr_pack       233.2/f 6.550/7.370/11.450 us 1.69%
  isr_dma_submit 233.2/f 0.615/0.916/6.665 us 0.21%
```

- Packing averaged 7.370 us per call and DMA submission averaged 0.916 us per call.
- Pack and DMA submit counters measure CPU time within the wake ISR.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.27%, corresponding to 60.456 ms foreground capacity per display window. Wake includes pack and submit. DMA completion has separate interrupt work outside this snapshot.

The 62.5 ms target would require a 1.615× reduction in the peak window's mean render cost. The individual peak is 1.643× the target. These inclusive costs already account for interrupts.

## Summary ranking

1. `fx_shader_draw` — 92.3% of the frame, 93.285 ms/f.
2. `fx_prepare_frame` — 5.9% of the frame, 5.972 ms/f.
3. `fx_prepare_column` — 3.6% of the frame, 3.674 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` takes its parent from the first caller. Counters with MIXED-PARENT or DUPLICATE-NAME tags may share time with siblings.
- The capture uses the standard profiling counters, and render timings include their overhead.
- The shipping selective-O3 regions along this path are `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. The frame completion bounds put frame 305 at or before 29.955365 s and frame 306 at or after 30.050268 s from profile instance construction. Handoff timing requires a full-playlist capture.
- The capture used a clean source tree. The build attestations identify the exact profile and shipping images.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
