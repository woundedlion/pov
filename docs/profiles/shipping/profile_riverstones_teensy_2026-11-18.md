# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\shared-mod-evidence\paired-bands.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `e2e5353c810d11b962dbaa945932da6a9af06fce` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1–434 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:111904, data:156720, headers:8880; free:1754112 / RAM1: variables:314080, code:26264, padding:6504; free:177440 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1172064, data:430244, headers:8436; free:420872 / RAM1: variables:313760, code:175528, padding:21080; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 2307773.757 us; measured wall sum 2307776 us, delta 0.97 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 11:18, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `52425b480c160eff9d986cacce3426b20fa0d67f279a91a066988860b86efa7f`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Revised artwork and shared renderer

The desktop revision doubles the burst rate from 1/128 to 1/64 and adds a translating affine warp between the projection and the hexagonal cell sampler. Both columns below use this artwork on Teensy 4.0 serial 6362390 (COM4), with the shipping 600 MHz segmented driver and the first 30 seconds of the opening slot. The handwritten build is `cb0e7e45fab4` and the shared build is `e2e5353c810d`.

| Measurement | Handwritten renderer | Shared renderer |
|---|---:|---:|
| Peak live render (ms) | 79.542 | 75.015 |
| Live frames over 62.5 ms | 367/406 | 381/433 |
| Setup (ms) | 94.713 | 91.667 |
| Full show flash code (bytes) | 1,174,376 | 1,172,064 |
| Full show flash data (bytes) | 430,136 | 430,244 |
| Marginal RiverStones code (bytes) | 17,976 | 16,664 |
| Marginal RiverStones data (bytes) | 372 | 404 |
| RAM1 code (bytes) | 175,512 | 175,528 |
| RAM1 variables (bytes) | 313,760 | 313,760 |

Marginal flash compares each full shipping image with the same source built with only RiverStones excluded. The previous artwork took 65.540 ms at peak with 113/471 live frames over budget on this board. The 480-frame revised artwork matched the handwritten renderer across all 4,976,640 output pixels (RGB16 SHA-256 `4435faa3cdf266170785fb49847511d65f0dbbfc21b8af9da8c2fe22b00679ad`). Separate tests compare straight alpha, edited controls and temporal endpoints with the independent stage evaluator.

During development, an ARM program-counter capture exposed duplicate hue-noise sampling. A serial diagnostic found four differing padding bytes in otherwise equivalent generators. The shared renderer compares the sampling inputs and keeps authored layer order. Both qualification captures used PJRC Teensy Loader CLI 2.3.

The outer camera retains zero spin and 0.1 wander. The published shader has no RiverStones-specific color loop. Its three authored noise layers use the shared column renderer. The affine coordinate transform and animated ripples remain active. ARM code generation removed the closest-lattice search and distance calculation for this fixed-spatial color stack. A displacement-driven layer still receives the measured distance.

During the 30 s slot, 52 of 433 live frames met the 62.5 ms budget, compared with 39 of 406 in the handwritten build. The spill rates were 88.0% and 90.4% respectively. Iridescence, which uses column-phase modulation, went from 48.106 to 47.242 ms peak on COM3, with zero spills in both 479-frame comparisons. Its setup fell from 91.760 to 88.278 ms.

The shipping code total fell by 2,312 bytes. Flash data grew by 108 bytes. RAM1 code grew by 16 bytes, with unchanged 13,920-byte stack headroom and unchanged RAM2 use. The six-effect native comparison covered 1,216 quadrant frames.

### Remaining cost

In a separate COM4 program-counter run on this revision, grouped ripples accounted for 22.8% of 253,229 render samples, the three-layer color kernel for 15.9%, color conversion for 9.8%, and noise-cube baking for 8.6%. The sampler reported zero dropped samples. The live frame times above come from the ordinary cycle-counter run. In a deterministic CPU workload covering 3,932,160 input directions, the new per-center bounds rejected 2,936,702 additional groups before the angle lookup.

In the full-versus-excluded shipping ELFs the largest additional code symbols were column-time frame interpolation (1,940 bytes), the hue/lift/opacity kernel (1,546), effect initialization (1,400), and column drawing (1,344). This marginal build includes the renderer and scheduling code that becomes reachable when RiverStones is in the playlist. The shipping HEX for the profiled commit and the final source commit `59f892bf2` matched byte for byte (SHA-256 `06889ce5068426beeba7a61fdc363cf5c4b1383f876d042da572ab0e0424069c`). The RiverStones profiling HEX also matched.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **91.667 ms**, spilled **382/434**. Setup frame 1 is **91.667 ms** and remains separate. Live show peak is **75.015 ms**, spilled **381/433**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–434. The raw log retains 477 frame records; the 43 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🔴 75.015, spilled 🔴 381/433 (87.99%).

README setup frame 1: 91.667 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    63.085 ms  37.851 Mcyc 100.0%
  pov_preserve_half      131.49 us   78.89 kcyc   0.2% x0.96875 136 us/call
  fx_shader_draw         48.611 ms  29.166 Mcyc  77.1%
    fx_prepare_column     2.028 ms   1.217 Mcyc   3.2% x144 14 us/call
  fx_prepare_frame        6.069 ms   3.641 Mcyc   9.6% x1 6069 us/call
  fx_advance              1.531 ms  918.81 kcyc   2.4%
    palette_cycle_step    1.475 ms  884.86 kcyc   2.3% x1 1475 us/call
  fx_timeline_step        92.64 us   55.59 kcyc   0.1% x1 93 us/call
  canvas_clear           100.37 us   60.22 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait      6.358 ms   3.815 Mcyc  10.1% x1 6358 us/call
```

Wall min/avg/max = 51.218/63.084/91.668 ms. `fx_shader_draw` averages 48.611 ms/f; total render averages 56.726 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    72.118 ms  43.271 Mcyc 100.0%
  pov_preserve_half      137.33 us   82.40 kcyc   0.2% x1 137 us/call
  fx_shader_draw         64.050 ms  38.430 Mcyc  88.8%
    fx_prepare_column     2.516 ms   1.510 Mcyc   3.5% x144 17 us/call
  fx_prepare_frame        5.986 ms   3.591 Mcyc   8.3% x1 5986 us/call
  fx_advance              1.453 ms  871.65 kcyc   2.0%
    palette_cycle_step    1.402 ms  841.29 kcyc   1.9% x1 1402 us/call
  fx_timeline_step        96.26 us   57.76 kcyc   0.1% x1 96 us/call
  canvas_clear           106.80 us   64.08 kcyc   0.1% x1 107 us/call
  canvas_buffer_wait      96.85 us   58.11 kcyc   0.1% x1 97 us/call
```

Wall min/avg/max = 69.274/72.118/75.320 ms. `fx_shader_draw` averages 64.050 ms/f; total render averages 72.021 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 385–416)

```
frame                    72.286 ms  43.371 Mcyc 100.0%
  pov_preserve_half      137.13 us   82.28 kcyc   0.2% x1 137 us/call
  fx_shader_draw         64.245 ms  38.547 Mcyc  88.9%
    fx_prepare_column     2.554 ms   1.532 Mcyc   3.5% x144 18 us/call
  fx_prepare_frame        5.983 ms   3.590 Mcyc   8.3% x1 5983 us/call
  fx_advance              1.445 ms  866.91 kcyc   2.0%
    palette_cycle_step    1.395 ms  836.75 kcyc   1.9% x1 1395 us/call
  fx_timeline_step        99.81 us   59.89 kcyc   0.1% x1 100 us/call
  canvas_clear           109.27 us   65.56 kcyc   0.2% x1 109 us/call
  canvas_buffer_wait      73.12 us   43.87 kcyc   0.1% x1 73 us/call
```

Wall min/avg/max = 70.284/72.285/73.723 ms. `fx_shader_draw` averages 64.245 ms/f; total render averages 72.213 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 2312371 us.

```
isr_wake       1332.0/f 0.491/1.766/23.885 us 3.25%
  isr_pack       166.5/f 6.516/7.398/11.228 us 1.70%
  isr_dma_submit 166.5/f 0.611/0.908/5.660 us 0.20%
```

- Mean pack 7.398 us/call versus submit 0.908 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.25%, corresponding to 60.469 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.908×; Peak live show render: 1.152×; Show end: 1.155×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 88.8% of the frame, 64.050 ms/f.
2. `fx_prepare_frame` — 8.3% of the frame, 5.986 ms/f.
3. `fx_prepare_column` — 3.5% of the frame, 2.516 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 434 completed by 29.968216 s from profile instance construction. A separate full-playlist telemetry run on COM4 confirmed the transition from RiverStones (index 0) to LiquidWindows (index 1). Initial construction and priming took 101.426 ms. The next effect took 73.862 ms to build, using 150.429 ms of the 249.999 ms handoff window and leaving 99.570 ms of margin. See the adjacent raw capture, bands-handoff.log.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
