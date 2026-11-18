# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river-opt-evidence\qualified.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `2ab3e4ef51a43caa2a9c0194cbcec237a88189bc` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1–472 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:111336, data:156200, headers:8944; free:1755136 / RAM1: variables:314080, code:28552, padding:4216; free:177440 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1173048, data:430136, headers:8592; free:419840 / RAM1: variables:313760, code:175512, padding:21096; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 289–320 root cycles / 600 MHz = 2027231.765 us; measured wall sum 2027235 us, delta 1.60 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 09:41, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `f476a6dae8f23ee71b32232f956409d697ff958a2d2e178095635e94227b6741`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Optimization measurements

The same desktop Teensy 4.0 (COM4, serial 6362390) ran the original promotion `3af89d143576` and the optimized source `2ab3e4ef51a4`. Both captures used the synchronized segment 0 driver at 600 MHz and the first 30 seconds of the authored slot.

| Measurement | Original | Optimized |
|---|---:|---:|
| Setup frame (ms) | 120.413 | 89.669 |
| Peak live frame (ms) | 102.660 | 65.479 |
| Live frames over 62.5 ms | 304/304 | 116/471 |
| Full show flash code (bytes) | 1,198,776 | 1,173,048 |
| Full show flash data (bytes) | 430,764 | 430,136 |
| Marginal RiverStones code (bytes) | 18,672 | 16,648 |
| Marginal RiverStones data (bytes) | 368 | 372 |
| Full show RAM1 code (bytes) | 175,960 | 175,512 |

Flash attribution compares each full shipping build with its corresponding build with only RiverStones excluded. The optimized excluded build contains 1,156,400 bytes of code and 429,764 bytes of data. The complete change saves 25,728 code bytes across the show. Within RiverStones, the shared preset initialization reduced its constructor from 1,330 to 112 bytes, and its preset-application entry shrank from 794 to 46 bytes. RAM1 variables use 313,760 bytes, with 13,920 bytes available for locals. RAM2 has 4,224 bytes free.

In the original 10 kHz program-counter sample, paired modulation evaluation accounted for 17.4% of render samples and OpenSimplex2 for 11.1%. The original 3,286-byte pixel shader had 33 call sites, including five calls to `memcpy`. RiverStones now shares the static hue-noise sample between temporal endpoints, computes one cube address for the ordered lift and opacity layers, and rejects pixels outside a same-center ripple group's complete support. The authored shader binds identity palette controls once per column and uses inline signed-byte cube interpolation. The candidate with these changes and the table-free ripple specialization reduced the measured peak from 70.379 to 65.479 ms. Edited modulation settings retain their fallback.

The final sampler recorded 243,788 render samples over 26.059 seconds, with zero dropped. The color tail accounted for 19.1%, grouped ripple geometry 15.5%, the outer geometry shader 10.9%, and OpenSimplex2 9.6%. The full-show ELF adds 2,112 bytes for the color tail, 1,912 for temporal frame interpolation, and 1,196 for the column scan. Raw samples and symbol attribution are in `verified_pc.log`, `verified_pc_report.txt`, and `qualified-symbols.txt` beside the timing log.

The original and optimized host renders produced the same RGB16 bytes for 480 deterministic quadrant frames (4,976,640 pixels, SHA-256 `a6f0e7e66db3c7ed760fac2da24fc97c4c6f12a93de781c4ab04cc408a0780a7`). Focused tests also compared straight alpha and edited palettes and mechanics against the pipeline evaluator.

The bench used PJRC Teensy Loader CLI under the per-board lock. The packaged 2.2 uploader failed on an early HID write. The subsequent full-show upload used PJRC Loader 2.3. Both timing captures reported one power-on reset.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **89.669 ms**, spilled **117/472**. Setup frame 1 is **89.669 ms** and remains separate. Live show peak is **65.479 ms**, spilled **116/471**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–472. The raw log retains 520 frame records; the 48 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🔴 65.479, spilled 🔴 116/471 (24.63%).

README setup frame 1: 89.669 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    63.005 ms  37.803 Mcyc 100.0%
  pov_preserve_half      131.75 us   79.05 kcyc   0.2% x0.96875 136 us/call
  fx_shader_draw         47.681 ms  28.609 Mcyc  75.7%
    fx_prepare_column     1.206 ms  723.31 kcyc   1.9% x144 8 us/call
  fx_prepare_frame        6.054 ms   3.632 Mcyc   9.6% x1 6054 us/call
  fx_advance              1.201 ms  720.63 kcyc   1.9%
    palette_cycle_step    1.140 ms  684.22 kcyc   1.8% x1 1140 us/call
  fx_timeline_step        64.39 us   38.63 kcyc   0.1% x1 64 us/call
  canvas_clear            92.79 us   55.67 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait      7.601 ms   4.561 Mcyc  12.1% x1 7601 us/call
```

Wall min/avg/max = 50.454/63.004/89.671 ms. `fx_shader_draw` averages 47.681 ms/f; total render averages 55.404 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 289–320)

```
frame                    63.351 ms  38.011 Mcyc 100.0%
  pov_preserve_half      138.40 us   83.04 kcyc   0.2% x1 138 us/call
  fx_shader_draw         54.958 ms  32.975 Mcyc  86.8%
    fx_prepare_column     1.372 ms  823.06 kcyc   2.2% x144 10 us/call
  fx_prepare_frame        5.963 ms   3.578 Mcyc   9.4% x1 5963 us/call
  fx_advance              1.183 ms  709.96 kcyc   1.9%
    palette_cycle_step    1.126 ms  675.48 kcyc   1.8% x1 1126 us/call
  fx_timeline_step        61.83 us   37.10 kcyc   0.1% x1 62 us/call
  canvas_clear            93.58 us   56.15 kcyc   0.1% x1 94 us/call
  canvas_buffer_wait     775.43 us  465.26 kcyc   1.2% x1 775 us/call
```

Wall min/avg/max = 59.568/63.351/65.781 ms. `fx_shader_draw` averages 54.958 ms/f; total render averages 62.576 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    63.306 ms  37.983 Mcyc 100.0%
  pov_preserve_half      138.67 us   83.20 kcyc   0.2% x1 139 us/call
  fx_shader_draw         55.361 ms  33.216 Mcyc  87.4%
    fx_prepare_column     1.368 ms  821.04 kcyc   2.2% x144 10 us/call
  fx_prepare_frame        5.966 ms   3.580 Mcyc   9.4% x1 5966 us/call
  fx_advance              1.218 ms  730.94 kcyc   1.9%
    palette_cycle_step    1.162 ms  697.37 kcyc   1.8% x1 1162 us/call
  fx_timeline_step        60.10 us   36.06 kcyc   0.1% x1 60 us/call
  canvas_clear           101.35 us   60.81 kcyc   0.2% x1 101 us/call
  canvas_buffer_wait     284.29 us  170.57 kcyc   0.4% x1 284 us/call
```

Wall min/avg/max = 61.165/63.305/65.388 ms. `fx_shader_draw` averages 55.361 ms/f; total render averages 63.022 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 289–320; CPU shares use its separately snapshotted ISR interval of 2031472 us.

```
isr_wake       1170.2/f 0.513/1.785/25.773 us 3.29%
  isr_pack       146.2/f 6.503/7.418/11.255 us 1.70%
  isr_dma_submit 146.2/f 0.565/0.918/6.765 us 0.21%
```

- Mean pack 7.418 us/call versus submit 0.918 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.29%, corresponding to 60.444 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.886×; Peak live show render: 1.001×; Show end: 1.008×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 86.8% of the frame, 54.958 ms/f.
2. `fx_prepare_frame` — 9.4% of the frame, 5.963 ms/f.
3. `fx_prepare_column` — 2.2% of the frame, 1.372 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 472 completed by 29.940544 s from profile instance construction. A separate full-show telemetry build recorded index 0 (RiverStones) followed by index 1 (LiquidWindows). The second build window used 138.864 ms of its 249.999 ms budget, leaving 111.135 ms. See `show-handoff.log` beside the timing log.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
