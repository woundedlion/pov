# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\jewelmelt_ship.log`. Replaces `profile_jewelmelt_teensy_2026-11-13.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | JewelMelt 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh JewelMelt profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:106360, data:152016, headers:8888; free:1764352 / RAM1: variables:314144, code:29704, padding:3064; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 289–320 root cycles / 600 MHz = 1997251.888 us; measured wall sum 1997252 us, delta 0.06 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:45, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `3363155fb37cb0cade9b791dbc77c5c955f9e8cce335ffae633c219d327c2171`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **79.906 ms**, spilled **1/960**. Setup frame 1 is **79.906 ms** and remains separate. Live show peak is **52.547 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 52.547, spilled 🟢 0/959 (0.00%).

README setup frame 1: 79.906 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.257 ms  37.354 Mcyc 100.0%
  pov_preserve_half      134.98 us   80.99 kcyc   0.2% x0.96875 139 us/call
  fx_shader_draw         38.218 ms  22.931 Mcyc  61.4%
    fx_prepare_column    741.96 us  445.17 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake   841.38 us  504.83 kcyc   1.4%
    fx_palette_rgb_bake  714.12 us  428.47 kcyc   1.1% x1 714 us/call
  fx_prepare_frame        3.810 ms   2.286 Mcyc   6.1% x1 3810 us/call
  fx_advance             959.76 us  575.86 kcyc   1.5%
    palette_cycle_step   919.41 us  551.64 kcyc   1.5% x1 919 us/call
  fx_timeline_step        91.86 us   55.12 kcyc   0.1% x1 92 us/call
  canvas_clear            91.91 us   55.15 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     17.981 ms  10.789 Mcyc  28.9% x1 17981 us/call
```

Wall min/avg/max = 43.064/62.257/79.907 ms. `fx_shader_draw` averages 38.218 ms/f; total render averages 44.276 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 289–320)

```
frame                    62.414 ms  37.448 Mcyc 100.0%
  pov_preserve_half      140.78 us   84.47 kcyc   0.2% x1 141 us/call
  fx_shader_draw         38.583 ms  23.150 Mcyc  61.8%
    fx_prepare_column    737.15 us  442.29 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    2.411 ms   1.447 Mcyc   3.9%
    fx_palette_rgb_bake   2.281 ms   1.369 Mcyc   3.7% x1 2281 us/call
  fx_prepare_frame        3.870 ms   2.322 Mcyc   6.2% x1 3870 us/call
  fx_advance             925.94 us  555.56 kcyc   1.5%
    palette_cycle_step   887.95 us  532.77 kcyc   1.4% x1 888 us/call
  fx_timeline_step        91.66 us   54.99 kcyc   0.1% x1 92 us/call
  canvas_clear            88.61 us   53.16 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     16.168 ms   9.701 Mcyc  25.9% x1 16168 us/call
```

Wall min/avg/max = 52.710/62.414/72.259 ms. `fx_shader_draw` averages 38.583 ms/f; total render averages 46.246 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.307 ms  37.384 Mcyc 100.0%
  pov_preserve_half      140.04 us   84.02 kcyc   0.2% x1 140 us/call
  fx_shader_draw         37.042 ms  22.225 Mcyc  59.5%
    fx_prepare_column    774.53 us  464.72 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    1.366 ms  819.42 kcyc   2.2%
    fx_palette_rgb_bake   1.235 ms  741.01 kcyc   2.0% x1 1235 us/call
  fx_prepare_frame        3.787 ms   2.272 Mcyc   6.1% x1 3787 us/call
  fx_advance             888.98 us  533.39 kcyc   1.4%
    palette_cycle_step   850.95 us  510.57 kcyc   1.4% x1 851 us/call
  fx_timeline_step        94.39 us   56.64 kcyc   0.2% x1 94 us/call
  canvas_clear            87.45 us   52.47 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     18.773 ms  11.264 Mcyc  30.1% x1 18773 us/call
```

Wall min/avg/max = 56.373/62.307/67.961 ms. `fx_shader_draw` averages 37.042 ms/f; total render averages 43.535 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 289–320; CPU shares use its separately snapshotted ISR interval of 2001326 us.

```
isr_wake       1152.8/f 0.355/1.745/26.333 us 3.21%
  isr_pack       144.1/f 7.830/8.531/11.033 us 1.96%
  isr_dma_submit 144.1/f 0.595/0.918/8.453 us 0.21%
```

- Mean pack 8.531 us/call versus submit 0.918 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.21%, corresponding to 60.494 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.708×; Peak live show render: 0.740×; Show end: 0.697×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 61.8% of the frame, 38.583 ms/f.
2. `fx_prepare_frame` — 6.2% of the frame, 3.870 ms/f.
3. `fx_palette_base_bake` — 3.9% of the frame, 2.411 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=JewelMelt`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
