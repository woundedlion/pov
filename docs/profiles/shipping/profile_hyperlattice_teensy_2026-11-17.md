# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\hyperlattice_ship.log`. Replaces `profile_hyperlattice_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 300 s / 4800 nominal frames at 16 fps; window 16. Raw capture 330 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 303 16`; analyze only the 300 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:136640, data:148760, headers:8488; free:1737728 / RAM1: variables:313024, code:13784, padding:18984; free:178496 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 3201–3216 root cycles / 600 MHz = 993082.925 us; measured wall sum 993081 us, delta 1.94 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:08, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `d4b604567e2128f4c9f9dfa4932b3d36dc60aec2cf9adf4b5fe6f2a70e66942e`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 300 s (4800 nominal frames at 16 fps). Show-frame peak render is **53.311 ms**, spilled **0/4800**. Setup frame 1 is **45.039 ms** and remains separate. Live show peak is **53.311 ms**, spilled **0/4799**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–4800. The raw log retains 5257 frame records; the 457 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 53.311 (9), spilled 🟢 0/4799 (0.00%).

README setup frame 1: 45.039 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    59.066 ms  35.440 Mcyc 100.0%
  pov_preserve_half      133.81 us   80.29 kcyc   0.2% x0.9375 143 us/call
  hl_shader_draw         22.522 ms  13.513 Mcyc  38.1% x0.9375 24023 us/call
  hl_shader_draw          2.777 ms   1.666 Mcyc   4.7% x0.0625 44424 us/call
  palette_cycle_step      1.562 ms  937.11 kcyc   2.6% x1 1562 us/call
  hl_timeline_step        16.91 us   10.14 kcyc   0.0% x1 17 us/call
  canvas_clear            94.12 us   56.47 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     31.652 ms  18.991 Mcyc  53.6% x1 31652 us/call
```

Wall min/avg/max = 25.615/59.066/64.327 ms. `hl_shader_draw` averages 2.777 ms/f; total render averages 27.414 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 3201–3216)

```
frame                    62.068 ms  37.241 Mcyc 100.0%
  pov_preserve_half      137.45 us   82.47 kcyc   0.2% x1 137 us/call
  hl_shader_draw         42.400 ms  25.440 Mcyc  68.3% x1 42400 us/call
  palette_cycle_step      1.696 ms   1.017 Mcyc   2.7% x1 1696 us/call
  hl_timeline_step        23.65 us   14.19 kcyc   0.0% x1 24 us/call
  canvas_clear            96.31 us   57.78 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     17.366 ms  10.420 Mcyc  28.0% x1 17366 us/call
```

Wall min/avg/max = 51.899/62.067/73.785 ms. `hl_shader_draw` averages 42.400 ms/f; total render averages 44.701 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 3089–3104)

```
frame                    62.477 ms  37.486 Mcyc 100.0%
  pov_preserve_half      136.87 us   82.12 kcyc   0.2% x1 137 us/call
  hl_shader_draw         43.149 ms  25.889 Mcyc  69.1% x1 43149 us/call
  palette_cycle_step      1.700 ms   1.020 Mcyc   2.7% x1 1700 us/call
  hl_timeline_step        14.91 us    8.95 kcyc   0.0% x1 15 us/call
  canvas_clear            93.87 us   56.32 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     17.049 ms  10.229 Mcyc  27.3% x1 17049 us/call
```

Wall min/avg/max = 59.884/62.477/64.985 ms. `hl_shader_draw` averages 43.149 ms/f; total render averages 45.429 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 3105–3120)

```
frame                    62.423 ms  37.454 Mcyc 100.0%
  pov_preserve_half      138.25 us   82.95 kcyc   0.2% x1 138 us/call
  hl_shader_draw         42.618 ms  25.571 Mcyc  68.3% x1 42618 us/call
  palette_cycle_step      1.709 ms   1.026 Mcyc   2.7% x1 1709 us/call
  hl_timeline_step        18.94 us   11.36 kcyc   0.0% x1 19 us/call
  canvas_clear            95.33 us   57.20 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     17.503 ms  10.502 Mcyc  28.0% x1 17503 us/call
```

Wall min/avg/max = 59.038/62.423/65.989 ms. `hl_shader_draw` averages 42.618 ms/f; total render averages 44.920 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 4785–4800)

```
frame                    62.474 ms  37.484 Mcyc 100.0%
  pov_preserve_half      141.68 us   85.01 kcyc   0.2% x1 142 us/call
  hl_shader_draw         31.439 ms  18.863 Mcyc  50.3% x1 31439 us/call
  palette_cycle_step      1.654 ms  992.23 kcyc   2.6% x1 1654 us/call
  hl_timeline_step        18.85 us   11.31 kcyc   0.0% x1 19 us/call
  canvas_clear            84.69 us   50.81 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     28.784 ms  17.270 Mcyc  46.1% x1 28784 us/call
```

Wall min/avg/max = 61.464/62.473/63.486 ms. `hl_shader_draw` averages 31.439 ms/f; total render averages 33.690 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 9 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 6 · octet-4d-flight | 19 | 3089–3104 | held | — | — | 43.149 | 45.429 | 16 | 53.311 | 0/559 |
| 9 · shell-4d-flight | 20 | 4593–4608 | held | — | — | 33.286 | 35.525 | 16 | 36.882 | 0/448 |
| 3 · hypercube-flight | 19 | 1185–1200 | held | — | — | 29.166 | 31.407 | 16 | 39.169 | 0/559 |
| 5 · octet-wide-flight | 19 | 2273–2288 | held | — | — | 29.138 | 31.348 | 16 | 31.945 | 0/679 |
| 8 · shell-close-flight | 19 | 4145–4160 | held | — | — | 29.086 | 31.340 | 16 | 34.329 | 0/679 |
| 4 · octet-flight | 19 | 1809–1824 | held | — | — | 27.778 | 29.985 | 16 | 32.844 | 0/439 |
| 2 · cubic-wide-flight | 19 | 769–784 | held | — | — | 25.769 | 27.964 | 16 | 30.181 | 0/679 |
| 7 · shell-flight | 19 | 3553–3568 | held | — | — | 25.069 | 27.321 | 16 | 30.192 | 0/439 |
| 1 · cubic-flight | 18 | 97–112 | held | — | — | 24.803 | 27.029 | 16 | 27.915 | 0/318 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 3201–3216; CPU shares use its separately snapshotted ISR interval of 994451 us.

```
isr_wake       1145.7/f 0.361/1.560/17.806 us 2.87%
  isr_pack       143.2/f 5.511/7.118/10.053 us 1.64%
  isr_dma_submit 143.2/f 0.756/0.915/0.980 us 0.21%
```

- Mean pack 7.118 us/call versus submit 0.915 us/call (7.8×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.87%, corresponding to 60.706 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.439×; Peak live show render: 0.715×; Costliest qualified preset: 0.727×; Transition / non-hold regime: 0.719×; Show end: 0.539×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `hl_shader_draw` — 68.3% of the frame, 42.400 ms/f.
2. `palette_cycle_step` — 2.7% of the frame, 1.696 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.137 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
