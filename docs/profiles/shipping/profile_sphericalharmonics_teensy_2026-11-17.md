# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\sphericalharmonics_ship.log`. Replaces `profile_sphericalharmonics_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 96 s / 1536 nominal frames at 16 fps; window 16. Raw capture 111 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh SphericalHarmonics profile 99 16`; analyze only the 96 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:43416, data:145172, headers:9044; free:1833984 / RAM1: variables:312896, code:12184, padding:20584; free:178624 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 1105–1120 root cycles / 600 MHz = 999221.250 us; measured wall sum 999222 us, delta 0.75 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:31, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `9a13dc4de5c2d4b80cf66f77eb441cb4b449bb215995278ae96c930b213b60dd`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 96 s (1536 nominal frames at 16 fps). Show-frame peak render is **21.305 ms**, spilled **0/1536**. Setup frame 1 is **21.305 ms** and remains separate. Live show peak is **12.841 ms**, spilled **0/1535**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1536. The raw log retains 1753 frame records; the 217 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 12.841 (24), spilled 🟢 0/1535 (0.00%).

README setup frame 1: 21.305 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    56.708 ms  34.025 Mcyc 100.0%
  pov_preserve_half      137.79 us   82.68 kcyc   0.2% x0.9375 147 us/call
  sh_rasterize           11.691 ms   7.015 Mcyc  20.6% x1 11691 us/call
  sh_timeline_step        25.56 us   15.34 kcyc   0.0% x1 26 us/call
  canvas_clear            92.34 us   55.40 kcyc   0.2% x1 92 us/call
  canvas_buffer_wait     44.742 ms  26.845 Mcyc  78.9% x1 44742 us/call
```

Wall min/avg/max = 11.348/56.708/62.895 ms. `sh_rasterize` averages 11.691 ms/f; total render averages 11.966 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 1105–1120)

```
frame                    62.451 ms  37.471 Mcyc 100.0%
  pov_preserve_half      146.78 us   88.07 kcyc   0.2% x1 147 us/call
  sh_rasterize           12.507 ms   7.504 Mcyc  20.0% x1 12507 us/call
  sh_timeline_step        20.23 us   12.14 kcyc   0.0% x1 20 us/call
  canvas_clear            82.34 us   49.40 kcyc   0.1% x1 82 us/call
  canvas_buffer_wait     49.676 ms  29.806 Mcyc  79.5% x1 49676 us/call
```

Wall min/avg/max = 62.224/62.451/62.568 ms. `sh_rasterize` averages 12.507 ms/f; total render averages 12.775 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 1121–1136)

```
frame                    62.440 ms  37.464 Mcyc 100.0%
  pov_preserve_half      145.45 us   87.27 kcyc   0.2% x1 145 us/call
  sh_rasterize           12.493 ms   7.496 Mcyc  20.0% x1 12493 us/call
  sh_timeline_step        20.24 us   12.15 kcyc   0.0% x1 20 us/call
  canvas_clear            82.37 us   49.42 kcyc   0.1% x1 82 us/call
  canvas_buffer_wait     49.681 ms  29.808 Mcyc  79.6% x1 49681 us/call
```

Wall min/avg/max = 62.224/62.440/62.524 ms. `sh_rasterize` averages 12.493 ms/f; total render averages 12.759 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 1137–1152)

```
frame                    62.452 ms  37.471 Mcyc 100.0%
  pov_preserve_half      153.34 us   92.00 kcyc   0.2% x1 153 us/call
  sh_rasterize           12.416 ms   7.449 Mcyc  19.9% x1 12416 us/call
  sh_timeline_step        22.17 us   13.30 kcyc   0.0% x1 22 us/call
  canvas_clear            82.56 us   49.54 kcyc   0.1% x1 83 us/call
  canvas_buffer_wait     49.759 ms  29.856 Mcyc  79.7% x1 49759 us/call
```

Wall min/avg/max = 62.351/62.452/62.584 ms. `sh_rasterize` averages 12.416 ms/f; total render averages 12.693 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1521–1536)

```
frame                    62.451 ms  37.470 Mcyc 100.0%
  pov_preserve_half      145.34 us   87.20 kcyc   0.2% x1 145 us/call
  sh_rasterize           11.186 ms   6.712 Mcyc  17.9% x1 11186 us/call
  sh_timeline_step        27.58 us   16.55 kcyc   0.0% x1 28 us/call
  canvas_clear            83.13 us   49.88 kcyc   0.1% x1 83 us/call
  canvas_buffer_wait     50.989 ms  30.594 Mcyc  81.6% x1 50989 us/call
```

Wall min/avg/max = 61.900/62.450/62.880 ms. `sh_rasterize` averages 11.186 ms/f; total render averages 11.461 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Mode:n logs the target mathematical harmonic index n=l*l+l+m. It is not the PRESET_IDS index. There is no full 16-frame static hold; all intervals are64-frame morphs. Per-owner peak telemetry describes the complete incoming morph, not an isolated harmonic.

The show-scoped data covers 24 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 24 · l=4, m=4 | 3 | 1121–1136 | near-target continuous-morph anchor | — | — | 12.493 | 12.759 | 16 | 12.841 | 0/64 |
| 17 · l=4, m=-3 | 3 | 673–688 | near-target continuous-morph anchor | — | — | 12.471 | 12.740 | 16 | 12.783 | 0/64 |
| 16 · l=4, m=-4 | 3 | 609–624 | near-target continuous-morph anchor | — | — | 12.316 | 12.586 | 16 | 12.676 | 0/64 |
| 18 · l=4, m=-2 | 3 | 737–752 | near-target continuous-morph anchor | — | — | 12.316 | 12.584 | 16 | 12.720 | 0/64 |
| 23 · l=4, m=3 | 3 | 1057–1072 | near-target continuous-morph anchor | — | — | 12.267 | 12.531 | 16 | 12.599 | 0/64 |
| 22 · l=4, m=2 | 3 | 993–1008 | near-target continuous-morph anchor | — | — | 12.106 | 12.371 | 16 | 12.448 | 0/64 |
| 19 · l=4, m=-1 | 3 | 801–816 | near-target continuous-morph anchor | — | — | 12.081 | 12.345 | 16 | 12.456 | 0/64 |
| 10 · l=3, m=-2 | 3 | 225–240 | near-target continuous-morph anchor | — | — | 12.001 | 12.274 | 16 | 12.308 | 0/64 |
| 15 · l=3, m=3 | 3 | 545–560 | near-target continuous-morph anchor | — | — | 11.999 | 12.273 | 16 | 12.324 | 0/64 |
| 9 · l=3, m=-3 | 3 | 161–176 | near-target continuous-morph anchor | — | — | 11.932 | 12.200 | 16 | 12.304 | 0/64 |
| 1 · l=1, m=-1 | 3 | 1185–1200 | near-target continuous-morph anchor | — | — | 11.887 | 12.159 | 16 | 12.709 | 0/64 |
| 14 · l=3, m=2 | 3 | 481–496 | near-target continuous-morph anchor | — | — | 11.822 | 12.086 | 16 | 12.149 | 0/64 |
| 11 · l=3, m=-1 | 3 | 289–304 | near-target continuous-morph anchor | — | — | 11.804 | 12.067 | 16 | 12.171 | 0/64 |
| 20 · l=4, m=0 | 3 | 865–880 | near-target continuous-morph anchor | — | — | 11.773 | 12.036 | 16 | 12.439 | 0/64 |
| 21 · l=4, m=1 | 3 | 929–944 | near-target continuous-morph anchor | — | — | 11.712 | 11.979 | 16 | 12.121 | 0/64 |
| 12 · l=3, m=0 | 3 | 353–368 | near-target continuous-morph anchor | — | — | 11.566 | 11.834 | 16 | 12.070 | 0/64 |
| 8 · l=2, m=2 | 3 | 97–112 | near-target continuous-morph anchor | — | — | 11.560 | 11.831 | 16 | 11.895 | 0/64 |
| 5 · l=2, m=-1 | 3 | 1441–1456 | near-target continuous-morph anchor | — | — | 11.518 | 11.790 | 16 | 11.828 | 0/64 |
| 13 · l=3, m=1 | 3 | 417–432 | near-target continuous-morph anchor | — | — | 11.514 | 11.786 | 16 | 11.949 | 0/64 |
| 4 · l=2, m=-2 | 3 | 1377–1392 | near-target continuous-morph anchor | — | — | 11.472 | 11.744 | 16 | 11.853 | 0/64 |
| 7 · l=2, m=1 | 2 | 33–48 | near-target continuous-morph anchor | — | — | 11.190 | 11.453 | 16 | 11.686 | 0/63 |
| 6 · l=2, m=0 | 3 | 1505–1520 | near-target continuous-morph anchor | — | — | 11.139 | 11.410 | 16 | 11.659 | 0/64 |
| 3 · l=1, m=1 | 3 | 1313–1328 | near-target continuous-morph anchor | — | — | 10.984 | 11.254 | 16 | 11.373 | 0/64 |
| 2 · l=1, m=0 | 3 | 1249–1264 | near-target continuous-morph anchor | — | — | 10.958 | 11.230 | 16 | 12.023 | 0/64 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 1105–1120; CPU shares use its separately snapshotted ISR interval of 1000035 us.

```
isr_wake       1152.1/f 0.315/1.671/12.873 us 3.08%
  isr_pack       144.0/f 7.831/8.445/10.651 us 1.94%
  isr_dma_submit 144.0/f 0.755/0.910/1.023 us 0.20%
```

- Mean pack 8.445 us/call versus submit 0.910 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.08%, corresponding to 60.575 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.191×; Peak live show render: 0.204×; Costliest qualified preset: 0.204×; Transition / non-hold regime: 0.203×; Show end: 0.183×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `sh_rasterize` — 20.0% of the frame, 12.507 ms/f.
2. `pov_preserve_half` — 0.2% of the frame, 0.147 ms/f.
3. `canvas_clear` — 0.1% of the frame, 0.082 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
