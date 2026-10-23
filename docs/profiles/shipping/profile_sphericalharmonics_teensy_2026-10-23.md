# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/sphericalharmonics_ship.log`, captured 2026-10-23 13:45 America/Los_Angeles on COM4.
Replaces `profile_sphericalharmonics_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | SphericalHarmonics, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 96000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh SphericalHarmonics profile 96 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:43008, data:147376, headers:8272   free for files:1832960` / `RAM1: variables:312896, code:11864, padding:20904   free for local variables:178624` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `725e9796997b6c7229cec95541ab27e2a99e46acd90b4633dc39d1855200e044`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **13.296 ms** (frame 924), mean **10.301 ms**, **0/1503 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–1504 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **17.983 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1513 raw frame rows, including that trailing portion, the observed peak is 17.983 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `sh_rasterize` averages 10.019 ms/frame; its worst window averages 12.971 ms/frame (frames 913–928).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 56.374 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 0.94 | 0.25% |
| ↳ `sh_rasterize` | 9.898 | 1.00 | 17.56% |
| ↳ `sh_timeline_step` | 0.032 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 46.203 | 1.00 | 81.96% |

Wall min/mean/max: 9.679/56.394/62.674 ms.

### Worst dominant-scope window (frames 913–928)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.435 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `sh_rasterize` | 12.971 | 1.00 | 20.78% |
| ↳ `sh_timeline_step` | 0.029 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.083 | 1.00 | 0.13% |
| ↳ `canvas_buffer_wait` | 49.203 | 1.00 | 78.81% |

Wall min/mean/max: 62.199/62.454/62.562 ms.

### Ownership-change window (frames 881–896)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.431 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.149 | 1.00 | 0.24% |
| ↳ `sh_rasterize` | 12.899 | 1.00 | 20.66% |
| ↳ `sh_timeline_step` | 0.031 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.083 | 1.00 | 0.13% |
| ↳ `canvas_buffer_wait` | 49.263 | 1.00 | 78.91% |

Wall min/mean/max: 62.354/62.450/62.547 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `sh_rasterize` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 21 | — | 13.296 | 0/64 | 12.971 | 3/4 |
| 20 | — | 13.284 | 0/64 | 12.970 | 3/4 |
| 22 | — | 13.258 | 0/64 | 11.843 | 3/4 |
| 19 | — | 12.136 | 0/64 | 11.822 | 3/4 |
| 23 | — | 12.079 | 0/64 | 10.499 | 3/4 |
| 13 | — | 11.619 | 0/64 | 11.259 | 3/4 |
| 12 | — | 11.618 | 0/64 | 11.295 | 3/4 |
| 14 | — | 11.535 | 0/64 | 10.053 | 3/4 |
| 18 | — | 10.814 | 0/64 | 10.504 | 3/4 |
| 24 | — | 10.746 | 0/64 | 9.812 | 3/4 |
| 11 | — | 10.356 | 0/64 | 10.017 | 3/4 |
| 15 | — | 10.255 | 0/64 | 9.340 | 3/4 |
| 17 | — | 10.127 | 0/64 | 9.804 | 3/4 |
| 1 | — | 10.052 | 0/64 | 9.018 | 3/4 |
| 16 | — | 9.823 | 0/64 | 9.487 | 3/4 |
| 7 | — | 9.822 | 0/62 | 9.898 | 3/4 |
| 6 | — | 9.808 | 0/33 | 9.448 | 2/2 |
| 8 | — | 9.767 | 0/64 | 8.885 | 3/4 |
| 10 | — | 9.685 | 0/64 | 9.350 | 3/4 |
| 9 | — | 9.361 | 0/64 | 9.048 | 3/4 |
| 5 | — | 9.210 | 0/64 | 8.870 | 3/4 |
| 2 | — | 9.208 | 0/64 | 8.325 | 3/4 |
| 4 | — | 8.902 | 0/64 | 8.580 | 3/4 |
| 3 | — | 8.673 | 0/64 | 8.337 | 3/4 |

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 1383.689 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1594.12 | 0.000428 | 0.001450 | 0.017696 | 2.67% |
| `isr_pack` | 127.44 | 0.008766 | 0.009228 | 0.011550 | 1.35% |
| `isr_dma_submit` | 127.44 | 0.000683 | 0.000931 | 0.001053 | 0.13% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `sh_rasterize` — 9.898 ms/frame, 17.56% of this root window.
2. `pov_preserve_half` — 0.143 ms/frame, 0.25% of this root window.
3. `canvas_clear` — 0.092 ms/frame, 0.16% of this root window.

README cells: peak 🟢 13.296 (24), spilled 🟢 0/1503 (0.00%).

README setup frame 1: 17.983 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=SphericalHarmonics  windows=94
  [PASS] captured >=3 windows (94)
  [PASS] every window names one effect and resolution (SphericalHarmonics 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 24, 24 distinct)
  [PASS] all 24 presets visited (24)
  [PASS] per-frame render is render, not wall (0 of 94 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 94 windows are short)
  [PASS] scope 'sh_rasterize' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (321.5 ppm, frames 913-928)

INVALID: SphericalHarmonics
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (321.5 ppm, frames 913-928). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
