# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/ashcloud_ship.log`, captured 2026-10-23 14:21 America/Los_Angeles on COM4.
Replaces `profile_ashcloud_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AshCloud profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88496, data:154660, headers:8748   free for files:1779712` / `RAM1: variables:313024, code:23304, padding:9464   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `e6170093f89e0fefd1b8f8a37433c85486220b5c43d3266e6acf4905cc42d200`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **47.180 ms** (frame 868), mean **44.551 ms**, **0/927 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–928 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **73.108 ms**, 1/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 936 raw frame rows, including that trailing portion, the observed peak is 73.108 ms and 1 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 37.983 ms/frame; its worst window averages 39.928 ms/frame (frames 545–576).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.119 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.97 | 0.22% |
| ↳ `fx_shader_draw` | 37.382 | 1.00 | 60.18% |
| ↳ ↳ `fx_prepare_column` | 0.592 | 144.00 | 0.95% |
| ↳ `fx_prepare_frame` | 3.684 | 1.00 | 5.93% |
| ↳ `fx_advance` | 2.179 | 1.00 | 3.51% |
| ↳ `fx_timeline_step` | 0.128 | 1.00 | 0.21% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 18.154 | 1.00 | 29.22% |

Wall min/mean/max: 43.063/62.136/73.108 ms.

### Worst dominant-scope window (frames 545–576)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.395 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 39.928 | 1.00 | 63.99% |
| ↳ ↳ `fx_prepare_column` | 0.604 | 144.00 | 0.97% |
| ↳ `fx_prepare_frame` | 3.701 | 1.00 | 5.93% |
| ↳ `fx_advance` | 2.184 | 1.00 | 3.50% |
| ↳ `fx_timeline_step` | 0.144 | 1.00 | 0.23% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 15.867 | 1.00 | 25.43% |

Wall min/mean/max: 60.980/62.404/63.633 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2416.615 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1392.06 | 0.000450 | 0.001703 | 0.017948 | 3.14% |
| `isr_pack` | 138.09 | 0.008873 | 0.009548 | 0.012203 | 1.74% |
| `isr_dma_submit` | 138.09 | 0.000715 | 0.000921 | 0.004616 | 0.16% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 37.382 ms/frame, 60.18% of this root window.
2. `fx_prepare_frame` — 3.684 ms/frame, 5.93% of this root window.
3. `fx_advance` — 2.179 ms/frame, 3.51% of this root window.

README cells: peak 🟢 47.180, spilled 🟢 0/927 (0.00%).

README setup frame 1: 73.108 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=AshCloud  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (AshCloud 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (157.7 ppm, frames 545-576)

INVALID: AshCloud
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (157.7 ppm, frames 545-576). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
