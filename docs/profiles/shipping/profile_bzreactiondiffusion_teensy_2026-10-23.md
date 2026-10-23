# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/bzreactiondiffusion_ship.log`, captured 2026-10-23 13:19 America/Los_Angeles on COM4.
Replaces `profile_bzreactiondiffusion_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | BZReactionDiffusion, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh BZReactionDiffusion profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:68320, data:241624, headers:8520   free for files:1713152` / `RAM1: variables:313024, code:16168, padding:16600   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `ac2808e65d7f5c92b1b0e56c7eb7dc861de32f793a2c55d8fc59c64f008b1e80`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **78.649 ms** (frame 1), mean **46.557 ms**, **1/928 spilled frames** (0.11%).
The frame scope is 1–928; no warmup or setup frame is excluded.
Frame 1 renders 78.649 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 936 raw frame rows, including that trailing portion, the observed peak is 78.649 ms and 1 frames spill.

`bz_render` averages 46.270 ms/frame; its worst window averages 47.139 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.476 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.131 | 0.97 | 0.21% |
| ↳ `bz_render` | 47.139 | 1.00 | 75.45% |
| ↳ ↳ `bz_raster` | 42.483 | 2.00 | 68.00% |
| ↳ ↳ `bz_physics` | 4.562 | 2.00 | 7.30% |
| ↳ ↳ `rd_lattice_basis` | 0.091 | 1.00 | 0.15% |
| ↳ `rd_timeline_step` | 0.043 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.068 | 1.00 | 24.12% |

Wall min/mean/max: 44.214/62.500/78.650 ms.

### Worst post-setup scope window (frames 33–64)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.417 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.135 | 1.00 | 0.22% |
| ↳ `bz_render` | 46.938 | 1.00 | 75.20% |
| ↳ ↳ `bz_raster` | 42.314 | 2.00 | 67.79% |
| ↳ ↳ `bz_physics` | 4.530 | 2.00 | 7.26% |
| ↳ ↳ `rd_lattice_basis` | 0.090 | 1.00 | 0.14% |
| ↳ `rd_timeline_step` | 0.039 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.086 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.215 | 1.00 | 24.38% |

Wall min/mean/max: 61.365/62.441/63.540 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2421.532 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1394.88 | 0.000471 | 0.001573 | 0.017423 | 2.90% |
| `isr_pack` | 138.44 | 0.005990 | 0.008360 | 0.012590 | 1.52% |
| `isr_dma_submit` | 138.44 | 0.000630 | 0.000933 | 0.002856 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `bz_render` — 47.139 ms/frame, 75.45% of this root window.
2. `pov_preserve_half` — 0.131 ms/frame, 0.21% of this root window.
3. `canvas_clear` — 0.090 ms/frame, 0.14% of this root window.

README cells: peak 🔴 78.649, spilled 🔴 1/928 (0.11%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=BZReactionDiffusion  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (BZReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'bz_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (389.6 ppm, frames 1-32)

INVALID: BZReactionDiffusion
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (389.6 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
