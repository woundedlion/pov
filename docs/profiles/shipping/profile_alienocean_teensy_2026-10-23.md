# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/alienocean_ship.log`, captured 2026-10-23 14:07 America/Los_Angeles on COM4.
Replaces `profile_alienocean_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienOcean, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienOcean profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:86800, data:154480, headers:8576   free for files:1781760` / `RAM1: variables:313024, code:22680, padding:10088   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `b65989cbc9fa37565c6b62d6a5632fe982ed64a11d13c69da9b1c9ede7b6e017`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **44.614 ms** (frame 1), mean **23.856 ms**, **0/448 spilled frames** (0.00%).
The frame scope is 1–448; no warmup or setup frame is excluded.
Frame 1 renders 44.614 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 457 raw frame rows, including that trailing portion, the observed peak is 44.614 ms and 0 frames spill.

`fx_shader_draw` averages 20.334 ms/frame; its worst window averages 20.811 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.637 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 20.811 | 1.00 | 34.32% |
| ↳ ↳ `fx_prepare_column` | 0.577 | 144.00 | 0.95% |
| ↳ `fx_prepare_frame` | 0.699 | 1.00 | 1.15% |
| ↳ `fx_advance` | 2.234 | 1.00 | 3.68% |
| ↳ `fx_timeline_step` | 0.115 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.100 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 36.174 | 1.00 | 59.66% |

Wall min/mean/max: 23.590/60.640/62.757 ms.

### Worst post-setup scope window (frames 161–192)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.388 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 20.363 | 1.00 | 32.64% |
| ↳ ↳ `fx_prepare_column` | 0.580 | 144.00 | 0.93% |
| ↳ `fx_prepare_frame` | 0.605 | 1.00 | 0.97% |
| ↳ `fx_advance` | 2.165 | 1.00 | 3.47% |
| ↳ `fx_timeline_step` | 0.112 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 38.548 | 1.00 | 61.79% |

Wall min/mean/max: 62.197/62.391/62.543 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2397.673 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1381.12 | 0.000543 | 0.001726 | 0.030136 | 3.18% |
| `isr_pack` | 136.72 | 0.008988 | 0.009560 | 0.012071 | 1.74% |
| `isr_dma_submit` | 136.72 | 0.000696 | 0.000923 | 0.006636 | 0.16% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 20.811 ms/frame, 34.32% of this root window.
2. `fx_advance` — 2.234 ms/frame, 3.68% of this root window.
3. `fx_prepare_frame` — 0.699 ms/frame, 1.15% of this root window.

README cells: peak 🟢 44.614, spilled 🟢 0/448 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=AlienOcean  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienOcean 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (58.9 ppm, frames 1-32)

INVALID: AlienOcean
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (58.9 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
