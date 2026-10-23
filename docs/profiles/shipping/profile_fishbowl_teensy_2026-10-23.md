# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/fishbowl_ship.log`, captured 2026-10-23 13:20 America/Los_Angeles on COM4.
Replaces `profile_fishbowl_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Fishbowl, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Fishbowl profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:136960, data:152936, headers:9112   free for files:1732608` / `RAM1: variables:313024, code:33000, padding:32536   free for local variables:145728` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `f3aa1b8c0521d0e029a93d22cc19bff5c59e37be628300770fd7387a8d7f9173`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **31.232 ms** (frame 197), mean **23.101 ms**, **0/928 spilled frames** (0.00%).
The frame scope is 1–928; no warmup or setup frame is excluded.
Frame 1 renders 1.811 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 937 raw frame rows, including that trailing portion, the observed peak is 31.232 ms and 0 frames spill.

`fish_build_vertices` averages 12.835 ms/frame; its worst window averages 13.691 ms/frame (frames 737–768).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 58.835 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 0.97 | 0.24% |
| ↳ `fish_multiline_draw` | 2.081 | 1.00 | 3.54% |
| ↳ `fish_build_vertices` | 2.133 | 1.00 | 3.63% |
| ↳ `fish_noise_prepare` | 0.004 | 1.00 | 0.01% |
| ↳ `fish_timeline_step` | 0.229 | 1.00 | 0.39% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 53.997 | 1.00 | 91.78% |

Wall min/mean/max: 1.766/58.838/64.052 ms.

### Peak-render window (frames 193–224)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.427 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 1.00 | 0.22% |
| ↳ `fish_multiline_draw` | 12.056 | 1.00 | 19.31% |
| ↳ `fish_build_vertices` | 13.637 | 1.00 | 21.84% |
| ↳ `fish_noise_prepare` | 0.004 | 1.00 | 0.01% |
| ↳ `fish_timeline_step` | 0.240 | 1.00 | 0.38% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 35.604 | 1.00 | 57.03% |

Wall min/mean/max: 54.979/62.432/69.836 ms.

### Worst dominant-scope window (frames 737–768)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.415 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.22% |
| ↳ `fish_multiline_draw` | 10.055 | 1.00 | 16.11% |
| ↳ `fish_build_vertices` | 13.691 | 1.00 | 21.94% |
| ↳ `fish_noise_prepare` | 0.004 | 1.00 | 0.01% |
| ↳ `fish_timeline_step` | 0.230 | 1.00 | 0.37% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 37.547 | 1.00 | 60.16% |

Wall min/mean/max: 59.386/62.418/65.345 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2001.211 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1152.75 | 0.000556 | 0.001945 | 0.027708 | 3.58% |
| `isr_pack` | 144.09 | 0.006023 | 0.006748 | 0.009211 | 1.55% |
| `isr_dma_submit` | 144.09 | 0.000618 | 0.000938 | 0.004166 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fish_build_vertices` — 13.637 ms/frame, 21.84% of this root window.
2. `fish_multiline_draw` — 12.056 ms/frame, 19.31% of this root window.
3. `fish_timeline_step` — 0.240 ms/frame, 0.38% of this root window.

README cells: peak 🟢 31.232, spilled 🟢 0/928 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=Fishbowl  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Fishbowl 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fish_build_vertices' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (66.0 ppm, frames 737-768)

INVALID: Fishbowl
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (66.0 ppm, frames 737-768). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
