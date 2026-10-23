# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/ringspin_ship.log`, captured 2026-10-23 13:32 America/Los_Angeles on COM4.
Replaces `profile_ringspin_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | RingSpin, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh RingSpin profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:91592, data:150960, headers:8328   free for files:1780736` / `RAM1: variables:312960, code:27448, padding:5320   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `5780221f2f740c9f0e9d202816fd6daeb8269581b0c572f19c22c774146665b6`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **53.743 ms** (frame 761), mean **38.705 ms**, **0/928 spilled frames** (0.00%).
The frame scope is 1–928; no warmup or setup frame is excluded.
Frame 1 renders 3.396 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 936 raw frame rows, including that trailing portion, the observed peak is 53.743 ms and 0 frames spill.

`rs_draw_rings` averages 38.369 ms/frame; its worst window averages 41.480 ms/frame (frames 737–768).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.494 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.134 | 0.97 | 0.23% |
| ↳ `rs_draw_rings` | 25.090 | 1.00 | 42.17% |
| ↳ ↳ `rs_ring_scan` | 24.297 | 54.62 | 40.84% |
| ↳ `rs_timeline_step` | 0.094 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 34.075 | 1.00 | 57.28% |

Wall min/mean/max: 3.396/59.496/89.088 ms.

### Peak-render window (frames 737–768)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.599 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.22% |
| ↳ `rs_draw_rings` | 41.480 | 1.00 | 66.26% |
| ↳ ↳ `rs_ring_scan` | 40.317 | 76.00 | 64.40% |
| ↳ `rs_timeline_step` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_clear` | 0.088 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 20.793 | 1.00 | 33.22% |

Wall min/mean/max: 43.638/62.601/81.166 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2005.237 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1155.09 | 0.000598 | 0.001895 | 0.020488 | 3.49% |
| `isr_pack` | 144.38 | 0.006110 | 0.008633 | 0.011740 | 1.98% |
| `isr_dma_submit` | 144.38 | 0.000603 | 0.000930 | 0.001128 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `rs_draw_rings` — 41.480 ms/frame, 66.26% of this root window.
2. `pov_preserve_half` — 0.141 ms/frame, 0.22% of this root window.
3. `rs_timeline_step` — 0.090 ms/frame, 0.14% of this root window.

README cells: peak 🟢 53.743, spilled 🟢 0/928 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=RingSpin  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (RingSpin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'rs_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (35.8 ppm, frames 737-768)

INVALID: RingSpin
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (35.8 ppm, frames 737-768). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
