# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/gridspace_ship.log`, captured 2026-10-23 14:12 America/Los_Angeles on COM4.
Replaces `profile_gridspace_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridSpace, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GridSpace profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:87544, data:155216, headers:9144   free for files:1779712` / `RAM1: variables:313024, code:22008, padding:10760   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `6011d9e9825b412ce91faefdac2bdf18a5e5af96511622884b97b5f6a8c3f18e`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **38.937 ms** (frame 1), mean **24.870 ms**, **0/928 spilled frames** (0.00%).
The frame scope is 1–928; no warmup or setup frame is excluded.
Frame 1 renders 38.937 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 937 raw frame rows, including that trailing portion, the observed peak is 38.937 ms and 0 frames spill.

`fx_shader_draw` averages 18.120 ms/frame; its worst window averages 18.631 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.491 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 18.631 | 1.00 | 30.80% |
| ↳ ↳ `fx_prepare_column` | 0.498 | 144.00 | 0.82% |
| ↳ `fx_prepare_frame` | 3.701 | 1.00 | 6.12% |
| ↳ `fx_advance` | 2.222 | 1.00 | 3.67% |
| ↳ `fx_timeline_step` | 0.123 | 1.00 | 0.20% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 35.226 | 1.00 | 58.23% |

Wall min/mean/max: 24.743/60.500/63.066 ms.

### Worst post-setup scope window (frames 577–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.387 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 18.245 | 1.00 | 29.24% |
| ↳ ↳ `fx_prepare_column` | 0.495 | 144.00 | 0.79% |
| ↳ `fx_prepare_frame` | 3.902 | 1.00 | 6.25% |
| ↳ `fx_advance` | 2.128 | 1.00 | 3.41% |
| ↳ `fx_timeline_step` | 0.132 | 1.00 | 0.21% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 37.398 | 1.00 | 59.95% |

Wall min/mean/max: 59.255/62.394/65.363 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2398.872 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1381.84 | 0.000451 | 0.001658 | 0.027728 | 3.05% |
| `isr_pack` | 136.81 | 0.006808 | 0.009116 | 0.012158 | 1.66% |
| `isr_dma_submit` | 136.81 | 0.000606 | 0.000938 | 0.006333 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 18.631 ms/frame, 30.80% of this root window.
2. `fx_prepare_frame` — 3.701 ms/frame, 6.12% of this root window.
3. `fx_advance` — 2.222 ms/frame, 3.67% of this root window.

README cells: peak 🟢 38.937, spilled 🟢 0/928 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=GridSpace  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (GridSpace 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (156.1 ppm, frames 1-32)

INVALID: GridSpace
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (156.1 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
