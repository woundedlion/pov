# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/chromaticlichen_ship.log`, captured 2026-10-23 14:19 America/Los_Angeles on COM4.
Replaces `profile_chromaticlichen_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ChromaticLichen, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ChromaticLichen profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88816, data:154524, headers:8564   free for files:1779712` / `RAM1: variables:313024, code:24328, padding:8440   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `a8151b1ae7f9c21ab3e050d539e5f6a116266601d0187aa3c63d2aa65724a032`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **65.897 ms** (frame 1), mean **39.428 ms**, **1/448 spilled frames** (0.22%).
The frame scope is 1–448; no warmup or setup frame is excluded.
Frame 1 renders 65.897 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 456 raw frame rows, including that trailing portion, the observed peak is 65.897 ms and 1 frames spill.

`fx_shader_draw` averages 32.628 ms/frame; its worst window averages 33.458 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 61.790 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 33.458 | 1.00 | 54.15% |
| ↳ ↳ `fx_prepare_column` | 0.798 | 144.00 | 1.29% |
| ↳ `fx_prepare_frame` | 4.039 | 1.00 | 6.54% |
| ↳ `fx_advance` | 1.950 | 1.00 | 3.16% |
| ↳ `fx_timeline_step` | 0.116 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 21.641 | 1.00 | 35.02% |

Wall min/mean/max: 39.587/61.812/65.899 ms.

### Worst post-setup scope window (frames 193–224)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.395 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 32.620 | 1.00 | 52.28% |
| ↳ ↳ `fx_prepare_column` | 0.809 | 144.00 | 1.30% |
| ↳ `fx_prepare_frame` | 4.324 | 1.00 | 6.93% |
| ↳ `fx_advance` | 1.855 | 1.00 | 2.97% |
| ↳ `fx_timeline_step` | 0.121 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 22.893 | 1.00 | 36.69% |

Wall min/mean/max: 61.765/62.416/62.973 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2413.299 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1390.12 | 0.000483 | 0.001676 | 0.021635 | 3.09% |
| `isr_pack` | 137.84 | 0.008198 | 0.009070 | 0.012286 | 1.65% |
| `isr_dma_submit` | 137.84 | 0.000616 | 0.000933 | 0.008006 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 33.458 ms/frame, 54.15% of this root window.
2. `fx_prepare_frame` — 4.039 ms/frame, 6.54% of this root window.
3. `fx_advance` — 1.950 ms/frame, 3.16% of this root window.

README cells: peak 🔴 65.897, spilled 🔴 1/448 (0.22%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=ChromaticLichen  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (ChromaticLichen 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (358.9 ppm, frames 1-32)

INVALID: ChromaticLichen
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (358.9 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
