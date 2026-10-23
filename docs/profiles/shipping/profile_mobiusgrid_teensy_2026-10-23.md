# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/mobiusgrid_ship.log`, captured 2026-10-23 13:28 America/Los_Angeles on COM4.
Replaces `profile_mobiusgrid_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MobiusGrid, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 105000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MobiusGrid profile 105 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:89568, data:155664, headers:8720   free for files:1777664` / `RAM1: variables:313024, code:21832, padding:10936   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `96955b39912cf70d67afde32c672314ada19819e5a4c4818e5e7361526223abe`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **42.938 ms** (frame 1), mean **26.714 ms**, **0/1648 spilled frames** (0.00%).
The frame scope is 1–1648; no warmup or setup frame is excluded.
Frame 1 renders 42.938 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1657 raw frame rows, including that trailing portion, the observed peak is 42.938 ms and 0 frames spill.

`fx_shader_draw` averages 23.305 ms/frame; its worst window averages 24.054 ms/frame (frames 1–16).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 58.946 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.94 | 0.23% |
| ↳ `fx_shader_draw` | 24.054 | 1.00 | 40.81% |
| ↳ ↳ `fx_prepare_column` | 0.686 | 144.00 | 1.16% |
| ↳ `fx_prepare_frame` | 0.603 | 1.00 | 1.02% |
| ↳ `fx_advance` | 2.261 | 1.00 | 3.83% |
| ↳ `fx_timeline_step` | 0.118 | 1.00 | 0.20% |
| ↳ `canvas_clear` | 0.102 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 31.473 | 1.00 | 53.39% |

Wall min/mean/max: 26.564/58.952/62.794 ms.

### Worst post-setup scope window (frames 1041–1056)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.386 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 23.989 | 1.00 | 38.45% |
| ↳ ↳ `fx_prepare_column` | 1.338 | 144.00 | 2.14% |
| ↳ `fx_prepare_frame` | 0.606 | 1.00 | 0.97% |
| ↳ `fx_advance` | 2.209 | 1.00 | 3.54% |
| ↳ `fx_timeline_step` | 0.196 | 1.00 | 0.31% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.993 | 1.00 | 56.09% |

Wall min/mean/max: 61.944/62.391/62.695 ms.

### Ownership-change window (frames 593–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.436 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 23.341 | 1.00 | 37.38% |
| ↳ ↳ `fx_prepare_column` | 1.026 | 144.00 | 1.64% |
| ↳ `fx_prepare_frame` | 0.647 | 1.00 | 1.04% |
| ↳ `fx_advance` | 2.105 | 1.00 | 3.37% |
| ↳ `fx_timeline_step` | 0.163 | 1.00 | 0.26% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 35.769 | 1.00 | 57.29% |

Wall min/mean/max: 57.975/62.442/66.670 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `fx_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 42.938 | 0/599 | 24.054 | 37/38 |
| 2 | — | 27.761 | 0/1049 | 23.989 | 65/65 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 1400.498 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1613.50 | 0.000511 | 0.001481 | 0.016355 | 2.73% |
| `isr_pack` | 129.88 | 0.008245 | 0.009365 | 0.011783 | 1.38% |
| `isr_dma_submit` | 129.88 | 0.000801 | 0.000933 | 0.004393 | 0.13% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 24.054 ms/frame, 40.81% of this root window.
2. `fx_advance` — 2.261 ms/frame, 3.83% of this root window.
3. `fx_prepare_frame` — 0.603 ms/frame, 1.02% of this root window.

README cells: peak 🟢 42.938 (2), spilled 🟢 0/1648 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=MobiusGrid  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (MobiusGrid 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (109.0 ppm, frames 1-16)

INVALID: MobiusGrid
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (109.0 ppm, frames 1-16). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
