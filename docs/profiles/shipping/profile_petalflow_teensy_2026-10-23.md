# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/petalflow_ship.log`, captured 2026-10-23 13:29 America/Los_Angeles on COM4.
Replaces `profile_petalflow_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | PetalFlow, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh PetalFlow profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:109896, data:148140, headers:8204   free for files:1765376` / `RAM1: variables:312992, code:28728, padding:4040   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `ad521a5df16c74489fe00195b82edd113e7f9eaa07396381bd71a7cbb8b8ac8f`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **30.286 ms** (frame 307), mean **24.655 ms**, **0/448 spilled frames** (0.00%).
The frame scope is 1–448; no warmup or setup frame is excluded.
Frame 1 renders 21.618 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 457 raw frame rows, including that trailing portion, the observed peak is 30.286 ms and 0 frames spill.

`pf_draw_rings` averages 24.323 ms/frame; its worst window averages 27.896 ms/frame (frames 321–352).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.847 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.133 | 0.97 | 0.22% |
| ↳ `pf_draw_rings` | 21.413 | 1.00 | 35.78% |
| ↳ ↳ `pf_ring_scan` | 20.357 | 23.34 | 34.02% |
| ↳ ↳ `pf_ring_build` | 0.887 | 23.34 | 1.48% |
| ↳ `pf_timeline_step` | 0.049 | 1.00 | 0.08% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 38.110 | 1.00 | 63.68% |

Wall min/mean/max: 20.327/59.850/66.496 ms.

### Peak-render window (frames 289–320)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.545 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.135 | 1.00 | 0.22% |
| ↳ `pf_draw_rings` | 27.853 | 1.00 | 44.53% |
| ↳ ↳ `pf_ring_scan` | 26.787 | 23.34 | 42.83% |
| ↳ ↳ `pf_ring_build` | 0.890 | 23.34 | 1.42% |
| ↳ `pf_timeline_step` | 0.054 | 1.00 | 0.09% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.363 | 1.00 | 54.94% |

Wall min/mean/max: 59.099/62.547/65.087 ms.

### Worst dominant-scope window (frames 321–352)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.405 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `pf_draw_rings` | 27.896 | 1.00 | 44.70% |
| ↳ ↳ `pf_ring_scan` | 26.819 | 23.34 | 42.97% |
| ↳ ↳ `pf_ring_build` | 0.892 | 23.34 | 1.43% |
| ↳ `pf_timeline_step` | 0.059 | 1.00 | 0.09% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.178 | 1.00 | 54.77% |

Wall min/mean/max: 59.637/62.408/64.257 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2003.165 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1153.91 | 0.000555 | 0.001690 | 0.023248 | 3.11% |
| `isr_pack` | 144.22 | 0.006081 | 0.006978 | 0.009766 | 1.60% |
| `isr_dma_submit` | 144.22 | 0.000681 | 0.000935 | 0.002841 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `pf_draw_rings` — 27.853 ms/frame, 44.53% of this root window.
2. `pov_preserve_half` — 0.135 ms/frame, 0.22% of this root window.
3. `canvas_clear` — 0.094 ms/frame, 0.15% of this root window.

README cells: peak 🟢 30.286, spilled 🟢 0/448 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=PetalFlow  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (PetalFlow 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'pf_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (50.0 ppm, frames 321-352)

INVALID: PetalFlow
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (50.0 ppm, frames 321-352). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
