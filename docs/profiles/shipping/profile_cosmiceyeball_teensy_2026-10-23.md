# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/cosmiceyeball_ship.log`, captured 2026-10-23 14:36 America/Los_Angeles on COM4.
Replaces `profile_cosmiceyeball_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | CosmicEyeball, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh CosmicEyeball profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:85792, data:154328, headers:8712   free for files:1782784` / `RAM1: variables:313024, code:22344, padding:10424   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `df5b352955154a8dee606015d46134dd6178995d47e356cfaa75e5aca9bb0588`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **26.390 ms** (frame 456), mean **26.074 ms**, **0/927 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–928 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **48.097 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 937 raw frame rows, including that trailing portion, the observed peak is 48.097 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 23.072 ms/frame; its worst window averages 23.761 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.864 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.137 | 0.97 | 0.22% |
| ↳ `fx_shader_draw` | 23.761 | 1.00 | 39.04% |
| ↳ ↳ `fx_prepare_column` | 0.564 | 144.00 | 0.93% |
| ↳ `fx_prepare_frame` | 0.592 | 1.00 | 0.97% |
| ↳ `fx_advance` | 1.952 | 1.00 | 3.21% |
| ↳ `fx_timeline_step` | 0.090 | 1.00 | 0.15% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.056 | 1.00 | 55.95% |

Wall min/mean/max: 26.097/60.869/62.872 ms.

### Worst post-setup scope window (frames 577–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.440 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 23.081 | 1.00 | 36.97% |
| ↳ ↳ `fx_prepare_column` | 0.563 | 144.00 | 0.90% |
| ↳ `fx_prepare_frame` | 0.581 | 1.00 | 0.93% |
| ↳ `fx_advance` | 1.854 | 1.00 | 2.97% |
| ↳ `fx_timeline_step` | 0.088 | 1.00 | 0.14% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 36.424 | 1.00 | 58.34% |

Wall min/mean/max: 60.870/62.444/63.871 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2400.113 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1382.56 | 0.000511 | 0.001710 | 0.031788 | 3.15% |
| `isr_pack` | 136.91 | 0.008778 | 0.009545 | 0.011933 | 1.74% |
| `isr_dma_submit` | 136.91 | 0.000660 | 0.000933 | 0.005950 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 23.761 ms/frame, 39.04% of this root window.
2. `fx_advance` — 1.952 ms/frame, 3.21% of this root window.
3. `fx_prepare_frame` — 0.592 ms/frame, 0.97% of this root window.

README cells: peak 🟢 26.390, spilled 🟢 0/927 (0.00%).

README setup frame 1: 48.097 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=CosmicEyeball  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (CosmicEyeball 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (86.3 ppm, frames 1-32)

INVALID: CosmicEyeball
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (86.3 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
