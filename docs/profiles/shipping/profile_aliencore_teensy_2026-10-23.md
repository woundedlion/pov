# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/aliencore_ship.log`, captured 2026-10-23 14:08 America/Los_Angeles on COM4.
Replaces `profile_aliencore_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienCore, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienCore profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:86760, data:154488, headers:8608   free for files:1781760` / `RAM1: variables:313024, code:22584, padding:10184   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `bb5eed718f93c62ce8e68c1a8e8d2fc739871fe5fa8b114cddc01bc1202c5429`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **21.483 ms** (frame 368), mean **21.134 ms**, **0/447 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–448 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **38.170 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 457 raw frame rows, including that trailing portion, the observed peak is 38.170 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 17.465 ms/frame; its worst window averages 17.895 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.347 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 17.895 | 1.00 | 29.65% |
| ↳ ↳ `fx_prepare_column` | 0.561 | 144.00 | 0.93% |
| ↳ `fx_prepare_frame` | 0.822 | 1.00 | 1.36% |
| ↳ `fx_advance` | 2.140 | 1.00 | 3.55% |
| ↳ `fx_timeline_step` | 0.094 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.100 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 38.761 | 1.00 | 64.23% |

Wall min/mean/max: 20.999/60.370/62.701 ms.

### Worst post-setup scope window (frames 321–352)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.384 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 17.453 | 1.00 | 27.98% |
| ↳ ↳ `fx_prepare_column` | 0.568 | 144.00 | 0.91% |
| ↳ `fx_prepare_frame` | 0.953 | 1.00 | 1.53% |
| ↳ `fx_advance` | 2.092 | 1.00 | 3.35% |
| ↳ `fx_timeline_step` | 0.090 | 1.00 | 0.15% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 41.180 | 1.00 | 66.01% |

Wall min/mean/max: 62.128/62.408/62.678 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2394.945 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1379.56 | 0.000543 | 0.001721 | 0.024030 | 3.17% |
| `isr_pack` | 136.53 | 0.008975 | 0.009506 | 0.012016 | 1.73% |
| `isr_dma_submit` | 136.53 | 0.000663 | 0.000938 | 0.006603 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 17.895 ms/frame, 29.65% of this root window.
2. `fx_advance` — 2.140 ms/frame, 3.55% of this root window.
3. `fx_prepare_frame` — 0.822 ms/frame, 1.36% of this root window.

README cells: peak 🟢 21.483, spilled 🟢 0/447 (0.00%).

README setup frame 1: 38.170 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=AlienCore  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienCore 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (387.4 ppm, frames 1-32)

INVALID: AlienCore
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (387.4 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
