# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/kaleidoscopepentbright_ship.log`, captured 2026-10-23 14:23 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopepentbright_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopePentBright, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopePentBright profile 60 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:86936, data:155040, headers:8904   free for files:1780736` / `RAM1: variables:313024, code:21768, padding:11000   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `252e869b841b601fd7c5e12e11122d44c1c610586a49a4d436afb6ed55fda6e6`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Column strobe is off in this shipping configuration.

## Frame cadence

**Steady render, frame 1 excluded:** peak **27.973 ms** (frame 874), mean **26.174 ms**, **0/927 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–928 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **48.788 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 937 raw frame rows, including that trailing portion, the observed peak is 48.788 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 22.596 ms/frame; its worst window averages 23.908 ms/frame (frames 417–448).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.834 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 22.586 | 1.00 | 37.13% |
| ↳ ↳ `fx_prepare_column` | 0.497 | 144.00 | 0.82% |
| ↳ `fx_prepare_frame` | 1.194 | 1.00 | 1.96% |
| ↳ `fx_advance` | 1.867 | 1.00 | 3.07% |
| ↳ `fx_timeline_step` | 0.110 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 34.465 | 1.00 | 56.65% |

Wall min/mean/max: 26.234/60.842/63.321 ms.

### Worst dominant-scope window (frames 417–448)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.390 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.145 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 23.908 | 1.00 | 38.32% |
| ↳ ↳ `fx_prepare_column` | 0.495 | 144.00 | 0.79% |
| ↳ `fx_prepare_frame` | 0.886 | 1.00 | 1.42% |
| ↳ `fx_advance` | 1.940 | 1.00 | 3.11% |
| ↳ `fx_timeline_step` | 0.116 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.937 | 1.00 | 56.00% |

Wall min/mean/max: 61.712/62.397/63.137 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2399.598 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1382.28 | 0.000451 | 0.001695 | 0.025621 | 3.12% |
| `isr_pack` | 136.88 | 0.008790 | 0.009566 | 0.017390 | 1.74% |
| `isr_dma_submit` | 136.88 | 0.000750 | 0.000923 | 0.004870 | 0.16% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 22.586 ms/frame, 37.13% of this root window.
2. `fx_advance` — 1.867 ms/frame, 3.07% of this root window.
3. `fx_prepare_frame` — 1.194 ms/frame, 1.96% of this root window.

README cells: peak 🟢 27.973, spilled 🟢 0/927 (0.00%).

README setup frame 1: 48.788 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=KaleidoscopePentBright  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (KaleidoscopePentBright 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (110.7 ppm, frames 417-448)

INVALID: KaleidoscopePentBright
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (110.7 ppm, frames 417-448). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopePentBright`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
