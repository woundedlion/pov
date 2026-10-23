# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/raymarch_ship.log`, captured 2026-10-23 13:30 America/Los_Angeles on COM4.
Replaces `profile_raymarch_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Raymarch, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Raymarch profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:112296, data:192288, headers:8760   free for files:1718272` / `RAM1: variables:312896, code:25032, padding:7736   free for local variables:178624` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `702a3dcff83d0e5453935233dd28d20f09c67278f8e990cea0f1b43cdc925e9f`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **57.751 ms** (frame 10), mean **52.159 ms**, **0/447 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–448 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **103.890 ms**, 1/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 456 raw frame rows, including that trailing portion, the observed peak is 103.890 ms and 1 frames spill.

Across all complete counter windows, including the setup-containing window, `rm_shader_draw` averages 48.704 ms/frame; its worst window averages 50.559 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 63.290 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.129 | 0.97 | 0.20% |
| ↳ `rm_shader_draw` | 50.559 | 1.00 | 79.89% |
| ↳ `rm_timeline_step` | 0.421 | 1.00 | 0.67% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 9.178 | 1.00 | 14.50% |

Wall min/mean/max: 53.132/63.291/103.891 ms.

### Worst post-setup scope window (frames 321–352)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.337 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.133 | 1.00 | 0.21% |
| ↳ `rm_shader_draw` | 49.649 | 1.00 | 79.65% |
| ↳ `rm_timeline_step` | 0.419 | 1.00 | 0.67% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 9.125 | 1.00 | 14.64% |

Wall min/mean/max: 54.381/62.339/70.951 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2421.911 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1395.09 | 0.000463 | 0.001566 | 0.018028 | 2.88% |
| `isr_pack` | 138.50 | 0.005990 | 0.007901 | 0.011105 | 1.44% |
| `isr_dma_submit` | 138.50 | 0.000650 | 0.000936 | 0.002830 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `rm_shader_draw` — 50.559 ms/frame, 79.89% of this root window.
2. `rm_timeline_step` — 0.421 ms/frame, 0.67% of this root window.
3. `pov_preserve_half` — 0.129 ms/frame, 0.20% of this root window.

README cells: peak 🟢 57.751, spilled 🟢 0/447 (0.00%).

README setup frame 1: 103.890 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=Raymarch  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (Raymarch 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'rm_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (36.3 ppm, frames 1-32)

INVALID: Raymarch
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (36.3 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
