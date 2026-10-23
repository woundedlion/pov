# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `C:/work/temp/astra/quest-verification/final-mindsplatter.log`, captured 2026-10-23 08:36
America/Los_Angeles on COM4. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, measured source `a6222f23b579ecdddb320cc991115d4acef8aaf3` |
| Method | 83 s capture, 16-frame windows. Shipping particle/trail configuration. Preset-adoption markers establish observed ownership; a complete wrap must be checked independently. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-quest-verify HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MindSplatter profile 83 16` |

Image size (`profile`, this effect only): `FLASH: code:151976, data:550400, headers:8280   free for files:1320960` / `RAM1: variables:313152, code:48472, padding:17064   free for local variables:145600` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `a5256b9163f5b5b39b0f40d51161d02fc505fcc82277716c5efd495adf15d31d`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 993–1008
have root cycles/600 = 996425.320 us
against 996488 us measured wall sum (**62.9 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **59.311 ms**
(frame 1003), mean live render **28.073 ms**,
**0/1295 spilled frames**. Setup frame 1 is excluded.
`msp_draw_particles` averages 25.448 ms/frame across complete
post-setup windows; its worst window averages
51.001 ms/frame
(frames 993–1008).

Setup frame 1 renders 0.407 ms, reported separately.
Complete counter windows end at frame 1296;
8 trailing raw frame rows without a complete
counter window are excluded from the published figures. Preset/shape wrap is observed. There are 8 parser-owned entries in complete-window live telemetry.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: Shipping particle/trail configuration. Preset-adoption markers establish observed ownership; a complete wrap must be checked independently.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 993–1008)

```text
frame                       62.277 ms 37.366 Mcyc 100.0%
  pov_preserve_half          137.0 us   82.2 kcyc  0.2%
  msp_draw_particles        51.001 ms 30.601 Mcyc 81.9%
  msp_particle_step          2.975 ms  1.785 Mcyc  4.8%
  msp_timeline_step          129.7 us   77.8 kcyc  0.2%
  canvas_clear                96.5 us   57.9 kcyc  0.2%
  canvas_buffer_wait         7.894 ms  4.737 Mcyc 12.7%
```

Wall min/avg/max = 58.053/62.280/66.115 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Preset-adoption boundary (frames 945–960)

```text
frame                       62.005 ms 37.203 Mcyc 100.0%
  pov_preserve_half          136.4 us   81.8 kcyc  0.2%
  msp_draw_particles        37.166 ms 22.299 Mcyc 59.9%
  msp_particle_step          7.781 ms  4.669 Mcyc 12.5%
  msp_timeline_step          116.6 us   69.9 kcyc  0.2%
  canvas_clear                92.7 us   55.6 kcyc  0.1%
  canvas_buffer_wait        16.663 ms  9.998 Mcyc 26.9%
```

Wall min/avg/max = 54.329/62.007/68.154 ms.
This window mixes owners across a preset adoption and is excluded from the single-owner window figures. Its live frames remain in aggregate and ownership peak/spill counts.

### Per-preset table

Preset/shape wrap is observed. There are 8 parser-owned entries in complete-window live telemetry.
Ownership includes the transition following an adoption. The window figures
use the highest-scope post-setup window with one owner; they are not certified
full-opacity clean holds. Where adoption timing gives no separately identified
hold, that distinction remains explicit. Window counts exclude advance-straddling
windows; per-entry peak/spill figures include every owned live frame.
No mesh V/E/F/I metadata is logged for these entries.

| Entry | Blended px/f | `msp_draw_particles` ms/f | Render ms/f | Wall fps | Windows | Owned peak ms | Spilled/frames |
|---|--:|--:|--:|--:|--:|--:|--:|
| 7 | n/a | 51.001 | 54.382 | 16.06 | 9 | 59.311 | 0/159 |
| 2 | n/a | 41.590 | 44.255 | 16.30 | 9 | 51.136 | 0/159 |
| 6 | n/a | 39.623 | 47.782 | 16.06 | 9 | 53.740 | 0/159 |
| 4 | n/a | 33.330 | 36.758 | 16.03 | 9 | 50.348 | 0/159 |
| 1 | n/a | 31.536 | 34.011 | 16.02 | 9 | 35.987 | 0/182 |
| 3 | n/a | 29.422 | 31.920 | 16.00 | 9 | 33.116 | 0/159 |
| 8 | n/a | 28.216 | 30.713 | 15.98 | 9 | 34.321 | 0/159 |
| 5 | n/a | 27.212 | 29.743 | 16.15 | 9 | 37.956 | 0/159 |

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1150.6/f  0.615/1.983/27.890 us  3.65%
isr_pack          143.8/f  6.436/8.503/11.601 us  1.95%
isr_dma_submit    143.8/f  0.755/0.943/6.813 us  0.21%
```

- Packing costs 1.95% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.65%, leaving
  60.219 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `msp_draw_particles` — 81.9% of the peak window, 51.001 ms/frame.
2. `msp_particle_step` — 4.8% of the peak window, 2.975 ms/frame.
3. `pov_preserve_half` — 0.2% of the peak window, 0.137 ms/frame.
4. `msp_timeline_step` — 0.2% of the peak window, 0.130 ms/frame.

The preceding canonical capture recorded 58.709 ms peak.
This capture records +0.602 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 59.311 (8), spilled 🟢 0/1295 (0.00%).

## Caveats

- CYCCNT free-runs, so every cycle scope includes ISR time.
- Counter scopes retain measurement overhead. This is not a per-pixel breakdown
  image; enabling per-pixel scopes changes the measured workload.
- `filter_blend` uses first-entry parenting where enabled and can be hidden
  when its recorded parent has no calls; absence is not a zero-blend count.
- Shipping selective-O3 Plot, raster and effect hot regions remain active.
  These are not global-O3 reference images.
- Exact selection, transition-speed and epoch flags are in Method and envdump.
  Epoch stretching prevents teardown; unobserved phases remain unmeasured.
- Measured source is `a6222f23b579ecdddb320cc991115d4acef8aaf3`. Provenance identifies the
  actual firmware; later parity/test-only commits are not attributed to it.
- Preset/shape wrap is observed. There are 8 parser-owned entries in complete-window live telemetry. Strict validation has limitations shown below.

Parser validation:

```text
effect=MindSplatter  windows=81
  [PASS] captured >=3 windows (81)
  [PASS] every window names one effect and resolution (MindSplatter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 8, 8 distinct)
  [PASS] all 8 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 81 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 81 windows are short)
  [PASS] scope 'msp_draw_particles' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (62.9 ppm, frames 993-1008)

INVALID: MindSplatter
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter` and
`HS_PROFILE_WINDOW=16`; `just profile MindSplatter` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
