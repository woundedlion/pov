# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Fishbowl`).
Raw capture: `C:/work/temp/astra/quest-verification/final-fishbowl.log`, captured 2026-10-23 08:39
America/Los_Angeles on COM4. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, measured source `a6222f23b579ecdddb320cc991115d4acef8aaf3` |
| Method | 60 s capture, 16-frame windows. Shipping Fishbowl configuration, live driver, no fixed-preset override. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-quest-verify HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Fishbowl profile 60 16` |

Image size (`profile`, this effect only): `FLASH: code:136848, data:152936, headers:8200   free for files:1733632` / `RAM1: variables:313024, code:33000, padding:32536   free for local variables:145728` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `ff0454ec37157150c633e376201f81fdbc2ea291dd04d2e13741d60b641eceb3`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 497–512
have root cycles/600 = 999527.407 us
against 999603 us measured wall sum (**75.6 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **31.100 ms**
(frame 197), mean live render **23.075 ms**,
**0/927 spilled frames**. Setup frame 1 is excluded.
`fish_build_vertices` averages 12.994 ms/frame across complete
post-setup windows; its worst window averages
13.646 ms/frame
(frames 497–512).

Setup frame 1 renders 1.799 ms, reported separately.
Complete counter windows end at frame 928;
9 trailing raw frame rows without a complete
counter window are excluded from the published figures. No preset markers are emitted; the capture describes the observed show segment.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: Shipping Fishbowl configuration, live driver, no fixed-preset override.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 193–208)

```text
frame                       62.446 ms 37.468 Mcyc 100.0%
  pov_preserve_half          139.1 us   83.5 kcyc  0.2%
  fish_multiline_draw       11.962 ms  7.177 Mcyc 19.2%
  fish_build_vertices       13.613 ms  8.168 Mcyc 21.8%
  fish_noise_prepare           4.1 us    2.5 kcyc  0.0%
  fish_timeline_step         249.2 us  149.5 kcyc  0.4%
    animation_noise_clock      0.3 us    0.2 kcyc  0.0%
  canvas_clear                90.4 us   54.3 kcyc  0.1%
  canvas_buffer_wait        35.723 ms 21.434 Mcyc 57.2%
```

Wall min/avg/max = 54.978/62.451/69.830 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Highest dominant-scope window (frames 497–512)

```text
frame                       62.470 ms 37.482 Mcyc 100.0%
  pov_preserve_half          140.5 us   84.3 kcyc  0.2%
  fish_multiline_draw       10.092 ms  6.055 Mcyc 16.2%
  fish_build_vertices       13.646 ms  8.188 Mcyc 21.8%
  fish_noise_prepare           4.3 us    2.6 kcyc  0.0%
  fish_timeline_step         238.8 us  143.3 kcyc  0.4%
    animation_noise_clock      0.2 us    0.1 kcyc  0.0%
  canvas_clear                90.3 us   54.2 kcyc  0.1%
  canvas_buffer_wait        37.594 ms 22.556 Mcyc 60.2%
```

Wall min/avg/max = 60.468/62.475/63.480 ms.
This window sets the highest dominant-scope average; it does not replace the higher individual render used in the headline.

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1152.6/f  0.560/2.008/29.361 us  3.70%
isr_pack          144.1/f  6.075/6.746/9.445 us  1.55%
isr_dma_submit    144.1/f  0.641/0.941/7.711 us  0.21%
```

- Packing costs 1.55% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.70%, leaving
  60.188 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `fish_build_vertices` — 21.8% of the peak window, 13.613 ms/frame.
2. `fish_multiline_draw` — 19.2% of the peak window, 11.962 ms/frame.
3. `fish_timeline_step` — 0.4% of the peak window, 0.249 ms/frame.
4. `pov_preserve_half` — 0.2% of the peak window, 0.139 ms/frame.

The preceding canonical capture recorded 28.480 ms peak.
This capture records +2.620 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 31.100, spilled 🟢 0/927 (0.00%).

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
- No preset markers are emitted; the capture describes the observed show segment. Strict validation has limitations shown below.

Parser validation:

```text
effect=Fishbowl  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (Fishbowl 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'fish_build_vertices' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (75.6 ppm, frames 497-512)

INVALID: Fishbowl
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Fishbowl` and
`HS_PROFILE_WINDOW=16`; `just profile Fishbowl` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
