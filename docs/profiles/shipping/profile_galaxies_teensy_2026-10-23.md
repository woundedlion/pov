# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Galaxies`).
Raw capture: `C:/work/temp/astra/quest-verification/final-galaxies.log`, captured 2026-10-23 08:41
America/Los_Angeles on COM4. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Galaxies 288×144, single-entry playlist, measured source `a6222f23b579ecdddb320cc991115d4acef8aaf3` |
| Method | 60 s capture, 16-frame windows. Shipping Galaxies configuration, live driver, no fixed-preset override. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-quest-verify HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Galaxies profile 60 16` |

Image size (`profile`, this effect only): `FLASH: code:91064, data:152588, headers:8252   free for files:1779712` / `RAM1: variables:312960, code:23544, padding:9224   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `55eda3f68ebc2fa79881fef7b0a4e2e5bad350be33785cdb36da77b9b498f7a0`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 865–880
have root cycles/600 = 998667.338 us
against 998690 us measured wall sum (**22.7 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **45.496 ms**
(frame 883), mean live render **25.109 ms**,
**0/927 spilled frames**. Setup frame 1 is excluded.
`gx_particle_step` averages 15.776 ms/frame across complete
post-setup windows; its worst window averages
27.105 ms/frame
(frames 865–880).

Setup frame 1 renders 1.483 ms, reported separately.
Complete counter windows end at frame 928;
8 trailing raw frame rows without a complete
counter window are excluded from the published figures. No preset markers are emitted; the capture describes the observed show segment.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: Shipping Galaxies configuration, live driver, no fixed-preset override.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 881–896)

```text
frame                       62.300 ms 37.380 Mcyc 100.0%
  pov_preserve_half          146.3 us   87.8 kcyc  0.2%
  gx_draw_particles         16.716 ms 10.030 Mcyc 26.8%
  gx_particle_step          27.094 ms 16.256 Mcyc 43.5%
  gx_timeline_step           116.8 us   70.1 kcyc  0.2%
  canvas_clear                93.7 us   56.2 kcyc  0.2%
  canvas_buffer_wait        18.084 ms 10.850 Mcyc 29.0%
```

Wall min/avg/max = 60.579/62.303/63.907 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Highest dominant-scope window (frames 865–880)

```text
frame                       62.417 ms 37.450 Mcyc 100.0%
  pov_preserve_half          142.8 us   85.7 kcyc  0.2%
  gx_draw_particles         17.107 ms 10.264 Mcyc 27.4%
  gx_particle_step          27.105 ms 16.263 Mcyc 43.4%
  gx_timeline_step           108.7 us   65.2 kcyc  0.2%
  canvas_clear                93.5 us   56.1 kcyc  0.1%
  canvas_buffer_wait        17.798 ms 10.679 Mcyc 28.5%
```

Wall min/avg/max = 61.802/62.418/63.245 ms.
This window sets the highest dominant-scope average; it does not replace the higher individual render used in the headline.

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1150.1/f  0.493/2.116/22.970 us  3.90%
isr_pack          143.8/f  6.145/7.788/10.321 us  1.79%
isr_dma_submit    143.8/f  0.770/0.935/5.160 us  0.21%
```

- Packing costs 1.79% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.90%, leaving
  60.062 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `gx_particle_step` — 43.5% of the peak window, 27.094 ms/frame.
2. `gx_draw_particles` — 26.8% of the peak window, 16.716 ms/frame.
3. `pov_preserve_half` — 0.2% of the peak window, 0.146 ms/frame.
4. `gx_timeline_step` — 0.2% of the peak window, 0.117 ms/frame.

The preceding canonical capture recorded 44.932 ms peak.
This capture records +0.564 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 45.496, spilled 🟢 0/927 (0.00%).

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
effect=Galaxies  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (Galaxies 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'gx_particle_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (22.7 ppm, frames 865-880)

INVALID: Galaxies
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Galaxies` and
`HS_PROFILE_WINDOW=16`; `just profile Galaxies` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
