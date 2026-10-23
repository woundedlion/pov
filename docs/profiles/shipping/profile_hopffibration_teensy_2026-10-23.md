# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HopfFibration`).
Raw capture: `C:/work/temp/Holosphere-p1polar/build/prof/frozen-hopf.log`, captured 2026-10-23 08:44
America/Los_Angeles on COM3. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, measured source `70fe9492b7e9c96c33791f0d156a8030884da030` |
| Method | 60 s capture, 16-frame windows. Trail nodes remain fixed in source space within each frame; world orientation supplies column-time motion. The measured short-edge path uses known-rigid transport. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-p1polar HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HopfFibration profile 60 16` |

Image size (`profile`, this effect only): `FLASH: code:124680, data:148236, headers:8684   free for files:1750016` / `RAM1: variables:312992, code:28296, padding:4472   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `a3b0a3fe584a5fc7f10d2aad0dbec49438f815ec5e631467ae95b23f2293c5bd`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 193–208
have root cycles/600 = 1015623.480 us
against 1015714 us measured wall sum (**89.1 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **56.986 ms**
(frame 824), mean live render **28.302 ms**,
**0/927 spilled frames**. Setup frame 1 is excluded.
`hf_render_trails` averages 28.148 ms/frame across complete
post-setup windows; its worst window averages
49.251 ms/frame
(frames 193–208).

Setup frame 1 renders 0.702 ms, reported separately.
Complete counter windows end at frame 928;
8 trailing raw frame rows without a complete
counter window are excluded from the published figures. No preset markers are emitted; the capture describes the observed show segment.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: Trail nodes remain fixed in source space within each frame; world orientation supplies column-time motion. The measured short-edge path uses known-rigid transport.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 817–832)

```text
frame                       62.338 ms 37.403 Mcyc 100.0%
  pov_preserve_half          140.8 us   84.5 kcyc  0.2%
  hf_render_trails          47.101 ms 28.260 Mcyc 75.6%
  hf_project_record          171.5 us  102.9 kcyc  0.3%
  hf_advance_tumble           10.8 us    6.5 kcyc  0.0%
  hf_timeline_step            57.9 us   34.8 kcyc  0.1%
  canvas_clear                95.8 us   57.5 kcyc  0.2%
  canvas_buffer_wait        14.736 ms  8.842 Mcyc 23.6%
```

Wall min/avg/max = 38.073/62.343/87.479 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Highest dominant-scope window (frames 193–208)

```text
frame                       63.476 ms 38.086 Mcyc 100.0%
  pov_preserve_half          142.1 us   85.3 kcyc  0.2%
  hf_render_trails          49.251 ms 29.551 Mcyc 77.6%
  hf_project_record          171.7 us  103.0 kcyc  0.3%
  hf_advance_tumble           10.8 us    6.5 kcyc  0.0%
  hf_timeline_step            63.0 us   37.8 kcyc  0.1%
  canvas_clear                95.4 us   57.2 kcyc  0.2%
  canvas_buffer_wait        13.718 ms  8.231 Mcyc 21.6%
```

Wall min/avg/max = 51.044/63.482/76.583 ms.
This window sets the highest dominant-scope average; it does not replace the higher individual render used in the headline.

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1150.5/f  0.598/3.733/28.956 us  6.88%
isr_pack          143.8/f  6.375/8.375/11.753 us  1.92%
isr_dma_submit    143.8/f  0.651/0.938/1.118 us  0.21%
```

- Packing costs 1.92% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 6.88%, leaving
  58.200 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `hf_render_trails` — 75.6% of the peak window, 47.101 ms/frame.
2. `hf_project_record` — 0.3% of the peak window, 0.171 ms/frame.
3. `pov_preserve_half` — 0.2% of the peak window, 0.141 ms/frame.
4. `canvas_clear` — 0.2% of the peak window, 0.096 ms/frame.

The preceding canonical capture recorded 57.358 ms peak.
This capture records -0.372 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 56.986, spilled 🟢 0/927 (0.00%).

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
- Measured source is `70fe9492b7e9c96c33791f0d156a8030884da030`. Provenance identifies the
  actual firmware; later parity/test-only commits are not attributed to it.
- No preset markers are emitted; the capture describes the observed show segment. Strict validation has limitations shown below.

Parser validation:

```text
effect=HopfFibration  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (HopfFibration 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'hf_render_trails' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (89.1 ppm, frames 193-208)

INVALID: HopfFibration
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HopfFibration` and
`HS_PROFILE_WINDOW=16`; `just profile HopfFibration` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
