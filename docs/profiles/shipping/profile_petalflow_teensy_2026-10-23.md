# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PetalFlow`).
Raw capture: `C:/work/temp/Holosphere-p1polar/build/prof/frozen-petalflow.log`, captured 2026-10-23 08:46
America/Los_Angeles on COM3. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, measured source `70fe9492b7e9c96c33791f0d156a8030884da030` |
| Method | 30 s capture, 16-frame windows. Shipping PetalFlow configuration, live driver, no fixed-preset override. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-p1polar HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh PetalFlow profile 30 16` |

Image size (`profile`, this effect only): `FLASH: code:109896, data:148140, headers:8204   free for files:1765376` / `RAM1: variables:312992, code:28728, padding:4040   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `66b03a7c9adf53ead29f25cfc26380b57e88990d55ee793003092e4a34de21d5`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 305–320
have root cycles/600 = 1000890.855 us
against 1000938 us measured wall sum (**47.1 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **30.159 ms**
(frame 307), mean live render **24.663 ms**,
**0/447 spilled frames**. Setup frame 1 is excluded.
`pf_draw_rings` averages 24.430 ms/frame across complete
post-setup windows; its worst window averages
28.305 ms/frame
(frames 305–320).

Setup frame 1 renders 21.570 ms, reported separately.
Complete counter windows end at frame 448;
9 trailing raw frame rows without a complete
counter window are excluded from the published figures. No preset markers are emitted; the capture describes the observed show segment.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: Shipping PetalFlow configuration, live driver, no fixed-preset override.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 305–320)

```text
frame                       62.556 ms 37.533 Mcyc 100.0%
  pov_preserve_half          134.1 us   80.4 kcyc  0.2%
  pf_draw_rings             28.305 ms 16.983 Mcyc 45.2%
    pf_ring_scan            27.239 ms 16.344 Mcyc 43.5%
    pf_ring_build            888.0 us  532.8 kcyc  1.4%
  pf_timeline_step            59.6 us   35.7 kcyc  0.1%
  canvas_clear                93.7 us   56.2 kcyc  0.1%
  canvas_buffer_wait        33.919 ms 20.351 Mcyc 54.2%
```

Wall min/avg/max = 60.691/62.558/64.456 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1154.2/f  0.555/1.688/19.201 us  3.11%
isr_pack          144.2/f  6.081/6.970/9.566 us  1.60%
isr_dma_submit    144.2/f  0.621/0.936/3.558 us  0.21%
```

- Packing costs 1.60% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.11%, leaving
  60.556 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `pf_draw_rings` — 45.2% of the peak window, 28.305 ms/frame.
2. `pf_ring_scan` — 43.5% of the peak window, 27.239 ms/frame.
3. `pf_ring_build` — 1.4% of the peak window, 0.888 ms/frame.
4. `pov_preserve_half` — 0.2% of the peak window, 0.134 ms/frame.

The preceding canonical capture recorded 30.379 ms peak.
This capture records -0.220 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 30.159, spilled 🟢 0/447 (0.00%).

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
effect=PetalFlow  windows=28
  [PASS] captured >=3 windows (28)
  [PASS] every window names one effect and resolution (PetalFlow 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 28 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 28 windows are short)
  [PASS] scope 'pf_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (47.1 ppm, frames 305-320)

INVALID: PetalFlow
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=PetalFlow` and
`HS_PROFILE_WINDOW=16`; `just profile PetalFlow` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
