# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `C:/work/temp/Holosphere-p1polar/build/prof/frozen-shapeshifter.log`, captured 2026-10-23 08:42
America/Los_Angeles on COM3. Replaces the October 22 canonical
shipping report; this capture’s coverage is stated independently below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` with selective `HS_O3` Plot/raster and effect hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, measured source `a6222f23b579ecdddb320cc991115d4acef8aaf3` |
| Method | 270 s capture, 16-frame windows. P1 retains its still route; P2 through P9 use column-time contours. The log reaches P9 without completing its hold or returning to P1. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-p1polar HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ShapeShifter profile 270 16 '-D HS_PROFILE_EPOCH_REVS=2560'` |

Image size (`profile`, this effect only): `FLASH: code:192728, data:155844, headers:8804   free for files:1674240` / `RAM1: variables:312992, code:53832, padding:11704   free for local variables:145760` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:993104, data:851736, headers:8600   free for files:178176` / `RAM1: variables:312736, code:171432, padding:25176   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
Frozen combined-production `.text.itcm` is 171424 B; RAM1 code is 171432 B,
variables 312736 B, with six code banks and 25176 B padding. Relative to the
landed IslamicStars image, RAM1 code rises 912 B (0.891 KiB), variables are
unchanged, FLASH data rises 68 B (0.066 KiB), and FLASH code rises 66864 B
(65.297 KiB). Relative to the earlier Mermaid/Hex master, RAM1 code rises
2096 B (2.047 KiB). No new discretionary pins are
included; the new uniform-color and polar-walk bodies remain in FLASH.
Image size and layout describe the captured source, not an assumed newer tip.

Profile ELF SHA-256: `dea7fecb4dcf52564dc38aaa242cac23bfda58a1e8ac10ce26f94754af70d5b9`.
The adjacent `.provenance`, envdump and build logs record compiler/source
attestation, exact flags and immutable image artifacts. Later commits do not
replace this captured source identity.

Engine revision `10df63d5432710f2ab22f1f726b7127ad0d6242e` additionally
corrects the PointRoute equator-cell speed upper bound. This capture predates
that correction. Timing and size figures remain captured-source measurements;
no device timing or performance claim is available for the correction.

Exactness cross-check: parser-selected frames 3553–3568
have root cycles/600 = 995967.372 us
against 996021 us measured wall sum (**53.8 ppm**).
The strict 5 ppm gate fails.
The headline uses logged per-frame render telemetry; all validator diagnostics
are retained in Caveats.

## Frame cadence

**Pass aggregate:** peak frame render **56.941 ms**
(frame 410), mean live render **35.022 ms**,
**0/4287 spilled frames**. Setup frame 1 is excluded.
`ss_draw_all` averages 34.641 ms/frame across complete
post-setup windows; its worst window averages
53.496 ms/frame
(frames 3553–3568).

Setup frame 1 renders 47.868 ms, reported separately.
Complete counter windows end at frame 4288;
8 trailing raw frame rows without a complete
counter window are excluded from the published figures. The log explicitly adopts P2 through P9; startup P1 ownership is inferred. P2 through P8 each contribute 479 complete-window live frames; P9 contributes a 464-frame partial tail. The capture does not return to P1.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. Every retained live render fits the 16 fps budget.
Wall time includes synchronization and scheduling jitter.
`canvas_buffer_wait` is idle waiting for the next display flip, excluded from
render peak and spill classification.

## Phase-by-phase readout

Phase schedule: P1 retains its still route; P2 through P9 use column-time contours. The log reaches P9 without completing its hold or returning to P1.
Trees preserve logged nesting and inclusive per-frame averages. Mixed-parent
and duplicate-name rows are not exclusive costs of the displayed branch.

### Peak-render regime (frames 401–416)

```text
frame                       62.793 ms 37.676 Mcyc 100.0%
  pov_preserve_half          138.1 us   82.9 kcyc  0.2%
  ss_draw_all               44.242 ms 26.545 Mcyc 70.5%
    ss_plot_dispatch        43.682 ms 26.209 Mcyc 69.6%
  ss_timeline_step            84.2 us   50.5 kcyc  0.1%
  ss_buffer_wait            18.318 ms 10.991 Mcyc 29.2%
    canvas_clear              85.0 us   51.0 kcyc  0.1%
    canvas_buffer_wait      18.226 ms 10.936 Mcyc 29.0%
```

Wall min/avg/max = 39.082/62.798/85.741 ms.
The tree shows inclusive averages over the window containing the highest individual live render. The peak frame and this window’s average are distinct measurements.

### Highest dominant-scope window (frames 3553–3568)

```text
frame                       62.248 ms 37.349 Mcyc 100.0%
  pov_preserve_half          138.9 us   83.4 kcyc  0.2%
  ss_draw_all               53.496 ms 32.098 Mcyc 85.9%
  ss_timeline_step            96.9 us   58.1 kcyc  0.2%
  ss_buffer_wait             8.498 ms  5.099 Mcyc 13.7%
    canvas_clear              94.3 us   56.6 kcyc  0.2%
    canvas_buffer_wait       8.390 ms  5.034 Mcyc 13.5%
```

Wall min/avg/max = 56.850/62.251/67.526 ms.
This window sets the highest dominant-scope average; it does not replace the higher individual render used in the headline.

### Preset-adoption boundary (frames 3345–3360)

```text
frame                       62.626 ms 37.576 Mcyc 100.0%
  pov_preserve_half          139.4 us   83.7 kcyc  0.2%
  ss_draw_all               45.858 ms 27.515 Mcyc 73.2%
  ss_timeline_step           115.6 us   69.4 kcyc  0.2%
  ss_buffer_wait            16.495 ms  9.897 Mcyc 26.3%
    canvas_clear              92.9 us   55.8 kcyc  0.1%
    canvas_buffer_wait      16.387 ms  9.832 Mcyc 26.2%
```

Wall min/avg/max = 53.481/62.628/71.713 ms.
This window mixes owners across a preset adoption and is excluded from the single-owner window figures. Its live frames remain in aggregate and ownership peak/spill counts.

### Per-preset table

The log explicitly adopts P2 through P9; startup P1 ownership is inferred. P2 through P8 each contribute 479 complete-window live frames; P9 contributes a 464-frame partial tail. The capture does not return to P1.
Ownership includes the transition following an adoption. The window figures
use the highest-scope post-setup window with one owner; they are not certified
full-opacity clean holds. Where adoption timing gives no separately identified
hold, that distinction remains explicit. Window counts exclude advance-straddling
windows; per-entry peak/spill figures include every owned live frame.
No mesh V/E/F/I metadata is logged for these entries.

| Entry | Blended px/f | `ss_draw_all` ms/f | Render ms/f | Wall fps | Windows | Owned peak ms | Spilled/frames |
|---|--:|--:|--:|--:|--:|--:|--:|
| P8 | n/a | 53.496 | 53.750 | 16.06 | 29 | 56.529 | 0/479 |
| P7 | n/a | 50.941 | 51.197 | 15.94 | 29 | 53.531 | 0/479 |
| P1 (startup, inferred) | n/a | 47.120 | 47.366 | 16.08 | 28 | 56.941 | 0/470 |
| P6 | n/a | 43.406 | 43.663 | 16.01 | 29 | 45.033 | 0/479 |
| P4 | n/a | 40.578 | 40.874 | 16.02 | 29 | 42.373 | 0/479 |
| P9 (partial tail) | n/a | 34.302 | 34.607 | 16.19 | 29 | 36.373 | 0/464 |
| P5 | n/a | 31.890 | 32.182 | 16.10 | 29 | 38.066 | 0/479 |
| P2 | n/a | 25.299 | 25.552 | 15.96 | 29 | 29.771 | 0/479 |
| P3 | n/a | 25.255 | 25.545 | 16.06 | 29 | 29.958 | 0/479 |

### Per-pixel figures

Blended-pixel telemetry is absent. Blends per frame, cycles per blend and scan cycles per blended pixel are unavailable; 10368 pixels is the geometric quadrant reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-window rates per rendered frame and per-call min/average/maximum:

```text
isr_wake         1158.8/f  0.546/1.918/13.590 us  3.53%
isr_pack          144.8/f  7.048/8.795/11.355 us  2.02%
isr_dma_submit    144.8/f  0.671/0.940/1.201 us  0.21%
```

- Packing costs 2.02% and submission 0.21%
  of the separately timestamped ISR interval. LED wire transfer is asynchronous;
  wire duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.53%, leaving
  60.294 ms per display window after that interrupt load.
  Pack and submit are nested inside wake; never add their shares to it.
  Render scopes already include interruptions, so do not subtract them twice.
- No observed live render requires a speedup for the display budget.
  DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `ss_draw_all` — 70.5% of the peak window, 44.242 ms/frame.
2. `ss_plot_dispatch` — 69.6% of the peak window, 43.682 ms/frame.
3. `ss_buffer_wait` — 29.2% of the peak window, 18.318 ms/frame.
4. `pov_preserve_half` — 0.2% of the peak window, 0.138 ms/frame.

The preceding canonical capture recorded 55.606 ms peak.
This capture records +1.335 ms against that
observation; sources and sampled motion differ, so this is not a matched
baseline experiment. Native/WASM rates do not establish on-device cadence.

README cells: peak 🟢 56.941 (9 observed; no wrap), spilled 🟢 0/4287 (0.00%).

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
- The log explicitly adopts P2 through P9; startup P1 ownership is inferred. P2 through P8 each contribute 479 complete-window live frames; P9 contributes a 464-frame partial tail. The capture does not return to P1. Strict validation has limitations shown below.

Parser validation:

```text
effect=ShapeShifter  windows=268
  [PASS] captured >=3 windows (268)
  [PASS] every window names one effect and resolution (ShapeShifter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 268 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 268 windows are short)
  [PASS] scope 'ss_draw_all' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (53.8 ppm, frames 3553-3568)

INVALID: ShapeShifter
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; `just profile ShapeShifter` invokes the standard
wrapper. Reproduce above names the actual capture duration and flags under
the shared device lock.
