# Comets on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/comets_ship.log`, captured 2026-10-22 19:30 America/Los_Angeles on COM4.
Replaces `profile_comets_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Comets 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 120 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Comets profile 120 16` |

Image size (`profile`, this effect only): `FLASH: code:90096, data:153016, headers:8792   free for files:1779712` / `RAM1: variables:312960, code:26728, padding:6040   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/comets_ship_aa61e787859b_1f67af43f716`. Profile ELF SHA-256: `1f67af43f716d852fea31a6d4cb35c4e4ba0defcc2412cac8fcd5ec8359526c4`.

## Frame cadence

**Pass aggregate**: peak frame render is **30.550 ms** (frame 647), mean live-frame render is **15.420 ms**, and **0/1887** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`cm_draw_trail` averages 14.579 ms/f across complete post-setup windows; its worst window is 20.657 ms/f (frames 513–528).

Setup frame 1 rendered 0.888 ms. Published live-frame figures use complete counter windows through frame 1888; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:07) recorded peak 🟢 31.92 (12) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 641–656)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.080 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `cm_draw_trail` | 20.092 | 1.00 | 32.36% |
| ↳ `cm_wipe_rebake` | 1.979 | 1.00 | 3.19% |
| ↳ `cm_timeline_step` | 0.169 | 1.00 | 0.27% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 39.597 | 1.00 | 63.79% |

Wall min/avg/max = 45.579/62.083/77.494 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `cm_draw_trail` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 5 | — | 30.550 | 0/160 | 20.092 | 9/10 |
| 11 | — | 27.557 | 0/160 | 18.121 | 9/10 |
| 10 | — | 27.221 | 0/160 | 20.301 | 9/10 |
| 4 | — | 26.451 | 0/160 | 20.657 | 9/10 |
| 9 | — | 25.965 | 0/160 | 19.080 | 9/10 |
| 2 | — | 25.278 | 0/160 | 19.790 | 9/10 |
| 12 | — | 23.437 | 0/129 | 20.224 | 8/8 |
| 8 | — | 22.577 | 0/160 | 19.585 | 9/10 |
| 6 | — | 21.102 | 0/160 | 15.101 | 9/10 |
| 3 | — | 19.144 | 0/160 | 15.043 | 9/10 |
| 7 | — | 15.904 | 0/160 | 12.626 | 9/10 |
| 1 | — | 14.944 | 0/158 | 12.310 | 8/9 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=994306 us
isr_wake       n=18328 cyc min/avg/max=292/1009/10262 ns min/avg/max=486/1681/17103 total=30841 us cpu=3.10%
isr_pack       n=2291 cyc min/avg/max=3594/4253/6354 ns min/avg/max=5990/7088/10590 total=16242 us cpu=1.63%
isr_dma_submit n=2291 cyc min/avg/max=393/563/5703 ns min/avg/max=655/938/9505 total=2151 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `cm_draw_trail` — 20.092 ms/f, 32.36% of the peak window.
2. `cm_wipe_rebake` — 1.979 ms/f, 3.19% of the peak window.
3. `cm_timeline_step` — 0.169 ms/f, 0.27% of the peak window.
4. `pov_preserve_half` — 0.136 ms/f, 0.22% of the peak window.
5. `canvas_clear` — 0.089 ms/f, 0.14% of the peak window.

README cells: peak 🟢 30.550 (12), spilled 🟢 0/1887 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=Comets  windows=118
  [PASS] captured >=3 windows (118)
  [PASS] every window names one effect and resolution (Comets 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 118 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 118 windows are short)
  [PASS] scope 'cm_draw_trail' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (69.0 ppm, frames 513-528)

INVALID: Comets
```

Validation limitations: cycle returns to its first index (peak idx 12, 11 distinct); all 12 presets visited (11); root cyc/600 vs wall sum within 5 ppm (69.0 ppm, frames 513-528). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
