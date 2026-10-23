# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/petalflow_ship.log`, captured 2026-10-22 19:09 America/Los_Angeles on COM4.
Replaces `profile_petalflow_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh PetalFlow profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:89472, data:148132, headers:9180   free for files:1784832` / `RAM1: variables:312992, code:29752, padding:3016   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/petalflow_ship_aa61e787859b_08d250477d40`. Profile ELF SHA-256: `08d250477d40898e23c52e4d5752d5cfd8ec44ac2faa3a9bd3e7de89db49076a`.

## Frame cadence

**Pass aggregate**: peak frame render is **30.379 ms** (frame 309), mean live-frame render is **24.809 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`pf_draw_rings` averages 24.727 ms/f across complete post-setup windows; its worst window is 27.957 ms/f (frames 321–352).

Setup frame 1 rendered 21.974 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:12) recorded peak 🟢 13.29 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 289–320)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.521 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `pf_draw_rings` | 27.901 | 1.00 | 44.63% |
| ↳ ↳ `pf_ring_scan` | 26.751 | 23.34 | 42.79% |
| ↳ ↳ `pf_ring_build` | 0.919 | 23.34 | 1.47% |
| ↳ `pf_timeline_step` | 0.032 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.348 | 1.00 | 54.94% |

Wall min/avg/max = 59.034/62.537/65.297 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2003490 us
isr_wake       n=36931 cyc min/avg/max=323/1131/12902 ns min/avg/max=538/1885/21503 total=69619 us cpu=3.47%
isr_pack       n=4616 cyc min/avg/max=3664/4202/6050 ns min/avg/max=6106/7003/10083 total=32332 us cpu=1.61%
isr_dma_submit n=4616 cyc min/avg/max=469/564/610 ns min/avg/max=781/940/1016 total=4341 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `pf_draw_rings` — 27.901 ms/f, 44.63% of the peak window.
2. `pov_preserve_half` — 0.136 ms/f, 0.22% of the peak window.
3. `canvas_clear` — 0.091 ms/f, 0.15% of the peak window.
4. `pf_timeline_step` — 0.032 ms/f, 0.05% of the peak window.

README cells: peak 🟢 30.379, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=PetalFlow  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (PetalFlow 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'pf_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (265.9 ppm, frames 321-352)

INVALID: PetalFlow
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (265.9 ppm, frames 321-352). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
