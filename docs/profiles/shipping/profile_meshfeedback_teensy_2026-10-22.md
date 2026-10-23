# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-10-22 19:34 America/Los_Angeles on COM4.
Replaces `profile_meshfeedback_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 181 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MeshFeedback profile 181 16` |

Image size (`profile`, this effect only): `FLASH: code:165960, data:196620, headers:9132   free for files:1659904` / `RAM1: variables:313120, code:42584, padding:22952   free for local variables:145632` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/meshfeedback_ship_aa61e787859b_16639e75a3e2`. Profile ELF SHA-256: `16639e75a3e23696898e510474144ddb533688001f66e99e69881d4634f9f357`.

## Frame cadence

**Pass aggregate**: peak frame render is **53.043 ms** (frame 2651), mean live-frame render is **43.817 ms**, and **0/2863** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`mf_feedback_flush` averages 36.465 ms/f across complete post-setup windows; its worst window is 40.306 ms/f (frames 513–528).

Setup frame 1 rendered 5.767 ms. Published live-frame figures use complete counter windows through frame 2864; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:18) recorded peak 🟢 48.53 (12) and spilled 🟢 0/6687 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 2641–2656)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.617 | 1.00 | 100.00% |
| ↳ `mf_mesh_draw` | 9.544 | 1.00 | 15.24% |
| ↳ `mf_feedback_flush` | 32.510 | 1.00 | 51.92% |
| ↳ ↳ `feedback_composite` | 31.899 | 1.00 | 50.94% |
| ↳ ↳ `feedback_populate` | 0.500 | 1.00 | 0.80% |
| ↳ ↳ `feedback_litscan` | 0.005 | 1.00 | 0.01% |
| ↳ `mf_timeline_step` | 0.091 | 1.00 | 0.15% |
| ↳ `mf_apply_params` | 0.046 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.394 | 1.00 | 0.63% |
| ↳ `canvas_buffer_wait` | 19.991 | 1.00 | 31.93% |

Wall min/avg/max = 53.717/62.626/74.200 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `mf_feedback_flush` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 53.043 | 0/214 | 32.943 | 13/13 |
| 6 | — | 50.587 | 0/241 | 39.163 | 14/15 |
| 3 | — | 50.363 | 0/241 | 40.306 | 14/15 |
| 8 | — | 50.066 | 0/241 | 38.751 | 14/15 |
| 11 | — | 49.897 | 0/241 | 32.697 | 14/15 |
| 9 | — | 49.458 | 0/241 | 40.260 | 14/15 |
| 5 | — | 49.350 | 0/241 | 40.197 | 14/15 |
| 10 | — | 48.723 | 0/241 | 33.651 | 14/15 |
| 2 | — | 45.195 | 0/241 | 39.201 | 15/16 |
| 7 | — | 44.844 | 0/241 | 36.376 | 14/15 |
| 1 | — | 44.274 | 0/239 | 35.080 | 14/14 |
| 4 | — | 43.790 | 0/241 | 39.576 | 14/15 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1003318 us
isr_wake       n=18494 cyc min/avg/max=303/1227/10222 ns min/avg/max=505/2045/17036 total=37823 us cpu=3.76%
isr_pack       n=2311 cyc min/avg/max=4938/5827/7195 ns min/avg/max=8230/9711/11991 total=22447 us cpu=2.23%
isr_dma_submit n=2311 cyc min/avg/max=464/560/617 ns min/avg/max=773/933/1028 total=2158 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `mf_feedback_flush` — 32.510 ms/f, 51.92% of the peak window.
2. `mf_mesh_draw` — 9.544 ms/f, 15.24% of the peak window.
3. `canvas_clear` — 0.394 ms/f, 0.63% of the peak window.
4. `mf_timeline_step` — 0.091 ms/f, 0.15% of the peak window.
5. `mf_apply_params` — 0.046 ms/f, 0.07% of the peak window.

README cells: peak 🟢 53.043 (12), spilled 🟢 0/2863 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=MeshFeedback  windows=179
  [PASS] captured >=3 windows (179)
  [PASS] every window names one effect and resolution (MeshFeedback 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 179 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 179 windows are short)
  [PASS] scope 'mf_feedback_flush' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (128.3 ppm, frames 513-528)

INVALID: MeshFeedback
```

Validation limitations: cycle returns to its first index (peak idx 12, 11 distinct); all 12 presets visited (11); root cyc/600 vs wall sum within 5 ppm (128.3 ppm, frames 513-528). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
