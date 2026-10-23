# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/dreamballs_ship.log`, captured 2026-10-22 19:43 America/Los_Angeles on COM4.
Replaces `profile_dreamballs_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 200 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh DreamBalls profile 200 16` |

Image size (`profile`, this effect only): `FLASH: code:139736, data:196704, headers:8648   free for files:1686528` / `RAM1: variables:313120, code:32952, padding:32584   free for local variables:145632` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/dreamballs_ship_aa61e787859b_db638fe284bf`. Profile ELF SHA-256: `db638fe284bf3d2305eef70ab003c73ce2791fa6712938879dfcac8dc4546458`.

## Frame cadence

**Pass aggregate**: peak frame render is **54.443 ms** (frame 2773), mean live-frame render is **29.625 ms**, and **0/3167** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`db_timeline_step` averages 29.308 ms/f across complete post-setup windows; its worst window is 51.342 ms/f (frames 2769–2784).

Setup frame 1 rendered 47.919 ms. Published live-frame figures use complete counter windows through frame 3168; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:28) recorded peak 🟢 41.44 (10) and spilled 🟢 0/3647 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 2769–2784)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.350 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.23% |
| ↳ `db_timeline_step` | 51.342 | 1.00 | 82.34% |
| ↳ ↳ `db_draw` | 51.217 | 1.00 | 82.14% |
| ↳ ↳ ↳ `db_draw_scene` | 51.199 | 1.00 | 82.11% |
| ↳ ↳ ↳ ↳ `db_mesh_plot` | 49.811 | 6.00 | 79.89% |
| ↳ ↳ ↳ ↳ `db_displace` | 1.121 | 6.00 | 1.80% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 10.740 | 1.00 | 17.23% |

Wall min/avg/max = 57.946/62.353/65.450 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `db_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 9 | — | 54.443 | 0/320 | 51.342 | 19/20 |
| 1 | — | 47.539 | 0/318 | 44.563 | 18/19 |
| 10 | — | 41.164 | 0/289 | 37.447 | 18/18 |
| 8 | — | 39.908 | 0/320 | 37.273 | 19/20 |
| 4 | — | 31.962 | 0/320 | 28.732 | 19/20 |
| 3 | — | 30.838 | 0/320 | 27.792 | 19/20 |
| 2 | — | 29.369 | 0/320 | 25.796 | 19/20 |
| 7 | — | 27.185 | 0/320 | 25.151 | 19/20 |
| 6 | — | 20.898 | 0/320 | 18.798 | 19/20 |
| 5 | — | 19.732 | 0/320 | 17.653 | 19/20 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999470 us
isr_wake       n=18424 cyc min/avg/max=355/1115/17521 ns min/avg/max=591/1858/29201 total=34246 us cpu=3.42%
isr_pack       n=2303 cyc min/avg/max=3954/4829/6645 ns min/avg/max=6590/8048/11075 total=18538 us cpu=1.85%
isr_dma_submit n=2303 cyc min/avg/max=434/559/616 ns min/avg/max=723/931/1026 total=2148 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `db_timeline_step` — 51.342 ms/f, 82.34% of the peak window.
2. `pov_preserve_half` — 0.144 ms/f, 0.23% of the peak window.
3. `canvas_clear` — 0.089 ms/f, 0.14% of the peak window.

README cells: peak 🟢 54.443 (10), spilled 🟢 0/3167 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=DreamBalls  windows=198
  [PASS] captured >=3 windows (198)
  [PASS] every window names one effect and resolution (DreamBalls 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 10, 9 distinct)
  [FAIL] all 10 presets visited (9)
  [PASS] per-frame render is render, not wall (0 of 198 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 198 windows are short)
  [PASS] scope 'db_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (54.3 ppm, frames 2769-2784)

INVALID: DreamBalls
```

Validation limitations: cycle returns to its first index (peak idx 10, 9 distinct); all 10 presets visited (9); root cyc/600 vs wall sum within 5 ppm (54.3 ppm, frames 2769-2784). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
