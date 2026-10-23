# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mindsplatter_ship.log`, captured 2026-10-22 19:20 America/Los_Angeles on COM4.
Replaces `profile_mindsplatter_teensy_2026-10-16.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 80 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MindSplatter profile 80 16` |

Image size (`profile`, this effect only): `FLASH: code:129432, data:550396, headers:8300   free for files:1343488` / `RAM1: variables:313152, code:49208, padding:16328   free for local variables:145600` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/mindsplatter_ship_aa61e787859b_7939fc906843`. Profile ELF SHA-256: `7939fc906843616b92e9fc985eacdcfc3bc12020c6e4aabed50da30ef1da4713`.

## Frame cadence

**Pass aggregate**: peak frame render is **58.709 ms** (frame 1003), mean live-frame render is **27.614 ms**, and **0/1247** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`msp_draw_particles` averages 24.982 ms/f across complete post-setup windows; its worst window is 50.813 ms/f (frames 993–1008).

Setup frame 1 rendered 0.390 ms. Published live-frame figures use complete counter windows through frame 1248; a trailing partial window is excluded.

The previous shipping report (2026-10-16 08:55) recorded peak 🟢 55.71 (8) and spilled 🟢 0/1727 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 993–1008)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.269 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `msp_draw_particles` | 50.813 | 1.00 | 81.60% |
| ↳ `msp_particle_step` | 2.985 | 1.00 | 4.79% |
| ↳ `msp_timeline_step` | 0.137 | 1.00 | 0.22% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 8.059 | 1.00 | 12.94% |

Wall min/avg/max = 58.433/62.274/65.953 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `msp_draw_particles` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 7 | — | 58.709 | 0/159 | 50.813 | 9/10 |
| 6 | — | 53.320 | 0/159 | 39.297 | 9/10 |
| 2 | — | 48.555 | 0/159 | 39.550 | 9/10 |
| 4 | — | 46.248 | 0/159 | 33.227 | 9/10 |
| 5 | — | 36.569 | 0/159 | 25.864 | 9/10 |
| 3 | — | 34.057 | 0/159 | 29.328 | 9/10 |
| 8 | — | 32.097 | 0/135 | 28.082 | 8/8 |
| 1 | — | 25.829 | 0/158 | 19.518 | 8/9 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999073 us
isr_wake       n=18416 cyc min/avg/max=357/1192/16537 ns min/avg/max=595/1986/27561 total=36598 us cpu=3.66%
isr_pack       n=2302 cyc min/avg/max=3818/5129/6811 ns min/avg/max=6363/8548/11351 total=19679 us cpu=1.96%
isr_dma_submit n=2302 cyc min/avg/max=389/560/2754 ns min/avg/max=648/933/4590 total=2149 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `msp_draw_particles` — 50.813 ms/f, 81.60% of the peak window.
2. `msp_particle_step` — 2.985 ms/f, 4.79% of the peak window.
3. `msp_timeline_step` — 0.137 ms/f, 0.22% of the peak window.
4. `pov_preserve_half` — 0.136 ms/f, 0.22% of the peak window.
5. `canvas_clear` — 0.093 ms/f, 0.15% of the peak window.

README cells: peak 🟢 58.709 (8), spilled 🟢 0/1247 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=MindSplatter  windows=78
  [PASS] captured >=3 windows (78)
  [PASS] every window names one effect and resolution (MindSplatter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 8, 7 distinct)
  [FAIL] all 8 presets visited (7)
  [PASS] per-frame render is render, not wall (0 of 78 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 78 windows are short)
  [PASS] scope 'msp_draw_particles' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (84.4 ppm, frames 993-1008)

INVALID: MindSplatter
```

Validation limitations: cycle returns to its first index (peak idx 8, 7 distinct); all 8 presets visited (7); root cyc/600 vs wall sum within 5 ppm (84.4 ppm, frames 993-1008). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
