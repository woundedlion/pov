# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-10-22 19:18 America/Los_Angeles on COM4.
Replaces `profile_shapeshifter_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 269 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ShapeShifter profile 269 16` |

Image size (`profile`, this effect only): `FLASH: code:108592, data:155712, headers:9104   free for files:1758208` / `RAM1: variables:312992, code:43080, padding:22456   free for local variables:145760` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/shapeshifter_ship_aa61e787859b_6da263e729f8`. Profile ELF SHA-256: `6da263e729f8703f4431dc66e11e1c8ec6bc8deba0c5eceb7497c852eb8d7f9f`.

## Frame cadence

**Pass aggregate**: peak frame render is **55.606 ms** (frame 410), mean live-frame render is **18.065 ms**, and **0/4271** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`ss_draw_all` averages 17.716 ms/f across complete post-setup windows; its worst window is 46.074 ms/f (frames 449–464).

Setup frame 1 rendered 47.697 ms. Published live-frame figures use complete counter windows through frame 4272; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:29) recorded peak 🟢 54.27 (9) and spilled 🟢 0/2447 (0.00%).

ShapeShifter uses its shipping still route.

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 401–416)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.764 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `ss_draw_all` | 43.180 | 1.00 | 68.80% |
| ↳ ↳ `ss_plot_dispatch` | 42.702 | 214.94 | 68.03% |
| ↳ `ss_timeline_step` | 0.054 | 1.00 | 0.09% |
| ↳ `ss_buffer_wait` | 19.379 | 1.00 | 30.88% |
| ↳ ↳ `canvas_clear` | 0.086 | 1.00 | 0.14% |
| ↳ ↳ `canvas_buffer_wait` | 19.285 | 1.00 | 30.73% |

Wall min/avg/max = 39.368/62.766/85.383 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `ss_draw_all` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 55.606 | 0/470 | 46.074 | 28/29 |
| 9 | — | 24.867 | 0/448 | 24.218 | 28/28 |
| 8 | — | 24.723 | 0/479 | 22.976 | 29/29 |
| 4 | — | 24.599 | 0/479 | 23.355 | 29/30 |
| 7 | — | 22.831 | 0/479 | 20.503 | 29/30 |
| 6 | — | 18.250 | 0/479 | 17.089 | 29/30 |
| 2 | — | 14.377 | 0/479 | 11.979 | 29/30 |
| 5 | — | 9.598 | 0/479 | 8.061 | 29/30 |
| 3 | — | 7.435 | 0/479 | 6.060 | 29/30 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1005333 us
isr_wake       n=18532 cyc min/avg/max=309/1149/16046 ns min/avg/max=515/1915/26743 total=35516 us cpu=3.53%
isr_pack       n=2316 cyc min/avg/max=4184/5280/7043 ns min/avg/max=6973/8800/11738 total=20384 us cpu=2.02%
isr_dma_submit n=2316 cyc min/avg/max=493/562/625 ns min/avg/max=821/936/1041 total=2171 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `ss_draw_all` — 43.180 ms/f, 68.80% of the peak window.
2. `pov_preserve_half` — 0.139 ms/f, 0.22% of the peak window.
3. `ss_timeline_step` — 0.054 ms/f, 0.09% of the peak window.

README cells: peak 🟢 55.606 (9), spilled 🟢 0/4271 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=ShapeShifter  windows=267
  [PASS] captured >=3 windows (267)
  [PASS] every window names one effect and resolution (ShapeShifter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 267 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 267 windows are short)
  [PASS] scope 'ss_draw_all' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (33.1 ppm, frames 449-464)

INVALID: ShapeShifter
```

Validation limitations: cycle returns to its first index (peak idx 9, 8 distinct); all 9 presets visited (8); root cyc/600 vs wall sum within 5 ppm (33.1 ppm, frames 449-464). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
