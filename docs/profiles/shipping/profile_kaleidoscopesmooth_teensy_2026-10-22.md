# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopesmooth_ship.log`, captured 2026-10-22 20:14 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopesmooth_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 240 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeSmooth profile 240 16` |

Image size (`profile`, this effect only): `FLASH: code:88496, data:155244, headers:9188   free for files:1778688` / `RAM1: variables:313024, code:17080, padding:15688   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopesmooth_ship_aa61e787859b_a0c66b4913c7`. Profile ELF SHA-256: `a0c66b4913c7320c7baa39ed659c4dfca15649ac8b96f177b9efd7c7f98542cc`.

## Frame cadence

**Pass aggregate**: peak frame render is **33.333 ms** (frame 1058), mean live-frame render is **28.034 ms**, and **0/3807** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 24.431 ms/f across complete post-setup windows; its worst window is 29.493 ms/f (frames 1057–1072).

Setup frame 1 rendered 51.144 ms. Published live-frame figures use complete counter windows through frame 3808; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:55) recorded peak 🟢 32.02 (4) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 1057–1072)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.346 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 29.493 | 1.00 | 47.31% |
| ↳ ↳ `fx_prepare_column` | 0.874 | 144.00 | 1.40% |
| ↳ `fx_prepare_frame` | 0.948 | 1.00 | 1.52% |
| ↳ `fx_advance` | 1.714 | 1.00 | 2.75% |
| ↳ `fx_timeline_step` | 0.209 | 1.00 | 0.34% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 29.344 | 1.00 | 47.07% |

Wall min/avg/max = 61.462/62.368/63.061 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 33.333 | 0/1079 | 29.493 | 66/67 |
| 4 | — | 31.369 | 0/1051 | 26.803 | 65/65 |
| 3 | — | 31.086 | 0/1079 | 26.268 | 67/68 |
| 1 | — | 31.037 | 0/598 | 27.251 | 36/37 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999762 us
isr_wake       n=18429 cyc min/avg/max=335/1204/15589 ns min/avg/max=558/2006/25981 total=36996 us cpu=3.70%
isr_pack       n=2304 cyc min/avg/max=5191/5594/7094 ns min/avg/max=8651/9323/11823 total=21481 us cpu=2.14%
isr_dma_submit n=2304 cyc min/avg/max=466/563/1313 ns min/avg/max=776/938/2188 total=2165 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 29.493 ms/f, 47.31% of the peak window.
2. `fx_advance` — 1.714 ms/f, 2.75% of the peak window.
3. `fx_prepare_frame` — 0.948 ms/f, 1.52% of the peak window.
4. `fx_timeline_step` — 0.209 ms/f, 0.34% of the peak window.
5. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.093 ms/f, 0.15% of the peak window.

README cells: peak 🟢 33.333 (4), spilled 🟢 0/3807 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopeSmooth  windows=238
  [PASS] captured >=3 windows (238)
  [PASS] every window names one effect and resolution (KaleidoscopeSmooth 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 4, 3 distinct)
  [FAIL] all 4 presets visited (3)
  [PASS] per-frame render is render, not wall (0 of 238 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 238 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (355.8 ppm, frames 1057-1072)

INVALID: KaleidoscopeSmooth
```

Validation limitations: cycle returns to its first index (peak idx 4, 3 distinct); all 4 presets visited (3); root cyc/600 vs wall sum within 5 ppm (355.8 ppm, frames 1057-1072). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
