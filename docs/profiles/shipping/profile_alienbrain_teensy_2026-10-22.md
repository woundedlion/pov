# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/alienbrain_ship.log`, captured 2026-10-22 19:48 America/Los_Angeles on COM4.
Replaces `profile_alienbrain_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 240 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienBrain profile 240 16` |

Image size (`profile`, this effect only): `FLASH: code:88240, data:155080, headers:8584   free for files:1779712` / `RAM1: variables:313024, code:17080, padding:15688   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/alienbrain_ship_aa61e787859b_6751c3b271de`. Profile ELF SHA-256: `6751c3b271def64e3d1aa120fd436d2775d7706c2590ae7eee8f231ff62a0693`.

## Frame cadence

**Pass aggregate**: peak frame render is **26.382 ms** (frame 2096), mean live-frame render is **25.285 ms**, and **0/3807** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 21.421 ms/f across complete post-setup windows; its worst window is 22.097 ms/f (frames 2913–2928).

Setup frame 1 rendered 46.851 ms. Published live-frame figures use complete counter windows through frame 3808; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:35) recorded peak 🟢 25.11 (4) and spilled 🟢 0/4767 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 2081–2096)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.374 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 22.057 | 1.00 | 35.36% |
| ↳ ↳ `fx_prepare_column` | 1.130 | 144.00 | 1.81% |
| ↳ `fx_prepare_frame` | 1.105 | 1.00 | 1.77% |
| ↳ `fx_advance` | 2.261 | 1.00 | 3.63% |
| ↳ `fx_timeline_step` | 0.208 | 1.00 | 0.33% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 36.151 | 1.00 | 57.96% |

Wall min/avg/max = 61.937/62.387/62.745 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 3 | — | 26.382 | 0/1079 | 22.096 | 67/68 |
| 2 | — | 26.327 | 0/1079 | 22.084 | 66/67 |
| 4 | — | 26.136 | 0/1051 | 22.097 | 65/65 |
| 1 | — | 25.671 | 0/598 | 21.414 | 36/37 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999986 us
isr_wake       n=18433 cyc min/avg/max=322/1165/16457 ns min/avg/max=536/1941/27428 total=35819 us cpu=3.58%
isr_pack       n=2304 cyc min/avg/max=4550/5305/6892 ns min/avg/max=7583/8841/11486 total=20373 us cpu=2.03%
isr_dma_submit n=2304 cyc min/avg/max=465/561/5738 ns min/avg/max=775/935/9563 total=2154 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 22.057 ms/f, 35.36% of the peak window.
2. `fx_advance` — 2.261 ms/f, 3.63% of the peak window.
3. `fx_prepare_frame` — 1.105 ms/f, 1.77% of the peak window.
4. `fx_timeline_step` — 0.208 ms/f, 0.33% of the peak window.
5. `pov_preserve_half` — 0.143 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.097 ms/f, 0.16% of the peak window.

README cells: peak 🟢 26.382 (4), spilled 🟢 0/3807 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=AlienBrain  windows=238
  [PASS] captured >=3 windows (238)
  [PASS] every window names one effect and resolution (AlienBrain 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 4, 3 distinct)
  [FAIL] all 4 presets visited (3)
  [PASS] per-frame render is render, not wall (0 of 238 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 238 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (222.2 ppm, frames 1-16)

INVALID: AlienBrain
```

Validation limitations: cycle returns to its first index (peak idx 4, 3 distinct); all 4 presets visited (3); root cyc/600 vs wall sum within 5 ppm (222.2 ppm, frames 1-16). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
