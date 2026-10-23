# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mobiusgrid_ship.log`, captured 2026-10-22 19:08 America/Los_Angeles on COM4.
Replaces `profile_mobiusgrid_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 105 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MobiusGrid profile 105 16` |

Image size (`profile`, this effect only): `FLASH: code:89560, data:155664, headers:8728   free for files:1777664` / `RAM1: variables:313024, code:17048, padding:15720   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/mobiusgrid_ship_aa61e787859b_921dd9a1e018`. Profile ELF SHA-256: `921dd9a1e0183772afaf553ede5e7933471cf86863a251cb7583fd1eb25571f4`.

## Frame cadence

**Pass aggregate**: peak frame render is **30.419 ms** (frame 792), mean live-frame render is **27.704 ms**, and **0/1647** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 24.280 ms/f across complete post-setup windows; its worst window is 26.115 ms/f (frames 849–864).

Setup frame 1 rendered 45.274 ms. Published live-frame figures use complete counter windows through frame 1648; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:08) recorded peak 🟢 21.77 (2) and spilled 🟢 0/2687 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 785–800)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.378 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 26.106 | 1.00 | 41.85% |
| ↳ ↳ `fx_prepare_column` | 2.910 | 144.00 | 4.67% |
| ↳ `fx_prepare_frame` | 0.684 | 1.00 | 1.10% |
| ↳ `fx_advance` | 2.219 | 1.00 | 3.56% |
| ↳ `fx_timeline_step` | 0.220 | 1.00 | 0.35% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 32.732 | 1.00 | 52.47% |

Wall min/avg/max = 60.983/62.392/63.780 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 30.419 | 0/1049 | 26.115 | 65/65 |
| 1 | — | 27.998 | 0/598 | 23.754 | 36/37 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1000245 us
isr_wake       n=18437 cyc min/avg/max=369/1200/14761 ns min/avg/max=615/2000/24601 total=36880 us cpu=3.68%
isr_pack       n=2304 cyc min/avg/max=4972/5568/7093 ns min/avg/max=8286/9280/11821 total=21382 us cpu=2.13%
isr_dma_submit n=2304 cyc min/avg/max=465/560/1071 ns min/avg/max=775/933/1785 total=2150 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 26.106 ms/f, 41.85% of the peak window.
2. `fx_advance` — 2.219 ms/f, 3.56% of the peak window.
3. `fx_prepare_frame` — 0.684 ms/f, 1.10% of the peak window.
4. `fx_timeline_step` — 0.220 ms/f, 0.35% of the peak window.
5. `pov_preserve_half` — 0.142 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.

README cells: peak 🟢 30.419 (2), spilled 🟢 0/1647 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=MobiusGrid  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (MobiusGrid 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (209.9 ppm, frames 849-864)

INVALID: MobiusGrid
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (209.9 ppm, frames 849-864). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
