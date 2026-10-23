# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/chromaticlichen_ship.log`, captured 2026-10-22 20:02 America/Los_Angeles on COM4.
Replaces `profile_chromaticlichen_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ChromaticLichen profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:88488, data:154524, headers:8892   free for files:1779712` / `RAM1: variables:313024, code:16984, padding:15784   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/chromaticlichen_ship_aa61e787859b_60cc2e4ddd62`. Profile ELF SHA-256: `60cc2e4ddd621a80564f547e193a98854941ade1703e87cbd9ee1075921c0d29`.

## Frame cadence

**Pass aggregate**: peak frame render is **40.095 ms** (frame 272), mean live-frame render is **39.722 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 33.003 ms/f across complete post-setup windows; its worst window is 33.045 ms/f (frames 193–224).

Setup frame 1 rendered 65.439 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:57) recorded peak 🟢 35.53 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 257–288)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.364 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 33.020 | 1.00 | 52.95% |
| ↳ ↳ `fx_prepare_column` | 1.717 | 144.00 | 2.75% |
| ↳ `fx_prepare_frame` | 4.347 | 1.00 | 6.97% |
| ↳ `fx_advance` | 1.817 | 1.00 | 2.91% |
| ↳ `fx_timeline_step` | 0.135 | 1.00 | 0.22% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 22.453 | 1.00 | 36.00% |

Wall min/avg/max = 62.130/62.384/62.648 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000022 us
isr_wake       n=36867 cyc min/avg/max=292/1173/14481 ns min/avg/max=486/1955/24135 total=72086 us cpu=3.60%
isr_pack       n=4608 cyc min/avg/max=4963/5461/6886 ns min/avg/max=8271/9101/11476 total=41945 us cpu=2.09%
isr_dma_submit n=4608 cyc min/avg/max=384/558/665 ns min/avg/max=640/930/1108 total=4288 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 33.020 ms/f, 52.95% of the peak window.
2. `fx_prepare_frame` — 4.347 ms/f, 6.97% of the peak window.
3. `fx_advance` — 1.817 ms/f, 2.91% of the peak window.
4. `pov_preserve_half` — 0.140 ms/f, 0.22% of the peak window.
5. `fx_timeline_step` — 0.135 ms/f, 0.22% of the peak window.
6. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.

README cells: peak 🟢 40.095, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=ChromaticLichen  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (ChromaticLichen 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (283.5 ppm, frames 1-32)

INVALID: ChromaticLichen
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (283.5 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
