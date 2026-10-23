# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gridspace_ship.log`, captured 2026-10-22 19:55 America/Los_Angeles on COM4.
Current shipping timing; the [October 17 comparison](profile_gridspace_teensy_2026-10-17.md) retains the matched temporal-on/off measurements.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GridSpace profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:87520, data:155216, headers:9168   free for files:1779712` / `RAM1: variables:313024, code:17208, padding:15560   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/gridspace_ship_aa61e787859b_3574a0d3e833`. Profile ELF SHA-256: `3574a0d3e8337de1fcf7b4e22b72f40ee353b0556870d9b7ebb741c275c44eeb`.

## Frame cadence

**Pass aggregate**: peak frame render is **25.957 ms** (frame 421), mean live-frame render is **25.260 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 18.253 ms/f across complete post-setup windows; its worst window is 18.405 ms/f (frames 577–608).

Setup frame 1 rendered 39.297 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-17 16:58) recorded peak 🟢 25.41 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 417–448)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.391 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.146 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 18.364 | 1.00 | 29.43% |
| ↳ ↳ `fx_prepare_column` | 0.775 | 144.00 | 1.24% |
| ↳ `fx_prepare_frame` | 3.990 | 1.00 | 6.39% |
| ↳ `fx_advance` | 2.502 | 1.00 | 4.01% |
| ↳ `fx_timeline_step` | 0.121 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 36.818 | 1.00 | 59.01% |

Wall min/avg/max = 61.618/62.406/63.367 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000683 us
isr_wake       n=36879 cyc min/avg/max=323/1173/10542 ns min/avg/max=538/1955/17570 total=72138 us cpu=3.60%
isr_pack       n=4610 cyc min/avg/max=3919/5450/7018 ns min/avg/max=6531/9083/11696 total=41876 us cpu=2.09%
isr_dma_submit n=4610 cyc min/avg/max=392/563/4751 ns min/avg/max=653/938/7918 total=4328 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 18.364 ms/f, 29.43% of the peak window.
2. `fx_prepare_frame` — 3.990 ms/f, 6.39% of the peak window.
3. `fx_advance` — 2.502 ms/f, 4.01% of the peak window.
4. `pov_preserve_half` — 0.146 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.121 ms/f, 0.19% of the peak window.
6. `canvas_clear` — 0.098 ms/f, 0.16% of the peak window.

README cells: peak 🟢 25.957, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=GridSpace  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (GridSpace 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (230.5 ppm, frames 1-32)

INVALID: GridSpace
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (230.5 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
