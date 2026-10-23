# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopepentbright_ship.log`, captured 2026-10-22 20:07 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopepentbright_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopePentBright 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopePentBright profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:86904, data:155052, headers:8924   free for files:1780736` / `RAM1: variables:313024, code:16968, padding:15800   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopepentbright_ship_aa61e787859b_cd3b1e466d82`. Profile ELF SHA-256: `cd3b1e466d82a91d1d3464e5794b78a479a5528e7bd5102025dc586c88076594`.

## Frame cadence

**Pass aggregate**: peak frame render is **28.723 ms** (frame 870), mean live-frame render is **26.889 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 23.248 ms/f across complete post-setup windows; its worst window is 24.575 ms/f (frames 417–448).

Setup frame 1 rendered 49.297 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:44) recorded peak 🟢 27.04 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 865–896)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.357 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 24.196 | 1.00 | 38.80% |
| ↳ ↳ `fx_prepare_column` | 1.164 | 144.00 | 1.87% |
| ↳ `fx_prepare_frame` | 0.923 | 1.00 | 1.48% |
| ↳ `fx_advance` | 2.032 | 1.00 | 3.26% |
| ↳ `fx_timeline_step` | 0.123 | 1.00 | 0.20% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.451 | 1.00 | 55.25% |

Wall min/avg/max = 60.824/62.379/63.987 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999819 us
isr_wake       n=36863 cyc min/avg/max=330/1201/18245 ns min/avg/max=550/2001/30408 total=73792 us cpu=3.68%
isr_pack       n=4607 cyc min/avg/max=5185/5738/7275 ns min/avg/max=8641/9563/12125 total=44065 us cpu=2.20%
isr_dma_submit n=4607 cyc min/avg/max=408/563/1950 ns min/avg/max=680/938/3250 total=4329 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 24.196 ms/f, 38.80% of the peak window.
2. `fx_advance` — 2.032 ms/f, 3.26% of the peak window.
3. `fx_prepare_frame` — 0.923 ms/f, 1.48% of the peak window.
4. `pov_preserve_half` — 0.144 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.123 ms/f, 0.20% of the peak window.
6. `canvas_clear` — 0.092 ms/f, 0.15% of the peak window.

README cells: peak 🟢 28.723, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopePentBright  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (KaleidoscopePentBright 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (350.7 ppm, frames 417-448)

INVALID: KaleidoscopePentBright
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (350.7 ppm, frames 417-448). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopePentBright`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
