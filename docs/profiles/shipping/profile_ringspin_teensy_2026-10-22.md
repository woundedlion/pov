# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/ringspin_ship.log`, captured 2026-10-22 19:12 America/Los_Angeles on COM4.
Replaces `profile_ringspin_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh RingSpin profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:91568, data:150960, headers:8352   free for files:1780736` / `RAM1: variables:312960, code:27448, padding:5320   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/ringspin_ship_aa61e787859b_b2463e2e54ad`. Profile ELF SHA-256: `b2463e2e54ad64c0f2a6209a32e742ba284e743562a72e189f049688255fbb92`.

## Frame cadence

**Pass aggregate**: peak frame render is **53.562 ms** (frame 761), mean live-frame render is **38.838 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`rs_draw_rings` averages 38.940 ms/f across complete post-setup windows; its worst window is 41.441 ms/f (frames 737–768).

Setup frame 1 rendered 3.379 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:22) recorded peak 🟢 47.16 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 737–768)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.611 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `rs_draw_rings` | 41.441 | 1.00 | 66.19% |
| ↳ ↳ `rs_ring_scan` | 40.265 | 76.00 | 64.31% |
| ↳ `rs_timeline_step` | 0.100 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.088 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 20.833 | 1.00 | 33.27% |

Wall min/avg/max = 43.921/62.612/80.940 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2005153 us
isr_wake       n=36961 cyc min/avg/max=335/1138/12729 ns min/avg/max=558/1896/21215 total=70111 us cpu=3.49%
isr_pack       n=4620 cyc min/avg/max=3666/5182/6967 ns min/avg/max=6110/8636/11611 total=39908 us cpu=1.99%
isr_dma_submit n=4620 cyc min/avg/max=380/561/654 ns min/avg/max=633/935/1090 total=4325 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `rs_draw_rings` — 41.441 ms/f, 66.19% of the peak window.
2. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
3. `rs_timeline_step` — 0.100 ms/f, 0.16% of the peak window.
4. `canvas_clear` — 0.088 ms/f, 0.14% of the peak window.

README cells: peak 🟢 53.562, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=RingSpin  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (RingSpin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'rs_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (36.6 ppm, frames 737-768)

INVALID: RingSpin
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (36.6 ppm, frames 737-768). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
