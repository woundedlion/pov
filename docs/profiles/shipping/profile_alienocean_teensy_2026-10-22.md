# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/alienocean_ship.log`, captured 2026-10-22 19:51 America/Los_Angeles on COM4.
Replaces `profile_alienocean_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienOcean profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:86496, data:154492, headers:8868   free for files:1781760` / `RAM1: variables:313024, code:16872, padding:15896   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/alienocean_ship_aa61e787859b_6c5ddc6dbdd7`. Profile ELF SHA-256: `6c5ddc6dbdd7e0f8135e95a4e8210d7c8e6f259667af297920940708d49eb797`.

## Frame cadence

**Pass aggregate**: peak frame render is **25.517 ms** (frame 55), mean live-frame render is **25.172 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 20.864 ms/f across complete post-setup windows; its worst window is 20.936 ms/f (frames 161–192).

Setup frame 1 rendered 45.694 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:39) recorded peak 🟢 22.37 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 33–64)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.395 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 20.880 | 1.00 | 33.46% |
| ↳ ↳ `fx_prepare_column` | 1.238 | 144.00 | 1.98% |
| ↳ `fx_prepare_frame` | 0.606 | 1.00 | 0.97% |
| ↳ `fx_advance` | 3.031 | 1.00 | 4.86% |
| ↳ `fx_timeline_step` | 0.118 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 37.126 | 1.00 | 59.50% |

Wall min/avg/max = 61.919/62.400/62.914 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000165 us
isr_wake       n=36869 cyc min/avg/max=318/1207/14482 ns min/avg/max=530/2011/24136 total=74195 us cpu=3.70%
isr_pack       n=4608 cyc min/avg/max=5295/5730/7280 ns min/avg/max=8825/9550/12133 total=44013 us cpu=2.20%
isr_dma_submit n=4608 cyc min/avg/max=434/569/3537 ns min/avg/max=723/948/5895 total=4370 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 20.880 ms/f, 33.46% of the peak window.
2. `fx_advance` — 3.031 ms/f, 4.86% of the peak window.
3. `fx_prepare_frame` — 0.606 ms/f, 0.97% of the peak window.
4. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.118 ms/f, 0.19% of the peak window.
6. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.

README cells: peak 🟢 25.517, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=AlienOcean  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienOcean 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (89.7 ppm, frames 1-32)

INVALID: AlienOcean
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (89.7 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
