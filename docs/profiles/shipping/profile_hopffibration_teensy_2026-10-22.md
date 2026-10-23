# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/hopffibration_ship.log`, captured 2026-10-22 19:05 America/Los_Angeles on COM4.
Replaces `profile_hopffibration_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HopfFibration profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:99752, data:148224, headers:9048   free for files:1774592` / `RAM1: variables:312992, code:29304, padding:3464   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/hopffibration_ship_aa61e787859b_71611fe90212`. Profile ELF SHA-256: `71611fe90212d95693d67987a1ef52589b55dc0c1bdc74902ebb15a05db1c571`.

## Frame cadence

**Pass aggregate**: peak frame render is **57.358 ms** (frame 824), mean live-frame render is **28.413 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`hf_render_trails` averages 28.436 ms/f across complete post-setup windows; its worst window is 48.511 ms/f (frames 193–224).

Setup frame 1 rendered 0.671 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:03) recorded peak 🟢 48.61 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 801–832)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.733 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 1.00 | 0.22% |
| ↳ `hf_render_trails` | 41.055 | 1.00 | 65.44% |
| ↳ `hf_project_record` | 0.186 | 1.00 | 0.30% |
| ↳ `hf_advance_tumble` | 0.002 | 1.00 | 0.00% |
| ↳ `hf_timeline_step` | 0.074 | 1.00 | 0.12% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 21.159 | 1.00 | 33.73% |

Wall min/avg/max = 31.353/62.746/95.503 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2010811 us
isr_wake       n=37065 cyc min/avg/max=364/1098/13725 ns min/avg/max=606/1830/22875 total=67836 us cpu=3.37%
isr_pack       n=4633 cyc min/avg/max=3661/4852/7024 ns min/avg/max=6101/8086/11706 total=37465 us cpu=1.86%
isr_dma_submit n=4633 cyc min/avg/max=410/563/643 ns min/avg/max=683/938/1071 total=4350 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `hf_render_trails` — 41.055 ms/f, 65.44% of the peak window.
2. `hf_project_record` — 0.186 ms/f, 0.30% of the peak window.
3. `pov_preserve_half` — 0.138 ms/f, 0.22% of the peak window.
4. `canvas_clear` — 0.087 ms/f, 0.14% of the peak window.
5. `hf_timeline_step` — 0.074 ms/f, 0.12% of the peak window.
6. `hf_advance_tumble` — 0.002 ms/f, 0.00% of the peak window.

README cells: peak 🟢 57.358, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=HopfFibration  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (HopfFibration 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'hf_render_trails' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (194.8 ppm, frames 193-224)

INVALID: HopfFibration
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (194.8 ppm, frames 193-224). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HopfFibration`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
