# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/bzreactiondiffusion_ship.log`, captured 2026-10-22 18:57 America/Los_Angeles on COM4.
Replaces `profile_bzreactiondiffusion_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh BZReactionDiffusion profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:68320, data:241624, headers:8520   free for files:1713152` / `RAM1: variables:313024, code:16168, padding:16600   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/bzreactiondiffusion_ship_aa61e787859b_ef9791eb28cb`. Profile ELF SHA-256: `ef9791eb28cb947ef7111b9cdf480b9a3ccd3fe7d1f0dac64c1a42e0946b31f8`.

## Frame cadence

**Pass aggregate**: peak frame render is **48.271 ms** (frame 643), mean live-frame render is **46.525 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`bz_render` averages 46.242 ms/f across complete post-setup windows; its worst window is 46.914 ms/f (frames 33–64).

Setup frame 1 rendered 78.475 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:05) recorded peak 🟢 48.34 and spilled 🟢 0/2047 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 641–672)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.420 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.135 | 1.00 | 0.22% |
| ↳ `bz_render` | 46.829 | 1.00 | 75.02% |
| ↳ ↳ `bz_raster` | 42.206 | 2.00 | 67.62% |
| ↳ ↳ `bz_physics` | 4.531 | 2.00 | 7.26% |
| ↳ ↳ `rd_lattice_basis` | 0.090 | 1.00 | 0.14% |
| ↳ `rd_timeline_step` | 0.038 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.086 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.328 | 1.00 | 24.56% |

Wall min/avg/max = 60.685/62.444/64.095 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999782 us
isr_wake       n=36862 cyc min/avg/max=289/1203/13602 ns min/avg/max=481/2005/22670 total=73935 us cpu=3.69%
isr_pack       n=4607 cyc min/avg/max=5364/5785/7700 ns min/avg/max=8940/9641/12833 total=44421 us cpu=2.22%
isr_dma_submit n=4607 cyc min/avg/max=391/559/688 ns min/avg/max=651/931/1146 total=4295 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `bz_render` — 46.829 ms/f, 75.02% of the peak window.
2. `pov_preserve_half` — 0.135 ms/f, 0.22% of the peak window.
3. `canvas_clear` — 0.086 ms/f, 0.14% of the peak window.
4. `rd_timeline_step` — 0.038 ms/f, 0.06% of the peak window.

README cells: peak 🟢 48.271, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=BZReactionDiffusion  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (BZReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'bz_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (372.9 ppm, frames 1-32)

INVALID: BZReactionDiffusion
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (372.9 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
