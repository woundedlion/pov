# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, captured 2026-10-22 19:04 America/Los_Angeles on COM4.
Replaces `profile_gsreactiondiffusion_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GSReactionDiffusion profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:79944, data:343828, headers:8356   free for files:1599488` / `RAM1: variables:313248, code:17224, padding:15544   free for local variables:178272` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/gsreactiondiffusion_ship_aa61e787859b_f3c5fec88b69`. Profile ELF SHA-256: `f3c5fec88b696ae0141e971c4caa752313d7cb4bb43dede239c79ef9a9b3f2d0`.

## Frame cadence

**Pass aggregate**: peak frame render is **48.119 ms** (frame 395), mean live-frame render is **40.898 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`grd_render` averages 41.804 ms/f across complete post-setup windows; its worst window is 46.968 ms/f (frames 385–416).

Setup frame 1 rendered 24.850 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:15) recorded peak 🟢 39.43 and spilled 🟢 0/2047 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 385–416)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.378 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.135 | 1.00 | 0.22% |
| ↳ `grd_render` | 46.968 | 1.00 | 75.30% |
| ↳ ↳ `grd_shader_draw` | 28.332 | 6.00 | 45.42% |
| ↳ ↳ `grd_cull_flags` | 0.701 | 6.00 | 1.12% |
| ↳ ↳ `grd_publish` | 4.130 | 6.00 | 6.62% |
| ↳ ↳ `grd_physics` | 5.170 | 6.00 | 8.29% |
| ↳ ↳ `grd_pigment` | 4.242 | 1.00 | 6.80% |
| ↳ ↳ `rd_lattice_basis` | 0.142 | 1.00 | 0.23% |
| ↳ ↳ `grd_color_noise` | 2.835 | 1.00 | 4.55% |
| ↳ `rd_timeline_step` | 0.051 | 1.00 | 0.08% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.127 | 1.00 | 24.25% |

Wall min/avg/max = 61.786/62.379/63.073 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1998630 us
isr_wake       n=36841 cyc min/avg/max=301/1219/20882 ns min/avg/max=501/2031/34803 total=74889 us cpu=3.74%
isr_pack       n=4605 cyc min/avg/max=4730/5538/8488 ns min/avg/max=7883/9230/14146 total=42505 us cpu=2.12%
isr_dma_submit n=4605 cyc min/avg/max=365/565/2472 ns min/avg/max=608/941/4120 total=4344 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `grd_render` — 46.968 ms/f, 75.30% of the peak window.
2. `pov_preserve_half` — 0.135 ms/f, 0.22% of the peak window.
3. `canvas_clear` — 0.087 ms/f, 0.14% of the peak window.
4. `rd_timeline_step` — 0.051 ms/f, 0.08% of the peak window.

README cells: peak 🟢 48.119, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=GSReactionDiffusion  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GSReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'grd_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (28.0 ppm, frames 385-416)

INVALID: GSReactionDiffusion
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (28.0 ppm, frames 385-416). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
