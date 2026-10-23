# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gnomonicstars_ship.log`, captured 2026-10-22 19:03 America/Los_Angeles on COM4.
Replaces `profile_gnomonicstars_teensy_2026-10-07.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GnomonicStars profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:75408, data:152040, headers:9096   free for files:1795072` / `RAM1: variables:313024, code:20776, padding:11992   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/gnomonicstars_ship_aa61e787859b_bf4194573a08`. Profile ELF SHA-256: `bf4194573a0830a74cdfabb9c7139c105a7facf86a1aab1783165104de567e48`.

## Frame cadence

**Pass aggregate**: peak frame render is **26.284 ms** (frame 206), mean live-frame render is **12.357 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`gn_draw_stars` averages 12.147 ms/f across complete post-setup windows; its worst window is 18.360 ms/f (frames 193–224).

Setup frame 1 rendered 15.394 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-07 13:42) recorded peak 🟢 21.84 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 193–224)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.480 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `gn_draw_stars` | 18.360 | 1.00 | 29.39% |
| ↳ ↳ `gn_star_scan` | 16.708 | 238.00 | 26.74% |
| ↳ `gn_timeline_step` | 0.098 | 1.00 | 0.16% |
| ↳ ↳ `animation_mobius_step` | 0.035 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.085 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 43.785 | 1.00 | 70.08% |

Wall min/avg/max = 48.204/62.484/75.703 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2002514 us
isr_wake       n=36913 cyc min/avg/max=329/1038/13849 ns min/avg/max=548/1730/23081 total=63895 us cpu=3.19%
isr_pack       n=4614 cyc min/avg/max=3676/4416/5968 ns min/avg/max=6126/7360/9946 total=33965 us cpu=1.69%
isr_dma_submit n=4614 cyc min/avg/max=361/557/697 ns min/avg/max=601/928/1161 total=4290 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `gn_draw_stars` — 18.360 ms/f, 29.39% of the peak window.
2. `pov_preserve_half` — 0.142 ms/f, 0.23% of the peak window.
3. `gn_timeline_step` — 0.098 ms/f, 0.16% of the peak window.
4. `canvas_clear` — 0.085 ms/f, 0.14% of the peak window.

README cells: peak 🟢 26.284, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=GnomonicStars  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GnomonicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'gn_draw_stars' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (62.3 ppm, frames 193-224)

INVALID: GnomonicStars
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (62.3 ppm, frames 193-224). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
