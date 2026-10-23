# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/islamicstars_ship.log`, captured 2026-10-22 19:39 America/Los_Angeles on COM4.
Replaces `profile_islamicstars_teensy_2026-10-07.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 306 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh IslamicStars profile 306 16` |

Image size (`profile`, this effect only): `FLASH: code:186968, data:207868, headers:8620   free for files:1628160` / `RAM1: variables:313632, code:31544, padding:1224   free for local variables:177888` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/islamicstars_ship_aa61e787859b_d89ed98e0888`. Profile ELF SHA-256: `d89ed98e0888316d7ad966a747f7342b092b004572143cc7347bf91f7b6e757b`.

## Frame cadence

**Pass aggregate**: peak frame render is **67.393 ms** (frame 2056), mean live-frame render is **27.035 ms**, and **33/4799** frames spilled (0.69%). Setup frame 1 is excluded from those figures.
`is_timeline_step` averages 26.771 ms/f across complete post-setup windows; its worst window is 64.701 ms/f (frames 2049–2064).

Setup frame 1 rendered 14.636 ms. Published live-frame figures use complete counter windows through frame 4800; a trailing partial window is excluded.

The previous shipping report (2026-10-07 15:16) recorded peak 🟢 48.510 (23) and spilled 🟢 0/3327 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 2049–2064)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| ↳ ↳ ↳ ↳ `scan_face_setup` DUPLICATE-NAME | 7.818 | 1149.44 | 4.36% |
| `mesh_copy_capture` MIXED-PARENT | 0.210 | 22.19 | 0.12% |
| `mesh_capture` MIXED-PARENT | 2.670 | 5.00 | 1.49% |
| `mesh_classify` MIXED-PARENT | 1.528 | 1.00 | 0.85% |
| `frame` | 179.466 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.08% |
| ↳ `is_timeline_step` | 64.701 | 1.00 | 36.05% |
| ↳ ↳ `is_draw_shape` | 64.505 | 1.00 | 35.94% |
| ↳ ↳ ↳ `is_mesh_scan` | 53.528 | 1.00 | 29.83% |
| ↳ ↳ ↳ ↳ `scan_face_sweep` | 5.117 | 1.00 | 2.85% |
| ↳ ↳ ↳ ↳ ↳ `scan_sweep_column` | 1.973 | 143.25 | 1.10% |
| ↳ ↳ ↳ ↳ `mesh_sweep_paths` | 0.116 | 1.00 | 0.06% |
| ↳ ↳ ↳ ↳ `scan_mesh_raster` | 34.161 | 1097.69 | 19.03% |
| ↳ ↳ ↳ `is_face_offsets` | 0.560 | 1.00 | 0.31% |
| ↳ ↳ ↳ `is_mesh_transform` | 4.110 | 1.00 | 2.29% |
| ↳ `is_ripple_prepare` | 0.019 | 1.00 | 0.01% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.05% |
| ↳ `canvas_buffer_wait` | 114.442 | 1.00 | 63.77% |

Wall min/avg/max = 66.137/179.483/189.517 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `is_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 67.393 | 33/352 | 64.701 | 21/22 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 57.235 | 0/288 | 54.333 | 17/18 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 54.216 | 0/312 | 47.321 | 18/19 |
| truncatedCube_hk90_meta | V=24 E=36 F=14 I=72 | 51.169 | 0/368 | 48.157 | 22/23 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | V=12 E=30 F=20 I=60 | 48.031 | 0/225 | 45.052 | 14/14 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 42.278 | 0/280 | 25.601 | 17/18 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 40.215 | 0/352 | 36.714 | 21/22 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 39.951 | 0/288 | 37.904 | 17/18 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 39.924 | 0/312 | 37.486 | 18/19 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 38.245 | 0/248 | 36.525 | 14/15 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 38.055 | 0/272 | 29.938 | 16/17 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 37.607 | 0/352 | 34.057 | 21/22 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 36.522 | 0/280 | 33.985 | 17/18 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 30.234 | 0/312 | 25.985 | 19/20 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 29.477 | 0/280 | 28.111 | 16/17 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 28.847 | 0/278 | 27.870 | 16/17 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2873995 us
isr_wake       n=52977 cyc min/avg/max=343/1208/16717 ns min/avg/max=571/2013/27861 total=106711 us cpu=3.71%
isr_pack       n=6621 cyc min/avg/max=5177/5581/7899 ns min/avg/max=8628/9301/13165 total=61590 us cpu=2.14%
isr_dma_submit n=6621 cyc min/avg/max=434/565/3769 ns min/avg/max=723/941/6281 total=6234 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `is_timeline_step` — 64.701 ms/f, 36.05% of the peak window.
2. `pov_preserve_half` — 0.140 ms/f, 0.08% of the peak window.
3. `canvas_clear` — 0.097 ms/f, 0.05% of the peak window.
4. `is_ripple_prepare` — 0.019 ms/f, 0.01% of the peak window.

README cells: peak 🔴 67.393 (1)<br>🟢 57.235 (15), spilled 🔴 33/352 (9.38%)<br>🟢 0/4447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=IslamicStars  windows=300
  [PASS] captured >=3 windows (300)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] shape markers advance and return to the first (cycle closed): 16 distinct
  [PASS] per-frame render is render, not wall (0 of 300 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 300 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): is_mesh_transform, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (96.7 ppm, frames 2049-2064)

INVALID: IslamicStars
```

Validation limitations: shape markers advance and return to the first (cycle closed): 16 distinct; root cyc/600 vs wall sum within 5 ppm (96.7 ppm, frames 2049-2064). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
