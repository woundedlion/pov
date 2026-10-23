# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/hankinsolids_ship.log`, captured 2026-10-22 19:23 America/Los_Angeles on COM4.
Replaces `profile_hankinsolids_teensy_2026-10-07.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 170 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HankinSolids profile 170 16 "-D HS_PROFILE_ORDERED_CYCLE"` |

Image size (`profile`, this effect only): `FLASH: code:154240, data:171688, headers:8920   free for files:1696768` / `RAM1: variables:313632, code:29352, padding:3416   free for local variables:177888` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/hankinsolids_ship_aa61e787859b_ee13957368bf`. Profile ELF SHA-256: `ee13957368bfd092118659419f50f365aa3c6359a171f5a469beb5d1684be6c2`.

## Frame cadence

**Pass aggregate**: peak frame render is **33.483 ms** (frame 2094), mean live-frame render is **16.354 ms**, and **0/2687** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`hk_timeline_step` averages 16.089 ms/f across complete post-setup windows; its worst window is 27.439 ms/f (frames 2097–2112).

Setup frame 1 rendered 14.951 ms. Published live-frame figures use complete counter windows through frame 2688; a trailing partial window is excluded.

The previous shipping report (2026-10-07 14:18) recorded peak 🟢 34.681 (19) and spilled 🟢 0/4447 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 2081–2096)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.297 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.147 | 1.00 | 0.24% |
| ↳ `hk_timeline_step` | 27.220 | 1.00 | 43.69% |
| ↳ ↳ `hk_update_hankin` DUPLICATE-NAME | 0.671 | 1.00 | 1.08% |
| ↳ ↳ `hk_draw_mesh` | 25.766 | 1.00 | 41.36% |
| ↳ ↳ ↳ `hk_mesh_scan` | 25.038 | 1.00 | 40.19% |
| ↳ ↳ ↳ ↳ `scan_face_setup` DUPLICATE-NAME | 1.681 | 197.38 | 2.70% |
| ↳ ↳ ↳ ↳ `mesh_copy_capture` | 0.061 | 8.38 | 0.10% |
| ↳ ↳ ↳ ↳ `mesh_capture` | 0.160 | 4.38 | 0.26% |
| ↳ ↳ ↳ ↳ `mesh_classify` | 0.224 | 0.88 | 0.36% |
| ↳ ↳ ↳ ↳ `scan_mesh_raster` | 22.215 | 198.62 | 35.66% |
| ↳ ↳ ↳ `hk_mesh_transform` | 0.509 | 1.00 | 0.82% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 34.814 | 1.00 | 55.88% |

Wall min/avg/max = 55.209/62.301/64.438 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `hk_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron | truncatedIcosidodecahedron | 33.483 | 0/117 | 27.439 | 6/7 |
| snubDodecahedron | snubDodecahedron | 30.340 | 0/61 | 23.884 | 3/3 |
| snubCube | snubCube | 26.757 | 0/117 | 22.245 | 7/8 |
| truncatedCuboctahedron | truncatedCuboctahedron | 26.603 | 0/117 | 23.032 | 6/7 |
| truncatedCube | truncatedCube | 26.217 | 0/105 | 22.563 | 6/7 |
| rhombicosidodecahedron | rhombicosidodecahedron | 26.063 | 0/117 | 23.245 | 6/7 |
| truncatedIcosahedron | truncatedIcosahedron | 25.316 | 0/105 | 22.313 | 6/7 |
| truncatedTetrahedron | truncatedTetrahedron | 22.890 | 0/105 | 18.348 | 6/7 |
| truncatedOctahedron | truncatedOctahedron | 22.860 | 0/105 | 19.668 | 6/7 |
| icosidodecahedron | icosidodecahedron | 22.239 | 0/222 | 20.268 | 12/14 |
| rhombicuboctahedron | rhombicuboctahedron | 22.025 | 0/105 | 19.369 | 6/7 |
| dodecahedron | dodecahedron | 21.189 | 0/234 | 18.911 | 13/15 |
| icosahedron | icosahedron | 21.036 | 0/105 | 18.721 | 6/7 |
| octahedron | octahedron | 20.876 | 0/210 | 17.103 | 11/13 |
| cuboctahedron | cuboctahedron | 20.041 | 0/327 | 18.019 | 17/20 |
| cube | cube | 19.500 | 0/327 | 16.396 | 16/19 |
| tetrahedron | tetrahedron | 16.351 | 0/105 | 13.586 | 5/6 |
| 0 | — | 16.335 | 0/103 | 13.584 | 5/6 |

Startup bucket 0 is the initial tetrahedron state. Named buckets merge repeated visits to the same shape, including bridge endpoints.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=998753 us
isr_wake       n=18410 cyc min/avg/max=375/1196/13605 ns min/avg/max=625/1993/22675 total=36721 us cpu=3.67%
isr_pack       n=2301 cyc min/avg/max=5229/5670/7072 ns min/avg/max=8715/9450/11786 total=21746 us cpu=2.17%
isr_dma_submit n=2301 cyc min/avg/max=442/564/2480 ns min/avg/max=736/940/4133 total=2164 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `hk_timeline_step` — 27.220 ms/f, 43.69% of the peak window.
2. `pov_preserve_half` — 0.147 ms/f, 0.24% of the peak window.
3. `canvas_clear` — 0.093 ms/f, 0.15% of the peak window.

README cells: peak 🟢 33.483 (18), spilled 🟢 0/2687 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=HankinSolids  windows=168
  [PASS] captured >=3 windows (168)
  [PASS] every window names one effect and resolution (HankinSolids 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] shape markers advance and return to the first (cycle closed): 17 distinct
  [PASS] per-frame render is render, not wall (0 of 168 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 168 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hk_draw_mesh, hk_mesh_scan, hk_mesh_transform, hk_update_hankin, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'hk_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (45.6 ppm, frames 2097-2112)

INVALID: HankinSolids
```

Validation limitations: shape markers advance and return to the first (cycle closed): 17 distinct; root cyc/600 vs wall sum within 5 ppm (45.6 ppm, frames 2097-2112). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
