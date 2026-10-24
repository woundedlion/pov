# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/hankinsolids_ship.log`, captured 2026-10-23 17:16 America/Los_Angeles on COM4.
Replaces `profile_hankinsolids_teensy_2026-10-23.md` (the earlier `49bb199f4` capture, a different source and frame scope).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids, 288×144, single-entry playlist, source `bda58b65d3fc317f5d6ccbe29888014738a0ee25` |
| Method | 170000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HankinSolids profile 170 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:161776, data:171904, headers:8336   free for files:1689600` / `RAM1: variables:313632, code:31432, padding:1336   free for local variables:177888` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995632, data:851824, headers:9056   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `ed8d30cbeebcaa01298de47814dd4c8120227c606b198c818633d616cb696a2d`. Phantasm ELF SHA-256: `5fc8212484e6c1fb01b5305253408e170a8226b46fb0cbbaff78393b0a20b80f`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Capture provenance:

```text
profile provenance: version=1
profile provenance: source_sha=bda58b65d3fc317f5d6ccbe29888014738a0ee25
profile provenance: compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile provenance: profile_elf_sha256=ed8d30cbeebcaa01298de47814dd4c8120227c606b198c818633d616cb696a2d
profile provenance: phantasm_elf_sha256=5fc8212484e6c1fb01b5305253408e170a8226b46fb0cbbaff78393b0a20b80f
profile provenance: profile_envdump_sha256=74ccd5579e1674e0727a55d4e1a9730b8e93d92b1c5d48c193c3e837d9988e3d
profile provenance: phantasm_envdump_sha256=f5b172fbaa6ad2b8d5af681df338fa70985682fe91967d36b794d3b3f973c9ea
profile provenance: artifact_profile_elf=build/prof/artifacts/hankinsolids_ship_bda58b65d3fc_ed8d30cbeebc/profile.elf
profile provenance: artifact_phantasm_elf=build/prof/artifacts/hankinsolids_ship_bda58b65d3fc_ed8d30cbeebc/phantasm.elf
profile provenance: artifact_dir=build/prof/artifacts/hankinsolids_ship_bda58b65d3fc_ed8d30cbeebc
profile provenance: profile_image=built:ec35cf6b835353580239e47d866d394d
profile provenance: phantasm_image=built:de16925e8d887db981c5aec1e331c698
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **31.644 ms** (frame 2379), mean **17.134 ms**, **0/2687 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–2688 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **19.276 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 2697 raw frame rows, including that trailing portion, the observed peak is 31.644 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `hk_timeline_step` averages 16.863 ms/frame; its worst window averages 27.581 ms/frame (frames 2369–2384).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 56.752 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.94 | 0.24% |
| ↳ `hk_timeline_step` | 15.515 | 1.00 | 27.34% |
| ↳ ↳ `hk_update_hankin` DUPLICATE-NAME | 0.072 | 1.00 | 0.13% |
| ↳ ↳ `hk_draw_mesh` | 15.229 | 1.00 | 26.84% |
| ↳ ↳ ↳ `hk_mesh_scan` | 14.921 | 1.00 | 26.29% |
| ↳ ↳ ↳ `hk_mesh_transform` | 0.101 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.100 | 1.00 | 0.18% |
| ↳ `canvas_buffer_wait` | 40.962 | 1.00 | 72.18% |

Wall min/mean/max: 16.727/56.755/65.036 ms.

### Peak-render window (frames 2369–2384)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.828 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.23% |
| ↳ `hk_timeline_step` | 27.581 | 1.00 | 43.90% |
| ↳ ↳ `hk_update_hankin` DUPLICATE-NAME | 0.705 | 1.00 | 1.12% |
| ↳ ↳ `hk_draw_mesh` | 26.087 | 1.00 | 41.52% |
| ↳ ↳ ↳ `hk_mesh_scan` | 25.390 | 1.00 | 40.41% |
| ↳ ↳ ↳ `hk_mesh_transform` | 0.444 | 1.00 | 0.71% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 34.987 | 1.00 | 55.69% |

Wall min/mean/max: 53.011/62.832/71.535 ms.

### Ownership-change window (frames 2161–2176)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.986 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.147 | 1.00 | 0.23% |
| ↳ `hk_timeline_step` | 21.148 | 1.00 | 33.58% |
| ↳ ↳ `hk_conway_previous` | 0.021 | 0.19 | 0.03% |
| ↳ ↳ `hk_draw_mesh` DUPLICATE-NAME | 17.990 | 0.81 | 28.56% |
| ↳ ↳ ↳ `hk_mesh_scan` DUPLICATE-NAME | 17.784 | 0.81 | 28.23% |
| ↳ ↳ ↳ `hk_mesh_transform` DUPLICATE-NAME | 0.028 | 0.81 | 0.04% |
| ↳ ↳ `hk_conway_compile` | 0.011 | 0.19 | 0.02% |
| ↳ ↳ `hk_conway_sweep` | 0.029 | 0.19 | 0.05% |
| ↳ ↳ `hk_update_hankin` DUPLICATE-NAME | 0.365 | 0.81 | 0.58% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 41.581 | 1.00 | 66.02% |

Wall min/mean/max: 60.338/62.990/71.033 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `hk_timeline_step` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| truncatedIcosidodecahedron | truncatedIcosidodecahedron | 31.644 | 0/117 | 27.581 | 6/7 |
| truncatedDodecahedron | truncatedDodecahedron | 30.537 | 0/105 | 26.590 | 6/7 |
| snubDodecahedron | snubDodecahedron | 30.296 | 0/117 | 26.838 | 6/7 |
| truncatedCuboctahedron | truncatedCuboctahedron | 30.225 | 0/117 | 27.573 | 7/8 |
| truncatedCube | truncatedCube | 29.544 | 0/105 | 25.963 | 6/7 |
| snubCube | snubCube | 28.938 | 0/117 | 23.003 | 7/8 |
| rhombicosidodecahedron | rhombicosidodecahedron | 28.863 | 0/117 | 26.144 | 6/7 |
| truncatedIcosahedron | truncatedIcosahedron | 27.797 | 0/105 | 25.478 | 6/7 |
| icosidodecahedron | icosidodecahedron | 25.791 | 0/285 | 22.157 | 10/17 |
| truncatedTetrahedron | truncatedTetrahedron | 25.549 | 0/105 | 20.046 | 5/6 |
| truncatedOctahedron | truncatedOctahedron | 25.255 | 0/105 | 22.990 | 6/7 |
| rhombicuboctahedron | rhombicuboctahedron | 24.154 | 0/105 | 20.867 | 6/7 |
| icosahedron | icosahedron | 22.555 | 0/160 | 18.328 | 6/9 |
| dodecahedron | dodecahedron | 22.533 | 0/225 | 19.584 | 10/15 |
| cube | cube | 22.231 | 0/213 | 18.073 | 8/13 |
| octahedron | octahedron | 21.673 | 0/177 | 17.693 | 8/11 |
| cuboctahedron | cuboctahedron | 21.670 | 0/285 | 18.673 | 9/16 |
| 0 | — | 18.640 | 0/103 | 15.515 | 6/7 |
| tetrahedron | tetrahedron | 15.162 | 0/24 | 11.927 | 1/2 |

The strict marker check confirms 18 distinct named solids and a return to the first solid. The additional startup owner is unmarked and appears separately in the ownership table.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 1007.050 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1160.19 | 0.000578 | 0.002018 | 0.030755 | 3.72% |
| `isr_pack` | 145.00 | 0.008541 | 0.009315 | 0.011675 | 2.14% |
| `isr_dma_submit` | 145.00 | 0.000781 | 0.000935 | 0.006073 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `hk_timeline_step` — 27.581 ms/frame, 43.90% of this root window.
2. `pov_preserve_half` — 0.144 ms/frame, 0.23% of this root window.
3. `canvas_clear` — 0.090 ms/frame, 0.14% of this root window.

README cells: peak 🟢 31.644 (19), spilled 🟢 0/2687 (0.00%).

README setup frame 1: 19.276 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

Previous `49bb199f4` capture (different source and capture) trapped during an arena allocation; no complete peak, mean or spill aggregate was available. This capture includes the scratch-arena correction.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=HankinSolids  windows=168
  [PASS] captured >=3 windows (168)
  [PASS] every window names one effect and resolution (HankinSolids 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] shape markers advance and return to the first (cycle closed): 18 distinct
  [PASS] per-frame render is render, not wall (0 of 168 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 168 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hk_draw_mesh, hk_mesh_scan, hk_mesh_transform, hk_update_hankin, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'hk_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (65.8 ppm, frames 2369-2384)

INVALID: HankinSolids
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (65.8 ppm, frames 2369-2384). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
