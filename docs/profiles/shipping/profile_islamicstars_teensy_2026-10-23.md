# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw captures: `C:/work/temp/Holosphere-isopt/build/prof/isopt_final_carousel.log`
(2026-10-23 08:04 America/Los_Angeles) and
`C:/work/temp/Holosphere-isopt/build/prof/isopt_final_heavy6_replay.log`
(2026-10-23 08:07, COM4). Replaces the October 22 shipping report.

The retained worst observed render is **61.337 ms**, with **0/535 spills** in
the fixed-shape stress replay. The complete carousel separately reaches
**53.742 ms**, with **0/3327 spills**. The stress result is the headline and
index ranking value; it exceeds the 58 ms optimization goal by **3.337 ms**
and retains **1.163 ms** below the 62.5 ms display budget.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` face construction, SDF distance and scan regions; no new discretionary ITCM pins |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, captured source `12ba767003ae3b5a24590e0b5879309bf678b357` |
| Method | Carousel: 210 s, 16-frame windows, `HS_PROFILE_TRANS_SPEED=4`, `HS_PROFILE_EPOCH_REVS=1920`. Stress: 35 s, 8-frame windows, zero-based `HS_PROFILE_PRESET=6`, `HS_PROFILE_TRANS_SPEED=1`, `HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-isopt HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh IslamicStars profile 210 16 '-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920'` |

Stress reproduction uses the same wrapper with arguments
`IslamicStars profile 35 8 '-D HS_PROFILE_PRESET=6 -D HS_PROFILE_TRANS_SPEED=1 -D HS_PROFILE_EPOCH_REVS=1200'`.
Transition speed changes the sampled
motion schedule; the two captures are distinct scenarios, not interchangeable
baselines. No master baseline was captured.

Image size (`profile`, carousel): `FLASH: code:196912, data:207880, headers:8904` /
`RAM1: variables:313664, code:30904, padding:1864, free:177856` /
`RAM2: variables:520064, free:4224`. The fixed-shape profile has FLASH code
196992 B, data 207984 B and headers 8720 B; its RAM figures match.

Full Phantasm image: `FLASH: code:926240, data:851668, headers:8972` /
`RAM1: variables:312736, code:170520, padding:26088, free:14944` /
`RAM2: variables:520064, free:4224`. Against the preceding master build,
RAM1 code rises **1.184 KB** (169336 to 170520 B), variables remain 312736 B,
FLASH code rises **9.600 KB**, and FLASH data rises **0.004 KB**.
The `.text.itcm` section rises **1.184 KB** (169328 to 170512 B); this is code
growth without additional pinned functions. KB here means 1000 B.

Provenance records retain compiler attestation, envdumps, build logs and ELF
artifacts under the captures' matching `build/prof/artifacts/` directories.
Profile ELF SHA-256: carousel `5b54e23ecdab251b1fdfb1ea8b781341c608a5a66ecc4b9977fa98ab1a5217f1`;
stress `edef53b47b1ff2066f292b3ca6bf222bc1bb6944c802796d2c7a4942150e1410`. Both use Phantasm ELF SHA-256
`0417a4c1cba80b4cce1aaa2de146c9fd458a463b683ace8cb01c647360d19b31`.

Exactness cross-check: parser-selected carousel frames 1473–1488 have root
cycles/600 = 979778.033 us
against a 980067 us wall sum (**294.8 ppm**).
Stress frames 513–520 measure 506484.827 us
against 506627 us (**280.6 ppm**). Both exceed the strict
5 ppm limit. Logged per-frame render telemetry supports the stated peaks;
neither capture is claimed fully cycle-validated.

## Frame cadence

| Scenario | Peak frame render | Mean live render | Spilled/live frames | Complete windows |
|---|---:|---:|---:|---:|
| Fixed-shape stress | 61.337 ms, frame 520 | 38.383 ms | 0/535 | 67 |
| All-shape carousel | 53.742 ms, frame 525 | 25.041 ms | 0/3327 | 208 |

Setup frame 1 renders 16.942 ms in stress and
18.793 ms in the carousel; it is excluded from live
peaks and spill counts. The carousel has nine additional trailing frame rows
without a complete counter window; published aggregates end at frame 3328.
The stress capture ends at complete frame 536.

`is_timeline_step` averages 38.452 ms/frame in stress
and 24.718 ms/frame in the carousel across complete
post-setup windows. Its worst window is 57.233 ms/frame (stress 513–520) and
41.202 ms/frame (carousel 1473–1488).

A display window is 62.5 ms at 480 RPM. IslamicStars renders one quadrant of
the 288×144 canvas, about 10368 pixels; all observed live renders fit one
display window (16 fps render budget). Wall telemetry includes scheduling,
instrumentation and synchronization jitter. `canvas_buffer_wait` is the idle
wait for the next display flip.

## Phase-by-phase readout

The carousel constructs each shape through its Conway operations, then draws
the shape through its segue envelope and ripple motion before advancing.
The stress capture holds the selected heavy shape and exercises its motion.
Counter trees are inclusive per-frame averages; tagged duplicate or mixed
parents are not exclusive costs of the displayed branch.

### Fixed heavy shape under deformation (frames 513–520)

```text
frame                       63.311 ms 37.986 Mcyc 100.0%
  pov_preserve_half          140.1 us   84.1 kcyc  0.2%
  is_timeline_step          57.233 ms 34.340 Mcyc 90.4%
    is_draw_shape           57.068 ms 34.241 Mcyc 90.1%
      is_mesh_scan          47.510 ms 28.506 Mcyc 75.0%
        scan_face_sweep      2.100 ms  1.260 Mcyc  3.3%
          scan_sweep_column   486.1 us  291.6 kcyc  0.8% x33.4 14.6us/call
        mesh_sweep_paths     109.8 us   65.9 kcyc  0.2%
        scan_face_setup      5.638 ms  3.383 Mcyc  8.9% x419.5 13.4us/call
        mesh_copy_capture    564.2 us  338.5 kcyc  0.9%
        mesh_capture         2.996 ms  1.797 Mcyc  4.7%
        mesh_classify        1.657 ms  994.1 kcyc  2.6%
        scan_mesh_raster    29.882 ms 17.929 Mcyc 47.2% x389.9 76.6us/call
      is_face_offsets         27.1 us   16.3 kcyc  0.0%
      is_mesh_transform      2.667 ms  1.600 Mcyc  4.2%
      mesh_visibility        1.228 ms  737.0 kcyc  1.9%
  is_ripple_prepare           15.5 us    9.3 kcyc  0.0%
  canvas_clear                97.6 us   58.6 kcyc  0.2%
  canvas_buffer_wait         5.784 ms  3.470 Mcyc  9.1%
```

Wall min/avg/max = 61.756/63.328/64.082 ms.
Ordinary face rasterization dominates this stress window. Face construction, capture and endpoint transformation remain separate substantial costs; the peak live frame is 61.337 ms, above this window’s mean.

### Carousel shape draw (frames 513–528)

```text
        scan_face_setup      5.029 ms  3.017 Mcyc  8.1% x385.1 13.1us/call
mesh_copy_capture             48.4 us   29.0 kcyc  0.1%
mesh_capture                 581.8 us  349.1 kcyc  0.9%
scan_face_sweep              1.484 ms  890.7 kcyc  2.4%
  scan_sweep_column          537.0 us  322.2 kcyc  0.9% x40.3 13.3us/call
mesh_sweep_paths              22.3 us   13.4 kcyc  0.0%
mesh_classify                324.7 us  194.8 kcyc  0.5%
frame                       62.293 ms 37.376 Mcyc 100.0%
  pov_preserve_half          143.9 us   86.3 kcyc  0.2%
  is_timeline_step          41.111 ms 24.667 Mcyc 66.0%
    is_draw_shape           40.963 ms 24.578 Mcyc 65.8%
      is_mesh_scan          36.513 ms 21.908 Mcyc 58.6%
        scan_mesh_raster    26.829 ms 16.097 Mcyc 43.1% x375.5 71.4us/call
      is_face_offsets         26.4 us   15.8 kcyc  0.0%
      is_mesh_transform      2.000 ms  1.200 Mcyc  3.2%
      mesh_visibility        1.210 ms  725.9 kcyc  1.9%
  is_ripple_prepare            9.5 us    5.7 kcyc  0.0%
  canvas_clear                92.7 us   55.6 kcyc  0.1%
  canvas_buffer_wait        20.889 ms 12.533 Mcyc 33.5%
```

Wall min/avg/max = 47.736/62.307/72.291 ms.
This window owns the heavy shape but uses the accelerated carousel schedule. Its lower capture work and raster mean do not replace the fixed-shape stress result.

### Carousel construction and transition (frames 1473–1488)

```text
scan_mesh_raster            27.140 ms 16.284 Mcyc 44.3% x271.7 99.9us/call
        scan_face_setup      428.1 us  256.9 kcyc  0.7% x15.9 27.0us/call
mesh_visibility              578.9 us  347.4 kcyc  0.9%
frame                       61.236 ms 36.742 Mcyc 100.0%
  pov_preserve_half          137.8 us   82.7 kcyc  0.2%
  is_timeline_step          41.202 ms 24.721 Mcyc 67.3%
    hk_conway_previous       637.5 us  382.5 kcyc  1.0%
    is_build_draw           38.365 ms 23.019 Mcyc 62.7%
      is_build_scan         36.954 ms 22.173 Mcyc 60.3%
        mesh_copy_capture    226.6 us  136.0 kcyc  0.4%
        mesh_capture         865.6 us  519.3 kcyc  1.4%
        scan_face_sweep      1.590 ms  953.9 kcyc  2.6%
          scan_sweep_column   542.9 us  325.8 kcyc  0.9% x121.9 4.5us/call
        mesh_sweep_paths      35.0 us   21.0 kcyc  0.1%
        mesh_classify        509.4 us  305.6 kcyc  0.8%
      is_mesh_transform       73.2 us   43.9 kcyc  0.1%
    hk_conway_compile        397.8 us  238.7 kcyc  0.6%
    hk_conway_sweep          895.0 us  537.0 kcyc  1.5%
  is_ripple_prepare            2.3 us    1.4 kcyc  0.0%
  canvas_clear                94.4 us   56.7 kcyc  0.2%
  canvas_buffer_wait        19.572 ms 11.743 Mcyc 32.0%
```

Wall min/avg/max = 40.142/61.254/72.686 ms.
Conway construction and its visible scan overlap the segue. This is the carousel’s highest timeline-window mean, while its highest individual render occurs in the shape-draw window above.

### Per-preset table

All 16 shape names occur and the capture returns to its first shape.
V/E/F/I are the final `Built Shape` geometry, not the smaller seed geometry
from `Spawning Shape`. Rows rank by the selected window’s dominant-scope cost.

Exact full-opacity clean holds cannot be certified from face-call counts:
visibility culling changes those counts with pose and ripple. This table uses
post-setup, single-owner, single-mesh drawing windows, excludes construction
windows and advance-straddling windows, and chooses each shape’s highest
timeline mean. The window count describes that selection; opacity transitions
may remain. Peak and spill columns cover all owned frames, transitions included.

| Shape | V/E/F/I | Blended px/f | Timeline ms/f | Render ms/f | Wall fps | Windows | Owned peak ms | Spilled/frames |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | n/a | 41.111 | 41.404 | 16.05 | 6 | 53.742 | 0/264 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | n/a | 35.522 | 35.809 | 15.88 | 4 | 44.712 | 0/216 |
| truncatedCube_hk90_meta | 218/648/432/1296 | n/a | 29.299 | 29.595 | 16.03 | 5 | 53.050 | 0/276 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 1620/2160/542/4320 | n/a | 29.047 | 29.332 | 15.95 | 3 | 41.374 | 0/160 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | n/a | 27.923 | 28.215 | 15.96 | 4 | 51.186 | 0/234 |
| icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | n/a | 27.559 | 27.849 | 16.02 | 3 | 37.616 | 0/144 |
| dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | n/a | 25.682 | 25.974 | 16.05 | 3 | 34.613 | 0/176 |
| rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | n/a | 25.138 | 25.435 | 16.03 | 4 | 36.421 | 0/177 |
| truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | n/a | 25.108 | 25.395 | 16.04 | 6 | 34.634 | 0/186 |
| dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | n/a | 23.655 | 23.950 | 16.07 | 3 | 30.873 | 0/140 |
| icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | n/a | 22.057 | 22.350 | 16.01 | 6 | 36.189 | 0/204 |
| dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | n/a | 22.013 | 22.303 | 16.02 | 5 | 31.212 | 0/208 |
| snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | n/a | 21.173 | 21.461 | 16.05 | 6 | 35.908 | 0/234 |
| octahedron_hk17_ambo_hk73 | 216/288/74/576 | n/a | 20.548 | 20.844 | 16.04 | 6 | 37.372 | 0/210 |
| dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | n/a | 19.935 | 20.223 | 16.03 | 6 | 34.803 | 0/264 |
| icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | n/a | 16.762 | 17.052 | 15.96 | 4 | 34.501 | 0/234 |

### Per-pixel figures

Blended-pixel telemetry is not enabled in these captures. Blended pixels per
frame, cycles per blend and scan cycles per blended pixel are unavailable;
the 10368-pixel quadrant is a geometric reference, not a measured blend count.

## Column-ISR / DMA marshaling cost

Stress peak-window ISR snapshot, rates per rendered frame and per-call
min/average/maximum times:

```text
isr_wake         1170.5/f  0.561/2.008/23.481 us  3.70%
isr_pack          146.2/f  8.713/9.348/11.781 us  2.15%
isr_dma_submit    146.2/f  0.726/0.940/1.016 us  0.21%
```

Carousel peak-window comparison:

```text
isr_wake         1150.7/f  0.563/2.003/27.075 us  3.69%
isr_pack          143.8/f  8.661/9.335/12.055 us  2.15%
isr_dma_submit    143.8/f  0.726/0.940/2.166 us  0.21%
```

- Packing costs 2.15% of its separately timestamped ISR interval;
  submission costs 0.21%. LED wire transfer runs asynchronously;
  transfer duration is not measured by these submit counters.
- Inclusive `isr_wake` costs 3.70%, leaving
  60.188 ms of a 62.5 ms window after that interrupt load.
  Render scopes already include these interruptions; do not subtract them twice.
  Pack and submit are nested inside wake and must not be added to it.
- No measured render spills require a speedup to meet the display budget.
  The independent optimization goal still requires removing 3.337 ms from
  the observed stress peak. DMA-completion and other interrupts are unmeasured.

## Summary ranking

1. `scan_mesh_raster` — 47.2% of the stress window, 29.882 ms/frame.
2. `scan_face_setup` — 8.9% of the stress window, 5.638 ms/frame.
3. `mesh_capture` — 4.7% of the stress window, 2.996 ms/frame.
4. `is_mesh_transform` — 4.2% of the stress window, 2.667 ms/frame.

The preceding shipping capture recorded 67.393 ms and 33/4799 spills.
The retained stress peak is 6.056 ms lower; the changed capture schedule and
source make this a comparison of observed scenarios, not a matched baseline.
Native/WASM rates are not substitutes for these on-device measurements.
Index ranking: **61.337 ms**, zero spills in both scenarios; the carousel’s
16 zero-spill ownership entries total 3327 frames and the stress replay adds
535 separately identified frames.

README cells: peak 🟢 61.337 (16 + stress), spilled 🟢 0/3862 (0.00%).

## Caveats

- All cycle scopes absorb ISR time because CYCCNT free-runs.
- Counter trees retain instrumentation overhead. No per-pixel breakdown image
  was used; adding per-pixel scopes changes the measured workload.
- `filter_blend` is absent. Where enabled, its parent is the first scope that
  enters it, and its subtree can be hidden when that parent has no calls.
- These are shipping selective-O3 images. Face construction, SDF distance and
  scan kernels use their existing hot regions; no global-O3 image or newly
  pinned capture/planning/bounds function is represented.
- Carousel transition speed is 4; stress transition speed is 1. Epoch stretching
  prevents mid-capture teardown. These knobs and selected schedules are explicit;
  a lower carousel peak does not establish the fixed replay’s stress cost.
- Captured firmware identifies source `12ba76700`. Later test-only fixture
  corrections landed with production-identical source on local master
  `965f7d04c`; they are not attributed to the captured image. Build provenance
  identifies the actual captured variants.
- Strict validation fails only root-cycle/wall agreement in each capture.
  Identity, frame monotonicity and complete-window telemetry checks pass;
  carousel shape wrap passes and the stress preset selection passes.

Carousel validator:

```text
effect=IslamicStars  windows=208
  [PASS] captured >=3 windows (208)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] shape markers advance and return to the first (cycle closed): 16 distinct
  [PASS] per-frame render is render, not wall (0 of 208 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 208 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): is_mesh_transform, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, mesh_visibility, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (294.8 ppm, frames 1473-1488)

INVALID: IslamicStars
```

Stress validator:

```text
effect=IslamicStars  windows=67
  [PASS] captured >=3 windows (67)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] fixed profile preset 6/16
  [PASS] per-frame render is render, not wall (0 of 67 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 67 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): mesh_visibility, scan_face_setup, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (280.6 ppm, frames 513-520)

INVALID: IslamicStars
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars` and
`HS_PROFILE_WINDOW=16` (carousel) or 8 (stress). `just profile IslamicStars`
uses the standard wrapper; the explicit commands above reproduce these two
scenario configurations under the COM4 device lock.
