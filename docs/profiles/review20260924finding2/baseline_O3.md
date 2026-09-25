# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-24, **-O3**)

Point-in-time snapshot (regenerate with the Harness command below). Baseline with only finding 2 reverted. Runtime excludes setup frame 1.

[Finding 2 baseline/candidate comparison](comparison.md).

[Shipping selective-O3 sibling](baseline_shipping.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @600 MHz, COM4, segment 0 master; DMA and flywheel ISRs live |
| Image | `profile_o3`; -O3 -ffast-math globally |
| Driver | POVSegmented<288, 4, 480>; IslamicStars 288×144; single-entry playlist |
| Effect | IslamicStars, source `b4a9848ff`; finding 2 baseline |
| Method | HS_PROFILE, 16-frame windows, 210 seconds, epoch 1920 revolutions, TRANS_SPEED=4 |
| Reproduce | `tools/profile_one.sh IslamicStars profile_o3 210 16` with the two defines and COM4 environment shown in Harness |
| Ranges | Runtime frames 2–3328; scopes/ISR windows 17–3328 |
| Captured | 2026-09-24 11:06 PDT |

Raw capture: [raw capture](data/review125_before_o3.log.txt). Replaces the 2026-09-20 snapshot in the standard candidate report locations; the matched baseline is retained with this experiment.

Image size:

```text
FLASH: code:152296, data:195032, headers:9024   free for files:1675264
RAM1: variables:315488, code:57416, padding:8120   free for local variables:143264
RAM2: variables:520064  free for malloc/new:4224
```

Selective optimization in the shipping hot path covers transforms/draw, SDF Face setup/distance, and the face-specialized scan loop.

Artifact provenance (recorded paths describe the original local build; ELF binaries are not published):

```text
source_sha=b4a9848ff44e8a63c7a35b744d34f695a8e37487
compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile_elf_sha256=75fd3d70a67d87ad9c08df43d8771a2f348d48b82c31b6b630ff049e917af76d
phantasm_elf_sha256=9d38ced6c9f92a8176df7e3f1248f507b36e4ce6ddef99a404fcf27e3ca28e8f
profile_envdump_sha256=45bcfada9a09f9d8df275f98de049092fb1275998177ef52234b4dfd3110b915
phantasm_envdump_sha256=1f6ba2d9835a97353750c7bb3ddaa3f3a6a7e41ba7cad8f75f820030a45462f6
artifact_profile_elf=build/prof/artifacts/review125_before_o3_b4a9848ff44e_75fd3d70a67d/profile.elf
artifact_phantasm_elf=build/prof/artifacts/review125_before_o3_b4a9848ff44e_75fd3d70a67d/phantasm.elf
artifact_dir=build/prof/artifacts/review125_before_o3_b4a9848ff44e_75fd3d70a67d
```

Both profile and Phantasm ELF SHA256 values were verified against preserved artifacts. The counter window with highest is_mesh_scan cost passes the root cycle/wall cross-check:

frames 2577–2592: 601,012,893 cycles ÷ 600 = 1001688.155 µs vs wall sum 1,001,689 µs, 0.844 ppm. All closed-window frame rows are contiguous and complete.

Source: `b4a9848ff`; log: [raw capture](data/review125_before_o3.log.txt).
Capture SHA256: `a8feea53c3fd0f460bd94b1a4e9047d576467ef98de071efbcd41b2871fcf82b`.

Teensy 4.0, 600 MHz, segmented 288×144, four segments, 480 RPM, TS=4.
TS=4 changes build/ripple sampling as well as hold duration; this is a matched TS=4 workload.

## Frame cadence

3327 runtime frames; mean render 22.800 ms; peak 50.353 ms; spilled 0/3327 (0.000%).
Wall min/mean/max: 8.725/62.431/79.949 ms.


Setup frame 1 render: 15.590 ms (excluded). The 62.5 ms budget leaves 12.147 ms at the worst observed runtime frame. All 23 presets have zero spills including transitions; the full ordered cycle wraps at frame 1776. One quadrant is approximately 10,368 pixels. Sync wait fills the remaining time until the next display flip.

Window summaries exclude frames 1?16: worst window mean render is 46.547 ms/frame (frames 2577?2592).

## Phase-by-phase readout

The ordered 23-shape cycle alternates incremental Conway/Hankin construction with completed-mesh ripple draws. The most expensive build window and the completed-mesh hold window are shown below; if the global peak falls in that build window it is identified without repeating the tree. Counter times are ms/us per frame. M/k/c = million/thousand/single cycles per frame; x = calls/frame; final us = µs/call. Percentages use root frame cycles.

### Build: frames 2801–2816

```text
frame                         64.93ms  38.96M 100.0%
  pov_preserve_half          143.43us  86.06k   0.2% x1.0 143.43us
  is_timeline_step            28.09ms  16.85M  43.3%
    is_build_draw             23.39ms  14.03M  36.0%
      is_build_scan           23.36ms  14.01M  36.0% x0.8 31143.78us
      is_mesh_transform       31.29us  18.78k   0.0% x0.8 41.72us
    hk_conway_compile        186.25us 111.75k   0.3% x0.8 248.33us
    hk_conway_sweep          467.61us 280.57k   0.7% x0.8 623.48us
    is_draw_shape              2.69ms   1.61M   4.1%
      is_mesh_scan             2.68ms   1.61M   4.1%
        scan_mesh_raster      23.05ms  13.83M  35.5%
          filter_blend         1.25ms 752.67k   1.9% x16367.7 0.08us
        scan_face_setup        2.84ms   1.71M   4.4% x287.0 9.91us
      is_face_offsets          6.83us   4.10k   0.0% x0.2 27.32us
      is_mesh_transform        2.17us   1.30k   0.0% x0.2 8.70us
  is_ripple_prepare            0.85us    512c   0.0% x1.0 0.85us
  canvas_clear                84.42us  50.65k   0.1% x1.0 84.42us
  canvas_buffer_wait          36.61ms  21.97M  56.4% x1.0 36611.09us
```

Shared/nonexclusive counters: is_mesh_transform (DUPLICATE-NAME), scan_face_setup (MIXED-PARENT), scan_mesh_raster (MIXED-PARENT).

Wall min/mean/max: 49.641/64.932/75.629 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

### Finished/ripple: frames 2577–2592

```text
frame                         62.61ms  37.56M 100.0%
  pov_preserve_half          140.21us  84.13k   0.2% x1.0 140.21us
  is_timeline_step            46.31ms  27.79M  74.0%
    is_draw_shape             46.25ms  27.75M  73.9%
      is_mesh_scan            42.25ms  25.35M  67.5%
        scan_mesh_raster      28.93ms  17.36M  46.2%
          filter_blend         1.47ms 883.85k   2.4% x18621.5 0.08us
        scan_face_setup       12.91ms   7.75M  20.6% x1082.0 11.93us
      is_face_offsets        526.66us 315.99k   0.8% x1.0 526.66us
      is_mesh_transform        3.48ms   2.09M   5.6% x1.0 3476.63us
  is_ripple_prepare            8.97us   5.38k   0.0% x1.0 8.97us
  canvas_clear                84.62us  50.77k   0.1% x1.0 84.62us
  canvas_buffer_wait          16.06ms   9.64M  25.7% x1.0 16059.64us
```

Wall min/mean/max: 59.205/62.606/65.886 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

The peak runtime window is the already shown frames 2801–2816; its tree is not repeated.

### Per-preset table

Ranked by worst clean-hold is_mesh_scan window. Window counts use the strict finished-geometry predicate below. Source indices follow the first 23 spawn markers; the cycle wraps back to its first entry. Window fps reflects measured wall time and boundary jitter; zero-spill shapes retain the observed 16 fps cadence.

| # | Shape | V/E/F/I | Windows | Scan ms/frame | Blends/frame | Render ms/frame | Window fps |
|---:|---|---|---:|---:|---:|---:|---:|
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4 | 42.249 | 18621.5 | 46.547 | 15.97 |
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4 | 37.802 | 18755.4 | 39.864 | 15.89 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4 | 36.271 | 18229.4 | 38.342 | 15.98 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2 | 36.198 | 16070.2 | 38.421 | 16.01 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4 | 35.015 | 17492.6 | 37.072 | 15.96 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4 | 31.340 | 17847.2 | 33.448 | 16.10 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4 | 29.836 | 17061.6 | 31.449 | 15.95 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1 | 28.671 | 16389.8 | 30.647 | 16.06 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4 | 26.464 | 16524.4 | 27.994 | 16.01 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4 | 26.219 | 15972.4 | 27.847 | 15.98 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2 | 24.777 | 18723.0 | 25.856 | 16.06 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4 | 24.042 | 15156.8 | 25.280 | 16.00 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2 | 24.025 | 14702.4 | 24.758 | 15.99 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4 | 22.918 | 14587.6 | 23.901 | 16.10 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4 | 21.211 | 14171.7 | 22.091 | 15.99 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4 | 20.687 | 13933.4 | 21.426 | 15.96 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4 | 19.612 | 13333.4 | 20.116 | 16.02 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4 | 19.452 | 17193.6 | 20.224 | 16.01 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2 | 18.436 | 15953.2 | 19.054 | 15.98 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4 | 18.037 | 15808.1 | 18.642 | 16.02 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2 | 17.654 | 13040.4 | 18.220 | 16.01 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2 | 15.156 | 14420.3 | 15.825 | 16.05 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4 | 14.049 | 14401.7 | 14.601 | 16.02 |

No strict clean-hold window: none.

### Per-pixel figures

Filter blend: 46.043324862940125 cycles/call, including instrumentation.


14,878.598 blended pixels/frame, or 1.435× quadrant coverage. Filter cost 46.043 cycles/blend; inclusive raster cost 721.861 cycles/blend. Raster includes build paths outside is_mesh_scan; their costs are not disjoint.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.08/frame  cpu 2.936%
  min/avg/max 0.438/1.592/26.743 us
  isr_pack       144.00/frame  cpu 1.631%
    min/avg/max 5.988/7.076/10.241 us
  isr_dma_submit 144.00/frame  cpu 0.215%
    min/avg/max 0.580/0.935/5.670 us
```


Wake ISR rate is 1152.08 calls/draw frame. Its inclusive CPU share consumes approximately 1.835 ms per 62.5 ms display period, leaving 60.665 ms for foreground work on average. Pack and submit are nested inside wake. DMA submit CPU cost excludes the asynchronous wire transfer. Measured render already includes ISR interruption; no second subtraction is used for headroom.

## Summary ranking

Inclusive averages over the post-setup windows; these scopes overlap:

1. `is_timeline_step`: 22.608 ms/frame (36.20% of wall).
2. `scan_mesh_raster`: 17.900 ms/frame (28.67% of wall).
3. `is_mesh_scan`: 14.385 ms/frame (23.04% of wall).
4. `is_build_draw`: 7.207 ms/frame (11.54% of wall).
5. `filter_blend`: 1.142 ms/frame (1.83% of wall).

## Caveats

`filter_blend` is registered under `scan_mesh_raster`; latched/shared parent relationships can hide a subtree when its recorded parent has zero calls. This is standard HS_PROFILE instrumentation, including per-blend scope overhead, on both images. No scan-metrics or deep probe-breakdown instrumentation was enabled. These instrumented cycle counts are not an overhead-free shipping timing estimate.

Scope times include interrupts. Wake includes pack and DMA submit; these are not additive. Counter trees preserve duplicate labels and mixed parents; shared scopes do not support exclusive phase attribution. ISR means use rounded window totals and header elapsed time. Asynchronous 600-byte SPI wire transfer is 200 µs at 24 MHz.

Clean holds require one finished geometry, no build draw, one shape draw per frame and exactly F raster calls per frame. Geometry comes from Built Shape, not seed geometry from Spawning Shape.

Trailing unclosed frames excluded: 9. Validation errors: [].

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1920`, `HS_PROFILE_TRANS_SPEED=4`; PlatformIO `profile_o3` and a 210-second capture. The parent capture logs and envdump files preserve the exact build and capture command. One run per configuration does not establish repeat-run uncertainty.

Both source commits were clean when built. Equivalent Git Bash invocation:

```bash
HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4 \
HS_PROFILE_OUT='build/prof/review125_before_o3.log' \
HS_SESSION=review125-profile2 \
bash tools/profile_one.sh IslamicStars \
  profile_o3 210 16 \
  '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'
```

## Global -O3 vs selective -O3

Mean render 23.070 → 22.800 ms; peak 50.859 → 50.353 ms. Global -O3 minus shipping: FLASH code +23,584 B; ITCM code +12,464 B. Global O3 is a single-effect reference, not the shipping full roster.
