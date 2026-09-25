# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-24, **-O3**)

Point-in-time snapshot (regenerate with the Harness command below). Unlanded finding 2 candidate. Runtime excludes setup frame 1.

[Finding 2 baseline/candidate comparison](../review20260924finding2/comparison.md).

[Shipping selective-O3 sibling](../shipping/profile_islamicstars_teensy_2026-09-24.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @600 MHz, COM4, segment 0 master; DMA and flywheel ISRs live |
| Image | `profile_o3`; -O3 -ffast-math globally |
| Driver | POVSegmented<288, 4, 480>; IslamicStars 288×144; single-entry playlist |
| Effect | IslamicStars, source `fa92838f9`; unlanded candidate |
| Method | HS_PROFILE, 16-frame windows, 210 seconds, epoch 1920 revolutions, TRANS_SPEED=4 |
| Reproduce | `tools/profile_one.sh IslamicStars profile_o3 210 16` with the two defines and COM4 environment shown in Harness |
| Ranges | Runtime frames 2–3328; scopes/ISR windows 17–3328 |
| Captured | 2026-09-24 11:12 PDT |

Raw capture: [raw capture](../review20260924finding2/data/review125_after_o3.log.txt). Replaces the 2026-09-20 snapshot in the standard candidate report locations; the matched baseline is retained with this experiment.

Image size:

```text
FLASH: code:152328, data:195032, headers:8992   free for files:1675264
RAM1: variables:315488, code:57448, padding:8088   free for local variables:143264
RAM2: variables:520064  free for malloc/new:4224
```

Selective optimization in the shipping hot path covers transforms/draw, SDF Face setup/distance, and the face-specialized scan loop.

Artifact provenance (recorded paths describe the original local build; ELF binaries are not published):

```text
source_sha=fa92838f91063daa49caf606e715f9748674a05b
compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile_elf_sha256=c31f26d6587fc5b1eb72118d2c2d179766a02571e7f23070a2b3857393eb282f
phantasm_elf_sha256=236a0c13712bfe293ea20dde67c53aef766f9cf66232958159d1533b19e9acf1
profile_envdump_sha256=1959e355b7f913c17641588f94690d18634432329e5eb2586bd8c880f54d813b
phantasm_envdump_sha256=e719295ac909cb49c55e27eb85c9ea4447491670ff889e4adae0b754f122d840
artifact_profile_elf=build/prof/artifacts/review125_after_o3_fa92838f9106_c31f26d6587f/profile.elf
artifact_phantasm_elf=build/prof/artifacts/review125_after_o3_fa92838f9106_c31f26d6587f/phantasm.elf
artifact_dir=build/prof/artifacts/review125_after_o3_fa92838f9106_c31f26d6587f
```

Both profile and Phantasm ELF SHA256 values were verified against preserved artifacts. The counter window with highest is_mesh_scan cost passes the root cycle/wall cross-check:

frames 2577–2592: 600,921,079 cycles ÷ 600 = 1001535.132 µs vs wall sum 1,001,535 µs, 0.131 ppm. All closed-window frame rows are contiguous and complete.

Source: `fa92838f9`; log: [raw capture](../review20260924finding2/data/review125_after_o3.log.txt).
Capture SHA256: `b18c5dccffeb61f10e1970198183216915206ef87bcb0aa707f4dda295c1545c`.

Teensy 4.0, 600 MHz, segmented 288×144, four segments, 480 RPM, TS=4.
TS=4 changes build/ripple sampling as well as hold duration; this is a matched TS=4 workload.

## Frame cadence

3327 runtime frames; mean render 24.816 ms; peak 55.674 ms; spilled 0/3327 (0.000%).
Wall min/mean/max: 8.951/62.431/96.107 ms.


Setup frame 1 render: 15.938 ms (excluded). The 62.5 ms budget leaves 6.826 ms at the worst observed runtime frame. All 23 presets have zero spills including transitions; the full ordered cycle wraps at frame 1776. One quadrant is approximately 10,368 pixels. Sync wait fills the remaining time until the next display flip.

Window summaries exclude frames 1?16: worst window mean render is 47.634 ms/frame (frames 2577?2592).

## Phase-by-phase readout

The ordered 23-shape cycle alternates incremental Conway/Hankin construction with completed-mesh ripple draws. The most expensive build window and the completed-mesh hold window are shown below; if the global peak falls in that build window it is identified without repeating the tree. Counter times are ms/us per frame. M/k/c = million/thousand/single cycles per frame; x = calls/frame; final us = µs/call. Percentages use root frame cycles.

### Build: frames 1041–1056

```text
scan_mesh_raster              27.56ms  16.53M  44.0%
  filter_blend                 1.45ms 872.44k   2.3% x18771.2 0.08us
scan_face_setup               10.49ms   6.29M  16.8% x1082.0 9.69us
frame                         62.56ms  37.54M 100.0%
  pov_preserve_half          140.95us  84.57k   0.2% x1.0 140.95us
  is_timeline_step            41.79ms  25.08M  66.8%
    is_build_draw             38.51ms  23.11M  61.6%
      is_build_scan           38.41ms  23.04M  61.4% x1.0 38407.30us
      is_mesh_transform      103.67us  62.20k   0.2% x1.0 103.67us
    hk_conway_compile        646.58us 387.95k   1.0% x1.0 646.58us
    hk_conway_sweep            1.41ms 844.31k   2.2% x1.0 1407.18us
  is_ripple_prepare            1.43us    859c   0.0% x1.0 1.43us
  canvas_clear                84.88us  50.93k   0.1% x1.0 84.88us
  canvas_buffer_wait          20.54ms  12.32M  32.8% x1.0 20540.89us
```

Shared/nonexclusive counters: scan_face_setup (MIXED-PARENT), scan_mesh_raster (MIXED-PARENT).

Wall min/mean/max: 52.815/62.563/77.660 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

### Finished/ripple: frames 2577–2592

```text
frame                         62.60ms  37.56M 100.0%
  pov_preserve_half          142.99us  85.79k   0.2% x1.0 142.99us
  is_timeline_step            47.39ms  28.44M  75.7%
    is_draw_shape             47.34ms  28.40M  75.6%
      is_mesh_scan            43.34ms  26.00M  69.2%
        scan_mesh_raster      29.58ms  17.75M  47.3%
          filter_blend         1.45ms 869.57k   2.3% x18626.0 0.08us
        scan_face_setup       13.38ms   8.03M  21.4% x1082.0 12.36us
      is_face_offsets        524.84us 314.90k   0.8% x1.0 524.84us
      is_mesh_transform        3.48ms   2.09M   5.6% x1.0 3476.09us
  is_ripple_prepare            8.93us   5.36k   0.0% x1.0 8.93us
  canvas_clear                84.60us  50.76k   0.1% x1.0 84.60us
  canvas_buffer_wait          14.96ms   8.98M  23.9% x1.0 14962.83us
```

Wall min/mean/max: 59.004/62.596/66.011 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

The peak runtime window is the already shown frames 1041–1056; its tree is not repeated.

### Per-preset table

Ranked by worst clean-hold is_mesh_scan window. Window counts use the strict finished-geometry predicate below. Source indices follow the first 23 spawn markers; the cycle wraps back to its first entry. Window fps reflects measured wall time and boundary jitter; zero-spill shapes retain the observed 16 fps cadence.

| # | Shape | V/E/F/I | Windows | Scan ms/frame | Blends/frame | Render ms/frame | Window fps |
|---:|---|---|---:|---:|---:|---:|---:|
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4 | 43.339 | 18626.0 | 47.634 | 15.98 |
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4 | 40.143 | 18761.5 | 42.195 | 15.88 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4 | 38.672 | 18233.2 | 40.743 | 15.97 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4 | 37.491 | 17499.1 | 39.543 | 15.95 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2 | 37.271 | 18892.2 | 38.344 | 16.04 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2 | 36.221 | 16070.2 | 38.441 | 16.01 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4 | 33.806 | 17851.1 | 35.906 | 16.10 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4 | 31.583 | 17063.8 | 33.192 | 15.95 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1 | 30.064 | 16390.9 | 32.032 | 16.07 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4 | 28.211 | 16526.7 | 29.737 | 16.01 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4 | 27.642 | 15973.8 | 29.264 | 15.98 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4 | 25.655 | 17277.8 | 26.414 | 16.00 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4 | 24.977 | 15157.1 | 26.211 | 16.00 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2 | 24.961 | 14702.6 | 25.694 | 15.98 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4 | 23.930 | 14588.9 | 24.905 | 16.10 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2 | 23.699 | 16005.6 | 24.308 | 15.99 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4 | 21.924 | 14171.9 | 22.798 | 15.98 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4 | 21.824 | 15762.1 | 22.839 | 16.00 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4 | 21.354 | 13934.8 | 22.082 | 15.96 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4 | 20.112 | 13333.5 | 20.608 | 16.02 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2 | 18.195 | 14430.6 | 18.859 | 16.04 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2 | 18.144 | 13040.4 | 18.707 | 16.01 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4 | 18.025 | 14434.9 | 18.570 | 16.02 |

No strict clean-hold window: none.

### Per-pixel figures

Filter blend: 45.40118841688525 cycles/call, including instrumentation.


14,894.526 blended pixels/frame, or 1.437× quadrant coverage. Filter cost 45.401 cycles/blend; inclusive raster cost 791.646 cycles/blend. Raster includes build paths outside is_mesh_scan; their costs are not disjoint.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.08/frame  cpu 2.834%
  min/avg/max 0.338/1.537/21.920 us
  isr_pack       144.00/frame  cpu 1.625%
    min/avg/max 5.988/7.051/19.130 us
  isr_dma_submit 144.00/frame  cpu 0.218%
    min/avg/max 0.588/0.945/5.893 us
```


Wake ISR rate is 1152.08 calls/draw frame. Its inclusive CPU share consumes approximately 1.771 ms per 62.5 ms display period, leaving 60.729 ms for foreground work on average. Pack and submit are nested inside wake. DMA submit CPU cost excludes the asynchronous wire transfer. Measured render already includes ISR interruption; no second subtraction is used for headroom.

## Summary ranking

Inclusive averages over the post-setup windows; these scopes overlap:

1. `is_timeline_step`: 24.632 ms/frame (39.45% of wall).
2. `scan_mesh_raster`: 19.652 ms/frame (31.47% of wall).
3. `is_mesh_scan`: 15.764 ms/frame (25.24% of wall).
4. `is_build_draw`: 7.851 ms/frame (12.57% of wall).
5. `filter_blend`: 1.127 ms/frame (1.80% of wall).

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
HS_PROFILE_OUT='build/prof/review125_after_o3.log' \
HS_SESSION=review125-profile2 \
bash tools/profile_one.sh IslamicStars \
  profile_o3 210 16 \
  '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'
```

## Global -O3 vs selective -O3

Mean render 25.263 → 24.816 ms; peak 60.695 → 55.674 ms. Global -O3 minus shipping: FLASH code +23,472 B; ITCM code +12,240 B. Global O3 is a single-effect reference, not the shipping full roster.
