# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-24, **selective -O3**)

Point-in-time snapshot (regenerate with the Harness command below). Unlanded finding 2 candidate. Runtime excludes setup frame 1.

The finding 2 baseline/candidate comparison and raw capture archive are no longer retained.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @600 MHz, COM4, segment 0 master; DMA and flywheel ISRs live |
| Image | `profile`; -Os plus selective HS_O3 regions |
| Driver | POVSegmented<288, 4, 480>; IslamicStars 288×144; single-entry playlist |
| Effect | IslamicStars, source `fa92838f9`; unlanded candidate |
| Method | HS_PROFILE, 16-frame windows, 210 seconds, epoch 1920 revolutions, TRANS_SPEED=4 |
| Reproduce | `tools/profile_one.sh IslamicStars profile 210 16` with the two defines and COM4 environment shown in Harness |
| Ranges | Runtime frames 2–3328; scopes/ISR windows 17–3328 |
| Captured | 2026-09-24 11:01 PDT |

Raw capture: `review125_after_ship.log.txt` (archive removed). Replaces the 2026-09-20 snapshot in the standard candidate report locations; the matched baseline archive has been removed.

Image size:

```text
FLASH: code:128856, data:194732, headers:9212   free for files:1698816
RAM1: variables:315488, code:45208, padding:20328   free for local variables:143264
RAM2: variables:520064  free for malloc/new:4224
```

Selective optimization in the shipping hot path covers transforms/draw, SDF Face setup/distance, and the face-specialized scan loop.

Artifact provenance (recorded paths describe the original local build; ELF binaries are not published):

```text
source_sha=fa92838f91063daa49caf606e715f9748674a05b
compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile_elf_sha256=a6cf467a653754981ea0ee30c237bcabf36a7092f6346a9b01f4016118b66699
phantasm_elf_sha256=4c4898a81df08f19dbe2d7cb6263d92a491b943a02834dedd109dd06c6612577
profile_envdump_sha256=016b1af4780daa0876532be1b6f5b14517e0605bc1472035df81eef803f8ca2f
phantasm_envdump_sha256=b481c20f81c8c86c996226ca38a467e3d77e321fe5412f428d2363107c856e1c
artifact_profile_elf=build/prof/artifacts/review125_after_ship_fa92838f9106_a6cf467a6537/profile.elf
artifact_phantasm_elf=build/prof/artifacts/review125_after_ship_fa92838f9106_a6cf467a6537/phantasm.elf
artifact_dir=build/prof/artifacts/review125_after_ship_fa92838f9106_a6cf467a6537
```

Both profile and Phantasm ELF SHA256 values were verified against preserved artifacts. The counter window with highest is_mesh_scan cost passes the root cycle/wall cross-check:

frames 2577–2592: 600,769,606 cycles ÷ 600 = 1001282.677 µs vs wall sum 1,001,283 µs, 0.323 ppm. All closed-window frame rows are contiguous and complete.

Source: `fa92838f9`; log: `review125_after_ship.log.txt` (archive removed).
Capture SHA256: `54a178e36fe3c9999b4f0b8070a58b5124acb7df6e1eef237438a6fd654854a3`.

Teensy 4.0, 600 MHz, segmented 288×144, four segments, 480 RPM, TS=4.
TS=4 changes build/ripple sampling as well as hold duration; this is a matched TS=4 workload.

## Frame cadence

3327 runtime frames; mean render 25.263 ms; peak 60.695 ms; spilled 0/3327 (0.000%).
Wall min/mean/max: 8.947/62.427/98.260 ms.


Setup frame 1 render: 15.820 ms (excluded). The 62.5 ms budget leaves 1.805 ms at the worst observed runtime frame. All 23 presets have zero spills including transitions; the full ordered cycle wraps at frame 1776. One quadrant is approximately 10,368 pixels. Sync wait fills the remaining time until the next display flip.

Window summaries exclude frames 1?16: worst window mean render is 46.069 ms/frame (frames 2577?2592).

## Phase-by-phase readout

The ordered 23-shape cycle alternates incremental Conway/Hankin construction with completed-mesh ripple draws. The most expensive build window and the completed-mesh hold window are shown below; if the global peak falls in that build window it is identified without repeating the tree. Counter times are ms/us per frame. M/k/c = million/thousand/single cycles per frame; x = calls/frame; final us = µs/call. Percentages use root frame cycles.

### Build: frames 2817–2832

```text
scan_mesh_raster              30.23ms  18.14M  48.5%
  filter_blend                 1.48ms 887.46k   2.4% x20161.9 0.07us
scan_face_setup               10.70ms   6.42M  17.2% x1082.0 9.89us
frame                         62.35ms  37.41M 100.0%
  pov_preserve_half          140.72us  84.43k   0.2% x1.0 140.72us
  is_timeline_step            45.51ms  27.31M  73.0%
    is_build_draw             41.42ms  24.85M  66.4%
      is_build_scan           41.32ms  24.79M  66.3% x1.0 41319.34us
      is_mesh_transform      104.39us  62.63k   0.2% x1.0 104.39us
    hk_conway_compile        851.91us 511.14k   1.4% x1.0 851.91us
    hk_conway_sweep            1.73ms   1.04M   2.8% x1.0 1734.63us
  is_ripple_prepare            0.19us    114c   0.0% x1.0 0.19us
  canvas_clear                86.64us  51.99k   0.1% x1.0 86.64us
  canvas_buffer_wait          16.62ms   9.97M  26.6% x1.0 16615.21us
```

Shared/nonexclusive counters: scan_face_setup (MIXED-PARENT), scan_mesh_raster (MIXED-PARENT).

Wall min/mean/max: 53.229/62.352/79.489 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

### Finished/ripple: frames 2577–2592

```text
frame                         62.58ms  37.55M 100.0%
  pov_preserve_half          140.06us  84.03k   0.2% x1.0 140.06us
  is_timeline_step            45.83ms  27.50M  73.2%
    is_draw_shape             45.78ms  27.47M  73.2%
      is_mesh_scan            41.75ms  25.05M  66.7%
        scan_mesh_raster      27.63ms  16.58M  44.1%
          filter_blend         1.33ms 796.85k   2.1% x18215.4 0.07us
        scan_face_setup       13.71ms   8.22M  21.9% x1082.0 12.67us
      is_face_offsets        517.65us 310.59k   0.8% x1.0 517.65us
      is_mesh_transform        3.51ms   2.11M   5.6% x1.0 3512.72us
  is_ripple_prepare            8.80us   5.28k   0.0% x1.0 8.80us
  canvas_clear                84.53us  50.72k   0.1% x1.0 84.53us
  canvas_buffer_wait          16.51ms   9.91M  26.4% x1.0 16511.37us
```

Wall min/mean/max: 60.556/62.580/64.502 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

The peak runtime window is the already shown frames 2817–2832; its tree is not repeated.

### Per-preset table

Ranked by worst clean-hold is_mesh_scan window. Window counts use the strict finished-geometry predicate below. Source indices follow the first 23 spawn markers; the cycle wraps back to its first entry. Window fps reflects measured wall time and boundary jitter; zero-spill shapes retain the observed 16 fps cadence.

| # | Shape | V/E/F/I | Windows | Scan ms/frame | Blends/frame | Render ms/frame | Window fps |
|---:|---|---|---:|---:|---:|---:|---:|
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4 | 41.748 | 18215.4 | 46.069 | 15.98 |
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4 | 40.765 | 18761.8 | 42.830 | 15.89 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4 | 39.202 | 18236.8 | 41.283 | 15.97 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2 | 38.661 | 18940.8 | 39.743 | 16.09 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4 | 37.918 | 17505.1 | 39.988 | 15.96 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2 | 36.483 | 16073.2 | 38.732 | 16.01 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4 | 34.557 | 17942.4 | 35.940 | 15.99 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4 | 31.948 | 16978.7 | 33.871 | 16.04 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1 | 29.545 | 16531.2 | 31.530 | 16.04 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4 | 28.811 | 16564.4 | 30.342 | 15.99 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4 | 27.953 | 15973.8 | 29.594 | 15.98 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4 | 26.673 | 15326.9 | 27.917 | 15.98 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4 | 26.172 | 17119.1 | 26.945 | 15.98 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2 | 24.473 | 14698.1 | 25.315 | 16.02 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4 | 24.410 | 14635.2 | 25.379 | 16.08 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2 | 23.556 | 15908.7 | 24.188 | 16.00 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4 | 22.457 | 15890.1 | 23.494 | 16.01 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4 | 22.061 | 14171.9 | 22.953 | 15.98 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4 | 21.098 | 13951.8 | 21.935 | 16.07 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4 | 19.690 | 13223.2 | 20.178 | 16.02 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2 | 19.091 | 14475.5 | 19.771 | 16.03 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4 | 18.808 | 14744.8 | 19.364 | 16.02 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2 | 18.347 | 13040.2 | 18.916 | 16.02 |

No strict clean-hold window: none.

### Per-pixel figures

Filter blend: 42.2369109434298 cycles/call, including instrumentation.


14,911.995 blended pixels/frame, or 1.438× quadrant coverage. Filter cost 42.237 cycles/blend; inclusive raster cost 800.730 cycles/blend. Raster includes build paths outside is_mesh_scan; their costs are not disjoint.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.08/frame  cpu 3.063%
  min/avg/max 0.485/1.661/27.885 us
  isr_pack       144.00/frame  cpu 1.643%
    min/avg/max 6.230/7.128/21.408 us
  isr_dma_submit 144.00/frame  cpu 0.216%
    min/avg/max 0.583/0.935/10.198 us
```


Wake ISR rate is 1152.08 calls/draw frame. Its inclusive CPU share consumes approximately 1.915 ms per 62.5 ms display period, leaving 60.585 ms for foreground work on average. Pack and submit are nested inside wake. DMA submit CPU cost excludes the asynchronous wire transfer. Measured render already includes ISR interruption; no second subtraction is used for headroom.

## Summary ranking

Inclusive averages over the post-setup windows; these scopes overlap:

1. `is_timeline_step`: 25.081 ms/frame (40.17% of wall).
2. `scan_mesh_raster`: 19.901 ms/frame (31.87% of wall).
3. `is_mesh_scan`: 16.014 ms/frame (25.65% of wall).
4. `is_build_draw`: 7.943 ms/frame (12.72% of wall).
5. `filter_blend`: 1.050 ms/frame (1.68% of wall).

## Caveats

`filter_blend` is registered under `scan_mesh_raster`; latched/shared parent relationships can hide a subtree when its recorded parent has zero calls. This is standard HS_PROFILE instrumentation, including per-blend scope overhead, on both images. No scan-metrics or deep probe-breakdown instrumentation was enabled. These instrumented cycle counts are not an overhead-free shipping timing estimate.

Scope times include interrupts. Wake includes pack and DMA submit; these are not additive. Counter trees preserve duplicate labels and mixed parents; shared scopes do not support exclusive phase attribution. ISR means use rounded window totals and header elapsed time. Asynchronous 600-byte SPI wire transfer is 200 µs at 24 MHz.

Clean holds require one finished geometry, no build draw, one shape draw per frame and exactly F raster calls per frame. Geometry comes from Built Shape, not seed geometry from Spawning Shape.

Trailing unclosed frames excluded: 9. Validation errors: [].

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1920`, `HS_PROFILE_TRANS_SPEED=4`; PlatformIO `profile` and a 210-second capture. The parent capture logs and envdump files preserve the exact build and capture command. One run per configuration does not establish repeat-run uncertainty.

Both source commits were clean when built. Equivalent Git Bash invocation:

```bash
HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4 \
HS_PROFILE_OUT='build/prof/review125_after_ship.log' \
HS_SESSION=review125-profile2 \
bash tools/profile_one.sh IslamicStars \
  profile 210 16 \
  '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'
```
