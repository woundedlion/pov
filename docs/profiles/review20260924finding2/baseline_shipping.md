# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-24, **selective -O3**)

Point-in-time snapshot (regenerate with the Harness command below). Baseline with only finding 2 reverted. Runtime excludes setup frame 1.

[Finding 2 baseline/candidate comparison](comparison.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @600 MHz, COM4, segment 0 master; DMA and flywheel ISRs live |
| Image | `profile`; -Os plus selective HS_O3 regions |
| Driver | POVSegmented<288, 4, 480>; IslamicStars 288×144; single-entry playlist |
| Effect | IslamicStars, source `b4a9848ff`; finding 2 baseline |
| Method | HS_PROFILE, 16-frame windows, 210 seconds, epoch 1920 revolutions, TRANS_SPEED=4 |
| Reproduce | `tools/profile_one.sh IslamicStars profile 210 16` with the two defines and COM4 environment shown in Harness |
| Ranges | Runtime frames 2–3328; scopes/ISR windows 17–3328 |
| Captured | 2026-09-24 10:56 PDT |

Raw capture: [raw capture](data/review125_before_ship.log.txt). Replaces the 2026-09-20 snapshot in the standard candidate report locations; the matched baseline is retained with this experiment.

Image size:

```text
FLASH: code:128712, data:194732, headers:8332   free for files:1699840
RAM1: variables:315488, code:44952, padding:20584   free for local variables:143264
RAM2: variables:520064  free for malloc/new:4224
```

Selective optimization in the shipping hot path covers transforms/draw, SDF Face setup/distance, and the face-specialized scan loop.

Artifact provenance (recorded paths describe the original local build; ELF binaries are not published):

```text
source_sha=b4a9848ff44e8a63c7a35b744d34f695a8e37487
compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile_elf_sha256=09b56f35589783f7704a11e68bc5b95be8f15c7d87457cc49f379028266b8688
phantasm_elf_sha256=1e82baf1d15747db9eb55eddb5eec4d17cb2ab8d0b543d67b5fecf172c765947
profile_envdump_sha256=7291fae26ee3409c2a42dd44457d5484d460f684d17403fbe20db946548da634
phantasm_envdump_sha256=c46560e6a6092197538a35f3e4f869a32a9bb340c50391a2578b676b23d9eb63
artifact_profile_elf=build/prof/artifacts/review125_before_ship_b4a9848ff44e_09b56f355897/profile.elf
artifact_phantasm_elf=build/prof/artifacts/review125_before_ship_b4a9848ff44e_09b56f355897/phantasm.elf
artifact_dir=build/prof/artifacts/review125_before_ship_b4a9848ff44e_09b56f355897
```

Both profile and Phantasm ELF SHA256 values were verified against preserved artifacts. The counter window with highest is_mesh_scan cost passes the root cycle/wall cross-check:

frames 2577–2592: 600,782,930 cycles ÷ 600 = 1001304.883 µs vs wall sum 1,001,306 µs, 1.115 ppm. All closed-window frame rows are contiguous and complete.

Source: `b4a9848ff`; log: [raw capture](data/review125_before_ship.log.txt).
Capture SHA256: `7d39a24b65f3b523f6e2842c9b395434df4b4d8555278cbd27820dee3ede08a0`.

Teensy 4.0, 600 MHz, segmented 288×144, four segments, 480 RPM, TS=4.
TS=4 changes build/ripple sampling as well as hold duration; this is a matched TS=4 workload.

## Frame cadence

3327 runtime frames; mean render 23.070 ms; peak 50.859 ms; spilled 0/3327 (0.000%).
Wall min/mean/max: 8.691/62.430/83.628 ms.


Setup frame 1 render: 15.471 ms (excluded). The 62.5 ms budget leaves 11.641 ms at the worst observed runtime frame. All 23 presets have zero spills including transitions; the full ordered cycle wraps at frame 1776. One quadrant is approximately 10,368 pixels. Sync wait fills the remaining time until the next display flip.

Window summaries exclude frames 1?16: worst window mean render is 45.090 ms/frame (frames 2577?2592).

## Phase-by-phase readout

The ordered 23-shape cycle alternates incremental Conway/Hankin construction with completed-mesh ripple draws. The most expensive build window and the completed-mesh hold window are shown below; if the global peak falls in that build window it is identified without repeating the tree. Counter times are ms/us per frame. M/k/c = million/thousand/single cycles per frame; x = calls/frame; final us = µs/call. Percentages use root frame cycles.

### Build: frames 2801–2816

```text
frame                         65.16ms  39.10M 100.0%
  pov_preserve_half          146.96us  88.18k   0.2% x1.0 146.96us
  is_timeline_step            28.55ms  17.13M  43.8%
    is_build_draw             23.42ms  14.05M  35.9%
      is_build_scan           23.39ms  14.03M  35.9% x0.8 31187.17us
      is_mesh_transform       31.39us  18.83k   0.0% x0.8 41.85us
    hk_conway_compile        245.68us 147.41k   0.4% x0.8 327.57us
    hk_conway_sweep          569.05us 341.43k   0.9% x0.8 758.73us
    is_draw_shape              2.65ms   1.59M   4.1%
      is_mesh_scan             2.65ms   1.59M   4.1%
        scan_mesh_raster      22.98ms  13.79M  35.3%
          filter_blend         1.12ms 674.33k   1.7% x16232.8 0.07us
        scan_face_setup        2.91ms   1.75M   4.5% x287.0 10.16us
      is_face_offsets          6.80us   4.08k   0.0% x0.2 27.21us
      is_mesh_transform        2.19us   1.32k   0.0% x0.2 8.77us
  is_ripple_prepare            0.40us    240c   0.0% x1.0 0.40us
  canvas_clear                86.23us  51.74k   0.1% x1.0 86.23us
  canvas_buffer_wait          36.38ms  21.83M  55.8% x1.0 36384.45us
```

Shared/nonexclusive counters: is_mesh_transform (DUPLICATE-NAME), scan_face_setup (MIXED-PARENT), scan_mesh_raster (MIXED-PARENT).

Wall min/mean/max: 50.811/65.165/74.840 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

### Finished/ripple: frames 2577–2592

```text
frame                         62.58ms  37.55M 100.0%
  pov_preserve_half          141.99us  85.20k   0.2% x1.0 141.99us
  is_timeline_step            44.85ms  26.91M  71.7%
    is_draw_shape             44.81ms  26.88M  71.6%
      is_mesh_scan            40.77ms  24.46M  65.2%
        scan_mesh_raster      27.12ms  16.27M  43.3%
          filter_blend         1.32ms 789.12k   2.1% x18214.2 0.07us
        scan_face_setup       13.24ms   7.95M  21.2% x1082.0 12.24us
      is_face_offsets        517.37us 310.42k   0.8% x1.0 517.37us
      is_mesh_transform        3.51ms   2.11M   5.6% x1.0 3514.55us
  is_ripple_prepare            9.07us   5.44k   0.0% x1.0 9.07us
  canvas_clear                84.38us  50.63k   0.1% x1.0 84.38us
  canvas_buffer_wait          17.49ms  10.50M  28.0% x1.0 17492.09us
```

Wall min/mean/max: 60.556/62.582/64.553 ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.

The peak runtime window is the already shown frames 2801–2816; its tree is not repeated.

### Per-preset table

Ranked by worst clean-hold is_mesh_scan window. Window counts use the strict finished-geometry predicate below. Source indices follow the first 23 spawn markers; the cycle wraps back to its first entry. Window fps reflects measured wall time and boundary jitter; zero-spill shapes retain the observed 16 fps cadence.

| # | Shape | V/E/F/I | Windows | Scan ms/frame | Blends/frame | Render ms/frame | Window fps |
|---:|---|---|---:|---:|---:|---:|---:|
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4 | 40.774 | 18214.2 | 45.090 | 15.98 |
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4 | 38.236 | 18755.6 | 40.299 | 15.89 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4 | 36.621 | 18233.1 | 38.687 | 15.98 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2 | 36.470 | 16073.2 | 38.708 | 16.01 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4 | 35.269 | 17499.0 | 37.338 | 15.96 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4 | 31.863 | 17937.0 | 33.244 | 16.00 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4 | 30.066 | 16978.2 | 31.985 | 16.04 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1 | 28.145 | 16530.4 | 30.123 | 16.04 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4 | 26.908 | 16558.1 | 28.434 | 15.99 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4 | 26.389 | 15972.4 | 28.027 | 15.98 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4 | 25.476 | 15325.3 | 26.717 | 15.98 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2 | 25.311 | 18775.1 | 26.403 | 16.08 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2 | 23.392 | 14697.9 | 24.227 | 16.02 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4 | 23.255 | 14633.1 | 24.221 | 16.08 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4 | 21.258 | 14171.7 | 22.143 | 15.99 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4 | 20.370 | 13951.7 | 21.209 | 16.06 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4 | 19.432 | 17006.5 | 20.291 | 16.02 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4 | 18.998 | 13223.2 | 19.480 | 16.03 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4 | 18.246 | 15878.1 | 19.274 | 16.01 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2 | 17.761 | 15855.2 | 18.384 | 16.00 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2 | 17.748 | 13040.0 | 18.309 | 16.02 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2 | 15.762 | 14465.1 | 16.438 | 16.03 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4 | 14.462 | 14700.6 | 15.004 | 16.01 |

No strict clean-hold window: none.

### Per-pixel figures

Filter blend: 41.612189239818335 cycles/call, including instrumentation.


14,895.963 blended pixels/frame, or 1.437× quadrant coverage. Filter cost 41.612 cycles/blend; inclusive raster cost 724.479 cycles/blend. Raster includes build paths outside is_mesh_scan; their costs are not disjoint.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.08/frame  cpu 3.117%
  min/avg/max 0.508/1.690/21.791 us
  isr_pack       144.00/frame  cpu 1.642%
    min/avg/max 6.230/7.124/16.628 us
  isr_dma_submit 144.00/frame  cpu 0.216%
    min/avg/max 0.581/0.935/5.881 us
```


Wake ISR rate is 1152.08 calls/draw frame. Its inclusive CPU share consumes approximately 1.948 ms per 62.5 ms display period, leaving 60.552 ms for foreground work on average. Pack and submit are nested inside wake. DMA submit CPU cost excludes the asynchronous wire transfer. Measured render already includes ISR interruption; no second subtraction is used for headroom.

## Summary ranking

Inclusive averages over the post-setup windows; these scopes overlap:

1. `is_timeline_step`: 22.879 ms/frame (36.64% of wall).
2. `scan_mesh_raster`: 17.986 ms/frame (28.80% of wall).
3. `is_mesh_scan`: 14.516 ms/frame (23.25% of wall).
4. `is_build_draw`: 7.242 ms/frame (11.60% of wall).
5. `filter_blend`: 1.033 ms/frame (1.65% of wall).

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
HS_PROFILE_OUT='build/prof/review125_before_ship.log' \
HS_SESSION=local \
bash tools/profile_one.sh IslamicStars \
  profile 210 16 \
  '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'
```
