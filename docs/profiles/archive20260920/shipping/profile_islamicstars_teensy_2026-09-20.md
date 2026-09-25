# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-20, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/review_20260920/islamic_baseline_ship.log`;
captured **2026-09-20 00:16 Pacific Daylight Time**, board **COM3**. This replaces the
2026-08-26 baseline report for this configuration. The roster currently contains
**23** shapes, not the 24 in the older skill table. See the paired
[review comparison](../../review_2026-09-20.md) for the experimental corrections.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile`; `-Os -ffast-math` base, newlib-nano, shipping selective-O3 regions |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, DMA LEDs |
| Effect | IslamicStars 288×144, one-entry playlist; source `d9a39abd8cd69fc6ed91175a6e0398e0d71be961` |
| Method | `HS_PROFILE`, 16-frame windows, 210-second capture, Trans Speed 4, epoch 1920 revolutions (240 s); all 23 shapes and wrap observed |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<source-tree> bash tools/profile_one.sh IslamicStars profile 210 16 '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'` |

Image size: `FLASH: code:128880, data:193044, headers:8828   free for files:1700864 / RAM1: variables:315520, code:43864, padding:21672   free for local variables:143232 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 1–16, root
`542439148` cycles ÷ 600 MHz versus measured wall sum
`904066` µs differ by **0.83 ppm**.
`tools/parse_profile.py ... validate` passes marker wrap, monotonic frames,
render/wall separation, and this counter check. Source/compiler/ELF attestations
are appended to the raw log; the source was clean baseline `d9a39abd`.

## Frame cadence

**Pass aggregate**: render mean **22.437 ms/frame**;
worst window mean **43.424 ms/frame**
(frames 2577–2592); peak frame render
**50.828 ms** (frames 2801–2816);
spilled **0/3328 frames**
(0.0%). All 23 shape buckets are green.

At 480 RPM each half-revolution lasts 62.5 ms, yielding 16 fps. One board renders
a 144×72 quadrant, about 10,368 pixels. Every captured build, finished-mesh,
ripple and fade phase fits that window; the measured peak leaves
**11.672 ms** of render headroom. `canvas_buffer_wait`
is deliberate idle until the next buffer handoff, not rendering work.

The aggregate covers 208 complete windows. The raw log's final
9 frame lines have no closed counter window and are excluded by the
canonical parser. Both the per-frame lines and window maxima give the same
peak render above. Matched-frame regression comparisons use the `f` lines.

## Phase-by-phase readout

Each entry sweeps in its seed, performs its lowered recipe build, holds and
ripples the finished mesh, then sweeps out. Trans Speed 4 compresses all stages:
the measured native schedule wraps after 1,776 frames (111 seconds at 16 fps).
It also changes build interpolation and ripple sampling; it is accelerated
choreography, not just shorter static holds.

### Recipe-build window (frames 2817–2832)

```text
scan_mesh_raster               23.91 ms  14.34 Mcyc  38.4%
  filter_blend                  1.46 ms 874.74 kcyc   2.3% x20032.7 44cyc
scan_face_setup                 9.57 ms   5.74 Mcyc  15.4% x1082.0 8.8us
frame                          62.30 ms  37.38 Mcyc 100.0%
  pov_preserve_half           144.69 us  86.84 kcyc   0.2% x1.0 144.7us
  is_timeline_step             37.70 ms  22.62 Mcyc  60.5%
    is_build_draw              33.98 ms  20.39 Mcyc  54.5%
      is_build_scan            33.87 ms  20.32 Mcyc  54.4% x1.0 33873.8us
      is_mesh_transform       106.50 us  63.94 kcyc   0.2% x1.0 106.5us
    hk_conway_compile         489.12 us 293.48 kcyc   0.8% x1.0 489.1us
    hk_conway_sweep             1.77 ms   1.06 Mcyc   2.8% x1.0 1771.5us
  is_ripple_prepare             0.19 us    119 cyc   0.0% x1.0 0.2us
  canvas_clear                 85.12 us  51.10 kcyc   0.1% x1.0 85.1us
  canvas_buffer_wait           24.36 ms  14.62 Mcyc  39.1%
```

Wall min/avg/max = 52.665/62.297/78.385 ms. All 16 draws are recipe-build draws. The finished-mesh draw scope is inactive. Shared face counters with that inactive parent are emitted as standalone roots; the build-scan parent still includes their raster work.

### Finished mesh / ripple window (frames 401–416)

```text
frame                          62.95 ms  37.77 Mcyc 100.0%
  pov_preserve_half           145.75 us  87.49 kcyc   0.2% x1.0 145.8us
  is_timeline_step             39.58 ms  23.75 Mcyc  62.9%
    is_draw_shape              39.52 ms  23.71 Mcyc  62.8%
      is_mesh_scan             37.75 ms  22.65 Mcyc  60.0%
        scan_mesh_raster       28.73 ms  17.24 Mcyc  45.6%
          filter_blend          1.37 ms 821.29 kcyc   2.2% x18755.6 44cyc
        scan_face_setup         8.73 ms   5.24 Mcyc  13.9% x722.0 12.1us
      is_face_offsets         348.81 us 209.30 kcyc   0.6% x1.0 348.8us
      is_mesh_transform         1.43 ms 856.47 kcyc   2.3% x1.0 1427.4us
  is_ripple_prepare             5.38 us   3.26 kcyc   0.0% x1.0 5.4us
  canvas_clear                 84.50 us  50.74 kcyc   0.1% x1.0 84.5us
  canvas_buffer_wait           23.13 ms  13.88 Mcyc  36.7%
```

Wall min/avg/max = 60.383/62.950/67.041 ms. All 16 frames belong to one shape, have zero recipe-build draws, and issue exactly 16 times its finished face count. This is a full finished-mesh window; ripple deformation may remain active.

### Window containing the peak frame (frames 2801–2816)

```text
frame                          65.05 ms  39.03 Mcyc 100.0%
  pov_preserve_half           147.44 us  88.49 kcyc   0.2% x1.0 147.4us
  is_timeline_step             28.36 ms  17.01 Mcyc  43.6%
    is_build_draw              23.30 ms  13.98 Mcyc  35.8%
      is_build_scan            23.27 ms  13.96 Mcyc  35.8% x0.8 31026.2us
      is_mesh_transform        33.75 us  20.29 kcyc   0.1% x0.8 45.0us
    hk_conway_compile         143.25 us  85.99 kcyc   0.2% x0.8 191.0us
    hk_conway_sweep           590.19 us 354.13 kcyc   0.9% x0.8 786.9us
    is_draw_shape               2.67 ms   1.60 Mcyc   4.1%
      is_mesh_scan              2.66 ms   1.60 Mcyc   4.1%
        scan_mesh_raster       22.90 ms  13.74 Mcyc  35.2%
          filter_blend          1.13 ms 678.94 kcyc   1.7% x16198.3 42cyc
        scan_face_setup         2.89 ms   1.73 Mcyc   4.4% x287.0 10.1us
      is_face_offsets           6.81 us   4.11 kcyc   0.0% x0.2 27.2us
      is_mesh_transform         2.12 us   1.31 kcyc   0.0% x0.2 8.5us
  is_ripple_prepare             0.19 us    119 cyc   0.0% x1.0 0.2us
  canvas_clear                 84.81 us  50.90 kcyc   0.1% x1.0 84.8us
  canvas_buffer_wait           36.46 ms  21.88 Mcyc  56.1%
```

Wall min/avg/max = 49.633/65.053/74.643 ms. This window locates the pass peak. Its averaged tree describes the whole 16-frame window and must not be mistaken for the individual peak frame. Build, fade, or ripple work can coexist within it.

Counter percentages above are relative to the root frame. Nesting follows the
raw counter registry, including its limitations: `scan_mesh_raster` and
`scan_face_setup` are marked `MIXED-PARENT`, and `is_mesh_transform` appears as
`DUPLICATE-NAME` under build and finished-mesh parents. Shared child totals can
exceed the displayed parent; do not sum those branches or interpret their
placement as exclusive phase attribution. Duplicate records are preserved in
these trees rather than overwritten by a label-keyed dictionary. The external
report reader accepts these annotation suffixes; the stock counter-line regex
ignores them. Exact per-frame peak/spill telemetry needs no such adaptation.

### Per-preset table

All 23 shapes are present and the log returns to the initial
`dodecahedron_hk62_ambo_hk62`. Geometry below comes from **Built Shape** records;
**Spawning Shape** records contain the seed geometry. The stock parser's modal
call-count rows alone do not prove a clean finished hold. Each row here passes
the stronger check: all 16 frame owners match, timeline and finished-mesh draw
each have 16 calls, `is_build_draw` has zero calls,
and `scan_mesh_raster` has exactly 16×finished-F calls. Among those windows,
the row selects the greatest raster time; a finished mesh can still be rippling.
Rows are ranked by raster time. `n` is the number of qualifying windows, and
all-phase peak includes the shape's builds/transitions as assigned by the parser.

| # | Shape | Finished V/E/F/I | n; frames | Blended px/f | Raster ms/f | Render ms/f | fps | All-phase peak ms |
|---:|---|---|---|---:|---:|---:|---:|---:|
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4; 401–416 | 18755.6 | 28.726 | 39.819 | 16 | 43.525 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4; 497–512 | 18233.1 | 27.592 | 36.548 | 16 | 44.472 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2; 577–592 | 16127.0 | 27.125 | 36.180 | 16 | 39.133 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4; 641–656 | 17496.3 | 26.288 | 35.252 | 16 | 38.288 |
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4; 2577–2592 | 18214.2 | 25.520 | 43.424 | 16 | 49.109 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4; 1425–1440 | 17937.0 | 23.733 | 31.981 | 16 | 37.576 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4; 2033–2048 | 16978.2 | 21.913 | 30.700 | 16 | 33.704 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4; 1873–1888 | 16558.1 | 21.865 | 28.162 | 16 | 29.361 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2; 2849–2864 | 18775.1 | 21.308 | 26.573 | 16 | 50.828 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1; 1745–1760 | 16530.4 | 20.827 | 28.582 | 16 | 29.996 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4; 177–192 | 15972.4 | 20.549 | 26.786 | 16 | 28.873 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4; 3137–3152 | 15325.3 | 20.459 | 25.395 | 16 | 28.840 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2; 1601–1616 | 14697.9 | 20.059 | 23.226 | 16 | 28.375 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4; 3281–3296 | 14633.1 | 19.850 | 23.171 | 16 | 27.067 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4; 33–48 | 14171.7 | 17.715 | 20.954 | 16 | 22.617 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4; 881–896 | 13223.2 | 17.109 | 18.576 | 16 | 20.186 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4; 2497–2512 | 13951.7 | 16.957 | 20.139 | 16 | 23.877 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4; 3009–3024 | 17006.5 | 16.066 | 20.425 | 16 | 41.464 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2; 1665–1680 | 15855.2 | 15.760 | 18.537 | 16 | 39.939 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2; 3073–3088 | 13040.0 | 15.758 | 17.313 | 16 | 18.480 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4; 2113–2128 | 15878.1 | 15.111 | 19.397 | 16 | 41.066 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2; 2769–2784 | 14465.1 | 14.224 | 16.557 | 16 | 29.003 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4; 1153–1168 | 14700.6 | 13.144 | 15.130 | 16 | 23.365 |

### Per-pixel figures

The selected heaviest finished-mesh window blends **18755.6 pixels/frame**,
1.809× quadrant coverage. `filter_blend` costs
**43.8 cycles/blend**;
`scan_mesh_raster` costs **919.0 cycles/blended pixel**.
Those scan cycles include probes that were tested and rejected; blended calls
are not a count of all probes or unique pixels. Across the complete capture,
`filter_blend` records 14884.7 blends/frame at
42.1 cycles/blend, and the shared face
raster scope records 701.2 cycles/blend. The annotation-aware reader
includes standalone roots when their latched parent is inactive.

## Column-ISR / DMA marshaling cost

```text
isr_wake            1154.2/f 0.44/1.71/23.00us cpu 3.16%
  isr_pack           143.9/f 6.24/7.16/10.39us cpu 1.65%
  isr_dma_submit     143.9/f 0.59/0.94/7.24us cpu 0.22%
```

The per-call columns are min/weighted-mean/max. CPU percentages use total
captured ISR time divided by the sum of capture-window wall durations.

- Pack plus submit consumes 1.165 ms of CPU per rendered frame; submit alone averages 0.938 µs/call.
- The 72-pixel strobe composite is 600 bytes at 24 MHz: **200 µs** of asynchronous wire time, distinct from the submit CPU cost.
- Pack and submit are nested in `isr_wake`. Its **3.16%** inclusive share must not be added to its children. That share implies approximately **60.53 ms** of foreground opportunity per display window.
- The measured render already includes ISR interruptions. Compare it with the full 62.5 ms window, avoiding a second ISR subtraction: its peak is 0.813× budget, so this captured workload needs no speedup to hold 16 fps.

## Summary ranking

1. `is_timeline_step` — 35.6% of aggregate frame wall time, 22.202 ms/frame.
2. `is_mesh_scan` — 22.4% of aggregate frame wall time, 13.967 ms/frame.
3. `is_build_scan` — 11.5% of aggregate frame wall time, 7.149 ms/frame.

The timeline parent includes both scan children and scheduling/mesh work; these
ranking entries overlap. No matched WASM/native timing capture is used. The
separate native applicability run measures geometry coverage, not device speed.

## Caveats

- All cycle scopes include ISR interruptions because CYCCNT free-runs.
- Shared face counters retain their first registered parent, and the duplicate transform label is ambiguous in the stock parser dictionary. Trees preserve raw records; exclusive per-phase totals cannot be reconstructed from those shared labels.
- `filter_blend` is a per-pixel scope and adds measurement overhead. Both configurations use the same instrumented harness; a counter whose latched parent is inactive is emitted as a standalone root.
- Shipping selective-O3 regions cover IslamicStars mesh transforms/draws, SDF Face setup/distance, and the face-specialized scan loop. The global-O3 image changes all eligible code generation and is a reference, not the full-roster shipping build.
- Trans Speed 4, window 16 and epoch 1920 are capture knobs. The current build/ripple choreography means TS4 changes temporal sampling as well as duration. Do not treat it as the default TS1 frame distribution.
- The native TS1 and TS4 full-cycle geometry surveys saw maximum cull radius about 1.514 and zero activations of the 0.01 cosine floor, whose threshold is about 99.995. Ordinary IslamicStars captures cannot measure the extra extreme-geometry work admitted by finding 126. See the [review comparison](../../review_2026-09-20.md).
- Baseline source is the clean attested commit above. This report does not include the candidate 126 or 24 corrections; before/after conclusions belong in the comparison report.

## Harness

`targets/Profile/Profile.ino` supplies `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1920`, and
`HS_PROFILE_TRANS_SPEED=4`. The device lock, image verification and capture
provenance are provided by `tools/profile_one.sh`; `just profile IslamicStars`
is its normal entry point. The explicit reproduction command above preserves
this cycling capture's nondefault knobs.
