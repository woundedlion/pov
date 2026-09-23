# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-20, **selective -O3**) - experimental vertical coverage

Point-in-time snapshot of an **unlanded experimental vertical coverage correction**.
Regenerate with the explicit worktree command below; the normal
`just profile IslamicStars` entry point does not select this experimental source.
Raw capture: `build/prof/review_20260920/islamic_24_ship.log`;
captured **2026-09-20 00:31 Pacific Daylight Time**, board **COM3**. This supplements, and does not replace, the
separate baseline report for this configuration. The roster currently contains
**23** shapes, not the 24 in the older skill table. See the paired
[review comparison](../../../review_2026-09-20.md) for the experimental corrections.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile`; `-Os -ffast-math` base, newlib-nano, shipping selective-O3 regions |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, DMA LEDs |
| Effect | IslamicStars 288×144, one-entry playlist; source `09cdc400921e40b63120d423582c0e3186b3edb8` |
| Method | `HS_PROFILE`, 16-frame windows, 150-second capture, Trans Speed 4, epoch 1920 revolutions (240 s); all 23 shapes and wrap observed |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<source-tree> bash tools/profile_one.sh IslamicStars profile 150 16 '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'` |

Image size: `FLASH: code:128848, data:193044, headers:8860   free for files:1700864 / RAM1: variables:315520, code:43832, padding:21704   free for local variables:143232 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 1–16, root
`542455928` cycles ÷ 600 MHz versus measured wall sum
`904095` µs differ by **1.98 ppm**.
`tools/parse_profile.py ... validate` passes marker wrap, monotonic frames,
render/wall separation, and this counter check. Source/compiler/ELF attestations
are appended to the raw log; the source was the clean experimental commit listed above.

## Frame cadence

**Pass aggregate**: render mean **22.843 ms/frame**;
worst window mean **42.980 ms/frame**
(frames 801–816); peak frame render
**50.001 ms** (frames 1025–1040);
spilled **0/2368 frames**
(0.0%). All 23 shape buckets are green.

At 480 RPM each half-revolution lasts 62.5 ms, yielding 16 fps. One board renders
a 144×72 quadrant, about 10,368 pixels. Every captured build, finished-mesh,
ripple and fade phase fits that window; the measured peak leaves
**12.499 ms** of render headroom. `canvas_buffer_wait`
is deliberate idle until the next buffer handoff, not rendering work.

The aggregate covers 148 complete windows. The raw log's final
9 frame lines have no closed counter window and are excluded by the
canonical parser. Both the per-frame lines and window maxima give the same
peak render above. Matched-frame regression comparisons use the `f` lines.

## Phase-by-phase readout

Each entry sweeps in its seed, performs its lowered recipe build, holds and
ripples the finished mesh, then sweeps out. Trans Speed 4 compresses all stages:
the measured native schedule wraps after 1,776 frames (111 seconds at 16 fps).
It also changes build interpolation and ripple sampling; it is accelerated
choreography, not just shorter static holds.

### Recipe-build window (frames 1041–1056)

```text
scan_mesh_raster               23.38 ms  14.03 Mcyc  37.5%
  filter_blend                  1.43 ms 855.14 kcyc   2.3% x19517.6 44cyc
scan_face_setup                 9.59 ms   5.76 Mcyc  15.4% x1082.0 8.9us
frame                          62.38 ms  37.43 Mcyc 100.0%
  pov_preserve_half           147.44 us  88.48 kcyc   0.2% x1.0 147.4us
  is_timeline_step             37.21 ms  22.33 Mcyc  59.7%
    is_build_draw              33.47 ms  20.08 Mcyc  53.7%
      is_build_scan            33.36 ms  20.02 Mcyc  53.5% x1.0 33363.6us
      is_mesh_transform       105.94 us  63.57 kcyc   0.2% x1.0 105.9us
    hk_conway_compile         488.81 us 293.29 kcyc   0.8% x1.0 488.8us
    hk_conway_sweep             1.78 ms   1.07 Mcyc   2.8% x1.0 1775.2us
  is_ripple_prepare             0.19 us    119 cyc   0.0% x1.0 0.2us
  canvas_clear                 84.62 us  50.79 kcyc   0.1% x1.0 84.6us
  canvas_buffer_wait           24.94 ms  14.96 Mcyc  40.0%
```

Wall min/avg/max = 52.310/62.381/78.984 ms. All 16 draws are recipe-build draws. The finished-mesh draw scope is inactive. Shared face counters with that inactive parent are emitted as standalone roots; the build-scan parent still includes their raster work.

### Finished mesh / ripple window (frames 401–416)

```text
frame                          62.95 ms  37.77 Mcyc 100.0%
  pov_preserve_half           146.69 us  88.03 kcyc   0.2% x1.0 146.7us
  is_timeline_step             39.61 ms  23.76 Mcyc  62.9%
    is_draw_shape              39.55 ms  23.73 Mcyc  62.8%
      is_mesh_scan             37.77 ms  22.66 Mcyc  60.0%
        scan_mesh_raster       28.74 ms  17.24 Mcyc  45.6%
          filter_blend          1.37 ms 820.26 kcyc   2.2% x18755.6 44cyc
        scan_face_setup         8.75 ms   5.25 Mcyc  13.9% x722.0 12.1us
      is_face_offsets         348.19 us 208.92 kcyc   0.6% x1.0 348.2us
      is_mesh_transform         1.43 ms 856.51 kcyc   2.3% x1.0 1427.5us
  is_ripple_prepare             5.38 us   3.26 kcyc   0.0% x1.0 5.4us
  canvas_clear                 84.38 us  50.66 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           23.11 ms  13.87 Mcyc  36.7%
```

Wall min/avg/max = 60.375/62.951/67.051 ms. All 16 frames belong to one shape, have zero recipe-build draws, and issue exactly 16 times its finished face count. This is a full finished-mesh window; ripple deformation may remain active.

### Window containing the peak frame (frames 1025–1040)

```text
frame                          65.00 ms  39.00 Mcyc 100.0%
  pov_preserve_half           147.94 us  88.77 kcyc   0.2% x1.0 147.9us
  is_timeline_step             28.22 ms  16.93 Mcyc  43.4%
    is_build_draw              23.18 ms  13.91 Mcyc  35.7%
      is_build_scan            23.15 ms  13.89 Mcyc  35.6% x0.8 30862.1us
      is_mesh_transform        33.62 us  20.21 kcyc   0.1% x0.8 44.8us
    hk_conway_compile         143.00 us  85.81 kcyc   0.2% x0.8 190.7us
    hk_conway_sweep           591.31 us 354.79 kcyc   0.9% x0.8 788.4us
    is_draw_shape               2.66 ms   1.60 Mcyc   4.1%
      is_mesh_scan              2.65 ms   1.59 Mcyc   4.1%
        scan_mesh_raster       22.76 ms  13.66 Mcyc  35.0%
          filter_blend          1.11 ms 666.19 kcyc   1.7% x15992.5 42cyc
        scan_face_setup         2.89 ms   1.74 Mcyc   4.5% x287.0 10.1us
      is_face_offsets           6.81 us   4.11 kcyc   0.0% x0.2 27.2us
      is_mesh_transform         2.12 us   1.30 kcyc   0.0% x0.2 8.5us
  is_ripple_prepare             0.19 us    119 cyc   0.0% x1.0 0.2us
  canvas_clear                 84.56 us  50.76 kcyc   0.1% x1.0 84.6us
  canvas_buffer_wait           36.54 ms  21.92 Mcyc  56.2%
```

Wall min/avg/max = 47.854/64.998/74.400 ms. This window locates the pass peak. Its averaged tree describes the whole 16-frame window and must not be mistaken for the individual peak frame. Build, fade, or ripple work can coexist within it.

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
call-count rows alone do not prove a clean finished hold. Each observed finished-window row here passes
the stronger check: all 16 frame owners match, timeline and finished-mesh draw
each have 16 calls, `is_build_draw` has zero calls,
and `scan_mesh_raster` has exactly 16×finished-F calls. Among those windows,
the row selects the greatest raster time; a finished mesh can still be rippling.
Rows without a qualifying finished window are labelled not observed; their
all-phase peak still comes from actual frame telemetry. Rows are ranked by
observed raster time. `n` is the number of qualifying windows, and
all-phase peak includes the shape's builds/transitions as assigned by the parser.

| # | Shape | Finished V/E/F/I | n; frames | Blended px/f | Raster ms/f | Render ms/f | fps | All-phase peak ms |
|---:|---|---|---|---:|---:|---:|---:|---:|
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4; 401–416 | 18755.6 | 28.736 | 39.843 | 16 | 43.548 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4; 497–512 | 18233.1 | 27.591 | 36.569 | 16 | 44.512 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2; 577–592 | 16127.0 | 27.123 | 36.201 | 16 | 39.156 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 2; 641–656 | 17496.3 | 26.272 | 35.271 | 16 | 38.296 |
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 2; 801–816 | 18154.7 | 25.044 | 42.980 | 16 | 48.966 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 2; 1425–1440 | 17937.0 | 23.731 | 31.997 | 16 | 37.604 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4; 2033–2048 | 16978.2 | 21.917 | 30.720 | 16 | 33.713 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4; 1873–1888 | 16558.1 | 21.858 | 28.174 | 16 | 29.371 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 1; 1073–1088 | 18893.4 | 21.188 | 26.461 | 16 | 50.001 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1; 1745–1760 | 16530.4 | 20.827 | 28.597 | 16 | 30.011 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4; 177–192 | 15972.4 | 20.544 | 26.798 | 16 | 28.889 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 2; 1361–1376 | 15104.6 | 20.387 | 25.361 | 16 | 28.852 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2; 1601–1616 | 14697.9 | 20.055 | 23.233 | 16 | 28.379 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 2; 1505–1520 | 14548.8 | 19.118 | 22.453 | 16 | 25.150 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4; 33–48 | 14171.7 | 17.712 | 20.960 | 16 | 22.629 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 2; 881–896 | 13223.2 | 17.108 | 18.578 | 16 | 20.195 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 2; 721–736 | 13700.8 | 16.405 | 19.616 | 16 | 23.026 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 2; 1233–1248 | 16943.1 | 15.931 | 20.295 | 16 | 41.302 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2; 1665–1680 | 15855.2 | 15.766 | 18.541 | 16 | 39.916 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 1; 1297–1312 | 12949.1 | 15.563 | 17.117 | 16 | 18.486 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4; 2113–2128 | 15878.1 | 15.108 | 19.405 | 16 | 41.068 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 1; 993–1008 | 14536.0 | 13.848 | 16.179 | 16 | 29.027 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 2; 1153–1168 | 14700.6 | 13.147 | 15.129 | 16 | 23.217 |

### Per-pixel figures

The selected heaviest finished-mesh window blends **18755.6 pixels/frame**,
1.809× quadrant coverage. `filter_blend` costs
**43.7 cycles/blend**;
`scan_mesh_raster` costs **919.3 cycles/blended pixel**.
Those scan cycles include probes that were tested and rejected; blended calls
are not a count of all probes or unique pixels. Across the complete capture,
`filter_blend` records 14993.6 blends/frame at
42.1 cycles/blend, and the shared face
raster scope records 705.5 cycles/blend. The annotation-aware reader
includes standalone roots when their latched parent is inactive.

## Column-ISR / DMA marshaling cost

```text
isr_wake            1155.3/f 0.45/1.71/13.90us cpu 3.16%
  isr_pack           143.9/f 6.24/7.16/10.50us cpu 1.64%
  isr_dma_submit     143.9/f 0.58/0.94/6.29us cpu 0.22%
```

The per-call columns are min/weighted-mean/max. CPU percentages use total
captured ISR time divided by the sum of capture-window wall durations.

- Pack plus submit consumes 1.166 ms of CPU per rendered frame; submit alone averages 0.938 µs/call.
- The 72-pixel strobe composite is 600 bytes at 24 MHz: **200 µs** of asynchronous wire time, distinct from the submit CPU cost.
- Pack and submit are nested in `isr_wake`. Its **3.16%** inclusive share must not be added to its children. That share implies approximately **60.53 ms** of foreground opportunity per display window.
- The measured render already includes ISR interruptions. Compare it with the full 62.5 ms window, avoiding a second ISR subtraction: its peak is 0.800× budget, so this captured workload needs no speedup to hold 16 fps.

## Summary ranking

1. `is_timeline_step` — 36.2% of aggregate frame wall time, 22.607 ms/frame.
2. `is_mesh_scan` — 23.0% of aggregate frame wall time, 14.341 ms/frame.
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
- The native TS1 and TS4 full-cycle geometry surveys saw maximum cull radius about 1.514 and zero activations of the 0.01 cosine floor, whose threshold is about 99.995. Ordinary IslamicStars captures cannot measure the extra extreme-geometry work admitted by the horizon interior correction. See the [review comparison](../../../review_2026-09-20.md).
- This is the unlanded variant at the attested commit above, separate from the baseline and the other experimental correction. Before/after conclusions belong in the comparison report.

## Harness

`targets/Profile/Profile.ino` supplies `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1920`, and
`HS_PROFILE_TRANS_SPEED=4`. The device lock, image verification and capture
provenance are provided by `tools/profile_one.sh`; `just profile IslamicStars`
is its normal entry point. The explicit reproduction command above preserves
this cycling capture's nondefault knobs.
