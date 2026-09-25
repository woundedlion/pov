# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-20, **-O3**) - experimental vertical coverage

Point-in-time snapshot of an **unlanded experimental vertical coverage correction**.
Regenerate with the explicit worktree command below; the normal
`just profile IslamicStars` entry point does not select this experimental source.
Raw capture: `build/prof/review_20260920/islamic_24_o3.log`;
captured **2026-09-20 00:35 Pacific Daylight Time**, board **COM3**. This supplements, and does not replace, the
separate baseline report for this configuration. The roster currently contains
**23** shapes, not the 24 in the older skill table. See the paired
[review comparison](../../../review_2026-09-20.md) for the experimental corrections.

Shipping sibling: [selective-O3 report](../shipping/profile_islamicstars_teensy_2026-09-20.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile_o3`; global `-O3 -ffast-math`, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, DMA LEDs |
| Effect | IslamicStars 288×144, one-entry playlist; source `09cdc400921e40b63120d423582c0e3186b3edb8` |
| Method | `HS_PROFILE`, 16-frame windows, 150-second capture, Trans Speed 4, epoch 1920 revolutions (240 s); all 23 shapes and wrap observed |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<source-tree> bash tools/profile_one.sh IslamicStars profile_o3 150 16 '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'` |

Image size: `FLASH: code:153480, data:193280, headers:8568   free for files:1676288 / RAM1: variables:315552, code:57000, padding:8536   free for local variables:143200 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 1–16, root
`542253242` cycles ÷ 600 MHz versus measured wall sum
`903760` µs differ by **5.09 ppm**.
`tools/parse_profile.py ... validate` passes marker wrap, monotonic frames,
render/wall separation, and this counter check. Source/compiler/ELF attestations
are appended to the raw log; the source was the clean experimental commit listed above.

## Frame cadence

**Pass aggregate**: render mean **22.497 ms/frame**;
worst window mean **42.373 ms/frame**
(frames 801–816); peak frame render
**48.866 ms** (frames 769–784);
spilled **0/2368 frames**
(0.0%). All 23 shape buckets are green.

At 480 RPM each half-revolution lasts 62.5 ms, yielding 16 fps. One board renders
a 144×72 quadrant, about 10,368 pixels. Every captured build, finished-mesh,
ripple and fade phase fits that window; the measured peak leaves
**13.634 ms** of render headroom. `canvas_buffer_wait`
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
scan_mesh_raster               21.70 ms  13.02 Mcyc  34.7%
  filter_blend                  1.46 ms 877.25 kcyc   2.3% x18659.9 47cyc
scan_face_setup                 9.41 ms   5.64 Mcyc  15.1% x1082.0 8.7us
frame                          62.49 ms  37.50 Mcyc 100.0%
  pov_preserve_half           145.44 us  87.28 kcyc   0.2% x1.0 145.4us
  is_timeline_step             34.90 ms  20.94 Mcyc  55.8%
    is_build_draw              31.59 ms  18.96 Mcyc  50.6%
      is_build_scan            31.49 ms  18.89 Mcyc  50.4% x1.0 31486.2us
      is_mesh_transform       106.31 us  63.81 kcyc   0.2% x1.0 106.3us
    hk_conway_compile         634.12 us 380.48 kcyc   1.0% x1.0 634.1us
    hk_conway_sweep             1.44 ms 861.66 kcyc   2.3% x1.0 1436.1us
  is_ripple_prepare             0.12 us    106 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.12 us  50.49 kcyc   0.1% x1.0 84.1us
  canvas_buffer_wait           27.36 ms  16.42 Mcyc  43.8%
```

Wall min/avg/max = 52.227/62.494/76.846 ms. All 16 draws are recipe-build draws. The finished-mesh draw scope is inactive. Shared face counters with that inactive parent are emitted as standalone roots; the build-scan parent still includes their raster work.

### Finished mesh / ripple window (frames 401–416)

```text
frame                          62.93 ms  37.76 Mcyc 100.0%
  pov_preserve_half           141.50 us  84.93 kcyc   0.2% x1.0 141.5us
  is_timeline_step             39.05 ms  23.43 Mcyc  62.1%
    is_draw_shape              38.99 ms  23.39 Mcyc  62.0%
      is_mesh_scan             37.22 ms  22.33 Mcyc  59.2%
        scan_mesh_raster       28.35 ms  17.01 Mcyc  45.1%
          filter_blend          1.50 ms 897.80 kcyc   2.4% x18755.4 48cyc
        scan_face_setup         8.59 ms   5.16 Mcyc  13.7% x722.0 11.9us
      is_face_offsets         354.19 us 212.52 kcyc   0.6% x1.0 354.2us
      is_mesh_transform         1.41 ms 846.52 kcyc   2.2% x1.0 1410.8us
  is_ripple_prepare             6.12 us   3.71 kcyc   0.0% x1.0 6.1us
  canvas_clear                 84.38 us  50.64 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           23.65 ms  14.19 Mcyc  37.6%
```

Wall min/avg/max = 60.455/62.930/66.971 ms. All 16 frames belong to one shape, have zero recipe-build draws, and issue exactly 16 times its finished face count. This is a full finished-mesh window; ripple deformation may remain active.

### Window containing the peak frame (frames 769–784)

```text
scan_mesh_raster               17.22 ms  10.33 Mcyc  27.0%
  filter_blend                  1.10 ms 659.40 kcyc   1.7% x14452.5 46cyc
scan_face_setup                 5.04 ms   3.03 Mcyc   7.9% x407.0 12.4us
frame                          63.69 ms  38.21 Mcyc 100.0%
  pov_preserve_half           143.44 us  86.09 kcyc   0.2% x1.0 143.4us
  is_timeline_step             26.11 ms  15.67 Mcyc  41.0%
    is_build_draw              22.53 ms  13.52 Mcyc  35.4%
      is_build_scan            22.45 ms  13.47 Mcyc  35.2% x1.0 22450.0us
      is_mesh_transform        76.06 us  45.65 kcyc   0.1% x1.0 76.1us
    hk_conway_compile         353.00 us 211.80 kcyc   0.6% x1.0 353.0us
    hk_conway_sweep           809.81 us 485.91 kcyc   1.3% x1.0 809.8us
  is_ripple_prepare             0.12 us    107 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.44 us  50.69 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           37.35 ms  22.41 Mcyc  58.6%
```

Wall min/avg/max = 48.971/63.688/78.757 ms. This window locates the pass peak. Its averaged tree describes the whole 16-frame window and must not be mistaken for the individual peak frame. Build, fade, or ripple work can coexist within it.

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
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4; 401–416 | 18755.4 | 28.351 | 39.285 | 16 | 42.948 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4; 497–512 | 18229.4 | 27.244 | 36.084 | 16 | 43.609 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2; 577–592 | 16123.9 | 26.985 | 35.939 | 16 | 37.998 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 2; 641–656 | 17489.6 | 26.043 | 34.894 | 16 | 37.794 |
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 2; 801–816 | 18264.3 | 24.751 | 42.373 | 16 | 48.866 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 2; 1441–1456 | 17847.2 | 23.068 | 31.945 | 16 | 36.616 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4; 2017–2032 | 17061.5 | 21.743 | 30.124 | 16 | 33.287 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4; 97–112 | 16524.4 | 21.479 | 27.702 | 16 | 28.509 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1; 1745–1760 | 16389.8 | 21.272 | 28.940 | 16 | 32.534 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2; 1585–1600 | 14702.4 | 20.649 | 23.702 | 16 | 28.427 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 1; 1073–1088 | 18660.9 | 20.551 | 25.723 | 16 | 48.204 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4; 177–192 | 15972.4 | 20.341 | 26.516 | 16 | 28.607 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 2; 1361–1376 | 15102.8 | 18.798 | 23.690 | 16 | 26.552 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 2; 1505–1520 | 14565.8 | 17.985 | 21.257 | 16 | 22.120 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4; 33–48 | 14171.7 | 17.623 | 20.833 | 16 | 23.679 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 2; 865–880 | 13294.2 | 17.399 | 18.884 | 16 | 19.853 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2; 1665–1680 | 15953.2 | 16.265 | 18.968 | 16 | 38.656 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 2; 721–736 | 13632.2 | 16.237 | 19.408 | 16 | 22.688 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 2; 1217–1232 | 17193.6 | 15.935 | 20.127 | 16 | 37.381 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 1; 1297–1312 | 12923.9 | 15.480 | 17.021 | 16 | 18.186 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4; 353–368 | 15808.1 | 14.765 | 18.562 | 16 | 37.880 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 1; 993–1008 | 14420.3 | 13.462 | 15.759 | 16 | 28.228 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 2; 1153–1168 | 14388.9 | 12.394 | 14.363 | 16 | 22.697 |

### Per-pixel figures

The selected heaviest finished-mesh window blends **18755.4 pixels/frame**,
1.809× quadrant coverage. `filter_blend` costs
**47.9 cycles/blend**;
`scan_mesh_raster` costs **907.0 cycles/blended pixel**.
Those scan cycles include probes that were tested and rejected; blended calls
are not a count of all probes or unique pixels. Across the complete capture,
`filter_blend` records 14970.1 blends/frame at
46.1 cycles/blend, and the shared face
raster scope records 697.6 cycles/blend. The annotation-aware reader
includes standalone roots when their latched parent is inactive.

## Column-ISR / DMA marshaling cost

```text
isr_wake            1155.3/f 0.34/1.53/27.11us cpu 2.82%
  isr_pack           143.9/f 0.53/7.04/10.64us cpu 1.62%
  isr_dma_submit     143.9/f 0.58/0.93/7.33us cpu 0.21%
```

The per-call columns are min/weighted-mean/max. CPU percentages use total
captured ISR time divided by the sum of capture-window wall durations.

- Pack plus submit consumes 1.147 ms of CPU per rendered frame; submit alone averages 0.929 µs/call.
- The 72-pixel strobe composite is 600 bytes at 24 MHz: **200 µs** of asynchronous wire time, distinct from the submit CPU cost.
- Pack and submit are nested in `isr_wake`. Its **2.82%** inclusive share must not be added to its children. That share implies approximately **60.74 ms** of foreground opportunity per display window.
- The measured render already includes ISR interruptions. Compare it with the full 62.5 ms window, avoiding a second ISR subtraction: its peak is 0.782× budget, so this captured workload needs no speedup to hold 16 fps.

## Summary ranking

1. `is_timeline_step` — 35.7% of aggregate frame wall time, 22.263 ms/frame.
2. `is_mesh_scan` — 22.6% of aggregate frame wall time, 14.125 ms/frame.
3. `is_build_scan` — 11.3% of aggregate frame wall time, 7.069 ms/frame.

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

### Global -O3 vs selective -O3

Global-O3 peak render is 48.866 ms versus shipping
50.001 ms; both report zero spilled frames. Those are
pass maxima, not a matched-phase speedup ratio. The global-O3 image adds
**+24,632 B FLASH code** and **+13,168 B ITCM code**.
The exact per-frame matched comparison is recorded in the review report.
