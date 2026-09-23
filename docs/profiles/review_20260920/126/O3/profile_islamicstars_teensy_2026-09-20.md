# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-20, **-O3**) - experimental horizon interiors

Point-in-time snapshot of an **unlanded experimental horizon interior correction**.
Regenerate with the explicit worktree command below; the normal
`just profile IslamicStars` entry point does not select this experimental source.
Raw capture: `build/prof/review_20260920/islamic_126_o3.log`;
captured **2026-09-20 00:28 Pacific Daylight Time**, board **COM3**. This supplements, and does not replace, the
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
| Effect | IslamicStars 288×144, one-entry playlist; source `846c2e746d08aca7ef3d6b38fc07c9edc0d73726` |
| Method | `HS_PROFILE`, 16-frame windows, 150-second capture, Trans Speed 4, epoch 1920 revolutions (240 s); all 23 shapes and wrap observed |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<source-tree> bash tools/profile_one.sh IslamicStars profile_o3 150 16 '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'` |

Image size: `FLASH: code:153480, data:193280, headers:8568   free for files:1676288 / RAM1: variables:315552, code:57000, padding:8536   free for local variables:143200 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 1–16, root
`542206818` cycles ÷ 600 MHz versus measured wall sum
`903682` µs differ by **4.39 ppm**.
`tools/parse_profile.py ... validate` passes marker wrap, monotonic frames,
render/wall separation, and this counter check. Source/compiler/ELF attestations
are appended to the raw log; the source was the clean experimental commit listed above.

## Frame cadence

**Pass aggregate**: render mean **22.608 ms/frame**;
worst window mean **42.601 ms/frame**
(frames 801–816); peak frame render
**49.074 ms** (frames 769–784);
spilled **0/2368 frames**
(0.0%). All 23 shape buckets are green.

At 480 RPM each half-revolution lasts 62.5 ms, yielding 16 fps. One board renders
a 144×72 quadrant, about 10,368 pixels. Every captured build, finished-mesh,
ripple and fade phase fits that window; the measured peak leaves
**13.426 ms** of render headroom. `canvas_buffer_wait`
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
scan_mesh_raster               21.82 ms  13.09 Mcyc  34.9%
  filter_blend                  1.44 ms 864.23 kcyc   2.3% x18659.9 46cyc
scan_face_setup                 9.46 ms   5.68 Mcyc  15.1% x1082.0 8.7us
frame                          62.50 ms  37.50 Mcyc 100.0%
  pov_preserve_half           142.56 us  85.56 kcyc   0.2% x1.0 142.6us
  is_timeline_step             35.07 ms  21.04 Mcyc  56.1%
    is_build_draw              31.77 ms  19.06 Mcyc  50.8%
      is_build_scan            31.66 ms  18.99 Mcyc  50.7% x1.0 31658.0us
      is_mesh_transform       106.25 us  63.78 kcyc   0.2% x1.0 106.2us
    hk_conway_compile         633.56 us 380.16 kcyc   1.0% x1.0 633.6us
    hk_conway_sweep             1.43 ms 858.32 kcyc   2.3% x1.0 1430.5us
  is_ripple_prepare             0.12 us    105 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.44 us  50.70 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           27.20 ms  16.32 Mcyc  43.5%
```

Wall min/avg/max = 52.212/62.495/76.899 ms. All 16 draws are recipe-build draws. The finished-mesh draw scope is inactive. Shared face counters with that inactive parent are emitted as standalone roots; the build-scan parent still includes their raster work.

### Finished mesh / ripple window (frames 401–416)

```text
frame                          62.93 ms  37.76 Mcyc 100.0%
  pov_preserve_half           144.81 us  86.89 kcyc   0.2% x1.0 144.8us
  is_timeline_step             39.30 ms  23.58 Mcyc  62.4%
    is_draw_shape              39.24 ms  23.54 Mcyc  62.3%
      is_mesh_scan             37.47 ms  22.48 Mcyc  59.5%
        scan_mesh_raster       28.54 ms  17.12 Mcyc  45.3%
          filter_blend          1.46 ms 875.42 kcyc   2.3% x18755.4 47cyc
        scan_face_setup         8.66 ms   5.20 Mcyc  13.8% x722.0 12.0us
      is_face_offsets         354.19 us 212.53 kcyc   0.6% x1.0 354.2us
      is_mesh_transform         1.41 ms 846.23 kcyc   2.2% x1.0 1410.4us
  is_ripple_prepare             6.25 us   3.76 kcyc   0.0% x1.0 6.2us
  canvas_clear                 84.12 us  50.49 kcyc   0.1% x1.0 84.1us
  canvas_buffer_wait           23.40 ms  14.04 Mcyc  37.2%
```

Wall min/avg/max = 60.402/62.933/67.052 ms. All 16 frames belong to one shape, have zero recipe-build draws, and issue exactly 16 times its finished face count. This is a full finished-mesh window; ripple deformation may remain active.

### Window containing the peak frame (frames 769–784)

```text
scan_mesh_raster               17.30 ms  10.38 Mcyc  27.2%
  filter_blend                  1.07 ms 644.83 kcyc   1.7% x14452.5 45cyc
scan_face_setup                 5.07 ms   3.04 Mcyc   8.0% x407.0 12.5us
frame                          63.69 ms  38.21 Mcyc 100.0%
  pov_preserve_half           144.06 us  86.45 kcyc   0.2% x1.0 144.1us
  is_timeline_step             26.22 ms  15.73 Mcyc  41.2%
    is_build_draw              22.64 ms  13.58 Mcyc  35.5%
      is_build_scan            22.56 ms  13.54 Mcyc  35.4% x1.0 22558.4us
      is_mesh_transform        75.44 us  45.28 kcyc   0.1% x1.0 75.4us
    hk_conway_compile         354.56 us 212.77 kcyc   0.6% x1.0 354.6us
    hk_conway_sweep           808.25 us 484.98 kcyc   1.3% x1.0 808.2us
  is_ripple_prepare             0.12 us    105 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.38 us  50.66 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           37.24 ms  22.34 Mcyc  58.5%
```

Wall min/avg/max = 48.937/63.689/78.735 ms. This window locates the pass peak. Its averaged tree describes the whole 16-frame window and must not be mistaken for the individual peak frame. Build, fade, or ripple work can coexist within it.

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
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4; 401–416 | 18755.4 | 28.537 | 39.538 | 16 | 43.236 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4; 497–512 | 18229.4 | 27.397 | 36.290 | 16 | 43.894 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2; 577–592 | 16123.9 | 27.154 | 36.145 | 16 | 38.223 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 2; 641–656 | 17489.6 | 26.201 | 35.103 | 16 | 38.048 |
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 2; 801–816 | 18264.3 | 24.889 | 42.601 | 16 | 49.074 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 2; 1441–1456 | 17847.2 | 23.209 | 32.129 | 16 | 36.860 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4; 2017–2032 | 17061.5 | 21.857 | 30.291 | 16 | 33.396 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4; 97–112 | 16524.4 | 21.602 | 27.861 | 16 | 28.682 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1; 1745–1760 | 16389.8 | 21.396 | 29.096 | 16 | 32.744 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2; 1585–1600 | 14702.4 | 20.784 | 23.851 | 16 | 28.609 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 1; 1073–1088 | 18660.9 | 20.672 | 25.878 | 16 | 48.732 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4; 177–192 | 15972.4 | 20.468 | 26.675 | 16 | 28.799 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 2; 1361–1376 | 15102.8 | 18.917 | 23.830 | 16 | 26.613 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 2; 1505–1520 | 14565.8 | 18.079 | 21.373 | 16 | 22.216 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4; 33–48 | 14171.7 | 17.721 | 20.943 | 16 | 23.727 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 2; 865–880 | 13294.2 | 17.503 | 18.991 | 16 | 19.940 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2; 1665–1680 | 15953.2 | 16.359 | 19.074 | 16 | 39.004 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 2; 721–736 | 13632.2 | 16.331 | 19.512 | 16 | 22.764 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 2; 1217–1232 | 17193.6 | 16.035 | 20.250 | 16 | 37.568 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 1; 1297–1312 | 12923.9 | 15.569 | 17.121 | 16 | 18.297 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4; 353–368 | 15808.1 | 14.839 | 18.666 | 16 | 38.071 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 1; 993–1008 | 14420.3 | 13.538 | 15.849 | 16 | 28.327 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 2; 1153–1168 | 14388.9 | 12.479 | 14.455 | 16 | 22.778 |

### Per-pixel figures

The selected heaviest finished-mesh window blends **18755.4 pixels/frame**,
1.809× quadrant coverage. `filter_blend` costs
**46.7 cycles/blend**;
`scan_mesh_raster` costs **912.9 cycles/blended pixel**.
Those scan cycles include probes that were tested and rejected; blended calls
are not a count of all probes or unique pixels. Across the complete capture,
`filter_blend` records 14970.1 blends/frame at
45.3 cycles/blend, and the shared face
raster scope records 701.1 cycles/blend. The annotation-aware reader
includes standalone roots when their latched parent is inactive.

## Column-ISR / DMA marshaling cost

```text
isr_wake            1155.3/f 0.33/1.54/24.00us cpu 2.84%
  isr_pack           143.9/f 0.53/7.08/10.68us cpu 1.63%
  isr_dma_submit     143.9/f 0.58/0.93/11.04us cpu 0.21%
```

The per-call columns are min/weighted-mean/max. CPU percentages use total
captured ISR time divided by the sum of capture-window wall durations.

- Pack plus submit consumes 1.153 ms of CPU per rendered frame; submit alone averages 0.931 µs/call.
- The 72-pixel strobe composite is 600 bytes at 24 MHz: **200 µs** of asynchronous wire time, distinct from the submit CPU cost.
- Pack and submit are nested in `isr_wake`. Its **2.84%** inclusive share must not be added to its children. That share implies approximately **60.72 ms** of foreground opportunity per display window.
- The measured render already includes ISR interruptions. Compare it with the full 62.5 ms window, avoiding a second ISR subtraction: its peak is 0.785× budget, so this captured workload needs no speedup to hold 16 fps.

## Summary ranking

1. `is_timeline_step` — 35.9% of aggregate frame wall time, 22.374 ms/frame.
2. `is_mesh_scan` — 22.8% of aggregate frame wall time, 14.205 ms/frame.
3. `is_build_scan` — 11.4% of aggregate frame wall time, 7.099 ms/frame.

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

Global-O3 peak render is 49.074 ms versus shipping
49.878 ms; both report zero spilled frames. Those are
pass maxima, not a matched-phase speedup ratio. The global-O3 image adds
**+24,536 B FLASH code** and **+13,072 B ITCM code**.
The exact per-frame matched comparison is recorded in the review report.
