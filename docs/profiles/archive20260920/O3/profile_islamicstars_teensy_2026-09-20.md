# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-20, **-O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/review_20260920/islamic_baseline_o3.log`;
captured **2026-09-20 00:21 Pacific Daylight Time**, board **COM3**. This replaces the
2026-08-26 baseline report for this configuration. The roster currently contains
**23** shapes, not the 24 in the older skill table. See the paired
[review comparison](../../review_2026-09-20.md) for the experimental corrections.

Shipping sibling: [selective-O3 report](../shipping/profile_islamicstars_teensy_2026-09-20.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile_o3`; global `-O3 -ffast-math`, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, DMA LEDs |
| Effect | IslamicStars 288×144, one-entry playlist; source `d9a39abd8cd69fc6ed91175a6e0398e0d71be961` |
| Method | `HS_PROFILE`, 16-frame windows, 210-second capture, Trans Speed 4, epoch 1920 revolutions (240 s); all 23 shapes and wrap observed |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<source-tree> bash tools/profile_one.sh IslamicStars profile_o3 210 16 '-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4'` |

Image size: `FLASH: code:153496, data:193280, headers:8552   free for files:1676288 / RAM1: variables:315552, code:57016, padding:8520   free for local variables:143200 / RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 1–16, root
`542202851` cycles ÷ 600 MHz versus measured wall sum
`903672` µs differ by **0.64 ppm**.
`tools/parse_profile.py ... validate` passes marker wrap, monotonic frames,
render/wall separation, and this counter check. Source/compiler/ELF attestations
are appended to the raw log; the source was clean baseline `d9a39abd`.

## Frame cadence

**Pass aggregate**: render mean **22.143 ms/frame**;
worst window mean **44.843 ms/frame**
(frames 2577–2592); peak frame render
**50.412 ms** (frames 2801–2816);
spilled **0/3328 frames**
(0.0%). All 23 shape buckets are green.

At 480 RPM each half-revolution lasts 62.5 ms, yielding 16 fps. One board renders
a 144×72 quadrant, about 10,368 pixels. Every captured build, finished-mesh,
ripple and fade phase fits that window; the measured peak leaves
**12.088 ms** of render headroom. `canvas_buffer_wait`
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
scan_mesh_raster               22.82 ms  13.69 Mcyc  36.5%
  filter_blend                  1.55 ms 929.72 kcyc   2.5% x19515.8 48cyc
scan_face_setup                 9.46 ms   5.67 Mcyc  15.1% x1082.0 8.7us
frame                          62.44 ms  37.47 Mcyc 100.0%
  pov_preserve_half           142.38 us  85.43 kcyc   0.2% x1.0 142.4us
  is_timeline_step             36.07 ms  21.64 Mcyc  57.8%
    is_build_draw              32.77 ms  19.66 Mcyc  52.5%
      is_build_scan            32.66 ms  19.59 Mcyc  52.3% x1.0 32656.8us
      is_mesh_transform       106.88 us  64.16 kcyc   0.2% x1.0 106.9us
    hk_conway_compile         635.06 us 381.04 kcyc   1.0% x1.0 635.1us
    hk_conway_sweep             1.43 ms 859.36 kcyc   2.3% x1.0 1432.2us
  is_ripple_prepare             0.12 us    105 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.38 us  50.63 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           26.15 ms  15.69 Mcyc  41.9%
```

Wall min/avg/max = 55.155/62.442/75.529 ms. All 16 draws are recipe-build draws. The finished-mesh draw scope is inactive. Shared face counters with that inactive parent are emitted as standalone roots; the build-scan parent still includes their raster work.

### Finished mesh / ripple window (frames 401–416)

```text
frame                          62.93 ms  37.76 Mcyc 100.0%
  pov_preserve_half           142.75 us  85.67 kcyc   0.2% x1.0 142.8us
  is_timeline_step             39.12 ms  23.47 Mcyc  62.2%
    is_draw_shape              39.05 ms  23.43 Mcyc  62.1%
      is_mesh_scan             37.29 ms  22.37 Mcyc  59.3%
        scan_mesh_raster       28.35 ms  17.01 Mcyc  45.1%
          filter_blend          1.49 ms 892.40 kcyc   2.4% x18755.4 48cyc
        scan_face_setup         8.66 ms   5.20 Mcyc  13.8% x722.0 12.0us
      is_face_offsets         354.38 us 212.63 kcyc   0.6% x1.0 354.4us
      is_mesh_transform         1.41 ms 846.37 kcyc   2.2% x1.0 1410.6us
  is_ripple_prepare             6.31 us   3.79 kcyc   0.0% x1.0 6.3us
  canvas_clear                 84.44 us  50.69 kcyc   0.1% x1.0 84.4us
  canvas_buffer_wait           23.58 ms  14.15 Mcyc  37.5%
```

Wall min/avg/max = 60.437/62.928/66.974 ms. All 16 frames belong to one shape, have zero recipe-build draws, and issue exactly 16 times its finished face count. This is a full finished-mesh window; ripple deformation may remain active.

### Window containing the peak frame (frames 2801–2816)

```text
frame                          64.80 ms  38.88 Mcyc 100.0%
  pov_preserve_half           144.38 us  86.63 kcyc   0.2% x1.0 144.4us
  is_timeline_step             27.93 ms  16.76 Mcyc  43.1%
    is_build_draw              23.20 ms  13.92 Mcyc  35.8%
      is_build_scan            23.16 ms  13.90 Mcyc  35.7% x0.8 30885.9us
      is_mesh_transform        32.50 us  19.52 kcyc   0.1% x0.8 43.3us
    hk_conway_compile         181.69 us 109.03 kcyc   0.3% x0.8 242.2us
    hk_conway_sweep           477.75 us 286.68 kcyc   0.7% x0.8 637.0us
    is_draw_shape               2.69 ms   1.62 Mcyc   4.2%
      is_mesh_scan              2.68 ms   1.61 Mcyc   4.1%
        scan_mesh_raster       22.84 ms  13.70 Mcyc  35.2%
          filter_blend          1.25 ms 748.20 kcyc   1.9% x16324.2 46cyc
        scan_face_setup         2.86 ms   1.72 Mcyc   4.4% x287.0 10.0us
      is_face_offsets           6.94 us   4.17 kcyc   0.0% x0.2 27.8us
      is_mesh_transform         2.12 us   1.30 kcyc   0.0% x0.2 8.5us
  is_ripple_prepare             0.12 us    106 cyc   0.0% x1.0 0.1us
  canvas_clear                 84.50 us  50.70 kcyc   0.1% x1.0 84.5us
  canvas_buffer_wait           36.64 ms  21.99 Mcyc  56.5%
```

Wall min/avg/max = 48.244/64.802/75.555 ms. This window locates the pass peak. Its averaged tree describes the whole 16-frame window and must not be mistaken for the individual peak frame. Build, fade, or ripple work can coexist within it.

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
| 5 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 2160/2880/722/5760 | 4; 401–416 | 18755.4 | 28.350 | 39.351 | 16 | 43.020 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 1620/2160/542/4320 | 4; 497–512 | 18229.4 | 27.234 | 36.132 | 16 | 43.676 |
| 10 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 3240/4320/1082/8640 | 4; 2577–2592 | 18621.5 | 27.136 | 44.843 | 16 | 48.946 |
| 7 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 1620/2160/542/4320 | 2; 577–592 | 16123.9 | 26.989 | 35.985 | 16 | 38.052 |
| 8 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 1620/2160/542/4320 | 4; 641–656 | 17489.6 | 26.040 | 34.936 | 16 | 37.841 |
| 18 | truncatedIcosahedron_hk54_ambo_hk72 | 1620/2160/542/4320 | 4; 1441–1456 | 17847.2 | 23.069 | 31.994 | 16 | 36.665 |
| 3 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 1620/2160/542/4320 | 4; 2017–2032 | 17061.5 | 21.746 | 30.164 | 16 | 33.320 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 990/1440/452/2880 | 4; 97–112 | 16524.4 | 21.468 | 27.735 | 16 | 28.542 |
| 22 | icosahedron_snub_relax_truncate033_hankin62 | 1350/1800/452/3600 | 1; 1745–1760 | 16389.8 | 21.287 | 28.977 | 16 | 32.580 |
| 20 | dodecahedron_hk72_ambo_dual_hk20 | 540/720/182/1440 | 2; 1585–1600 | 14702.4 | 20.655 | 23.717 | 16 | 28.447 |
| 13 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 542/1080/540/2160 | 2; 2849–2864 | 18723.1 | 20.554 | 25.762 | 16 | 50.412 |
| 2 | dodecahedron_ambo_bevel33_relax_hk66 | 1080/1440/362/2880 | 4; 177–192 | 15972.4 | 20.355 | 26.550 | 16 | 28.669 |
| 19 | dodecahedron_hk54_ambo_hk72 | 540/720/182/1440 | 4; 3281–3296 | 14587.6 | 19.324 | 22.630 | 16 | 25.747 |
| 17 | rhombicuboctahedron_hk63_ambo_hk63 | 864/1152/290/2304 | 4; 3137–3152 | 15156.8 | 19.004 | 23.931 | 16 | 26.559 |
| 11 | octahedron_hk17_ambo_hk73 | 216/288/74/576 | 4; 2641–2656 | 13333.4 | 17.637 | 19.126 | 16 | 22.319 |
| 0 | dodecahedron_hk62_ambo_hk62 | 540/720/182/1440 | 4; 33–48 | 14171.7 | 17.625 | 20.847 | 16 | 23.671 |
| 9 | icosahedron_ambo_truncate033_hankin59 | 540/720/182/1440 | 4; 2481–2496 | 13933.4 | 17.166 | 20.240 | 16 | 22.701 |
| 21 | truncatedIcosahedron_truncate50d_ambo_dual | 272/540/270/1080 | 2; 1665–1680 | 15953.2 | 16.265 | 18.978 | 16 | 38.701 |
| 15 | snubDodecahedron_truncate5d_ambo_dual | 452/900/450/1800 | 4; 1217–1232 | 17193.6 | 15.935 | 20.151 | 16 | 37.570 |
| 16 | octahedron_hk34_ambo_hk72 | 216/288/74/576 | 2; 3073–3088 | 13040.4 | 15.630 | 17.196 | 16 | 18.190 |
| 4 | dodecahedron_bevel2_relax_gyro | 542/900/360/1800 | 4; 353–368 | 15808.1 | 14.760 | 18.581 | 16 | 37.945 |
| 12 | icosahedron_kis_gyro | 272/450/180/900 | 2; 993–1008 | 14420.3 | 13.456 | 15.769 | 16 | 28.225 |
| 14 | icosidodecahedron_truncate5d_ambo_dual | 182/360/180/720 | 4; 2929–2944 | 14401.7 | 12.567 | 14.540 | 16 | 23.707 |

### Per-pixel figures

The selected heaviest finished-mesh window blends **18755.4 pixels/frame**,
1.809× quadrant coverage. `filter_blend` costs
**47.6 cycles/blend**;
`scan_mesh_raster` costs **906.9 cycles/blended pixel**.
Those scan cycles include probes that were tested and rejected; blended calls
are not a count of all probes or unique pixels. Across the complete capture,
`filter_blend` records 14867.6 blends/frame at
46.1 cycles/blend, and the shared face
raster scope records 693.6 cycles/blend. The annotation-aware reader
includes standalone roots when their latched parent is inactive.

## Column-ISR / DMA marshaling cost

```text
isr_wake            1154.2/f 0.36/1.55/24.01us cpu 2.86%
  isr_pack           143.9/f 0.53/7.07/10.58us cpu 1.63%
  isr_dma_submit     143.9/f 0.58/0.93/9.76us cpu 0.21%
```

The per-call columns are min/weighted-mean/max. CPU percentages use total
captured ISR time divided by the sum of capture-window wall durations.

- Pack plus submit consumes 1.152 ms of CPU per rendered frame; submit alone averages 0.931 µs/call.
- The 72-pixel strobe composite is 600 bytes at 24 MHz: **200 µs** of asynchronous wire time, distinct from the submit CPU cost.
- Pack and submit are nested in `isr_wake`. Its **2.86%** inclusive share must not be added to its children. That share implies approximately **60.71 ms** of foreground opportunity per display window.
- The measured render already includes ISR interruptions. Compare it with the full 62.5 ms window, avoiding a second ISR subtraction: its peak is 0.807× budget, so this captured workload needs no speedup to hold 16 fps.

## Summary ranking

1. `is_timeline_step` — 35.1% of aggregate frame wall time, 21.909 ms/frame.
2. `is_mesh_scan` — 22.1% of aggregate frame wall time, 13.786 ms/frame.
3. `is_build_scan` — 11.4% of aggregate frame wall time, 7.084 ms/frame.

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

### Global -O3 vs selective -O3

Global-O3 peak render is 50.412 ms versus shipping
50.828 ms; both report zero spilled frames. Those are
pass maxima, not a matched-phase speedup ratio. The global-O3 image adds
**+24,616 B FLASH code** and **+13,152 B ITCM code**.
The exact per-frame matched comparison is recorded in the review report.
