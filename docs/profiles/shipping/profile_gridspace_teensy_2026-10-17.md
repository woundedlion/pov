# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-17, **selective -O3**)

Historical matched temporal A/B comparison at the source snapshots below. All measurements and implementation descriptions in this supplement refer to those snapshots. [The current canonical report](profile_gridspace_teensy_2026-10-23.md) supplies current shipping timing.
Final raw captures are `C:/work/temp/gridspace-temporal-final.log`,
`C:/work/temp/gridspace-optimized-off-final.log`, and
`C:/work/temp/gridspace-original-off-final.log`. Captured on COM4 at
2026-10-17 16:58, 16:55, and 16:32 America/Los_Angeles, respectively.
This replaces the earlier October 17 temporal report at the same path.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | Shipping `profile` environment: `-Os`, newlib-nano, selective hot pullback shading/color regions and hot flash column traversal |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, clean profiled source `095ad16d969a0c90f533b893acfdb77905f68cd1` |
| Matched baseline | Same optimized source with only `GridSpaceSpec::TEMPORAL_SAMPLING` set false; includes the shared geometry optimization |
| Method | One 70 s capture per image, 32-frame windows; no epoch stretch or dwell compression |
| Reproduce | Pin `HS_PROFILE_TREE` to the selected checkout and `HS_TEENSY_PORT=COM4`; run `bash tools/profile_one.sh GridSpace profile 70 32` |

Candidate image size (`profile`, GridSpace only):
`FLASH: code:85740, data:155280, headers:8832` /
`RAM1: variables:313024, code:22344, padding:10424, free:178496` /
`RAM2: variables:520064, free:4224`.

The full shipping `phantasm` size/layout gate passes. The memory baseline below
is the previous temporal implementation, not the temporal-off timing control.

| Memory | Previous temporal bytes | Optimized temporal bytes | Delta bytes |
|---|---:|---:|---:|
| RAM1 code | 179,816 | 181,704 | +1,888 |
| RAM1 variables | 312,736 | 312,736 | 0 |
| FLASH data | 851,388 | 851,388 | 0 |
| FLASH total | 1,498,108 | 1,511,420 | +13,312 |
| RAM2 variables | 520,064 | 520,064 | 0 |
| GridSpace persistent arena, host probe | 25,536 | 77,816 | +52,280 |

The final host arena charge, including common overhead, is 110,584 bytes;
measured host peak stack is 3,167 bytes. These are host-ABI measurements, not
device allocations. Retained endpoint and float lookup caches use the existing
reserved arena. The device's fixed RAM1 variable allocation and 14,944-byte
stack allowance are unchanged. The original non-temporal host persistent
high-water mark was 13,584 bytes.

All three final logs pass parser validation: 34 complete windows, 1,088 frames,
correct effect/resolution, monotonic numbering, and no epoch reset. For candidate
frames 1–32, root cycles 1,161,312,721 / 600 give 1,935,521.202 us versus wall
sum 1,935,525 us: **2.0 ppm**. Matched-off and original-off captures agree
within **1.9 ppm** and **1.4 ppm**, respectively.

The wrapper attested GCC 15.2.1 and matching profile/shipping compiler and ABI.
Candidate profile ELF SHA-256 is
`2563c69ec02ee1890fe9a12a11e8faa5f5983ac3f060a4ffca1f6de817462aef`;
matched-off ELF SHA-256 is
`225f34158fb27ccc9218cd2f56b141baa0c33dedeabf772a1d8b6a54b6eb809a`.
Raw logs, build logs, environment dumps, ELF files, and provenance are local
artifacts under `C:/work/temp`, with candidate artifacts in
`artifacts/gridspace-temporal-final_095ad16d969a_2563c69ec02e`.
The matched-off working tree was based on `47d2724f989975d9310b1de1a539a9cb167a7a8e`
plus the optimization patch. `gridspace-final-pair-source-check.txt` confirms
the only compiled difference from the final candidate was the temporal flag.
Its retained patch is `gridspace-off-final-archive.patch`. The candidate was
subsequently rebased to `efc5bc974`; only peer profile-index documentation differs
from the captured source.

## Frame cadence

**Initialization/priming draw: 38.694 ms.** The driver primes a full-width
288×72 band before publishing the effect. This is construction cost, excluded
from the live-frame statistics below. Matched-off initialization is 36.067 ms;
original-off initialization is 40.265 ms.

**Live frames 2–1,088:** render averages **24.754 ms**, peak **25.409 ms** at
frame 415, and **0/1,087 spilled (0.00%)**. A half-turn is 62.5 ms; each live
frame emits one 144×72 quadrant (10,368 pixels). The observed cadence remains
**16 fps**, with 37.091 ms between the worst live render and the deadline.
`canvas_buffer_wait` is the round-up idle until the display flip.

The unfiltered parser footer includes setup frame 1: peak 38.694 ms and
0/1,088 spills. The ranked indexes follow the archive convention of excluding
that setup draw. Over complete steady windows (frames 33–1,088),
`fx_shader_draw` averages **18.005 ms/frame**; its worst window is
**18.145 ms/frame**, frames 65–96.

### Matched optimized-code temporal overhead

| Work | Optimized off ms/frame | Optimized on ms/frame | Delta ms | Delta |
|---|---:|---:|---:|---:|
| Live render | 22.268 | 24.754 | +2.485 | +11.2% |
| Shading and column preparation | 15.648 | 18.005 | +2.357 | +15.1% |
| Frame preparation | 4.060 | 3.857 | −0.203 | −5.0% |
| Animation advancement | 2.236 | 2.308 | +0.072 | +3.2% |
| Live peak render | 22.825 | 25.409 | +2.584 | +11.3% |

Render means use individual live-frame rows. Scope means omit the complete
startup window because its counters cannot separate individual frames.
Both images include the shared factored-geometry path; disabling the compile-time
temporal flag also removes its retained endpoint/cache work. This is the
incremental temporal cost after shared optimization. **The under-5% incremental
overhead target was not reached:** the measured mean overhead is **11.2%**.
Both captures have zero live spills.

### Total rendering improvement

| Image | Live mean ms | Live peak ms | Spilled | Initialization ms |
|---|---:|---:|---:|---:|
| Original single-time renderer, final 70 s capture | 24.222 | 24.581 | 0/1087 (0.00%) | 40.265 |
| Previous temporal renderer, historical 70 s capture | 36.130 | 36.615 | 0/1087 (0.00%) | 62.977 |
| Optimized single-time renderer, final 70 s capture | 22.268 | 22.825 | 0/1087 (0.00%) | 36.067 |
| Optimized temporal renderer, final 70 s capture | 24.754 | 25.409 | 0/1087 (0.00%) | 38.694 |

The optimized temporal renderer is **31.5% faster** than the previous temporal
renderer, saving **11.376 ms/frame**. Against the original single-time renderer,
it costs **0.532 ms/frame (2.2%)** more. That smaller percentage includes shared
geometry/color improvements and is not the matched temporal-overhead figure.

The historical temporal row comes from the earlier 70 s COM4 capture at
`5cbdd7d91c493c5af7d84d5ecadf61d1bcd4cc22`, retained locally as
`grid-temporal-profile-20261017/candidate_ship.log`. The separate short
36.068 ms temporal measurement is not used in these 70 s comparisons.
The final original-off control uses `47d2724f989975d9310b1de1a539a9cb167a7a8e`
with temporal sampling disabled and no optimization patch.

## Phase-by-phase readout

Phase schedule: one preset with continuous camera, warp, and color animation.
Initialization is reported separately above.

### Steady peak window (frames 385–416)

```
frame                    62.392 ms  37.435 Mcyc 100.0%
  pov_preserve_half     148.06 us   88.85 kcyc   0.2% x1 148.1 us/call
  fx_shader_draw         17.907 ms  10.744 Mcyc  28.7% x1 17907.1 us/call
    fx_prepare_column   599.84 us  359.91 kcyc   1.0% x146 4.11 us/call
  fx_prepare_frame        4.005 ms   2.403 Mcyc   6.4% x1 4005.3 us/call
  fx_advance              2.310 ms   1.386 Mcyc   3.7% x1 2309.6 us/call
  fx_timeline_step       70.25 us   42.16 kcyc   0.1% x1 70.3 us/call
  canvas_clear           85.22 us   51.13 kcyc   0.1% x1 85.2 us/call
  canvas_buffer_wait     37.584 ms  22.551 Mcyc  60.2% x1 37584.3 us/call
```

Wall min/avg/max = 61.406/62.392/63.404 ms. Wall variation includes sync wait
and is not render spill. Column preparation is included in the shading scope;
the observed 146 preparations include the clip's border columns. Retaining
endpoint resources lies outside the named shading/preparation scopes and is
included in total render. Over all complete steady windows, column preparation
averages 0.596 ms/frame.

The generic meridian basis removes repeated rotation work from both on/off
paths. Temporal preparation shares quaternion arcs and periodic clock deltas.
Float copies of quantized hue/noise tables, shared endpoint addressing, and
prepared color kernels reduce repeated lookup conversion and branching.
Equal noise tables and constant endpoint alpha have shared sampling paths.

### Per-pixel figures

There is no `filter_blend` counter: this shader writes pixels directly.
Steady shading cost is **1,042 cycles per nominal output pixel**, versus
1,716 in the previous temporal capture, dividing the inclusive scope by
10,368 emitted pixels. These figures include column preparation, border
traversal, and ISR interference; they are not isolated shader instruction counts.

## Column-ISR / DMA marshaling cost

Peak-window rates are calls per 32 rendered frames. CPU shares use the separately
logged 1,999,960 us ISR interval:

```
counter          calls/frame  min/avg/max per call     CPU
isr_wake        1152.06/frame  0.50/1.95/20.48 us  3.60%
isr_pack         144.00/frame  8.39/9.17/12.17 us  2.11%
isr_dma_submit   144.00/frame  0.64/0.94/7.40 us   0.21%
```

- DMA submission averages 0.936 us/call; packing averages 9.170 us/call.
- The 72-LED image plus black strobe frame occupies 600 bytes. At the configured
  24 MHz SPI rate, the driver budgets 230 us including byte framing, within the
  approximately 434 us column period. The transfer is asynchronous.
- Inclusive `isr_wake` consumes 3.60%, leaving about 60.25 ms of foreground
  CPU time per 62.5 ms window. Pack and submit run inside wake and are not added
  to it. Render measurements already include ISR interference; no speedup is
  required to sustain this observed cadence. DMA-completion and other interrupts
  are not measured by these counters.

## Summary ranking

README cells: peak 🟢 25.41, spilled 🟢 0/1087 (0.00%).

1. `fx_shader_draw` — 28.7% of the peak window, 17.907 ms/frame, including
   0.600 ms of column preparation.
2. `fx_prepare_frame` — 6.4%, 4.005 ms/frame.
3. `fx_advance` — 3.7%, 2.310 ms/frame.
4. `pov_preserve_half` — 0.2%, 0.148 ms/frame.

Device render cost fell 31.5% from the previous temporal implementation.
The final result is 2.2% above the original renderer, but temporal sampling
still adds 11.2% relative to the optimized single-time path. No native or WASM
timing result is substituted for these device measurements.

The final WASM bundle passed actual-browser checks at 96×20 and 288×144:
single-engine and four-worker fidelity changes, paused refresh, visual response,
removed presentation/exposure controls, and legacy URL handling. These are
functional checks, not device performance or perceived-seam measurements.

## Caveats

- CYCCNT runs through interrupts; all foreground scopes include ISR time.
- Frame 1 primes the effect before publication and is initialization cost,
  not a missed live display frame.
- There are no per-pixel profiling scopes or `filter_blend` parenting effects.
  The column-preparation scope remains included in the measured image.
- The geometry optimization changes floating-point evaluation order. GridSpace's
  authored-preset tests bound the difference to 16 RGB16 codes maximum and
  0.05 mean; temporal color tests allow two RGB16 codes against independently
  rounded endpoints. Alpha continues to interpolate its endpoints.
- Selective-O3 shipping behavior includes hot pullback shading/color and column
  traversal. Global-O3 timing and historical code-size comparison cells were
  not refreshed.
- The final candidate capture used a clean commit. The two off controls use
  the profiling-only source states described above. No dwell compression or
  epoch override was used.
- One sequential capture per image measures this 70-second trajectory; it does
  not establish a worst-case bound for all GUI parameter settings.
- Optimization support is generic, but this capture validates GridSpace only.
  It does not establish temporal continuity or performance for other effects.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=GridSpace` and
`HS_PROFILE_WINDOW=32`. `just profile GridSpace 70` runs the supported locked
build, flash, and capture wrapper. Pin `HS_PROFILE_TREE` and `HS_TEENSY_PORT`
for the selected source and board. Validate with
`python tools/parse_profile.py <log> validate`.
