# ShapeShifter: finding 77 before/after device comparison

The only source difference is `799ed81d1`, the geodesic replay parameter clamp, compared with `ee19e9dac`. Both trees are clean. This measures the real shipping nine-preset ShapeShifter choreography, without a synthetic fixture or modified dwell. The stale four-shape profiling schedule was replaced by 155-second captures with 1600-revolution epochs, and all captures completed a nine-preset wrap without reinitialization.

Shipping pairs use COM3; global-O3 pairs use COM4. Capture order was before-shipping + after-O3, then after-shipping + before-O3. Every build, flash and capture used `profile_one.sh`, board/tree locks and ELF provenance attestation. Two initial 70-second captures covered only part of the cycle and are excluded from the published evidence and all final results.

Frame 1 is excluded. Comparison rows pair the common frame-number range and verify matching preset ownership; differing trailing capture lengths cannot bias the comparison. Scope and ISR summaries use complete windows starting at frame 17. Timings include interrupts and instrumentation; one paired capture per configuration does not establish a statistical confidence interval.

## Attribution limit

ShapeShifter sets `SAMPLED_RASTER_CONFIG.single_pass = true` in `effects/ShapeShifter.h:548` and passes it to `Plot::rasterize` at line 571. The `if constexpr (SINGLE_PASS)` branch in `core/render/plot/raster.h:648` returns at line 808, before the changed two-pass replay clamp at line 879. Consequently this effect never executes the changed clamp. These measurements describe ShapeShifter's overall image behavior, including run variability and possible compiler/code-layout differences; they do **not** measure the clamp's per-replay-sample execution cost. The hot raster symbols are not all byte-identical (`raster-symbols.json`), so no whole-code identity or statistical-zero claim is made.

No cadence regression was observed: both configurations had zero spilled live frames before and after. Shipping mean increased 0.1931% (+45.888 µs/frame); global O3 increased 0.0016% (+0.380 µs/frame). The paired spherical-polygon presets show shipping mean increases of roughly 0.55–0.61%; those remain observable image differences, not evidence that the bypassed clamp executed.

## Paired frame results

| Configuration | Common frames | Mean before ms | Mean after ms | Mean delta | Peak before ms | Peak after ms | Spilled before/after |
|---|---|---:|---:|---:|---:|---:|---:|
| ship | 2–2458 | 23.7669 | 23.8128 | +0.0459 ms (+0.193%) | 58.376 | 58.395 | 0/2457 → 0/2457 |
| o3 | 2–2458 | 23.8055 | 23.8059 | +0.0004 ms (+0.002%) | 60.028 | 60.047 | 0/2457 → 0/2457 |

## Paired preset means

| Build | Preset | Frames | Before ms | After ms | Delta ms | Delta % |
|---|---|---:|---:|---:|---:|---:|
| ship | 1: Planar star, 208 | 478 | 42.4442 | 42.4472 | +0.0031 | +0.007% |
| ship | 2: Spherical polygon, 74.645 | 299 | 11.4384 | 11.5086 | +0.0701 | +0.613% |
| ship | 3: Planar star, 43.328 | 240 | 7.6678 | 7.6696 | +0.0018 | +0.023% |
| ship | 4: Flower, 70 | 240 | 34.2683 | 34.2862 | +0.0179 | +0.052% |
| ship | 5: Planar star, 72 | 240 | 9.4043 | 9.4064 | +0.0021 | +0.022% |
| ship | 6: Spherical polygon, 128 | 240 | 17.0325 | 17.1255 | +0.0930 | +0.546% |
| ship | 7: Spherical polygon, 144 / 4 sides | 240 | 19.2793 | 19.3882 | +0.1088 | +0.565% |
| ship | 8: Spherical polygon, 144 / 3.195 sides | 240 | 22.2118 | 22.3471 | +0.1353 | +0.609% |
| ship | 9: Flower, 72 | 240 | 34.6644 | 34.6818 | +0.0174 | +0.050% |
| o3 | 1: Planar star, 208 | 478 | 44.2274 | 44.2281 | +0.0007 | +0.002% |
| o3 | 2: Spherical polygon, 74.645 | 299 | 11.5394 | 11.5383 | -0.0011 | -0.009% |
| o3 | 3: Planar star, 43.328 | 240 | 6.6014 | 6.6049 | +0.0035 | +0.053% |
| o3 | 4: Flower, 70 | 240 | 32.9610 | 32.9621 | +0.0011 | +0.003% |
| o3 | 5: Planar star, 72 | 240 | 8.4913 | 8.4911 | -0.0003 | -0.003% |
| o3 | 6: Spherical polygon, 128 | 240 | 17.2108 | 17.2102 | -0.0006 | -0.003% |
| o3 | 7: Spherical polygon, 144 / 4 sides | 240 | 19.5698 | 19.5696 | -0.0002 | -0.001% |
| o3 | 8: Spherical polygon, 144 / 3.195 sides | 240 | 22.4347 | 22.4347 | -0.0000 | -0.000% |
| o3 | 9: Flower, 72 | 240 | 33.9772 | 33.9775 | +0.0003 | +0.001% |

## Memory and layout

| Image | RAM1 code before/after | RAM1 variables before/after | FLASH data before/after | FLASH code before/after |
|---|---:|---:|---:|---:|
| Full shipping Phantasm | 194,776 → 194,840 (+64 B) | 314,816 → 314,816 (+0 B) | 725,384 → 725,384 (+0 B) | 501,016 → 501,080 (+64 B) |
| ShapeShifter shipping profile | 49,656 → 49,656 (+0 B) | 315,072 → 315,072 (+0 B) | 152,148 → 152,148 (+0 B) | 77,112 → 77,112 (+0 B) |
| ShapeShifter O3 profile | 74,152 → 74,152 (+0 B) | 315,072 → 315,072 (+0 B) | 152,028 → 152,028 (+0 B) | 106,400 → 106,400 (+0 B) |

Every full shipping and single-effect profile build passed its size/layout gates. Full shipping RAM1 padding and remaining bank-boundary headroom are 1768 B after the clamp.

## Reports and evidence

The evidence below is committed with this report. Original capture text and
provenance retain their recorded local paths; those paths are historical metadata,
not required downloads. ELF binaries, environment dumps, build logs, and temporary
README proposals remain local and are not published.

- [Shipping profile](../shipping/profile_shapeshifter_teensy_2026-09-25.md)
- [Global-O3 reference](../O3/profile_shapeshifter_teensy_2026-09-25.md)
- [Paired metrics](metrics.json), [code-section hashes](code-sections.json), and [raster-symbol hashes](raster-symbols.json).
- [Before shipping capture](data/before-ship.log.txt), [after shipping capture](data/after-ship.log.txt), [before O3 capture](data/before-o3.log.txt), and [after O3 capture](data/after-o3.log.txt).
- Matching `.provenance.txt` files in `data/` record compiler, source, and ELF hashes.

The candidate reports supersede the historical August26 ShapeShifter reports. Their timing change cannot be attributed solely to finding77; the paired comparison isolates the source revision while measuring the overall ShapeShifter image; it does not isolate replay-clamp execution cost.
