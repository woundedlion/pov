# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, copied with original mtime from `C:/work/temp/flash-20261030-gs-after.log`. This replaces `profile_gsreactiondiffusion_teensy_2026-10-23.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; flywheel and DMA ISRs live |
| Image | Shipping selective-O3 `profile`, effect render/physics and Scan hot regions active; no deep, replay, PC or count instrumentation |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GSReactionDiffusion 288x144, single-entry playlist, source `480cd8ac1b417fcbe381ac645f4b5761268ab8e5` |
| Method | 130 s, window 32; `HS_PROFILE_EPOCH_REVS=1200` stretches the instance to 150 s, avoiding an epoch crossing |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-flash-neighbors-20261030 HS_TEENSY_PORT=COM4 HS_TEENSY_TOOLS=C:/work/temp/teensy-tools-loopback-20261030 bash tools/profile_one.sh GSReactionDiffusion profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Profile image size: `FLASH: code:87976, data:248996, headers:9140   free for files:1685504 / RAM1: variables:313248, code:17864, padding:14904   free for local variables:178272 / RAM2: variables:520064  free for malloc/new:4224`.

The shipping image uses FLASH code 1,021,896 B, data 749,104 B and headers 8,712 B, total 1,779,712 B. Against parent `93def4f2a`, Flash falls 91,136 B: data -91,880 B, code +272 B, headers +472 B. RAM1 code changes 187,160 to 187,192 B (+32 B), with padding shrinking by 32 B; bank allocation and 14,944 B stack headroom are unchanged. RAM1 variables remain 312,736 B and RAM2 variables remain 520,064 B.

Exactness cross-check: frames 1729-1760 root cycles / 600 MHz = 1997514.362 us; measured wall sum 1997511 us, delta 1.68 ppm. `parse_profile.py ... validate` passes frame coverage and no-reset checks. Captured 2026-10-30 09:50, America/Los_Angeles.

The matched COM4 baseline `flash-20261030-gs-before.log` ran source `aaef299bd` with the same effect, flags, duration, driver and epoch stretch. The candidate includes the independently landed cold-registration change and pullback/Millefiori work; those do not change these effects' simulation equations or authored configuration. This comparison measures the complete candidate against the recorded baseline, not the old shorter canonical capture.

## Frame cadence

Complete-window pass aggregate (`parse_profile.py ... windows`): 53.504 ms peak, 0/2048 spilled. Setup frame 1 is reported separately in the ranked index.

| COM4 metric | Baseline | Candidate |
|---|--:|--:|
| Setup render ms | 25.709 | 25.722 |
| Steady peak render ms | 53.470 | 53.504 |
| Steady spilled frames | 0/2047 | 0/2047 |

Steady peak changes +0.06%. A display window is 62.5 ms, giving 16 fps for steady frames; the active quadrant is 10,368 pixels. `canvas_buffer_wait` is the round-up idle until the next display flip. The peak steady frame has 8.996 ms margin. Setup and every steady frame fit one display window.

All raw telemetry, including the trailing incomplete counter window, contains 2056 frames: peak 53.504 ms and 0 spills. The ranked index uses complete counter windows, frames 2-2048, matching its declared archive convention.

README cells: peak 🟢 53.504, spilled 🟢 0/2047 (0.00%).

## Phase-by-phase readout

The simulation evolves continuously. Startup, peak post-setup load and the late captured state sample its observed cost regimes; no preset-cycling marker is expected. Epoch stretching changes instance lifetime, not the per-frame physics or shading work.

### Startup (frames 1-32)

```
frame                    60.195 ms  36.117 Mcyc 100.0%
  pov_preserve_half      136.69 us   82.02 kcyc   0.2% x0.97 141 us/call
  grd_render             25.513 ms  15.308 Mcyc  42.4%
    grd_shader_draw       9.776 ms   5.866 Mcyc  16.2% x6 1629 us/call
    grd_cull_flags       691.43 us  414.86 kcyc   1.1% x6 115 us/call
    grd_publish           1.966 ms   1.180 Mcyc   3.3% x6 328 us/call
    grd_physics           5.185 ms   3.111 Mcyc   8.6% x6 864 us/call
    grd_pigment           3.102 ms   1.861 Mcyc   5.2% x1 3102 us/call
    rd_lattice_basis     114.70 us   68.82 kcyc   0.2% x1 115 us/call
    grd_color_noise       3.775 ms   2.265 Mcyc   6.3% x1 3775 us/call
  rd_timeline_step        94.21 us   56.53 kcyc   0.2% x1 94 us/call
  canvas_clear            95.60 us   57.36 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     34.331 ms  20.599 Mcyc  57.0% x1 34331 us/call
grd_color_palette        426.17 us  255.70 kcyc   0.7% x0.03 13638 us/call
grd_seed_reaction        422.88 us  253.73 kcyc   0.7% x0.03 13532 us/call
```

Wall min/avg/max = 22.213/60.194/65.105 ms. The first draw includes startup work; the other frames in this window already follow the normal display cadence. The individual first-frame value is separated from this window mean.

### Peak post-setup load (frames 1729-1760)

```
frame                    62.422 ms  37.453 Mcyc 100.0%
  pov_preserve_half      134.31 us   80.59 kcyc   0.2% x1 134 us/call
  grd_render             51.484 ms  30.890 Mcyc  82.5%
    grd_shader_draw      31.112 ms  18.667 Mcyc  49.8% x6 5185 us/call
    grd_cull_flags       709.16 us  425.49 kcyc   1.1% x6 118 us/call
    grd_publish           4.638 ms   2.783 Mcyc   7.4% x6 773 us/call
    grd_physics           5.212 ms   3.127 Mcyc   8.3% x6 869 us/call
    grd_pigment           4.259 ms   2.555 Mcyc   6.8% x1 4259 us/call
    rd_lattice_basis     112.39 us   67.43 kcyc   0.2% x1 112 us/call
    grd_color_noise       3.813 ms   2.288 Mcyc   6.1% x1 3813 us/call
  rd_timeline_step        92.39 us   55.44 kcyc   0.1% x1 92 us/call
  canvas_clear            95.40 us   57.24 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     10.592 ms   6.355 Mcyc  17.0% x1 10592 us/call
```

Wall min/avg/max = 59.566/62.422/65.355 ms. This window contains the largest post-setup render frame. Its mean render and buffer wait are distinct from the individual-frame peak used for cadence.

### Late simulation (frames 2017-2048)

```
frame                    62.411 ms  37.446 Mcyc 100.0%
  pov_preserve_half      136.02 us   81.61 kcyc   0.2% x1 136 us/call
  grd_render             43.279 ms  25.967 Mcyc  69.3%
    grd_shader_draw      26.344 ms  15.806 Mcyc  42.2% x6 4391 us/call
    grd_cull_flags       702.20 us  421.32 kcyc   1.1% x6 117 us/call
    grd_publish           1.947 ms   1.168 Mcyc   3.1% x6 325 us/call
    grd_physics           5.204 ms   3.123 Mcyc   8.3% x6 867 us/call
    grd_pigment           4.248 ms   2.549 Mcyc   6.8% x1 4248 us/call
    rd_lattice_basis     111.90 us   67.14 kcyc   0.2% x1 112 us/call
    grd_color_noise       3.809 ms   2.285 Mcyc   6.1% x1 3809 us/call
  rd_timeline_step        91.39 us   54.84 kcyc   0.1% x1 91 us/call
  canvas_clear            93.70 us   56.22 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     18.787 ms  11.272 Mcyc  30.1% x1 18787 us/call
```

Wall min/avg/max = 57.698/62.410/67.290 ms. This is the last complete counter window of the 130-second pass. It extends the observation beyond the previous shorter report and retains the changing simulation workload.

### Per-pixel figures

The quadrant holds 10,368 output pixels. This image does not emit a complete blended-pixel counter for the direct shader path, so blended coverage and cycles per blend cannot be inferred from draw-call counts. Counter instrumentation overhead remains included.

## Column-ISR / DMA marshaling cost

Representative peak post-setup window; rates, min/avg/max CPU times and measured ISR-window CPU shares:

```
isr_wake          1152.8/f 0.543/1.991/26.710 us 3.67%
  isr_pack         144.1/f 8.078/9.273/18.830 us 2.13%
  isr_dma_submit   144.1/f 0.641/0.948/6.016 us 0.21%
```

- DMA submission is CPU setup; pixel packing is separately measured inside wake.
- At 72 LEDs the image buffer is 300 B and its strobed composite is 600 B. At the configured SPI divider, `dma::transfer_us` budgets 115/230 us including per-byte gaps, outside submit CPU time.
- Inclusive wake CPU share is 3.67%, leaving about 60.206 ms foreground capacity per display window. The measured inclusive steady render already fits 62.5 ms, so no cadence speedup is required. Pack and submit nest inside wake and are not summed with it. DMA-completion and other interrupts are not measured by these rows.

## Summary ranking

1. `grd_render`: 51.484 ms/f, 82.5% of the frame.
2. `grd_shader_draw`: 31.112 ms/f, 49.8% of the frame.
3. `grd_physics`: 5.212 ms/f, 8.3% of the frame.

Native Release checks preserve the full, clipped and shipped render goldens. The exhaustive production traversal matches all 46,080 ordered neighbor slots against the independent full table. Hoisting the hill-climb run lookup fails with 153 changed cubemap texels; swapping two production slots fails exactly two ordered-traversal assertions. Debug and Release reaction-graph/BZ/GS/guard checks and 242 Python firmware-gate tests pass.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` is parented by its first caller and is absent from this direct-path capture; no complete blend count is inferred. MIXED-PARENT palette/seed counters, when present, are non-exclusive and must not be summed with their displayed siblings.
- No per-pixel scopes were added. Existing scope overhead is part of the measured image.
- Shipping selective-O3 physics/render/Scan regions remain enabled; the full neighbor array is excluded from production, and compact decoding preserves neighbor order. Node positions and visual quality settings are unchanged.
- The 150-second profiling epoch is a lifetime stretch only; this 130-second pass does not cross a reset. The older canonical report covered a shorter segment, so its peak is not the matched baseline for this change.
- The measured candidate worktree was clean. Source/ELF/compiler attestations and immutable profile/shipping artifacts are named in the raw log.
- The supported locked wrapper uses isolated uploader helpers from `C:/work/temp/teensy-tools-loopback-20261030` to correct false localhost-listener discovery. Installed vendor binaries and measured firmware are unchanged by that host workaround; hashes and exact helper edits are retained in `loopback-preflight.json`.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET` and `HS_PROFILE_WINDOW`; `just profile` invokes the supported build/flash/capture wrapper. Reproduce this extended steady-effect pass with the command in Setup, including its epoch stretch. Use `parse_profile.py` modes `validate` and `windows` on the retained raw log.
