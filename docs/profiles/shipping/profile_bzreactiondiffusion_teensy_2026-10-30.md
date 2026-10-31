# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile BZReactionDiffusion`).
Raw capture: `build/prof/bzreactiondiffusion_ship.log`, copied with original mtime from `C:/work/temp/flash-20261030-bz-after.log`. This replaces `profile_bzreactiondiffusion_teensy_2026-10-23.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; flywheel and DMA ISRs live |
| Image | Shipping selective-O3 `profile`, effect render/physics and Scan hot regions active; no deep, replay, PC or count instrumentation |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | BZReactionDiffusion 288x144, single-entry playlist, source `480cd8ac1b417fcbe381ac645f4b5761268ab8e5` |
| Method | 130 s, window 32; `HS_PROFILE_EPOCH_REVS=1200` stretches the instance to 150 s, avoiding an epoch crossing |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-flash-neighbors-20261030 HS_TEENSY_PORT=COM4 HS_TEENSY_TOOLS=C:/work/temp/teensy-tools-loopback-20261030 bash tools/profile_one.sh BZReactionDiffusion profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Profile image size: `FLASH: code:73736, data:156512, headers:8344   free for files:1793024 / RAM1: variables:313024, code:16200, padding:16568   free for local variables:178496 / RAM2: variables:520064  free for malloc/new:4224`.

The shipping image uses FLASH code 1,021,896 B, data 749,104 B and headers 8,712 B, total 1,779,712 B. Against parent `93def4f2a`, Flash falls 91,136 B: data -91,880 B, code +272 B, headers +472 B. RAM1 code changes 187,160 to 187,192 B (+32 B), with padding shrinking by 32 B; bank allocation and 14,944 B stack headroom are unchanged. RAM1 variables remain 312,736 B and RAM2 variables remain 520,064 B.

Exactness cross-check: frames 577-608 root cycles / 600 MHz = 1997880.350 us; measured wall sum 1997878 us, delta 1.18 ppm. `parse_profile.py ... validate` passes frame coverage and no-reset checks. Captured 2026-10-30 09:46, America/Los_Angeles.

The matched COM4 baseline `flash-20261030-bz-before.log` ran source `aaef299bd` with the same effect, flags, duration, driver and epoch stretch. The candidate includes the independently landed cold-registration change and pullback/Millefiori work; those do not change these effects' simulation equations or authored configuration. This comparison measures the complete candidate against the recorded baseline, not the old shorter canonical capture.

## Frame cadence

Complete-window pass aggregate (`parse_profile.py ... windows`): 76.549 ms peak, 1/2048 spilled. Setup frame 1 is reported separately in the ranked index.

| COM4 metric | Baseline | Candidate |
|---|--:|--:|
| Setup render ms | 78.582 | 76.549 |
| Steady peak render ms | 48.319 | 44.608 |
| Steady spilled frames | 0/2047 | 0/2047 |

Steady peak changes -7.68%. A display window is 62.5 ms, giving 16 fps for steady frames; the active quadrant is 10,368 pixels. `canvas_buffer_wait` is the round-up idle until the next display flip. The peak steady frame has 17.892 ms margin. The startup frame already spills in the baseline and remains the sole spill; no steady frame spills.

All raw telemetry, including the trailing incomplete counter window, contains 2056 frames: peak 76.549 ms and 1 spills. The ranked index uses complete counter windows, frames 2-2048, matching its declared archive convention.

README cells: peak 🟢 44.608, spilled 🟢 0/2047 (0.00%).

## Phase-by-phase readout

The simulation evolves continuously. Startup, peak post-setup load and the late captured state sample its observed cost regimes; no preset-cycling marker is expected. Epoch stretching changes instance lifetime, not the per-frame physics or shading work.

### Startup (frames 1-32)

```
frame                    62.310 ms  37.386 Mcyc 100.0%
  pov_preserve_half      137.32 us   82.39 kcyc   0.2% x0.97 142 us/call
  bz_render              44.265 ms  26.559 Mcyc  71.0%
    bz_raster            39.979 ms  23.988 Mcyc  64.2% x2 19990 us/call
    bz_physics            4.189 ms   2.513 Mcyc   6.7% x2 2094 us/call
    rd_lattice_basis      94.06 us   56.44 kcyc   0.2% x1 94 us/call
  rd_timeline_step        47.78 us   28.67 kcyc   0.1% x1 48 us/call
  canvas_clear            91.25 us   54.75 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     17.757 ms  10.654 Mcyc  28.5% x1 17757 us/call
```

Wall min/avg/max = 42.293/62.310/76.550 ms. The first draw includes startup work; the other frames in this window already follow the normal display cadence. The individual first-frame value is separated from this window mean.

### Peak post-setup load (frames 577-608)

```
frame                    62.434 ms  37.460 Mcyc 100.0%
  pov_preserve_half      140.95 us   84.57 kcyc   0.2% x1 141 us/call
  bz_render              43.970 ms  26.382 Mcyc  70.4%
    bz_raster            39.721 ms  23.833 Mcyc  63.6% x2 19861 us/call
    bz_physics            4.155 ms   2.493 Mcyc   6.7% x2 2078 us/call
    rd_lattice_basis      92.07 us   55.24 kcyc   0.1% x1 92 us/call
  rd_timeline_step        39.24 us   23.54 kcyc   0.1% x1 39 us/call
  canvas_clear            87.35 us   52.41 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     18.186 ms  10.912 Mcyc  29.1% x1 18186 us/call
```

Wall min/avg/max = 61.618/62.433/63.150 ms. This window contains the largest post-setup render frame. Its mean render and buffer wait are distinct from the individual-frame peak used for cadence.

### Late simulation (frames 2017-2048)

```
frame                    62.434 ms  37.461 Mcyc 100.0%
  pov_preserve_half      144.94 us   86.96 kcyc   0.2% x1 145 us/call
  bz_render              43.875 ms  26.325 Mcyc  70.3%
    bz_raster            39.618 ms  23.771 Mcyc  63.5% x2 19809 us/call
    bz_physics            4.163 ms   2.498 Mcyc   6.7% x2 2081 us/call
    rd_lattice_basis      91.60 us   54.96 kcyc   0.1% x1 92 us/call
  rd_timeline_step        48.28 us   28.97 kcyc   0.1% x1 48 us/call
  canvas_clear            86.76 us   52.06 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     18.269 ms  10.962 Mcyc  29.3% x1 18269 us/call
```

Wall min/avg/max = 61.821/62.434/62.939 ms. This is the last complete counter window of the 130-second pass. It extends the observation beyond the previous shorter report and retains the changing simulation workload.

### Per-pixel figures

The quadrant holds 10,368 output pixels. This image does not emit a complete blended-pixel counter for the direct shader path, so blended coverage and cycles per blend cannot be inferred from draw-call counts. Counter instrumentation overhead remains included.

## Column-ISR / DMA marshaling cost

Representative peak post-setup window; rates, min/avg/max CPU times and measured ISR-window CPU shares:

```
isr_wake          1152.1/f 0.485/1.990/14.825 us 3.66%
  isr_pack         144.0/f 8.976/9.551/12.208 us 2.20%
  isr_dma_submit   144.0/f 0.628/0.931/1.160 us 0.21%
```

- DMA submission is CPU setup; pixel packing is separately measured inside wake.
- At 72 LEDs the image buffer is 300 B and its strobed composite is 600 B. At the configured SPI divider, `dma::transfer_us` budgets 115/230 us including per-byte gaps, outside submit CPU time.
- Inclusive wake CPU share is 3.66%, leaving about 60.213 ms foreground capacity per display window. The measured inclusive steady render already fits 62.5 ms, so no cadence speedup is required. Pack and submit nest inside wake and are not summed with it. DMA-completion and other interrupts are not measured by these rows.

## Summary ranking

1. `bz_render`: 43.970 ms/f, 70.4% of the frame.
2. `bz_raster`: 39.721 ms/f, 63.6% of the frame.
3. `bz_physics`: 4.155 ms/f, 6.7% of the frame.

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
