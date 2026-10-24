# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-10-23 17:21 America/Los_Angeles on COM4.
Replaces `profile_shapeshifter_teensy_2026-10-23.md` (the earlier `49bb199f4` capture, a different source and frame scope).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter, 288×144, single-entry playlist, source `bda58b65d3fc317f5d6ccbe29888014738a0ee25` |
| Method | 269000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ShapeShifter profile 269 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:192856, data:155844, headers:8676   free for files:1674240` / `RAM1: variables:312992, code:53832, padding:11704   free for local variables:145760` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995632, data:851824, headers:9056   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `9aa090642d5524640a3774b19f90becf33126ff420cc5c760450daa19080032f`. Phantasm ELF SHA-256: `5fc8212484e6c1fb01b5305253408e170a8226b46fb0cbbaff78393b0a20b80f`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Capture provenance:

```text
profile provenance: version=1
profile provenance: source_sha=bda58b65d3fc317f5d6ccbe29888014738a0ee25
profile provenance: compiler=GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203
profile provenance: profile_elf_sha256=9aa090642d5524640a3774b19f90becf33126ff420cc5c760450daa19080032f
profile provenance: phantasm_elf_sha256=5fc8212484e6c1fb01b5305253408e170a8226b46fb0cbbaff78393b0a20b80f
profile provenance: profile_envdump_sha256=14b68df161aa3ab3b04113714664f4b211a61c141c4c56039efc03b1fb23c993
profile provenance: phantasm_envdump_sha256=f5b172fbaa6ad2b8d5af681df338fa70985682fe91967d36b794d3b3f973c9ea
profile provenance: artifact_profile_elf=build/prof/artifacts/shapeshifter_ship_bda58b65d3fc_9aa090642d55/profile.elf
profile provenance: artifact_phantasm_elf=build/prof/artifacts/shapeshifter_ship_bda58b65d3fc_9aa090642d55/phantasm.elf
profile provenance: artifact_dir=build/prof/artifacts/shapeshifter_ship_bda58b65d3fc_9aa090642d55
profile provenance: profile_image=built:f06da468b817c4da3e7dbe26e5b732cd
profile provenance: phantasm_image=cache:de16925e8d887db981c5aec1e331c698
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **57.537 ms** (frame 410), mean **35.048 ms**, **0/4271 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–4272 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **47.865 ms**, 0/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 4280 raw frame rows, including that trailing portion, the observed peak is 57.537 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `ss_draw_all` averages 34.667 ms/frame; its worst window averages 53.365 ms/frame (frames 3553–3568).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.668 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.128 | 0.94 | 0.21% |
| ↳ `ss_draw_all` | 34.082 | 1.00 | 57.12% |
| ↳ ↳ `ss_plot_dispatch` | 33.474 | 281.19 | 56.10% |
| ↳ `ss_timeline_step` | 0.073 | 1.00 | 0.12% |
| ↳ `ss_buffer_wait` | 25.372 | 1.00 | 42.52% |
| ↳ ↳ `canvas_clear` | 0.094 | 1.00 | 0.16% |
| ↳ ↳ `canvas_buffer_wait` | 25.272 | 1.00 | 42.35% |

Wall min/mean/max: 34.341/59.672/66.077 ms.

### Peak-render window (frames 401–416)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.809 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 1.00 | 0.22% |
| ↳ `ss_draw_all` | 44.833 | 1.00 | 71.38% |
| ↳ ↳ `ss_plot_dispatch` | 44.245 | 214.94 | 70.44% |
| ↳ `ss_timeline_step` | 0.078 | 1.00 | 0.12% |
| ↳ `ss_buffer_wait` | 17.746 | 1.00 | 28.25% |
| ↳ ↳ `canvas_clear` | 0.085 | 1.00 | 0.14% |
| ↳ ↳ `canvas_buffer_wait` | 17.656 | 1.00 | 28.11% |

Wall min/mean/max: 39.022/62.812/85.768 ms.

### Worst dominant-scope window (frames 3553–3568)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.249 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `ss_draw_all` | 53.365 | 1.00 | 85.73% |
| ↳ `ss_timeline_step` | 0.096 | 1.00 | 0.15% |
| ↳ `ss_buffer_wait` | 8.621 | 1.00 | 13.85% |
| ↳ ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ ↳ `canvas_buffer_wait` | 8.508 | 1.00 | 13.67% |

Wall min/mean/max: 57.029/62.252/67.366 ms.

### Ownership-change window (frames 3345–3360)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.631 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `ss_draw_all` | 45.745 | 1.00 | 73.04% |
| ↳ `ss_timeline_step` | 0.111 | 1.00 | 0.18% |
| ↳ `ss_buffer_wait` | 16.604 | 1.00 | 26.51% |
| ↳ ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ ↳ `canvas_buffer_wait` | 16.492 | 1.00 | 26.33% |

Wall min/mean/max: 53.668/62.634/71.528 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `ss_draw_all` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 57.537 | 0/470 | 47.783 | 29/30 |
| 8 | — | 56.329 | 0/479 | 53.365 | 29/29 |
| 7 | — | 53.384 | 0/479 | 50.832 | 29/30 |
| 6 | — | 44.793 | 0/479 | 43.268 | 29/30 |
| 4 | — | 42.417 | 0/479 | 40.802 | 29/30 |
| 5 | — | 37.098 | 0/479 | 31.785 | 29/30 |
| 9 | — | 35.807 | 0/448 | 34.444 | 28/28 |
| 2 | — | 29.667 | 0/479 | 25.233 | 29/30 |
| 3 | — | 29.660 | 0/479 | 25.061 | 29/30 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 1005.808 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1158.75 | 0.000568 | 0.001931 | 0.029736 | 3.56% |
| `isr_pack` | 144.88 | 0.006903 | 0.008803 | 0.011281 | 2.02% |
| `isr_dma_submit` | 144.88 | 0.000751 | 0.000948 | 0.001253 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `ss_draw_all` — 44.833 ms/frame, 71.38% of this root window.
2. `pov_preserve_half` — 0.140 ms/frame, 0.22% of this root window.
3. `ss_timeline_step` — 0.078 ms/frame, 0.12% of this root window.

README cells: peak 🟢 57.537 (9), spilled 🟢 0/4271 (0.00%).

README setup frame 1: 47.865 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

Previous `49bb199f4` capture (different source and capture): steady peak 70.224 ms, setup frame 1 47.913 ms, steady mean 39.203 ms, 134/3999 steady frames spilled. The still earlier `a6222f23b` capture peaked at 56.941 ms with 0/4287 spilled frames under its recorded scope.

The ShapeShifter regression came from flash image placement of the uniform column walker. Device-clip host counters were unchanged by restoring the older sweep-arc bound: two poses, equal sweep arc and sweep rate, and at sampled frame 3554 the same 14,463 columns, 262 prepared edges, 53 culled edges and 22,183 emitted samples. The walker retained the same 4,100-byte body; aligning it to a flash cache line restored device timing without changing its arithmetic or the correctness bounds. Shipping ITCM remains 171,616 bytes.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=ShapeShifter  windows=267
  [PASS] captured >=3 windows (267)
  [PASS] every window names one effect and resolution (ShapeShifter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 267 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 267 windows are short)
  [PASS] scope 'ss_draw_all' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (52.7 ppm, frames 3553-3568)

INVALID: ShapeShifter
```

Validation limitations: cycle returns to its first index (peak idx 9, 8 distinct); all 9 presets visited (8); root cyc/600 vs wall sum within 5 ppm (52.7 ppm, frames 3553-3568). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
