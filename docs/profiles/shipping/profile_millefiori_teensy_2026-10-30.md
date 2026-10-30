# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Millefiori`) of clean source `50e659aa535548d98e4d33258164aa5c0147f254`.
Raw capture: `C:/work/temp/millefiori-study-20261030/direct_fields_50e659aa5.log`. This replaces the
2026-10-29 Millefiori report.

**The user’s ≤58 ms ongoing target is met over all 2,047 fully reported ongoing
frames: 57.488 ms peak at frame 1621, 0 frames above 58 ms, and 0/2,047 display
spills.** Frame 1 is setup: 81.095 ms and one display spill. The complete pass
therefore has an 81.095 ms peak and 1/2,048 spills. Startup is excluded only from
the explicitly ongoing target; it remains included in full-pass reporting.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live; wrapper selected COM3 |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` functions, `-ffast-math -fno-finite-math-only`; no PC or deep profiling |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, 480 RPM |
| Effect | Millefiori 288×144, single-entry playlist, clean source `50e659aa5` |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture, `HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

`profile`: FLASH code 104,168, data 153,376, headers 8,696; RAM1 variables 313,056, code 32,200, padding 568, free local 178,464; RAM2 variables 520,064, free 4,224 bytes.

`phantasm`: FLASH code 1,037,536, data 749,244, headers 8,292; RAM1 variables 312,736, code 196,600, padding 8, free local 14,944; RAM2 variables 520,064, free 4,224 bytes.

Cold final `phantasm` and `holosphere` builds compiled all ten first-party
translation units across both environments, with zero first-party warnings. Shipping free local space is 14,944 B,
above the 12,288 B shipping floor. Only **8 B** remains below the 196,608 B
instruction-memory boundary. Linker free space is distinct from measured stack
high-water usage.

Exactness cross-check: canonical `parse_profile.py` reports **VALID**, with root
cycles ÷ 600 matching wall sum within **1.1 ppm**, frames 1953–1984, selected
scope `fx_shader_draw`. All 64 windows have complete per-frame telemetry and
render distinct from wall. Frames 2049–2056 have no completed window/footer and
are excluded; the final accepted window ends at frame 2048.

### Frozen provenance

| Artifact | SHA-256 / value |
|---|---|
| Source | `50e659aa535548d98e4d33258164aa5c0147f254` |
| Profile ELF | `3675721a068c291d70a712d996812b1ddd54091b1d188c0e790ebeeb04d7bc1d` |
| Shipping ELF | `c4a081fd136b54826be65ed5f21c0b8dcbc2189041129b1551266db7ce73ea18` |
| Raw capture | `32ac1f38f4d1316e86450fea3fc9ca5504f891a7ed7136875b59cb23b3282b8a` |
| Profile envdump | `17ac25fc4ddb93d59e36a1de42550c15e4b542da95d5e665d807d3bcf4bae55a` |
| Compiler | GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203 |
| Artifact state | clean; empty source status/diff; PC=false; deep=false |

Artifact directory: `C:/work/temp/millefiori-study-20261030/artifacts/direct_fields_50e659aa5_50e659aa5355_3675721a068c`.
Both ELF hashes and the profile envdump hash were checked against archived files;
sidecar and logged provenance agree. Machine-readable evidence is retained in
`C:/work/temp/millefiori-study-20261030/profile_millefiori_50e659aa5_data.json`
and `C:/work/temp/millefiori-study-20261030/compare_ongoing_report.json`.

## Frame cadence

| Scope | Frames | Peak render | >58 ms | Display spills (>62.5 ms) |
|---|---:|---:|---:|---:|
| Setup | 1 | 81.095 ms | 1 | 1/1 |
| Ongoing acceptance | 2–2048 | **57.488 ms**, frame 1621 | **0** | **0/2,047** |
| Full pass | 1–2048 | **81.095 ms**, frame 1 | 1 | **1/2,048 (0.049%)** |

README cells: peak 🟢 57.488, spilled 🟢 0/2047 (0.00%).
These cells use ongoing frames; setup is listed in its own README column.

Ongoing mean render is 50.325 ms. Full-pass mean
render is 50.340 ms. `fx_shader_draw` averages 47.905 ms/f;
its highest window mean is 52.897 ms/f
(frames 1953–1984). These means do not replace
the exact maximum used for acceptance. The ongoing peak leaves 0.512 ms against
the user’s 58 ms limit and 5.012 ms against the 62.5 ms display window.

The nominal display quadrant is 144×72 = 10,368 pixels. Owned clip margins make
the ongoing shader workload 146×73 = **10,658 pixels/frame** (1.028× nominal).
The prime draw uses full width, 288×73 = **21,024 pixels**. Both paths prepare
144 sweep-time contexts/frame; opposite arms share preparation. Source call
chain: `prime_segment_frame` sets the owned y band and full x clip before its
initial draw; steady segmented draws set the half-width clip with margins.
`canvas_buffer_wait` is display synchronization idle, removed from render.

## Phase-by-phase readout

Millefiori animates continuously; these are observed workload windows, not
separate presets. Counter trees show per-frame averages and root-relative cycle
percentages. Parents are inclusive; do not sum them with their children.

### Startup-containing window (frames 1–32)

```text
frame                      62.501 ms 37.500 Mcyc 100.0%
  pov_preserve_half         0.135 ms  0.081 Mcyc   0.2% x0.96875
  fx_shader_draw           45.349 ms 27.209 Mcyc  72.6% x1
    fx_prepare_column       1.858 ms  1.115 Mcyc   3.0% x144
  fx_prepare_frame          0.623 ms  0.374 Mcyc   1.0% x1
  fx_advance                1.157 ms  0.694 Mcyc   1.9% x1
    palette_cycle_step      1.060 ms  0.636 Mcyc   1.7% x1
  fx_timeline_step          0.097 ms  0.058 Mcyc   0.2% x1
  canvas_clear              0.093 ms  0.056 Mcyc   0.1% x1
  canvas_buffer_wait       14.663 ms  8.798 Mcyc  23.5% x1
```

Wall min/avg/max = 42.255/62.500/81.098 ms.
Mean render = 47.838 ms; the exact peak render in this window is
81.095 ms. Column preparation runs 144 times/frame,
averaging 12.903 µs/call. These wall figures include synchronization idle;
acceptance uses the per-frame render rows.

### Ongoing peak window (frames 1601–1632)

```text
frame                      62.449 ms 37.469 Mcyc 100.0%
  pov_preserve_half         0.138 ms  0.083 Mcyc   0.2% x1
  fx_shader_draw           52.751 ms 31.651 Mcyc  84.5% x1
    fx_prepare_column       3.180 ms  1.908 Mcyc   5.1% x144
  fx_prepare_frame          0.846 ms  0.507 Mcyc   1.4% x1
  fx_advance                1.114 ms  0.668 Mcyc   1.8% x1
    palette_cycle_step      1.087 ms  0.652 Mcyc   1.7% x1
  fx_timeline_step          0.095 ms  0.057 Mcyc   0.2% x1
  canvas_clear              0.089 ms  0.054 Mcyc   0.1% x1
  canvas_buffer_wait        7.038 ms  4.223 Mcyc  11.3% x1
```

Wall min/avg/max = 58.274/62.448/66.719 ms.
Mean render = 55.411 ms; the exact peak render in this window is
57.488 ms. Column preparation runs 144 times/frame,
averaging 22.081 µs/call. These wall figures include synchronization idle;
acceptance uses the per-frame render rows.

### Highest mean shader cost (frames 1953–1984)

```text
frame                      62.439 ms 37.463 Mcyc 100.0%
  pov_preserve_half         0.137 ms  0.082 Mcyc   0.2% x1
  fx_shader_draw           52.897 ms 31.738 Mcyc  84.7% x1
    fx_prepare_column       2.761 ms  1.657 Mcyc   4.4% x144
  fx_prepare_frame          0.827 ms  0.496 Mcyc   1.3% x1
  fx_advance                0.961 ms  0.577 Mcyc   1.5% x1
    palette_cycle_step      0.932 ms  0.559 Mcyc   1.5% x1
  fx_timeline_step          0.094 ms  0.056 Mcyc   0.2% x1
  canvas_clear              0.090 ms  0.054 Mcyc   0.1% x1
  canvas_buffer_wait        7.051 ms  4.230 Mcyc  11.3% x1
```

Wall min/avg/max = 60.978/62.438/63.822 ms.
Mean render = 55.388 ms; the exact peak render in this window is
56.172 ms. Column preparation runs 144 times/frame,
averaging 19.172 µs/call. These wall figures include synchronization idle;
acceptance uses the per-frame render rows.

### Lowest mean shader cost (frames 737–768)

```text
frame                      62.311 ms 37.387 Mcyc 100.0%
  pov_preserve_half         0.139 ms  0.083 Mcyc   0.2% x1
  fx_shader_draw           41.648 ms 24.989 Mcyc  66.8% x1
    fx_prepare_column       2.330 ms  1.398 Mcyc   3.7% x144
  fx_prepare_frame          0.577 ms  0.346 Mcyc   0.9% x1
  fx_advance                1.050 ms  0.630 Mcyc   1.7% x1
    palette_cycle_step      1.023 ms  0.614 Mcyc   1.6% x1
  fx_timeline_step          0.107 ms  0.064 Mcyc   0.2% x1
  canvas_clear              0.088 ms  0.053 Mcyc   0.1% x1
  canvas_buffer_wait       18.325 ms 10.995 Mcyc  29.4% x1
```

Wall min/avg/max = 59.131/62.311/65.402 ms.
Mean render = 43.987 ms; the exact peak render in this window is
46.985 ms. Column preparation runs 144 times/frame,
averaging 16.182 µs/call. These wall figures include synchronization idle;
acceptance uses the per-frame render rows.

### Per-pixel figures

The direct shader path has no `filter_blend` subtree. Figures divide the inclusive
shader counter by actual clipped output pixels, including ISR interruptions and
shared column preparation. They are not isolated sample-kernel costs.

| Window | Actual pixels/window | Shader cycles/pixel | Shader µs/pixel |
|---|---:|---:|---:|
| 1–32 | 351,422 | 2477.6 | 4.129 |
| 1601–1632 | 341,056 | 2969.7 | 4.949 |
| 1953–1984 | 341,056 | 2977.9 | 4.963 |
| 737–768 | 341,056 | 2344.6 | 3.908 |

## Column-ISR / DMA marshaling cost

Rates and CPU shares use each ISR reporting interval, which is a separate clock
from the frame wall sum. Per-call min/avg/max below come from the firmware’s
cycle-derived nanosecond telemetry. `isr_wake` includes nested pack and submit;
their CPU shares must not be added to wake.

```text
isr_wake
  isr_pack
  isr_dma_submit
```

### Ongoing peak window ISR (frames 1601–1632, interval 2,001,279 µs)

| Scope | Calls/frame | Calls/s | Min/avg/max µs/call | CPU share |
|---|---:|---:|---:|---:|
| `isr_wake` | 1152.81 | 18433.2 | 0.543/2.018/19.076 | 3.72% |
| `isr_pack` | 144.09 | 2304.0 | 8.756/9.503/12.305 | 2.18% |
| `isr_dma_submit` | 144.09 | 2304.0 | 0.611/0.945/6.851 | 0.21% |

Inclusive wake costs 3.721% CPU, equivalent to 2.326 ms
of a 62.5 ms display window. Render scopes already absorb interruptions; this
is not an additional deduction from the measured render result.

### Highest mean shader cost ISR (frames 1953–1984, interval 2,000,856 µs)

| Scope | Calls/frame | Calls/s | Min/avg/max µs/call | CPU share |
|---|---:|---:|---:|---:|
| `isr_wake` | 1152.56 | 18433.1 | 0.543/2.015/25.456 | 3.71% |
| `isr_pack` | 144.06 | 2304.0 | 8.671/9.460/12.066 | 2.17% |
| `isr_dma_submit` | 144.06 | 2304.0 | 0.623/0.940/1.091 | 0.21% |

Inclusive wake costs 3.716% CPU, equivalent to 2.322 ms
of a 62.5 ms display window. Render scopes already absorb interruptions; this
is not an additional deduction from the measured render result.

Pack is the larger nested marshaling cost; submit measures CPU register/setup
work. Async LED wire-transfer time and DMA completion are not timed here. Other
interrupts may also occur inside the render scopes.

## Summary ranking

- `fx_shader_draw`: 47.905 ms/frame, 76.75% of full-pass root time.
- `fx_prepare_column`: 2.570 ms/frame, 4.12% of full-pass root time (nested in shader).
- `fx_advance`: 1.009 ms/frame, 1.62% of full-pass root time.
- `fx_prepare_frame`: 0.722 ms/frame, 1.16% of full-pass root time.
- `canvas_buffer_wait`: 12.074 ms/frame, 19.35% of full-pass root time (display idle).

The largest measured region remains the shader. The clean equal-coverage
precursor `1921ee80b` peaks at 61.908 ms, frame 1621, over frames 2–2048, with
412 frames above 58 ms and 0 display spills. `50e659aa5` reduces that peak by
**4.420 ms (7.14%)** to 57.488 ms and removes all >58 ms frames. Ongoing mean
falls from 54.746 to
50.325 ms. This is a same-coverage comparison;
older shorter captures are retained in the external comparison JSON and do not
replace any candidate’s ALL-ongoing peak.

### Diagnostic PC sampling and direct carrier writes

The separate precursor PC run `palette_lab_pc_1921ee80b.log` recorded 620,600
ticks over 62,060 ms, including 553,232 render samples. Its report attributes
22,376 samples (4.0% of render samples) to `memcpy`, 8.6% to the generated
palette path, and 27.2% to the inlined Millefiori shader. Statistical PC bins
are not exclusive phase timing or proof that all memcpy belongs to one caller.
The PC capture **fails strict checksum validation** (144.4 ppm, frames 193–224)
and is diagnostic only; it is excluded from target acceptance.

Source change `50e659aa5` writes cached palette carrier fields directly in
`GeneratedPalette::apply_shared_base` instead of assigning a temporary aggregate,
then stores cached Lab and its validity flag. It retains the field, opacity,
palette coordinate, hue/lift defaults, and brightness gain. The clean no-PC
before/after captures establish the net speed improvement; sampling explains
the investigation without claiming an isolated memcpy speedup.

ARM disassembly confirms that the cached-base helper shrinks from 784 to
752 bytes and its local stack allocation from 92 to 12 bytes. Its 36-byte
`memset` and 69-byte `memcpy` calls disappear. Every carrier member is set
explicitly; a dirty-output test catches an omitted lift reset. The final
rendering change is bit-exact over 600 native frames. The palette cache itself
uses two 64-by-3 OKLab buffers (4,608 bytes of node payload) in the existing
arena. Compared with the preceding phase-only renderer, its 600-frame RGB16
difference has maximum 3,229 and mean 3.12 channel codes; ripple geometry is
unchanged.

### Validation results

- User explicitly approved the radial-table ripple appearance: “I think it’s ok, yes.”
- Final cold `phantasm`/`holosphere` compile: all ten first-party translation
  units across both environments, zero first-party warnings; code/RAM gates pass.
- Native full CTest: **124/125 pass, 0 skips**; the sole `unit_case_calls`
  static check (`test_composed_effect.h:verify_export reaches no assertion`)
  fails identically on clean accepted `53f24ce02`. All 175 goldens pass;
  `HS_SMOKE_FRAMES=120`.
- Debug: **1,362,939 + 11 assertions pass**. Four relevant Release modules
  pass; the final palette write change is **bit-exact over 600 frames**.
- Final WASM 120-frame smoke: **48 failures**, unchanged versus `1921ee80b`
  and accepted `53f24ce02`: 46 control-mapping failures and two Millefiori
  stack high-water checks. Measured 2,324/8,192 B exceeds the existing 2,048 B
  creep tripwire. This is 16 B above the 2,308 B observation at `53f24ce02`;
  no tolerance was increased and this remains a known failure.
- Installed final WASM SHA-256:
  `e5c9167e6d48d73b9d89004f0a19a8493787411228930e06a614cdc4423fd24e`.
  JavaScript SHA-256:
  `e2185c91652fa21f0631281a5e0f7231584f5ec82f9e313479881c1c9be2ee9a`.
- Versus accepted `53f24ce02`, shipping RAM1 code grows 368 B; RAM1 variables
  and FLASH data are unchanged. Firmware gates pass.

## Caveats

- This is a finite 130 s capture of one continuously animated default effect,
  not a proof over every control setting or all future animation frames.
- Exact ongoing acceptance covers every fully reported frame 2–2048, including
  later frame 1621; neither startup nor a matched prefix hides a later peak.
- CYCCNT free-runs: render and shader scopes include ISR and instrumentation.
  Parent scopes and children overlap. No per-pixel scope was added.
- The shader writes pixels directly; `filter_blend` counters are absent.
- Shipping uses an `-Os` base with selected `HS_O3_FN` shader/scan paths and hot
  color/math helpers; this is not an all-translation-unit `-O3` measurement.
- `HS_PROFILE_EPOCH_REVS=1200` prevents playlist rollover during capture; it
  changes the capture lifetime, not the per-frame effect workload.
- Final artifact source is clean, with PC/deep profiling disabled. Dirty or PC
  diagnostic captures remain labeled separately in the comparison report.
- Color fidelity to the old framebuffer was relaxed by the user; radial
  geometry and visible motion were qualified and the user approved the look.
  Rendering retains one composed color finish per final pixel.
- Native and WASM validation retains the explicitly listed baseline failures;
  the measured timing target does not waive those failures.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`,
`HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`.

```sh
bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'
python tools/parse_profile.py C:/work/temp/millefiori-study-20261030/direct_fields_50e659aa5.log validate --scope fx_shader_draw
python tools/parse_profile.py C:/work/temp/millefiori-study-20261030/direct_fields_50e659aa5.log windows --scope fx_shader_draw
python C:/work/temp/millefiori-study-20261030/compare_ongoing.py
```

The bench capture used the supported wrapper through the portable TyTools
adapter at `C:/work/temp/millefiori-uploader-20261030/capture.sh`, with
`HS_PROFILE_TREE` set explicitly to the clean source checkout. Both board locks
remained held through upload and capture. No loader, driver, or execution-policy
setting was changed.
