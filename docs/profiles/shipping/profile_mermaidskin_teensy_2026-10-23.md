# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MermaidSkin`).
Raw capture: `C:/work/temp/Holosphere-msopt/build/prof/final_tip_mermaid.log`, captured 2026-10-22 23:18
America/Los_Angeles on COM4. Replaces the October 22 shipping report
with a newer-source scoped measurement; the former capture’s coverage is
not inherited.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` pullback hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, source `4e75923c44d52467b3d6383740e6987c457ddd17` |
| Method | 30 s capture, 32-frame windows; first 32 frames excluded from the matched table. Scoped post-warmup check of the resident-group and stationary-grid implementation; PC sampling is disabled. This 30-second capture is not a whole-show claim. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-msopt HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MermaidSkin profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:88888, data:154520, headers:8496   free for files:1779712` / `RAM1: variables:313024, code:25256, padding:7512   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:916640, data:851664, headers:8336   free for files:254976` / `RAM1: variables:312736, code:169336, padding:27272   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
The approved placement ledger records `.text.itcm` 160048 to 169328 B,
an increase of **9280 B (9.063 KiB)**, with 27280 B below the 196608 B cap.
RAM1 code grows 160056 to 169336 B, variables remain 312736 B, FLASH code
grows 915824 to 916640 B, and FLASH data remains 851664 B. These sizes belong
to the captured source; later engine landings have their own image sizes.

Profile ELF SHA-256: `bba90c05915e99734dbb9dc8653d3dfc3dc28b193aab6e72c3875b02ebdb9d31`.
The matching `.provenance`, envdump and build logs retain source/compiler
attestation and artifact paths alongside the raw log. Displacement, gradient, palette and grid hot-path placements accompany the stationary-grid path and resident-group correction. Shader-entry pins are absent.

Exactness cross-check: parser-selected frames 1–32 have
root cycles/600 = 1984118.638 us against
1984272 us measured wall sum (**77.3 ppm**). This
exceeds the strict 5 ppm limit; direct per-frame render telemetry supports
the timing figures, but this is not a fully cycle-validated capture.

## Frame cadence

**Scoped pass aggregate:** peak frame render **42.116 ms**
(frame 254), mean **41.638 ms**, **0/416 spills**
over complete-window frames 33–448.
`fx_shader_draw` averages 34.936 ms/frame; its worst
retained window is 35.097 ms/frame
(frames 65–96).

Setup frame 1 renders **69.943 ms** and exceeds the
display budget. Frames 1–32 are warmup, not included
in the headline zero-spill count. 8 trailing raw frame rows
without a complete counter window are excluded. MermaidSkin emits no preset markers. This is a scoped 30-second steady-effect capture, with no full-show duration guarantee.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. All retained live renders fit the
16 fps render budget; wall readings include sync and scheduling jitter.
`canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: this capture is a short selected regime after warmup;
unobserved show phases are not represented. Counter trees are inclusive
per-frame averages and mirror the logged nesting.

### Retained regime, worst dominant-scope window (frames 65–96)

```text
frame                       62.411 ms 37.447 Mcyc 100.0%
  pov_preserve_half          143.2 us   85.9 kcyc  0.2%
  fx_shader_draw            35.097 ms 21.058 Mcyc 56.2%
    fx_prepare_column        1.314 ms  788.5 kcyc  2.1%
  fx_prepare_frame           4.079 ms  2.448 Mcyc  6.5%
  fx_advance                 1.908 ms  1.145 Mcyc  3.1%
  fx_timeline_step           102.5 us   61.5 kcyc  0.2%
  canvas_clear                93.5 us   56.1 kcyc  0.1%
  canvas_buffer_wait        20.669 ms 12.402 Mcyc 33.1%
```

Wall min/avg/max = 62.140/62.416/62.969 ms.
The draw kernel dominates retained foreground cost. Buffer wait accounts for synchronization idle; it is excluded from render peak and mean figures.

### Per-pixel figures

Blended-pixel counts are unavailable: `filter_blend` and scan-count telemetry
are absent. Cycles per blend and scan cycles per blended pixel cannot be
derived; the 10368-pixel quadrant is only the geometric footprint.

## Column-ISR / DMA marshaling cost

Worst retained scope-window rates per rendered frame and per-call
min/average/maximum times:

```text
isr_wake         1152.0/f  0.558/2.028/28.850 us  3.73%
isr_pack          144.0/f  8.891/9.558/27.048 us  2.20%
isr_dma_submit    144.0/f  0.625/0.941/8.606 us  0.21%
```

- Packing costs 2.20% and submission 0.21%
  of the independently timestamped ISR interval. Wire transfer runs
  asynchronously; its duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.73%, leaving
  60.169 ms after that interrupt load per
  62.5 ms window. Render scopes already include interruptions. Pack and
  submit run inside wake and must not be added to its share.
- Retained frames require no speedup for the display budget. DMA-completion
  and other interrupts are unmeasured.

## Summary ranking

1. `fx_shader_draw` — 56.2% of the window, 35.097 ms/frame.
2. `fx_prepare_frame` — 6.5% of the window, 4.079 ms/frame.
3. `fx_advance` — 3.1% of the window, 1.908 ms/frame.
4. `fx_prepare_column` — 2.1% of the window, 1.314 ms/frame.

The approved matched table records peak 61.164 to 42.116 ms, a 19.048 ms reduction, and mean 59.723 to 41.638 ms, an 18.085 ms reduction. These are frames 33–448 of the existing matched captures; no new baseline was captured for this refresh.
Native/WASM rates do not establish on-device cadence. The index row uses
this scoped peak and identifies the limited scope; it does not establish
unobserved preset or full-show performance.

README cells: peak 🟢 42.116 (scoped), spilled 🟢 0/416 (0.00%).

## Caveats

- All cycle scopes include ISR time because CYCCNT free-runs.
- Instrumentation overhead remains. No per-pixel scope breakdown was used;
  per-pixel instrumentation would change the measured workload.
- `filter_blend` is absent. When enabled, first-entry parenting can hide its
  subtree in windows where the recorded parent has no calls.
- These are shipping selective-O3 measurements of the pullback paths,
  including the landed shared placements and stationary-grid implementation.
  Global-O3 results remain historical in the separate report set.
- There is no dwell-compression or epoch-stretch flag in these envdumps;
  the optional fixed-preset flag and observed markers are stated above.
- Captured source is `4e75923c44d52467b3d6383740e6987c457ddd17`. No newer engine state is
  substituted for this provenance; this refresh only parses existing logs.
- The strict validator fails root-cycle/wall agreement. Identity, frame
  monotonicity and complete-window telemetry pass. Its fixed-preset check
  is not a full-cycle check; zero spills applies only to retained frames.

Parser validation:

```text
effect=MermaidSkin  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (MermaidSkin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (77.3 ppm, frames 1-32)

INVALID: MermaidSkin
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MermaidSkin`,
`HS_PROFILE_WINDOW=32`; `just profile MermaidSkin` invokes the
standard build/flash/capture wrapper. The explicit Reproduce command names
this capture’s duration, window and selection flags under the device lock.
