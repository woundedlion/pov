# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexOil`).
Raw capture: `C:/work/temp/Holosphere-msopt/build/prof/final_tip_hex.log`, captured 2026-10-22 23:19
America/Los_Angeles on COM4. Replaces the October 22 shipping report
with a newer-source scoped measurement; the former capture’s coverage is
not inherited.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` pullback hot paths; PC sampling disabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, source `4e75923c44d52467b3d6383740e6987c457ddd17` |
| Method | 30 s capture, 32-frame windows; first 32 frames excluded from the matched table. Scoped post-warmup check with fixed profile selection 0/2; the effect immediately adopts and logs Preset 2/2. PC sampling is disabled. No subsequent preset advance or wrap occurs. |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-msopt HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeHexOil profile 30 32 '-D HS_PROFILE_PRESET=0'` |

Image size (`profile`, this effect only): `FLASH: code:88280, data:158900, headers:8820   free for files:1775616` / `RAM1: variables:313024, code:22728, padding:10040   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image: `FLASH: code:916640, data:851664, headers:8336   free for files:254976` / `RAM1: variables:312736, code:169336, padding:27272   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.
The approved placement ledger records `.text.itcm` 160048 to 169328 B,
an increase of **9280 B (9.063 KiB)**, with 27280 B below the 196608 B cap.
RAM1 code grows 160056 to 169336 B, variables remain 312736 B, FLASH code
grows 915824 to 916640 B, and FLASH data remains 851664 B. These sizes belong
to the captured source; later engine landings have their own image sizes.

Profile ELF SHA-256: `8913ba703a843eba0e87c24e56e46e7b4dbb50bb77130c0c1f456abae7a215ab`.
The matching `.provenance`, envdump and build logs retain source/compiler
attestation and artifact paths alongside the raw log. This shader crosses the shared displacement, gradient, palette and grid hot paths, including the stationary-grid specialization. The profile retains the landed selective placements; shader-entry pins are absent.

Exactness cross-check: parser-selected frames 385–416 have
root cycles/600 = 1997037.338 us against
1997454 us measured wall sum (**208.6 ppm**). This
exceeds the strict 5 ppm limit; direct per-frame render telemetry supports
the timing figures, but this is not a fully cycle-validated capture.

## Frame cadence

**Scoped pass aggregate:** peak frame render **43.790 ms**
(frame 447), mean **41.873 ms**, **0/416 spills**
over complete-window frames 33–448.
`fx_shader_draw` averages 38.151 ms/frame; its worst
retained window is 39.474 ms/frame
(frames 385–416).

Setup frame 1 renders **69.802 ms** and exceeds the
display budget. Frames 1–32 are warmup, not included
in the headline zero-spill count. 8 trailing raw frame rows
without a complete counter window are excluded. Only the setup marker Profile preset 0/2 and the adopted Preset 2/2 appear. No full two-preset cycle or wrap is demonstrated; the prior October 22 full-segment report's coverage cannot be transferred to this source.

A display window is 62.5 ms at 480 RPM. The effect renders one quadrant of
the 288×144 canvas, about 10368 pixels. All retained live renders fit the
16 fps render budget; wall readings include sync and scheduling jitter.
`canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: this capture is a short selected regime after warmup;
unobserved show phases are not represented. Counter trees are inclusive
per-frame averages and mirror the logged nesting.

### Retained regime, worst dominant-scope window (frames 385–416)

```text
frame                       62.407 ms 37.444 Mcyc 100.0%
  pov_preserve_half          139.3 us   83.6 kcyc  0.2%
  fx_shader_draw            39.474 ms 23.684 Mcyc 63.3%
    fx_prepare_column        2.323 ms  1.394 Mcyc  3.7%
  fx_prepare_frame           1.097 ms  658.1 kcyc  1.8%
  fx_advance                 2.186 ms  1.312 Mcyc  3.5%
  fx_timeline_step           217.3 us  130.4 kcyc  0.3%
  canvas_clear                96.8 us   58.1 kcyc  0.2%
  canvas_buffer_wait        19.071 ms 11.442 Mcyc 30.6%
```

Wall min/avg/max = 62.174/62.420/62.651 ms.
The draw kernel dominates retained foreground cost. Buffer wait accounts for synchronization idle; it is excluded from render peak and mean figures.

### Per-preset table

Only the setup marker Profile preset 0/2 and the adopted Preset 2/2 appear. No full two-preset cycle or wrap is demonstrated; the prior October 22 full-segment report's coverage cannot be transferred to this source.
One scoped ownership entry is represented; no independent clean-hold or
full-cycle guarantee is made. Geometry and blend counts are not logged.

| Observed entry | Blended px/f | `fx_shader_draw` ms/f | Render ms/f | Wall fps | Windows | Peak ms | Spilled/frames |
|---|--:|--:|--:|--:|--:|--:|--:|
| Adopted 2/2; setup selection 0/2 | n/a | 39.474 | 43.350 | 16.02 | 13 | 43.790 | 0/416 |

### Per-pixel figures

Blended-pixel counts are unavailable: `filter_blend` and scan-count telemetry
are absent. Cycles per blend and scan cycles per blended pixel cannot be
derived; the 10368-pixel quadrant is only the geometric footprint.

## Column-ISR / DMA marshaling cost

Worst retained scope-window rates per rendered frame and per-call
min/average/maximum times:

```text
isr_wake         1152.3/f  0.495/1.975/29.950 us  3.64%
isr_pack          144.0/f  8.401/9.348/12.165 us  2.15%
isr_dma_submit    144.0/f  0.696/0.936/7.091 us  0.21%
```

- Packing costs 2.15% and submission 0.21%
  of the independently timestamped ISR interval. Wire transfer runs
  asynchronously; its duration is not measured by submit cost.
- Inclusive `isr_wake` costs 3.64%, leaving
  60.225 ms after that interrupt load per
  62.5 ms window. Render scopes already include interruptions. Pack and
  submit run inside wake and must not be added to its share.
- Retained frames require no speedup for the display budget. DMA-completion
  and other interrupts are unmeasured.

## Summary ranking

1. `fx_shader_draw` — 63.3% of the window, 39.474 ms/frame.
2. `fx_prepare_column` — 3.7% of the window, 2.323 ms/frame.
3. `fx_advance` — 3.5% of the window, 2.186 ms/frame.
4. `fx_prepare_frame` — 1.8% of the window, 1.097 ms/frame.

The approved matched table records peak 50.312 to 43.790 ms, a 6.522 ms reduction, and mean 48.103 to 41.873 ms, a 6.230 ms reduction. These are frames 33–448 of the existing matched captures; no new baseline was captured for this refresh.
Native/WASM rates do not establish on-device cadence. The index row uses
this scoped peak and identifies the limited scope; it does not establish
unobserved preset or full-show performance.

README cells: peak 🟢 43.790 (scoped), spilled 🟢 0/416 (0.00%).

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
effect=KaleidoscopeHexOil  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (KaleidoscopeHexOil 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] fixed profile preset 0/2
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (208.6 ppm, frames 385-416)

INVALID: KaleidoscopeHexOil
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexOil`,
`HS_PROFILE_WINDOW=32`; `just profile KaleidoscopeHexOil` invokes the
standard build/flash/capture wrapper. The explicit Reproduce command names
this capture’s duration, window and selection flags under the device lock.
