# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/dreamballs_ship.log`. Captured 2026-11-04 23:43 America/Los_Angeles on COM4. Replaces profile_dreamballs_teensy_2026-11-01.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DreamBalls, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 200 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh DreamBalls profile 200 16` |

Image size:

```text
teensy_size:   FLASH: code:157968, data:190564, headers:8844   free for files:1674240
teensy_size:    RAM1: variables:314144, code:37048, padding:28488   free for local variables:144608
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `db8318abbb3358b6ddc233279aaeb1d78c8865b6c4457eb9618421241abb5ca3`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 2769–2784, root cyc/600 versus measured wall sum differs by **1.095 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   DreamBalls               worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 38,308 B  full-set worst 4
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **52.024 ms** (frame 2752), **0/3167 spilled** (0.00%). Setup frame 1: 52.869 ms. Complete windows cover frames 1–3168; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 28.867 ms/frame. `db_timeline_step` averages 28.599 ms/frame over complete windows; its worst window is 48.917 ms/frame (frames 2769–2784).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 0.802–52.024 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     60.512 ms 36.307 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.082 Mcyc   0.2%
    x0.94/frame; 145.223 us/call
  db_timeline_step                        39.791 ms 23.875 Mcyc  65.8%
    db_draw                               39.707 ms 23.824 Mcyc  65.6%
      db_draw_scene                       39.689 ms 23.813 Mcyc  65.6%
        db_mesh_plot                      39.153 ms 23.492 Mcyc  64.7%
          x18.00/frame; 2175.158 us/call
        db_displace                        0.466 ms  0.280 Mcyc   0.8%
          x18.00/frame; 25.883 us/call
  canvas_clear                             0.099 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.710 us/call
  canvas_buffer_wait                      20.448 ms 12.269 Mcyc  33.8%
    x1.00/frame; 20448.194 us/call
```

Wall min/avg/max: 37.558/60.512/66.666 ms. Window mean render is 40.064 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 2737–2752)

```text
frame                                     62.778 ms 37.667 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.412 us/call
  db_timeline_step                        48.541 ms 29.125 Mcyc  77.3%
    db_draw                               48.378 ms 29.027 Mcyc  77.1%
      db_draw_scene                       48.361 ms 29.017 Mcyc  77.0%
        db_mesh_plot                      47.001 ms 28.201 Mcyc  74.9%
          x6.00/frame; 7833.528 us/call
        db_displace                        1.118 ms  0.671 Mcyc   1.8%
          x6.00/frame; 186.326 us/call
  canvas_clear                             0.093 ms  0.056 Mcyc   0.1%
    x1.00/frame; 93.372 us/call
  canvas_buffer_wait                      13.960 ms  8.376 Mcyc  22.2%
    x1.00/frame; 13959.983 us/call
```

Wall min/avg/max: 59.746/62.778/65.199 ms. Window mean render is 48.818 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 2769–2784)

```text
frame                                     62.393 ms 37.436 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.172 us/call
  db_timeline_step                        48.917 ms 29.350 Mcyc  78.4%
    db_draw                               48.752 ms 29.251 Mcyc  78.1%
      db_draw_scene                       48.735 ms 29.241 Mcyc  78.1%
        db_mesh_plot                      47.364 ms 28.418 Mcyc  75.9%
          x6.00/frame; 7894.008 us/call
        db_displace                        1.124 ms  0.675 Mcyc   1.8%
          x6.00/frame; 187.396 us/call
  canvas_clear                             0.096 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.693 us/call
  canvas_buffer_wait                      13.193 ms  7.916 Mcyc  21.1%
    x1.00/frame; 13193.121 us/call
```

Wall min/avg/max: 59.198/62.393/65.423 ms. Window mean render is 49.200 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 2865–2880)

```text
frame                                     59.565 ms 35.739 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.971 us/call
  db_timeline_step                        44.057 ms 26.434 Mcyc  74.0%
    db_draw                               43.895 ms 26.337 Mcyc  73.7%
      db_draw_scene                       43.878 ms 26.327 Mcyc  73.7%
        db_mesh_plot                      42.601 ms 25.561 Mcyc  71.5%
          x5.62/frame; 7573.510 us/call
        db_displace                        1.053 ms  0.632 Mcyc   1.8%
          x5.62/frame; 187.260 us/call
  canvas_clear                             0.091 ms  0.055 Mcyc   0.2%
    x1.00/frame; 91.435 us/call
  canvas_buffer_wait                      15.230 ms  9.138 Mcyc  25.6%
    x1.00/frame; 15230.324 us/call
```

Wall min/avg/max: 16.938/59.565/66.485 ms. Window mean render is 44.335 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 9  | 20/19 | unavailable | 48.917 | 49.200 | 16.03 | modal-call-count hold |
| 1  | 20/19 | unavailable | 45.431 | 45.711 | 15.97 | modal-call-count hold |
| 10  | 18/18 | unavailable | 36.451 | 36.729 | 16.00 | modal-call-count hold |
| 8  | 20/19 | unavailable | 35.964 | 36.248 | 16.02 | modal-call-count hold |
| 4  | 20/19 | unavailable | 28.129 | 28.408 | 16.03 | modal-call-count hold |
| 3  | 20/19 | unavailable | 26.175 | 26.455 | 16.01 | modal-call-count hold |
| 2  | 20/19 | unavailable | 25.274 | 25.550 | 16.04 | modal-call-count hold |
| 7  | 20/19 | unavailable | 23.864 | 24.142 | 15.97 | modal-call-count hold |
| 6  | 20/19 | unavailable | 17.492 | 17.772 | 16.06 | modal-call-count hold |
| 5  | 20/19 | unavailable | 16.936 | 17.211 | 16.00 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 9 | 52.024 | 0/320 | 16.00 |
| 1 | 48.520 | 0/318 | 16.00 |
| 10 | 40.414 | 0/289 | 16.00 |
| 8 | 38.804 | 0/320 | 16.00 |
| 4 | 32.220 | 0/320 | 16.00 |
| 3 | 28.918 | 0/320 | 16.00 |
| 2 | 28.707 | 0/320 | 16.00 |
| 7 | 25.907 | 0/320 | 16.00 |
| 6 | 20.302 | 0/320 | 16.00 |
| 5 | 18.518 | 0/320 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1153.06/frame
  min/avg/max 0.525/1.793/19.673 us; cpu 3.30%
  isr_pack         144.12/frame
    min/avg/max 6.500/7.688/10.733 us; cpu 1.77%
  isr_dma_submit   144.12/frame
    min/avg/max 0.761/0.935/1.378 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1000836 us
isr_wake       n=18449 cyc min/avg/max=315/1076/11804 ns min/avg/max=525/1793/19673 total=33088 us cpu=3.30%
isr_pack       n=2306 cyc min/avg/max=3900/4613/6440 ns min/avg/max=6500/7688/10733 total=17730 us cpu=1.77%
isr_dma_submit n=2306 cyc min/avg/max=457/561/827 ns min/avg/max=761/935/1378 total=2159 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.30% leaves approximately 60.438 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `db_timeline_step` — 48.541 ms/frame, 77.32% of this root window.
2. `db_draw` — 48.378 ms/frame, 77.06% of this root window.
3. `db_draw_scene` — 48.361 ms/frame, 77.04% of this root window.
4. `db_mesh_plot` — 47.001 ms/frame, 74.87% of this root window.
5. `db_displace` — 1.118 ms/frame, 1.78% of this root window.

README cells: peak 🟢 52.024 (10), spilled 🟢 0/3167 (0.00%).

README setup frame 1: 52.869 ms. Captured 2026-11-04 23:43 America/Los_Angeles.

This capture measures the generated hot-code layout against the 11-04 matched sweep; see the sweep scorecard for the delta.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=DreamBalls  windows=198
  [PASS] captured >=3 windows (198)
  [PASS] every window names one effect and resolution (DreamBalls 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 10, 9 distinct)
  [FAIL] all 10 presets visited (9)
  [PASS] per-frame render is render, not wall (0 of 198 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 198 windows are short)
  [PASS] scope 'db_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.1 ppm, frames 2769-2784)

INVALID: DreamBalls
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 10, 9 distinct); [FAIL] all 10 presets visited (9). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`. `just profile DreamBalls` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
