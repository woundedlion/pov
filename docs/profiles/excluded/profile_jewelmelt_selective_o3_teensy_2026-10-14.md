# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-10-14, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile JewelMelt`).
Raw capture: `build/prof/jewelmelt_ship.log`, captured 2026-10-14 20:16 PDT on COM4.
Replaces `profile_jewelmelt_teensy_2026-10-08.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os`, newlib-nano; shipping `HS_O3` shader scan, composed shading, surface displacement and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | JewelMelt 288×144, softness 0.5, single-entry playlist, clean source `fe0b439cb846241dd28476ec956a08f8e5cd9e91` |
| Method | `HS_PROFILE`, 32-frame windows, 70 s capture; no cadence overrides |
| Reproduce | `bash tools/profile_one.sh JewelMelt profile 70 32` with `HS_PROFILE_TREE` set explicitly |

Image size: `FLASH: code:69908, data:153360, headers:9176` / `RAM1: variables:312992, code:19208, padding:13560, free:178528` / `RAM2: variables:520064, free:4224`.

Full Phantasm passed memory/layout gates: FLASH code 656,004 B, FLASH data 853,176 B, RAM1 code 179,512 B, RAM1 variables 312,736 B, padding 17,096 B and locals free 14,944 B. No firmware source was changed for this capture.

Exactness cross-check: frames 865–896, root 1,198,550,614 cycles ÷ 600 = 1,997,584.357 us versus wall sum 1,997,588 us: **1.8 ppm**. Parser validation passed all 34 complete windows, with no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `fx_shader_draw` averages 37.48 ms/f, worst window 39.61 ms/f (865–896); peak render 77.45 ms, spilled 1/1088 frames, including initialization.

**Display-frame aggregate** (parser rows 2–1088): peak **46.35 ms** at frame 854; **0/1087 frames (0.00%) spilled**. The 62.5 ms window gives 16 fps and 16.15 ms peak render margin. Nine trailing rows outside complete counter windows are excluded.

Setup frame 1 renders the full-width 288×72 segment band before display publication, taking 77.45 ms. Live frames render one 144×72 quadrant (10,368 pixels). The first window has 31 preservation calls for 32 draws. Setup is retained in the pass aggregate and excluded from ranked display figures, matching the archive convention. `canvas_buffer_wait` is round-up idle until the next display flip.

## Phase-by-phase readout

Phase schedule: initialization followed by one continuously animated preset.

### Startup (window frames 1–32)

```
frame                     62.30 ms  37.38 Mcyc 100.0%
  pov_preserve_half       134.5 us   80.7 kcyc   0.2% x1.0 138.9 us/call
  fx_shader_draw          37.64 ms  22.58 Mcyc  60.4% x1.0 37641.4 us/call
  fx_prepare_frame         3.98 ms   2.39 Mcyc   6.4% x1.0 3976.7 us/call
  fx_advance               1.97 ms   1.18 Mcyc   3.2% x1.0 1969.8 us/call
  fx_timeline_step         67.6 us   40.5 kcyc   0.1% x1.0 67.6 us/call
  canvas_clear             91.0 us   54.6 kcyc   0.1% x1.0 91.0 us/call
  canvas_buffer_wait      18.38 ms  11.03 Mcyc  29.5% x1.0 18381.2 us/call
```

Wall min/avg/max = 44.43/62.30/77.45 ms. The wider setup draw accounts for the sole over-budget render. This mixed window does not represent steady-state cost.

### Steady peak shader cost (window frames 865–896)

```
frame                     62.42 ms  37.45 Mcyc 100.0%
  pov_preserve_half       140.0 us   84.0 kcyc   0.2% x1.0 140.0 us/call
  fx_shader_draw          39.61 ms  23.76 Mcyc  63.4% x1.0 39605.9 us/call
  fx_prepare_frame         3.88 ms   2.33 Mcyc   6.2% x1.0 3882.2 us/call
  fx_advance               2.02 ms   1.21 Mcyc   3.2% x1.0 2020.5 us/call
  fx_timeline_step         73.4 us   44.1 kcyc   0.1% x1.0 73.4 us/call
  canvas_clear             87.2 us   52.3 kcyc   0.1% x1.0 87.2 us/call
  canvas_buffer_wait      16.58 ms   9.95 Mcyc  26.6% x1.0 16577.5 us/call
```

Wall min/avg/max = 61.42/62.42/63.28 ms. Shader work dominates, with display synchronization filling the remaining interval. Every render meets the deadline; the highest individual display render occurs in the preceding window.

### Per-pixel figures

No `filter_blend` counter is emitted, so blend counts and cycles per blend are unavailable. The steady peak shader window costs 2292.0 cycles per quadrant pixel, using 10,368 pixels rather than measured blend calls.

## Column-ISR / DMA marshaling cost

Steady peak shader window; shares use the separately logged 1,999,856 us ISR interval, including report output.

```
isr_wake          1152.0/frame  0.55/1.54/17.42 us  cpu 2.83%
  isr_pack         144.0/frame  5.26/5.68/8.58 us   cpu 1.30%
  isr_dma_submit   144.0/frame  0.67/0.94/8.34 us   cpu 0.21%
```

- Submission costs 0.94 us/call versus 5.68 us/call for packing.
- At 24 MHz SPI, the 600-byte image-plus-black composite for 72 LEDs takes 200 us asynchronously; submission timing does not measure completion latency.
- Inclusive `isr_wake` consumes 2.83% CPU, leaving about 60.73 ms foreground CPU per 62.5 ms window. Pack and submit are nested inside wake, not added to it. No steady-state speedup is required. DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 63.4% of the steady peak window, 39.61 ms/f.
2. `fx_prepare_frame` — 6.2%, 3.88 ms/f.
3. `fx_advance` — 3.2%, 2.02 ms/f.
4. `pov_preserve_half` — 0.2%, 0.14 ms/f.

Display peak was 46.39 ms on October 8, with zero display spills in both captures. Different boards and source snapshots prevent attributing this small difference to a specific change. No same-workload WASM/native measurement accompanies this capture.

README cells: peak 🟢 46.35, spilled 🟢 0/1087 (0.00%).

## Caveats

- All scopes include ISR time because CYCCNT free-runs.
- `filter_blend` normally parents under its first caller and can be hidden when that parent is inactive; this capture emits no such counter.
- Scopes are per-frame; per-pixel scope overhead is not present.
- Shipping selective-O3 covers shader scan, composed shading, surface displacement and color evaluation; this is not global-O3.
- No dwell compression, ordered cycle or epoch stretch was used; capture stayed within the default 120 s epoch.
- The main checkout's unrelated uncommitted playlist edit failed Phantasm assertions. A clean detached checkout at the recorded revision supplied this capture; that edit was preserved and excluded.
- Compiler, flags, clean source status and ELF/map artifacts are retained under `build/prof/artifacts/jewelmelt_ship_fe0b439cb846_85d0a5f040ce/`. Profile ELF SHA-256: `85d0a5f040ce1aa3180a5a50e59349719e01f22c8b231daa080424898c70b94a`.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=JewelMelt`, `HS_PROFILE_WINDOW=32`; `just profile JewelMelt` builds, flashes and captures under the shared device lock.
