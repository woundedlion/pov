# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/aliencore_ship.log`, captured 2026-10-22 19:52 America/Los_Angeles on COM4.
Replaces `profile_aliencore_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienCore profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:86448, data:154488, headers:8920   free for files:1781760` / `RAM1: variables:313024, code:16792, padding:15976   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/aliencore_ship_aa61e787859b_e48233d0c565`. Profile ELF SHA-256: `e48233d0c5657bc6ccb6009d09794c8cf12b3059bfcac05c83cce073fd95cdce`.

## Frame cadence

**Pass aggregate**: peak frame render is **23.648 ms** (frame 314), mean live-frame render is **23.239 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 18.936 ms/f across complete post-setup windows; its worst window is 18.959 ms/f (frames 321–352).

Setup frame 1 rendered 39.672 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:42) recorded peak 🟢 17.86 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 289–320)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.396 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 18.950 | 1.00 | 30.37% |
| ↳ ↳ `fx_prepare_column` | 1.979 | 144.00 | 3.17% |
| ↳ `fx_prepare_frame` | 0.902 | 1.00 | 1.45% |
| ↳ `fx_advance` | 2.765 | 1.00 | 4.43% |
| ↳ `fx_timeline_step` | 0.106 | 1.00 | 0.17% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 39.032 | 1.00 | 62.55% |

Wall min/avg/max = 62.127/62.407/62.643 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999968 us
isr_wake       n=36866 cyc min/avg/max=304/1190/15682 ns min/avg/max=506/1983/26136 total=73152 us cpu=3.65%
isr_pack       n=4608 cyc min/avg/max=5344/5651/7330 ns min/avg/max=8906/9418/12216 total=43401 us cpu=2.17%
isr_dma_submit n=4608 cyc min/avg/max=409/562/1534 ns min/avg/max=681/936/2556 total=4318 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 18.950 ms/f, 30.37% of the peak window.
2. `fx_advance` — 2.765 ms/f, 4.43% of the peak window.
3. `fx_prepare_frame` — 0.902 ms/f, 1.45% of the peak window.
4. `pov_preserve_half` — 0.142 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.106 ms/f, 0.17% of the peak window.
6. `canvas_clear` — 0.093 ms/f, 0.15% of the peak window.

README cells: peak 🟢 23.648, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=AlienCore  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienCore 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (185.3 ppm, frames 1-32)

INVALID: AlienCore
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (185.3 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
