# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopemandala_ship.log`, captured 2026-10-22 19:54 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopemandala_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 105 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeMandala profile 105 32` |

Image size (`profile`, this effect only): `FLASH: code:88528, data:155308, headers:9092   free for files:1778688` / `RAM1: variables:313024, code:16872, padding:15896   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopemandala_ship_aa61e787859b_6819471f1a67`. Profile ELF SHA-256: `6819471f1a6748b9e97dd4c783b1244ce451f9186b9d222639a800ff6b73e531`.

## Frame cadence

**Pass aggregate**: peak frame render is **39.963 ms** (frame 1030), mean live-frame render is **32.902 ms**, and **0/1631** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 28.324 ms/f across complete post-setup windows; its worst window is 32.201 ms/f (frames 1025–1056).

Setup frame 1 rendered 54.019 ms. Published live-frame figures use complete counter windows through frame 1632; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:45) recorded peak 🟢 34.97 (2) and spilled 🟢 0/2367 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 1025–1056)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.362 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 32.201 | 1.00 | 51.64% |
| ↳ ↳ `fx_prepare_column` | 1.940 | 144.00 | 3.11% |
| ↳ `fx_prepare_frame` | 4.083 | 1.00 | 6.55% |
| ↳ `fx_advance` | 1.974 | 1.00 | 3.16% |
| ↳ `fx_timeline_step` | 0.195 | 1.00 | 0.31% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 23.322 | 1.00 | 37.40% |

Wall min/avg/max = 60.248/62.374/64.284 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 39.963 | 0/1033 | 32.201 | 32/32 |
| 1 | — | 33.140 | 0/598 | 28.831 | 17/18 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999907 us
isr_wake       n=36864 cyc min/avg/max=298/1215/17666 ns min/avg/max=496/2025/29443 total=74710 us cpu=3.73%
isr_pack       n=4608 cyc min/avg/max=5280/5686/7204 ns min/avg/max=8800/9476/12006 total=43675 us cpu=2.18%
isr_dma_submit n=4608 cyc min/avg/max=393/565/4111 ns min/avg/max=655/941/6851 total=4346 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 32.201 ms/f, 51.64% of the peak window.
2. `fx_prepare_frame` — 4.083 ms/f, 6.55% of the peak window.
3. `fx_advance` — 1.974 ms/f, 3.16% of the peak window.
4. `fx_timeline_step` — 0.195 ms/f, 0.31% of the peak window.
5. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.

README cells: peak 🟢 39.963 (2), spilled 🟢 0/1631 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopeMandala  windows=51
  [PASS] captured >=3 windows (51)
  [PASS] every window names one effect and resolution (KaleidoscopeMandala 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 51 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 51 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (198.4 ppm, frames 1025-1056)

INVALID: KaleidoscopeMandala
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (198.4 ppm, frames 1025-1056). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
