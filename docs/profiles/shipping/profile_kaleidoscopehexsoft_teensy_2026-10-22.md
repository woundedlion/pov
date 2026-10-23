# KaleidoscopeHexSoft on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopehexsoft_ship.log`, captured 2026-10-22 19:50 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopehexsoft_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexSoft 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeHexSoft profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:86296, data:154656, headers:8904   free for files:1781760` / `RAM1: variables:313024, code:17064, padding:15704   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopehexsoft_ship_aa61e787859b_817ad00a2735`. Profile ELF SHA-256: `817ad00a27352a5ca7fc5bfc87859577448b4fbac173278eacb8c621df89a8a1`.

## Frame cadence

**Pass aggregate**: peak frame render is **27.376 ms** (frame 223), mean live-frame render is **26.995 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 20.399 ms/f across complete post-setup windows; its worst window is 20.503 ms/f (frames 193–224).

Setup frame 1 rendered 42.224 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:37) recorded peak 🟢 23.88 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 193–224)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.377 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.147 | 1.00 | 0.24% |
| ↳ `fx_shader_draw` | 20.503 | 1.00 | 32.87% |
| ↳ ↳ `fx_prepare_column` | 0.876 | 144.00 | 1.41% |
| ↳ `fx_prepare_frame` | 3.658 | 1.00 | 5.86% |
| ↳ `fx_advance` | 2.127 | 1.00 | 3.41% |
| ↳ `fx_timeline_step` | 0.145 | 1.00 | 0.23% |
| ↳ `canvas_clear` | 0.096 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 35.362 | 1.00 | 56.69% |

Wall min/avg/max = 61.816/62.392/62.708 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999980 us
isr_wake       n=36866 cyc min/avg/max=295/1207/15213 ns min/avg/max=491/2011/25355 total=74217 us cpu=3.71%
isr_pack       n=4608 cyc min/avg/max=5340/5756/10626 ns min/avg/max=8900/9593/17710 total=44211 us cpu=2.21%
isr_dma_submit n=4608 cyc min/avg/max=391/559/2858 ns min/avg/max=651/931/4763 total=4300 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 20.503 ms/f, 32.87% of the peak window.
2. `fx_prepare_frame` — 3.658 ms/f, 5.86% of the peak window.
3. `fx_advance` — 2.127 ms/f, 3.41% of the peak window.
4. `pov_preserve_half` — 0.147 ms/f, 0.24% of the peak window.
5. `fx_timeline_step` — 0.145 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.096 ms/f, 0.15% of the peak window.

README cells: peak 🟢 27.376, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopeHexSoft  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (KaleidoscopeHexSoft 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (210.4 ppm, frames 1-32)

INVALID: KaleidoscopeHexSoft
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (210.4 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexSoft`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
