# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `build/prof/hyperlattice_ship.log` (standard cycle) and
`build/prof/hyperlattice_experimental_ship.log` (opt-in experimental cycle).
Replaces `profile_hyperlattice_teensy_2026-09-28.md`. Every per-change capture
of the optimization campaign is kept under
`build/prof/hyperlattice_campaign_2026-09-29/`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live; standard cycle on COM4, experimental cycle and pinned presets on COM3 (paired captures of one image on both boards agree within 0.06%) |
| Image | `profile` env; the traced paths run from cached flash (`HS_HOT_FLASH_MEMBER` scans) and cross no `HS_O3` region |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, landed tip `4d1bb7368` (captures from the pre-rebase branch; the rebase added only other effects' peer commits and a palette-rebake placement fix) |
| Method | `HS_PROFILE` cycle scopes, window 16. Standard cycle: 120 s, `-D HS_PROFILE_EPOCH_REVS=1200`. Experimental cycle: 345 s, `-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_EPOCH_REVS=2800`. Per-preset A/B: 30 s pinned captures, `-D HS_PROFILE_PRESET=<i>` |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 120 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size (standard profile image): `FLASH: code:77880, data:150016,
headers:8648` / `RAM1: variables:314944, code:13384, padding:19384, free:176576`
/ `RAM2: variables:520064, free:4224`.

Default Phantasm image against master `67f935538`: RAM1 code 168,856 →
168,936 B (+80 B, all inlining drift in other effects; HyperLattice's own ITCM
fell 240 B), FLASH code 512,264 → 529,344 B (+17,080 B). Experimental image:
RAM1 code 172,472 → 179,752 B (+7,280 B, again other effects' functions whose
inlining moved as the unit grew), FLASH code +88,720 B. RAM1 variables and the
12,896 B stack floor are unchanged in both.

Exactness cross-check: window frames 1185–1200 root counter cyc ÷ 600 MHz
matches the measured wall sum within **0.7 ppm**.

## Frame cadence

**Pass aggregate** (standard cycle): `hl_shader_draw` avg 23.86 ms/f, worst
window 28.09 ms/f (frames 1185–1200), peak frame render **37.46 ms** (frame
1198), spilled **0/1887**. Setup frame 1 (43.32 ms) is excluded.

A display window is 62.5 ms; the effect renders one quadrant, 146×73 = 10,658
rays per frame. Every frame of the standard cycle holds 16 fps with at least
25 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: three presets, each a 320-frame hold followed by a 240-frame
segue that now interpolates every continuous control; the pattern, view and
shell count switch at the segue midpoint.

### Segue into the hypercube (frames 1185–1200, worst window)

```
frame                  62.79 ms  37.67 Mcyc  100%
  hl_shader_draw       28.09 ms  16.86 Mcyc   45%  x1.0
  pov_preserve_half     0.14 ms   0.09 Mcyc    0%
  canvas_clear          0.09 ms   0.05 Mcyc    0%
  hl_timeline_step      6.4 us    3.9 kcyc     0%
  canvas_buffer_wait   32.68 ms  19.61 Mcyc   52%
```

Wall min/avg/max = 58.18/62.79/67.42 ms. The segue from Cubic Wide Flight into
Hypercube Flight carries the pass peak: its first half grows the cubic cell
toward the hypercube's while the far distance shrinks, and after the midpoint
the 4D slice starts with the lerped wider wires. Holds run 25–27 ms mean.

### Per-preset table

Pinned 30 s captures at the final code, frame 1 excluded, and the full
experimental cycle's buckets. A bucket opens at its preset's marker, so it
holds the segue into that preset and then its hold.

| # | Preset | Pinned peak ms | Pinned mean ms | Cycle bucket peak ms | Cycle spilled |
|---|---|--:|--:|--:|--:|
| 6 | Shell Flight | 27.88 | 23.55 | 🔴 65.49 | 9/559 |
| 5 | Octet 4D Flight | 56.12 | 47.00 | 🟢 49.24 | 0/559 |
| 8 | Shell 4D Flight | 51.17 | 35.04 | 🟢 51.00 | 0/559 |
| 3 | Octet Flight | 31.92 | 28.12 | 🟢 42.33 | 0/559 |
| 2 | Hypercube Flight | 34.36 | 26.90 | 🟢 37.34 | 0/559 |
| 7 | Shell Close Flight | 31.64 | 27.20 | 🟢 31.97 | 0/559 |
| 4 | Octet Wide Flight | 31.26 | 29.77 | 🟢 31.65 | 0/559 |
| 1 | Cubic Wide Flight | 29.15 | 25.59 | 🟢 29.58 | 0/681 |
| 0 | Cubic Flight | 26.37 | 24.72 | 🟢 43.30 | 0/878 |

The experimental cycle wrapped to preset 1 (nine `Preset:` markers). Cubic
Flight's bucket peak is setup frame 1. The only spilling frames are 3162–3234,
the first half of the segue from Octet 4D Flight into Shell Flight: the lerped
cell size shrinks toward the shell preset's, and each ray crosses about 36%
more octet planes before the pattern switches, 9 of 5,471 live frames.

### Per-pixel figures

The shader writes premultiplied pixels directly; there is no `filter_blend`
counter. Standard cycle: 23.86 ms/f mean over 10,658 rays is about 1,340
cycles per ray.

## Column-ISR / DMA marshaling cost

```
isr_wake         1158/frame  min/avg/max 0.6/1.6/11.0 us  cpu 2.99%
isr_pack          145/frame  min/avg/max 6.2/6.7/9.2 us   cpu 1.54%
isr_dma_submit    145/frame  min/avg/max 0.7/0.9/1.0 us   cpu 0.21%
```

- Packing marshals the LED data on the CPU; submission only launches the DMA.
- A 600-byte transfer takes about 400 µs at 12 MHz, asynchronously.
- The ISRs take 4.7% of the CPU, leaving about 59.5 ms of render budget per
  62.5 ms window: the standard cycle needs no further speedup; the Octet 4D
  segue needs about 9% on its worst frames.

## Summary ranking

1. `hl_shader_draw` — 45% of the worst window, 28.09 ms/f.
2. `pov_preserve_half` — 0.2%, 0.14 ms/f.
3. `canvas_clear` — 0.1%, 0.09 ms/f.

Against the previous shipping report the standard cycle's peak fell from 56.12
to 37.46 ms, and every experimental preset now holds 16 fps pinned.

## Optimization ledger

Each change was profiled on device before landing. Peak ms are pinned 30 s
captures unless noted; "ITCM" is the profile image's RAM1 code delta, zero
throughout except where listed (the traced paths all run from cached flash).

| Change | Presets | Peak ms before → after | ITCM |
|---|---|---|--:|
| Shell layer march (3D, single owner) | Shell / Shell Close | 123.5 → 35.0 / 126.0 → 40.3 | 0 |
| Dedicated cubic merge loop | Hypercube / Cubic Wide / Cubic | 57.7 → 46.5 / 48.9 → 39.6 / 45.5 → 36.2 | 0 |
| Per-axis crossings, sort covered only, one division | same | 46.5 → 40.2 / 39.6 → 34.4 / 36.2 → 30.9 | 0 |
| Inline cubic scan lambda (reverted) | Hypercube / Cubic | ±0 | +12,880 |
| Octet per-stream walks, sort covered only | Octet / Octet Wide | 46.1 → 41.3 / 45.0 → 39.6 | 0 |
| Float palette composite | Hypercube / Octet / Shell Close | 40.2 → 37.0 / 41.3 → 39.2 / 40.3 → 36.0 | 0 |
| Direct shell scan | Shell / Shell Close | 35.0 → 29.4 / 36.0 → 33.3 | 0 |
| Direct 3D octet scan | Octet / Octet Wide | 39.2 → 36.6 / 39.6 → 35.9 | 0 |
| Fixed two-shell slice scan | Hypercube | 37.0 → 35.5 | 0 |
| Fixed two-shell cubic scan (with the palette composite) | Cubic Wide / Cubic | 34.4 → 30.0 / 30.9 → 27.1 | 0 |
| Static octet owners in registers | Octet / Octet Wide | 36.6 → 32.6 / 35.9 → 31.6 | 0 |
| Canonical-frame 4D octet | Octet 4D (old slice params) | 137.7 → 63.5 | 0 |
| Tighter 4D bound plus lazy parity (reverted) | Octet 4D (old params) | 63.5 → 67.0 | 0 |
| 4D shell neighbor march | Shell 4D | 179.8 → 59.9 | 0 |
| Stateless 4D octet plane bound | Octet 4D | 74.5 → 72.7 | 0 |
| Paired crossings (reverted) | Octet 4D | 72.7 → 76.6 | 0 |
| 16-bit fixed-point plane bound | Octet 4D | 72.7 → 67.8 | 0 |
| Packed halfword adds (reverted) | Octet 4D | 67.8 → 69.3 | 0 |
| Counted crossings, fixed-point class bound | Octet 4D | 67.8 → 65.7 | 0 |
| Across-residual bound | Octet 4D | 65.7 → 63.8 | 0 |
| Two-stage bound (reverted) | Octet 4D | 63.8 → 66.3 | 0 |
| Integer forward-differenced threshold (reverted) | Octet 4D | 63.8 → 67.5 | 0 |
| Shared family setup, one division | Octet 4D | 63.8 → 62.7 | 0 |
| Integer magnitude ranking | Octet 4D | 62.7 → 62.0 | 0 |
| Counted shell layers | Shell 4D / Shell Close / Shell | 59.8 → 56.6 / 33.5 → 32.2 / 29.5 → 28.5 | 0 |
| Lazy class bounds (reverted) | Octet 4D | 62.0 → 64.1 | 0 |
| `fminf` threshold clamp | Octet 4D | 62.0 → 61.0 | 0 |
| `fmaxf` shell offset | Shell 4D | 56.6 → 50.6 | 0 |
| `fminf`/`fmaxf` sweep of traced paths | all | e.g. Hypercube 35.5 → 34.4, Octet 4D 61.0 → 60.0; Shell 4D 50.6 → 51.2 | 0 |
| Branch-free class mask and across shift | Octet 4D | 60.0 → 56.1 | 0 |
| 32-bit fixed-point fractions (reverted) | Octet 4D | 56.1 → 61.4 | 0 |
| 3D neighbor march for wide lerped shells | Shell Close → Shell 4D segue | 142.7 → 51.0 (cycle bucket) | 0 |

Octet 4D Flight and Shell 4D Flight joined the cycle with new parameters
during the campaign; each row measures the parameters current at its time.

## Caveats

- CYCCNT includes ISR time (4.7% of the CPU); it is part of every scope.
- No per-pixel scopes run in these captures; the campaign's deep captures
  (`HS_PROFILE_DEEP`) were attribution-only and are not timing figures.
- The traced paths run from cached flash; no `HS_O3` region is on them.
- The experimental presets require `HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`;
  epoch stretches only lengthen the effect instance, never a frame's cost.
- The pinned captures and both cycle captures ran the pre-rebase branch; the
  rebase onto master added peer commits outside HyperLattice.

## Harness

`targets/Profile/Profile.ino` knobs: `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS`, `HS_PROFILE_PRESET`, and
`HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1` for the experimental cycle. One-liner:
`just profile HyperLattice 120`.
