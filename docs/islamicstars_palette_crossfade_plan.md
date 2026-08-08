# IslamicStars: per-leg fresh palette classification with smooth crossfade

**Status: LANDED 2026-07-26 (60733687 + fix 56c102e8 — the structural dual
handoff no longer requires an immutable handoff, unblocking the smooth dual
bridge's closing leg).** Measured: worst leg carries 19/36
distinct blend pairs; scratch_b peak unchanged at 70,228 / 73,728 B; phantasm
ITCM 1,400 B free after the cohort-plumbing deletion. Fragment path unchanged.
The section-2 hazard was resolved by deriving MAX_BLEND_PAIRS = N^2 (bank is
N = 6 after Blood Stream's removal); the leg-start pair-cap fallback was never
needed.

Goal: replace the birth-colour ("stable face colour") model in the IslamicStars
build tour with per-leg fresh palette assignment — each leg's arrival mesh is
classified on its own and mapped to a freshly shuffled palette set — and a
smooth per-face crossfade from the departed palette to the target over the leg.
Hard constraint: no performance regression.

## 1. Shipped behavior

Each leg classifies its arrival mesh and derives a fresh target palette. The
handoff carries the landed palette between legs, and distinct source/target
pairs are pre-blended into scratch ramps once per frame. The deleted immutable,
pinned-target and birth-counter modes are not configuration surfaces.

## 2. Blend-pair bound

`OpLeg::MAX_BLEND_PAIRS` is derived as `PALETTES * PALETTES`, covering all 36
pairs for the six-palette bank. The measured roster maximum is 19 pairs, so the
leg-start fallback considered below was not needed.

Scratch budget: each blended pair costs one 256-entry Color4 LUT = 3 KB per
frame in scratch_arena_b. Device splits (`IslamicStars.h:163-167`) give
scratch_b 72-75 KB, with a measured Needle peak of 70,228 B — but that peak is
at generation/classification time, not during swept-frame draw. Worst case
36 x 3 KB = 108 KB does **not** fit; measured-realistic counts (reconcile
precedent: 18 pairs = 54 KB) likely do.

## 3. How it landed

**Measurement first (host).** A host-only counter in `intern_palette_ramp`
replayed the full recipe roster (test_conway_morph.h build replay, test_solids
sweep paths) under a prototype non-immutable flip, logging the per-leg
distinct-pair maximum and the swept-frame scratch_b high-water (ramps plus
concurrent draw content). The maximum came in at 19 pairs, comfortably inside
headroom, so `MAX_BLEND_PAIRS` took the `PALETTES * PALETTES` ceiling and the
fallback leg-start pair cap — sort pairs by face population, keep the top B,
merge the tail onto its `to` endpoint — was never built. A per-fragment
two-LUT lerp was rejected outright: it regresses the memory-bound shade path.

**The behaviour flip**, all in `effects/IslamicStars.h`. The handoff sites
`seed_handoff`, `landing_handoff` and dual-bridge leg 3 stopped pinning a
target, freezing the palette, and keying newborns off a birth counter. Legs
carry the **landed** palette rather than `from_palette` (valid only when
to == from), through a `Landing::landed_palette(f)` helper reading
`to_palette[wrap(topology[f], PALETTES)]` — the HankinSolids pattern
(`effects/HankinSolids.h`). The standing-display adoption reads the same
helper, so post-leg standing frames match the w = 1 frame exactly with no pop.
Newborn faces mid-leg key by `wrap(class, PALETTES)` into the leg's fresh
`to_palette`. The dead birth-cohort plumbing went with the flip; the
engine-side immutable machinery (`PaletteHandoff::immutable` / `birth_counter`
/ `pinned_to` and everything they gated) has since been deleted from
`core/animation/opleg.h`, as no shipped effect, target, or tool ever selected
it. Spawn colouring is untouched — birth colours of the spawned mesh are still
class + shuffle.

**Tests.** The intent-pinned expectations moved once: build replay
(test_conway_morph.h), high-water sweeps (test_solids.h), and the arena survey
(test_opchain_arena_survey.h, where swept-frame scratch_b rises). The per-leg
shuffle consumes extra Pcg32 draws, so the downstream stream shifted
deterministically. The measurement probe became a permanent CI-wired pin:
roster-wide max distinct pairs <= `MAX_BLEND_PAIRS`, guarding against the
HS_CHECK trap.

**Perf and size.** The pre-blend upper bound is 36 pairs x 256 RampBlend
samples (~2 LUT lerps each), and only mid-leg frames pay it — endpoint plateaus
alias for free. The before/after on-device IslamicStars profile on identical
builds stayed binary green: 0 spilled frames, peak column green. All new code
is HS_COLD_MEMBER / cold-path and the cohort-plumbing deletion offsets it,
leaving 1,400 B of phantasm ITCM free under `pio run -e phantasm` (the only
build that catches the granule cliff) and the CI size gate.

**Visual QA.** Framebuffer dumps at w ~ 0.5 show the blend, the arrival frame
and the first standing frame are identical, and leg boundaries show no pop.
The daydream/wasm build was checked after install.

## 4. Configuration surface

No halfway immutable or pinned-target mode was retained. Changing the palette
policy now requires a new handoff design rather than setting the deleted
`immutable`, `pinned_to` or `birth_counter` fields.
