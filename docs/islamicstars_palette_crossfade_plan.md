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

## 1. What already exists

The requested behaviour is `OpLeg`'s *default* colour model; IslamicStars opts
out of it via `PaletteHandoff::immutable = true`:

- `immutable = false` makes every carried face crossfade `from -> to` under the
  leg's `BlendWeightFn` (`core/animation/opleg.h:195`). IslamicStars already passes
  `classic_blend` at all three leg-construction sites
  (`effects/IslamicStars.h:1077,1187,1254`), so the whole-leg smooth lerp comes
  for free the moment `immutable` flips.
- `pinned_to = nullptr` gives the per-leg fresh shuffle of target palettes
  (`opleg.h:1771-1772`, via `hs::shuffle` — determinism preserved).
- Fresh classification per arrival is *already computed* — every leg site runs
  `classify_faces_by_topology` on the arrival seed and passes it as
  `BookendClasses`. Colour targets key on it (`opleg.h:1805`), so the w = 1
  frame equals the standing display by construction.
- The fragment path is identical in both modes: `ramps[face_ramp[f]].get(t)`,
  one baked-LUT lookup. Crossfade cost is a once-per-frame pre-blend of the
  distinct (from, to) pairs into `scratch_arena_b` (`opleg.h:1441-1455`), with a
  zero-cost aliasing fast path at w <= 0 / w >= 1
  (`composition.h:1093-1102`).

So per-fragment cost cannot regress; the entire perf question is the per-frame
pre-blend (time + scratch_b bytes) on mid-leg frames.

## 2. The one real hazard: MAX_BLEND_PAIRS is stale

`OpLeg::MAX_BLEND_PAIRS = 25` (`opleg.h:42`) with a comment claiming it covers
`BakedPaletteBank::N^2 = 25`. **N is now 6** (`composition.h:1129`), so the
true pair space is 36. `intern_palette_ramp` traps via `HS_CHECK` on overflow
(`opleg.h:1715`) — a fail-fast device crash. Today immutable legs stay tiny and
the worst measured leg (gyro_kis reconcile) hit 18; fresh-classified legs raise
the ceiling everywhere. This must be resolved by measurement + a bound, not
hope.

Scratch budget: each blended pair costs one 256-entry Color4 LUT = 3 KB per
frame in scratch_arena_b. Device splits (`IslamicStars.h:163-167`) give
scratch_b 72-75 KB, with a measured Needle peak of 70,228 B — but that peak is
at generation/classification time, not during swept-frame draw. Worst case
36 x 3 KB = 108 KB does **not** fit; measured-realistic counts (reconcile
precedent: 18 pairs = 54 KB) likely do.

## 3. Plan

### Phase A — measure before touching behaviour (host)

1. Instrument `intern_palette_ramp` (host-only counter) and replay the full
   recipe roster (test_conway_morph.h build replay / test_solids sweep paths)
   with a prototype non-immutable flip, logging per-leg distinct-pair max and
   the swept-frame scratch_b high-water (ramps + concurrent draw content).
2. Decision from the numbers:
   - max pairs x 3 KB fits swept-frame scratch_b headroom -> raise
     `MAX_BLEND_PAIRS` to that max (ceiling 36) and fix the stale comment.
   - it does not fit -> add a leg-start pair cap: sort pairs by face
     population, keep top B (budget from scratch_b headroom), merge tail pairs
     onto their `to` endpoint (those few faces snap at leg start —
     deterministic, bounded, visually minor). Per-fragment two-LUT lerp is
     rejected: it regresses the memory-bound shade path.

### Phase B — the behaviour flip

All in `effects/IslamicStars.h`:

1. Handoff sites `seed_handoff` (:989), `landing_handoff` (:1031), dual-bridge
   leg 3 (:1169): `pinned_to = nullptr`, `immutable = false`,
   `birth_counter = nullptr`.
2. Carry the **landed** palette between legs, not `from_palette` (valid only
   when to == from). Landed palette per face f is
   `to_palette[wrap(topology[f], PALETTES)]` — the HankinSolids pattern
   (`effects/HankinSolids.h:842`). Sites: `landing_handoff` pal fill (:1029),
   `carry_landing_to_seed` (:1211), leg-3 pal snapshot (:1140), and the
   standing-display adoption `slot_face_palette[front]` (:1493) so post-leg
   standing frames match the w = 1 frame exactly (no pop).
   Prefer adding a small `Landing::landed_palette(f)` helper over repeating the
   expression.
3. Delete the now-dead birth-cohort plumbing: `build_palette_order_`,
   `build_birth_counter_`, and the immutable-comment blocks (spawn colouring at
   :692-708 stays — birth colours of the spawned mesh are still class + shuffle).
   The engine-side immutable machinery (`PaletteHandoff::immutable` /
   `birth_counter` / `pinned_to` and everything they gated) has since been
   deleted from `core/animation/opleg.h`; no shipped effect, target, or tool
   ever selected it.
4. Newborn faces mid-leg: with a null birth counter the engine keys newborns by
   `wrap(class, PALETTES)` into the leg's fresh `to_palette` — accept, verify
   visually.

### Phase C — tests

- Update intent-pinned expectations: build replay (test_conway_morph.h),
  high-water sweeps (test_solids.h), arena survey
  (test_opchain_arena_survey.h — swept-frame scratch_b rises).
- New pin: roster-wide max distinct pairs <= MAX_BLEND_PAIRS (turns the
  Phase A probe into a permanent guard against the HS_CHECK trap). Wire into
  the CI runner.
- RNG note: the per-leg shuffle consumes extra Pcg32 draws, shifting the
  downstream stream — replay/golden expectations move once, deterministically.

### Phase D — perf + size gates

- Per-frame pre-blend time upper bound: 36 pairs x 256 RampBlend samples
  (~2 LUT lerps each) — sub-millisecond expected, but only mid-leg frames pay
  it (endpoint plateaus alias for free), and mid-leg frames are also the
  geometry-heavy ones. Gate: on-device profile (teensy-profile) of IslamicStars
  before/after on identical builds; ship criterion is the binary green —
  0 spilled frames, peak column green.
- ITCM: phantasm sits at ~40 B headroom. All new code is HS_COLD_MEMBER /
  cold-path; deletion of the cohort plumbing offsets. Verify with
  `pio run -e phantasm` (the only build that catches the granule cliff) and the
  CI size gate.

### Phase E — visual QA

- Framebuffer dump harness: mid-leg frames (w ~ 0.5) show the blend; the
  arrival frame vs. the first standing frame must be identical (carry fix
  check); leg boundaries show no pop.
- daydream/wasm check after install (hard-refresh; verify disk sha256 first).

## 4. Design knob kept open

If the per-leg re-shuffle reads as too busy (every leg re-rolls every colour),
the halfway point is `pinned_to = &build_palette_order_` with
`immutable = false`: fresh classification per leg, shared target set — faces
whose class survives keep their colour, only reclassified faces fade. Fewer
blend pairs too. Worth an A/B before committing the shuffle.
