# Phantasm Synchronization Architecture — Design Spec

*Status: DRAFT for review. Not yet implemented. This is a full-stack redesign of
Phantasm's cross-board timing, of which code-review finding #1
(`pov_segmented.h:420-451`, the frame-sync/`advance_display` desync) is one
symptom at the frame layer. Decisions taken so far: (1) spec the hybrid fully
before implementing; (2) **collapse to a single inter-board wire** — delete the
shared column clock and generate columns from a local, time-derived flywheel on
every board; (3) flip ownership = reconcile + self-describing sync pulses (A+C).*

*Review pass folded in (feasibility + correctness): the architecture is feasible
and the §4.5 drift math holds. Five findings were incorporated rather than left
implicit — (i) the symbol must be **count/burst-coded, not width-coded**, so the
very `FastLED.show()` masking the flywheel defeats can't corrupt the one wire all
layers ride (§5.2, §11.3); (ii) "extra flip is benign" is a **Layer-2-only**
statement — content `t` is flip-paced through the `buffer_free()` gate, so a
spurious flip offsets `t` until the next epoch (§8.4, §6.1); (iii) the §9
**master-dead** row is graceful precession, not a watchdog trap, now that
downstream owns its columns; (iv) the epoch **absolute index** is promoted toward
baseline (trap, don't assume index 0) (§6.3); (v) the non-preemption invariant is
a **downgrade** of today's same-vector guarantee (§8.2). Plus a best-possible
caveat: this is the best *open-loop* design; a rotor index would supersede it
(§13).*

---

## 1. Goal

Make cross-board divergence on Phantasm **impossible to latch and fast to
self-correct at every layer** — column position, buffer-flip phase, *and*
displayed content (animation frame + effect) — without adding time-measurement
or branching to the per-column hot path, and **over a single wire** between
boards.

The insight driving the redesign: aligning the *display* (column + flip) is
necessary but **not sufficient** for a coherent image. Each board renders its
own copy of the effect in a free-running loop, so two boards can flip in perfect
time yet show different animation frames or even different effects. Coherence
requires a third, content layer that the current design does not synchronize.

The second insight (this revision): once each board runs a disciplined local
flywheel, the continuous per-column clock wire is **redundant**. A flywheel that
derives column position from time and snaps to a boundary mark twice per
revolution holds sub-column alignment on its own (§4.5). So the two-wire design
collapses to one — the low-rate **sync-symbol wire** — and the column clock wire
is deleted outright.

---

## 2. System model

Hardware constants (`targets/Phantasm/Phantasm.ino`, `core/platform.h`):

| Quantity            | Value    | Derivation                         |
|---------------------|----------|------------------------------------|
| RPM                 | 480      | 8 rev/s                            |
| `CANVAS_W` (W)      | 288      | columns per full image (360°)      |
| Column rate         | 2304 Hz  | `(RPM/60)·W`                       |
| Column period (T0)  | ~434 µs  | `1 / 2304`                         |
| Revolution          | 125 ms   | `60 / 480`                         |
| Half-revolution     | 62.5 ms  | one full image (two arms, 180°)    |
| Frame rate          | 16 fps   | 2 flips/rev                        |
| Boards (N)          | 4        | seg 0 = clock master / conductor   |
| Effect duration     | 120 s    | = 960 revolutions                  |

**Rotor model — open-loop (load-bearing).** There is no hall/index/encoder
input anywhere; both drivers synthesize the column clock at a *fixed* frequency
from nominal RPM (`pov_single.h:9`, `pov_segmented.h:259`). Therefore:

- "Sync" means *boards agreeing with master*, never *locking to the physical
  rotor*. Master's flywheel **is** the reference by definition.
- If the motor's true speed wanders, the whole image precesses uniformly. That
  is acceptable because every segment precesses together → still coherent.
- Absolute angular stability (image pinned to physical rotation) is **out of
  scope**; it would require a rotor index into master. See §13.

**Content clock — frame-counter, not wall-clock (load-bearing, good news).**
On-device animation advances off a per-frame counter `AnimationBase::t`,
incremented once per `step()` (`animation.h:160,226`), *not* `millis()`. So
content sync reduces to: (a) a shared `t` origin, (b) no dropped frames, (c)
deterministic RNG — and (c) is already done. The remaining gap is that the
effect **epoch** (playlist advance + `t` reset) is driven by per-board
`millis()` (`Phantasm.ino:47-53` → `run()`'s `while (millis()-start<…)`), so
boards change effects at different real moments.

**The single wire (role in the new design):**

- **Sync wire (the only inter-board connection).** Master emits width/count-coded
  symbols on `PIN_FRAME_SYNC_OUT` 3 → all boards decode on `PIN_FRAME_SYNC_IN` 4.
  Carries the **boundary** marks (Layer 2, 2/rev) and the **epoch** mark
  (Layer 3, 1/effect). These same boundary marks also discipline each board's
  Layer-1 flywheel (phase snap, optional frequency trim).
- **Deleted: the column clock wire.** The former `PIN_CLOCK_OUT` 5 →
  `PIN_COLUMN_SYNC` 2 link and master's fixed-frequency PWM are removed. Master
  no longer streams a per-column clock; it runs the same flywheel as every other
  board and only emits the low-rate sync symbols. Those two pins are freed.

---

## 3. Architecture: one flywheel, three disciplined layers, one wire

```
                     sync-symbol wire  (master → all, 2/rev + 1/effect)
                                    │  phase snap  (+ optional freq trim)
                                    ▼
   ┌──────────────────────────────────────────────────────────────┐
   │  LOCAL FLYWHEEL TIMEBASE  (per board, master included)         │
   │  column position derived from a FREE-RUNNING HARDWARE CLOCK    │
   │  (cycle counter), NOT from counting timer interrupts           │
   └───────────────┬───────────────┬───────────────┬───────────────┘
        derives     │      derives  │       advances │
                    ▼               ▼                ▼
   LAYER 1: column x_     LAYER 2: boundary     LAYER 3: content t
   (which column)         → advance_display     (which frame / effect)
                    ▲               ▲                ▲
                    │ snapped       │ flipped        │ epoch-reset
                    └── boundary ───┴── boundary ────┴── epoch ──── sync wire
                        symbol          symbol           symbol
```

Every layer reads the **same** disciplined timebase, so correcting the timebase
corrects all three coherently. Each layer also has an **absolute reference**
delivered on the one wire that pulls it back if it ever drifts:

| Layer | Local advance source | Absolute reference | Resync granularity |
|-------|----------------------|--------------------|--------------------|
| 1 Column/phase | flywheel (time-derived) | boundary-symbol snap (+ opt. freq trim) | every half-rev (sub-column) |
| 2 Frame/flip | flywheel boundary crossing | boundary symbol (width-coded) | every half-rev |
| 3 Content | `t++` per synced flip | epoch symbol + deterministic playlist | every effect (960 revs) |

Master's flywheel is the conductor: it free-runs and *defines* the reference;
every symbol it emits is timed from its own timebase. Downstream boards snap to
those symbols. There is no longer any board-to-board continuous clock — only the
master's flywheel as the common reference and the 2/rev + 1/effect symbols that
broadcast it.

---

## 4. Layer 1 — Column / phase: free-running flywheel, snapped by the sync wire

With the column clock wire gone, columns must be generated **locally on every
board**. The flywheel does it, and the sync wire's boundary symbols snap it back
into phase twice per revolution. The drift budget (§4.5) proves a Teensy crystal
holds sub-column alignment across the 62.5 ms between snaps, so a continuous
per-column reference buys nothing the snap doesn't already deliver.

This replaces today's design, where every column is an interrupt *slaved* to a
shared wire (one missed edge = one dropped column, permanent ±1 until the next
boundary snap). The flywheel inverts the dependency: **columns come from a local
timebase; the wire only snaps phase.**

### 4.1 The flywheel timebase — position from time, not from interrupt count (load-bearing)

This is the crux that makes everything else work. The flywheel must derive the
current column from a **free-running hardware clock** (ARM `DWT->CYCCNT` cycle
counter, or `micros()`), *not* by counting how many times its `IntervalTimer`
ISR has fired:

```
x_target = ((now_cycles - epoch_cycles) / cycles_per_column) mod W
```

Why this matters: `FastLED.show()` / the WS2801 bit-bang path masks interrupts
for windows **longer than T0** (`FastLED.show()` masks IRQs, §10). If `x_` were
incremented once per ISR fire, a masked window would coalesce several pending
timer interrupts into one fire and **lose** columns — exactly the dropped-column
bug we are eliminating. With position derived from time, the ISR that finally
runs after a mask simply reads the clock, computes the time-correct `x_target`,
and resumes there. The masked columns were not displayable anyway (the strip was
being clocked out), so jumping to the time-correct column is not just lossless —
it is the *correct* behavior. Counting ISR fires would instead resume one column
behind real angular position and smear.

The `IntervalTimer` at `T0` is therefore only a *wake-up*, not the source of
truth: it paces ISR entry near each column, but the cycle counter decides which
column it is. A late, early, or coalesced wake-up cannot inject drift.

### 4.2 Boundary-snap discipline (baseline)

Each boundary symbol (2/rev) names an absolute column (`x==0` or `x==W/2`,
§5.2). On decode, the board re-bases the flywheel so `x_` equals that boundary,
compensated for columns elapsed during the decode window (§5.2). Between snaps
the flywheel free-runs at *its own* nominal `T0` (its own crystal). The phase
error this accumulates over one half-revolution is the crystal's relative offset
times 62.5 ms — **0.006 column at a worst-case 40 ppm** (§4.5) — and the next
snap zeroes it.

No PI loop, no per-column phase math, no tuning. A free-running flywheel plus a
hard snap twice per revolution is the entire baseline, and it is already
sub-column.

### 4.3 Optional frequency trim (refinement, not required for correctness)

The snap corrects *phase*; it leaves a ~0.006-col sawtooth seam at each boundary
(the accumulated inter-snap drift, reset each half-rev). If that seam ever proves
visible — it will not until frequency error exceeds ~0.7 % (≈100× a Teensy
crystal, §4.5) — a board can additionally trim its column period to match
master's actual frequency, learned from the *interval between consecutive
boundary symbols* measured on its own cycle counter:

```
measured_half_rev = t_boundary_k - t_boundary_{k-1}
cycles_per_column ← (1-α)·cycles_per_column + α·(measured_half_rev / (W/2))
```

A slow first-order trim (small `α`) suffices because the quantity being tracked
— relative crystal offset — drifts only on thermal timescales (seconds to
minutes), not per-frame. This is a smoothness upgrade layered *on top of* the
snap, not a replacement for it and not a prerequisite. Ship snap-only first
(§11.1).

> This is the disciplined descendant of the full DPLL that earlier drafts
> proposed. With only one wire it collapses from a fast PI phase-loop (locking
> per-column edges) to a slow scalar frequency estimate (one update per
> half-rev), because phase is already handled by the snap. Simpler, and provably
> sufficient.

### 4.4 Layer-1 behavior summary

| Event | Baseline (snap-only) | With frequency trim |
|-------|----------------------|---------------------|
| Normal | flywheel free-runs at own T0; snapped 2/rev | period trimmed to master; snap nulls residual |
| Masked-IRQ window (`FastLED.show()`) | ISR resumes at time-correct column; no drift | same |
| 1 dropped boundary symbol | coasts ≤1 rev on own crystal (~0.01 col), re-snaps next | same, with even less coast error |
| 1 spurious symbol | width-decode + `try_flip` identity reject it | trim ignores out-of-window intervals |
| Sync wire dead | free-runs at nominal T0, precesses on own crystal | precesses on *trimmed* (last-known) frequency |

### 4.5 Clock drift budget

Crystal drift only accumulates while a clock free-runs *between* snaps. A board's
column counter ticks at `f = 2304·(1+δ)` Hz; two boards diverge at `2304·δ_rel`
columns/sec where `δ_rel` is the **relative** offset between their crystals. Time
to slip one column:

```
τ = T0 / δ_rel = 434 µs / δ_rel
```

| Relative offset δ_rel | τ (1 column) | drift rate |
|-----------------------|--------------|------------|
| 1 ppm  | 434 s (7.2 min) | 0.14 col/min |
| 10 ppm | 43 s            | 1.4 col/min  |
| 20 ppm | 22 s            | 2.8 col/min  |
| 30 ppm | 14 s            | 4.1 col/min  |
| 50 ppm | 8.7 s           | 5.5 col/min  |

Teensy 4's 24 MHz crystal is ~±10 ppm at room temperature, widening to ±30–50
ppm over its temperature range. Boards differ by the *relative* spread (up to the
sum), so realistically `δ_rel ≈ 20–40 ppm → one column in ~10–20 s`, a few
columns per minute *if uncorrected*. Shared enclosure/temperature pulls the
relative offset below the worst-case sum, but the design does not rely on that.

The whole point: the flywheel **is** corrected — snapped every half-rev — so it
never free-runs longer than 62.5 ms in normal operation:

| Regime | Snap interval | Drift before snap (@40 ppm) |
|--------|---------------|-----------------------------|
| Normal (snap-only) | every boundary (62.5 ms) | 2.5 µs ≈ **0.006 col** — sub-column |
| Normal (+ freq trim) | frequency matched between snaps | ~0 |
| 1 dropped symbol (coast 1 full rev) | 125 ms | 5 µs ≈ **0.01 col** — still sub-column |
| Sync wire dead (degradation) | never re-snapped | ≥1 col visible in ~10–20 s — slow smear, not a break |

Three consequences for the design decisions:

- **Snap-only is sub-column (0.006 col).** This is the quantitative
  justification for deleting the column clock wire: a continuous per-column
  reference would drive accumulated drift to ~0, but 0.006 col is already far
  below the visible threshold, so the wire earns nothing.
- **The seam is invisible until ~0.7 % frequency error.** One column of drift per
  half-rev needs `δ_rel = T0 / 62.5 ms ≈ 0.7 %` — about 100× a Teensy crystal's
  worst case. This is why §4.3's frequency trim is optional, not required.
- **A dead sync wire degrades gracefully over ~10–20 s**, not catastrophically.
  With the wire gone there is no re-snap, so unlike the old two-wire design the
  slip is no longer capped at ½ rev — it becomes a slow unbounded precession. The
  §9 "sync wire dead" row is "slow visible smear that needs the wire restored,"
  not "instant break."

This is purely **board-vs-board** (differential) crystal drift. It is distinct
from **image-vs-rotor** drift (motor holding RPM to ±%, not ±ppm — far larger),
which is *common-mode*: all boards follow master's flywheel, so the image
precesses uniformly and stays coherent (§2, rotor model). Crystal drift is the
only differential error between boards, and is exactly what the snap (and
optional trim) nulls.

---

## 5. Layer 2 — Frame / flip: exactly-once, self-describing boundaries (A + C)

The boundary (`x==0`, `x==W/2`) is derived from the disciplined `x_`. The flip
must fire **exactly twice per revolution on every board** regardless of column
drift, and a single sync glitch must not latch.

### 5.1 Exactly-once flip via boundary identity (A)

Boundaries strictly alternate (0, W/2, 0, W/2 …). A single deduplicated
primitive keyed on *which* boundary makes the flip exactly-once across both the
flywheel-crossing path and the sync-symbol path:

```
enum Boundary { NONE, ZERO, HALF };
static Boundary last_flipped_;            // shared, see §8

inline void try_flip(Boundary b) {
  if (b != last_flipped_) { effect_->advance_display(); last_flipped_ = b; }
}
```

`advance_display` is idempotent at `prev_==next_` (`canvas.h:387`), so a
redundant call is harmless; `try_flip` makes the common case exact.

Flip paths per board:

| Board | Flip trigger (normal) | Backstop |
|-------|-----------------------|----------|
| Master (0) | flywheel boundary crossing | — (its flywheel is the reference) |
| Downstream | flywheel boundary crossing (fires first) | boundary **symbol** |

A note on which path actually fires, because §5.2's "authoritative symbol" is
about the *snap*, not the flip: by timing the snapped flywheel crosses the
boundary **before** the symbol's falling-edge decode completes (~K cols later,
§11.6), so in normal operation the crossing triggers the flip and the symbol's
flip is the backstop. The symbol is authoritative for the **phase snap** (Layer
1, §4.2), not for triggering the flip. Only a drifted flywheel lets the symbol
flip first.

Downstream thus has two flip paths, deduped by `try_flip`: a flywheel that has
drifted slightly → the symbol still flips on time; a dropped *symbol* → the
flywheel crossing still flips (and is the normal trigger anyway); only losing
both in one half-rev glitches, self-healing next half-rev.

### 5.2 Self-describing boundary symbols (C)

Master codes the boundary into the sync symbol so each one names its boundary
*absolutely* — deleting the parity counter (`sync_seeded_`/`sync_at_zero_`,
`pov_segmented.h:444-450`) that today inverts permanently on one dropped pulse.
The symbols below are written as widths for concreteness, but the encoding is not
settled: the callout after the list recommends count/burst coding instead (§11.3).

- **NARROW** (≈2 columns high) = boundary HALF (`x==W/2`).
- **WIDE** (≈4 columns high) = boundary ZERO (`x==0`).
- **EXTRA-WIDE** (≈6–8 columns, or a coded burst) = boundary ZERO **+ EPOCH**
  (Layer 3). Epoch only ever lands on a ZERO boundary.

> **Encoding choice — count/burst over pure width (load-bearing, resolve in §11.3).**
> Width-coding is the obvious first cut, but it carries the symbol's meaning in
> the *spacing of two edges*, and the whole redesign exists because
> `FastLED.show()` masks IRQs for windows approaching T0 (§4.1, §10). A mask
> landing on a pulse's falling edge mis-measures its width → misclassified
> boundary → a bad snap *and* a bad flip on the **one** wire all three layers
> ride. The old edge-only frame-sync was immune to this (one edge, no width).
> Recommendation: encode each symbol as a **count of narrow pulses** (HALF = 1,
> ZERO = 2, ZERO+EPOCH = 3+, each pulse just an edge to count), so decode needs
> only the edge *count*, not accurate inter-edge timing — degrading gracefully
> when a mask window swallows part of a symbol. This generalizes the "coded
> burst" already proposed for EXTRA-WIDE to every symbol. Keep widths only if
> §11.3's margin analysis shows the worst-case mask window stays well clear of
> every threshold.

Generation (no blocking): master raises the line at the boundary in its
timebase ISR and lowers it via a column-counted countdown — width measured in
column ticks, no extra peripheral. Decode (cold ISR, `micros()` allowed): attach
the sync wire on **CHANGE**; on rising record `t_rise` + `x_at_rise`; on falling
classify width → Boundary (+epoch), snap `x_` to the boundary compensated for
columns elapsed during the measurement window, then `try_flip`. This same snap is
the Layer-1 phase discipline (§4.2) — one decode serves both layers.

> Decision pending (§11): elapsed-column-compensated snap vs. hard-snap with
> ≤K-column rewind. Recommendation: compensated (cheap, removes a visible seam
> under chronic drift).

---

## 6. Layer 3 — Content: epoch-synchronized animation + playlist

The gap finding #1's fix alone leaves open. Two divergence sources:

1. **Frame-index drift** — a board that drops a render falls behind in `t`
   permanently within an effect.
2. **Effect-transition drift** — the playlist advances on per-board `millis()`,
   so boards switch effects at different times and drift over the sequence.

### 6.1 Epoch symbol drives playlist + `t`

Master is the **conductor**. It counts revolutions on its own timebase and, when
an effect's 960 revolutions elapse, emits the **EPOCH** boundary symbol. On
epoch, every board (master included) atomically: advances to the next effect in
the deterministic roster, re-inits it, and resets `t=0`. Between epochs `t`
increments once per synced flip (Layer 2) — paced through the `buffer_free()`
gate, so flip count *is* `t` — and stays equal across boards **so long as each
keeps up rendering and sees no spurious flip**. A dropped render (§6.2) or a
spurious extra flip (§8.4) offsets that board's `t` by one until the next epoch;
on adjacent *spatial* segments (e.g. seg 0 / seg 1, the top and bottom of one
arm) a one-frame `t` offset is a visible junction seam on fast stateless
effects, so the epoch reset is what bounds it, not per-flip equality alone. This
replaces the `millis()`-gated `show<E>(120)` sequencing with epoch-counted
sequencing.

Because all boards iterate the **same** `HS_EFFECT_LIST` order, the epoch mark
need only say "advance," not "advance to N" — *provided* no epoch is ever missed
and all boards start at index 0 together.

### 6.2 Stateless vs. stateful (inherent regimes)

- **Stateless / parametric effects** (most of `animation.h`, driven directly by
  `t`): a dropped (or spurious-extra) frame leaves `t` permanently offset by one
  *within the effect* — the visible *hitch* clears next frame as motion resumes,
  but the absolute phase is one frame off until the epoch reset re-zeros `t`. On
  a single board that is imperceptible; the case that matters is two adjacent
  *spatial* segments drifting one frame apart (§6.1, junction seam). The
  in-principle distinction from stateful still holds — a stateless effect *would*
  resync instantly if a correct `t` were delivered, whereas a stateful one cannot
  be retro-synced even given `t`. But this design delivers no per-frame shared
  `t` mid-effect, so in practice both are bounded by the epoch reset; stateless
  just carries a far smaller (1-frame) artifact until then. If the junction seam
  ever proves visible, the cheapest fix is broadcasting low bits of `t` as a
  multi-bit content symbol (out of scope today).
- **Stateful integrators** (BZ/GS reaction-diffusion, trails, comets): cannot be
  retro-synced mid-effect (history is unrecoverable). The time-derived flywheel
  (Layer 1) minimizes dropped frames; the **epoch re-init bounds any stateful
  divergence to ≤ one effect**. This is inherent, not a flaw.

### 6.3 Epoch reliability

Epoch is rare (1 / 120 s) but a *missed* epoch is very visible (one segment
stuck on the previous effect for 120 s). So the epoch symbol must be more robust
than per-column timing:

> Decisions pending (§11):
> - **Redundancy:** repeat the EPOCH symbol on the next R ZERO-boundaries so a
>   single glitch can't drop it (effect changes are 960 revs apart — free).
> - **Absolute index (recommend baseline, not optional):** encode the effect
>   index (coded burst) so a late-booting or mid-show-rebooted board joins at the
>   *correct* effect rather than assuming index 0. Redundancy above only guards a
>   *dropped* symbol — it does nothing for a board that powers up at 0 while peers
>   are at K, which then never re-syncs (one segment on the wrong effect
>   indefinitely). The "all boot together at 0" assumption is the weakest link in
>   Layer 3. Under the project's fail-fast doctrine the alternative to encoding
>   the index is to **trap** when a board cannot establish it — *not* to silently
>   assume 0. Costs a real multi-bit symbol at 1/120 s; cheap insurance.

---

## 7. Sync-wire symbol alphabet (summary)

All inter-board information now rides this one wire:

| Symbol | Encoding (proposed) | Meaning | Rate |
|--------|---------------------|---------|------|
| NARROW | ~2 col high | boundary HALF (`x==W/2`) | 1/rev |
| WIDE | ~4 col high | boundary ZERO (`x==0`) | 1/rev |
| EXTRA-WIDE | ~6–8 col high (or coded burst, optionally repeated R×) | boundary ZERO + EPOCH (advance effect, reset `t`) | 1/effect |

Thresholds tunable (§11). All widths ≪ 62.5 ms half-rev, so consecutive symbols
never overlap. There is no separate column-clock wire; the boundary symbols both
flip the frame (Layer 2) and snap the column phase (Layer 1).

---

## 8. Concurrency & ISR model

ISRs on each board: **(a)** flywheel `IntervalTimer` (Layer-1 advance + pixel
pack — the hot path; reads the cycle counter for position, §4.1), **(b)** sync-wire
CHANGE (Layer-2/3 decode + snap + flip + epoch). **The former column-wire edge
ISR is gone** — there is no per-column input to service, which is both one fewer
ISR and the removal of the design's only per-column EMI surface (§10).

Invariants:

1. **Hot path stays branchless / time-light.** The per-column ISR does one cycle-
   counter read to compute `x_target` and packs pixels; no `digitalRead`, no
   symbol logic. Width measurement and epoch decode live only in the cold sync-wire
   ISR (≤ 2/rev + 1/effect).
2. **Non-preemption (load-bearing, *weakens* `pov_segmented.h:292-306`).** All
   shared state (`x_`, `last_flipped_`, the flywheel epoch/period registers,
   `t_rise`, `x_at_rise`, epoch flags) is touched only by these handlers, which
   must remain mutually non-preempting. Note this is a **downgrade** of today's
   guarantee, not an extension: today both handlers are GPIO pins on Teensy 4's
   single shared port vector, so they *literally cannot overlap regardless of
   priority*. The flywheel is now a PIT/`IntervalTimer` interrupt on a different
   vector, so that strong same-vector guarantee is gone — correctness now rests
   on the weaker "equal NVIC priority ⇒ no preemption" (a Cortex-M exception
   cannot preempt one of equal-or-lower priority). That holds only while both
   install at the *same* priority (`attachInterrupt` and `IntervalTimer` both
   default to 128) and SysTick stays demoted. This is now a cross-peripheral
   coupling of two default priorities — re-state it as a trap-worthy invariant at
   both attach sites; raising the flywheel timer above the sync-wire handler
   breaks it. Any new shared field carries the same note.
3. **Flip exactly-once** via `try_flip`'s identity check + idempotent
   `advance_display` backstop.
4. **Failure asymmetry honored — but only at Layer 2.** A *missed* flip (stale
   frame) is the only harmful **display** outcome; an *extra* flip is benign
   *for the display*, and the protocol biases toward flipping (two redundant
   paths) at every layer. Caveat for Layer 3: because content `t` is paced
   one-for-one by `advance_display` through the `buffer_free()` gate (§6.1, the
   render loop blocks until the ISR flips, then `step()`s once), a *spurious*
   flip that consumes an already-queued frame advances that board's `t` by one
   relative to its peers — a persistent content offset until the next epoch, not
   a benign no-op. `try_flip`'s identity dedup means the two *legitimate* paths
   never double-consume, and `advance_display` is idempotent when the buffer is
   already free, so this bites only on a genuine spurious boundary event landing
   while a fresh frame is queued (rare, epoch-bounded). Still: "extra flip is
   benign" is a Layer-2 statement, not a Layer-3 one.
5. **Seeding:** `run()` / epoch resets `x_=0`, `last_flipped_=NONE`,
   `t=0`, flywheel `epoch_cycles=now`, period = `T0`. First boundary is HALF
   (`HALF≠NONE` ⇒ flips).

---

## 9. Failure analysis (post-redesign)

| Event | Layer 1 (column) | Layer 2 (flip) | Layer 3 (content) |
|-------|------------------|----------------|-------------------|
| Masked-IRQ window (`FastLED.show()`) | flywheel resumes at time-correct column; no drift | unaffected | unaffected |
| 1 dropped boundary symbol | coasts ≤1 rev (~0.01 col); re-snaps next | crossing fallback flips | unaffected |
| 1 spurious symbol | width-decode + `try_flip` reject | identity check no-ops it | epoch redundancy guards it |
| 1 board renders slow (drops a frame) | — | shows prior frame 1 period | stateless: heals next frame; stateful: heals next epoch |
| 1 dropped epoch symbol | — | — | bounded by redundancy (R repeats); else 1 segment stale ≤120 s |
| Sync wire dead | boards free-run at T0, precess on own crystal (≥1 col in ~10–20 s) | crossing still flips 2/rev (no longer phase-snapped) | playlist drifts → degrades to per-board timing |
| Master dead | downstream flywheels free-run at T0, precess on own crystal (same as "sync wire dead" — master is just the symbol source) | crossing still flips 2/rev (no re-snap) | playlist drifts → per-board timing |

No single-glitch event latches a permanent error at any layer. The only visible
artifacts require either two coincident losses in one half-rev (self-heals) or a
dropped epoch (mitigated by redundancy). Losing the one wire is now a single
point of failure for *all three* layers at once — the cost of collapsing to one
wire — but it degrades as a slow smear (§4.5), not an instant break, and is a
wiring/connector concern addressed physically, not in protocol.

One consequence of the flywheel for the watchdog: the `buffer_free()` trap
(`canvas.h`) no longer fires on **master death**, because each downstream board
now generates its own columns and keeps flipping locally — master is only the
symbol source, not the column source. The watchdog's role narrows to catching a
board's *own* flywheel stalling (detached `IntervalTimer`, priority inversion),
which is the only thing that now stops `advance_display`. This is why the table
row above is graceful precession, not a trap — a correction from the current
design, where columns come from master's PWM and its death *does* starve the
column ISR into the watchdog.

---

## 10. Likelihood (why the flywheel matters)

A missed column today requires one board to coalesce two edges or latch a
spurious one; ranked causes: (1) interrupt-masked foreground windows > 434 µs —
dominant on the FastLED/WS2801 bit-bang path (`FastLED.show()` masks IRQs),
largely designed out on DMA; (2) worst-case ISR overrun; (3) EMI on a motorized
spinner (spurious edges). At ~1.1 M column edges/min across 4 boards, even a
1e-6 per-edge rate is a visible glitch every minute or two.

The time-derived flywheel **removes causes (1)/(2) as column-drop sources
outright** — a masked window or a long ISR just means the next ISR reads the
clock and resumes at the time-correct column (§4.1). And deleting the column
clock wire **removes cause (3) entirely from the column path**: there is no
per-column input left to pick up EMI. The only remaining EMI surface is spurious
*symbols* on the sync wire — 2/rev instead of 2304/rev — guarded by the symbol
encoding (count/burst recommended over width, §5.2) and epoch redundancy.
Collapsing to one disciplined flywheel is the core
robustness win of the redesign, and dropping the clock wire makes the hot path
strictly cleaner, not weaker.

---

## 11. Open decisions (resolve before implementing)

1. **Frequency trim:** ship snap-only first (§4.2 — simplest, already sub-column)
   and add the §4.3 frequency trim only if the ~0.006-col seam ever shows.
   Recommendation: **snap-only first.** If/when trimming, choose `α` and the
   out-of-window reject guard for the interval estimate.
2. **Free-running clock source:** `DWT->CYCCNT` (32-bit cycle counter, highest
   resolution) vs. `micros()` (simpler, 1 µs ≈ 0.0023 col quantization). Affects
   `x_target` precision and rollover handling.
3. **Symbol encoding — count/burst vs. width (see §5.2 callout):** recommend
   **count of narrow pulses** (HALF 1 / ZERO 2 / EPOCH 3+) so decode needs only
   an edge count, not accurate inter-edge timing, and survives a `FastLED.show()`
   mask window landing mid-symbol. Keep the proposed widths (NARROW 2 / WIDE 4 /
   EXTRA-WIDE 6–8 columns) only if the margin analysis shows the worst-case mask
   window *and* sync-ISR entry latency stay well clear of every width threshold.
4. **Snap compensation:** elapsed-column-compensated vs. hard-snap rewind.
   Recommend compensated.
5. **Epoch robustness:** repeat count R, **and** encode an absolute effect index
   (recommend baseline, §6.3) — redundancy guards a dropped symbol but not a
   late-boot/mid-show-reboot board joining at the wrong index; without the index,
   trap rather than assume 0.
6. **Inter-board flip skew:** in normal operation downstream flips on its **own
   flywheel crossing** (sub-column phase, §5.1), not the symbol — so skew is
   sub-column, not ~K cols. The ~K-col / ≤1.7 ms figure is the worst case where a
   drifted flywheel lets the symbol's falling edge flip first; still ≪ 62.5 ms
   frame and believed invisible. Confirm the crossing reliably precedes the
   symbol decode under normal drift.
7. **Share the flywheel with `pov_single`?** The single-board driver already has
   a local `IntervalTimer`; now that *master itself* runs the time-derived
   flywheel (no PWM clock-out), factoring a common flywheel core covers all three
   roles (single, master, downstream) and is more natural than before — but
   expands scope.

---

## 12. Test plan (host-testable where possible)

Following the `pov_segment_map.h` precedent (pure, host-tested index math):

- **Pure functions:** factor symbol→Boundary/epoch classification (whatever the
  chosen encoding, §11.3) and the `try_flip` state machine into free functions;
  assert exactly-once across interleaved crossing/symbol arrivals, and that a
  mask window swallowing part of a symbol degrades to a *missed* symbol (one
  crossing fallback), never a *misclassified* boundary.
- **Flywheel sim:** drive a mock time source (cycle-counter advances with
  injected masked-IRQ windows and a ppm frequency offset) plus a boundary-symbol
  stream with injected drops/dups/jitter; assert (a) no column lost across a
  masked window (position is time-correct on resume), (b) phase error bounded to
  the §4.5 budget between snaps, (c) phase re-acquired after a burst of glitches.
  For the optional frequency trim, assert the estimate converges and does not
  oscillate across the `α` range.
- **Layer-2 invariants:** `advance_display` count == 2/rev/board under arbitrary
  column drift; no latched inversion after a single symbol glitch.
- **Layer-3 content:** simulate 4 boards over a multi-effect playlist with
  injected per-board render-time jitter, a *spurious* extra flip (§8.4), a
  dropped epoch, and a mid-show board reboot; assert all boards on the same
  `(effect, t)` except within the proven bounds (stateless/spurious-flip: ≤1
  frame until epoch; stateful: ≤ 1 effect; dropped epoch: ≤ R-bounded; rebooted
  board: rejoins the correct effect index, or traps — §6.3, not "assumes 0").
- **Determinism harness** stays green (device-only ISR changes; host renders
  unaffected). The daydream `segment_controller`/`segment_worker` reproduces the
  4-board partition and is the end-to-end coherence check before fabrication.

---

## 13. Out of scope

- **Absolute angular lock to the rotor** — needs a hall/index sensor into
  master; the system is open-loop on rotation by design (§2). Uniform precession
  under motor-speed drift is accepted.

  *Best-possible caveat.* This is a **hardware** constraint (no sensor on the
  board), not a design-optimal choice — and it is the one axis on which this spec
  is not the best solution achievable, only the best *open-loop* one. A single
  index sensor into master would dominate the flywheel design on three fronts at
  once: (1) pin the image to physical rotation, removing the accepted precession;
  (2) give a crystal-drift-free revolution count, removing the §6.3 "never miss
  an epoch across 960 revs of crystal drift" fragility; and (3) provide a
  hardware boundary reference more robust than any wire symbol (§5.2). The whole
  flywheel-plus-snap apparatus exists to *approximate*, open-loop, what one
  sensor would give directly. If a sensor is ever fitted, this design is its
  graceful-degradation fallback, not a competitor.
- **Duplicate-master detection** (a peer holding the same hardware ID) — already
  out of scope (`pov_segmented.h:199-203`); needs an open-drain arbitration line.
- **A second (column-clock) wire / per-column genlock** — *removed*, not merely
  unused. The flywheel + 2/rev snap is proven sub-column (§4.5), so the
  continuous clock wire earns nothing and is deleted. Re-adding it would only buy
  the difference between 0.006 col and ~0 — below the visible threshold either
  way.
- **`pov_single` content/rotor changes** beyond an optional shared flywheel core.
