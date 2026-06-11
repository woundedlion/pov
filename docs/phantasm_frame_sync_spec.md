# Phantasm Synchronization Architecture — Design Spec

*Status: IMPLEMENTED, and this document describes the implementation. The
protocol core lives in `hardware/pov_sync.h` (pure, host-tested — flywheel
timebase, symbol alphabet, edge mailbox, acceptance gate, beacon codec,
content tracker, emitter); the device shell is `hardware/pov_segmented.h`
(two ISRs, pixel packing, effect handoff); `targets/Phantasm/Phantasm.ino`
supplies the roster factory table; the §12 test plan is
`tests/test_pov_sync.h` (pure units plus a 4-board simulator with crystal
offsets, single-latch mask windows, EMI, symbol drops, missed epochs, and
mid-show reboots). Refinements that emerged during implementation are folded
into their sections: oversampled wake-ups (§4.1), the suspect-burst
demarcation guard (§5.3), the join grid and per-effect RNG reseed (§6.1,
§6.5), the beacon schedule (§6.4), and the effect-handoff protocol (§8.5).
One validation item remains open on hardware: measuring the worst-case
mask window M on the shipped LED path (§11.3).*

*Original status notes: this is a full-stack redesign of Phantasm's
cross-board timing, of which code-review finding #1 (the
frame-sync/`advance_display` desync) was one symptom at the frame layer.
Decisions taken: (1) spec the hybrid fully before implementing; (2)
**collapse to a single inter-board wire** — delete the shared column clock
and generate columns from a local, time-derived flywheel on every board;
(3) flip ownership = reconcile + self-describing sync pulses (A+C).*

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
a **downgrade** of the old design's same-vector guarantee (§8.2). Plus a best-possible
caveat: this is the best *open-loop* design; a rotor index would supersede it
(§13).*

*Reliability review folded in (this revision): (a) the ISR model is now
**single-writer** — the sync-wire ISR is a pure edge publisher and the flywheel
ISR the sole consumer/owner of all sync state, deleting the fragile
cross-peripheral equal-priority invariant outright (§8.2); (b) symbol
acceptance gets a **plausibility gate + ACQUIRE/LOCKED acquisition states**
(§5.3), closing the residual misclassification and spurious-flip holes; (c)
master **self-censors late emissions** and emits pin-first (§5.2), so emission
jitter — previously unbudgeted, and potentially ~170× the §4.5 drift budget —
cannot poison downstream phase; (d) the flywheel gains a **rebase rule and
64-bit position math** (§4.1), making cycle-counter wrap structurally
impossible rather than handled; (e) the epoch becomes a **deadline-scheduled
commit** (§6.1) — "atomic re-init" hid a multi-ms foreground operation with
per-board skew; (f) a mid-revolution **index beacon** (§6.4) lets a rebooted
board rejoin at the correct effect within ~2 s instead of trapping or
staying wrong for 120 s. Open decisions §11.2/3/4/5/7 are resolved. The sync
wire itself is assumed physically reliable — a hard, soldered line by
construction — so wire-dead is **accepted as out of scope**, not designed for.*

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
requires a third, content layer that the previous design did not synchronize.

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
input anywhere; both drivers synthesize column timing at a *fixed* frequency
from nominal RPM (`pov_single.h`'s IntervalTimer period; the nominal
`cycles_per_half_rev` timebase constant in `pov_sync.h`). Therefore:

- "Sync" means *boards agreeing with master*, never *locking to the physical
  rotor*. Master's flywheel **is** the reference by definition.
- If the motor's true speed wanders, the whole image precesses uniformly. That
  is acceptable because every segment precesses together → still coherent.
- Absolute angular stability (image pinned to physical rotation) is **out of
  scope**; it would require a rotor index into master. See §13.

**Content clock — frame-counter, not wall-clock (load-bearing, good news).**
On-device animation advances off a per-frame counter `AnimationBase::t`,
incremented once per `step()` (`animation.h`), *not* `millis()`. So content
sync reduces to: (a) a shared `t` origin, (b) no dropped frames, (c)
deterministic RNG. For (c) the driver reseeds `hs::random()` (1337) at every
effect construction, so all boards render identical streams regardless of
boot/join history — the same per-effect determinism contract the simulator
uses. (a) is the gap the previous design left open: the effect **epoch**
(playlist advance + `t` reset) was driven by per-board `millis()`
(`show<E>(120)`'s `while (millis()-start<…)`), so boards changed effects at
different real moments; §6 replaces it with epoch-counted sequencing.

**The single wire (role in the new design):**

- **Sync wire (the only inter-board connection).** Master emits count-coded
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
| 2 Frame/flip | flywheel boundary crossing | boundary symbol (count-coded) | every half-rev |
| 3 Content | `t++` per synced flip | epoch symbol + deterministic playlist + index beacon (§6.4) | every effect; index re-verified every 16 revs |

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

This replaces the previous design, where every column was an interrupt
*slaved* to a shared wire (one missed edge = one dropped column, permanent ±1
until the next boundary snap). The flywheel inverts the dependency: **columns
come from a local timebase; the wire only snaps phase.**

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

The `IntervalTimer` is therefore only a *wake-up*, not the source of truth:
it paces ISR entry, but the cycle counter decides which column it is. A late,
early, or coalesced wake-up cannot inject drift.

**Wake-ups are oversampled 8× per column** (`kOversample` in
`pov_segmented.h`; the timer runs at T0/8 ≈ 54 µs). The timer grid is not
phase-locked to the flywheel — a snap shifts the flywheel's boundary instants
but not the timer — so at one wake per column, both the rendered column and
the master's first emitted pulse would lag their true instants by up to a
full, slowly beat-drifting column. Waking 8× per column bounds both
quantizations to ⅛ column, keeping the §5.2 self-censor budget meaningful
and the inter-board render skew far below a visible seam. The cost is ~18 kHz
of near-empty entries (one cycle-counter read, the 64-bit position
computation, a compare — ≈1 % CPU at 600 MHz); the pixel pack and DMA submit
still run only on a column change, at the same 2304 Hz as before.

**Wake-up contract.** Because wake-ups are advisory, the ISR must be
*idempotent* when `x_target` equals the column it last rendered (early or
beat-frequency double wake: do nothing) and *skip-tolerant* when `x_target` has
jumped (late or coalesced wake: render the time-correct column, never
back-fill the skipped ones — they were masked precisely because the strip was
busy).

**Clock source and rebase rule (load-bearing — resolves §11.2).** The clock is
`DWT->CYCCNT`: highest resolution (1.67 ns), zero cost to read. It is 32-bit at
600 MHz and wraps every **~7.16 s** — long enough to forget, short enough to
corrupt any elapsed-time computation that coasts. Rather than *handling* wrap,
make it structurally impossible: at every locally-crossed boundary the flywheel
folds its own epoch forward,

```
epoch_cycles += cycles_per_half_rev      // exact integer add, no drift
```

so `now_cycles - epoch_cycles` never exceeds ~63 ms regardless of how long the
board coasts between snaps. A snap re-bases `epoch_cycles` absolutely; the fold
keeps it fresh in between. (The same rule covers `micros()` and its 71.6-min
wrap if the cycle counter is ever unavailable.)

**Position math is 64-bit; the stored period is per-half-rev, not per-column.**
`cycles_per_column` ≈ 260,416.67 is **non-integer** — a truncated per-column
divisor bakes in a ~2.6 ppm systematic error and makes the §4.3 trim
inexpressible. Store `cycles_per_half_rev` (= 37,500,000 exactly at nominal
480 RPM) as the single, optionally-trimmed timebase constant and compute

```
x_target = (x_boundary + ((now_cycles - epoch_cycles) · (W/2))
                          / cycles_per_half_rev) mod W
```

with a 64-bit intermediate (the product peaks ≈ 5.4 × 10⁹). One long multiply
and divide per column ISR at 2304 Hz — negligible on the M7.

### 4.2 Boundary-snap discipline (baseline)

Each boundary symbol (2/rev) names an absolute column (`x==0` or `x==W/2`,
§5.2). On decode — and acceptance through the §5.3 gate — the board re-bases
the flywheel so `x_` equals that boundary,
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
snap, not a replacement for it and not a prerequisite. Shipped snap-only
(§11.1): the implementation exposes the trim hook
(`Flywheel::set_cycles_per_half_rev`, exercised at ±40 ppm by the tests) but
no estimator.

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
| 1 spurious symbol | count-decode + `try_flip` identity reject it | trim ignores out-of-window intervals |
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
- **A dead sync wire is out of scope (hard line by construction)** — but the
  behavior is still well-defined, not undefined: with no re-snap the boards
  precess apart at crystal rate (≥1 col in ~10–20 s, a slow smear, never an
  instant break), and the §4.1 rebase rule keeps the flywheel arithmetic valid
  indefinitely (no cycle-counter-wrap jump at 7.16 s). Nothing beyond that is
  designed for; the wire is a soldered, hard connection.

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

`advance_display` is idempotent at `prev_==next_` (`canvas.h`), so a
redundant call is harmless; `try_flip` makes the common case exact.

Flip paths per board:

| Board | Flip trigger (normal) | Backstop |
|-------|-----------------------|----------|
| Master (0) | flywheel boundary crossing | — (its flywheel is the reference) |
| Downstream | flywheel boundary crossing (fires first) | boundary **symbol** |

A note on which path actually fires, because §5.2's "authoritative symbol" is
about the *snap*, not the flip: by timing the snapped flywheel crosses the
boundary **before** the symbol's burst classification completes (burst + gap
timeout later, §5.2/§11.6), so in normal operation the crossing triggers the flip and the symbol's
flip is the backstop. The symbol is authoritative for the **phase snap** (Layer
1, §4.2), not for triggering the flip. Only a drifted flywheel lets the symbol
flip first.

Downstream thus has two flip paths, deduped by `try_flip`: a flywheel that has
drifted slightly → the symbol still flips on time; a dropped *symbol* → the
flywheel crossing still flips (and is the normal trigger anyway); only losing
both in one half-rev glitches, self-healing next half-rev.

### 5.2 Self-describing boundary symbols (C) — pulse-count encoding (FINAL)

Master codes the boundary into the sync symbol so each one names its boundary
*absolutely* — this deleted the previous design's parity counter
(`sync_seeded_`/`sync_at_zero_`), which inverted permanently on one dropped
pulse.

**Encoding (decision §11.3, RESOLVED): every symbol is a burst of short pulses
at fixed pitch; the meaning is the count of rising edges. Pulse widths carry no
information.**

- **1 pulse** = boundary HALF (`x==W/2`).
- **3 pulses** = boundary ZERO (`x==0`).
- **5 pulses** = boundary ZERO **+ EPOCH** (Layer 3; "advance" only — the
  absolute effect index, decision §11.5, does not ride the boundary symbol).
  Epoch only ever lands on a ZERO boundary.
- **Any other count (even, or >5) is INVALID: discard the whole symbol — no
  snap, no flip.** The flywheel crossing already covers the flip (§5.1) and the
  next boundary re-snaps 62.5 ms later, so a discarded symbol costs only a
  ~0.01-col coast (§4.5). *Fail to "missed," never to "wrong."*

Why count beats width, precisely: on i.MX RT each pin has **one** latched
interrupt flag, so an IRQ-mask window *delays* an edge's ISR but cannot lose
the edge — unless **two** edges land inside one mask window and share the
single latch. Therefore, with pulse pitch chosen **greater than the worst-case
mask window M**, the decoder's edge *count* is exact even when `FastLED.show()`
masks IRQs mid-symbol — whereas a width measurement is wrong whenever a mask
touches *any* edge, because the delayed ISR timestamps it late. Count decoding
moves the failure mode from "misclassified boundary" (a wrong-by-W/2 snap) to
"discarded symbol" (a free coast).

The odd-only, distance-2 alphabet {1, 3, 5} then guards the residual cases: a
single lost edge (two pulses inside one mask window, if M ever exceeds the
pitch margin) or a single spurious EMI edge each move the count by one — onto
an even, invalid value — and are discarded rather than misread.
Misclassification now requires **two coincident edge errors in one burst**; the
locked-mode snap plausibility gate (§11.7) rejects even those when the implied
correction or boundary identity contradicts the disciplined flywheel.

Parameters — all derived from the worst-case mask window M. **Phantasm ships
on the DMA LED path (`USE_DMA_LEDS`), where M ≈ 0** — confirm by measurement,
after which the margins below hold with order-of-magnitude headroom. (The
bit-bang fallback would require re-deriving them, and risks chronic master
self-censoring — see §9.1, row "lost boundary symbol.")

| Parameter | Value | Rule |
|-----------|-------|------|
| Pulse high time | one ISR body (pin HIGH at entry, LOW at exit; tens of µs) | width carries no information — only the rising edge registers; the glitch filter constrains edge *spacing*, not width |
| Pulse pitch | 2 columns (~868 µs) | **pitch > M** ⇒ no edge ever lost to the single latch |
| Burst gap timeout | 4 columns (~1.7 ms) | **timeout > pitch + M** ⇒ a mask-stretched gap cannot split one burst into two |
| Glitch filter | reject rising edges <100 µs apart | an EMI spike adds an isolated count → invalid → discard |

Generation (no blocking, no extra peripheral): the master's `SymbolEmitter`
schedules pulses in cycle time — the first due at the exact boundary instant,
the rest at 2-column pitch — and the flywheel ISR emits any due pulse at
wake-up, **pin write first, LED work after**, so emission timing carries only
the ⅛-column wake quantization (§4.1) plus ISR-entry jitter. If the master is
late at a boundary (the first pulse would start > ~½ column after the
boundary instant), it **self-censors**: it skips that boundary's symbol
entirely rather than emit a late one (downstream coasts one half-rev, ~0.01
col, §4.5); lateness detected *mid-burst* stops the remaining pulses. A
truncated count is usually invalid (discarded downstream); the two valid
truncations are harmless — 5→3 keeps the right boundary and only drops the
epoch flag (the §6.3 repeats cover it), and 3→1 names the wrong boundary,
which the §5.3 gate rejects on identity. Never emit a lie.

Decode (split across the two ISRs — the §8 single-writer model): the sync-wire
**RISING** ISR is a pure *publisher* — it applies the glitch filter, increments
the burst's edge count, and records the cycle-counter timestamp of the burst's
**first** edge into a small mailbox (`EdgeMailbox`); it touches no flywheel,
flip, or epoch state. The flywheel ISR — which wakes ~8× per column anyway —
is the sole *consumer*: when it observes the line quiet past the gap timeout
(compared on its own clock; no extra timer needed), it claims the burst under
a brief IRQ-off copy, classifies the count, applies the §5.3 acceptance gate,
and on acceptance snaps `x_` so the first-edge timestamp corresponds to the
named boundary — the elapsed-column compensation falls out for free, because
position is always derived as "time since epoch" — then `try_flip`s. This
same snap is the Layer-1 phase discipline (§4.2); one decode serves both
layers. (Classification *completes* up to burst + timeout after the boundary
— ≈13 columns for EPOCH, plus ≤⅛ column for the flywheel ISR to consume it.
That is irrelevant for the snap, which is timestamp-compensated, and sets the
§11.6 worst-case skew for the backstop flip.)

> Decision §11.4 — RESOLVED: **elapsed-column-compensated snap.** Cheap, removes
> a visible seam under chronic drift, and under count decoding it is *required*,
> since classification completes columns after the boundary instant.

### 5.3 Symbol acceptance — plausibility gate + acquisition states

The §5.2 alphabet makes single edge errors harmless; the gate handles
everything the alphabet cannot, by exploiting the fact that a disciplined
flywheel is provably sub-column accurate (§4.5) — so any symbol that *disagrees
substantially with the flywheel* is overwhelmingly more likely to be wrong than
the flywheel is. Two states per board (resolves §11.7's mechanism; the
constants are tunables):

- **ACQUIRE** (boot, mid-show reboot, or fallback): accept any *valid* symbol
  unconditionally — hard snap, no gate. The board displays **black** until it
  has both (a) phase, from one accepted boundary symbol, and (b) content
  identity, from an epoch or beacon (§6.4) — it never renders a guessed effect
  (fail-fast doctrine: show nothing rather than the wrong thing). First
  accepted snap → LOCKED.
- **LOCKED** (steady state): a valid symbol is accepted only if its implied
  phase correction is **≤ G columns** (G = 4) *and* its boundary identity
  matches the flywheel's nearest predicted boundary. Anything else is
  *rejected*: counted in telemetry (§8.6), no snap, no flip. After **R
  consecutive rejections** (R = 4, ≈2 revolutions) the board concludes its own
  timebase — not the wire — is at fault and falls back to ACQUIRE, hard-snapping
  to the next valid symbol. The fallback is mandatory: a gate without an escape
  deadlocks a genuinely-lost board into rejecting good symbols forever.

Two demarcation guards complete the mechanism (both fell out of host
simulation, not the original design):

- **ACQUIRE quiet-before guard:** a hard snap is taken only on a burst
  preceded by ≥ 16 columns of wire silence. Boundary symbols are isolated
  (half a revolution apart); beacon digits follow each other within ~12
  columns — so a digit train cannot capture a just-rebooted board mid-frame.
  (The train's *first* digit can still be mistaken once; the §9.1 mis-snap
  row bounds the recovery.)
- **Suspect bursts:** in LOCKED, a lone valid-count burst far (> G) from
  every predicted boundary cannot be told apart, at decode time, from a
  beacon's first digit — so it is fed to the beacon parser *and* held as a
  suspect until the beacon interdigit window (24 columns) passes. If another
  burst follows inside the window, it was beacon data; if the wire stays
  silent, it is counted as a gate rejection toward the R-fallback. Without
  this, a board with a corrupted timebase would route every REAL boundary
  symbol to the beacon parser (> G from its broken predictions) and never
  accumulate the R rejections — exactly the deadlock the fallback exists to
  prevent. Healthy boards pay nothing: a beacon's odd-count first digit is
  cleared by its own train, and the 2/rev accepted symbols reset the counter
  anyway.

What the gate buys on top of the alphabet: it rejects (a) the
two-coincident-edge-error residual — a misclassified boundary implies a ~W/2
correction, wildly implausible against a sub-column flywheel; (b) a late
emission that escaped master's self-censor (implied correction > G); (c)
spurious EMI bursts that happen to form a valid count (wrong position *and*
identity). Consequence for Layer 3: a spurious flip now requires a forged burst
that is simultaneously valid, plausible, and boundary-consistent — closing the
§8.4 spurious-flip `t`-offset hole in practice, and making "extra flip is
benign" true at Layer 3 as well as Layer 2.

---

## 6. Layer 3 — Content: epoch-synchronized animation + playlist

The gap finding #1's fix alone leaves open. Two divergence sources:

1. **Frame-index drift** — a board that drops a render falls behind in `t`
   permanently within an effect.
2. **Effect-transition drift** — the playlist advances on per-board `millis()`,
   so boards switch effects at different times and drift over the sequence.

### 6.1 Epoch symbol drives playlist + `t`

Master is the **conductor**. It counts revolutions on its own timebase and, when
an effect's 960 revolutions elapse, emits the **EPOCH** boundary symbol at
boundary B.

**The epoch commit is deadline-scheduled, not "atomic."** Effect re-init is a
multi-millisecond *foreground* operation — arena allocation, mesh builds — not
an ISR action, and per-board init times differ; an immediate switch would have
boards starting the new effect on different revolutions, a whole-sphere
mismatch every 120 s. Instead, EPOCH at boundary B means **commit at the
absolute boundary B+R+K** (R = redundancy repeats, K = construction window,
both fixed): the countdown runs in two phases. Through the **announce phase**
(B to B+R — the revolutions carrying the repeats) every board keeps playing
the outgoing effect. At B+R the **construction window** opens on every board
at once: each foreground tears down the old effect and constructs the next
one in the deterministic roster, displaying **black** for exactly K
revolutions (deterministic and identical on all boards — never a stale frame
on some and black on others); at B+R+K every board flips to the new effect's
frame 0 and resets `t=0` simultaneously. A board that accepted a *repeat*
rather than the primary copy counts down to the **same** boundary — see
§6.3.1 — which is why construction cannot start before B+R: only then is the
window's start common knowledge regardless of which copy each board heard,
and only then does every hearer get the full K-revolution build budget.
`HS_CHECK(init_complete)` at the deadline: an effect that cannot construct
inside K revolutions is an invariant violation and traps (fail-fast), rather
than silently skewing the show — and because the budget is K for every
hearer, the trap can only mean a firmware bug, never a missed symbol. K is
chosen comfortably above the slowest measured effect init. Construction
itself is deterministic across boards: the driver reseeds `hs::random()`
(1337) for every effect build, so the new instance is bit-identical no
matter what a board rendered — or whether it even existed — before the
epoch.

**Epoch dedup:** an accepted EPOCH opens a refractory window (~16 revs) in
which further EPOCH symbols are ignored — this is what makes the §6.3
redundancy repeats idempotent. Epochs are 960 revs apart, so the window is two
orders of magnitude clear of a real successor.

Between epochs `t`
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
need only say "advance," not "advance to N." The absolute index rides the §6.4
beacon instead, which both verifies each advance after the fact and lets a
late-booting board establish it without assuming "everyone started at 0
together."

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
  just carries a far smaller (1-frame) artifact until then. The §6.4 beacon
  carries the revolution count, so a slipped board *detects* the offset at the
  next beacon (`beacon_rev_mismatches` telemetry, §8.6) — but it does not
  retro-correct `t`: animation counters live inside each effect's animations
  and there is no frame fast-forward API, so correction stays epoch-bounded
  by design.
- **Stateful integrators** (BZ/GS reaction-diffusion, trails, comets): cannot be
  retro-synced mid-effect (history is unrecoverable). The time-derived flywheel
  (Layer 1) minimizes dropped frames; the **epoch re-init bounds any stateful
  divergence to ≤ one effect**. This is inherent, not a flaw.

### 6.3 Epoch reliability (resolved — three stacked mechanisms)

Epoch is rare (1 / 120 s) but a *missed* epoch is very visible (one segment
stuck on the previous effect). Three mechanisms stack, each catching what the
previous one cannot:

1. **Redundancy:** the EPOCH symbol is repeated on the next R ZERO-boundaries
   (R = 3); the §6.1 refractory window makes the repeats idempotent. A board
   must lose all R+1 in ~R revolutions to miss the advance itself. The
   repeats are **lockstep-safe**: the symbol carries no "which repeat am I"
   payload, but none is needed — the master starts the train exactly when
   `rev_in_effect` reaches 960, so a board hearing copy j infers j from its
   own (crossing-exact) revolution count and counts down to the same
   absolute B+R+K boundary as everyone else (`ContentTracker::
   on_epoch_symbol`). This also covers the master self-censoring its own
   primary emission (§5.2): downstream first hears the B+1 repeat and still
   commits with the master. The one case that cannot infer j — a board that
   beacon-joined mid-effect, whose revolution count is mod-64, lands outside
   the train window — falls back to j = 0 and commits up to j revolutions
   late for that one effect; its first commit re-zeros `rev_in_effect` in
   step with the master, so the next epoch is lockstep again.
2. **Absolute index on the beacon (§6.4):** a board that *does* miss every
   repeat — or that booted late, or rebooted mid-show — corrects at the next
   beacon, ≤16 revs (~2 s) later, instead of staying on the wrong effect for up
   to 120 s. No board ever *assumes* index 0; the "all boot together at 0"
   assumption is gone. (A beacon-corrected joiner starts the effect's history
   fresh mid-flight — stateless and stateful alike, since there is no frame
   fast-forward; full coherence restores at the next epoch.)
3. **Fail-dark, not fail-wrong:** until a board has established the index from
   an epoch or beacon it displays black (ACQUIRE, §5.3). If the index can never
   be established the segment stays dark — with the wire hard by construction
   (§9), that is the correct terminal fail-state, and it is visible at a glance
   rather than subtly wrong.

### 6.4 The index beacon — a data symbol off the boundary

Boundary symbols are timing-critical and must stay short; payload does not
belong on them. The beacon is a **data** symbol placed mid-revolution
(first burst starting at x ≈ W/4), where the wire is otherwise quiet and ±ms of
emission timing is irrelevant — the timing channel and the data channel are
separated *in time* on the same wire.

- **Frame:** five base-8 digits, each a burst of `digit + 1` pulses at
  1-column pitch, bursts separated by five quiet columns (one past the gap
  timeout, so the decoder reliably terminates each digit):
  `[index_hi, index_lo, rev_hi, rev_lo, checksum]` — effect index (0–63),
  revolution-count-within-effect mod 64, and a mod-8 sum checksum. Fixed digit
  count = end-of-frame detection; total worst-case length ≈ 26 ms at 1-column
  pitch, comfortably inside the half-rev.
- **Integrity model:** unlike boundary symbols (exactness via pitch > M, §5.2),
  the beacon tolerates corruption by *rejection*: any checksum mismatch, wrong
  digit count, or out-of-range digit drops the whole frame — the next beacon
  arrives ≤2 s later. This is why it may use the tighter 1-column pitch. A
  partial frame staler than the 24-column interdigit timeout is likewise
  dropped before the next burst starts a fresh frame.
- **Demarcation from boundary symbols:** a burst whose first edge lands far
  (> G columns, §5.3) from a predicted boundary is routed to the beacon
  parser; the LOCKED gate already excludes such bursts from snap/flip
  consideration, so the two parsers cannot claim each other's symbols. The
  converse hazard — a lost board's real boundary symbols all landing "far" —
  is what the §5.3 suspect-burst rule exists for.
- **Schedule:** revolution 1 of every 16 (rev ≡ 1 mod 16 — never revolution
  0, so a just-powered board meets clean, isolated boundary symbols before
  any data train), plus revolutions 1..R of a fresh effect (confirming the
  post-commit index immediately). Suppressed while a commit is pending: a
  beacon there would broadcast the outgoing index to joiners.
- **Consumers:** a LOCKED board cross-checks `(effect, rev)` — a mismatched
  index corrects a missed epoch (§6.3.2, triggering a rebuild and a §6.5
  grid-aligned rejoin); a mismatched rev is recorded as telemetry (§6.2). An
  ACQUIRE/joining board adopts `(effect, rev)` as its content identity and
  goes live at the next join-grid boundary (§6.5) with the effect at frame 0
  — there is no frame fast-forward, so only the grid-aligned boot realizes
  the ideal `t = 2·rev + parity` join; a mid-show rejoiner runs `t`-offset
  until the next epoch (§6.2).

### 6.5 The join grid — going live in lockstep

Boards take a constructed effect live only at ZERO boundaries where
`rev_in_effect ≡ 0 (mod 4)` (`join_grid_revs`; epoch commits are anchored by
B+R+K and unaffected). This is what makes **boot** coherent: without it, the
master would go live at its first boundary while downstream boards waited
~1–2 revolutions for beacon identity, leaving the master's `t` a few frames
ahead of its neighbors for the entire first effect — a visible junction
shear on fast effects. With the grid, every board (master included) holds
dark until the same crossing — typically revolution 4, ~0.5 s after power-on
— and starts frame 0 together. The grid length divides 64, so a
beacon-joined board's mod-64 revolution count lands on the same grid as the
master's true count; a mid-show rejoiner waits ≤ 4 revolutions (500 ms),
well inside the §9.1 ~2 s rejoin budget.

---

## 7. Sync-wire symbol alphabet (summary)

All inter-board information now rides this one wire:

| Symbol | Encoding (final, §5.2) | Meaning | Rate |
|--------|------------------------|---------|------|
| HALF | burst of 1 pulse | boundary HALF (`x==W/2`) | 1/rev |
| ZERO | burst of 3 pulses | boundary ZERO (`x==0`) | 1/rev |
| ZERO+EPOCH | burst of 5 pulses (repeated R×, §6.3) | boundary ZERO + EPOCH (advance at B+R+K, §6.1) | 1/effect |
| BEACON | 5 base-8 count digits at x≈W/4 (§6.4) | absolute effect index + rev count, checksummed | rev ≡ 1 (mod 16) + first R revs of an effect; silent during commits |
| *(invalid)* | any other count / failed checksum | discarded — no snap, no flip, no advance | — |

Boundary pulse pitch 2 columns, gap timeout 4 columns (§5.2); beacon digits at
1-column pitch under checksum (§6.4). All bursts ≪ 62.5 ms
half-rev, so consecutive symbols never overlap. There is no separate column-clock wire; the boundary symbols both
flip the frame (Layer 2) and snap the column phase (Layer 1).

---

## 8. Concurrency & ISR model

ISRs on each board: **(a)** flywheel `IntervalTimer` at T0/8 (§4.1) — the hot
path (one cycle-counter read for position, plus pixel pack on a column
change) *and* the sole owner of all sync state: `SyncBoard::tick()` consumes
the edge mailbox, runs gap-timeout classification, the §5.3 gate, the snap,
`try_flip`, and epoch scheduling; **(b)** sync-wire **RISING** (downstream
boards only — the master does not listen to its own emissions) — a pure
*publisher*: glitch filter, edge count, first-edge timestamp, written into a
small mailbox and nothing else. **The former column-wire edge ISR is gone** —
there is no per-column input to service, which is both one fewer ISR and the
removal of the design's only per-column EMI surface (§10).

Invariants:

1. **Hot path stays branchless / time-light.** Each wake-up does one
   cycle-counter read and the 64-bit position computation; ~7 of 8 entries
   end there (≈1 % CPU at 600 MHz — the foreground keeps the rest for
   rendering), and only a column change packs pixels and submits DMA
   (2304 Hz, same as the previous per-column interrupt). No `digitalRead`,
   no decode logic beyond a single mailbox check; the classify/snap/flip
   branch is taken ≤ 2/rev. All edge handling lives in the cold sync-wire
   ISR (≤ 2/rev of boundary edges + the occasional beacon). There is no
   busy-waiting anywhere in the protocol — emission alignment comes from the
   oversampled wake grid (§4.1), not from spinning to hit boundary instants.
2. **Single-writer ownership (replaces the equal-priority invariant).** Every
   piece of sync state — `x_`, `epoch_cycles`, `cycles_per_half_rev`,
   `last_flipped_`, lock state, epoch schedule, telemetry counters — has exactly
   **one writer: the flywheel ISR**. The sync-wire ISR writes only the mailbox
   (`edge_count`, `first_edge_cycles`, `last_edge_cycles`); the foreground only
   reads published flags. The mailbox handoff is a brief
   IRQ-off copy in the consumer — nanoseconds, not a masked window. With
   single-writer ownership the NVIC priority relationship between the two ISRs
   is **free**: the previous draft's cross-peripheral equal-priority invariant —
   itself a fragile downgrade of the old design's same-vector guarantee —
   is *deleted*, not restated. The sync ISR may
   preempt the flywheel ISR mid-column with no correctness consequence; worst
   case, a completed burst is consumed one column (434 µs) later, against a
   62.5 ms cadence.
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
   benign" is a Layer-2 statement, not a Layer-3 one — though the §5.3 gate
   closes the gap in practice: a spurious flip now requires a forged burst that
   is valid, plausible, *and* boundary-consistent at once.
5. **Lifecycle & effect handoff.** The flywheel and sync-wire ISRs attach
   **once per show and persist across epochs** — the flywheel *is* the
   timebase and is never detached between effects (the previous per-effect
   attach/detach did not carry over). The effect handoff keeps single-writer
   ownership: the foreground constructs and deletes instances; the ISR owns
   the live pointer. Protocol: the foreground polls `SyncBoard::build_word()`
   (a single aligned word, `generation << 8 | index`); on a generation change
   it bumps `release_req_`, the ISR drops its live pointer and acks within
   one wake-up, the foreground deletes the old instance, reseeds the RNG,
   constructs the new effect, draws its frame 0 (fresh buffers never block),
   and publishes it to the pending slot under a brief interrupts-off bracket.
   The ISR swaps pending → live only at a ZERO boundary: the §6.1 commit
   (`HS_CHECK` that the pending instance exists — the deadline trap) or a
   §6.5 join-grid crossing. During the dark window the ISR **keeps flipping
   the live effect** until released — the foreground may be blocked in the
   Canvas `buffer_free()` gate on its final frame of the outgoing effect, and
   `advance_display()` is what releases it to go tear the effect down. Boot
   seeding: ACQUIRE state (§5.3), display black, `last_flipped_=NONE` (the
   first accepted boundary flips, `HALF≠NONE`), `epoch_cycles=now`, period
   nominal; master is born LOCKED with identity (effect 0, rev 0) — it *is*
   the reference.
6. **Health telemetry (foreground-polled).** The flywheel ISR maintains
   counters — symbols accepted / gate-rejected / discarded (invalid count),
   beacons decoded / rejected, beacon index corrections and rev mismatches,
   epochs ignored by the refractory window, lock-state transitions, flips,
   emissions self-censored / aborted, and the longest coast (half-revs
   without a snap) — and the foreground render loop reports changes behind
   `hs::debug` at ≤1 Hz, exactly like the existing DMA-overrun counter
   pattern. Nothing in any ISR formats or prints. Degradation that the
   protocol absorbs silently (a discarded symbol, a rejected snap) must still
   be *visible* in one glance of debug output, or field diagnosis is
   guesswork.

---

## 9. Failure analysis (post-redesign)

| Event | Layer 1 (column) | Layer 2 (flip) | Layer 3 (content) |
|-------|------------------|----------------|-------------------|
| Masked-IRQ window (`FastLED.show()`) | flywheel resumes at time-correct column; no drift | unaffected | unaffected |
| 1 dropped boundary symbol | coasts ≤1 rev (~0.01 col); re-snaps next | crossing fallback flips | unaffected |
| 1 spurious symbol | count alphabet discards (even count) or §5.3 gate rejects | identity check no-ops it | epoch refractory + gate guard it |
| Late-emitted symbol (master masked) | master self-censors (§5.2); residual rejected by gate (§5.3) | crossing flips on time regardless | unaffected |
| 1 board renders slow (drops a frame) | — | shows prior frame 1 period | stateless: heals next frame (or next beacon, §6.4); stateful: heals next epoch |
| 1 dropped epoch symbol | — | — | R repeats; missed-all-R corrected by next beacon ≤16 revs (~2 s) |
| Board reboots mid-show | ACQUIRE: hard-snaps to first valid symbol | flips resume on first accepted boundary | black until index from beacon (≤2 s), then rejoins at the correct effect (frame 0, §6.5 grid; `t` offset until the next epoch) |
| Sync wire dead *(out of scope — hard line)* | free-runs at T0, precesses on own crystal (≥1 col in ~10–20 s); rebase rule keeps arithmetic valid (§4.1) | crossing still flips 2/rev | playlist freezes on current effect (epoch never arrives); ACQUIRE boards stay dark |
| Master dead | downstream flywheels free-run at T0, precess on own crystal (same as "sync wire dead" — master is just the symbol source) | crossing still flips 2/rev (no re-snap) | playlist freezes on current effect |

No single-glitch event latches a permanent error at any layer, and with the
§5.3 gate even *multi*-error symbol corruption fails to "discarded," not
"wrong." The only visible artifacts require either two coincident losses in one
half-rev (self-heals) or a missed epoch (R repeats, then beacon-corrected
within ~2 s). Losing the one wire is a single point of failure for all three
layers at once — the accepted cost of collapsing to one wire. It is **out of
scope by construction** (a hard, soldered line, not a connector); the rows
above exist to show the degradation is a well-defined slow smear rather than
undefined behavior, not because it is designed against.

### 9.1 Failure-mode budget (DMA path)

Worst-case recovery and expected frequency per failure mode, on the shipped
DMA LED path (mask window M ≈ 0, so all masked-IRQ modes are non-events).
Constants: gate G = 4 col, fallback R = 4 rejections, construction window
K = 2 revs (commit at B+R+K, §6.1), beacon every 16 revs, EPOCH ×(1+3). EMI rate anchor: λ ≈ 1 induced event/min —
the §10 old-design glitch estimate, deliberately pessimistic for a terminated
hard line. Time anchors: 1 col = 434 µs; 144 col = ½ rev = 62.5 ms;
4,608 col = 16 revs = 2 s; 276,480 col = 960 revs = 120 s.

| Failure mode | Worst-case artifact | Worst-case recovery | Expected frequency |
|---|---|---|---|
| Crystal drift between snaps (normal operation) | 0.006 col phase error | 144 col (next snap) | continuous; ~40× below visibility |
| Lost boundary symbol (discarded burst, glitch-filtered EMI, master self-censor) | coast error 0.006 → 0.01 col | 288 col (1-rev coast) | ≈0 on DMA; harmless at any plausible rate |
| EMI on the sync wire | binding case: an edge within G of the matching predicted boundary → ≤G col (≈5°) seam on one board for ≤½ rev; all other cases rejected with no artifact | ≤144 col (next real symbol re-snaps) | accepted-case ≈ λ·2G/288 ≈ **1.7/hr**; rejected ≈ λ ≈ 1/min (telemetry only); misclassification (2 coincident errors) ≈ 1/2 yrs *and* gate-rejected |
| Mis-snap despite the gate / corrupted timebase (incl. forged burst during ACQUIRE) | one board off by up to W/2 | ≤720 col ≈ 313 ms (R rejections = 576 col → ACQUIRE → re-snap ≤144) | effectively never — needs a 2-coincident-error burst *during* a ~2 s ACQUIRE window, or a firmware bug; the fallback bounds it either way |
| Dropped render (effect misses the 62.5 ms budget) | stale frame for 1 period; 1-frame `t` seam vs neighbors | display 144 col; `t`: ≤4,608 col via beacon (stateless) / ≤276,480 col via epoch (stateful) | ≈0 within budget; watched by the overrun/`ft` telemetry |
| Missed epoch (all R+1 copies) or corrupted beacon frame | one segment on the old effect ≤2 s; a dropped beacon alone is consequence-free redundancy | ≤4,608 col (~2 s, next beacon) | ≈0 — requires 4 independent symbol losses; beacon bounds it regardless |
| Board reboot mid-show | one segment dark (fail-dark, never wrong) | ≤4,608 col (~2 s): phase ≤144 col, index at next beacon, rejoins at the correct effect on the §6.5 grid | per external reboot event |
| Firmware invariant violation (init > K, flywheel stall) | trap (`HS_CHECK` / `buffer_free()` watchdog) | none — fail-fast by design | 0 in correct firmware; a caught bug class, not a runtime mode |
| Sync wire / master dead | uniform slow smear ~1 col per 10–20 s; playlist freezes; arithmetic stays valid (§4.1 rebase) | physical repair | out of scope — hard line by construction |

Reading by tier: everything the wire can plausibly throw at the design recovers
sub-column within ≤2 revolutions; the only in-principle-visible stochastic
artifact is the accepted-EMI case (~5° one-board seam for one half-rev,
~1.7/hr at the pessimistic λ — halve G to 2 to halve both rate and magnitude
if it ever shows); content-layer slips are beacon-bounded to 2 s; firmware
defects trap rather than recover.

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

In the old slaved-clock design, a missed column required one board to
coalesce two edges or latch a spurious one; ranked causes: (1)
interrupt-masked foreground windows > 434 µs —
dominant on the FastLED/WS2801 bit-bang path (`FastLED.show()` masks IRQs),
largely designed out on DMA; (2) worst-case ISR overrun; (3) EMI on a motorized
spinner (spurious edges). At ~1.1 M column edges/min across 4 boards, even a
1e-6 per-edge rate is a visible glitch every minute or two.

The time-derived flywheel **removes causes (1)/(2) as column-drop sources
outright** — a masked window or a long ISR just means the next ISR reads the
clock and resumes at the time-correct column (§4.1). And deleting the column
clock wire **removes cause (3) entirely from the column path**: there is no
per-column input left to pick up EMI. The only remaining EMI surface is spurious
*symbols* on the sync wire — 2/rev instead of 2304/rev — guarded by the
count-coded symbol alphabet (§5.2) and epoch redundancy.
Collapsing to one disciplined flywheel is the core
robustness win of the redesign, and dropping the clock wire makes the hot path
strictly cleaner, not weaker.

---

## 11. Decision log (all resolved; shipped values noted)

1. **Frequency trim — SHIPPED snap-only** (§4.2 — simplest, already
   sub-column). The trim hook exists (`Flywheel::set_cycles_per_half_rev`,
   exercised at ±40 ppm by the tests) with no estimator; add §4.3 only if the
   ~0.006-col seam ever shows, choosing `α` and the out-of-window reject
   guard then.
2. **Free-running clock source — SHIPPED: `DWT->CYCCNT` + the §4.1 rebase
   rule.** Highest resolution, free to read; the 7.16 s wrap is made
   structurally impossible by folding `epoch_cycles` forward every boundary
   (the same rule covers `micros()` if CYCCNT is ever unavailable). Position
   math is 64-bit with `cycles_per_half_rev` (= 37,500,000 nominal) as the
   single stored, optionally-trimmed timebase constant (§4.1).
3. **Symbol encoding — SHIPPED: pulse-count bursts (§5.2).** Alphabet HALF=1 /
   ZERO=3 / ZERO+EPOCH=5 rising edges at 2-column pitch, gap-timeout-terminated;
   even/other counts are invalid and discarded whole. Chosen over width because
   the pin's single latched IRQ flag makes edge *counts* exact whenever pitch >
   worst-case mask window M (a masked ISR is delayed, not lost), while widths
   are corrupted by any mask touching any edge; the odd-only distance-2 alphabet
   turns single edge errors — lost *or* spurious — into discards instead of
   misclassifications, and master self-censors late/interrupted emissions.
   **Still open on hardware: measure M on the shipped LED path** and confirm
   the pitch/timeout margins in the §5.2 table (expected ≈ 0 on the DMA
   path; the host simulator already covers M up to several columns).
4. **Snap compensation — SHIPPED: elapsed-column-compensated** (§5.2).
   Inherent to position-from-time (a snap re-bases the epoch to the
   first-edge timestamp and position is always "time since epoch"), required
   under count decoding, and it removes the visible seam under chronic drift.
5. **Epoch robustness — SHIPPED: three stacked mechanisms (§6.3/§6.4).**
   R repeats made idempotent by the §6.1 refractory window; absolute effect
   index + revolution count on the mid-rev beacon (≤2 s correction for any
   missed epoch or late-booting board); fail-dark in ACQUIRE rather than
   assume index 0. Repeats are lockstep-safe: every heard copy counts down
   to the same absolute B+R+K boundary (§6.3.1). Shipped constants: R = 3,
   beacon period 16 revs, refractory window 16 revs, construction window
   K = 2 revs (HS_CHECK-trapped; confirm against the slowest measured effect
   init on hardware), join grid 4 revs (§6.5).
6. **Inter-board flip skew — RESOLVED by construction.** In normal operation
   downstream flips on its **own flywheel crossing** (within ⅛ column of the
   boundary under the oversampled wake grid, §4.1), and symbol
   classification cannot complete before burst + gap timeout (≥ ~4.5
   columns) — so the crossing always leads and skew is sub-column. The
   backstop case (a drifted flywheel letting the symbol flip first) lands at
   ≈13 columns / ~5.6 ms worst case for EPOCH; still ≪ 62.5 ms frame. The
   simulator asserts exactly 2 flips/rev/board across drift, masks, EMI, and
   drops.
7. **Snap plausibility gate + acquisition states — SHIPPED (§5.3).** ACQUIRE
   (hard snap behind the quiet-before guard, fail-dark) ⇄ LOCKED (correction
   ≤ G columns; with G < W/4 the distance gate subsumes the
   boundary-identity check; R consecutive rejections — including
   suspect-burst timeouts — fall back to ACQUIRE so the gate can never
   deadlock a lost board). Shipped constants: G = 4 columns, R = 4.
8. **Share the flywheel with `pov_single`? — NOT DONE (future option).** The
   single-board driver keeps its per-column IntervalTimer ISR. Now that
   master itself runs the time-derived flywheel, factoring a common flywheel
   core covering all three roles (single, master, downstream) is natural —
   but it expands scope and Holosphere has no sync problem to solve.

---

## 12. Test plan (host-testable where possible)

Implemented as `tests/test_pov_sync.h`: pure units for every protocol piece
plus a 4-board event-driven simulator — per-board crystal ppm offsets, a
single-latch masked-IRQ model (edges during a mask merge and arrive late,
the i.MX RT pin-flag behavior the count coding is designed around), EMI
injection, symbol-drop windows, foreground build delays, mid-show reboots,
and the commit-deadline trap. The simulator earned its keep before any
hardware existed: it caught a master emitter that went permanently idle
after its first beacon frame, and the §5.3 demarcation deadlock that the
suspect-burst rule now closes. Two deltas against the plan below: the
*spurious-flip* case is realized as forged-burst rejection (the gate makes a
deliverable spurious flip unconstructible short of a valid + plausible +
boundary-consistent forgery), and rejoiners are asserted on `(effect, rev)`
— `t` equality is asserted for the grid-aligned boot and post-epoch, since
there is no frame fast-forward (§6.4).

Following the `pov_segment_map.h` precedent (pure, host-tested index math):

- **Pure functions:** factor symbol→Boundary/epoch classification (pulse-count
  encoding, §5.2) and the `try_flip` state machine into free functions;
  assert exactly-once across interleaved crossing/symbol arrivals, and that a
  mask window swallowing part of a symbol degrades to a *missed* symbol (one
  crossing fallback), never a *misclassified* boundary.
- **Flywheel sim:** drive a mock time source (cycle-counter advances with
  injected masked-IRQ windows and a ppm frequency offset) plus a boundary-symbol
  stream with injected drops/dups/jitter; assert (a) no column lost across a
  masked window (position is time-correct on resume), (b) phase error bounded to
  the §4.5 budget between snaps, (c) phase re-acquired after a burst of glitches.
  (The shipped tests exercise the trim *hook* at ±40 ppm; estimator
  convergence tests come with §4.3 if it ever ships.)
- **Layer-2 invariants:** `advance_display` count == 2/rev/board under arbitrary
  column drift; no latched inversion after a single symbol glitch.
- **Timebase arithmetic:** run the flywheel sim past a 32-bit cycle-counter
  wrap (>7.16 s of mock time) with sparse snaps and assert the §4.1 rebase rule
  keeps `x_target` continuous across the wrap; exercise the 64-bit position
  math at the trim extremes (`cycles_per_half_rev` ± worst-case ppm) and assert
  zero accumulated truncation drift over thousands of revolutions.
- **Emission self-censor:** inject master-side ISR lateness (mask windows over
  scheduled emission instants); assert the symbol is suppressed (not emitted
  late), a mid-burst interruption truncates to an invalid count, and no
  downstream snap correction ever exceeds the §5.3 gate as a result.
- **Acceptance gate:** assert a forged burst implying a ~W/2 correction (the
  two-coincident-edge-error residual) is rejected in LOCKED; assert a board
  seeded with a corrupted timebase re-acquires within R symbols via the
  ACQUIRE fallback (no rejection deadlock); assert mid-rev beacon bursts are
  never consumed as boundary symbols and vice versa (§6.4 demarcation).
- **Epoch commit:** simulate per-board init-time spread inside the K-rev
  window; assert every board flips to the new effect's frame 0 at exactly
  B+R+K with black during the construction window (and NOT during the
  announce phase), and that an init exceeding K traps (`HS_CHECK`), never
  silently skews. Assert repeat-lockstep (§6.3.1): a board deafened for just
  the primary copy, and the whole downstream side when the master
  self-censors its primary emission, commit at the same boundary as their
  peers with equal frame counters after.
- **Beacon:** corrupt single digits, the checksum, and the digit count; assert
  every corrupted frame is dropped whole (no partial application) and that a
  rebooted board joins at the correct `(effect, rev)` from the next good
  beacon (display from the next §6.5 grid boundary).
- **Layer-3 content:** simulate 4 boards over a multi-effect playlist with
  injected per-board build delays, a forged plausible burst (§8.4's
  spurious-flip vector, gate-rejected), a dropped epoch train, and a mid-show
  board reboot; assert all boards on the same `(effect, t)` except within the
  proven bounds (dropped epoch: beacon-corrected within one beacon period +
  join grid; rebooted board: dark through ACQUIRE, then rejoins at the
  correct effect — §6.3/§6.4, never "assumes 0" — with `t` offset until the
  next epoch).
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
  out of scope (`pov_segmented.h` `read_id()`); needs an open-drain arbitration
  line.
- **A second (column-clock) wire / per-column genlock** — *removed*, not merely
  unused. The flywheel + 2/rev snap is proven sub-column (§4.5), so the
  continuous clock wire earns nothing and is deleted. Re-adding it would only buy
  the difference between 0.006 col and ~0 — below the visible threshold either
  way.
- **`pov_single` content/rotor changes** beyond an optional shared flywheel core.
