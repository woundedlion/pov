# Phantasm Synchronization Architecture — Design Spec

*Status: DRAFT for review. Not yet implemented. This is a full-stack redesign of
Phantasm's cross-board timing, of which code-review finding #1
(`pov_segmented.h:420-451`, the frame-sync/`advance_display` desync) is one
symptom at the frame layer. Decisions taken so far: (1) spec the hybrid fully
before implementing; (2) keep the shared column-clock wire as the genlock
reference; (3) flip ownership = reconcile + self-describing sync pulses (A+C).*

---

## 1. Goal

Make cross-board divergence on Phantasm **impossible to latch and fast to
self-correct at every layer** — column position, buffer-flip phase, *and*
displayed content (animation frame + effect) — without adding time-measurement
or branching to the per-column hot path.

The insight driving the redesign: aligning the *display* (column + flip) is
necessary but **not sufficient** for a coherent image. Each board renders its
own copy of the effect in a free-running loop, so two boards can flip in perfect
time yet show different animation frames or even different effects. Coherence
requires a third, content layer that the current design does not synchronize.

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
  rotor*. Master's clock **is** the reference by definition.
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

**Wires (roles in the new design):**

- **Wire 1 — column clock (KEPT, role changes).** Master's fixed-freq PWM
  (`PIN_CLOCK_OUT` 5 → all boards' `PIN_COLUMN_SYNC` 2). No longer the *primary*
  clock that each edge directly advances; now the **genlock phase reference**
  for each board's local flywheel timer (Layer 1).
- **Wire 2 — frame sync (ENRICHED).** Master emits width/count-coded symbols on
  `PIN_FRAME_SYNC_OUT` 3 → all boards decode on `PIN_FRAME_SYNC_IN` 4. Carries
  the **boundary** marks (Layer 2) and the **epoch** mark (Layer 3).

---

## 3. Architecture: one flywheel, three disciplined layers

```
                       shared column wire (genlock ref, ~2304 Hz)
                                    │  phase/freq correction
                                    ▼
   ┌──────────────────────────────────────────────────────────────┐
   │  LOCAL FLYWHEEL TIMEBASE  (per board: IntervalTimer @ ~T0)     │
   │  free-runs; coasts through missed reference edges              │
   └───────────────┬───────────────┬───────────────┬───────────────┘
        advances    │      derives  │       advances │
                    ▼               ▼                ▼
   LAYER 1: column x_     LAYER 2: boundary     LAYER 3: content t
   (which column)         → advance_display     (which frame / effect)
                    ▲               ▲                ▲
                    │ disciplined   │ corrected      │ epoch-reset
                    └── column wire │── boundary ────│── epoch ──── sync wire
                                    │   symbol        │   symbol
```

Every layer reads the **same** disciplined timebase, so correcting the timebase
corrects all three coherently. Each layer also has an **absolute reference**
delivered on a wire that pulls it back if it ever drifts:

| Layer | Local advance source | Absolute reference | Resync granularity |
|-------|----------------------|--------------------|--------------------|
| 1 Column/phase | flywheel timer tick | column-wire edge (phase+freq) | every column (sub-column lock) |
| 2 Frame/flip | flywheel boundary crossing | boundary symbol (width-coded) | every half-rev |
| 3 Content | `t++` per synced flip | epoch symbol + deterministic playlist | every effect (960 revs) |

---

## 4. Layer 1 — Column / phase: flywheel genlocked to the column wire

Today every column is an interrupt *slaved* to the shared wire, so one missed
edge = one dropped column (permanent ±1 until the next boundary snap). The
flywheel inverts the dependency: **columns come from a local `IntervalTimer`;
the wire only disciplines it.** A missed/spurious reference edge → the flywheel
coasts → no dropped column.

Two formulations; the spec recommends the simpler **X** because the wire is kept
and is reliable-when-present, then notes **Y** for smoothness / future wire-drop.

### 4.1 Formulation X — wire drives, timer backstops (recommended baseline)

- Column-wire edge ISR advances `x_` and packs pixels (today's hot path,
  essentially unchanged).
- A local `IntervalTimer` runs at `T0` and is **re-armed by every edge**, so in
  steady state it never fires (the edge always beats it). If an expected edge is
  late by > a guard (e.g. `1.5·T0`), the timer fires and advances `x_` itself —
  filling exactly the missed column — then continues until the next real edge
  re-arms it.
- Net: the wire provides timing and frequency when present; the timer fills
  gaps. No per-column phase math, no PLL. Minimal change to the hot path.

Cost: one timer re-arm per edge (cheap on Teensy). Caveat: a *spurious* extra
edge still over-advances `x_` by one — corrected at the next boundary symbol
(Layer 2), exactly as today, but now bounded and content-safe via Layers 2–3.

### 4.2 Formulation Y — flywheel drives, wire disciplines (option, full DPLL)

- The `IntervalTimer` ISR is the *only* thing that advances `x_` (uniform column
  pacing, immune to edge-arrival jitter).
- The column-wire edge ISR runs a **phase detector**: `e_k = t_edge −
  t_predicted`. A PI loop adjusts the timer period:
  `T_{k+1} = T0 − Kp·e_k − Ki·Σe`, slew-limited for stability (start
  `Kp≈1/16`, `Ki` smaller). This locks both phase and *frequency* (this board's
  crystal vs master's).
- With frequency locked, residual drift between *any* two references is
  ≈ `ppm · interval` (≈10–20 ppm Teensy crystals → sub-column even at one
  reference per half-rev). This is what would make dropping Wire 1 safe later;
  with Wire 1 kept it is a smoothness upgrade, not a necessity.

> Decision pending (§11): ship X first (simplest, fixes dropped columns) and
> treat Y as a follow-up, or build Y directly. Recommendation: **X**.

### 4.3 Layer-1 behavior summary

| Event | Formulation X | Formulation Y |
|-------|---------------|---------------|
| Normal | edge advances; timer idle | timer advances; edge trims phase |
| 1 missed edge | timer fills the column | timer coasts; next edge re-trims |
| 1 spurious edge | +1 column; fixed at next boundary | phase detector rejects/averages it out |
| Wire 1 dead | board free-runs at nominal T0 (precesses on own crystal) | same, but frequency-trimmed to last lock |

### 4.4 Clock drift budget

Crystal drift only accumulates while a clock free-runs. A board's column counter
ticks at `f = 2304·(1+δ)` Hz; two boards diverge at `2304·δ_rel` columns/sec
where `δ_rel` is the **relative** offset between their crystals. Time to slip one
column:

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
columns per minute. Shared enclosure/temperature pulls the relative offset below
the worst-case sum, but the design does not rely on that.

What this means per regime — note how little it bites until Wire 1 is gone:

| Regime | Correction interval | Drift before correction (@40 ppm) |
|--------|---------------------|-----------------------------------|
| Wire 1 alive (normal) | every column (434 µs) | ~0 (re-locked each edge) |
| Missed-edge gap (flywheel coast) | a few columns (~1 ms) | ~40 ns ≈ 0.0001 col — invisible |
| Wire-2-only (half-rev = 62.5 ms) | per boundary symbol | 2.5 µs ≈ **0.006 col** — sub-column |
| Wire 1 dead (degradation) | only Wire-2 boundary snaps | ≥1 col visible in ~10–20 s, but Wire-2 still re-snaps `x_` 2/rev → slip bounded to ≤ ½ rev, not unbounded |

Three consequences for the design decisions:

- **Keeping Wire 1 (§11.1) gives per-column genlock → zero accumulation.** This
  is the quantitative justification for the "keep the wire" decision.
- **Sub-column drift between half-rev syncs (0.006 col)** is the proof that
  Formulation Y's frequency lock + sparse sync *could* hold alignment without
  Wire 1 — i.e. why Y is a viable future option, not a present necessity.
- **A dead column wire degrades gracefully over ~10–20 s**, not catastrophically,
  and Wire-2 boundary snaps cap the slip at ½ rev. The §9 "Wire 1 dead" row is
  "slow visible smear," not "breaks."

This is purely **board-vs-board** (differential) crystal drift. It is distinct
from **image-vs-rotor** drift (motor holding RPM to ±%, not ±ppm — far larger),
which is *common-mode*: all boards share master's clock, so the image precesses
uniformly and stays coherent (§2, rotor model). Crystal drift is the only
differential error between boards, and is exactly what the sync layers null.

---

## 5. Layer 2 — Frame / flip: exactly-once, self-describing boundaries (A + C)

The boundary (`x==0`, `x==W/2`) is derived from the disciplined `x_`. The flip
must fire **exactly twice per revolution on every board** regardless of column
drift, and a single sync glitch must not latch.

### 5.1 Exactly-once flip via boundary identity (A)

Boundaries strictly alternate (0, W/2, 0, W/2 …). A single deduplicated
primitive keyed on *which* boundary makes the flip exactly-once across both the
counter-crossing path and the sync-symbol path:

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

| Board | Primary | Fallback |
|-------|---------|----------|
| Master (0) | flywheel boundary crossing | — (its counter is the reference) |
| Downstream | boundary **symbol** (authoritative) | flywheel boundary crossing |

Downstream has two flip paths, deduped by `try_flip`: a dropped *column* → the
symbol still flips; a dropped *symbol* → the crossing still flips; only losing
both in one half-rev glitches, self-healing next half-rev.

### 5.2 Self-describing boundary symbols (C)

Master width-codes the boundary on Wire 2 so each pulse names its boundary
*absolutely* — deleting the parity counter (`sync_seeded_`/`sync_at_zero_`,
`pov_segmented.h:444-450`) that today inverts permanently on one dropped pulse.

- **NARROW** (≈2 columns high) = boundary HALF (`x==W/2`).
- **WIDE** (≈4 columns high) = boundary ZERO (`x==0`).
- **EXTRA-WIDE** (≈6–8 columns, or a coded burst) = boundary ZERO **+ EPOCH**
  (Layer 3). Epoch only ever lands on a ZERO boundary.

Generation (no blocking): master raises the line at the boundary in its
timebase ISR and lowers it via a column-counted countdown — width measured in
column ticks, no extra peripheral. Decode (cold ISR, `micros()` allowed): attach
Wire 2 on **CHANGE**; on rising record `t_rise` + `x_at_rise`; on falling
classify width → Boundary (+epoch), snap `x_` to the boundary compensated for
columns elapsed during the measurement window, then `try_flip`.

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
increments once per synced flip (Layer 2), identically on all boards. This
replaces the `millis()`-gated `show<E>(120)` sequencing with epoch-counted
sequencing.

Because all boards iterate the **same** `HS_EFFECT_LIST` order, the epoch mark
need only say "advance," not "advance to N" — *provided* no epoch is ever missed
and all boards start at index 0 together.

### 6.2 Stateless vs. stateful (inherent regimes)

- **Stateless / parametric effects** (most of `animation.h`, driven directly by
  `t`): resync is *instantaneous* — set `t` from the shared counter and content
  matches. A dropped frame self-heals next frame.
- **Stateful integrators** (BZ/GS reaction-diffusion, trails, comets): cannot be
  retro-synced mid-effect (history is unrecoverable). Frequency discipline
  (Layer 1) minimizes dropped frames; the **epoch re-init bounds any stateful
  divergence to ≤ one effect**. This is inherent, not a flaw.

### 6.3 Epoch reliability

Epoch is rare (1 / 120 s) but a *missed* epoch is very visible (one segment
stuck on the previous effect for 120 s). So the epoch symbol must be more robust
than per-column timing:

> Decisions pending (§11):
> - **Redundancy:** repeat the EPOCH symbol on the next R ZERO-boundaries so a
>   single glitch can't drop it (effect changes are 960 revs apart — free).
> - **Absolute index:** optionally encode the effect index (coded burst) so a
>   late-booting or desynced board can join at the correct effect rather than
>   assuming index 0. Costs a real multi-bit symbol; decide if worth it vs. the
>   "all boot together at 0" assumption.

---

## 7. Wire 2 symbol alphabet (summary)

| Symbol | Encoding (proposed) | Meaning | Rate |
|--------|---------------------|---------|------|
| NARROW | ~2 col high | boundary HALF (`x==W/2`) | 1/rev |
| WIDE | ~4 col high | boundary ZERO (`x==0`) | 1/rev |
| EXTRA-WIDE | ~6–8 col high (or coded burst, optionally repeated R×) | boundary ZERO + EPOCH (advance effect, reset `t`) | 1/effect |

Thresholds tunable (§11). All widths ≪ 62.5 ms half-rev, so consecutive symbols
never overlap. Wire 1 continues to carry the per-column phase reference.

---

## 8. Concurrency & ISR model

ISRs on each board: **(a)** flywheel `IntervalTimer` (Layer-1 advance + pixel
pack — the hot path), **(b)** column-wire edge (Layer-1 discipline / backstop),
**(c)** Wire-2 CHANGE (Layer-2/3 decode + snap + flip + epoch).

Invariants:

1. **Hot path stays branchless / time-free.** The per-column ISR gains no
   `micros()` and no per-column `digitalRead`. Width measurement and epoch
   decode live only in the cold Wire-2 ISR (≤ 2/rev + 1/effect).
2. **Non-preemption (load-bearing, extends `pov_segmented.h:292-306`).** All
   shared state (`x_`, `last_flipped_`, the flywheel period/phase registers,
   `t_rise`, `x_at_rise`, epoch flags) is touched only by these handlers, which
   must remain mutually non-preempting (equal NVIC priority; SysTick already
   demoted). Moving the flywheel timer to a higher priority than the GPIO
   handlers would break this — call it out at the attach site. Any new shared
   field carries the same note.
3. **Flip exactly-once** via `try_flip`'s identity check + idempotent
   `advance_display` backstop.
4. **Failure asymmetry honored:** a *missed* flip (stale frame) is the only
   harmful display outcome; an *extra* flip is benign. The protocol biases
   toward flipping (two redundant paths) at every layer.
5. **Seeding:** `run()` / epoch resets `x_=0`, `last_flipped_=NONE`,
   `t=0`, flywheel period = `T0`. First boundary is HALF (`HALF≠NONE` ⇒ flips).

---

## 9. Failure analysis (post-redesign)

| Event | Layer 1 (column) | Layer 2 (flip) | Layer 3 (content) |
|-------|------------------|----------------|-------------------|
| 1 missed column edge | flywheel fills it; no drift | unaffected | unaffected |
| 1 spurious column edge | +1 col → snapped at next boundary | symbol flip unaffected | unaffected |
| 1 dropped boundary symbol | — | crossing fallback flips | — |
| 1 dropped epoch symbol | — | — | bounded by redundancy (R repeats); else 1 segment stale ≤120 s |
| 1 board renders slow (drops a frame) | — | shows prior frame 1 period | stateless: heals next frame; stateful: heals next epoch |
| Wire 1 dead | board free-runs at T0, precesses on crystal | boundary symbol still aligns flip 2/rev | epoch still aligns content |
| Wire 2 dead | flywheel still genlocked to Wire 1 | crossing fallback flips 2/rev | playlist drifts (no epoch) → degrades to per-board timing |
| Master dead | no references; boards coast then `buffer_free()` watchdog traps (`canvas.h`) | — | — |

No single-glitch event latches a permanent error at any layer. The only visible
artifacts require either two coincident losses in one half-rev (self-heals) or a
dropped epoch (mitigated by redundancy).

---

## 10. Likelihood (why the flywheel matters)

A missed column today requires one board to coalesce two edges or latch a
spurious one; ranked causes: (1) interrupt-masked foreground windows > 434 µs —
dominant on the FastLED/WS2801 bit-bang path (`FastLED.show()` masks IRQs),
largely designed out on DMA; (2) worst-case ISR overrun; (3) EMI on a motorized
spinner (spurious edges). At ~1.1 M column edges/min across 4 boards, even a
1e-6 per-edge rate is a visible glitch every minute or two. **The flywheel
removes causes (1)/(2) as column-drop sources outright** (a masked window or a
long ISR just lets the timer coast), and the boundary/epoch symbols bound the
residual from (3). This is the core robustness win of the redesign.

---

## 11. Open decisions (resolve before implementing)

1. **Layer-1 formulation:** X (wire drives, timer backstops — recommended) vs. Y
   (flywheel drives, PI discipline). Affects hot-path structure and tuning.
2. **PI gains / guard time** (if Y): `Kp`, `Ki`, slew limit, and (X) the
   missed-edge guard (`~1.5·T0`).
3. **Boundary symbol widths & threshold:** proposed NARROW 2 / WIDE 4 /
   EXTRA-WIDE 6–8 columns; confirm margin vs. worst-case sync-ISR entry latency.
4. **Snap compensation:** elapsed-column-compensated vs. hard-snap rewind.
   Recommend compensated.
5. **Epoch robustness:** repeat count R; whether to encode an absolute effect
   index for late-boot recovery.
6. **Inter-board flip skew:** master flips at the boundary; downstream on the
   symbol's falling edge (~K cols later, uniform). ≤ ~1.7 ms ≪ 62.5 ms frame —
   believed invisible; confirm or align master to its own falling edge.
7. **Share the flywheel with `pov_single`?** The single-board driver already has
   a local `IntervalTimer`; factoring a common flywheel/genlock core is natural
   but expands scope.

---

## 12. Test plan (host-testable where possible)

Following the `pov_segment_map.h` precedent (pure, host-tested index math):

- **Pure functions:** factor (width → Boundary/epoch) classification and the
  `try_flip` state machine into free functions; assert exactly-once across
  interleaved counter/symbol arrivals.
- **Flywheel/genlock sim:** drive a mock event stream (timer ticks + reference
  edges with injected drops/dups/jitter and a ppm frequency offset); assert
  (a) no column dropped under missed reference edges, (b) phase error bounded,
  (c) lock re-acquired after a burst of glitches. For Y, assert PI stability
  (no oscillation) across the gain range.
- **Layer-2 invariants:** `advance_display` count == 2/rev/board under arbitrary
  column drift; no latched inversion after a single symbol glitch.
- **Layer-3 content:** simulate 4 boards over a multi-effect playlist with
  injected per-board render-time jitter and a dropped epoch; assert all boards
  on the same `(effect, t)` except within the proven bounds (stateless: 1 frame;
  stateful: ≤ 1 effect; dropped epoch: ≤ R-bounded).
- **Determinism harness** stays green (device-only ISR changes; host renders
  unaffected). The daydream `segment_controller`/`segment_worker` reproduces the
  4-board partition and is the end-to-end coherence check before fabrication.

---

## 13. Out of scope

- **Absolute angular lock to the rotor** — needs a hall/index sensor into
  master; the system is open-loop on rotation by design (§2). Uniform precession
  under motor-speed drift is accepted.
- **Duplicate-master detection** (a peer holding the same hardware ID) — already
  out of scope (`pov_segmented.h:199-203`); needs an open-drain arbitration line.
- **Dropping Wire 1** — explicitly *not* done; kept as the genlock reference.
  Formulation Y's frequency lock is what would make it safe to revisit later.
- **`pov_single` content/rotor changes** beyond an optional shared flywheel core.
