# PHANTASM Segment Board — PCB Design Specification

Source of truth for the KiCad schematic + layout of the per-segment carrier board.
One identical PCB is built **×4**; a 2-bit solder strap selects each board's role
(segment 0 = master/conductor, segments 1–3 = flywheel slaves).

Companion documents:
- [phantasm_circuit.svg](phantasm_circuit.svg) — the wiring schematic this spec formalizes.
- Firmware contracts: [pov_segmented.h](../hardware/pov_segmented.h), [pov_sync.h](../hardware/pov_sync.h), [dma_led.h](../hardware/dma_led.h).

---

## 1. Scope & system context

- **Assembly:** 4 boards co-rotate on the arm at **480 RPM**. Mechanical robustness
  and rotor balance are first-class requirements, not afterthoughts.
- **Per board:** the PCB is a **logic controller card** — 1× Teensy 4.0 (3.3 V logic) + 1× 74AHCT125
  level shifter @ 5 V + sync conditioning. It drives one HD107S 72-px strip's **data (DI/CI)** and
  carries only **logic current (~0.15 A)** — the strip's 4.3 A power does **not** flow through the card.
- **LED power is injected off-board (§2.3).** Heavy +5 V/GND runs from the rotor busbar **directly to
  each strip's injection point** (bulk cap there), bypassing the card. The card provides DI/CI plus a
  **dedicated low-current signal-GND** that ties to the strip's ground at the **strip (load) end** —
  keeping the 4.3 A return's IR drop out of the data reference (R-SI-2).
- **Crosses the slip ring:** **+5 V** and **GND** only, on **two separate runs** — a heavy LED-power
  feed to the strips and a light logic feed to the cards. SYNC and SPI stay entirely on-rotor.
- **Shared between boards:** one **SYNC** wire (multidrop bus) + the power rails.
  **Nothing high-speed is shared** — each board's 24 MHz SPI data is local to its own strip.
- **Inter-board sync cable:** **Belden 8451** (shielded twisted pair, 2× 22 AWG + foil +
  drain), daisy-chained board→board. One conductor = **SYNC**, the other = its **GND return**
  (tight pair, small loop); the overall **shield is grounded at a single point — the master**
  (§4.4, §7). 5 V/GND power is distributed on the separate heavy rotor rail/busbar, **not** on
  this pair.
- **Assembly:** **partial PCBA** — low-mass SMD parts reflow-placed by the house; heavy/mechanical
  through-hole parts hand-soldered. See §11.

### 1.1 Teensy 4.0 pin map (fixed by firmware — do not reassign)

| Teensy pin | Net | Direction | Notes |
|---|---|---|---|
| 11 (MOSI) | LED **DATA** → DI | out | → '125 ch A → 33 Ω → strip DI |
| 13 (SCK)  | LED **CLK** → CI  | out | → '125 ch B → 33 Ω → strip CI |
| 3  | **FRAME_SYNC** | in/out | OUTPUT on master, INPUT on slaves (mutually exclusive); drive via '125 ch C, receive via divider |
| 5  | **MASTER_EN** | out | LOW on master; gates '125 ch C `/OE`; **R_MEN 10 kΩ pull-up → 3V3** (boot-safe disable, R-LS-5) ([pov_segmented.h:112](../hardware/pov_segmented.h#L112)) |
| 21 | **ID0** | in (PULLUP) | strap bit 0; ground = bit set |
| 22 | **ID1** | in (PULLUP) | strap bit 1; ground = bit set |
| 23 | **ID2** | in (PULLUP) | strap bit 2 — **reserved for 8-seg (N=8)**; unread at N=4 |
| VIN | **+5 V** | power in | board + strip rail |
| GND | **GND** | power | common return |
| 3V3 | **+3.3 V** | power out | on-board regulator; sources **R_MEN** (mandatory) + R_ID0 pull-up (opt) + J4 |

SPI clock = **24 MHz** ([pov_segmented.h:151](../hardware/pov_segmented.h#L151)). Lay out the fast nets to be
clean to **≥30 MHz** so headroom exists.

---

## 2. Power

### 2.1 Budget

**Controller card (on-board):**

| Load | Peak | Notes |
|---|---|---|
| Teensy 4.0 | ~0.1 A | |
| 74AHCT125 + dividers | < 5 mA | |
| **Card total (peak)** | **≈ 0.15 A** | low-current board: 1 oz copper, small connectors |

**LED power (off-board, §2.3):** HD107S 72 px ≈ **4.3 A** peak (72 × 60 mA, full white), delivered to
the strip by the separate power harness — **never through the card**. Aggregate ≈ **17 A** across 4
strips at the slip ring (R-PWR-9).

### 2.2 Requirements — controller card

- **R-PWR-1** Logic power enters on a **small** TH connector (J1), `+5 V`/`GND`, **~0.5 A** rating
  (card draws ~0.15 A). This is the **light logic feed**, separate from the LED-power harness (§2.3).
- **R-PWR-2** Card copper carries only ~0.15 A → **1 oz is sufficient**, no heavy pours. Strip 5 V/GND
  sizing lives off-board with the power harness (§2.3, R-PWR-10).
- **R-PWR-3 — Bulk moves off-board.** The **1000 µF strip bulk is at the off-board injection point**
  (§2.3, R-PWR-11), **not** on the card. The card keeps only a small **C_IN (≥100 µF)** on its logic feed.
- **R-PWR-4** **0.1 µF (C_DEC)** decoupling at the Teensy VIN pin and the U1 Vcc pin, within 3 mm,
  short via to the ground plane; plus **C_LF** on the logic rail (R-PWR-5).
- **R-PWR-5 — Logic rail filter.** Two nodes now (no strip rail on the card): **+5V_RAW** (J1 → F1 →
  Q_REV) → ferrite bead (≈600 Ω @ 100 MHz) → **+5V_LOGIC** (Teensy VIN + U1 Vcc + C_DEC + C_LF). The
  bead carries only the **~0.15 A logic branch**, isolating it from conducted noise on the shared rotor
  rail. **Damp the bead-LC** — a bead into a low-ESR ceramic is a high-Q tank that *peaks* noise at f₀:
  a **small series R (default 1–2 Ω, R_LF) ahead of C_LF**, or a few µF of **tantalum/ESR cap** in
  parallel, or a **lossy bead**. Keep R_LF low (≤0.3 V at 0.15 A). **C_LF DC-bias:** a 1206 10 µF X5R
  derates to ~6 µF at 5 V → use a **22 µF** part (or the tantalum option, which also damps).
- **R-PWR-6** With the 1000 µF off-board (§2.3), the card's only electrolytic is **C_IN (≥100 µF)** —
  retain it with RTV and place it (and the connectors/Teensy) toward the hub for balance (R-MECH-2).
- **R-PWR-7 — Reverse protection collapses to logic-only.** The card carries ~0.15 A, so a **small
  series Schottky or P-FET (Q_REV) on the logic feed** is enough — the 4.4 A ideal-diode problem (and
  its body-diode-orientation / SOA subtleties) is **gone with the LED current**. Still **key/polarize
  J1** (a reversed feed destroys the Teensy + '125). The **LED-power harness gets its own
  protection/keying off-board** (§2.3).
- **R-PWR-8 — Per-card overcurrent.** The card's ~0.15 A logic feed wants only a **small fuse/PTC (F1,
  ~0.5–1 A)** at J1, or documentation that it's covered upstream. (Strip overcurrent lives with the
  power harness — §2.3, R-PWR-12.)
- **R-PWR-9 — Aggregate + inrush (upstream).** The slip ring + upstream feed still carry **≈17 A peak**
  (4 × 4.3 A) on the **heavy LED-power run**; **soft-start / inrush limiting lives upstream** at the
  PSU / slip-ring distribution. The big ~1 mF bulk is now at the strip injection points (§2.3,
  R-PWR-11), off the cards, so per-card inrush is just the small C_IN.

### 2.3 Power injection (off-board)

The strips' 4.3 A is delivered by a **separate heavy harness from the rotor busbar straight to each
strip**, bypassing the controller cards. This is what makes the card low-current and electrically
quiet (R-SI-2).

- **R-PWR-10 — Heavy feed sizing & cable.** Size the busbar/harness 5 V/GND for **4.3 A per strip** at
  < 250 mV drop and acceptable rise (the old on-board 2 oz pour requirement, relocated). Use **16–18 AWG
  or a busbar** — **not** the 22 AWG Belden 8451 sync pair (≈0.28 V drop over ~2 ft at 4.3 A, and it's a
  signal cable). The heavy feed **may run alongside the shielded sync 8451** (the sync's shield + twist +
  RC + deglitch tolerate it), but it stays a **separate conductor** — never share the sync pair's
  conductors or shield for power, and cross power/sync at ~90° where they meet. Its own
  reverse-protection/keying lives here.
- **R-PWR-11 — Bulk at the injection point.** Put the **1000 µF (2200 µF for margin)** bulk
  electrolytic **at each strip's 5 V/GND injection point**, `+` to 5 V — right where the 4.3 A transient
  is sourced. This point is also the **ground star**: the card's signal-GND meets the strip ground here
  (R-SI-2, §6).
- **R-PWR-12 — Per-strip fusing.** State whether each strip's 4.3 A feed is fused at injection or
  protected upstream — don't leave it implicit.

---

## 3. Level shifter — 74AHCT125 (U1, SOIC-14, Vcc = 5 V)

Teensy 4.0 is **3.3 V and NOT 5 V tolerant**. AHCT reads 3.3 V as a valid TTL HIGH and
drives a clean 5 V output — the correct in-spec 3.3 → 5 V up-shifter.

| Ch | Input (3.3 V) | Output (5 V) | `/OE` | Series R | To |
|---|---|---|---|---|---|
| A | Teensy 11 DATA | DI | tied LOW (always on) | **33 Ω** | strip DI (J2) |
| B | Teensy 13 CLK | CI | tied LOW (always on) | **33 Ω** | strip CI (J2) |
| C | Teensy 3 SYNC-OUT | SYNC bus | **Teensy 5 (MASTER_EN)** | **100 Ω** | SYNC bus (J3) |
| D | **tie input → GND** | NC | **tie /OE → Vcc** | — | unused (disabled) |

- **R-LS-1** Series terminations sit **at the '125 output pin** (source termination): 33 Ω on
  DATA/CLK, 100 Ω on SYNC-OUT.
- **R-LS-2** Channel C `/OE` is driven by Teensy **pin 5**. Only the master (pin 5 = LOW) enables
  SYNC-OUT; DATA/CLK channels have `/OE` tied LOW so they are always active on every board.
- **R-LS-3** Tie the **unused channel D input to GND** *and* its **`/OE` to Vcc** (output disabled).
  Never float a CMOS input or control pin — `/OE` is itself a logic input.
- **R-LS-4** **U1 is SMD (SOIC-14), reflow-placed by the PCBA house** (§11). SOIC is strictly better
  than DIP here — lower profile, low mass, no socket/shunt to sling off at 480 RPM, so no RTV needed.
  Part = **SN74AHCT125** (e.g. SN74AHCT125DR). The 5 V-tolerance caveat below is about signal
  *direction*, not package.
- **R-LS-5 — Default the sync driver disabled at boot.** Ch C `/OE` is MASTER_EN (Teensy pin 5),
  which **floats from power-on until `run_show()` drives it** ([pov_segmented.h:236-237](../hardware/pov_segmented.h#L236-L237));
  a floating `/OE` that settles LOW briefly enables a slave's bus driver — the transient phantom-master
  hazard. Fit a **pull-up R_MEN (10 kΩ) on pin 5 → 3V3** (not 5 V — keeps pin 5 safe; 3.3 V is a solid
  AHCT TTL HIGH = `/OE` disabled). Every board then boots with its sync driver **off**, enabled only
  when firmware actively asserts master — consistent with the DNP ID0 pull-up (R-ID-3). *Companion
  (optional):* weak pulls on the ch A/B inputs (Teensy 11/13) to define DATA/CLK during the same
  boot window, since those '125 inputs float too — **intentionally not in the base BOM** (APA102-class
  strips latch no garbage from a brief boot transient; add only if a startup flash is objectionable).

> **Why a divider for sync receive (not a 5 V-tolerant buffer)?** Pin 3 isn't 5 V-tolerant and the
> bus swings to 5 V, so the receive path must down-shift. Now that SMD is in play a 74LVC125A @3.3 V
> *could* do it, but that spends another buffer + channel; the passive R1/R2 divider in §4 is simpler
> and already required (its R2 is pin 3's pull-down). **Keep the divider.**

---

## 4. SYNC bus (the only shared signal)

Sync is a **low-rate symbol stream**, not a clock: 2 boundary marks/revolution + rare epoch/beacon
([pov_sync.h:20](../hardware/pov_sync.h#L20)). Pulse pitch ≈ 868 µs; edges are ≥100 µs apart and pass a
~100 µs firmware glitch filter ([pov_sync.h:105](../hardware/pov_sync.h#L105)).

**Topology: single source-terminated multidrop bus.** The master's one '125 channel drives the
shared wire; all four boards (incl. master) tap it through high-impedance dividers (~25 kΩ each);
with the 10 kΩ idle pull-down the total bus load is **≈1.3 mA**. **No star, no second buffer** — see §4.3.

### 4.1 Drive path (master only)
- Teensy 3 → '125 ch C → **100 Ω** series → SYNC bus (J3). Gated by `/OE` = pin 5.
- **Source termination is valid because the master sits at a bus _end_** (daisy master→1→2→3): the
  100 Ω damps the outgoing edge and absorbs the reflection back at the driver. If the daisy order ever
  changes so the master is mid-bus, single-source-termination no longer applies — re-evaluate.

### 4.2 Receive path (every board)
- SYNC bus → **R1 = 10 kΩ** → node → **R2 = 15 kΩ** → GND. Node ≈ **3.0 V** → Teensy pin 3.
  **R2 = 15 kΩ is the default** (not 18 kΩ): pin 3 is **not** 5 V-tolerant, so this divider is the
  hard clamp. 10k/15k gives 3.0 V nom and **3.15 V at a hot 5.25 V rail** — under the 3.3 V max
  *independent of V_OH*, still ≫ the ~2.3 V Teensy V_IH. (The legacy 18 kΩ from the SVG yields
  5.25 × 18/28 ≈ **3.38 V**, over the 3.3 V max and rescued only by V_OH droop — don't ship it on a
  not-5V-tolerant pin.)
- **R-SYNC-1** R2 doubles as pin 3's defined-state pull-down. **Do not add any other pull-down.**
- **R-SYNC-2** One **10 kΩ idle pull-down** on the bus itself (keeps the bus defined when no
  master is driving / during boot).
- **R-SYNC-3 — Node RC filter (populate by default).** Fit **C_SYNC ≈ 220 pF** at the pin-3 divider
  node. With R_th = R1‖R2 ≈ 6.0 kΩ this gives **RC ≈ 1.3 µs** — negligible against the 434 µs column
  and 868 µs pulse pitch, but it attenuates sub-µs BLDC/LED spikes *before* they cross the GPIO
  threshold. The firmware glitch filter ([pov_sync.h:105](../hardware/pov_sync.h#L105)) gates edge *spacing*
  (≥100 µs), not *amplitude* — a single fast spike that crosses threshold still registers as a real
  edge and can corrupt a burst count. In this known-noisy environment, **default-populate**; keep the
  value (100 pF–1 nF) tunable. Pair with receiver hysteresis (R-SYNC-7).
- **R-SYNC-4** Master ignores its own echo in firmware; no hardware action needed.
- **R-SYNC-7 — Receiver hysteresis (implemented in firmware).** A slowed edge (R-SYNC-3) into a
  non-hysteresis i.MX RT pad can double-trigger near threshold, so pad **HYS** is enabled on pin 3
  after `pinMode` in [pov_segmented.h](../hardware/pov_segmented.h) (`portControlRegister(PIN_FRAME_SYNC) |=
  IOMUXC_SW_PAD_CTL_PAD_HYS`). No PCB action — listed so layout/firmware stay paired with C_SYNC.

### 4.3 Decision: bus vs. star (resolved — bus)
The 74AHCT125 sources **±8 mA** (datasheet, AHCT class) into the bus load — four ~25 kΩ receive
dividers plus the 10 kΩ idle pull-down, **≈1.3 mA total** — so DC fan-out is a non-issue (≈6× margin).
A star would (a) require a **second
'125 on the master** — its quad already spends 3
channels on DATA/CLK/SYNC-OUT — and (b) break the identical-board design, all to clean up edges on
a signal that is 100 µs-paced and software-deglitched. Inter-board skew on a ~1–2 m bus is a few ns
versus a 434 µs column — negligible. **Keep the single bus.** Revisit only if sync rate is ever
raised by orders of magnitude or false triggers are observed on a scope.

### 4.4 Inter-board cable & shield (Belden 8451)
The bus runs board→board on **Belden 8451** (STP, 2× 22 AWG + foil shield + drain). The two
conductors carry the signal pair; the drain handles the screen:
- **R-SYNC-5 — Pair assignment.** Land one conductor on `SYNC` and the other on `GND` at each
  tap (J3A/J3B, §6) so signal and return travel as the twisted pair — small loop area, low
  pickup/emission next to the BLDC and the 4 A LED return (§7). Land each cable's `GND` conductor
  on the plane *adjacent to its own* `SYNC` pin.
  **Accepted DC ground loop:** this cable-GND conductor *and* the heavy rotor busbar both common all
  board GNDs (master→cable→slave→busbar→back), so the twisted pair is **not** an isolated return —
  for a 100 µs-paced, divided, RC-filtered, software-deglitched signal that is acceptable and is the
  chosen tradeoff. If it ever proves noisy, break it with a **small series R or ferrite in the
  cable-GND tap**. (A stated choice, like R-PWR-8/9 — not an oversight.)
- **R-SYNC-6 — Single-point shield.** The cable shield/drain is grounded at **exactly one node —
  the master board** — through jumper **JP_SHLD** (§6). The drain passes through every board's
  connector (bridged in↔out) for continuity, but **JP_SHLD is DNP on slaves**, so the continuous
  shield has one ground reference and forms **no shorted turn** around the rotor. Heatshrink any
  unused drain tail; never ground the shield at two boards.
- **R-SYNC-8 — Optional bus transient clamp.** Pin 3 is well protected (R1 = 10 kΩ limits its clamp
  current), but the '125 output sits directly on a 1–2 m cable in a BLDC field with no clamp. A small
  **TVS / clamp diode (D_BUS) at the bus tap** is cheap insurance against inductive/ESD transients on
  that exposed low-impedance node. Optional (DNP footprint).

---

## 5. Hardware ID strap (pins 21 / 22 / 23)

Straps are `INPUT_PULLUP`; firmware reads **log2(N) of them**, inverts, and masks:
`segment_id = (~raw) & (N−1)` ([pov_segmented.h:411](../hardware/pov_segmented.h#L411)). **Grounding a strap
sets its bit; all straps open = ID 0 = master.** The build is **N = 4** (reads ID0/ID1, pins 21/22); the
board also breaks out **ID2 (pin 23)** so a recompiled **N = 8** build decodes segments 4–7 by the same
scheme. N must be a power of two ≤ 8.

### 5.1 Truth table (built config, N = 4)

| Segment | Role | pin 21 (ID0) | pin 22 (ID1) |
|---|---|---|---|
| **0** | master / conductor | open | open |
| **1** | slave | **→ GND** | open |
| **2** | slave | open | **→ GND** |
| **3** | slave | **→ GND** | **→ GND** |

At N = 4, **ID2 (pin 23) stays open** and is unread (`ID_STRAPS = log2(4) = 2`). For an **N = 8** build,
ground **ID2** to add bit 2 → segments 4–7 follow the same "ground to set the bit" pattern; firmware
widens the mask to `& 7` automatically.

### 5.2 Implementation
- **R-ID-1** Break out **pins 21, 22 and 23 (ID0/ID1/ID2)** to pads, each with a **GND pad/via within
  ~5 mm** — ID2 reserved for the 8-segment build. The role is set by a **short soldered wire-link** (or
  a stuffed 0 Ω / closed solder-jumper) from the pin pad to its GND pad, per the table. "Open" straps
  get nothing — the on-die pull-up holds HIGH.
- **R-ID-2** **The grounded link is the load-bearing connection.** Its failure mode is an *open*,
  which reads HIGH → inverts toward ID 0 → **elects a phantom second master** and causes sync-bus
  contention ([pov_segmented.h:379-385](../hardware/pov_segmented.h#L379-L385)). Use a fully soldered link
  (not a removable header/shunt — it can sling off at speed). Keep the link short and tack with RTV.
- **R-ID-3** **No external pull-downs** on the straps (would invert the decode). **Pull-ups are
  safe**: optionally provide a **DNP 10 kΩ pull-up footprint on ID0 → 3.3 V** to harden the
  open/master state against a cold strap (firmware author's suggestion, [pov_segmented.h:379-381](../hardware/pov_segmented.h#L379-L381)).
- **R-ID-4** **Silkscreen the truth table** beside the straps and a large **board-ID digit field**,
  so each board's role is readable during balancing/assembly. Mark "MASTER = both OPEN."

---

## 6. Connectors (all through-hole, for rotor mechanical strength)

| Ref | Function | Pins | Type / rating | Pinout |
|---|---|---|---|---|
| **J1** | Logic power in (light feed) | 2 | 0.1″ TH header or small JST, **~1 A** | `+5 V`, `GND` |
| **J2** | Strip **signal** out | 3 | 0.1″ TH header | `DI` (DATA, post-33 Ω), `CI` (CLK, post-33 Ω), `SIG_GND` |
| **J3A** | SYNC daisy — **in** | 3 | 0.1″ TH header | `SYNC`, `GND`, `SHLD` (one Belden 8451) |
| **J3B** | SYNC daisy — **out** | 3 | 0.1″ TH header | `SYNC`, `GND`, `SHLD` (one Belden 8451) |
| **JP_SHLD** | shield ground point | — | solder jumper / 0 Ω | drain net → GND; **stuff master only** |
| **J4** (opt) | Debug/serial breakout | 3–4 | 0.1″ TH header | `3V3`, `GND`, `pin 5`, TX/RX |

- **R-CON-1 — J2 is signal-only; tie SIG_GND at the strip (load) end.** J2 carries **DI/CI + a
  dedicated SIG_GND** at logic current — a plain 0.1″ header is fine (no power pins). **SIG_GND is the
  card's logic GND, run alongside DI/CI and landed on the strip's GND pin at the strip input** — *not*
  at the slip ring. The strip's GND pin is the **ground star** where this low-current reference meets
  the heavy LED return from the power harness (R-PWR-11); the 4.3 A return's IR drop then sits upstream
  of the junction and never appears across the data pair (R-SI-2/3). The **1000 µF bulk lives at that
  injection point off-board** (R-PWR-11), not at J2.
- **R-CON-2** The 33 Ω series resistors sit **between the '125 outputs and J2** DI/CI pins (so the
  termination is at the source, ahead of the cable to the strip).
- **R-CON-3** **J3A/J3B each terminate one Belden 8451** (SYNC + its GND return + drain). On-board:
  the two `SYNC` pins are **bridged** (the bus tap that feeds the divider / master source), each
  `GND` lands on the plane next to its own `SYNC` pin so each pair's loop stays small (§7), and the
  two `SHLD` pins are **bridged to a shared drain net**. **End boards (master, seg-3) populate only
  one connector;** mid-boards populate both for through-daisy.
- **R-CON-4** All connectors TH with adequate annular ring/strain relief; orient so the mating
  cables route toward the hub, not radially outward into the spin.
- **R-CON-5** **Shield single-point ground:** the bridged drain net reaches GND only through
  **JP_SHLD**, stuffed **on the master board only**; slaves leave it open (§4.4, §7 R-SI-8).
- **R-CON-6** The 8451 daisy is co-rotating at 480 RPM: **strain-relieve it within ~2 cm of
  J3A/J3B** and either **solder the tails or use a latching connector** — a bare friction 0.1″
  shunt can walk off at speed (cf. R-ID-2).

---

## 7. Signal integrity & grounding (noisy environment: motor + 24–30 MHz LED SPI)

This is the section that matters most for a board spinning next to a BLDC. With LED power now off the
card (§2.3), the on-board ground is quiet by construction (R-SI-2); the remaining work is the fast SPI
nets, the sync line, and the load-end signal-ground tie.

- **R-SI-1 — Reference plane.** Minimum **2-layer with a continuous bottom GND pour** under all
  signals; **4-layer (SIG / GND / PWR / SIG) preferred**. Never route a fast net over a plane split.
- **R-SI-2 — LED return is off-board (the #1 item, now solved by construction).** Because the 4.3 A
  strip 5 V/GND never touch the card (§2.3), the pulsed LED di/dt is **inherently off** the
  Teensy/'125/sync ground reference — the on-board ground is quiet by design, no on-board star-pour
  needed. The **one remaining rule:** the card's **SIG_GND** ties to the strip ground at the **strip
  (load) end** (R-CON-1), so the heavy return's IR drop stays upstream of the data reference. Do **not**
  also bond logic GND to the heavy power-GND at the card / slip-ring end — that reintroduces the drop
  across the data pair.
- **R-SI-3 — Fast nets short & referenced.** DATA (pin 11) and CLK (pin 13) Teensy→'125→J2 runs:
  keep **short, direct, over unbroken ground**, no stubs, roughly length-matched to each other.
  33 Ω source termination at the '125 (R-LS-1). At 30 MHz keep these under ~6 cm where practical.
  **Off-board J2→strip harness:** carries **DI/CI + SIG_GND only** (no power). Keep it **short
  (≤ ~15 cm)** with **DI/CI paired against SIG_GND**, which is the data return and lands on the strip's
  GND pin (R-CON-1). Source termination tames the near end; a long lead into the strip's load still
  degrades the 24 MHz edges regardless of the series R.
- **R-SI-4 — Decoupling discipline.** 0.1 µF caps (R-PWR-4) on the **same side** as their IC, < 3 mm
  from the pin, with a dedicated via to the plane (minimize the Vcc→cap→GND loop).
- **R-SI-5 — Sync conditioning.** 100 Ω source series (R-LS-1) + 10 k idle pull-down (R-SYNC-2) +
  optional RC (R-SYNC-3). Route **SYNC adjacent to a ground return** (small loop), keep it **away from
  and orthogonal to** the CLK net and any motor/ESC leads.
- **R-SI-6 — Logic-rail filtering.** The ferrite-bead filter (R-PWR-5) isolates the card's logic supply
  from conducted noise on the shared rotor rail; the heavy LED feed is a separate run (§2.3).
- **R-SI-7 — EMI hygiene.** Minimize all current-loop areas; if 4-layer, stitch the perimeter and add
  ground vias around the fast nets. Keep motor wiring entirely off this board.
- **R-SI-8 — Shield single-point ground.** The Belden 8451 overall shield is bonded to GND at **one
  node only — the master (JP_SHLD, §6).** Grounding it at both ends closes a shorted turn around the
  rotor (circulating current, noise antenna); leaving it fully floating loses the screen. The drain is
  bridged through every daisy connector and referenced exactly once.

---

## 8. Mechanical / rotor (480 RPM)

- **R-MECH-1 — Mass-based package split (§11).** Heavy parts (electrolytics, connectors, Teensy) are
  **through-hole** for joint strength + RTV bonding under centrifugal load; low-mass parts (U1 SOIC,
  all passives) are **SMD** — their pad-only joints carry negligible load at 480 RPM.
- **R-MECH-2** Mounting holes sized to the rotor hardware; place heavy parts (electrolytics,
  connectors, Teensy) **symmetrically and as near the hub as routing allows** for balance.
- **R-MECH-3** Provide pad area / clearance to **bond C_IN and any tall part** with RTV (the 1000 µF
  bulk is off-board now, §2.3).
- **R-MECH-4** Keep board outline and component height within the arm's swept envelope; chamfer/relieve
  the leading edge if it sees airflow.
- **R-MECH-5** **1 oz copper suffices** — the card carries only ~0.15 A (the 4.3 A LED current is
  off-board, §2.3).
- **R-MECH-6 — Board width ≤ 35 mm.** The card co-mounts along the rotor arm; the outline is a narrow
  strip with **width hard-capped at 35 mm** (raised from 30 mm so the small SMD parts pack *beside* the
  Teensy and the board stays short), and length is the free dimension to be minimised. The floor part is
  the **Teensy 4.0 (17.8 mm wide)**; at 35 mm there is ~17 mm of width left alongside it for the SMD
  parts (the D16 mm bulk that used to set the width is off-board, §2.3). 2-D-pack the parts to the
  shortest length within the width cap (generated board ≈ **54 × 35 mm**); keep heavy parts (C_IN,
  connectors, Teensy) toward the **hub end** for balance (R-MECH-2).

---

## 9. Bill of materials (per board)

Assembly column: **SMD** = reflow-placed by the PCBA house (top side, §11); **TH** = through-hole,
hand-soldered by you.

| Ref | Part | Value / PN | Pkg | Asm |
|---|---|---|---|---|
| U_MCU | Teensy 4.0 | — | 2× 0.1″ header rows | TH |
| U1 | Level shifter | SN74AHCT125 (…DR) | SOIC-14 | **SMD** |
| C_IN | Electrolytic | ≥100 µF | radial TH | TH (RTV) |
| C_LF | Ceramic | 22 µF (≥10 µF *effective* after DC-bias derate) | 1206 | **SMD** |
| R_LF (opt) | Bead-LC damping | 1–2 Ω default (≤4.7 Ω; or lossy bead / ESR cap instead) | 0805 | **SMD** |
| C_DEC1,2 | Ceramic | 0.1 µF | 0603/0805 | **SMD** |
| C_SYNC | Ceramic | 220 pF (default-populated; 100 pF–1 nF tunable) | 0603 | **SMD** |
| R_D1, R_D2 | Series term | 33 Ω | 0603/0805 | **SMD** |
| R_S | Sync source | 100 Ω | 0603/0805 | **SMD** |
| R1 | Divider top | 10 kΩ | 0603 | **SMD** |
| R2 | Divider btm | **15 kΩ** (legacy 18 kΩ — see §4.2) | 0603 | **SMD** |
| R_PD | Bus idle pull-down | 10 kΩ | 0603 | **SMD** |
| R_ID0 (opt) | ID0 pull-up | 10 kΩ → 3V3 | 0603 | **SMD** DNP |
| R_MEN | MASTER_EN boot pull-up | 10 kΩ → 3V3 | 0603 | **SMD** |
| FB | Ferrite bead | ≈600 Ω @ 100 MHz, logic branch (~0.15 A) | 1206 | **SMD** |
| Q_REV | Reverse protect (logic) | small Schottky or P-FET, ~1 A | SOD-123 / SOT | **SMD** |
| D_BUS (opt) | Bus transient clamp | TVS / clamp diode | 0603 | **SMD** DNP |
| F1 | Fuse / PTC (logic) | ~0.5–1 A | 1206 / TH | **SMD** |
| J1 | Logic power in | 2-pin ~1 A (0.1″ / JST) | TH | TH |
| J2 | Strip signal out | 3-pin 0.1″ (DI/CI/SIG_GND) | TH | TH |
| J3A, J3B | SYNC daisy in / out | 2× 3-pin 0.1″ | TH (one Belden 8451 each) | TH |
| JP_SHLD | Shield ground jumper | 0 Ω / solder jumper | 0603 or SJ pad | hand, **master only** |
| J4 (opt) | Debug | 3–4-pin 0.1″ | TH | TH |
| — | ID strap links | wire / 0 Ω | — | hand, per §5 |
| — | Inter-board cable | Belden 8451 (STP 2×22 AWG + drain) | per arm/hub run | hand |
| — | **LED power harness** (off-board) | heavy 5 V/GND, **4.3 A/strip** (§2.3) | busbar / wire | hand |
| C_BULK | **Injection bulk** (off-board) | 1000 µF / ≥10 V (2200 µF margin) | radial TH | at strip, hand |

---

## 10. Net summary

| Net | Members |
|---|---|
| **+5V_RAW** | J1 → F1 → Q_REV (logic reverse-protect) → FB in |
| **+5V_LOGIC** (post-bead) | FB out, R_LF, C_LF, C_IN+, Teensy VIN, U1 Vcc, C_DEC1/2 |
| **3V3** | Teensy 3V3 pin → R_ID0 top (opt), R_MEN top, J4 |
| **GND** | single quiet logic plane: J1 GND, Teensy GND, U1 GND, C_IN/C_LF/C_DEC −, J3A/J3B GND, R2/R_PD bottoms, strap GND, JP_SHLD (master), SIG_GND → J2 |
| **DATA** | Teensy 11 → U1 chA in; U1 chA out → R_D1(33 Ω) → J2 DI |
| **CLK** | Teensy 13 → U1 chB in; U1 chB out → R_D2(33 Ω) → J2 CI |
| **SIG_GND** | logic GND → J2 pin 3 → strip GND pin (**load-end star**, R-CON-1) |
| **SYNC_OUT** | Teensy 3 → U1 chC in; U1 chC out → R_S(100 Ω) → SYNC bus |
| **SYNC bus** | J3A/J3B SYNC (bridged), R1 top, R_PD top, D_BUS (opt), U1 chC out (via R_S) |
| **SYNC_RX** | SYNC bus → R1 → node (≈3.0 V; R2 to GND, C_SYNC) → Teensy 3 |
| **SHIELD** | J3A SHLD, J3B SHLD (bridged drain) → JP_SHLD → GND (**master only**) |
| **MASTER_EN** | Teensy 5 → U1 chC `/OE`; R_MEN (10 kΩ) pull-up → 3V3 |
| **ID0** | Teensy 21 → strap pad (→GND per role); opt R_ID0 to 3V3 |
| **ID1** | Teensy 22 → strap pad (→GND per role) |

> **Pin 3 is one physical node.** `SYNC_OUT` and `SYNC_RX` are the two roles of the *same* Teensy
> pin-3 net: the master drives it (pin 3 = OUTPUT → ch C), slaves sense the bus through the divider
> (pin 3 = INPUT). The roles are **mutually exclusive per board** ([pov_segmented.h:239-242](../hardware/pov_segmented.h#L239-L242));
> on the master the divider is bypassed because pin 3 is a driven output.

> **Strip +5 V / GND are _not_ on the card.** They come from the off-board power harness
> (R-PWR-10/11). The card's only tie to the strip ground is **SIG_GND at the strip's GND pin (load
> end)** — that single junction is the ground star (R-SI-2). Never also common logic GND to the heavy
> power-GND at the card end.

---

## 11. Assembly strategy — partial PCBA (SMD pre-assembly + TH hand-build)

The board is built as a **partial PCBA**: the assembly house reflow-places the low-mass SMD parts;
you hand-solder the heavy / mechanical through-hole parts. This keeps the rotor-balance intent of
§6/§8 (heavy parts TH + RTV-bonded) while removing the tedium and cold-joint risk of hand-soldering
a dozen tiny passives ×4 boards. The original "all through-hole" framing (old §3 note) was about
enabling *hand* assembly, not a hard electrical or mechanical constraint — partial PCBA supersedes it.

### 11.1 Package split

| Pre-assembled by house (SMD, top side) | Hand-soldered by you (TH) |
|---|---|
| U1 74AHCT125 (SOIC-14) | Teensy 4.0 (header rows) |
| R_D1/R_D2, R_S, R1, R2, R_PD, R_LF (0603/0805) | J1, J2, J3A/J3B, J4 connectors |
| C_DEC1,2 (0.1 µF), C_LF (22 µF) | C_IN (≥100 µF) electrolytic — RTV-bonded |
| FB ferrite (1206), Q_REV (Schottky/P-FET), F1 (fuse) | ID strap links, JP_SHLD — manual |
| C_SYNC (220 pF, **populated** — noise filter, R-SYNC-3) | C_BULK + LED power harness — **off-board** (§2.3) |
| R_ID0, D_BUS — **DNP** (footprints only) | — |

Rationale: SMD pad joints carry negligible centrifugal load at these masses; electrolytics rely on
pads alone would be a bad trade at 480 RPM, so they stay TH + RTV (R-PWR-6 / R-MECH-3).

### 11.2 Order/layout action items (do before fab)

- **R-ASM-1 — SMD footprints + part numbers.** Swap every pre-assembled symbol from the TH packages
  named in §3/§9 to its SMD footprint and assign a **manufacturer / LCSC part number** in the
  schematic (the PCBA house places from these). Passives are all JLC **basic parts**.
- **R-ASM-2 — All SMD on one side.** Place every reflow part on the **top** side; double-sided
  assembly carries a surcharge. Teensy, electrolytics, and connectors go on whichever side suits hand
  assembly (typically bottom or the same top edge), but **no SMD on the second side**.
- **R-ASM-3 — Check U1 part class.** SN74AHCT125 may be an **extended part** on JLCPCB (one-time
  feeder/setup fee); confirm at quote time. The passives are basic (no fee).
- **R-ASM-4 — DNP discipline.** **C_SYNC IS populated** (220 pF — the BLDC/LED noise filter, R-SYNC-3);
  do *not* omit it. **R_ID0 and D_BUS are Do-Not-Populate** (footprints present, parts omitted).
  JP_SHLD is hand-stuffed on the master only.
- **R-ASM-5 — Tight SMD placement.** With SMD, honor R-PWR-4 / R-SI-4 easily: C_DEC within 3 mm of
  U1 Vcc and Teensy VIN; C_LF at the bead output. Keep U1 and its decoupling clustered.
- **R-ASM-6 — Hand-solder thermals.** TH pads on the GND pour (connectors, electrolytics) need
  **thermal-relief spokes** so a hand iron can wet the joint against the copper ground plane.
- **R-ASM-7 — Sever Teensy VUSB↔VIN.** Each board powers Teensy VIN from the rotor rail and J4
  exposes serial for USB debug. **Cut the VIN/VUSB pad on every Teensy 4.0** so a live rail can't
  back-feed a USB host's VBUS (or vice-versa) during bring-up/debug/balancing; VIN is then fed only
  from the board rail. Known Teensy gotcha — **mandatory** given J4.
