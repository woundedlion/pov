/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_segmented.h
 * @brief Multi-Teensy segmented POV display driver for Phantasm.
 *
 * Phantasm uses N Teensys (typically 4), each controlling a contiguous
 * segment of LEDs on a single arm of the POV spinner.  Two Teensys sit
 * on each arm — one driving the top half (N pole toward junction) and
 * one driving the bottom half (S pole toward junction).
 *
 * Physical strip layout (N=4, S=288, H=144):
 *
 *   Arm A:
 *     Segment 0 (top):    LED 0 at N pole (y=0)   → LED 71 at junction (y=71)
 *     Segment 1 (bottom):  LED 0 at S pole (y=143) → LED 71 at junction (y=72)  ← reversed
 *
 *   Arm B (x offset by W/2):
 *     Segment 2 (top):    LED 0 at N pole (y=0)   → LED 71 at junction (y=71)
 *     Segment 3 (bottom):  LED 0 at S pole (y=143) → LED 71 at junction (y=72)  ← reversed
 *
 * Each Teensy reads a hardware ID from GPIO pins at boot to determine
 * which segment it owns.  One wire connects all boards
 * (docs/phantasm_frame_sync_spec.md):
 *
 *   Sync wire: segment 0 (the master/conductor) emits count-coded symbol
 *   bursts — two boundary marks per revolution, an epoch mark once per
 *   effect, and a mid-revolution index beacon.  Every board generates its
 *   own columns from a local flywheel timebase (position derived from the
 *   free-running cycle counter, never from counting timer interrupts); the
 *   symbols snap each flywheel's phase and synchronize buffer flips and
 *   the effect playlist.  All protocol logic lives in pov_sync.h
 *   (host-tested); this file is the device shell: it reads the cycle
 *   counter, services two ISRs, packs pixels, and toggles one pin.
 *
 * Effects render the full CANVAS_W × CANVAS_H canvas.  Segmentation is
 * handled entirely in the ISR, which packs only this segment's pixels
 * into a local DMALEDController frame.
 *
 * Include directly from target .ino files:
 *   #include "../../hardware/pov_segmented.h"
 */
#pragma once
#include "led.h"
#include "pov_segment_map.h" // pure index math (host-testable; see that file)
#include "pov_sync.h"        // pure sync protocol (host-testable; see that file)

#ifdef ARDUINO
#include <Arduino.h>

#ifdef USE_DMA_LEDS
  #include "dma_led.h"
#else
  // A segmented build ships one uniform firmware to every board; the master
  // (segment 0) is elected at runtime from the ID pins, so ANY board may become
  // the conductor that drives the sync wire. The master's pulse-width contract
  // (spec §5.2: a "tens of µs" HIGH bracketing the LED work, dropped right after)
  // cannot be honored over the blocking FastLED transport — FastLED.show() stalls
  // the flywheel ISR for the full strip transfer, widening the pulse far past spec
  // and skewing column timing. That combination is untested and unshippable, so
  // fail the build rather than silently degrade the only multi-board transport.
  // (No real target hits this: targets/Phantasm defines USE_DMA_LEDS; the host
  // unit test never includes this Arduino-only driver.)
  #error "POVSegmented requires USE_DMA_LEDS (the Phantasm DMA LED transport): the FastLED fallback cannot honor the master sync pulse-width contract (spec §5.2)."
#endif

#include "canvas.h"
#include "geometry.h"
#include "memory.h"

#include <atomic>

/**
 * @brief Multi-Teensy segmented POV display driver.
 * @tparam S   Total number of physical LEDs across the full strip (both arms).
 * @tparam N   Number of Teensy segments (must be even; N/2 per arm).
 * @tparam RPM Rotations per minute of the spinner.
 *
 * Each segment drives S/N LEDs on a single arm.  The top segment on each
 * arm starts at the N pole (y=0) and counts upward.  The bottom segment
 * starts at the S pole (y=H-1) and counts downward — its strip is
 * physically reversed.
 */
template <int S, int N, int RPM>
class POVSegmented {

  // ── Compile-time geometry ───────────────────────────────────────────

  static constexpr int PPS          = S / N;     /**< Pixels per segment.       */
  static constexpr int ROWS         = S / 2;     /**< Canvas rows (height).     */
  static constexpr int SEGS_PER_ARM = N / 2;     /**< Segments on each arm.     */

  static_assert(S % N == 0,
      "Total pixel count must be evenly divisible by segment count");
  static_assert(N % 2 == 0,
      "Segment count must be even (equal split across two arms)");
  static_assert(S >= N,
      "Must have at least one pixel per segment");
  static_assert((N & (N - 1)) == 0 && N <= 4,
      "Segment ID is decoded from 2 GPIO pins as (raw & (N-1)); that requires "
      "N to be a power of two and N <= 4 (2 ID bits)");

  // ── Pin assignments ─────────────────────────────────────────────────

  /** @brief GPIO pins for hardware ID (active-low with internal pull-up). */
  static constexpr int PIN_ID0 = 21;
  static constexpr int PIN_ID1 = 22;

  /** @brief Sync-symbol output (segment 0 / master only). */
  static constexpr int PIN_FRAME_SYNC_OUT = 3;

  /** @brief Sync-symbol input (downstream boards decode symbol bursts). */
  static constexpr int PIN_FRAME_SYNC_IN = 4;

  // ── Flywheel timing ─────────────────────────────────────────────────

  /** @brief Nominal column period in µs (T0 ≈ 434 µs at 480 RPM × 288). */
  static constexpr float kColumnUs =
      1000000.0f * 60.0f / (float(RPM) * float(CANVAS_W));

  /**
   * @brief Flywheel wake-up oversampling factor.
   *
   * Wake-ups are advisory (position comes from the cycle counter), so the
   * only cost of a coarse wake grid is quantization: a board renders column
   * k at its first wake after k's instant, and the master emits a symbol's
   * first pulse at its first wake after the boundary instant. The timer grid
   * is not phase-locked to the flywheel, so at 1 wake per column both lags
   * reach a full column. Waking 8× per column bounds both to ⅛ column — far
   * inside the §5.2 self-censor budget and below any visible seam — for
   * ~18 kHz of near-empty ISR entries (a position computation and compare).
   */
  static constexpr int kOversample = 8;

  /**
   * @brief HD107S SPI clock for the Phantasm DMA path, in Hz.
   *
   * 24 MHz (vs the TeensySPIDMA 12 MHz default) roughly halves the per-column
   * transfer time, restoring the ~2× headroom the README §1 Phantasm table
   * documents over the ~434 µs column period. HD107S (APA102-compatible) parts
   * clock well past this; the column ISR, not the strip, is the binding budget.
   */
  static constexpr uint32_t SPI_CLOCK_HZ = 24000000;

public:

  /** @brief Foreground effect constructor: builds, arena-configures, and
   *  init()s one roster entry, ready to draw its first frame. */
  using EffectFactory = Effect *(*)();

  /**
   * @brief Initializes hardware: reads segment ID and configures the LED
   *        driver.
   * @details CONTRACT — construct only from setup(), never as a file-scope
   *          global. This constructor performs hardware I/O directly: read_id()
   *          does pinMode/delay/Serial/DWT-enable, and the body then brings up
   *          the DMA/FastLED driver — all valid only once the Arduino core is
   *          initialized. dma_led.h deliberately keeps hardware bring-up out of
   *          its constructor (an explicit begin()) and warns against
   *          constructor-time I/O; this class diverges on purpose because its
   *          sole instantiation site is the Phantasm setup() (a function-local
   *          singleton), so the "after core init" precondition always holds. Do
   *          not promote this object to a global or construct it before setup().
   */
  POVSegmented() {
    read_id();
    configure_segment();

#ifdef USE_DMA_LEDS
    ledController_.begin();
    ledController_.setCorrection(255, 176, 240);   // TypicalLEDStrip
    ledController_.setTemperature(255, 147, 41);    // Candle
    ledController_.setBrightness(255);
#else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(
        leds_, PPS);
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
#endif

    // SysTick (millis) runs at NVIC priority 32 — the Teensy 4 boot default —
    // and so may briefly preempt the flywheel ISR, whose IntervalTimer sits
    // at the default 128. That is fine: the handler is short, shares no state
    // with the ISR, and the flywheel tolerates wake-up jitter by construction
    // (position is time-derived).

    // The flywheel timebase reads DWT->CYCCNT. Teensyduino normally enables it
    // during startup, but enable it explicitly here rather than silently
    // depending on that: TRCENA gates the DWT block on, then CYCCNTENA starts
    // the cycle counter. (Enabling unconditionally, not asserting — the counter
    // is running either way once these two writes land.)
    ARM_DEMCR |= ARM_DEMCR_TRCENA;
    ARM_DWT_CTRL |= ARM_DWT_CTRL_CYCCNTENA;

    Serial.print("[Phantasm] Segment ");
    Serial.print(segment_id_);
    Serial.print(arm_b_ ? " arm-B" : " arm-A");
    Serial.print(y_step_ < 0 ? " (bottom, reversed)" : " (top)");
    Serial.print(segment_id_ == 0 ? " MASTER" : "");
    Serial.print(" | y_base=");
    Serial.print(y_base_);
    Serial.print(" y_step=");
    Serial.print(y_step_);
    Serial.print(" pixels=");
    Serial.println(PPS);
  }

  /**
   * @brief Total DMA frames dropped on overrun since boot.
   *
   * The flywheel ISR drops a frame (rather than blocking) when a prior DMA
   * transfer is still in flight, so a nonzero — and
   * especially a rising — value means the per-column ISR/transfer budget is
   * being exceeded. Always 0 on the FastLED build, which has no DMA
   * controller.
   * @return Cumulative count of DMA frames dropped on overrun since boot; 0 on
   *         the FastLED build.
   */
  uint32_t overrun_count() const {
#if defined(USE_DMA_LEDS)
    return ledController_.getOverrunCount();
#else
    return 0;
#endif
  }

  /**
   * @brief Runs the synchronized show forever (spec §6).
   *
   * The playlist is epoch-counted, not millis()-gated: the master counts
   * revolutions on its own timebase and broadcasts an EPOCH symbol when an
   * effect's revolutions elapse; every board (master included) then tears
   * down its current effect, constructs the next roster entry during the
   * K-revolution commit window (display black), and swaps to its frame 0 at
   * exactly the same boundary. Downstream boards join the running show via
   * the index beacon — they never assume the playlist position.
   *
   * @param factories     One constructor per roster entry (HS_EFFECT_LIST
   *                      order — identical on every board).
   * @param effect_count  Roster length.
   */
  [[noreturn]] void run_show(const EffectFactory *factories,
                             int effect_count) {
    factories_ = factories;

    const pov::sync::Config cfg = pov::sync::phantasm_config(
        F_CPU, RPM, CANVAS_W, effect_count);
    HS_CHECK(cfg.valid(), "pov::sync::Config invariants");
    sync_ = pov::sync::SyncBoard(cfg);

    const bool master = (segment_id_ == 0);
    if (master) {
      pinMode(PIN_FRAME_SYNC_OUT, OUTPUT);
      digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
    }
    pinMode(PIN_FRAME_SYNC_IN, INPUT);

    sync_.seed(ARM_DWT_CYCCNT, master);

    // Attach ONCE per show; the flywheel is the timebase and persists across
    // epochs (spec §8.5) — there is no per-effect attach/detach.
    if (!master) {
      attachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC_IN),
                      sync_edge_isr, RISING);
    }
    // begin() returns false if all four PIT channels are already taken; an
    // unstarted flywheel is a silent dead board (no timebase, no display), so
    // trap at the violation site rather than spin forever in the loop below.
    HS_CHECK(timer_.begin(flywheel_isr, kColumnUs / float(kOversample)),
             "flywheel IntervalTimer failed to start (no PIT channel)");

    // ── Foreground: construct effects on request, render, report ────────
    Effect *cur = nullptr;
    uint32_t built_gen = 0;
#if defined(USE_DMA_LEDS)
    uint32_t last_overrun = ledController_.getOverrunCount();
#endif
    pov::sync::Telemetry last_tm{};
    unsigned long last_report = millis();

    for (;;) {
      const uint32_t bw = sync_.build_word();
      const uint32_t gen = pov::sync::SyncBoard::build_gen_of(bw);
      if (gen != built_gen) {
        // The ISR owns the live-effect pointer; ask it to let go before the
        // old instance is destroyed (it acknowledges within one wake-up).
        release_req_.fetch_add(1, std::memory_order_relaxed);
        const unsigned long t0 = micros();
        while (release_ack_.load(std::memory_order_relaxed) !=
               release_req_.load(std::memory_order_relaxed)) {
          // 100 ms is a deliberately loose liveness bound, not a tight timing
          // budget: the flywheel ISR acks within one wake interval
          // (kColumnUs / kOversample ≈ 54 µs at 480 RPM × 288), so 100 ms is
          // ~1800 wake intervals of slack. It absorbs any plausible transient
          // (a long `__disable_irq` window, a stalled wake) and traps only a
          // genuinely wedged ISR — never a momentary handshake delay.
          HS_CHECK(micros() - t0 < 100000UL,
                   "flywheel ISR failed to release the live effect");
        }
        delete cur;
        // Deterministic content (spec §2): every board restarts the shared
        // RNG stream per effect — exactly what the simulator does for a
        // standalone effect — so frames are identical across boards
        // regardless of boot or join history.
        hs::random().seed(1337);
        cur = factories_[pov::sync::SyncBoard::build_index_of(bw)]();
        // ISR seam (project doctrine: trap at the cold construction site, not
        // the hot ISR): render_column() walks PPS pixels over canvas rows in
        // [0, ROWS) and indexes buf[y * width + x_col], in-bounds only when the
        // effect's canvas height equals ROWS (S/2) AND its width equals the
        // CANVAS_W the sync engine sweeps x over — x_col derives from x in
        // [0, CANVAS_W) but is bounds-checked against the row stride w =
        // width(). Trap a resolution mismatch on either axis here, before the
        // flywheel ISR can take this instance live.
        HS_CHECK(cur->height() == ROWS,
                 "POVSegmented: effect canvas height must equal S/2 (ROWS)");
        HS_CHECK(cur->width() == CANVAS_W,
                 "POVSegmented: effect canvas width must equal CANVAS_W");
        cur->draw_frame(); // frame 0, queued; fresh buffers never block
        hs::disable_interrupts();
        // Release store last so it publishes both the new generation and every
        // constructor/draw_frame() write into *cur; the bracket keeps the
        // (effect, gen) pair atomic against the ISR.
        pending_gen_.store(gen, std::memory_order_relaxed);
        pending_effect_.store(cur, std::memory_order_release);
        hs::enable_interrupts();
        built_gen = gen;
      }

      // Render only once the ISR has taken this instance live; before that
      // the queued frame 0 is waiting for the commit/join boundary. The
      // Canvas buffer_free() gate paces this to exactly one step() per flip.
      if (cur && consumed_gen_.load(std::memory_order_relaxed) == built_gen) {
        const unsigned long f0 = micros();
        cur->draw_frame();
        if (hs::debug) {
          Serial.print("ft ");
          Serial.println(micros() - f0);
        }
      }

      // Health telemetry (spec §8.6): foreground-polled, never ISR-printed,
      // emitted only on change — the existing DMA-overrun pattern. The
      // snapshot read races the ISR's counter increments; for debug counters
      // a torn read is benign — the worst case is a one-cycle-late or briefly
      // inconsistent debug print (behind hs::debug), never anything that feeds
      // control flow. This is a deliberate, documented data race. If telemetry
      // is ever promoted beyond debug output, snapshot it under a brief IRQ-off
      // bracket like the mailbox handoff.
      if (hs::debug && millis() - last_report >= 1000UL) {
        last_report = millis();
        const pov::sync::Telemetry tm = sync_.telemetry();
        if (memcmp(&tm, &last_tm, sizeof tm) != 0) {
          // hs::log (integer-only vsniprintf + Serial.println) rather than
          // Serial.printf: Teensy's Print::printf routes through vdprintf, which
          // unconditionally drags newlib's float formatter (_dtoa_r + bignum
          // helpers, ~5 KB of ITCM) into the image. These counters are all %lu.
          hs::log("sync acc=%lu rej=%lu inv=%lu cens=%lu abrt=%lu "
                  "bok=%lu brej=%lu fix=%lu rmis=%lu lock=%lu "
                  "flip=%lu coast=%lu epi=%lu",
                  (unsigned long)tm.symbols_accepted,
                  (unsigned long)tm.symbols_rejected_gate,
                  (unsigned long)tm.symbols_discarded_invalid,
                  (unsigned long)tm.emit_censored,
                  (unsigned long)tm.emit_aborted,
                  (unsigned long)tm.beacons_ok,
                  (unsigned long)tm.beacons_rejected,
                  (unsigned long)tm.beacon_index_corrections,
                  (unsigned long)tm.beacon_rev_mismatches,
                  (unsigned long)tm.lock_transitions,
                  (unsigned long)tm.flips,
                  (unsigned long)tm.max_coast_halves,
                  (unsigned long)tm.epochs_refractory_ignored);
          last_tm = tm;
        }
#if defined(USE_DMA_LEDS)
        const uint32_t overruns = ledController_.getOverrunCount();
        if (overruns != last_overrun) {
          Serial.print("overrun ");
          Serial.println(overruns);
          last_overrun = overruns;
        }
#endif
      }
    }
  }

private:

  // ── Hardware ID ─────────────────────────────────────────────────────

  /**
   * @brief Reads the 2-bit hardware segment ID from GPIO pins.
   *
   * Pins are configured as INPUT_PULLUP; grounding a pin sets its bit.
   * The result is inverted so that all-floating = ID 0 (master).
   *
   * A floating or cold-soldered strap reads HIGH, which inverts toward ID 0 —
   * silently promoting the board to a second master and driving the shared
   * sync wire into contention. To catch unstable straps we sample twice and
   * trap on disagreement: a stable strap is an invariant.
   *
   * Note: this validates *this* board's strap only. It cannot detect a *peer*
   * holding the same ID — the sync wire is push-pull, so a duplicate master
   * is not observable from a single board without a dedicated (open-drain)
   * arbitration line. That cross-check is intentionally out of scope; the
   * debounce above removes the most likely cause (a floating pin).
   */
  void read_id() {
    pinMode(PIN_ID0, INPUT_PULLUP);
    pinMode(PIN_ID1, INPUT_PULLUP);
    // Settle the pull-ups, then debounce by sampling three times spread across a
    // ~10 ms window and requiring all three to agree. The heuristic 10 ms settle
    // assumes the strap RC (internal pull-up x trace+pin capacitance) resolves
    // well within it; the resampling guards against reading a strap mid-
    // transition.
    //
    // Why three samples over ~10 ms, not one re-sample at 2 ms: a single 2 ms
    // re-sample only catches bounces shorter than that gap. A strap that settles
    // slower than 2 ms (long traces / added caps) could read identically at both
    // samples while still mid-transition and pass — mis-ID'ing the board into a
    // *second* master and driving the push-pull sync wire into bus contention, a
    // fault no single board can detect (see the class note on the absent peer
    // cross-check). Three samples spaced ~5 ms widen the debounce window an order
    // of magnitude past the expected sub-ms RC for a one-time ~20 ms boot cost.
    // Widen the spacing further if a future board adds long strap traces or caps.
    delay(10);  // settle time for pull-ups

    const int raw0 = digitalReadFast(PIN_ID0) | (digitalReadFast(PIN_ID1) << 1);
    for (int i = 0; i < 2; ++i) {
      delay(5);  // spread the resamples across the ~10 ms debounce window
      const int rawN = digitalReadFast(PIN_ID0) | (digitalReadFast(PIN_ID1) << 1);
      HS_CHECK(rawN == raw0);  // unstable strap → mis-ID → bus contention
    }

    // ~ binds tighter than &: invert the 2-bit reading (so all-floating
    // pull-ups => ID 0), then mask off the inverted high bits with & (N-1).
    // The mask is load-bearing — without it the full-width ~ would set every
    // upper bit. Parens make the (~...) & (N-1) grouping explicit.
    segment_id_ = (~raw0) & (N - 1);
  }

  // ── Segment mapping ─────────────────────────────────────────────────

  /**
   * @brief Computes the precomputed ISR mapping from hardware segment ID.
   *
   * Layout (N=4, SEGS_PER_ARM=2, ROWS=144):
   *
   *   Segments 0 .. SEGS_PER_ARM-1  →  arm A (x column)
   *   Segments SEGS_PER_ARM .. N-1   →  arm B (x + W/2 column)
   *
   *   Within each arm (arm_seg 0 = top, arm_seg 1 = bottom):
   *
   *     arm_seg 0 (top):    LED 0 at N pole (y=0), strip runs toward junction.
   *                         y_base = 0, y_step = +1  →  y ∈ [0, PPS-1]
   *
   *     arm_seg 1 (bottom): LED 0 at S pole (y=ROWS-1), strip runs toward junction.
   *                         y_base = ROWS-1, y_step = -1, descending over the PPS
   *                         rows  →  y ∈ [ROWS-PPS, ROWS-1] (strip physically
   *                         reversed). The lower bound is ROWS-PPS for any S/N;
   *                         it only coincides with PPS in the N=4 case (ROWS=2*PPS).
   */
  void configure_segment() {
    const pov::SegmentMap m = pov::segment_map(segment_id_, S, N);
    arm_b_  = m.arm_b;
    y_base_ = m.y_base;
    y_step_ = m.y_step;
  }

  // ── ISRs ────────────────────────────────────────────────────────────

  /**
   * @brief Sync-wire edge ISR (downstream boards only): a pure publisher.
   *
   * Applies the glitch filter and records the edge into the mailbox; touches
   * no flywheel, flip, or epoch state (spec §8.2 single-writer model). NVIC
   * priority relative to the flywheel ISR is therefore free — preemption has
   * no correctness consequence.
   */
  static FASTRUN void sync_edge_isr() {
    sync_.on_sync_edge(ARM_DWT_CYCCNT);
  }

  /**
   * @brief Flywheel ISR: the sole owner of all sync state (spec §8).
   *
   * Paced by an IntervalTimer at T0/kOversample as a wake-up only — the
   * cycle counter decides which column it is (spec §4.1), so a late, early,
   * or coalesced wake-up cannot inject drift: the ISR is idempotent when the
   * column is unchanged and skip-tolerant when it jumped (the skipped
   * columns were masked precisely because the strip was busy).
   */
  static FASTRUN void flywheel_isr() {
    const uint32_t now = ARM_DWT_CYCCNT;

    // Complete a deferred dark-path sync pulse from the previous wake. On the
    // dark-latched path the ISR body is only a few dozen instructions (no
    // render_black / render_column), so a same-wake HIGH→LOW pulse collapses to
    // ~50–200 ns — two to three orders below spec §5.2's "tens of µs," exactly
    // during the boot-acquisition window downstream boards lock onto. The dark
    // path instead holds the pin HIGH and drops it here at the next wake, giving
    // the pulse a full wake interval (~T0/kOversample) of width.
    if (sync_low_pending_) {
      digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
      sync_low_pending_ = false;
    }

    // Mailbox handoff (spec §8.2): a brief IRQ-off copy in the consumer.
    pov::sync::BurstSnapshot burst;
    const pov::sync::BurstSnapshot *bp = nullptr;
    if (segment_id_ != 0) {
      __disable_irq();
      if (sync_.mailbox().try_claim(now, sync_.gap_timeout_cycles(), &burst))
        bp = &burst;
      // Retire a stale glitch-filter reference so the cycle counter cannot wrap
      // out from under it during a long wire silence (spec §8).
      sync_.mailbox().age_prior(now, sync_.glitch_filter_cycles());
      __enable_irq();
    }

    const pov::sync::TickActions a = sync_.tick(now, bp);

    // Pin write first, LED work after (spec §5.2): emission timing carries
    // only ISR-entry jitter plus the tick() computation.
    if (a.pulse)
      digitalWriteFast(PIN_FRAME_SYNC_OUT, HIGH);

    // Release handshake: the foreground wants the live pointer dropped so it
    // can destroy the instance (epoch teardown / beacon rebuild).
    if (release_ack_.load(std::memory_order_relaxed) !=
        release_req_.load(std::memory_order_relaxed)) {
      live_effect_ = nullptr;
      release_ack_.store(release_req_.load(std::memory_order_relaxed),
                         std::memory_order_relaxed);
    }

    // Swap in the foreground-constructed pending effect, only ever at a ZERO
    // boundary (frame parity). Two paths:
    //   commit — the B+K epoch deadline (spec §6.1). The next effect MUST be
    //            ready; an init that outruns K revolutions is an invariant
    //            violation and traps rather than silently skewing the show.
    //   join   — first display (boot / beacon join / index correction): take
    //            the pending effect at the next join-grid boundary once it
    //            exists. The grid (rev ≡ 0 mod join_grid_revs) is what makes
    //            all four boards go live at the SAME crossing at boot, frame
    //            counters aligned.
    if (a.commit) {
      // Acquire load pairs with the foreground's release store, ordering the
      // pending instance's construction writes before any dereference below.
      Effect *p = pending_effect_.load(std::memory_order_acquire);
      // This equality check is sound only because pending_gen_ is stable across
      // the construction window: it is set once when the window opens (the
      // foreground's publish_build) and no beacon/symbol path bumps it again
      // while a commit is pending — SyncBoard::handle_beacon_burst suppresses
      // beacon index-corrections during commit_pending, and the master gates its
      // epoch publish on !commit_pending. So a mismatch here means only the
      // genuine fault in the message — init outran the K-revolution window — and
      // never a mid-window generation bump from a beacon correction.
      HS_CHECK(p && pending_gen_.load(std::memory_order_relaxed) ==
                        pov::sync::SyncBoard::build_gen_of(sync_.build_word()),
               "epoch commit: effect init exceeded the K-revolution window");
      live_effect_ = p;
      consumed_gen_.store(pending_gen_.load(std::memory_order_relaxed),
                          std::memory_order_relaxed);
    } else if (a.join_boundary && !a.dark && live_effect_ == nullptr) {
      // Re-read pending_gen_ against the current build_word() so a late joiner
      // only adopts an effect still matching the sync wire's advertised
      // generation. A momentary visibility lag between the foreground's
      // pending_gen_ publish and build_word() reflecting it can fail this match
      // and skip the join — benign: join_boundary recurs on the next grid
      // boundary, so the board simply joins one grid step later.
      Effect *p = pending_effect_.load(std::memory_order_acquire);
      const uint32_t pg = pending_gen_.load(std::memory_order_relaxed);
      if (p && pg != consumed_gen_.load(std::memory_order_relaxed) &&
          pg == pov::sync::SyncBoard::build_gen_of(sync_.build_word())) {
        live_effect_ = p;
        consumed_gen_.store(pg, std::memory_order_relaxed);
      }
    }

    // Flip whenever the effect is live — even during the dark commit window,
    // where the foreground may be blocked in the Canvas buffer_free() gate on
    // its final frame of the outgoing effect; advance_display() is what
    // releases it to go tear the effect down.
    Effect *e = live_effect_;
    if (a.flip && e)
      e->advance_display();

    // Tracks whether this wake did real LED work between the pulse edges. The
    // bright path (render_column) and the first dark frame (render_black) take
    // tens of µs and so width the pulse themselves; the dark-latched path skips
    // both and must defer the pin drop (below) to keep the pulse spec-wide.
    bool did_render = false;
    if (a.dark || e == nullptr) {
      // Fail-dark (spec §5.3/§6.3): show nothing rather than the wrong
      // thing. One black frame latches the strip; then stay idle.
      //
      // Latch only once the black frame is actually accepted: submitFrame()
      // drops on a DMA overrun (the prior bright column's transfer is still in
      // flight), and latching on a dropped frame would hold that stale bright
      // column for the entire dark window — inverting the design's
      // missed-never-wrong asymmetry. Leaving dark_latched_ false retries on
      // the next wake until the black frame lands, so a drop self-heals in one
      // column tick (~434 µs) instead of latching indefinitely.
      if (!dark_latched_ && render_black()) {
        dark_latched_ = true;
        did_render = true;
      }
    } else {
      dark_latched_ = false;
      if (a.render_column >= 0) {
        render_column(e, a.render_column);
        did_render = true;
      }
    }

    if (a.pulse) {
      if (did_render) {
        // Bright / first-dark path: the ISR body already gave the pulse its
        // spec width (tens of µs). Drop it now.
        digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
      } else {
        // Dark-latched (or render-less) path: the body is too short to width
        // the pulse. Hold the pin HIGH and drop it at the next wake — see the
        // deferred-drop at the top of the ISR.
        sync_low_pending_ = true;
      }
    }
  }

  // ── Pixel packing ───────────────────────────────────────────────────

  /**
   * @brief Packs this segment's pixels for canvas column @p x and submits.
   * @param e Live effect supplying the display buffer to sample.
   * @param x Canvas column index in [0, CANVAS_W); arm-B segments sample
   *          column x + W/2.
   * @details The loop is branchless — all per-segment decisions are resolved at
   *          boot in configure_segment().  Arm B segments read from x + W/2
   *          (opposite half of the image).
   */
  static FASTRUN void render_column(Effect *e, int x) {
    const int w = e->width();
    const int x_col = pov::segment_x_col(arm_b_, x, w);

    // ISR fast path: fetch the display buffer base once and index it
    // directly, dropping PPS per-pixel virtual get_pixel() dispatches per
    // column. Sound here because (1) prev_ is stable for this whole column —
    // advance_display() runs at the boundary, before rendering — and (2) no
    // effect reachable on this segmented path overrides get_pixel (the one
    // effect that does is out of scope here and never runs on this path).
    const Pixel *buf = e->display_buffer();

#if defined(USE_DMA_LEDS)
    auto &frame = ledController_.backFrame();
    int y = y_base_;
    for (int i = 0; i < PPS; ++i, y += y_step_) {
      frame.packPixel(i, buf[y * w + x_col]);
    }
    // Steady-state column path intentionally drops the accept/overrun result:
    // a dropped image column self-heals next tick (see submitFrame's doc and
    // render_black, which DOES gate on the return for the fail-dark latch).
    (void)ledController_.submitFrame(e->show_bg());
#else
    int y = y_base_;
    for (int i = 0; i < PPS; ++i, y += y_step_) {
      leds_[i] = static_cast<CRGB>(buf[y * w + x_col]);
    }
    FastLED.show();
    if (e->show_bg()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }
#endif
  }

  /**
   * @brief Submits one all-black frame (ACQUIRE / construction window).
   * @return true if the black frame was accepted by the LED transport; false
   *         if it was dropped on a DMA overrun (caller must retry, not latch).
   */
  static FASTRUN bool render_black() {
#if defined(USE_DMA_LEDS)
    auto &frame = ledController_.backFrame();
    for (int i = 0; i < PPS; ++i) {
      frame.packPixel(i, Pixel(0, 0, 0));
    }
    return ledController_.submitFrame(false);
#else
    FastLED.showColor(CRGB(0, 0, 0));
    return true; // synchronous path always lands
#endif
  }

  // ── Static state ────────────────────────────────────────────────────

#ifndef USE_DMA_LEDS
  static CRGB leds_[PPS]; /**< FastLED output buffer, one CRGB per segment pixel. */
#endif

  /**
   * @brief The sync engine: sole owner of all sync/flywheel state.
   * @details Written ONLY by the flywheel ISR (tick()) and the edge ISR
   *          (mailbox publisher) per the spec §8 single-writer model. Foreground
   *          reads are single aligned words (build_word) or debug telemetry.
   */
  static pov::sync::SyncBoard sync_;
  static IntervalTimer timer_;             /**< Flywheel wake-up timer (PIT channel).   */
  static const EffectFactory *factories_;  /**< Roster of effect constructors (HS_EFFECT_LIST order). */

  /**
   * @brief Effect handoff state between the foreground and the ISR.
   * @details Ownership: the foreground constructs and deletes; the ISR only ever
   *          dereferences the instance it has been handed. live_effect_ is
   *          ISR-written only; pending_* are foreground-written only (published
   *          under a brief interrupts-off bracket); the release_req_/release_ack_
   *          pair is the teardown handshake (foreground bumps req, ISR drops the
   *          live pointer and copies req to ack).
   *
   *          pending_effect_ is published with a release store and consumed with
   *          an acquire load: the release/acquire edge — not the IRQ bracket's
   *          compiler barrier alone — is what orders the freshly-constructed
   *          instance's member writes (vtable, buffers) before the ISR
   *          dereferences it. The interrupts-off bracket still makes the
   *          (effect, gen) publish atomic with respect to the ISR. The acquire
   *          load runs only at a commit/join boundary, never on the per-column
   *          hot path (which reads the plain, ISR-owned live_effect_), so the
   *          ordering carries no hot-path cost. Mirrors Canvas's std::atomic
   *          ISR handoff rather than the bare volatile it replaces.
   */
  static Effect *live_effect_;             /**< Effect the ISR renders; ISR-owned.       */
  static std::atomic<Effect *> pending_effect_; /**< Next effect awaiting commit; foreground-written (release), ISR-read (acquire). */
  // Cross-context counters as std::atomic (relaxed), not bare volatile: each
  // has a single writer (named below) but is read from the other context, so a
  // plain volatile RMW/read is a formal data race the way dma_led.h's overrun
  // counters were before they moved to std::atomic (:388-396). Relaxed ordering
  // adds no fence on the single-core Cortex-M — these carry no happens-before of
  // their own (the (effect, gen) handoff is ordered by pending_effect_'s
  // release/acquire), so they compile to the same load/store the volatile did,
  // now well-defined. Completes the volatile->atomic migration pending_effect_
  // already made.
  static std::atomic<uint32_t> pending_gen_;   /**< Build generation of pending_effect_; foreground-written. */
  static std::atomic<uint32_t> consumed_gen_;  /**< Build generation taken live; ISR-written.  */
  static std::atomic<uint32_t> release_req_;   /**< Teardown request counter; foreground-written. */
  static std::atomic<uint32_t> release_ack_;   /**< Teardown acknowledge counter; ISR-written.  */
  static bool dark_latched_;               /**< True once the black frame has latched; ISR-owned. */
  static bool sync_low_pending_;           /**< ISR-owned: dark-path pulse drop deferred to next wake. */

  static int segment_id_;                  /**< Decoded 2-bit hardware segment ID.       */
  static bool arm_b_;                      /**< True if this segment lives on arm B (x + W/2). */
  static int y_base_;                      /**< Canvas row of this segment's LED 0.      */
  static int y_step_;                      /**< Row stride per LED: +1 (top) or -1 (bottom, reversed). */
#if defined(USE_DMA_LEDS)
  static DMALEDController<PPS> ledController_; /**< DMA SPI LED controller for the segment strip. */
#endif
};

// ── Static member definitions ───────────────────────────────────────────

template <int S, int N, int RPM>
pov::sync::SyncBoard POVSegmented<S, N, RPM>::sync_{pov::sync::Config{}};

template <int S, int N, int RPM>
IntervalTimer POVSegmented<S, N, RPM>::timer_;

template <int S, int N, int RPM>
const typename POVSegmented<S, N, RPM>::EffectFactory
    *POVSegmented<S, N, RPM>::factories_ = nullptr;

template <int S, int N, int RPM>
Effect *POVSegmented<S, N, RPM>::live_effect_ = nullptr;

template <int S, int N, int RPM>
std::atomic<Effect *> POVSegmented<S, N, RPM>::pending_effect_{nullptr};

template <int S, int N, int RPM>
std::atomic<uint32_t> POVSegmented<S, N, RPM>::pending_gen_{0};

template <int S, int N, int RPM>
std::atomic<uint32_t> POVSegmented<S, N, RPM>::consumed_gen_{0};

template <int S, int N, int RPM>
std::atomic<uint32_t> POVSegmented<S, N, RPM>::release_req_{0};

template <int S, int N, int RPM>
std::atomic<uint32_t> POVSegmented<S, N, RPM>::release_ack_{0};

template <int S, int N, int RPM>
bool POVSegmented<S, N, RPM>::dark_latched_ = false;

template <int S, int N, int RPM>
bool POVSegmented<S, N, RPM>::sync_low_pending_ = false;

template <int S, int N, int RPM>
int POVSegmented<S, N, RPM>::segment_id_ = 0;

template <int S, int N, int RPM>
bool POVSegmented<S, N, RPM>::arm_b_ = false;

template <int S, int N, int RPM>
int POVSegmented<S, N, RPM>::y_base_ = 0;

template <int S, int N, int RPM>
int POVSegmented<S, N, RPM>::y_step_ = 1;

#ifndef USE_DMA_LEDS
template <int S, int N, int RPM>
CRGB POVSegmented<S, N, RPM>::leds_[POVSegmented<S, N, RPM>::PPS];
#endif

#if defined(USE_DMA_LEDS)
// DMAMEM (OCRAM): the controller's HD107SFrame buffers are the actual eDMA TX
// source, so they must live in DMA-reachable, cached OCRAM — which is exactly
// what HD107SFrame's arm_dcache_flush() assumes. Default placement is DTCM,
// where that flush is a dead no-op; here it does the required write-back.
template <int S, int N, int RPM>
DMAMEM DMALEDController<POVSegmented<S, N, RPM>::PPS>
    POVSegmented<S, N, RPM>::ledController_{POVSegmented<S, N, RPM>::SPI_CLOCK_HZ};
#endif

#endif // ARDUINO
