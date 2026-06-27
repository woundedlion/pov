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
   * reach a full column. Waking 8× per column drops the lag to ⅛ column on the
   * cheap position-only wakes. The bound is looser on the wakes that actually
   * render or emit: a column render (~96 µs, README:1265) overruns the ~54 µs
   * wake grid and the PIT ISR cannot re-enter, so the next 1–2 grid slots
   * coalesce onto it and that column's lag is closer to ~2/8 column. Still far
   * inside the §5.2 self-censor budget and below any visible seam, and the
   * design is unaffected — position is time-derived and tick() is idempotent /
   * skip-tolerant. The remaining ~18 kHz of wakes are near-empty (a position
   * computation and compare).
   */
  static constexpr int kOversample = 8;

  static_assert(kColumnUs / float(kOversample) >= 1.0f,
      "Flywheel wake period must be >= 1 us for IntervalTimer::begin");

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
   *          the DMA LED driver — all valid only once the Arduino core is
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

    ledController_.begin();
    ledController_.setCorrection(255, 176, 240);   // TypicalLEDStrip
    ledController_.setTemperature(255, 147, 41);    // Candle
    ledController_.setBrightness(255);

    // Enable the DWT cycle counter the flywheel timebase reads: TRCENA gates the
    // DWT block on, then CYCCNTENA starts the counter.
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
   * being exceeded.
   * @return Cumulative count of DMA frames dropped on overrun since boot.
   */
  uint32_t overrun_count() const { return ledController_.getOverrunCount(); }

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

    if (!master) {
      attachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC_IN),
                      sync_edge_isr, RISING);
    }
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
          HS_CHECK(micros() - t0 < 100000UL,
                   "flywheel ISR failed to release the live effect");
        }
        pending_effect_.store(nullptr, std::memory_order_release);
        delete cur;
        // Restart the shared RNG stream per effect so frames match the
        // simulator and across boards regardless of boot/join history (spec §2).
        hs::random().seed(1337);
        cur = factories_[pov::sync::SyncBoard::build_index_of(bw)]();
        HS_CHECK(cur->height() == ROWS,
                 "POVSegmented: effect canvas height must equal S/2 (ROWS)");
        HS_CHECK(cur->width() == CANVAS_W,
                 "POVSegmented: effect canvas width must equal CANVAS_W");
        cur->draw_frame();
        hs::disable_interrupts();
        // Release store last so the (effect, gen) pair publishes atomically and
        // orders every constructor/draw_frame() write before the ISR sees it.
        pending_gen_.store(gen, std::memory_order_relaxed);
        pending_effect_.store(cur, std::memory_order_release);
        hs::enable_interrupts();
        built_gen = gen;
      }

      if (cur && consumed_gen_.load(std::memory_order_relaxed) == built_gen) {
        const unsigned long f0 = micros();
        cur->draw_frame();
        if (hs::debug) {
          Serial.print("ft ");
          Serial.println(micros() - f0);
        }
      }

      // Health telemetry (spec §8.6): foreground-polled, emitted only on change.
      // telemetry() aliases the live, ISR-mutated block, so copy it under a brief
      // IRQ-off bracket or a torn read could mis-fire the comparison below.
      if (hs::debug && millis() - last_report >= 1000UL) {
        last_report = millis();
        __disable_irq();
        const pov::sync::Telemetry tm = sync_.telemetry();
        __enable_irq();
        if (memcmp(&tm, &last_tm, sizeof tm) != 0) {
          // hs::log, not Serial.printf: Teensy's printf drags in newlib's float
          // formatter (~5 KB ITCM); these counters are all %lu.
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
   * Pins are configured as INPUT_PULLUP, so a floating pin reads HIGH and a
   * grounded pin reads LOW. The raw 2-bit reading is inverted, so a grounded pin
   * contributes a 1 and all-floating = ID 0 (master).
   *
   * A floating or cold-soldered strap reads HIGH, which inverts toward ID 0 —
   * silently promoting the board to a second master and driving the shared
   * sync wire into contention. To catch unstable straps we sample three
   * times and trap on disagreement: a stable strap is an invariant.
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
    delay(10);  // settle time for pull-ups

    // Debounce: three samples ~5 ms apart must agree; an unstable strap reads as
    // a second master and drives the push-pull sync wire into bus contention.
    const int raw0 = digitalReadFast(PIN_ID0) | (digitalReadFast(PIN_ID1) << 1);
    for (int i = 0; i < 2; ++i) {
      delay(5);
      const int rawN = digitalReadFast(PIN_ID0) | (digitalReadFast(PIN_ID1) << 1);
      HS_CHECK(rawN == raw0, "unstable segment-ID strap (field/manufacturing fault)");
    }

    // Invert the 2-bit reading (all-floating pull-ups => ID 0), then mask the
    // inverted high bits; the mask is load-bearing.
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

    // Complete a deferred dark-path sync pulse from the previous wake: the
    // dark-latched path's body is too short to width a same-wake pulse to spec
    // §5.2's "tens of µs," so it holds the pin HIGH and drops it here instead.
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

    // Pin write first, LED work after (spec §5.2).
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
    // boundary. Two paths: commit (the B+K epoch deadline, spec §6.1) and join
    // (first display: boot / beacon join / index correction, taken at the next
    // join-grid boundary so all boards go live at the same crossing).
    if (a.commit) {
      Effect *p = pending_effect_.load(std::memory_order_acquire);
      HS_CHECK(p && pending_gen_.load(std::memory_order_relaxed) ==
                        pov::sync::SyncBoard::build_gen_of(sync_.build_word()),
               "epoch commit: effect init exceeded the K-revolution window");
      // render_column samples the raw display buffer, bypassing get_pixel(), so
      // a get_pixel-overriding effect must never go live here.
      HS_CHECK(!p->overrides_get_pixel(),
               "segmented render_column bypasses get_pixel(); a get_pixel-"
               "overriding effect cannot run on this path");
      live_effect_ = p;
      consumed_gen_.store(pending_gen_.load(std::memory_order_relaxed),
                          std::memory_order_relaxed);
    } else if (a.join_boundary && !a.dark && live_effect_ == nullptr) {
      // Adopt only an effect still matching the wire's advertised generation; a
      // visibility lag that fails the match simply joins one grid step later.
      Effect *p = pending_effect_.load(std::memory_order_acquire);
      const uint32_t pg = pending_gen_.load(std::memory_order_relaxed);
      if (p && pg != consumed_gen_.load(std::memory_order_relaxed) &&
          pg == pov::sync::SyncBoard::build_gen_of(sync_.build_word())) {
        HS_CHECK(!p->overrides_get_pixel(),
                 "segmented render_column bypasses get_pixel(); a get_pixel-"
                 "overriding effect cannot run on this path");
        live_effect_ = p;
        consumed_gen_.store(pg, std::memory_order_relaxed);
      }
    }

    // Flip whenever the effect is live, even during the dark commit window:
    // advance_display() is what releases a foreground blocked in buffer_free().
    Effect *e = live_effect_;
    if (a.flip && e)
      e->advance_display();

    bool did_render = false;
    if (a.dark || e == nullptr) {
      // Fail-dark (spec §5.3/§6.3). Latch only once the black frame is actually
      // accepted: latching on a DMA-overrun drop would hold the stale bright
      // column for the whole dark window. Leaving it false retries next wake.
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
        digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
      } else {
        // Body too short to width the pulse: hold HIGH, drop at the next wake.
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

    // ISR fast path: index the display buffer directly, dropping PPS per-pixel
    // virtual get_pixel() dispatches. No effect on this path overrides get_pixel.
    const Pixel *buf = e->display_buffer();

    auto &frame = ledController_.backFrame();
    int y = y_base_;
    for (int i = 0; i < PPS; ++i, y += y_step_) {
      frame.packPixel(i, buf[y * w + x_col]);
    }
    // Drop the accept/overrun result: a dropped image column self-heals next
    // tick (render_black, by contrast, gates on it for the fail-dark latch).
    (void)ledController_.submitFrame(e->strobe_columns());
  }

  /**
   * @brief Submits one all-black frame (ACQUIRE / construction window).
   * @return true if the black frame was accepted by the LED transport; false
   *         if it was dropped on a DMA overrun (caller must retry, not latch).
   */
  static FASTRUN bool render_black() {
    auto &frame = ledController_.backFrame();
    for (int i = 0; i < PPS; ++i) {
      frame.packPixel(i, Pixel(0, 0, 0));
    }
    return ledController_.submitFrame(false);
  }

  // ── Static state ────────────────────────────────────────────────────

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

#if defined(USE_DMA_LEDS)
// ledController_ is intentionally NOT defined out-of-line here. Its HD107SFrame
// buffers are the eDMA TX source and belong in cached, DMA-reachable OCRAM (where
// HD107SFrame's arm_dcache_flush() write-back is meaningful — in DTCM it is a
// dead no-op). DMAMEM (a section attribute) is silently dropped by GCC on a
// vague-linkage template static member, so a generic definition would land in
// DTCM regardless. Each instantiating target therefore defines it as an explicit
// specialization (ordinary strong linkage, so DMAMEM sticks) — see Phantasm.ino.
#endif

#endif // ARDUINO
