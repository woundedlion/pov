/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_segmented.h
 * @brief Multi-Teensy segmented POV display driver for Phantasm.
 *
 * Phantasm uses N Teensys (4 by default, up to 8), each controlling a
 * contiguous segment of LEDs on a single arm of the POV spinner. Northern
 * segments run toward the junction in increasing canvas-row order; southern
 * segments run toward it in decreasing order.
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
#include "render/led.h"
#include "pov_segment_map.h" // pure index math (host-testable; see that file)
#include "pov_sync.h"        // pure sync protocol (host-testable; see that file)
#include "pov_handoff.h"     // pure effect-handoff state machine (host-testable)

#ifdef ARDUINO
#include <Arduino.h>

#ifdef USE_DMA_LEDS
  #include "dma_led.h"
#else
  #error "POVSegmented requires USE_DMA_LEDS (the Phantasm DMA LED transport): the FastLED fallback cannot honor the master sync pulse-width contract (spec §5.2)."
#endif

#include "render/canvas.h"
#include "math/geometry.h"
#include "engine/memory.h"

#include <atomic>

#ifdef HS_PROFILE_ENABLE
namespace hs {
/** Column-ISR profiling accumulators (see IsrCycleStats): the whole flywheel
 *  wake, render_column's pixel pack, and the submitFrame DMA marshal+kick.
 *  ISR-written; read + reset from the foreground under IRQ-off. */
inline IsrCycleStats g_flywheel_wake_cycles;
inline IsrCycleStats g_column_pack_cycles;
inline IsrCycleStats g_dma_submit_cycles;
} // namespace hs
#endif

/**
 * @brief Multi-Teensy segmented POV display driver.
 * @tparam S   Total number of physical LEDs across the full strip (both arms).
 * @tparam N   Number of Teensy segments (must be even; N/2 per arm).
 * @tparam RPM Rotations per minute of the spinner.
 *
 * Each segment drives S/N LEDs on a single arm. Northern segments count
 * upward; southern segments count downward from the S pole toward the
 * junction. N=2 assigns one forward strip to each complete arm.
 */
template <int S, int N, int RPM>
class POVSegmented {

  // ── Compile-time geometry ───────────────────────────────────────────

  static constexpr int PPS          = S / N;     /**< Pixels per segment.       */
  static constexpr int ROWS         = S / 2;     /**< Canvas rows (height).     */
  static constexpr int SEGS_PER_ARM = N / 2;     /**< Segments on each arm.     */

  static_assert(RPM > 0, "RPM must be positive (COLUMN_US divides by RPM)");
  static_assert(S % N == 0,
      "Total pixel count must be evenly divisible by segment count");
  static_assert(N % 2 == 0,
      "Segment count must be even (equal split across two arms)");
  static_assert(S >= N,
      "Must have at least one pixel per segment");
  static_assert((N & (N - 1)) == 0 && N <= 8,
      "N must be a power of two and <= 8: ID is decoded from up to 3 GPIO "
      "straps as (~raw) & (N-1), pins 21/22/23");

  // ── Pin assignments ─────────────────────────────────────────────────

  /**
   * @brief GPIO straps for hardware ID (active-low with internal pull-up).
   * Up to three straps decode N <= 8 segments; the build reads log2(N) of them.
   */
  static constexpr int PIN_ID0 = 21;
  static constexpr int PIN_ID1 = 22;
  static constexpr int PIN_ID2 = 23;

  /** @brief Number of ID straps actually read = log2(N). */
  static constexpr int ID_STRAPS = pov::segment_id_strap_count(N);

  /**
   * @brief Shared sync wire: master drives it OUTPUT (symbol bursts),
   *        downstream boards read it INPUT (RISING edge ISR). One pin serves
   *        both roles since a board is master XOR downstream.
   */
  static constexpr int PIN_FRAME_SYNC = 3;

  /**
   * @brief Master-enable strap for the external sync-out level shifter.
   *        OUTPUT, driven LOW on the master (segment 0) and HIGH elsewhere, so
   *        only the master drives the shared sync bus.
   */
  static constexpr int PIN_MASTER_EN = 5;

  // ── Flywheel timing ─────────────────────────────────────────────────

  /** @brief Nominal column period in µs (T0 ≈ 434 µs at 480 RPM × 288). */
  static constexpr float COLUMN_US =
      1000000.0f * 60.0f / (float(RPM) * float(CANVAS_W));

  /**
   * @brief Flywheel wake-up oversampling factor.
   *
   * Wakes are advisory (position comes from the cycle counter); a coarse grid
   * only quantizes when a column renders/emits. 8× per column keeps that lag
   * well inside the §5.2 self-censor budget. position() is time-derived and
   * tick() is idempotent/skip-tolerant, so the exact grid does not affect
   * correctness.
   */
  static constexpr int OVERSAMPLE = 8;

  static_assert(COLUMN_US / float(OVERSAMPLE) >= 1.0f,
      "Flywheel wake period must be >= 1 us for IntervalTimer::begin");

  /**
   * @brief HD107S SPI clock for the Phantasm DMA path, in Hz.
   *
   * 24 MHz (vs the 12 MHz default) halves the per-column transfer time. The
   * column ISR, not the strip, is the binding budget.
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
  HS_COLD_MEMBER POVSegmented() {
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
    Serial.print(y_step_ < 0 ? " (-y, reversed)" : " (+y)");
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

    pov::sync::Config cfg = pov::sync::phantasm_config(
        F_CPU, RPM, CANVAS_W, effect_count);
#ifdef HS_PROFILE_EPOCH_REVS
    // Profiling knob: stretch the epoch so one effect instance covers a full
    // preset cycle in a single capture.
    cfg.revs_per_effect = HS_PROFILE_EPOCH_REVS;
#endif
    HS_CHECK(cfg.valid(), "pov::sync::Config invariants");
    sync_.reconstruct(cfg);

    const bool master = (segment_id_ == 0);
    pinMode(PIN_MASTER_EN, OUTPUT);
    digitalWriteFast(PIN_MASTER_EN, master ? LOW : HIGH);
    if (master) {
      pinMode(PIN_FRAME_SYNC, OUTPUT);
      digitalWriteFast(PIN_FRAME_SYNC, LOW);
    } else {
      pinMode(PIN_FRAME_SYNC, INPUT);
      // Schmitt-trigger the sync input. The on-board divider + C_SYNC RC slows the
      // edge to reject BLDC/LED spikes; pad hysteresis then gives exactly one clean
      // interrupt per edge instead of multiple threshold recrossings on the slow
      // ramp. pinMode rewrites the pad-control register, so enable HYS afterward.
      *(portControlRegister(PIN_FRAME_SYNC)) |= IOMUXC_PAD_HYS;
    }

    sync_.seed(ARM_DWT_CYCCNT, master);

    if (!master) {
      attachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC),
                      sync_edge_isr, RISING);
    }
    HS_CHECK(timer_.begin(flywheel_isr, COLUMN_US / float(OVERSAMPLE)),
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
        handoff_.request_release();
        const unsigned long t0 = micros();
        while (!handoff_.release_complete()) {
          HS_CHECK(micros() - t0 < 100000UL,
                   "flywheel ISR failed to release the live effect");
        }
        handoff_.clear_pending();
        delete cur;
        // Restart the shared RNG stream per effect, seeded from the beacon-
        // synchronized effect index (spec §2): every board derives the same
        // per-visit stream locally, regardless of boot/join history — a board
        // wrong about the index is already building the wrong effect.
        const int32_t effect_index = pov::sync::SyncBoard::build_index_of(bw);
        hs::random().seed(hs::epoch_seed(static_cast<uint32_t>(effect_index)));
        cur = factories_[effect_index]();
        HS_CHECK(cur->height() == ROWS,
                 "POVSegmented: effect canvas height must equal S/2 (ROWS)");
        HS_CHECK(cur->width() == CANVAS_W,
                 "POVSegmented: effect canvas width must equal CANVAS_W");
        HS_CHECK(!cur->overrides_get_pixel(),
                 "POVSegmented: effect must not override get_pixel(); the "
                 "segmented render_column path bypasses it");
        // The first frame commits on a ZERO boundary, an arm-A-left window.
        clip_to_segment(cur, /*arm_a_left=*/true);
        cur->draw_frame();
        hs::disable_interrupts();
        // Publish under IRQ-off so the (effect, gen) pair reaches the ISR
        // atomically; publish()'s release store orders every constructor/
        // draw_frame() write before the ISR's acquire load.
        handoff_.publish(cur, gen);
        hs::enable_interrupts();
        built_gen = gen;
      }

      if (cur && handoff_.consumed(built_gen)) {
        // Render the quadrant the next display window paints: the live window's
        // opposite half (windows alternate ZERO/HALF).
        clip_to_segment(cur, handoff_.window_left() == 0);
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
        hs::disable_interrupts();
        const pov::sync::Telemetry tm = sync_.telemetry();
        hs::enable_interrupts();
        if (memcmp(&tm, &last_tm, sizeof tm) != 0) {
          // hs::log, not Serial.printf: Teensy's printf drags in newlib's float
          // formatter (~5 KB ITCM); these counters are all %lu.
          hs::log("sync acc=%lu rej=%lu inv=%lu cens=%lu abrt=%lu bdrop=%lu "
                  "bok=%lu brej=%lu fix=%lu rmis=%lu lock=%lu "
                  "flip=%lu coast=%lu epi=%lu",
                  (unsigned long)tm.symbols_accepted,
                  (unsigned long)tm.symbols_rejected_gate,
                  (unsigned long)tm.symbols_discarded_invalid,
                  (unsigned long)tm.emit_censored,
                  (unsigned long)tm.emit_aborted,
                  (unsigned long)tm.beacons_overrun_dropped,
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
   * @brief Samples the raw ID straps (ID_STRAPS bits, LSB = ID0).
   * @return Raw reading; floating (HIGH) bits set, grounded bits clear.
   */
  int sample_strap_() const {
    int raw = digitalReadFast(PIN_ID0);
    if constexpr (ID_STRAPS >= 2) raw |= digitalReadFast(PIN_ID1) << 1;
    if constexpr (ID_STRAPS >= 3) raw |= digitalReadFast(PIN_ID2) << 2;
    return raw;
  }

  /**
   * @brief Reads the hardware segment ID from the GPIO straps (log2(N) bits).
   *
   * Straps are INPUT_PULLUP; the raw reading is inverted, so a grounded strap
   * contributes a 1 and all-floating = ID 0 (master). A floating/cold strap
   * therefore elects a phantom second master and drives the push-pull sync wire
   * into contention — the triple-sample debounce below traps on that instability.
   *
   * Validates *this* board's strap only. Two failure modes exist; only the
   * first is trappable on the push-pull sync bus (U1 74AHCT125, ch C):
   *
   *   1. Unstable/cold strap on this board — a flaky link momentarily reads as
   *      ID 0 and enables this board's bus driver. The triple-sample debounce
   *      below catches the instability and traps.
   *
   *   2. A *peer* stably strapped to the same ID 0 — undetectable here. Both
   *      boards read a clean ID 0, both assert MASTER_EN, and both push-pull
   *      ch-C drivers source-fight on the shared wire. The two 100 ohm source
   *      resistors (R_S) current-limit the fight to ~22 mA (no part damage),
   *      but the bus parks near mid-rail and sync is silently corrupted — no
   *      trap fires, because each board's own strap reads valid.
   *
   * Mode 2 is an assembly/wiring fault the firmware cannot observe on this
   * board revision: the push-pull driver source-fights instead of wired-ANDing
   * (an open-drain bus would make a duplicate master benign and read-back
   * detectable, but that is a board respin, not the shipped design). It is
   * prevented procedurally — one board strapped master (all ID straps open) and
   * unique soldered ID links elsewhere, per PCB rules R-ID-2 (soldered links)
   * and R-ID-4 (silkscreen truth table).
   */
  void read_id() {
    pinMode(PIN_ID0, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 2) pinMode(PIN_ID1, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 3) pinMode(PIN_ID2, INPUT_PULLUP);
    delay(10);  // settle time for pull-ups

    // Debounce: three samples ~5 ms apart must agree; an unstable strap reads as
    // a second master and drives the push-pull sync wire into bus contention.
    const int raw0 = sample_strap_();
    for (int i = 0; i < 2; ++i) {
      delay(5);
      HS_CHECK(sample_strap_() == raw0,
               "unstable segment-ID strap (field/manufacturing fault)");
    }

    // Invert the reading (all-floating pull-ups => ID 0), then mask to log2(N)
    // bits; the mask is load-bearing.
    segment_id_ = pov::decode_segment_id(raw0, N);
  }

  // ── Segment mapping ─────────────────────────────────────────────────

  /**
   * @brief Computes the precomputed ISR mapping from hardware segment ID.
   * @details IDs [0, N/2) map to arm A and [N/2, N) map to arm B. Each arm's
   * northern bands advance in +y; its southern bands advance in -y.
   */
  void configure_segment() {
    const pov::SegmentMap m = pov::segment_map(segment_id_, S, N);
    arm_b_  = m.arm_b;
    y_base_ = m.y_base;
    y_step_ = m.y_step;
  }

  /**
   * @brief Clip @p e to this segment's quadrant for the upcoming display window.
   * @param e Effect to clip.
   * @param arm_a_left True if the window this frame displays in sweeps arm-A
   *        columns [0, CANVAS_W/2); arm B paints the opposite half.
   * @details Left full-canvas for an effect that reads cross-segment or prior-
   *          frame state (needs_full_frame / persists_pixels), so trails and
   *          feedback stay correct under the per-frame arm-half alternation.
   */
  void clip_to_segment(Effect *e, bool arm_a_left) {
    if (e->needs_full_frame() || e->persists_pixels())
      return;
    const pov::SegmentMap m{arm_b_, y_base_, y_step_};
    const pov::SegmentClip c = pov::segment_clip(m, arm_a_left, S, N, CANVAS_W);
    e->set_clip(c.y0, c.y1, c.x0, c.x1);
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

  // render_column samples the raw display buffer, bypassing get_pixel(), so a
  // get_pixel-overriding effect must never go live on either adoption path.
  static FASTRUN void assert_render_column_safe(Effect *p) {
    HS_CHECK(!p->overrides_get_pixel(),
             "segmented render_column bypasses get_pixel(); a get_pixel-"
             "overriding effect cannot run on this path");
  }

  /**
   * @brief Flywheel ISR: the sole owner of all sync state (spec §8).
   *
   * Paced by an IntervalTimer at T0/OVERSAMPLE as a wake-up only — the
   * cycle counter decides which column it is (spec §4.1), so a late, early,
   * or coalesced wake-up cannot inject drift: the ISR is idempotent when the
   * column is unchanged and skip-tolerant when it jumped (the skipped
   * columns were masked precisely because the strip was busy).
   */
  static FASTRUN void flywheel_isr() {
    HS_ISR_PROFILE(hs::g_flywheel_wake_cycles);
    const uint32_t now = ARM_DWT_CYCCNT;

    // Complete a deferred dark-path sync pulse from the previous wake: the
    // dark-latched path's body is too short to width a same-wake pulse to spec
    // §5.2's "tens of µs," so it holds the pin HIGH and drops it here instead.
    if (sync_low_pending_) {
      digitalWriteFast(PIN_FRAME_SYNC, LOW);
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
      digitalWriteFast(PIN_FRAME_SYNC, HIGH);

    // Release handshake: the foreground wants the live pointer dropped so it
    // can destroy the instance (epoch teardown / beacon rebuild).
    handoff_.service_release();

    // Swap in the foreground-constructed pending effect, only ever at a ZERO
    // boundary. Two paths: commit (the B+K epoch deadline, spec §6.1) and join
    // (first display: boot / beacon join / index correction, taken at the next
    // join-grid boundary so all boards go live at the same crossing).
    if (a.commit) {
      const auto p = handoff_.pending_acquire();
      HS_CHECK(handoff_.committable(
                   p, pov::sync::SyncBoard::build_gen_of(sync_.build_word())),
               "epoch commit: effect init exceeded the K-revolution window");
      assert_render_column_safe(p.effect);
      handoff_.adopt(p.effect, p.gen);
    } else if (a.join_boundary && !a.dark && handoff_.live() == nullptr) {
      // Adopt only an effect still matching the wire's advertised generation; a
      // visibility lag that fails the match simply joins one grid step later.
      const auto p = handoff_.pending_acquire();
      if (handoff_.joinable(
              p, pov::sync::SyncBoard::build_gen_of(sync_.build_word()))) {
        assert_render_column_safe(p.effect);
        handoff_.adopt(p.effect, p.gen);
      }
    }

    // Publish which half the now-open display window sweeps so the foreground
    // clips the next frame to the quadrant this segment will paint: a ZERO flip
    // opens the arm-A-left [0,W/2) half-rev, a HALF flip opens [W/2,W).
    if (a.flip)
      handoff_.set_window_left(a.zero_crossing);

    // Flip whenever the effect is live, even during the dark commit window:
    // advance_display() is what releases a foreground blocked in buffer_free().
    Effect *e = handoff_.live();
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
        digitalWriteFast(PIN_FRAME_SYNC, LOW);
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
    {
      HS_ISR_PROFILE(hs::g_column_pack_cycles);
      int y = y_base_;
      for (int i = 0; i < PPS; ++i, y += y_step_) {
        frame.packPixel(i, buf[y * w + x_col]);
      }
    }
    // Drop the accept/overrun result: a dropped image column self-heals next
    // tick (render_black, by contrast, gates on it for the fail-dark latch).
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
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
   * @brief Effect handoff state machine between the foreground and the ISR.
   * @details The teardown handshake, the acquire/release publish/adopt of the
   *          pending effect, the consumed-generation gate, and the display-window
   *          alternation live in pov_handoff.h (host-tested). Ownership: the
   *          foreground constructs and deletes; the ISR only ever dereferences
   *          the instance it has been handed via live().
   */
  static pov::EffectHandoff<Effect> handoff_;
  static bool dark_latched_;               /**< True once the black frame has latched; ISR-owned. */
  static bool sync_low_pending_;           /**< ISR-owned: dark-path pulse drop deferred to next wake. */

  static int segment_id_;                  /**< Decoded hardware segment ID (up to 3 strap bits, 0..N-1). */
  static bool arm_b_;                      /**< True if this segment lives on arm B (x + W/2). */
  static int y_base_;                      /**< Canvas row of this segment's LED 0.      */
  static int y_step_;                      /**< Row stride per LED: +1 north band or -1 reversed south band. */
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
pov::EffectHandoff<Effect> POVSegmented<S, N, RPM>::handoff_;

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
// DTCM regardless. Each instantiating target must instead invoke
// HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM) once at file scope; it emits
// the required explicit specialization, whose ordinary strong linkage keeps the
// DMAMEM section attribute — see Phantasm.ino.
#define HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM)                       \
  template <>                                                                  \
  DMAMEM DMALEDController<(S) / (N)> POVSegmented<S, N, RPM>::ledController_{    \
      POVSegmented<S, N, RPM>::SPI_CLOCK_HZ}
#endif

#endif // ARDUINO
