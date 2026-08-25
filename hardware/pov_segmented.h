/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
 * (docs/specs/phantasm_frame_sync_spec.md):
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
#include "pov_single_map.h"  // pov::transfer_us (host-tested)
#include "pov_sync.h"    // pure sync protocol (host-testable; see that file)
#include "pov_handoff.h" // pure effect-handoff state machine (host-testable)
#include "pov_submit_gate.h" // pure LED-submit decision (host-testable)

#ifdef ARDUINO
#include <Arduino.h>

#ifdef USE_DMA_LEDS
#include "dma_led.h"
#else
#error                                                                         \
    "POVSegmented requires USE_DMA_LEDS (the Phantasm DMA LED transport): FastLED's bit-bang show() masks IRQs for windows that break the sync symbol margins (pitch > M, gap timeout > pitch + M; spec §5.2)."
#endif

#include "render/canvas.h"
#include "math/geometry.h"
#include "engine/memory.h"

#include <atomic>
#include <cstring>
#include <type_traits>

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
template <int S, int N, int RPM> class POVSegmented {

  // ── Compile-time geometry ───────────────────────────────────────────

  static constexpr int PPS = S / N;  /**< Pixels per segment. */
  static constexpr int ROWS = S / 2; /**< Canvas rows (height). */

  static_assert(RPM > 0, "RPM must be positive (COLUMN_US divides by RPM)");
  static_assert(S % N == 0,
                "Total pixel count must be evenly divisible by segment count");
  static_assert(N % 2 == 0,
                "Segment count must be even (equal split across two arms)");
  static_assert(S >= N, "Must have at least one pixel per segment");
  static_assert(ROWS == CANVAS_H,
                "S/2 must equal CANVAS_H: effects render the full CANVAS_W x "
                "CANVAS_H canvas and the ISR packs S/2 rows");
  static_assert(
      (N & (N - 1)) == 0 && N <= 8,
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

  static_assert(
      COLUMN_US / float(OVERSAMPLE) >= 1.0f,
      "Flywheel wake period must be >= 1 us for IntervalTimer::begin");

  /**
   * @brief NVIC priority for the sync-wire edge IRQ (Teensy 4 pin interrupts
   *        all share IRQ_GPIO6789).
   *
   * The edge ISR timestamps the edge with ARM_DWT_CYCCNT, and snap() re-bases
   * the flywheel epoch onto that stamp — so it must preempt the flywheel ISR or
   * the stamp becomes a service time up to one full column-ISR body late. The
   * Cortex-M7 NVIC implements 4 priority bits (levels in steps of 16, lower =
   * higher priority) and Teensy defaults every IRQ, IntervalTimer included, to
   * 128; 16 puts the edge above all of them while leaving 0 free.
   */
  static constexpr uint8_t SYNC_EDGE_IRQ_PRIORITY = 16;

  /**
   * @brief HD107S SPI clock for the Phantasm DMA path, in Hz.
   *
   * 24 MHz (vs the 12 MHz default) halves the per-column transfer time. The
   * column ISR, not the strip, is the binding budget.
   */
  static constexpr uint32_t SPI_CLOCK_HZ = 24000000;

  /**
   * @brief Worst-case duration of one column's LED transfer, in µs.
   * @details Image frame plus the trailing black frame strobe_columns() appends,
   * eight clocks per byte at SPI_CLOCK_HZ, rounded up so the overrun check
   * below never under-counts.
   */
  static constexpr unsigned long COLUMN_TRANSFER_US =
      pov::transfer_us(HD107SFrame<PPS>::COMPOSITE_SIZE, SPI_CLOCK_HZ);

  // A transfer wider than the column period overruns every column; the submit
  // gate absorbs the drops, so it surfaces as a dim image, not as a fault.
  static_assert(COLUMN_TRANSFER_US < COLUMN_US,
                "LED transfer outlasts the column period (S, N, RPM and "
                "SPI_CLOCK_HZ overrun the DMA every column)");

public:
  /** @brief Foreground effect constructor: builds, arena-configures, and
   *  init()s one roster entry, ready to draw its first frame. */
  using EffectFactory = Effect *(*)();

  /**
   * @brief Drives MASTER_EN to its disabled level, parking the external
   *        sync-out buffer.
   * @details Call as the first statement of setup(). MASTER_EN is the '125
   *          channel-C /OE and only takes its board-role level in run_show();
   *          until this runs, only the R_MEN pull-up (PCB rule R-LS-5) holds a
   *          board's sync driver off the shared bus, and a driver that is
   *          enabled while a peer drives the bus source-fights undetectably
   *          (see read_id()).
   */
  HS_COLD_MEMBER static void park_sync_out() {
    digitalWriteFast(PIN_MASTER_EN, HIGH);
    pinMode(PIN_MASTER_EN, OUTPUT);
  }

  /**
   * @brief Hardware segment ID decoded from the GPIO straps.
   * @return Segment index in [0, N); 0 is the master.
   * @details Meaningful only once an instance exists — the constructor's
   *          read_id() is what samples the straps.
   */
  static int segment_index() { return segment_id; }

  /**
   * @brief Initializes hardware: reads segment ID and configures the LED
   *        driver.
   * @details CONTRACT — construct only from setup(), never as a file-scope
   *          global. This constructor performs hardware I/O directly: read_id()
   *          does pinMode/delay, and the body then brings up the DMA LED
   *          driver, enables the DWT cycle counter, and prints to Serial — all
   *          valid only once the Arduino core is initialized. dma_led.h
   *          deliberately keeps hardware bring-up out of its constructor (an
   *          explicit begin()) and warns against constructor-time I/O; this
   *          class diverges on purpose because its sole instantiation site is
   *          the Phantasm setup() (a function-local singleton), so the "after
   *          core init" precondition always holds. Do not promote this object
   *          to a global or construct it before setup().
   */
  HS_COLD_MEMBER POVSegmented() {
    read_id();
    configure_segment();

    ledController.begin();
    ledController.setCorrection(hd107s::TYPICAL_LED_STRIP.r,
                                hd107s::TYPICAL_LED_STRIP.g,
                                hd107s::TYPICAL_LED_STRIP.b);
    ledController.setTemperature(hd107s::CANDLE.r, hd107s::CANDLE.g,
                                 hd107s::CANDLE.b);
    ledController.setBrightness(255);

    // Enable the DWT cycle counter the flywheel timebase reads: TRCENA gates the
    // DWT block on, then CYCCNTENA starts the counter.
    ARM_DEMCR |= ARM_DEMCR_TRCENA;
    ARM_DWT_CTRL |= ARM_DWT_CTRL_CYCCNTENA;

    Serial.print("[Phantasm] Segment ");
    Serial.print(segment_id);
    Serial.print(segment.arm_b ? " arm-B" : " arm-A");
    Serial.print(segment.y_step < 0 ? " (-y, reversed)" : " (+y)");
    Serial.print(segment_id == 0 ? " MASTER" : "");
    Serial.print(" | y_base=");
    Serial.print(segment.y_base);
    Serial.print(" y_step=");
    Serial.print(segment.y_step);
    Serial.print(" pixels=");
    Serial.println(PPS);
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
   * @tparam R            Roster length, deduced from `factories`.
   * @param factories     One constructor per roster entry (HS_EFFECT_LIST
   *                      order — identical on every board). Its length is the
   *                      roster length.
   * @param effect_revolutions Optional duration for each roster entry, in
   *                           revolutions.
   * @param stable_effect_seeds Optional RNG identity for each roster entry;
   *                           without it the per-visit stream is seeded from
   *                           the roster index alone.
   * @details Both optional tables are taken as pointers to arrays of R, so a
   * table that does not span the roster is a compile error rather than an
   * unguarded index: revolutions_for_effect() and Config::valid()'s own sweep
   * read effect_revolutions over [0, R), and the per-effect reseed below reads
   * stable_effect_seeds at the published index.
   */
  template <int R>
  [[noreturn]] void
  run_show(const EffectFactory (&factories)[R],
           const uint32_t (*effect_revolutions)[R] = nullptr,
           const uint64_t (*stable_effect_seeds)[R] = nullptr) {
    effect_factories = factories;
    effect_seed_identities =
        stable_effect_seeds ? &(*stable_effect_seeds)[0] : nullptr;

    // F_CPU_ACTUAL, not F_CPU: the flywheel timebase counts ARM_DWT_CYCCNT
    // ticks, which run at the clock the core actually booted to.
    pov::sync::Config cfg =
        pov::sync::phantasm_config(F_CPU_ACTUAL, RPM, CANVAS_W, R);
    cfg.effect_revolutions =
        effect_revolutions ? &(*effect_revolutions)[0] : nullptr;
#ifdef HS_PROFILE_EPOCH_REVS
    // Profiling knob: stretch the epoch so one effect instance covers a full
    // preset cycle in a single capture.
    cfg.revs_per_effect = HS_PROFILE_EPOCH_REVS;
    cfg.effect_revolutions = nullptr;
#endif
    const char *const bad_invariant = cfg.valid();
    HS_CHECK(bad_invariant == nullptr, "pov::sync::Config invariant: %s",
             bad_invariant);
    sync.configure(cfg);
    gap_timeout_cycles = sync.gap_timeout_cycles();
    max_burst_cycles = sync.max_burst_cycles();
    glitch_filter_cycles = sync.glitch_filter_cycles();

    const bool master = (segment_id == 0);
    if (master) {
      digitalWriteFast(PIN_FRAME_SYNC, LOW);
      pinMode(PIN_FRAME_SYNC, OUTPUT);
    } else {
      pinMode(PIN_FRAME_SYNC, INPUT);
      // Schmitt-trigger the sync input. The on-board divider + C_SYNC RC slows the
      // edge to reject BLDC/LED spikes; pad hysteresis then gives exactly one clean
      // interrupt per edge instead of multiple threshold recrossings on the slow
      // ramp. pinMode rewrites the pad-control register, so enable HYS afterward.
      *(portControlRegister(PIN_FRAME_SYNC)) |= IOMUXC_PAD_HYS;
    }
    // park_sync_out() already left MASTER_EN an output at its disabled level, so
    // this write is what enables the sync-bus driver: take the board-role level
    // only once PIN_FRAME_SYNC is driven, or a pad keeper puts one spurious edge
    // on the wire that downstream boards read as a symbol.
    digitalWriteFast(PIN_MASTER_EN, master ? LOW : HIGH);

    sync.seed(ARM_DWT_CYCCNT, master);

    if (!master) {
      attachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC), sync_edge_isr,
                      RISING);
      // attachInterrupt leaves the IRQ at the Teensy default (128), equal to
      // the flywheel's; raise it so the edge stamp is taken at the edge.
      NVIC_SET_PRIORITY(IRQ_GPIO6789, SYNC_EDGE_IRQ_PRIORITY);
    }
    HS_CHECK(timer.begin(flywheel_isr, COLUMN_US / float(OVERSAMPLE)),
             "flywheel IntervalTimer failed to start (no PIT channel)");

    // ── Foreground: construct effects on request, render, report ────────
    Effect *cur = nullptr;
    uint32_t built_gen = 0;
    uint32_t last_overrun = ledController.getOverrunCount();
    pov::sync::Telemetry last_tm{};
    unsigned long last_report = millis();
    // The K-revolution construction budget, in the flywheel's own timebase.
    const uint32_t commit_budget_cycles =
        cfg.col_cycles(static_cast<int32_t>(cfg.commit_revs) * CANVAS_W);
    const uint32_t cycles_per_us = F_CPU_ACTUAL / 1000000u;
    // ARM_DWT_CYCCNT, not micros(): this stamp is taken every spin and
    // micros() brackets its read with __disable_irq().
    uint32_t poll_prev_cycles = ARM_DWT_CYCCNT;

    for (;;) {
      const uint32_t poll_cycles = ARM_DWT_CYCCNT;
      const uint32_t bw = sync.build_word();
      const uint32_t gen = pov::sync::SyncBoard::build_gen_of(bw);
      if (gen != built_gen) {
        const uint32_t build_start_cycles = ARM_DWT_CYCCNT;
        // The ISR owns the live-effect pointer; ask it to let go before the
        // old instance is destroyed (it acknowledges within one wake-up).
        handoff.request_release();
        const unsigned long t0 = micros();
        while (!handoff.release_complete()) {
          HS_CHECK(micros() - t0 < 100000UL,
                   "flywheel ISR failed to release the live effect");
        }
        handoff.clear_pending();
        delete cur;
        // Restart the shared RNG stream per effect, seeded from the beacon-
        // synchronized effect index (spec §2): every board derives the same
        // per-visit stream locally, regardless of boot/join history — a board
        // wrong about the index is already building the wrong effect.
        const int32_t effect_index = pov::sync::SyncBoard::build_index_of(bw);
        HS_CHECK(effect_index >= 0 && effect_index < R,
                 "sync published an out-of-roster effect index");
        hs::random().seed(
            effect_seed_identities
                ? effect_seed_identities[effect_index]
                : hs::epoch_seed(static_cast<uint32_t>(effect_index)));
        cur = effect_factories[effect_index]();
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
        cur->set_buffer_ready_hook(prepare_segment_clip);
        hs::disable_interrupts();
        // Publish under IRQ-off so the (effect, gen) pair reaches the ISR
        // atomically; publish()'s release store orders every constructor/
        // draw_frame() write before the ISR's acquire load.
        handoff.publish(cur, gen);
        hs::enable_interrupts();
        built_gen = gen;
        // Per-build budget report. The ISR's commit trap is the only other
        // signal that this window ran tight, and it fires on every board at
        // once. The request cannot predate the previous poll, so `window`
        // covers the render that was in flight plus the build, and `build`
        // is the construction cost alone; `margin` is the headroom left.
        const uint32_t done_cycles = ARM_DWT_CYCCNT;
        const unsigned long window_us =
            (done_cycles - poll_prev_cycles) / cycles_per_us;
        const unsigned long budget_us = commit_budget_cycles / cycles_per_us;
        if (hs::debug)
          hs::log("build idx=%ld build=%lu window=%lu budget=%lu margin=%ld us",
                  (long)effect_index,
                  (unsigned long)((done_cycles - build_start_cycles) /
                                  cycles_per_us),
                  window_us, budget_us, (long)budget_us - (long)window_us);
      }
      poll_prev_cycles = poll_cycles;

      if (cur && handoff.consumed(built_gen)) {
        const unsigned long f0 = micros();
        cur->draw_frame();
        if (hs::debug) {
          Serial.print("ft ");
          Serial.println(micros() - f0);
        }
      }

      // Health telemetry (spec §8.6): foreground-polled, emitted only on change.
      if (hs::debug && millis() - last_report >= 1000UL) {
        last_report = millis();
        const pov::sync::Telemetry tm = sync.telemetry_snapshot();
        // Change detection is a byte compare, so Telemetry must have no
        // padding: a member of another width would make stale padding read as
        // a change and spam the log every poll.
        static_assert(
            std::has_unique_object_representations_v<pov::sync::Telemetry>,
            "Telemetry must be padding-free");
        if (std::memcmp(&tm, &last_tm, sizeof tm) != 0) {
          // hs::log, not Serial.printf: Teensy's printf drags in newlib's float
          // formatter (~5 KB ITCM); these counters are all %lu.
          // Two lines, not one: 18 saturating uint32 counters run to 10 digits
          // each, which overruns hs::log's fixed 256-byte buffer and truncates
          // the tail. Each line below fits its own worst case.
          hs::log("sync coast=%lu stall=%lu epi=%lu lock=%lu flip=%lu acc=%lu "
                  "rej=%lu inv=%lu",
                  (unsigned long)tm.max_coast_halves,
                  (unsigned long)tm.master_stalls,
                  (unsigned long)tm.epochs_refractory_ignored,
                  (unsigned long)tm.lock_transitions, (unsigned long)tm.flips,
                  (unsigned long)tm.symbols_accepted,
                  (unsigned long)tm.symbols_rejected_gate,
                  (unsigned long)tm.symbols_discarded_invalid);
          hs::log(
              "sync emit cens=%lu abrt=%lu bdrop=%lu bbusy=%lu blate=%lu "
              "sdrop=%lu bok=%lu brej=%lu fix=%lu rmis=%lu",
              (unsigned long)tm.emit_censored, (unsigned long)tm.emit_aborted,
              (unsigned long)tm.beacons_overrun_dropped,
              (unsigned long)tm.beacons_busy_dropped,
              (unsigned long)tm.beacons_late_dropped,
              (unsigned long)tm.boundary_bursts_dropped,
              (unsigned long)tm.beacons_ok, (unsigned long)tm.beacons_rejected,
              (unsigned long)tm.beacon_index_corrections,
              (unsigned long)tm.beacon_rev_mismatches);
          last_tm = tm;
        }
        const uint32_t overruns = ledController.getOverrunCount();
        if (overruns != last_overrun) {
          Serial.print("overrun ");
          Serial.println(overruns);
          last_overrun = overruns;
        }
      }
    }
  }

private:
  // ── Hardware ID ─────────────────────────────────────────────────────

  /**
   * @brief Samples the raw ID straps (ID_STRAPS bits, LSB = ID0).
   * @return Raw reading; floating (HIGH) bits set, grounded bits clear.
   */
  HS_COLD_MEMBER static int sample_strap() {
    int raw = digitalReadFast(PIN_ID0);
    if constexpr (ID_STRAPS >= 2)
      raw |= digitalReadFast(PIN_ID1) << 1;
    if constexpr (ID_STRAPS >= 3)
      raw |= digitalReadFast(PIN_ID2) << 2;
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
  HS_COLD_MEMBER static void read_id() {
    pinMode(PIN_ID0, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 2)
      pinMode(PIN_ID1, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 3)
      pinMode(PIN_ID2, INPUT_PULLUP);
    delay(10); // settle time for pull-ups

    // Debounce: three samples ~5 ms apart must agree; an unstable strap reads as
    // a second master and drives the push-pull sync wire into bus contention.
    const int raw0 = sample_strap();
    for (int i = 0; i < 2; ++i) {
      delay(5);
      HS_CHECK(sample_strap() == raw0,
               "unstable segment-ID strap (field/manufacturing fault)");
    }

    // Invert the reading (all-floating pull-ups => ID 0), then mask to log2(N)
    // bits; the mask is load-bearing.
    segment_id = pov::decode_segment_id(raw0, N);
  }

  // ── Segment mapping ─────────────────────────────────────────────────

  /**
   * @brief Computes the precomputed ISR mapping from hardware segment ID.
   * @details IDs [0, N/2) map to arm A and [N/2, N) map to arm B. Each arm's
   * northern bands advance in +y; its southern bands advance in -y.
   */
  static void configure_segment() {
    segment = pov::segment_map(segment_id, S, N);
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
  static void clip_to_segment(Effect *e, bool arm_a_left) {
    if (e->needs_full_frame() || e->persists_pixels())
      return;
    const pov::SegmentClip c =
        pov::segment_clip(segment, arm_a_left, S, N, CANVAS_W);
    e->set_clip(c.y0, c.y1, c.x0, c.x1);
  }

  static void prepare_segment_clip(Effect &e) {
    // One window ahead: the frame being drawn displays in the window after the
    // open one, which sweeps the opposite half.
    clip_to_segment(&e, handoff.window_left() == 0);
  }

  // ── ISRs ────────────────────────────────────────────────────────────

  /**
   * @brief Sync-wire edge ISR (downstream boards only): a pure publisher.
   *
   * Applies the glitch filter and records the edge into the mailbox; touches
   * no flywheel, flip, or epoch state (spec §8.2 single-writer model). It runs
   * above the flywheel ISR (SYNC_EDGE_IRQ_PRIORITY) so the captured cycle count
   * is the edge time rather than a service time; that preemption is also what
   * makes the consumer's interrupt-save bracket around the mailbox claim
   * load-bearing.
   */
  static void sync_edge_isr() { sync.on_sync_edge(ARM_DWT_CYCCNT); }

  /**
   * @brief Flywheel ISR: the sole owner of all sync state (spec §8).
   *
   * Paced by an IntervalTimer at T0/OVERSAMPLE as a wake-up only — the
   * cycle counter decides which column it is (spec §4.1), so a late, early,
   * or coalesced wake-up cannot inject drift: the ISR is idempotent when the
   * column is unchanged and skip-tolerant when it jumped (the skipped
   * columns were masked precisely because the strip was busy).
   */
  static void flywheel_isr() {
    HS_ISR_PROFILE(hs::g_flywheel_wake_cycles);
    const uint32_t now = ARM_DWT_CYCCNT;
    pov::run_wake_sequence(
        sync_pulse, submit_gate, handoff,
        [&] {
          pov::sync::BurstSnapshot burst;
          const pov::sync::BurstSnapshot *bp = nullptr;
          if (segment_id != 0) {
            const uint32_t primask = hs::save_disable_interrupts();
            if (sync.mailbox().try_claim(now, gap_timeout_cycles,
                                         max_burst_cycles, &burst))
              bp = &burst;
            sync.mailbox().age_prior(now, glitch_filter_cycles);
            hs::restore_interrupts(primask);
          }
          return sync.tick(now, bp);
        },
        [] { return pov::sync::SyncBoard::build_gen_of(sync.build_word()); },
        [](bool high) { digitalWriteFast(PIN_FRAME_SYNC, high ? HIGH : LOW); },
        [] {
          HS_CHECK(
              false,
              "epoch commit: effect init exceeded the K-revolution window");
        },
        [](Effect *e, int32_t column) {
          e->set_output_envelope(sync.effect_envelope(column, CANVAS_W));
        },
        [](pov::SubmitAction action, Effect *e, int32_t column) {
          switch (action) {
          case pov::SubmitAction::BLACK:
            return render_black();
          case pov::SubmitAction::COLUMN:
            return render_column(e, column);
          case pov::SubmitAction::RESUBMIT:
            return resubmit_frame(e->strobe_columns());
          case pov::SubmitAction::NONE:
            break;
          }
          return false;
        });
  }

  // ── Pixel packing ───────────────────────────────────────────────────

  /**
   * @brief Packs this segment's pixels for canvas column @p x and submits.
   * @param e Live effect supplying the display buffer to sample.
   * @param x Canvas column index in [0, CANVAS_W); arm-B segments sample
   *          column x + W/2.
   * @return true if the LED transport accepted the frame; false if it was
   *         dropped on a DMA overrun (caller retries via resubmit_frame()).
   * @details The loop is branchless — all per-segment decisions are resolved at
   *          boot in configure_segment().  Arm B segments read from x + W/2
   *          (opposite half of the image).
   */
  [[nodiscard]] HS_O3_FN static bool render_column(Effect *e, int x) {
    const int w = e->width();
    const int x_col = pov::segment_x_col(segment.arm_b, x, w);

    // ISR fast path: index the display buffer directly, dropping PPS per-pixel
    // virtual get_pixel() dispatches. No effect on this path overrides get_pixel.
    const Pixel *buf = e->display_buffer();

    auto &frame = ledController.backFrame();
    {
      HS_ISR_PROFILE(hs::g_column_pack_cycles);
      const int stride = pov::segment_row_stride(segment, w);
      int off = pov::segment_pixel_base(segment, x_col, w);
      if (e->output_envelope_u16() == 65535u)
        for (int i = 0; i < PPS; ++i, off += stride)
          frame.packPixel(i, buf[off]);
      else
        for (int i = 0; i < PPS; ++i, off += stride)
          frame.packPixel(i, e->apply_output_envelope(buf[off]));
    }
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submitFrame(e->strobe_columns());
  }

  /**
   * @brief Re-submits the frame a previous wake had dropped on overrun.
   * @param strobe Whether the live effect wants the trailing black frame.
   * @return true if the transport accepted it this time.
   * @details No repack: submitFrame() returns before swapping buffers on an
   *          overrun, so backFrame() still holds the dropped column's pixels.
   */
  [[nodiscard]] static bool resubmit_frame(bool strobe) {
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submitFrame(strobe);
  }

  /**
   * @brief Submits one all-black frame (ACQUIRE / construction window).
   * @return true if the black frame was accepted by the LED transport; false
   *         if it was dropped on a DMA overrun (caller must retry, not latch).
   */
  [[nodiscard]] static bool render_black() {
    auto &frame = ledController.backFrame();
    {
      HS_ISR_PROFILE(hs::g_column_pack_cycles);
      for (int i = 0; i < PPS; ++i) {
        frame.packPixel(i, Pixel(0, 0, 0));
      }
    }
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submitFrame(false);
  }

  // ── Static state ────────────────────────────────────────────────────

  /**
   * @brief The sync engine: sole owner of all sync/flywheel state.
   * @details Written ONLY by the flywheel ISR (tick()) and the edge ISR
   *          (mailbox publisher) per the spec §8 single-writer model. Foreground
   *          reads are single aligned words (build_word) or debug telemetry.
   */
  static pov::sync::SyncBoard sync;
  static uint32_t gap_timeout_cycles;
  static uint32_t max_burst_cycles;
  static uint32_t glitch_filter_cycles;
  static IntervalTimer timer; /**< Flywheel wake-up timer (PIT channel).   */
  /** @brief Roster of effect constructors (HS_EFFECT_LIST order). */
  static const EffectFactory *effect_factories;
  static const uint64_t *effect_seed_identities;

  /**
   * @brief Effect handoff state machine between the foreground and the ISR.
   * @details The teardown handshake, the acquire/release publish/adopt of the
   *          pending effect, the consumed-generation gate, and the display-window
   *          alternation live in pov_handoff.h (host-tested). Ownership: the
   *          foreground constructs and deletes; the ISR only ever dereferences
   *          the instance it has been handed via live().
   */
  static pov::EffectHandoff<Effect> handoff;
  /**
   * @brief Per-wake LED-submit decision and its overrun-retry latches.
   * @details The accept/drop bookkeeping lives in pov_submit_gate.h
   *          (host-tested); the ISR only performs the action it names.
   */
  static pov::SubmitGate submit_gate;
  /**
   * @brief Sync-pulse width decision and its deferred-drop latch.
   * @details Lives in pov_submit_gate.h (host-tested); the ISR only performs
   *          the pin writes it names.
   */
  static pov::SyncPulseGate sync_pulse;

  static int
      segment_id; /**< Decoded hardware segment ID (up to 3 strap bits, 0..N-1). */
  static pov::SegmentMap segment; /**< Precomputed canvas mapping. */
  static DMALEDController<PPS>
      ledController; /**< DMA SPI LED controller for the segment strip. */
};

// ── Static member definitions ───────────────────────────────────────────

template <int S, int N, int RPM>
pov::sync::SyncBoard POVSegmented<S, N, RPM>::sync{pov::sync::Config{}};

template <int S, int N, int RPM>
uint32_t POVSegmented<S, N, RPM>::gap_timeout_cycles = 0;

template <int S, int N, int RPM>
uint32_t POVSegmented<S, N, RPM>::max_burst_cycles = 0;

template <int S, int N, int RPM>
uint32_t POVSegmented<S, N, RPM>::glitch_filter_cycles = 0;

template <int S, int N, int RPM> IntervalTimer POVSegmented<S, N, RPM>::timer;

template <int S, int N, int RPM>
const typename POVSegmented<S, N, RPM>::EffectFactory
    *POVSegmented<S, N, RPM>::effect_factories = nullptr;

template <int S, int N, int RPM>
const uint64_t *POVSegmented<S, N, RPM>::effect_seed_identities = nullptr;

template <int S, int N, int RPM>
pov::EffectHandoff<Effect> POVSegmented<S, N, RPM>::handoff;

template <int S, int N, int RPM>
pov::SubmitGate POVSegmented<S, N, RPM>::submit_gate;

template <int S, int N, int RPM>
pov::SyncPulseGate POVSegmented<S, N, RPM>::sync_pulse;

template <int S, int N, int RPM> int POVSegmented<S, N, RPM>::segment_id = 0;

template <int S, int N, int RPM>
pov::SegmentMap POVSegmented<S, N, RPM>::segment{false, 0, 1};

// ledController has no out-of-line definition: DMAMEM survives only on an
// explicit specialization (see DMALEDController in dma_led_controller.h). Each
// instantiating target invokes HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM)
// once at file scope — see targets/Phantasm/phantasm_target.h.
#define HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM)                      \
  template <>                                                                  \
  DMAMEM DMALEDController<(S) / (N)> POVSegmented<S, N, RPM>::ledController {  \
    POVSegmented<S, N, RPM>::SPI_CLOCK_HZ                                      \
  }

#endif // ARDUINO
