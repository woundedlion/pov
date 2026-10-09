/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
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
 *     Segment 0 (top):    LED 0 at N end (y=0)   → LED 71 at junction (y=71)
 *     Segment 1 (bottom):  LED 0 at S end (y=143) → LED 71 at junction (y=72)  ← reversed
 *
 *   Arm B (x offset by W/2):
 *     Segment 2 (top):    LED 0 at N end (y=0)   → LED 71 at junction (y=71)
 *     Segment 3 (bottom):  LED 0 at S end (y=143) → LED 71 at junction (y=72)  ← reversed
 *
 * Each Teensy reads a hardware ID from GPIO straps at boot to determine
 * which segment it owns. Segment 0 (the master) emits count-coded symbol
 * bursts on one shared sync wire (docs/specs/phantasm_frame_sync_spec.md);
 * every board generates its columns from a local flywheel timebase derived
 * from the free-running cycle counter, and the symbols snap each flywheel's
 * phase and synchronize buffer flips and the effect playlist. This file is
 * the device shell around pov::sync::SyncBoard.
 *
 * Effects use full-canvas coordinates with rendering clipped per board.
 */
#pragma once
#include "core/platform/led.h"
#include "pov_segment_map.h"
#include "pov_segment_frame.h"
#include "pov_sync.h"
#include "pov_handoff.h"
#include "pov_submit_gate.h"

#ifdef ARDUINO
#include <Arduino.h>

#ifndef HS_PHANTASM_BOARD_REV
#error "Define HS_PHANTASM_BOARD_REV as 11 (rev 1.1) or 12 (rev 1.2)"
#elif HS_PHANTASM_BOARD_REV != 11 && HS_PHANTASM_BOARD_REV != 12
#error "Unsupported HS_PHANTASM_BOARD_REV; supported revisions are 11 and 12"
#endif

#ifdef USE_DMA_LEDS
#include "dma_led.h"
#else
#error                                                                         \
    "POVSegmented requires USE_DMA_LEDS (the Phantasm DMA LED transport): FastLED's bit-bang show() masks IRQs for windows that break the sync symbol margins (pitch > M, gap timeout > pitch + M; spec §5.2)."
#endif

#include "core/render/canvas.h"

#ifdef HS_PROFILE_ENABLE
namespace hs {
/** Column-ISR profiling accumulators (see IsrCycleStats): the whole flywheel
 *  wake, pack_column's pixel pack, and the submit_frame DMA marshal+kick.
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
 * junction. N=2 assigns one strip per arm only when the per-column DMA budget
 * permits it; at S=288 and 480 RPM, the supported counts are N=4 or N=8.
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
   * @brief Sync receive pin; shared with transmit on rev 1.1.
   */
  static constexpr int PIN_SYNC_RX = 3;
  /** @brief Sync transmit pin selected by the board revision. */
  static constexpr int PIN_SYNC_TX = HS_PHANTASM_BOARD_REV == 12 ? 4 : 3;

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
   * Wakes are advisory (position comes from the cycle counter); the grid only
   * quantizes when a column renders/emits, within the §5.2 self-censor budget.
   */
  static constexpr int OVERSAMPLE = 8;

  static_assert(COLUMN_US / float(OVERSAMPLE) >= 1.0f,
                "Flywheel wake period must be >= 1 us");

  /**
   * @brief NVIC priority for the sync-wire edge IRQ (Teensy 4 pin interrupts
   *        all share IRQ_GPIO6789).
   *
   * Must preempt the flywheel ISR so the ARM_DWT_CYCCNT stamp is the edge
   * time, not a service time. Cortex-M7 levels step by 16 (lower = higher
   * priority); Teensy defaults every IRQ to 128.
   */
  static constexpr uint8_t SYNC_EDGE_IRQ_PRIORITY = 16;

  /**
   * @brief HD107S SPI clock for the Phantasm DMA path, in Hz.
   */
  static constexpr uint32_t SPI_CLOCK_HZ = dma::SEGMENTED_CLOCK_HZ;

  /**
   * @brief Worst-case duration of one column's LED transfer, in µs.
   * @details Image frame plus the trailing strobe black frame, SPI data and
   * per-byte framing clocks, rounded up.
   */
  static constexpr unsigned long COLUMN_TRANSFER_US =
      dma::transfer_us(HD107SFrame<PPS>::COMPOSITE_SIZE, SPI_CLOCK_HZ);

  // An overrunning transfer surfaces as a dim image (dropped columns), not a
  // fault.
  static_assert(COLUMN_TRANSFER_US < COLUMN_US,
                "LED transfer outlasts the column period (S, N, RPM and "
                "SPI_CLOCK_HZ overrun the DMA every column)");

public:
  /**
   * @brief Foreground effect constructor: builds, arena-configures, and
   *  init()s one roster entry, ready to draw its first frame.
   * @details Must return a non-null, init()ed Effect; the result is
   *  dereferenced unguarded.
   */
  using EffectFactory = Effect *(*)();

  /**
   * @brief Drives MASTER_EN to its disabled level, parking the external
   *        sync-out buffer.
   * @details Call as the first statement of setup(). MASTER_EN is the '125
   *          channel-C/D /OE; until this runs, only the R_MEN pull-up (PCB rule
   *          R-LS-5) holds a board's sync driver off the shared bus.
   */
  HS_COLD_MEMBER static void park_sync_out() {
    digitalWriteFast(PIN_MASTER_EN, HIGH);
    pinMode(PIN_MASTER_EN, OUTPUT);
    if constexpr (HS_PHANTASM_BOARD_REV == 12) {
      digitalWriteFast(PIN_SYNC_TX, LOW);
      pinMode(PIN_SYNC_TX, OUTPUT);
    }
  }

  /**
   * @brief Hardware segment ID decoded from the GPIO straps.
   * @return Segment index in [0, N); 0 is the master.
   * @details Valid after begin() samples the straps with read_id().
   */
  static int segment_index() { return segment_id; }

  POVSegmented() = delete;
  /** Initializes the transport after Arduino core startup; call once. */
  HS_COLD_MEMBER static void begin() {
    read_id();
    configure_segment();

    ledController.begin();
    ledController.set_correction(hd107s::LINEAR_STRIP_GAIN.r,
                                 hd107s::LINEAR_STRIP_GAIN.g,
                                 hd107s::LINEAR_STRIP_GAIN.b);
    ledController.set_temperature(hd107s::LINEAR_WARM_GAIN.r,
                                  hd107s::LINEAR_WARM_GAIN.g,
                                  hd107s::LINEAR_WARM_GAIN.b);
    ledController.set_brightness(255);

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
   * The master broadcasts an EPOCH symbol when an effect's revolutions
   * elapse; every board then constructs the next roster entry during the
   * K-revolution commit window (display black) and swaps to its frame 0 at
   * the same boundary. Downstream boards join via the index beacon.
   *
   * @tparam R            Roster length, deduced from `factories`.
   * @param factories     One constructor per roster entry (HS_PHANTASM_EFFECT_LIST
   *                      order — identical on every board). Its length is the
   *                      roster length.
   * @param effect_revolutions Optional duration for each roster entry, in
   *                           revolutions.
   * @param stable_effect_seeds Optional RNG identity for each roster entry;
   *                           without it the per-visit stream is seeded from
   *                           the roster index alone.
   * @details Both optional tables must span the roster (arrays of R).
   */
  template <int R>
  [[noreturn]] static void
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
    if (effect_revolutions)
      cfg.set_effect_revolutions(*effect_revolutions);
#ifdef HS_PROFILE_EPOCH_REVS
    // Profiling: fixed epoch length in revolutions, replacing per-entry durations.
    cfg.revs_per_effect = HS_PROFILE_EPOCH_REVS;
    cfg.clear_effect_revolutions();
#endif
    const char *const bad_invariant = cfg.valid();
    HS_CHECK(bad_invariant == nullptr, "pov::sync::Config invariant: %s",
             bad_invariant);
    sync.configure(cfg);

    const bool master = (segment_id == 0);
    if (master || HS_PHANTASM_BOARD_REV == 12) {
      digitalWriteFast(PIN_SYNC_TX, LOW);
      pinMode(PIN_SYNC_TX, OUTPUT);
    }
    if (!master || HS_PHANTASM_BOARD_REV == 12) {
      pinMode(PIN_SYNC_RX, INPUT);
      // Pad hysteresis gives one interrupt per slow RC-filtered edge. pinMode
      // rewrites the pad-control register, so enable HYS afterward.
      *(portControlRegister(PIN_SYNC_RX)) |= IOMUXC_PAD_HYS;
    }
    // Enables the sync-bus driver; PIN_SYNC_TX must already be driven or a pad
    // keeper puts a spurious edge on the wire.
    digitalWriteFast(PIN_MASTER_EN, master ? LOW : HIGH);

    sync.seed(ARM_DWT_CYCCNT, master);

    if (!master) {
      attachInterrupt(digitalPinToInterrupt(PIN_SYNC_RX), sync_edge_isr,
                      RISING);
      // Above the flywheel's default 128 so the stamp is taken at the edge.
      NVIC_SET_PRIORITY(IRQ_GPIO6789, SYNC_EDGE_IRQ_PRIORITY);
    }
    HS_CHECK(timer.begin(flywheel_isr, COLUMN_US / float(OVERSAMPLE)),
             "flywheel IntervalTimer failed to start (no PIT channel)");

    // ── Foreground: construct effects on request, render, report ────────
    Effect *cur = nullptr;
    uint32_t built_gen = 0;
    uint32_t last_overrun = ledController.get_overrun_count();
    unsigned long last_report = millis();
    // The K-revolution construction budget, in the flywheel's own timebase.
    const uint32_t commit_budget_cycles =
        cfg.col_cycles(static_cast<int32_t>(cfg.commit_revs) * CANVAS_W);
    const uint32_t cycles_per_us = F_CPU_ACTUAL / 1000000u;
    uint32_t poll_prev_cycles = ARM_DWT_CYCCNT;

    for (;;) {
      const uint32_t poll_cycles = ARM_DWT_CYCCNT;
      const uint32_t bw = sync.build_word();
      const uint32_t gen = pov::sync::SyncBoard::build_gen_of(bw);
      if (gen != built_gen) {
        const uint32_t build_start_cycles = ARM_DWT_CYCCNT;
        const int32_t effect_index = pov::sync::SyncBoard::build_index_of(bw);
        const unsigned long t0 = micros();
        cur = pov::rebuild_effect(
            handoff, gen,
            [&] {
              HS_CHECK(micros() - t0 < 100000UL,
                       "flywheel ISR failed to release the live effect");
            },
            [&] { delete cur; },
            [&] {
              // Per-effect RNG seed from the synchronized effect index (spec §2).
              HS_CHECK(effect_index >= 0 && effect_index < R,
                       "sync published an out-of-roster effect index");
              hs::random().seed(
                  effect_seed_identities
                      ? effect_seed_identities[effect_index]
                      : hs::epoch_seed(static_cast<uint32_t>(effect_index)));
              cur = effect_factories[effect_index]();
              HS_CHECK(
                  cur->height() == ROWS,
                  "POVSegmented: effect canvas height must equal S/2 (ROWS)");
              HS_CHECK(cur->width() == CANVAS_W,
                       "POVSegmented: effect canvas width must equal CANVAS_W");
              HS_CHECK(
                  !cur->overrides_get_pixel(),
                  "POVSegmented: effect must not override get_pixel(); the "
                  "segmented pack_column path bypasses it");
              clip_to_segment(cur, /*arm_a_left=*/true);
              cur->set_clip_x(0, CANVAS_W);
              cur->draw_frame();
              cur->set_buffer_ready_hook(prepare_segment_clip);
              if (pov::segment_clip_applies(cur->needs_full_frame(),
                                            cur->persists_pixels()))
                cur->set_buffer_complete_hook(pov::preserve_segment_half);
              return cur;
            },
            [&](auto publish) {
              const uint32_t primask = hs::save_disable_interrupts();
              publish();
              hs::restore_interrupts(primask);
            });
        built_gen = gen;
        // `build` includes release, teardown, construction, and the first frame;
        // `window` starts at the previous poll, before the request arrived.
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

      // Health telemetry (spec §8.6): foreground-polled 1 Hz heartbeat.
      if (hs::debug && millis() - last_report >= 1000UL) {
        last_report = millis();
        const pov::sync::Telemetry tm = sync.telemetry_snapshot();
        // Each line fits hs::log's 256-byte buffer with saturated counters.
        hs::log("sync coast=%lu stall=%lu epi=%lu lock=%lu flip=%lu acc=%lu "
                "rej=%lu inv=%lu",
                (unsigned long)tm.max_coast_halves,
                (unsigned long)tm.master_stalls,
                (unsigned long)tm.epochs_refractory_ignored,
                (unsigned long)tm.lock_transitions, (unsigned long)tm.flips,
                (unsigned long)tm.symbols_accepted,
                (unsigned long)tm.symbols_rejected_gate,
                (unsigned long)tm.symbols_discarded_invalid);
        hs::log("sync emit cens=%lu abrt=%lu bdrop=%lu bbusy=%lu blate=%lu "
                "sdrop=%lu bok=%lu brej=%lu fix=%lu rmis=%lu",
                (unsigned long)tm.emit_censored, (unsigned long)tm.emit_aborted,
                (unsigned long)tm.beacons_overrun_dropped,
                (unsigned long)tm.beacons_busy_dropped,
                (unsigned long)tm.beacons_late_dropped,
                (unsigned long)tm.boundary_bursts_dropped,
                (unsigned long)tm.beacons_ok,
                (unsigned long)tm.beacons_rejected,
                (unsigned long)tm.beacon_index_corrections,
                (unsigned long)tm.beacon_rev_mismatches);
        const uint32_t overruns = ledController.get_overrun_count();
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
   * @details Grounded straps read LOW and decode_segment_id inverts them to set
   *          their ID bits; all-open straps decode to master ID 0. A duplicate
   *          master ID causes sync-bus contention. Duplicate peer IDs paint one
   *          segment twice and leave another dark. Assembly requires unique
   *          IDs and one master (R-ID-2/R-ID-4).
   */
  HS_COLD_MEMBER static void read_id() {
    pinMode(PIN_ID0, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 2)
      pinMode(PIN_ID1, INPUT_PULLUP);
    if constexpr (ID_STRAPS >= 3)
      pinMode(PIN_ID2, INPUT_PULLUP);
    delay(10); // settle time for pull-ups

    // Debounce: three samples ~5 ms apart must agree.
    const int raw0 = sample_strap();
    for (int i = 0; i < 2; ++i) {
      delay(5);
      HS_CHECK(sample_strap() == raw0,
               "unstable segment-ID strap (field/manufacturing fault)");
    }

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
   * @brief Clip @p e to this segment's display rectangle for the upcoming window.
   * @param e Effect to clip.
   * @param arm_a_left True if the window this frame displays in sweeps arm-A
   *        columns [0, CANVAS_W/2); arm B paints the opposite half.
   * @details No-op for an effect that needs_full_frame() or persists_pixels().
   */
  static void clip_to_segment(Effect *e, bool arm_a_left) {
    if (!pov::segment_clip_applies(e->needs_full_frame(), e->persists_pixels()))
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
   * no flywheel, flip, or epoch state (spec §8.2 single-writer model). Preempts
   * the flywheel ISR, so the mailbox claim must run with interrupts saved off.
   */
  static void sync_edge_isr() { sync.on_sync_edge(ARM_DWT_CYCCNT); }

  /**
   * @brief Flywheel ISR: the sole owner of all sync state (spec §8).
   *
   * Paced by an IntervalTimer at T0/OVERSAMPLE as a wake-up only; the cycle
   * counter decides the column (spec §4.1). Idempotent when the column is
   * unchanged, skip-tolerant when it jumped.
   */
  static void flywheel_isr() {
    HS_ISR_PROFILE(hs::g_flywheel_wake_cycles);
    const uint32_t now = ARM_DWT_CYCCNT;
    pov::run_wake_sequence(
        sync_pulse, submit_gate, handoff,
        [&] {
          pov::sync::BurstSnapshot burst;
          const pov::sync::BurstSnapshot *bp = nullptr;
          if (segment_id != 0 && sync.claim_sync_burst(now, &burst))
            bp = &burst;
          return sync.tick(now, bp);
        },
        [] { return pov::sync::SyncBoard::build_gen_of(sync.build_word()); },
        [](bool high) { digitalWriteFast(PIN_SYNC_TX, high ? HIGH : LOW); },
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
            return pack_black();
          case pov::SubmitAction::COLUMN:
            return pack_column(e, column);
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
   */
  [[nodiscard]] HS_O3_FN static bool pack_column(Effect *e, int x) {
    const int w = e->width();
    const int x_col = pov::segment_x_col(segment.arm_b, x, w);

    // Bypasses get_pixel(); effects are checked not to override it.
    const Pixel *buf = e->display_buffer();

    auto &frame = ledController.back_frame();
    {
      HS_ISR_PROFILE(hs::g_column_pack_cycles);
      const int stride = pov::segment_row_stride(segment, w);
      int off = pov::segment_pixel_base(segment, x_col, w);
      if (e->output_envelope_u16() == 65535u)
        for (int i = 0; i < PPS; ++i, off += stride)
          frame.pack_pixel(i, buf[off]);
      else
        for (int i = 0; i < PPS; ++i, off += stride)
          frame.pack_pixel(i, e->apply_output_envelope(buf[off]));
    }
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submit_frame(e->strobe_columns());
  }

  /**
   * @brief Re-submits the frame a previous wake had dropped on overrun.
   * @param strobe Whether the live effect wants the trailing black frame.
   * @return true if the transport accepted it this time.
   * @details No repack: submit_frame() returns before swapping buffers on an
   *          overrun, so back_frame() still holds the dropped column's pixels.
   */
  [[nodiscard]] static bool resubmit_frame(bool strobe) {
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submit_frame(strobe);
  }

  /**
   * @brief Submits one all-black frame (ACQUIRE / construction window).
   * @return true if the black frame was accepted by the LED transport; false
   *         if it was dropped on a DMA overrun (caller must retry, not latch).
   */
  [[nodiscard]] static bool pack_black() {
    auto &frame = ledController.back_frame();
    {
      HS_ISR_PROFILE(hs::g_column_pack_cycles);
      for (int i = 0; i < PPS; ++i) {
        frame.pack_pixel(i, Pixel(0, 0, 0));
      }
    }
    HS_ISR_PROFILE(hs::g_dma_submit_cycles);
    return ledController.submit_frame(false);
  }

  // ── Static state ────────────────────────────────────────────────────

  /**
   * @brief The sync engine: sole owner of all sync/flywheel state.
   * @details Written only by the flywheel ISR (tick()) and the edge ISR
   *          (mailbox publisher), spec §8. Foreground reads are single aligned
   *          words (build_word) or debug telemetry.
   */
  static pov::sync::SyncBoard sync;
  static IntervalTimer timer; /**< Flywheel wake-up timer (PIT channel).   */
  /** @brief Roster of effect constructors (HS_PHANTASM_EFFECT_LIST order). */
  static const EffectFactory *effect_factories;
  static const uint64_t *effect_seed_identities;

  /**
   * @brief Effect handoff state machine between the foreground and the ISR.
   * @details The foreground constructs and deletes; the ISR only dereferences
   *          the instance handed to it via live().
   */
  static pov::EffectHandoff<Effect> handoff;
  /**
   * @brief Per-wake LED-submit decision and its overrun-retry latches.
   */
  static pov::SubmitGate submit_gate;
  /**
   * @brief Sync-pulse width decision and its deferred-drop latch.
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

// DMAMEM survives only on an explicit specialization, so each instantiating
// target invokes HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM) once at file
// scope.
#define HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(S, N, RPM)                      \
  template <>                                                                  \
  DMAMEM DMALEDController<(S) / (N)> POVSegmented<S, N, RPM>::ledController {  \
    POVSegmented<S, N, RPM>::SPI_CLOCK_HZ                                      \
  }

#endif // ARDUINO
