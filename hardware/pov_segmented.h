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
 * which segment it owns.  Two shared wires connect all boards:
 *
 *   Wire 1 (column clock): Segment 0 generates PWM; all boards trigger
 *          their column ISR on each rising edge.
 *
 *   Wire 2 (frame sync):   Segment 0 pulses HIGH once per revolution
 *          (at x=0).  All boards reset their column counter on this
 *          edge, bounding any drift to ≤ 1 revolution.
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

#ifdef ARDUINO
#include <Arduino.h>

#ifdef USE_DMA_LEDS
  #include "dma_led.h"
#else
  #include <FastLED.h>
#endif

#include "canvas.h"
#include "memory.h"

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

  // ── Pin assignments ─────────────────────────────────────────────────

  /** @brief GPIO pins for hardware ID (active-low with internal pull-up). */
  static constexpr int PIN_ID0 = 21;
  static constexpr int PIN_ID1 = 22;

  /** @brief External interrupt input for shared column clock. */
  static constexpr int PIN_COLUMN_SYNC = 2;

  /** @brief Frame sync output (segment 0 pulses at x=0). */
  static constexpr int PIN_FRAME_SYNC_OUT = 3;

  /** @brief Frame sync input (all boards reset x on rising edge). */
  static constexpr int PIN_FRAME_SYNC_IN = 4;

  /** @brief PWM output for column clock (segment 0 / clock master only). */
  static constexpr int PIN_CLOCK_OUT = 5;

public:

  /**
   * @brief Initializes hardware: reads segment ID, configures LED driver,
   *        starts clock (master only), and demotes SysTick priority.
   */
  POVSegmented() {
    read_id();
    configure_segment();

#ifdef USE_DMA_LEDS
    ledController_.begin();
    ledController_.setCorrection(255, 176, 240);   // TypicalLEDStrip
    ledController_.setTemperature(255, 147, 41);    // Candle
    ledController_.setBrightness(255);
    // Slew-rate limiting for signal integrity
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_02 & ~IOMUXC_PAD_SRE;   // Pin 11
    IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 =
        IOMUXC_SW_PAD_CTL_PAD_GPIO_B0_03 & ~IOMUXC_PAD_SRE;   // Pin 13
#else
    FastLED.addLeds<WS2801, PIN_DATA, PIN_CLOCK, RGB, DATA_RATE_MHZ(6)>(
        leds_, PPS);
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
#endif

    // Demote SysTick so the column ISR is never preempted by millis().
    SCB_SHPR3 = 0x20200000;

    // Configure frame sync pins.
    if (segment_id_ == 0) {
      pinMode(PIN_FRAME_SYNC_OUT, OUTPUT);
      digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
    }
    pinMode(PIN_FRAME_SYNC_IN, INPUT);

    start_clock();

    Serial.print("[Phantasm] Segment ");
    Serial.print(segment_id_);
    Serial.print(arm_b_ ? " arm-B" : " arm-A");
    Serial.print(y_step_ < 0 ? " (bottom, reversed)" : " (top)");
    Serial.print(" | y_base=");
    Serial.print(y_base_);
    Serial.print(" y_step=");
    Serial.print(y_step_);
    Serial.print(" pixels=");
    Serial.println(PPS);
  }

  /**
   * @brief Runs a specific effect for a given duration.
   * @tparam E The Effect class to instantiate and run.
   * @param duration Time in seconds to display the effect.
   */
  template <typename E>
  void show(unsigned long duration) {
    E *e = new E();
    e->init();
    run(e, duration);
  }

private:

  // ── Hardware ID ─────────────────────────────────────────────────────

  /**
   * @brief Reads the 2-bit hardware segment ID from GPIO pins.
   *
   * Pins are configured as INPUT_PULLUP; grounding a pin sets its bit.
   * The result is inverted so that all-floating = ID 0 (clock master).
   */
  void read_id() {
    pinMode(PIN_ID0, INPUT_PULLUP);
    pinMode(PIN_ID1, INPUT_PULLUP);
    delay(10);  // settle time for pull-ups

    segment_id_ = ~(digitalReadFast(PIN_ID0)
                   | (digitalReadFast(PIN_ID1) << 1)) & (N - 1);
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
   *                         y_base = ROWS-1, y_step = -1  →  y ∈ [ROWS-1, PPS]
   *                         (strip physically reversed)
   */
  void configure_segment() {
    arm_b_   = (segment_id_ >= SEGS_PER_ARM);
    int arm_seg = segment_id_ % SEGS_PER_ARM;

    if (arm_seg == 0) {
      // Top strip: LED 0 at N pole, counting toward junction
      y_base_ = 0;
      y_step_ = 1;
    } else {
      // Bottom strip: LED 0 at S pole, counting toward junction (reversed)
      y_base_ = ROWS - 1;
      y_step_ = -1;
    }
  }

  // ── Clock ───────────────────────────────────────────────────────────

  /**
   * @brief Starts the column clock (segment 0 only).
   *
   * Generates a PWM signal at (RPM / 60) × CANVAS_W Hz on the clock
   * output pin.  All segments trigger their column ISR from this shared
   * clock line via external interrupt.
   */
  void start_clock() {
    if (segment_id_ == 0) {
      float col_freq = (RPM / 60.0f) * CANVAS_W;
      analogWriteFrequency(PIN_CLOCK_OUT, col_freq);
      analogWrite(PIN_CLOCK_OUT, 128);
    }
  }

  // ── Effect lifecycle ────────────────────────────────────────────────

  /**
   * @brief Non-template core of show(). Takes ownership of the effect.
   * @param e        Effect instance (deleted on exit).
   * @param duration Display time in seconds.
   */
  void run(Effect *e, unsigned long duration) {
    configure_arenas_default();

    long start = millis();
    effect_ = e;
    x_ = 0;

    // Attach column + frame sync ISRs.
    pinMode(PIN_COLUMN_SYNC, INPUT);
    attachInterrupt(digitalPinToInterrupt(PIN_COLUMN_SYNC),
                    show_col, RISING);
    attachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC_IN),
                    frame_sync_isr, RISING);

    while (millis() - start < duration * 1000) {
      e->draw_frame();
    }

    detachInterrupt(digitalPinToInterrupt(PIN_COLUMN_SYNC));
    detachInterrupt(digitalPinToInterrupt(PIN_FRAME_SYNC_IN));
    delete effect_;
    effect_ = nullptr;
  }

  // ── Column ISR ──────────────────────────────────────────────────────

  /**
   * @brief ISR called on every rising edge of the shared column clock.
   *
   * Packs this segment's pixels into the DMA frame using the precomputed
   * y_base_ / y_step_ mapping.  Arm B segments read from x + W/2
   * (opposite half of the image).
   *
   * The loop is branchless — all per-segment decisions are resolved at
   * boot time in configure_segment().
   *
   * When the column counter wraps to 0 or W/2, the display buffer is
   * advanced (double-buffer flip for the ISR read pointer).
   */
  static FASTRUN void show_col() {
    const int w = effect_->width();
    const int x_col = arm_b_ ? (x_ + w / 2) % w : x_;

#if defined(USE_DMA_LEDS)
    auto &frame = ledController_.backFrame();

    int y = y_base_;
    for (int i = 0; i < PPS; ++i, y += y_step_) {
      frame.packPixel(i, effect_->get_pixel(x_col, y));
    }
    ledController_.submitFrame(effect_->show_bg());

#else
    int y = y_base_;
    for (int i = 0; i < PPS; ++i, y += y_step_) {
      leds_[i] = effect_->get_pixel(x_col, y);
    }
    FastLED.show();
    if (effect_->show_bg()) {
      FastLED.showColor(CRGB(0, 0, 0));
    }
#endif

    x_ = (x_ + 1) % w;
    if (x_ == 0) {
      effect_->advance_display();
      // Segment 0: pulse the frame sync line so all boards re-align.
      if (segment_id_ == 0) {
        digitalWriteFast(PIN_FRAME_SYNC_OUT, HIGH);
        // Pulse will be cleared by the frame_sync_isr on this board too.
      }
    } else if (x_ == w / 2) {
      effect_->advance_display();
    }
  }

  /**
   * @brief ISR for the frame sync line.  Resets the column counter to 0,
   *        re-aligning this board if any column pulses were missed.
   *
   * Fired by segment 0's pulse at the start of each revolution.  On
   * segment 0 itself this also clears the sync output.
   */
  static FASTRUN void frame_sync_isr() {
    x_ = 0;
    if (segment_id_ == 0) {
      digitalWriteFast(PIN_FRAME_SYNC_OUT, LOW);
    }
  }

  // ── Static state ────────────────────────────────────────────────────

#ifndef USE_DMA_LEDS
  static CRGB leds_[PPS];
#endif
  static Effect *effect_;
  static int x_;
  static int segment_id_;
  static bool arm_b_;
  static int y_base_;
  static int y_step_;
#if defined(USE_DMA_LEDS)
  static DMALEDController<PPS> ledController_;
#endif
};

// ── Static member definitions ───────────────────────────────────────────

template <int S, int N, int RPM>
int POVSegmented<S, N, RPM>::x_ = 0;

template <int S, int N, int RPM>
Effect *POVSegmented<S, N, RPM>::effect_ = nullptr;

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
template <int S, int N, int RPM>
DMALEDController<POVSegmented<S, N, RPM>::PPS>
    POVSegmented<S, N, RPM>::ledController_;
#endif

#endif // ARDUINO
