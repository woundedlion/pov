/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the double-buffered DMA LED controller
 * (hardware/dma_led_controller.h). The controller is templated on its SPI/DMA
 * transport, so these run it against MockStrip — a recording double standing in
 * for the Arduino-only TeensySPIDMA — to cover the orchestration the pure-math
 * (test_dma_core.h) and wire-format (test_hd107s_frame.h) suites cannot reach:
 * the double-buffer flip, the overrun-drop path (and its watchdog consult), the
 * with_bg transfer-length select, and the end-to-end byte stream a real strip
 * would clock in. The register/eDMA/ISR internals of TeensySPIDMA stay
 * hardware-only and out of host-test scope (see the wedged-channel death case in
 * test_death.h for the overrun-drop watchdog trap).
 */
#pragma once

#include "hardware/dma_led_controller.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstddef>
#include <cstdint>

namespace hs_test {
namespace dma_controller {

constexpr int N = 8; // BUFFER_SIZE=37, COMPOSITE_SIZE=74

using Frame = HD107SFrame<N>;

/**
 * @brief Recording transport double for DMALEDController.
 * @details Mirrors TeensySPIDMA's transmit/completion/watchdog contract while
 *          capturing each transfer's pointer, length, and byte snapshot and
 *          letting a test drive completion. The controller owns its transport by
 *          value and constructs it from a clock alone, so the observable state
 *          lives in a single static block (the device likewise holds exactly one
 *          transport per image); reset() clears it between tests.
 */
class MockStrip {
public:
  static constexpr int CAPTURE_CAP = 128; // >= COMPOSITE_SIZE for this N

  /**
   * @brief State observed and driven by the test.
   */
  struct State {
    uint32_t clock = 0;      /**< Clock the controller forwarded at construction. */
    bool complete = true;    /**< Completion flag; a fresh channel is idle. */
    bool wedged = false;     /**< When set, checkStaleTransfer() traps. */
    int init_calls = 0;      /**< init() invocations. */
    int transmit_calls = 0;  /**< transmitAsync() invocations. */
    int check_stale_calls = 0; /**< checkStaleTransfer() invocations. */
    const uint8_t *last_data = nullptr; /**< Pointer handed to the last transmit. */
    size_t last_len = 0;     /**< Length handed to the last transmit. */
    uint8_t capture[CAPTURE_CAP] = {}; /**< Snapshot of the last transmit's bytes. */
  };

  /**
   * @brief Returns the shared observable state.
   * @return Reference to the single static State block.
   */
  static State &state() {
    static State s;
    return s;
  }

  /**
   * @brief Clears the shared state so a test starts from an idle channel.
   */
  static void reset() { state() = State{}; }

  /**
   * @brief Records the forwarded clock; matches the transport ctor contract.
   * @param clock SPI clock in Hz forwarded by the controller.
   */
  explicit MockStrip(uint32_t clock = 12000000) { state().clock = clock; }

  /**
   * @brief Counts a hardware-init call.
   */
  void init() { ++state().init_calls; }

  /**
   * @brief Reports whether the in-flight transfer has finished.
   * @return The test-driven completion flag.
   */
  bool isComplete() const { return state().complete; }

  /**
   * @brief Watchdog consult on the overrun-drop path; traps when wedged.
   */
  void checkStaleTransfer() {
    ++state().check_stale_calls;
    HS_CHECK(!state().wedged,
             "MockStrip: wedged channel — watchdog trap on the overrun path");
  }

  /**
   * @brief Records a transfer and snapshots its bytes.
   * @param data Pointer to the frame buffer being transmitted.
   * @param len Number of bytes in the transfer.
   * @details Marks the channel in-flight; the test completes it via complete().
   */
  void transmitAsync(const uint8_t *data, size_t len) {
    HS_CHECK(state().complete,
             "MockStrip: transmitAsync while a transfer is still in flight");
    State &s = state();
    s.last_data = data;
    s.last_len = len;
    ++s.transmit_calls;
    size_t n = len < static_cast<size_t>(CAPTURE_CAP) ? len : CAPTURE_CAP;
    for (size_t i = 0; i < n; ++i)
      s.capture[i] = data[i];
    s.complete = false;
  }
};

/**
 * @brief Restores HD107SFrame's shared static correction state to unity.
 */
inline void reset_correction() {
  Frame::setCorrection(255, 255, 255);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
}

/**
 * @brief begin() forwards to the transport's one-time init().
 */
inline void test_begin_inits() {
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;
  HS_EXPECT_EQ(MockStrip::state().init_calls, 0);
  ctl.begin();
  HS_EXPECT_EQ(MockStrip::state().init_calls, 1);
}

/**
 * @brief A clean submit hands the image-frame bytes to the transport and bumps
 * the transfer counter.
 */
inline void test_submit_happy_path() {
  reset_correction();
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;

  ctl.backFrame().packPixel(0, Pixel16(CRGB(255, 0, 0)));
  bool ok = ctl.submitFrame(/*withBg=*/false);

  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(ctl.getTransferCount(), 1u);
  HS_EXPECT_EQ(ctl.getOverrunCount(), 0u);
  HS_EXPECT_EQ(MockStrip::state().transmit_calls, 1);
  HS_EXPECT_EQ(MockStrip::state().last_len, static_cast<size_t>(Frame::BUFFER_SIZE));
  HS_EXPECT_TRUE(MockStrip::state().last_data != nullptr);
}

/**
 * @brief Successive submits alternate buffers, and the buffer just handed to DMA
 * is never the one the next backFrame() exposes for writing.
 */
inline void test_double_buffer_flip() {
  reset_correction();
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;

  const uint8_t *write0 =
      reinterpret_cast<const uint8_t *>(ctl.backFrame().data());
  HS_EXPECT_TRUE(ctl.submitFrame(false));
  const uint8_t *dma1 = MockStrip::state().last_data;
  HS_EXPECT_EQ(dma1, write0); // the frame just written is the one DMA'd

  // Front buffer is in flight; the exposed back buffer must be the other one.
  const uint8_t *write1 =
      reinterpret_cast<const uint8_t *>(ctl.backFrame().data());
  HS_EXPECT_NE(write1, dma1);

  MockStrip::state().complete = true; // signal the prior transfer done
  HS_EXPECT_TRUE(ctl.submitFrame(false));
  const uint8_t *dma2 = MockStrip::state().last_data;
  HS_EXPECT_EQ(dma2, write1);
  HS_EXPECT_NE(dma2, dma1); // buffers alternate

  MockStrip::state().complete = true;
  HS_EXPECT_TRUE(ctl.submitFrame(false));
  HS_EXPECT_EQ(MockStrip::state().last_data, dma1); // back to the first buffer
  HS_EXPECT_EQ(ctl.getTransferCount(), 3u);
}

/**
 * @brief An overrun drops the frame: no transmit, counters reflect the drop, the
 * watchdog is consulted, and the active buffer does not flip.
 */
inline void test_overrun_drop() {
  reset_correction();
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;

  HS_EXPECT_TRUE(ctl.submitFrame(false)); // first transfer, now in flight
  const uint8_t *back_before =
      reinterpret_cast<const uint8_t *>(ctl.backFrame().data());

  bool ok = ctl.submitFrame(false); // prior still in flight -> overrun
  HS_EXPECT_FALSE(ok);
  HS_EXPECT_EQ(ctl.getOverrunCount(), 1u);
  HS_EXPECT_EQ(ctl.getTransferCount(), 1u);      // unchanged
  HS_EXPECT_EQ(MockStrip::state().transmit_calls, 1); // no new transmit
  HS_EXPECT_EQ(MockStrip::state().check_stale_calls, 1); // watchdog consulted

  const uint8_t *back_after =
      reinterpret_cast<const uint8_t *>(ctl.backFrame().data());
  HS_EXPECT_EQ(back_after, back_before); // no buffer flip on a drop
}

/**
 * @brief withBg selects the composite (image + trailing black) transfer length.
 */
inline void test_withbg_length() {
  reset_correction();
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;

  HS_EXPECT_TRUE(ctl.submitFrame(/*withBg=*/true));
  HS_EXPECT_EQ(MockStrip::state().last_len,
               static_cast<size_t>(Frame::COMPOSITE_SIZE));

  MockStrip::state().complete = true;
  HS_EXPECT_TRUE(ctl.submitFrame(/*withBg=*/false));
  HS_EXPECT_EQ(MockStrip::state().last_len,
               static_cast<size_t>(Frame::BUFFER_SIZE));
}

/**
 * @brief The bytes the controller hands to DMA are exactly the HD107S wire image
 * an independently packed frame produces — pixels → correction → wire, verified
 * end-to-end through the controller.
 */
inline void test_end_to_end_wire_bytes() {
  reset_correction();
  MockStrip::reset();
  DMALEDController<N, MockStrip> ctl;

  const CRGB colors[N] = {CRGB(255, 0, 0),     CRGB(0, 255, 0),
                          CRGB(0, 0, 255),     CRGB(255, 255, 255),
                          CRGB(200, 100, 50),  CRGB(0, 0, 0),
                          CRGB(12, 240, 60),   CRGB(90, 30, 210)};

  Frame ref;
  for (int i = 0; i < N; ++i)
    ref.packPixel(i, Pixel16(colors[i]));

  for (int i = 0; i < N; ++i)
    ctl.backFrame().packPixel(i, Pixel16(colors[i]));
  HS_EXPECT_TRUE(ctl.submitFrame(false));

  HS_EXPECT_EQ(MockStrip::state().last_len, static_cast<size_t>(Frame::BUFFER_SIZE));
  for (int k = 0; k < Frame::BUFFER_SIZE; ++k)
    HS_EXPECT_EQ(MockStrip::state().capture[k], ref.data()[k]);
}

/**
 * @brief Runs the DMA controller suite.
 * @return The module's failure count from end_module().
 */
inline int run_dma_controller_tests() {
  hs_test::ModuleFixture fixture("dma_controller");
  test_begin_inits();
  test_submit_happy_path();
  test_double_buffer_flip();
  test_overrun_drop();
  test_withbg_length();
  test_end_to_end_wire_bytes();
  reset_correction(); // leave shared static state clean for any later module
  return fixture.result();
}

} // namespace dma_controller
} // namespace hs_test
