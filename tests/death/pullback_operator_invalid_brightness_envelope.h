/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/** @brief Death case: the generated-palette operator rejects an unknown
    brightness envelope. */
inline void case_pullback_operator_invalid_brightness_envelope() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.envelope_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}
