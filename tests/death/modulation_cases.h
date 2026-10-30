/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file modulation_cases.h
 * @brief Modulation configuration and prepared-resource admission traps. */

inline void case_pattern_direction_invalid() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(0.0f), 0.0f, 0.0f));
}
inline void case_modulation_hue_nonfinite() {
  (void)ColorModulation::HueRotation::prepare({0.5f});
  (void)ColorModulation::HueRotation::prepare(
      {opaque(std::numeric_limits<float>::infinity())});
}
inline void case_modulation_unit_amount_invalid() {
  (void)ColorModulation::Opacity::prepare({0.5f});
  (void)ColorModulation::Opacity::prepare({opaque(2.0f)});
}
inline void case_modulation_unknown_prepare() {
  ColorModulation::Params params;
  ColorModulation::validate(params);
  params.kind = static_cast<ColorModulation::Mechanic>(opaque(uint8_t{255}));
  ColorModulation::validate(params);
}
inline void case_modulation_invalid_initialization() {
  ColorModulation::Config config;
  (void)Pullback::Modulation::init(config);
  config.enabled = false;
  config.version = opaque(uint16_t{0});
  (void)Pullback::Modulation::init(config);
}
inline void case_modulation_invalid_synchronization() {
  ColorModulation::Config config;
  auto state = Pullback::Modulation::init(config);
  Pullback::Modulation::synchronize(state, config);
  config.version = opaque(uint16_t{0});
  Pullback::Modulation::synchronize(state, config);
}
inline void case_modulation_invalid_preparation() {
  ColorModulation::Config config;
  const auto STATE = Pullback::Modulation::init(config);
  Pullback::Modulation::Resources<> resources;
  (void)resources.prepare(config, STATE);
  config.version = opaque(uint16_t{0});
  (void)resources.prepare(config, STATE);
}
inline void case_modulation_missing_cube_capacity() {
  ColorModulation::Config config;
  config.mechanic.amount = 0.5f;
  const auto STATE = Pullback::Modulation::init(config);
  Pullback::Modulation::Resources<false> resources;
  (void)resources.prepare(config, STATE);
}

inline void case_pattern_cube_invalid_direction() {
  (void)Pattern::prepare_sphere_lut_sample(
      math::Vector(opaque(0.0f), 0.0f, 0.0f));
}
inline void case_modulation_unknown_kernel() {
  (void)ColorModulation::kernel(
      static_cast<ColorModulation::Mechanic>(opaque(uint8_t{255})));
}
inline void case_modulation_control_index_invalid() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  (void)params.extra_address(opaque(size_t{38}));
}
