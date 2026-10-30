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

/** @brief Field color stages require retained field metadata in both paths. */
inline void reject_field_modulation_missing_metadata(bool temporal) {
  using namespace Pullback::Interp;
  Op::ModulateField::Params params;
  params.modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  params.modulation.evaluation = Pattern::Evaluation::DIRECT;
  params.modulation.sine.frequency = 0.0f;
  params.modulation.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  Op::ModulateField::State state;
  Op::ModulateField::init(state, {"field", Op::ModulateField::ID, 1337});
  Pullback::FieldSample sample{.5f, 1.0f, math::Y_AXIS, 1.0f};
  FrameContext context;
  context.color_sample = &sample;
  const auto prepared = Op::ModulateField::prepare(context, params, state);
  const Color4 INPUT(Pixel(1000, 2000, 3000), .8f);
  const auto run = [&] {
    return temporal ? Op::ModulateField::run_temporal(INPUT, context, prepared,
                                                      prepared, .5f)
                    : Op::ModulateField::run(INPUT, context, params, prepared);
  };
  const Color4 valid = run();
  HS_CHECK(valid.color == INPUT.color && valid.alpha == .4f,
           "valid field metadata companion failed");
  context.color_sample = opaque<const Pullback::FieldSample *>(nullptr);
  (void)run();
}

inline void case_modulation_field_missing_metadata() {
  reject_field_modulation_missing_metadata(false);
}
inline void case_modulation_temporal_field_missing_metadata() {
  reject_field_modulation_missing_metadata(true);
}
inline void case_modulation_hue_invalid_mode() {
  using namespace Pullback::Interp;
  Op::ModulateHueShift::Params params;
  Op::ModulateHueShift::State state;
  FrameContext context;
  (void)Op::ModulateHueShift::prepare(context, params, state);
  params.hue_mode = opaque(uint8_t{255});
  (void)Op::ModulateHueShift::prepare(context, params, state);
}
inline void case_modulation_hue_missing_metadata() {
  using namespace Pullback::Interp;
  Op::ModulateHueShift::Params params;
  params.hue_shift_amount = 1.0f;
  Op::ModulateHueShift::State state;
  Pullback::FieldSample sample{.5f, 1.0f, math::Y_AXIS, 1.0f};
  FrameContext context;
  context.color_sample = &sample;
  const auto prepared = Op::ModulateHueShift::prepare(context, params, state);
  const Color4 INPUT(Pixel(1000, 2000, 3000), .8f);
  const Color4 valid =
      Op::ModulateHueShift::run(INPUT, context, params, prepared);
  HS_CHECK(valid.color == INPUT.color && valid.alpha == INPUT.alpha,
           "valid hue metadata companion failed");
  context.color_sample = opaque<const Pullback::FieldSample *>(nullptr);
  (void)Op::ModulateHueShift::run(INPUT, context, params, prepared);
}
