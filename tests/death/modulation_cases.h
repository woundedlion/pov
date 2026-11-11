/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Modulation configuration and prepared-resource admission traps.

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
  const auto state = Pullback::Modulation::init(config);
  Pullback::Modulation::Resources<> resources;
  (void)resources.prepare(config, state);
  config.version = opaque(uint16_t{0});
  (void)resources.prepare(config, state);
}
inline void case_modulation_missing_cube_capacity() {
  ColorModulation::Config config;
  config.mechanic.amount = 0.5f;
  const auto state = Pullback::Modulation::init(config);
  Pullback::Modulation::Resources<false> resources;
  (void)resources.prepare(config, state);
}

inline void case_modulation_retain_table_size() {
  static int8_t cube[Pattern::SphereLutView::SIZE] = {};
  static int8_t table[Pattern::SphereLutView::SIZE] = {};
  static Pullback::Modulation::Prepared source;
  static Pullback::Modulation::Prepared destination;
  source.cube = {cube, true};
  (void)Pullback::Modulation::retain(&source, destination, table);
  (void)Pullback::Modulation::retain(
      &source, destination, std::span<int8_t>(table, opaque(size_t{1})));
}
inline void case_modulation_retain_linear_cube() {
  static int8_t cube[Pattern::SphereLutView::SIZE] = {};
  static int8_t table[Pattern::SphereLutView::SIZE] = {};
  static float linear[1] = {};
  static Pullback::Modulation::Prepared source;
  static Pullback::Modulation::Prepared destination;
  source.cube = {cube, true};
  (void)Pullback::Modulation::retain(&source, destination, table);
  source.cube.linear_data = opaque(static_cast<const float *>(linear));
  (void)Pullback::Modulation::retain(&source, destination, table);
}

inline void case_modulation_unknown_kernel() {
  (void)ColorModulation::kernel(
      static_cast<ColorModulation::Mechanic>(opaque(uint8_t{255})));
}
inline void case_modulation_control_index_invalid() {
  Pullback::Interp::Op::ModulationParams params;
  (void)params.extra_address(
      opaque(Pullback::Interp::Op::ModulationParams::EXTRA_SCHEMA.size()));
}

/** @brief Active spatial modulation requires its instance's prepared endpoint. */
inline void case_modulation_missing_endpoint() {
  using Op = Pullback::Interp::Op::ModulateColorV1;
  Op::Params params;
  params.modulation.evaluation = Pattern::Evaluation::DIRECT;
  params.modulation.mechanic.amount = 0.5f;
  Op::State state;
  state.key = "layer";
  Pullback::Modulation::Resources<false> resources;
  const auto &source = resources.prepare(params.modulation, {});
  const Pullback::Modulation::Endpoint endpoint{state.key, &source};
  Pullback::Interp::FrameContext context;
  context.modulations = std::span(&endpoint, 1);
  (void)Op::prepare(context, params, state);
  context.modulations = {};
  (void)Op::prepare(context, params, state);
}
