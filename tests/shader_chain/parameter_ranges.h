/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Finds a numeric range and checks the emitted interval contract. */
inline std::array<float, 2> parameter_range(ShaderChainWhiteBox::FX &effect,
                                            const char *name) {
  std::array<float, 2> result{NAN, NAN};
  effect.parameter_ranges(
      [&](const char *candidate, float minimum, float maximum) {
        HS_CONTEXT(candidate);
        if (std::strcmp(candidate, name) == 0)
          result = {minimum, maximum};
        const auto *parameter = effect.getParameters().find(candidate);
        HS_EXPECT_TRUE(parameter != nullptr);
        if (parameter) {
          HS_EXPECT_LE(parameter->min, minimum);
          HS_EXPECT_GE(parameter->max, maximum);
          HS_EXPECT_LE(minimum, parameter->get_requested());
          HS_EXPECT_GE(maximum, parameter->get_requested());
        }
      });
  return result;
}

/** @brief Slider limits follow whole-chain growth without changing live state. */
inline void test_shader_chain_parameter_ranges() {
  AffineWarpChain fixture;
  constexpr float MIN_SCALE = AffineWarpChain::MIN_SCALE;
  const auto &ids = fixture.ids;
  const std::string refused_scale_x =
      ids[AffineWarpChain::ADMITTED] + ".scale-x";
  const auto &chain = fixture.chain;
  reset_globals();
  ShaderChainWhiteBox::FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  for (int index = 0; index < AffineWarpChain::ADMITTED; ++index) {
    const std::string x = ids[index] + ".scale-x";
    const std::string y = ids[index] + ".scale-y";
    const ShaderChainParameterWrite writes[] = {{x.c_str(), MIN_SCALE},
                                                {y.c_str(), MIN_SCALE}};
    HS_EXPECT_EQ(effect.update_parameters(writes), ParamSetResult::APPLIED);
  }
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(), MIN_SCALE),
               ParamSetResult::INADMISSIBLE);
  effect.setAnimationsPaused(false);
  const auto generation = effect.getParameterSchemaGeneration();
  const char *warning = effect.parameter_warning(refused_scale_x.c_str());
  auto &program = ShaderChainWhiteBox::program(effect);
  std::vector<std::vector<uint8_t>> params;
  std::vector<std::vector<uint8_t>> states;
  for (size_t index = 0; index < program.ops().size(); ++index) {
    const auto &runtime = program.ops()[index].op->runtime;
    const auto *p = static_cast<const uint8_t *>(program.param_block(index));
    const auto *s = static_cast<const uint8_t *>(program.state_block(index));
    params.emplace_back(p, p + runtime.param.size);
    states.emplace_back(s, s + runtime.state.size);
  }
  effect.parameter_ranges([&](const char *name, float, float) {
    const auto *parameter = effect.getParameters().find(name);
    HS_EXPECT_TRUE(parameter != nullptr);
    if (parameter)
      HS_EXPECT_FALSE(parameter->is_enum());
  });
  const auto limited = parameter_range(effect, refused_scale_x.c_str());
  HS_EXPECT_GT(limited[0], MIN_SCALE);
  HS_EXPECT_LE(limited[0], 1.0f);
  HS_EXPECT_GT(limited[1], 1.0f);
  HS_EXPECT_LT(limited[1], 64.0f);
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), generation);
  HS_EXPECT_TRUE(effect.parameter_warning(refused_scale_x.c_str()) == warning);
  for (size_t index = 0; index < program.ops().size(); ++index) {
    HS_EXPECT_EQ(std::memcmp(params[index].data(), program.param_block(index),
                             params[index].size()),
                 0);
    HS_EXPECT_EQ(std::memcmp(states[index].data(), program.state_block(index),
                             states[index].size()),
                 0);
  }
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(), limited[0]),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(),
                                      std::nextafter(limited[0], 0.0f)),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(), limited[1]),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(),
                                      std::nextafter(limited[1], INFINITY)),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(), 1.0f),
               ParamSetResult::APPLIED);
  const ShaderChainParameterWrite release[] = {{"warp0.scale-x", 1.0f},
                                               {"warp0.scale-y", 1.0f}};
  HS_EXPECT_EQ(effect.update_parameters(release), ParamSetResult::APPLIED);
  const auto expanded = parameter_range(effect, refused_scale_x.c_str());
  HS_EXPECT_LT(expanded[0], limited[0]);
  HS_EXPECT_EQ(effect.updateParameter(refused_scale_x.c_str(), expanded[0]),
               ParamSetResult::APPLIED);
}

/** @brief Axis sliders stay in a normalizable connected interval. */
inline void test_shader_chain_axis_ranges() {
  reset_globals();
  ShaderChainWhiteBox::FX effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"project", In::Op::ProjectStereographic::ID},
      {"sample", In::Op::SampleGridV4::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID},
      {"modulation", In::Op::ModulateColorV1::ID},
  };
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  const ShaderChainParameterWrite direction[] = {
      {"modulation.modulation-source", 1.0f},
      {"modulation.modulation-evaluation", 0.0f},
      {"modulation.sine-axis-x", 1.0f},
      {"modulation.sine-axis-y", 0.0f},
      {"modulation.sine-axis-z", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(direction), ParamSetResult::APPLIED);
  const ParamDef *axis = effect.getParameters().find("modulation.sine-axis-x");
  HS_EXPECT_TRUE(axis != nullptr);
  if (!axis)
    return;
  const float AXIS_MIN = axis->min, AXIS_MAX = axis->max;
  const auto positive = parameter_range(effect, "modulation.sine-axis-x");
  HS_EXPECT_GT(positive[0], 0.0f);
  HS_EXPECT_EQ(positive[1], AXIS_MAX);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x", positive[0]),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x",
                                      std::nextafter(positive[0], 0.0f)),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x", -1.0f),
               ParamSetResult::APPLIED);
  const auto negative = parameter_range(effect, "modulation.sine-axis-x");
  HS_EXPECT_EQ(negative[0], AXIS_MIN);
  HS_EXPECT_LT(negative[1], 0.0f);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x", negative[1]),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x",
                                      std::nextafter(negative[1], 0.0f)),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-y", 1.0f),
               ParamSetResult::APPLIED);
  const auto crossing = parameter_range(effect, "modulation.sine-axis-x");
  HS_EXPECT_EQ(crossing[0], AXIS_MIN);
  HS_EXPECT_EQ(crossing[1], AXIS_MAX);
  HS_EXPECT_EQ(effect.updateParameter("modulation.sine-axis-x", 0.0f),
               ParamSetResult::APPLIED);
}

/** @brief Harmonic choices follow the active term and its paired degree/order. */
inline void test_shader_chain_harmonic_ranges() {
  reset_globals();
  ShaderChainWhiteBox::FX effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"field", In::Op::DisplacePotential::ID},
      {"project", In::Op::ProjectStereographic::ID},
      {"sample", In::Op::SampleGridV4::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  const ShaderChainParameterWrite mode[] = {{"field.term-1-l", 2},
                                            {"field.term-1-m", 1}};
  HS_EXPECT_EQ(effect.update_parameters(mode), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-l")[0], 1.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-l")[1], 8.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-m")[0], 0.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-m")[1], 2.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-m", 3),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-l", 0),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-type", 1),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-l")[0], 0.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-l")[1], 8.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-m")[0], 0.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-m")[1], 8.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-m", 7),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-type")[0], 1.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-type")[1], 3.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-type", 0),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(requested_param(effect, "field.term-1-m"), 7.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-l", 8),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-type")[0], 0.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-type")[1], 3.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-type", 0),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-l")[0], 7.0f);
  HS_EXPECT_EQ(parameter_range(effect, "field.term-1-m")[1], 8.0f);
}
