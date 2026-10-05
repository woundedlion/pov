/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_effects.h.

/** @brief Chooses a distant writable value for a parameter lint. */
inline float parameter_probe_target(const ParamDef &def, float current) {
  if (def.option_values != nullptr) {
    float target = static_cast<float>(def.option_values[0]);
    for (int i = 1; i < def.option_count; ++i) {
      const float candidate = static_cast<float>(def.option_values[i]);
      if (fabsf(candidate - current) > fabsf(target - current))
        target = candidate;
    }
    return target;
  }
  return def.is_bool() ? (current > 0.5f ? 0.0f : 1.0f)
         : (current - def.min) > (def.max - current)
             ? def.min + 0.25f * (def.max - def.min)
             : def.min + 0.75f * (def.max - def.min);
}

inline void test_parameter_probe_targets() {
  static constexpr const int64_t IDS[] = {0, 1, 6};
  ParamDef sparse;
  sparse.min = 0;
  sparse.max = 6;
  sparse.option_count = std::size(IDS);
  sparse.option_values = IDS;
  sparse.target_type = ParamDef::TargetType::INT_U8;
  for (float current : {0.0f, 1.0f, 6.0f}) {
    float target = parameter_probe_target(sparse, current);
    HS_EXPECT_EQ(target, current == 6.0f ? 0.0f : 6.0f);
    HS_EXPECT_NE(target, current);
    HS_EXPECT_EQ(sparse.normalize(target), ParamSetResult::APPLIED);
  }
  ParamDef dense = sparse;
  dense.max = 2;
  dense.option_values = nullptr;
  HS_EXPECT_EQ(parameter_probe_target(dense, 0.0f), 1.5f);
  HS_EXPECT_EQ(parameter_probe_target(dense, 2.0f), 0.5f);
  ParamDef boolean;
  boolean.target_type = ParamDef::TargetType::BOOL;
  HS_EXPECT_EQ(parameter_probe_target(boolean, 0.0f), 1.0f);
  HS_EXPECT_EQ(parameter_probe_target(boolean, 1.0f), 0.0f);
  ParamDef slider;
  slider.min = -4;
  slider.max = 4;
  HS_EXPECT_EQ(parameter_probe_target(slider, 1.0f), -2.0f);
}
