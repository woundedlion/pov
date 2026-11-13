/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Sparse enum values and their parameter specification. */
namespace sparse_parameter_tests {
inline constexpr const char *LABELS[] = {"Zero", "One", "Six"};
inline constexpr const int64_t IDS[] = {0, 1, 6};
enum class SparseMode : uint8_t { ZERO = 0, ONE = 1, SIX = 6 };
inline constexpr ParamSpec<SparseMode> SPARSE_SPEC{.min = 0,
                                                   .max = 6,
                                                   .options = LABELS,
                                                   .option_count = 3,
                                                   .option_values = IDS};
} // namespace sparse_parameter_tests

/** @brief Sparse fields admit only the enumerated IDs. */
inline void test_sparse_parameter_fields() {
  using namespace sparse_parameter_tests;
  const auto spec = SPARSE_SPEC;
  struct SparseParams {
    SparseMode mode;
  };
  const auto fields = std::tuple{Control::Field<SparseParams, SparseMode>{
      "mode", &SparseParams::mode, "Mode", spec}};
  for (int64_t id : IDS)
    HS_EXPECT_TRUE(Control::valid_fields(
        SparseParams{static_cast<SparseMode>(id)}, fields));
  for (int64_t id : {2, 3, 4, 5, 7})
    HS_EXPECT_FALSE(Control::valid_fields(
        SparseParams{static_cast<SparseMode>(id)}, fields));
}

/** @brief Malformed sparse schemas reject their option values. */
inline void test_sparse_parameter_schema_rejections() {
  using namespace sparse_parameter_tests;
  const auto spec = SPARSE_SPEC;
  HS_EXPECT_TRUE(spec.valid_option_values(SparseMode::SIX));
  HS_EXPECT_FALSE(spec.valid_option_values(static_cast<SparseMode>(2)));
  static constexpr const int64_t DUPLICATE[] = {0, 6, 6};
  static constexpr const int64_t OUTSIDE[] = {0, 6, 7};
  static constexpr const int64_t MISSING_MAX[] = {0, 1, 5};
  {
    auto invalid = spec;
    invalid.option_count = 2;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  }
  {
    auto invalid = spec;
    invalid.options = nullptr;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  }
  {
    auto invalid = spec;
    invalid.option_values = DUPLICATE;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  }
  {
    auto invalid = spec;
    invalid.option_values = OUTSIDE;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  }
  static constexpr const int64_t BELOW[] = {-1, 0, 6};
  ParamSpec<int8_t> signed_spec{.min = 0,
                                .max = 6,
                                .options = LABELS,
                                .option_count = 3,
                                .option_values = IDS};
  HS_EXPECT_TRUE(signed_spec.valid_option_values(0));
  const ParamSpec<int8_t> below{.min = 0,
                                .max = 6,
                                .options = LABELS,
                                .option_count = 3,
                                .option_values = BELOW};
  HS_EXPECT_FALSE(below.valid_option_values(0));
  static constexpr const int64_t MISSING_MIN[] = {1, 2, 6};
  {
    auto invalid = spec;
    invalid.option_values = MISSING_MIN;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ONE));
  }
  {
    auto invalid = spec;
    invalid.option_values = MISSING_MAX;
    HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  }
  static constexpr const int64_t INEXACT[] = {0, 16777217, 33554432};
  ParamSpec<uint32_t> wide{.min = 0,
                           .max = 33554432,
                           .options = LABELS,
                           .option_count = 3,
                           .option_values = INEXACT};
  HS_EXPECT_FALSE(wide.valid_option_values(0));
  static constexpr const int64_t WIDE_VALID[] = {0, 1, 33554432};
  const ParamSpec<uint32_t> valid_wide{.min = 0,
                                       .max = 33554432,
                                       .options = LABELS,
                                       .option_count = 3,
                                       .option_values = WIDE_VALID};
  HS_EXPECT_TRUE(valid_wide.valid_option_values(0));
  static constexpr const int64_t TOO_WIDE[] = {0, 1, 4294967296LL};
  const ParamSpec<uint32_t> too_wide{.min = 0,
                                     .max = 4294967296LL,
                                     .options = LABELS,
                                     .option_count = 3,
                                     .option_values = TOO_WIDE};
  HS_EXPECT_FALSE(too_wide.valid_option_values(0));
  static constexpr const int64_t SIGNED_MIN_VALID[] = {INT32_MIN, 0, 6};
  ParamSpec<int32_t> signed_wide{.min = INT32_MIN,
                                 .max = 6,
                                 .options = LABELS,
                                 .option_count = 3,
                                 .option_values = SIGNED_MIN_VALID};
  HS_EXPECT_TRUE(signed_wide.valid_option_values(0));
  static constexpr const int64_t TOO_LOW[] = {INT32_MIN - 256LL, 0, 6};
  const ParamSpec<int32_t> too_low{.min = INT32_MIN - 256LL,
                                   .max = 6,
                                   .options = LABELS,
                                   .option_count = 3,
                                   .option_values = TOO_LOW};
  HS_EXPECT_FALSE(too_low.valid_option_values(0));
}

/** @brief Dense factories and float fields enforce integral option IDs. */
inline void test_sparse_parameter_factory_and_float_fields() {
  using namespace sparse_parameter_tests;
  const auto dense = ParamSpec<uint8_t>::enumerated(LABELS, 3);
  HS_EXPECT_TRUE(dense.option_values.empty());
  HS_EXPECT_EQ(dense.min, 0);
  HS_EXPECT_EQ(dense.max, 2);
  HS_EXPECT_EQ(dense.option_count, 3);
  HS_EXPECT_TRUE(dense.options == LABELS);

  struct FloatParams {
    float mode;
  };
  const ParamSpec<float> float_spec{.min = 0,
                                    .max = 6,
                                    .options = LABELS,
                                    .option_count = 3,
                                    .option_values = IDS};
  const Control::Field<FloatParams> float_field{"mode", &FloatParams::mode,
                                                "Mode", float_spec};
  for (float id : {0.0f, 1.0f, 6.0f})
    HS_EXPECT_TRUE(float_field.valid(FloatParams{id}));
  for (float id : {0.5f, 2.0f, 5.6f, 7.0f})
    HS_EXPECT_FALSE(float_field.valid(FloatParams{id}));
  auto dense_float_field = float_field;
  dense_float_field.spec.option_values = {};
  HS_EXPECT_TRUE(dense_float_field.valid(FloatParams{2.0f}));
  HS_EXPECT_FALSE(dense_float_field.valid(FloatParams{0.5f}));
}

/** @brief Sparse updates reject gaps and apply rounded, clamped IDs. */
inline void test_sparse_parameter_updates() {
  using namespace sparse_parameter_tests;
  const auto spec = SPARSE_SPEC;
  TestEffect fx(4, 4);
  SparseMode mode = SparseMode::SIX;
  float float_mode = 6;
  fx.register_param("Sparse", &mode, spec);
  fx.register_param("FloatSparse", &float_mode,
                    ParamSpec<float>{.min = 0,
                                     .max = 6,
                                     .options = LABELS,
                                     .option_count = 3,
                                     .option_values = IDS});
  const auto *parameter = required_param(fx, "Sparse");
  if (!parameter)
    return;
  HS_EXPECT_TRUE(parameter->option_values == IDS);
  for (float gap : {2.0f, 3.0f, 4.0f, 5.0f}) {
    HS_EXPECT_EQ(fx.updateParameter("Sparse", gap),
                 ParamSetResult::INADMISSIBLE);
    HS_EXPECT_EQ(fx.updateParameter("FloatSparse", gap),
                 ParamSetResult::INADMISSIBLE);
    HS_EXPECT_EQ(mode, SparseMode::SIX);
    HS_EXPECT_EQ(float_mode, 6.0f);
  }
  for (int64_t id : IDS) {
    HS_EXPECT_EQ(fx.updateParameter("Sparse", static_cast<float>(id)),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(mode, static_cast<SparseMode>(id));
  }
  HS_EXPECT_EQ(fx.updateParameter("Sparse", 5.6f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(mode, SparseMode::SIX);
  HS_EXPECT_EQ(fx.updateParameter("Sparse", 99.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(mode, SparseMode::SIX);
  HS_EXPECT_EQ(fx.updateParameter("Sparse", -99.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(mode, SparseMode::ZERO);
}

/** @brief Signed and wide enum registrations retain their descriptors and IDs. */
inline void test_sparse_integer_registration_metadata() {
  enum class SignedMode : int8_t { LOW = -4, HIGH = 7 };
  enum class WideMode : uint32_t { LOW = 0, HIGH = 4294967040u };
  static constexpr const char *LABELS[] = {"Low", "High"};
  static constexpr const char *EXPORTS[] = {"Mode::LOW", "Mode::HIGH"};
  static constexpr int64_t SIGNED_IDS[] = {-4, 7};
  static constexpr int64_t WIDE_IDS[] = {0, 4294967040LL};
  TestEffect fx(4, 4);
  SignedMode signed_mode = SignedMode::LOW;
  WideMode wide_mode = WideMode::HIGH;
  const auto generation = fx.getParameterSchemaGeneration();
  fx.register_param("Signed", &signed_mode,
                    ParamSpec<SignedMode>{.min = -4,
                                          .max = 7,
                                          .animated = true,
                                          .preset = false,
                                          .options = LABELS,
                                          .export_options = EXPORTS,
                                          .option_count = 2,
                                          .option_values = SIGNED_IDS});
  fx.register_param("Wide", &wide_mode,
                    ParamSpec<WideMode>{.min = 0,
                                        .max = 4294967040LL,
                                        .options = LABELS,
                                        .export_options = EXPORTS,
                                        .option_count = 2,
                                        .option_values = WIDE_IDS});
  HS_EXPECT_EQ(fx.getParameterSchemaGeneration(), generation + 2);
  const auto &parameters = fx.getParameters();
  HS_EXPECT_SIZE_OR_RETURN(parameters, size_t(2));
  HS_EXPECT_STREQ(parameters.begin()[0].name, "Signed");
  HS_EXPECT_STREQ(parameters.begin()[1].name, "Wide");
  HS_REQUIRED_PARAM(SIGNED, fx, "Signed");
  HS_REQUIRED_PARAM(WIDE, fx, "Wide");
  HS_EXPECT_EQ(SIGNED->target, &signed_mode);
  HS_EXPECT_EQ(SIGNED->target_type, ParamDef::TargetType::INT_I8);
  HS_EXPECT_EQ(SIGNED->min, -4.0f);
  HS_EXPECT_EQ(SIGNED->max, 7.0f);
  HS_EXPECT_EQ(SIGNED->options, LABELS);
  HS_EXPECT_EQ(SIGNED->export_options, EXPORTS);
  HS_EXPECT_EQ(SIGNED->option_values, SIGNED_IDS);
  HS_EXPECT_TRUE(SIGNED->animated);
  HS_EXPECT_FALSE(SIGNED->preset);
  HS_EXPECT_EQ(WIDE->target, &wide_mode);
  HS_EXPECT_EQ(WIDE->target_type, ParamDef::TargetType::INT_U32);
  HS_EXPECT_EQ(WIDE->max, 4294967040.0f);
  HS_EXPECT_EQ(WIDE->get(), 4294967040.0f);
  HS_EXPECT_EQ(WIDE->options, LABELS);
  HS_EXPECT_EQ(WIDE->export_options, EXPORTS);
  HS_EXPECT_EQ(WIDE->option_values, WIDE_IDS);
  HS_EXPECT_FALSE(WIDE->animated);
  HS_EXPECT_TRUE(WIDE->preset);
  HS_EXPECT_EQ(fx.updateParameter("Signed", 0.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(signed_mode, SignedMode::LOW);
  HS_EXPECT_EQ(fx.updateParameter("Signed", 7.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(signed_mode, SignedMode::HIGH);
  HS_EXPECT_EQ(fx.updateParameter("Wide", 0.0f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(wide_mode, WideMode::LOW);
}
