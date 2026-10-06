/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_canvas.h.

inline void test_sparse_parameter_specs() {
  static constexpr const char *LABELS[] = {"Zero", "One", "Six"};
  static constexpr const int64_t IDS[] = {0, 1, 6};
  enum class SparseMode : uint8_t { ZERO = 0, ONE = 1, SIX = 6 };
  ParamSpec<SparseMode> spec{.min = 0,
                             .max = 6,
                             .options = LABELS,
                             .option_count = 3,
                             .option_values = IDS};
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
  HS_EXPECT_TRUE(spec.valid_option_values(SparseMode::SIX));
  HS_EXPECT_FALSE(spec.valid_option_values(static_cast<SparseMode>(2)));
  auto invalid = spec;
  invalid.option_count = 2;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  invalid = spec;
  invalid.options = nullptr;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  static constexpr const int64_t DUPLICATE[] = {0, 6, 6};
  invalid = spec;
  invalid.option_values = DUPLICATE;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  static constexpr const int64_t OUTSIDE[] = {0, 6, 7};
  invalid.option_values = OUTSIDE;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  static constexpr const int64_t BELOW[] = {-1, 0, 6};
  ParamSpec<int8_t> signed_spec{.min = 0,
                                .max = 6,
                                .options = LABELS,
                                .option_count = 3,
                                .option_values = IDS};
  HS_EXPECT_TRUE(signed_spec.valid_option_values(0));
  signed_spec.option_values = BELOW;
  HS_EXPECT_FALSE(signed_spec.valid_option_values(0));
  static constexpr const int64_t MISSING_MIN[] = {1, 2, 6};
  invalid.option_values = MISSING_MIN;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ONE));
  static constexpr const int64_t MISSING_MAX[] = {0, 1, 5};
  invalid.option_values = MISSING_MAX;
  HS_EXPECT_FALSE(invalid.valid_option_values(SparseMode::ZERO));
  static constexpr const int64_t INEXACT[] = {0, 16777217, 33554432};
  ParamSpec<uint32_t> wide{.min = 0,
                           .max = 33554432,
                           .options = LABELS,
                           .option_count = 3,
                           .option_values = INEXACT};
  HS_EXPECT_FALSE(wide.valid_option_values(0));
  static constexpr const int64_t WIDE_VALID[] = {0, 1, 33554432};
  wide.option_values = WIDE_VALID;
  HS_EXPECT_TRUE(wide.valid_option_values(0));
  static constexpr const int64_t TOO_WIDE[] = {0, 1, 4294967296LL};
  wide.max = 4294967296LL;
  wide.option_values = TOO_WIDE;
  HS_EXPECT_FALSE(wide.valid_option_values(0));
  static constexpr const int64_t SIGNED_MIN_VALID[] = {INT32_MIN, 0, 6};
  ParamSpec<int32_t> signed_wide{.min = INT32_MIN,
                                 .max = 6,
                                 .options = LABELS,
                                 .option_count = 3,
                                 .option_values = SIGNED_MIN_VALID};
  HS_EXPECT_TRUE(signed_wide.valid_option_values(0));
  static constexpr const int64_t TOO_LOW[] = {INT32_MIN - 256LL, 0, 6};
  signed_wide.min = INT32_MIN - 256LL;
  signed_wide.option_values = TOO_LOW;
  HS_EXPECT_FALSE(signed_wide.valid_option_values(0));
  HS_EXPECT_TRUE(
      ParamSpec<uint8_t>::enumerated(LABELS, 3).valid_option_values(2));

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
  HS_EXPECT_TRUE(fx.getParameters().find("Sparse")->option_values == IDS);
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
