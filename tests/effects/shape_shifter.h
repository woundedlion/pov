/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// ShapeShifter slider contract and preset rows.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies ShapeShifter's slider contract and preset-row invariants.
 * @details Preset magnitudes are not pinned. Every value must sit inside its
 *          registered range, each row keeps its structural selections and
 *          shape/falloff pairing, the non-preset Alpha slider survives a
 *          selection, and every row is a distinct parameter vector.
 */
inline void test_shapeshifter_preset_defaults() {
  reset_effect_globals();
  ShapeShifter<DEFAULT_W, DEFAULT_H> ss;
  ss.init();

  auto value = [&](const char *name) {
    for (const auto &def : ss.getParameters())
      if (std::strcmp(def.name, name) == 0)
        return def.get();
    HS_EXPECT(false, "ShapeShifter parameter is missing");
    return -1.0f;
  };

  // Preset rows are assigned straight into params, bypassing register_param's
  // range check.
  auto expect_in_range = [&](const char *label) {
    HS_CONTEXT(label);
    for (const auto &def : ss.getParameters()) {
      HS_CONTEXT(def.name);
      const float v = def.get();
      HS_EXPECT_TRUE(std::isfinite(v));
      HS_EXPECT_GE(v, def.min);
      HS_EXPECT_LE(v, def.max);
      if (def.option_count > 0) {
        HS_EXPECT_EQ(v, std::floor(v));
        HS_EXPECT_LT(v, static_cast<float>(def.option_count));
      }
    }
  };

  expect_in_range("boot state");
  HS_EXPECT_EQ(value("Alpha"), 1.0f); // boots fully opaque
  HS_EXPECT_EQ(value("Shape"), 3.0f);
  HS_EXPECT_EQ(value("Spacing"), 1.0f);
  HS_EXPECT_EQ(value("Function"),
               static_cast<float>(
                   ShapeShifter<DEFAULT_W, DEFAULT_H>::PhaseFunction::SINE));
  HS_EXPECT_EQ(value("Opposite"), 0.0f);
  HS_EXPECT_EQ(value("Alpha Falloff"), 1.0f);

  const char *expected_export_order[] = {
      "Shape", "Count",    "Sides",         "Function", "Amplitude",
      "Speed", "Opposite", "Alpha Falloff", "Spacing"};
  size_t export_index = 0;
  for (const auto &def : ss.getParameters()) {
    if (!def.preset)
      continue;
    HS_EXPECT(export_index < std::size(expected_export_order),
              "ShapeShifter exports an unexpected parameter");
    if (export_index < std::size(expected_export_order))
      HS_EXPECT_EQ(std::string_view(def.name),
                   std::string_view(expected_export_order[export_index]));
    ++export_index;
  }
  HS_EXPECT_EQ(export_index, std::size(expected_export_order));

  const auto *alpha = ss.getParameters().find("Alpha");
  const auto *shape = ss.getParameters().find("Shape");
  const auto *falloff = ss.getParameters().find("Alpha Falloff");
  const auto *spacing = ss.getParameters().find("Spacing");
  const auto *count = ss.getParameters().find("Count");
  const auto *speed = ss.getParameters().find("Speed");
  HS_EXPECT(alpha != nullptr, "ShapeShifter Alpha parameter is missing");
  HS_EXPECT(shape != nullptr, "ShapeShifter Shape parameter is missing");
  HS_EXPECT(falloff != nullptr,
            "ShapeShifter Alpha Falloff parameter is missing");
  HS_EXPECT(spacing != nullptr, "ShapeShifter Spacing parameter is missing");
  HS_EXPECT(count != nullptr, "ShapeShifter Count parameter is missing");
  HS_EXPECT(speed != nullptr, "ShapeShifter Speed parameter is missing");
  if (count)
    HS_EXPECT_EQ(count->max, 288.0f);
  if (speed)
    HS_EXPECT_EQ(speed->max, 0.16f);
  if (alpha)
    HS_EXPECT_FALSE(alpha->preset);
  if (shape) {
    HS_EXPECT_TRUE(shape->is_enum());
    HS_EXPECT_EQ(std::string_view(shape->export_options[3]),
                 std::string_view("ShapeType::PLANAR_STAR"));
    HS_EXPECT_EQ(std::string_view(shape->export_options[4]),
                 std::string_view("ShapeType::SPHERICAL_STAR"));
  }
  if (falloff) {
    HS_EXPECT_TRUE(falloff->is_enum());
    HS_EXPECT_EQ(std::string_view(falloff->export_options[1]),
                 std::string_view("AlphaFalloff::TOWARD_EQUATOR"));
  }
  if (spacing) {
    HS_EXPECT_TRUE(spacing->is_enum());
    HS_EXPECT_EQ(std::string_view(spacing->export_options[1]),
                 std::string_view("RadiusSpacing::SCREEN_BALANCED"));
  }

  HS_EXPECT_EQ(ss.updateParameter("Alpha", 0.37f), ParamSetResult::APPLIED);

  // Structural selections: which primitive and falloff each row draws.
  const float expected_shapes[] = {3.0f, 1.0f, 3.0f, 2.0f, 3.0f,
                                   1.0f, 1.0f, 1.0f, 2.0f};
  const float expected_falloffs[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f,
                                     0.0f, 0.0f, 0.0f, 0.0f};
  HS_EXPECT_EQ(std::size(expected_shapes), ss.getPresetCount());
  HS_EXPECT_EQ(std::size(expected_falloffs), ss.getPresetCount());
  std::vector<std::vector<float>> rows;
  for (size_t i = 0; i < std::size(expected_shapes); ++i) {
    HS_CONTEXT("preset", static_cast<int>(i));
    ss.profile_select_preset(i);
    expect_in_range("preset row");
    HS_EXPECT_EQ(value("Alpha"), 0.37f); // a preset never writes a non-preset
    HS_EXPECT_TRUE(ss.animations_paused());
    HS_EXPECT_EQ(value("Shape"), expected_shapes[i]);
    HS_EXPECT_EQ(value("Alpha Falloff"), expected_falloffs[i]);
    HS_EXPECT_EQ(value("Shape") == 3.0f, value("Alpha Falloff") == 1.0f);

    std::vector<float> row;
    for (const auto &def : ss.getParameters())
      if (def.preset)
        row.push_back(def.get());
    rows.push_back(row);
  }

  // Identical parameter vectors would be one preset visited twice.
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = i + 1; j < rows.size(); ++j) {
      HS_CONTEXT("preset pair", static_cast<int>(i), static_cast<int>(j));
      HS_EXPECT(rows[i] != rows[j],
                "each ShapeShifter preset must be distinct");
    }
}

/**
 * @brief Renders every Shape and Function slider selection.
 * @details Each primitive is exercised at radii on both sides of the antipode
 * fold; two selections producing the same frame means the switch did not
 * dispatch on them.
 */
inline void test_shapeshifter_slider_selections_render() {
  using SS = ShapeShifter<SMALL_W, SMALL_H>;

  auto render = [](const char *slider, int selection) {
    reset_effect_globals();
    SS ss;
    ss.init();
    ss.setAnimationsPaused(true);
    HS_EXPECT_EQ(ss.updateParameter(slider, static_cast<float>(selection)),
                 ParamSetResult::APPLIED);
    ss.draw_frame();
    ss.advance_display();

    const uint64_t acc = frame_energy<SMALL_W, SMALL_H>(ss);
    uint64_t fold = hs_test::HASH_SEED64;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x) {
        const Pixel &pixel = ss.get_pixel(x, y);
        for (uint16_t channel : {pixel.r, pixel.g, pixel.b})
          fold = hs_test::fnv1a64_channel(fold, channel);
      }
    HS_EXPECT_GT(acc, 0u);
    return fold;
  };

  auto sweep = [&](const char *slider, int selections) {
    std::vector<uint64_t> folds;
    for (int selection = 0; selection < selections; ++selection)
      folds.push_back(render(slider, selection));
    for (int i = 0; i < selections; ++i)
      for (int j = i + 1; j < selections; ++j) {
        if (folds[i] == folds[j])
          std::printf("  SHAPESHIFTER %s selections %d and %d render the same "
                      "frame (fold %llu)\n",
                      slider, i, j, static_cast<unsigned long long>(folds[i]));
        HS_EXPECT(folds[i] != folds[j],
                  "each slider selection must render a distinct frame");
      }
  };

  sweep("Shape", SS::NUM_SHAPES);
  sweep("Function", SS::NUM_FUNCTIONS);
}
