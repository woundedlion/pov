/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// ShapeShifter slider contract and preset rows.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies ShapeShifter's slider contract and preset-row invariants.
 * @details Preset values are not pinned. Boot equals the first row, every
 *          value sits inside its registered range, the non-preset Alpha slider
 *          survives a selection, and every row is a distinct parameter vector.
 */
inline void test_shapeshifter_preset_defaults() {
  using SS = ShapeShifter<DEFAULT_W, DEFAULT_H>;
  reset_globals();
  SS ss;
  ss.init();

  auto value = [&](const char *name) {
    const auto *def = hs_test::required_param(ss, name);
    return def ? def->get() : -1.0f;
  };
  auto preset_row = [&] {
    std::vector<float> row;
    for (const auto &def : ss.getParameters())
      if (def.preset)
        row.push_back(def.get());
    return row;
  };

  hs_test::expect_params_in_range(ss, "boot state");
  HS_EXPECT_EQ(value("Alpha"), 1.0f); // boots fully opaque
  const std::vector<float> boot_row = preset_row();

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

  const auto *alpha = hs_test::required_param(ss, "Alpha");
  const auto *shape = hs_test::required_param(ss, "Shape");
  const auto *falloff = hs_test::required_param(ss, "Alpha Falloff");
  const auto *spacing = hs_test::required_param(ss, "Spacing");
  const auto *count = hs_test::required_param(ss, "Count");
  HS_EXPECT_TRUE(hs_test::required_param(ss, "Speed") != nullptr);
  if (count)
    HS_EXPECT_EQ(count->max, static_cast<float>(SS::DRAW_LIMIT));
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

  HS_EXPECT_GT(ss.getPresetCount(), size_t{0});
  std::vector<std::vector<float>> rows;
  for (size_t i = 0; i < ss.getPresetCount(); ++i) {
    HS_CONTEXT("preset", static_cast<int>(i));
    ss.profile_select_preset(i);
    hs_test::expect_params_in_range(ss, "preset row");
    HS_EXPECT_EQ(value("Alpha"), 0.37f); // a preset never writes a non-preset
    HS_EXPECT_TRUE(ss.animations_paused());
    if (shape)
      HS_EXPECT_LT(value("Shape"), static_cast<float>(shape->option_count));
    rows.push_back(preset_row());
  }
  HS_EXPECT_TRUE(rows.front() == boot_row);

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
    reset_globals();
    SS ss;
    ss.init();
    ss.setAnimationsPaused(true);
    HS_EXPECT_EQ(ss.updateParameter(slider, static_cast<float>(selection)),
                 ParamSetResult::APPLIED);
    ss.draw_frame();
    ss.advance_display();

    const uint64_t acc = frame_energy<SMALL_W, SMALL_H>(ss);
    const uint64_t fold = frame_hash<SMALL_W, SMALL_H>(ss);
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
