/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file dodecahedral_matrix_cache.h
 * @brief Runtime-owned dodecahedral matrix parity and lifetime tests. */

/** @brief Owned matrices retain their exact bytes through history and interval copies. */
void test_dodecahedral_matrix_cache_lifetime() {
  static_assert(!Pullback::Spec::CACHE_DODECAHEDRAL_MATRICES);
  static_assert(sizeof(Pullback::OptionalDodecahedralMatrices<false>) == 1);
  static_assert(sizeof(Pullback::OptionalDodecahedralMatrices<true>) == 4320);
  reset_globals();
  pin_frame_clock(0);
  AshCloud<96, 20> effect;
  effect.init();
  const auto first = ComposedFrameWhiteBox::frame(effect);
  const auto *matrices = first.dodecahedral_matrices;
  HS_EXPECT_TRUE(matrices != nullptr);
  if (!matrices)
    return;
  HS_EXPECT_TRUE(matrices != lenses::DodecahedralTable::GROUP.values.data());
  HS_EXPECT_EQ(std::memcmp(matrices,
                           lenses::DodecahedralTable::GROUP.values.data(),
                           4320),
               0);
  effect.draw_frame();
  effect.advance_display();
  const auto previous = ComposedFrameWhiteBox::current_frame(effect);
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto current = ComposedFrameWhiteBox::current_frame(effect);
  for (float fraction : {0.0f, .25f, .5f, .75f, 1.0f})
    HS_EXPECT_TRUE(Pullback::FrameInterval(previous, current)
                       .at(fraction)
                       .dodecahedral_matrices == matrices);
  HS_EXPECT_EQ(std::memcmp(matrices,
                           lenses::DodecahedralTable::GROUP.values.data(),
                           4320),
               0);
}

/** @brief One million folded directions match the flash-backed outputs byte for byte. */
void test_dodecahedral_matrix_cache_fold_parity() {
  std::array<lenses::DodecahedralTable::Matrix, 120> matrices;
  std::memcpy(matrices.data(), lenses::DodecahedralTable::GROUP.values.data(),
              sizeof(matrices));
  uint32_t random = 0x218bf41;
  const auto next = [&]() {
    random ^= random << 13;
    random ^= random >> 17;
    random ^= random << 5;
    return float(random & 0xffffff) * (2.0f / 0xffffff) - 1.0f;
  };
  bool identical = true;
  for (int sample = 0; sample < 1000000; ++sample) {
    const auto input = math::Vector(next(), next(), next()).normalized();
    const auto flash = lenses::dodecahedral_table_lens(input);
    const auto cached = lenses::dodecahedral_table_lens(input, matrices.data());
    if (std::memcmp(&flash, &cached, sizeof(flash)) != 0) {
      identical = false;
      break;
    }
  }
  HS_EXPECT_TRUE(identical);
}

/** @brief Interior lens samples use the supplied matrices and literal endpoints stay canonical. */
void test_dodecahedral_matrix_cache_lens_binding() {
  std::array<lenses::DodecahedralTable::Matrix, 120> doubled =
      lenses::DodecahedralTable::GROUP.values;
  for (auto &matrix : doubled)
    for (auto &row : matrix)
      for (float &value : row)
        value *= 2.0f;
  AshCloud<96, 20>::FrameState frame{};
  frame.dodecahedral_matrices = doubled.data();
  frame.temporal.fraction = .5f;
  const auto input = math::Vector(.1f, 1.0f, .2f).normalized();
  const auto actual = AshCloudDodecahedralLens::apply(input, frame);
  const auto expected = 2.0f * input;
  HS_EXPECT_EQ(std::memcmp(&actual, &expected, sizeof(actual)), 0);
  const auto canonical = lenses::dodecahedral_kaleidoscope_lens(input);
  for (float endpoint : {0.0f, 1.0f}) {
    frame.temporal.fraction = endpoint;
    const auto literal = AshCloudDodecahedralLens::apply(input, frame);
    HS_EXPECT_EQ(std::memcmp(&literal, &canonical, sizeof(literal)), 0);
  }
  frame.temporal.fraction = .5f;
  frame.dodecahedral_matrices = nullptr;
  const auto fallback = AshCloudDodecahedralLens::apply(input, frame);
  const auto flash = lenses::dodecahedral_table_lens(input);
  HS_EXPECT_EQ(std::memcmp(&fallback, &flash, sizeof(fallback)), 0);
}

/** @brief Face maps retain exact classifier bytes through frame history. */
void test_dodecahedral_face_cache_lifetime() {
  static_assert(!Pullback::Spec::CACHE_DODECAHEDRAL_FACES);
  static_assert(sizeof(Pullback::OptionalDodecahedralFaces<false>) == 1);
  static_assert(sizeof(Pullback::OptionalDodecahedralFaces<true>) == 24576);
  reset_globals();
  pin_frame_clock(0);
  AshCloud<96, 20> effect;
  effect.init();
  const auto first = ComposedFrameWhiteBox::frame(effect);
  const auto *faces = first.dodecahedral_faces;
  HS_EXPECT_TRUE(faces != nullptr);
  if (!faces)
    return;
  HS_EXPECT_TRUE(faces != lenses::DodecahedralTable::FACES[0]);
  for (int face = 0; face < 6; ++face)
    HS_EXPECT_EQ(std::memcmp(faces + face * 4096,
                             lenses::DodecahedralTable::FACES[face], 4096),
                 0);
  effect.draw_frame();
  effect.advance_display();
  const auto previous = ComposedFrameWhiteBox::current_frame(effect);
  pin_frame_clock(1);
  effect.draw_frame();
  effect.advance_display();
  const auto current = ComposedFrameWhiteBox::current_frame(effect);
  for (float fraction : {0.0f, .25f, .5f, .75f, 1.0f})
    HS_EXPECT_TRUE(Pullback::FrameInterval(previous, current)
                       .at(fraction)
                       .dodecahedral_faces == faces);
  for (int face = 0; face < 6; ++face)
    HS_EXPECT_EQ(std::memcmp(faces + face * 4096,
                             lenses::DodecahedralTable::FACES[face], 4096),
                 0);
}

/** @brief One million full-RAM folds match the original flash outputs byte for byte. */
void test_dodecahedral_face_cache_fold_parity() {
  const auto matrices = lenses::DodecahedralTable::GROUP.values;
  std::array<uint8_t, 24576> faces;
  for (int face = 0; face < 6; ++face)
    std::copy_n(lenses::DodecahedralTable::FACES[face], 4096,
                faces.data() + face * 4096);
  uint32_t random = 0x72fab315;
  const auto next = [&]() {
    random ^= random << 13;
    random ^= random >> 17;
    random ^= random << 5;
    return float(random & 0xffffff) * (2.0f / 0xffffff) - 1.0f;
  };
  bool identical = true;
  for (int sample = 0; sample < 1000000; ++sample) {
    const auto input = math::Vector(next(), next(), next()).normalized();
    const auto flash = lenses::dodecahedral_table_lens(input);
    const auto cached =
        lenses::dodecahedral_table_lens(input, matrices.data(), faces.data());
    if (std::memcmp(&flash, &cached, sizeof(flash))) {
      identical = false;
      break;
    }
  }
  HS_EXPECT_TRUE(identical);
}

/** @brief Face views select the RAM classifiers while null views retain flash lookup. */
void test_dodecahedral_face_cache_lens_binding() {
  auto doubled = lenses::DodecahedralTable::GROUP.values;
  for (auto &matrix : doubled)
    for (auto &row : matrix)
      for (float &value : row)
        value *= 2;
  std::array<uint8_t, 24576> rejected;
  rejected.fill(255);
  AshCloud<96, 20>::FrameState frame{};
  frame.temporal.fraction = .5f;
  frame.dodecahedral_matrices = doubled.data();
  frame.dodecahedral_faces = rejected.data();
  const auto input = math::Vector(.1f, 1.0f, .2f).normalized();
  const auto canonical = lenses::dodecahedral_kaleidoscope_lens(input);
  const auto fallback = AshCloudDodecahedralLens::apply(input, frame);
  HS_EXPECT_EQ(std::memcmp(&fallback, &canonical, sizeof(fallback)), 0);
  frame.dodecahedral_faces = nullptr;
  const auto flash = AshCloudDodecahedralLens::apply(input, frame);
  const auto expected = 2.0f * input;
  HS_EXPECT_EQ(std::memcmp(&flash, &expected, sizeof(flash)), 0);
  frame.dodecahedral_faces = rejected.data();
  for (float endpoint : {0.0f, 1.0f}) {
    frame.temporal.fraction = endpoint;
    const auto literal = AshCloudDodecahedralLens::apply(input, frame);
    HS_EXPECT_EQ(std::memcmp(&literal, &canonical, sizeof(literal)), 0);
  }
}
