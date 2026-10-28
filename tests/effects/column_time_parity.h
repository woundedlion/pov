/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Column-time scan effects: fidelity-off frames pinned to recorded hashes.
// ---------------------------------------------------------------------------

/**
 * @brief Folds an effect's frames into one hash.
 * @tparam E Effect class template.
 * @tparam W Render width.
 * @tparam H Render height.
 * @tparam FRAMES Frames drawn and folded.
 * @tparam ENABLED Temporal sampling mode the frames are drawn in.
 * @return Running frame hash over every displayed frame.
 */
template <template <int, int> class E, int W, int H, int FRAMES,
          bool ENABLED = false>
inline uint64_t column_time_parity_hash() {
  reset_effect_globals();
  hs::random().seed(1337u);
  pin_frame_clock(0);
  Spacetime::set_flat_surface(!ENABLED);
  uint64_t hash = HASH_SEED64;
  {
    E<W, H> effect;
    effect.init();
    for (int f = 0; f < FRAMES; ++f) {
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      hash = frame_hash<W, H>(effect, hash);
    }
  }
  Spacetime::set_flat_surface(false);
  return hash;
}

/** @brief One recorded fidelity-off frame hash. */
struct ColumnTimeParityCase {
  const char *name;   /**< Effect and resolution label. */
  uint64_t (*hash)(); /**< Renders and folds the frames. */
  uint64_t expected;  /**< Recorded fidelity-off hash. */
};

/** @brief Frames folded at 96x20. */
inline constexpr int PARITY_SMALL_FRAMES = 48;
/** @brief Frames folded at 288x144. */
inline constexpr int PARITY_FULL_FRAMES = 12;

/**
 * @brief Column-time scan effects match their recorded fidelity-off frames.
 * @details The hashes are this toolchain's: the Windows Clang preset. Set
 * `HS_PRINT_PARITY_HASHES=1` to print the current values.
 */
inline void test_column_time_scans_keep_fidelity_off_frames() {
  const ColumnTimeParityCase CASES[] = {
      {"HankinSolids 96x20",
       column_time_parity_hash<HankinSolids, 96, 20, PARITY_SMALL_FRAMES>,
       0xc625c70426485e58ull},
      {"HankinSolids 288x144",
       column_time_parity_hash<HankinSolids, 288, 144, PARITY_FULL_FRAMES>,
       0x2f9c16a4c36f174dull},
      {"IslamicStars 96x20",
       column_time_parity_hash<IslamicStars, 96, 20, PARITY_SMALL_FRAMES>,
       0x82d6d0bfb72493e9ull},
      {"IslamicStars 288x144",
       column_time_parity_hash<IslamicStars, 288, 144, PARITY_FULL_FRAMES>,
       0x87f42faa2cc8457dull},
      {"GnomonicStars 96x20",
       column_time_parity_hash<GnomonicStars, 96, 20, PARITY_SMALL_FRAMES>,
       0x0c36227e7f48d1b9ull},
      {"GnomonicStars 288x144",
       column_time_parity_hash<GnomonicStars, 288, 144, PARITY_FULL_FRAMES>,
       0x84af7ac582a6f7a1ull},
      {"RingSpin 96x20",
       column_time_parity_hash<RingSpin, 96, 20, PARITY_SMALL_FRAMES>,
       0x16d1807f9dab9933ull},
      {"RingSpin 288x144",
       column_time_parity_hash<RingSpin, 288, 144, PARITY_FULL_FRAMES>,
       0x9cd2c180ca327e22ull},
      {"ShapeShifter 96x20",
       column_time_parity_hash<ShapeShifter, 96, 20, PARITY_SMALL_FRAMES>,
       0xba30d52d7e7978aeull},
      {"ShapeShifter 288x144",
       column_time_parity_hash<ShapeShifter, 288, 144, PARITY_FULL_FRAMES>,
       0x0059215f82610085ull},
      {"DisplacementField 96x20",
       column_time_parity_hash<DisplacementField, 96, 20, PARITY_SMALL_FRAMES>,
       0xe8cb92684225dbedull},
      {"DisplacementField 288x144",
       column_time_parity_hash<DisplacementField, 288, 144, PARITY_FULL_FRAMES>,
       0x623a4c9b3f3e1fceull},
      {"Comets 96x20",
       column_time_parity_hash<Comets, 96, 20, PARITY_SMALL_FRAMES>,
       0xc28706a13cce2932ull},
      {"Comets 288x144",
       column_time_parity_hash<Comets, 288, 144, PARITY_FULL_FRAMES>,
       0xab2a9200c93c77bdull},
      {"Galaxies 96x20",
       column_time_parity_hash<Galaxies, 96, 20, PARITY_SMALL_FRAMES>,
       0x6d13b10c6b9019c3ull},
      {"Galaxies 288x144",
       column_time_parity_hash<Galaxies, 288, 144, PARITY_FULL_FRAMES>,
       0x2f7ad3256dd907ffull},
      {"Raymarch 96x20",
       column_time_parity_hash<Raymarch, 96, 20, PARITY_SMALL_FRAMES>,
       0x3e80b2c7e91bd6f5ull},
      {"Raymarch 288x144",
       column_time_parity_hash<Raymarch, 288, 144, PARITY_FULL_FRAMES>,
       0x87e2793afac1a37full},
  };
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *print = std::getenv("HS_PRINT_PARITY_HASHES");
#pragma clang diagnostic pop
  for (const auto &c : CASES) {
    HS_CONTEXT(c.name);
    const uint64_t hash = c.hash();
    if (print && std::atoi(print) != 0)
      std::printf("  parity %-28s 0x%016llxull\n", c.name,
                  static_cast<unsigned long long>(hash));
    HS_EXPECT_EQ(hash, c.expected);
  }
}

/**
 * @brief Checks that an effect's column-time scans change its frames when
 *        temporal sampling is on.
 * @tparam E Effect class template.
 * @param name Effect name for failure context.
 */
template <template <int, int> class E>
inline void check_column_time_scans_engage(const char *name) {
  HS_CONTEXT(name);
  constexpr int FRAMES = 12;
  HS_EXPECT_TRUE((column_time_parity_hash<E, 96, 20, FRAMES, true>()) !=
                 (column_time_parity_hash<E, 96, 20, FRAMES, false>()));
}

/**
 * @brief Effects migrated to column-time scans draw differently with temporal
 *        sampling on.
 */
inline void test_column_time_scans_engage_with_fidelity_on() {
  check_column_time_scans_engage<HankinSolids>("HankinSolids");
  check_column_time_scans_engage<IslamicStars>("IslamicStars");
  check_column_time_scans_engage<GnomonicStars>("GnomonicStars");
  check_column_time_scans_engage<ShapeShifter>("ShapeShifter");
  check_column_time_scans_engage<Comets>("Comets");
  check_column_time_scans_engage<Galaxies>("Galaxies");
  check_column_time_scans_engage<DisplacementField>("DisplacementField");
  check_column_time_scans_engage<Raymarch>("Raymarch");
  check_column_time_scans_engage<RingSpin>("RingSpin");
}
