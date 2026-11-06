/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Column-time scan effects: fidelity-off frames pinned to recorded hashes.
// ---------------------------------------------------------------------------

/** @brief Renders the selected fidelity mode through the shared golden fold. */
template <template <int, int> class E, int W, int H, int FRAMES,
          bool ENABLED = false>
inline PixelGolden column_time_parity_golden() {
  return untimed_frame_golden<E, W, H>(FRAMES, nullptr, ENABLED);
}

/** @brief One recorded fidelity-off render sequence. */
struct ColumnTimeParityCase {
  const char *name;        /**< Effect and resolution label. */
  PixelGolden (*render)(); /**< Renders and folds the frames. */
  RecordedGolden expected; /**< Recorded fidelity-off sequence. */
};

/** @brief Frames folded at 96x20. */
inline constexpr int PARITY_SMALL_FRAMES = 48;
/** @brief Frames folded at 288x144. */
inline constexpr int PARITY_FULL_FRAMES = 12;

/**
 * @brief Column-time scan effects match their recorded fidelity-off frames.
 * @details Fingerprints run across platforms; exact hashes use the recording environment.
 */
inline void test_column_time_scans_keep_fidelity_off_frames() {
  const ColumnTimeParityCase CASES[] = {
      {"HankinSolids 96x20",
       column_time_parity_golden<HankinSolids, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0xc625c70426485e58ull,
                 {{2405286481ull, 920543984ull, 991082918ull},
                  40097352575ull,
                  206105599699ull}},
        .fastmath = {}}},
      {"HankinSolids 288x144",
       column_time_parity_golden<HankinSolids, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x2f9c16a4c36f174dull,
                 {{2399661886ull, 9356570567ull, 8743403705ull},
                  1463302998150ull,
                  2940154522173ull}},
        .fastmath = {}}},
      {"IslamicStars 96x20",
       column_time_parity_golden<IslamicStars, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0x82d6d0bfb72493e9ull,
                 {{996645661ull, 635178966ull, 561335503ull},
                  21206809168ull,
                  104413817022ull}},
        .fastmath = {}}},
      {"IslamicStars 288x144",
       column_time_parity_golden<IslamicStars, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x87f42faa2cc8457dull,
                 {{54739243ull, 1270350040ull, 1363694935ull},
                  241118668929ull,
                  418037705963ull}},
        .fastmath = {}}},
      {"GnomonicStars 96x20",
       column_time_parity_golden<GnomonicStars, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0x0c36227e7f48d1b9ull,
                 {{2208959805ull, 1217740422ull, 390979467ull},
                  27163485468ull,
                  182919251671ull}},
        .fastmath = {}}},
      {"GnomonicStars 288x144",
       column_time_parity_golden<GnomonicStars, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x84af7ac582a6f7a1ull,
                 {{1310088287ull, 671996065ull, 256177046ull},
                  115244955543ull,
                  328110094610ull}},
        .fastmath = {}}},
      {"RingSpin 96x20",
       column_time_parity_golden<RingSpin, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0x16d1807f9dab9933ull,
                 {{804330226ull, 393048872ull, 122535173ull},
                  12539185369ull,
                  64921389436ull}},
        .fastmath = {}}},
      {"RingSpin 288x144",
       column_time_parity_golden<RingSpin, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x9cd2c180ca327e22ull,
                 {{644241967ull, 286983139ull, 105920220ull},
                  74155888695ull,
                  151266740548ull}},
        .fastmath = {}}},
      {"ShapeShifter 96x20",
       column_time_parity_golden<ShapeShifter, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0xa3519029bbec2b34ull,
                 {{471052378ull, 101295125ull, 18146353ull},
                  5604968973ull,
                  23316663868ull}},
        .fastmath = {}}},
      {"ShapeShifter 288x144",
       column_time_parity_golden<ShapeShifter, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x2c88e7346ec70f1eull,
                 {{2355375271ull, 519016794ull, 82298113ull},
                  211434200810ull,
                  388204021380ull}},
        .fastmath = {}}},
      {"DisplacementField 96x20",
       column_time_parity_golden<DisplacementField, 96, 20,
                                 PARITY_SMALL_FRAMES>,
       {.ieee = {0xe8cb92684225dbedull,
                 {{51572137ull, 48567174ull, 15133576ull},
                  1093337454ull,
                  5449820239ull}},
        .fastmath = {}}},
      {"DisplacementField 288x144",
       column_time_parity_golden<DisplacementField, 288, 144,
                                 PARITY_FULL_FRAMES>,
       {.ieee = {0x623a4c9b3f3e1fceull,
                 {{253982950ull, 240608237ull, 73720692ull},
                  40637933758ull,
                  81429583492ull}},
        .fastmath = {}}},
      {"Comets 96x20",
       column_time_parity_golden<Comets, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0xc28706a13cce2932ull,
                 {{34426585ull, 30933203ull, 50355801ull},
                  1046976537ull,
                  2786908034ull}},
        .fastmath = {}}},
      {"Comets 288x144",
       column_time_parity_golden<Comets, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0xab2a9200c93c77bdull,
                 {{27827681ull, 26436834ull, 43063049ull},
                  1779521037ull,
                  8551408899ull}},
        .fastmath = {}}},
      {"Galaxies 96x20",
       column_time_parity_golden<Galaxies, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0x6d13b10c6b9019c3ull,
                 {{356360678ull, 343562003ull, 350173699ull},
                  9767397711ull,
                  50337773639ull}},
        .fastmath = {}}},
      {"Galaxies 288x144",
       column_time_parity_golden<Galaxies, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x2f7ad3256dd907ffull,
                 {{234056023ull, 232552089ull, 233480476ull},
                  50010321876ull,
                  99664755720ull}},
        .fastmath = {}}},
      {"Raymarch 96x20",
       column_time_parity_golden<Raymarch, 96, 20, PARITY_SMALL_FRAMES>,
       {.ieee = {0x3e80b2c7e91bd6f5ull,
                 {{255434371ull, 182341233ull, 232767990ull},
                  6534981470ull,
                  30752386289ull}},
        .fastmath = {}}},
      {"Raymarch 288x144",
       column_time_parity_golden<Raymarch, 288, 144, PARITY_FULL_FRAMES>,
       {.ieee = {0x87e2793afac1a37full,
                 {{1794455107ull, 1235362061ull, 1619888466ull},
                  347030946155ull,
                  655972168544ull}},
        .fastmath = {}}},
  };
  for (const auto &c : CASES) {
    HS_CONTEXT(c.name);
    expect_pixel_golden(c.name, c.render(), c.expected);
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
  HS_EXPECT_TRUE((column_time_parity_golden<E, 96, 20, FRAMES, true>().hash) !=
                 (column_time_parity_golden<E, 96, 20, FRAMES, false>().hash));
}

/**
 * @brief Column-time scan effects draw differently with temporal sampling on.
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
