/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @brief Advances an orientation one frame through three sub-steps about
 *        changing axes, as a motion animation's upsampled history does.
 * @param orientation Orientation to advance; its history becomes the frame's.
 * @param frame Frame number, varying the axes.
 * @param angle Rotation per sub-step, radians.
 */
inline void advance_turning_history(math::Orientation<> &orientation, int frame,
                                    float angle) {
  orientation.collapse();
  for (int step = 1; step <= 3; ++step) {
    const float phase = 0.9f * static_cast<float>(frame) + 1.7f * step;
    const math::Vector axis =
        math::Vector(cosf(phase), 0.6f, sinf(phase)).normalized();
    orientation.push(math::make_rotation(axis, angle) * orientation.get());
  }
}

/**
 * @brief Largest elementwise difference between two rotation matrices.
 * @param a,b The matrices.
 * @return The largest absolute entry of `a - b`.
 */
inline float matrix_distance(const math::RotationMatrix &a,
                             const math::RotationMatrix &b) {
  float worst = 0.0f;
  for (const math::Vector d : {a.r0 - b.r0, a.r1 - b.r1, a.r2 - b.r2})
    worst = fmaxf(worst, fmaxf(fabsf(d.x), fmaxf(fabsf(d.y), fabsf(d.z))));
  return worst;
}

/**
 * @brief Whether two rotation matrices are bit-identical.
 * @param a,b The matrices.
 * @return True when every entry matches exactly.
 */
inline bool same_matrix(const math::RotationMatrix &a,
                        const math::RotationMatrix &b) {
  return std::memcmp(&a, &b, sizeof(a)) == 0;
}

/**
 * @brief Whether a table entry holds a pose's matrix: its stored first two
 *        rows bit-identical, its third, their cross product, within rounding.
 * @param entry Table entry.
 * @param pose Matrix expanded from the pose.
 * @return True when the entry holds the pose.
 */
inline bool same_pose(const math::RotationMatrix &entry,
                      const math::RotationMatrix &pose) {
  const math::Vector third = entry.r2 - pose.r2;
  return std::memcmp(&entry.r0, &pose.r0, sizeof(pose.r0)) == 0 &&
         std::memcmp(&entry.r1, &pose.r1, sizeof(pose.r1)) == 0 &&
         math::dot(third, third) <= 1e-12f;
}

/** @brief The pose table follows the sub-frame history at every column time
 *         and starts each frame where the last ended. */
inline void test_pose_table_follows_history_path() {
  for (const int width : {288, 96}) {
    HS_CONTEXT("width", width);
    const int half = width / 2;
    ScratchScope storage(scratch_arena_a);
    auto table = std::make_unique<Spacetime::PoseTable>();
    const size_t before = scratch_arena_a.get_offset();
    table->bind(scratch_arena_a, width);
    HS_EXPECT_LE(scratch_arena_a.get_offset() - before,
                 Spacetime::PoseTable::storage_bytes(width));
    math::Orientation<> orientation(math::make_rotation(math::X_AXIS, 0.3f));
    math::RotationMatrix previous_end;
    float history_error = 0.0f, endpoint_departure = 0.0f;
    for (int frame = 0; frame < 6; ++frame) {
      advance_turning_history(orientation, frame, 0.09f);
      table->build(orientation);
      HS_EXPECT_EQ(table->times(), half);
      HS_EXPECT_FALSE(table->stationary());
      const math::SlerpArc endpoints(orientation.get(0), orientation.get());
      for (int c = 0; c < half; ++c) {
        const float t = Spacetime::surface_time(c, width);
        const math::RotationMatrix want(orientation.get(t));
        history_error =
            fmaxf(history_error, matrix_distance(table->at(c), want));
        endpoint_departure = fmaxf(
            endpoint_departure,
            matrix_distance(math::RotationMatrix(endpoints.sample(t)), want));
        HS_EXPECT_TRUE(same_matrix(table->column(c), table->at(c)));
        HS_EXPECT_TRUE(same_matrix(table->column(c + half), table->at(c)));
      }
      HS_EXPECT_TRUE(
          same_pose(table->at(half), math::RotationMatrix(orientation.get())));
      HS_EXPECT_TRUE(
          same_pose(table->at(0), math::RotationMatrix(orientation.get(0))));
      if (frame > 0)
        HS_EXPECT_TRUE(same_matrix(table->at(0), previous_end));
      previous_end = table->at(half);
    }
    HS_EXPECT_LT(history_error, 2e-6f);
    // The history bends away from its endpoint slerp, so an endpoint-only
    // table fails the comparison above.
    HS_EXPECT_GT(endpoint_departure, 100.0f * 2e-6f);
  }
}

/** @brief A one-pose history fills the table with that pose, still. */
inline void test_pose_table_holds_a_still_history() {
  constexpr int W = 288;
  ScratchScope storage(scratch_arena_a);
  auto table = std::make_unique<Spacetime::PoseTable>();
  table->bind(scratch_arena_a, W);
  const math::Orientation<> still(math::make_rotation(math::Z_AXIS, 1.1f));
  table->build(still);
  HS_EXPECT_TRUE(table->stationary());
  const math::RotationMatrix pose(still.get());
  for (int c = 0; c <= W / 2; ++c)
    HS_EXPECT_TRUE(same_pose(table->at(c), pose));
}

/** @brief Rotation through the table is exact at column times and tracks the
 *         history path between them. */
inline void test_pose_table_rotates_between_column_times() {
  constexpr int W = 288;
  ScratchScope storage(scratch_arena_a);
  auto table = std::make_unique<Spacetime::PoseTable>();
  table->bind(scratch_arena_a, W);
  math::Orientation<> orientation;
  advance_turning_history(orientation, 2, 0.12f);
  table->build(orientation);
  const math::Vector v = math::Vector(0.3f, -0.5f, 0.81f).normalized();
  float worst = 0.0f;
  for (int c = 0; c <= W / 2; ++c) {
    const float t = static_cast<float>(c) / (W / 2);
    const math::Vector exact = table->at(c).apply(v);
    const math::Vector got = table->rotate(v, t);
    HS_EXPECT_EQ(got.x, exact.x);
    HS_EXPECT_EQ(got.y, exact.y);
    HS_EXPECT_EQ(got.z, exact.z);
    if (c < W / 2) {
      const float between = (static_cast<float>(c) + 0.37f) / (W / 2);
      worst = fmaxf(worst, (table->rotate(v, between) -
                            math::rotate(v, orientation.get(between)))
                               .magnitude());
    }
  }
  HS_EXPECT_LT(worst, 2e-6f);
}

/** @brief The column pose generator samples a two-pose history as
 *         math::SlerpArc does. */
inline void test_column_poses_match_endpoint_slerp() {
  constexpr int W = 288;
  const math::Quaternion ends[] = {
      math::make_rotation(math::Y_AXIS, 0.2f),
      -math::make_rotation(math::Vector(0.3f, 1.0f, 0.2f).normalized(), 0.5f)};
  const math::SlerpArc arc(ends[0], ends[1]);
  int visited = 0;
  float worst = 0.0f;
  Spacetime::for_each_time_pose(
      ends, 2, W, [&](int c, const math::Quaternion &pose) {
        HS_EXPECT_EQ(c, visited);
        ++visited;
        const math::Quaternion want =
            arc.sample(Spacetime::SamplingSurface<W>::time(c));
        worst = fmaxf(worst, sqrtf((pose - want).squared_magnitude()));
      });
  HS_EXPECT_EQ(visited, W / 2);
  HS_EXPECT_LT(worst, 1e-6f);
}
