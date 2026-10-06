/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// deep_tween (global_t span)
// ============================================================================

/**
 * @brief Pins deep_tween's admissible input: an OrientationTrail, not a bare
 * Orientation.
 * @details A bare Orientation's get() yields a Quaternion, which has no
 * sub-frame history to flatten.
 */
inline void test_tweenable_rejects_bare_orientation() {
  static_assert(Tweenable<Animation::OrientationTrail<math::Orientation<8>, 8>>,
                "OrientationTrail must satisfy Tweenable");
  static_assert(!Tweenable<math::Orientation<8>>,
                "a bare Orientation must not satisfy Tweenable");
}

/** @brief Verifies Orientation tween skips only a shared motion boundary. */
inline void test_tween_orientation_skips_shared_boundary() {
  math::Orientation<5> orientation;
  std::vector<float> ts;
  tween(orientation,
        [&](const math::Quaternion &, float t) { ts.push_back(t); });
  HS_EXPECT_SIZE_OR_RETURN(ts, 1);
  HS_EXPECT_NEAR(ts.front(), 1.0f, 1e-6f);

  orientation.push(math::make_rotation(math::Z_AXIS, 0.5f));
  orientation.upsample(5);
  ts.clear();
  tween(orientation,
        [&](const math::Quaternion &, float t) { ts.push_back(t); });
  HS_EXPECT_SIZE_OR_RETURN(ts, 4);
  HS_EXPECT_NEAR(ts.front(), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < ts.size(); ++i)
    HS_EXPECT_GE(ts[i], ts[i - 1]);
}

/**
 * @brief Verifies deep_tween emits a global t spanning [0,1] across the whole
 * trail, with the expected sample count and non-decreasing t.
 * @details The emitted count is M + (N-1)*(M-1) for N frames of M sub-frames,
 * since sub-frame 0 of every frame after the first is a skipped shared
 * boundary.
 */
inline void test_deep_tween_global_t_spans_unit_interval() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  const int N = 3, M = 3;
  for (int k = 0; k < N; ++k) {
    Ori o;
    o.push(math::make_rotation(math::Z_AXIS, 0.3f * (k + 1)));
    o.upsample(M);
    trail.record(o);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const math::Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_SIZE_OR_RETURN(gts, M + (N - 1) * (M - 1));
  HS_EXPECT_NEAR(gts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < gts.size(); ++i)
    HS_EXPECT_GE(gts[i], gts[i - 1] - 1e-6f);
}

/**
 * @brief Verifies deep_tween still reaches t=1.0 when the newest frame is a
 * collapsed (motionless) length-1 frame.
 * @details deep_tween normalizes against the newest contentful frame so global
 * t reaches 1.0 at the leading edge, without re-plotting the boundary the
 * collapsed frame shares with its predecessor.
 */
inline void test_deep_tween_collapsed_newest_frame_reaches_one() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  for (int k = 0; k < 2; ++k) {
    Ori o;
    o.push(math::make_rotation(math::Z_AXIS, 0.3f * (k + 1)));
    o.upsample(3);
    trail.record(o);
  }
  Ori still; // length 1 — no motion this frame
  trail.record(still);

  std::vector<float> gts;
  deep_tween(trail,
             [&](const math::Quaternion &, float gt) { gts.push_back(gt); });

  // The collapsed tail frame is dropped, so the count matches the two
  // contentful frames (M=3 sub-frames each): M + (contentful-1)*(M-1).
  const size_t M = 3, contentful = 2;
  HS_EXPECT_SIZE_OR_RETURN(gts, M + (contentful - 1) * (M - 1));
  HS_EXPECT_NEAR(gts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
  for (size_t i = 1; i < gts.size(); ++i)
    HS_EXPECT_GE(gts[i], gts[i - 1] - 1e-6f);
}

/**
 * @brief Verifies that when every frame is motionless the lone plotted
 * orientation (the trail head) reads t = 1.0.
 * @details Mirrors tween(Orientation) for a lone snapshot.
 */
inline void test_deep_tween_all_collapsed_reaches_one() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  for (int k = 0; k < 3; ++k) {
    Ori still; // length 1 — no motion any frame
    trail.record(still);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const math::Quaternion &, float gt) { gts.push_back(gt); });

  HS_EXPECT_SIZE_OR_RETURN(gts, 1);
  HS_EXPECT_NEAR(gts.back(), 1.0f, 1e-6f);
}

/**
 * @brief Verifies deep_tween_frames groups deep_tween's emission by frame.
 * @details Callback k receives frame k's contribution — M sub-positions for
 * the first frame, M-1 (shared boundary skipped) after — and the concatenated
 * (quaternion, t) stream equals deep_tween's exactly.
 */
inline void test_deep_tween_frames_groups_flat_emission() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  const int N = 3, M = 4;
  for (int k = 0; k < N; ++k) {
    Ori o;
    o.push(math::make_rotation(math::Z_AXIS, 0.3f * (k + 1)));
    o.upsample(M);
    trail.record(o);
  }

  std::vector<std::pair<math::Quaternion, float>> flat;
  deep_tween(trail, [&](const math::Quaternion &q, float t) {
    flat.emplace_back(q, t);
  });

  size_t idx = 0;
  int frames = 0;
  deep_tween_frames(
      trail, [&](const math::Quaternion *qs, const float *ts, int count) {
        HS_EXPECT_EQ(count, frames == 0 ? M : M - 1);
        for (int i = 0; i < count && idx < flat.size(); ++i, ++idx) {
          HS_EXPECT_EQ(qs[i].r, flat[idx].first.r);
          HS_EXPECT_EQ(qs[i].v.x, flat[idx].first.v.x);
          HS_EXPECT_EQ(qs[i].v.y, flat[idx].first.v.y);
          HS_EXPECT_EQ(qs[i].v.z, flat[idx].first.v.z);
          HS_EXPECT_EQ(ts[i], flat[idx].second);
        }
        ++frames;
      });
  HS_EXPECT_EQ(idx, flat.size());
  HS_EXPECT_EQ(frames, N);
}

/**
 * @brief Verifies a motionless interior frame leaves no hole in the age ramp.
 * @details A moving / motionless / moving trail: the length-1 interior frame
 * contributes no sample and is excluded from the span, so the remaining moving
 * frames stay evenly spaced across [0,1].
 */
inline void test_deep_tween_interior_motionless_frame_no_gap() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  {
    Ori o;
    o.push(math::make_rotation(math::Z_AXIS, 0.3f));
    o.upsample(3);
    trail.record(o);
  }
  {
    Ori still; // length 1 — no motion this interior frame
    trail.record(still);
  }
  {
    Ori o;
    o.push(math::make_rotation(math::Z_AXIS, 0.6f));
    o.upsample(3);
    trail.record(o);
  }

  std::vector<float> gts;
  deep_tween(trail,
             [&](const math::Quaternion &, float gt) { gts.push_back(gt); });

  const float expected[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  HS_EXPECT_SIZE_OR_RETURN(gts, 5);
  for (size_t i = 0; i < gts.size(); ++i)
    HS_EXPECT_NEAR(gts[i], expected[i], 1e-6f);
}

/** @brief Verifies a motionless oldest frame occupies its age slot endpoint. */
inline void test_deep_tween_oldest_motionless_frame_no_gap() {
  using Ori = math::Orientation<8>;
  Animation::OrientationTrail<Ori, 8> trail;
  trail.record(Ori());
  Ori moving;
  moving.push(math::make_rotation(math::Z_AXIS, 0.6f));
  moving.upsample(3);
  trail.record(moving);

  std::vector<float> gts;
  deep_tween(trail,
             [&](const math::Quaternion &, float gt) { gts.push_back(gt); });
  HS_EXPECT_SIZE_OR_RETURN(gts, 3);
  HS_EXPECT_NEAR(gts[0], 0.5f, 1e-6f);
  HS_EXPECT_NEAR(gts[1], 0.75f, 1e-6f);
  HS_EXPECT_NEAR(gts[2], 1.0f, 1e-6f);
}

/**
 * @brief Verifies a single-sample VectorTrail reads t = 1.0 (the lone trail
 * head), while a multi-sample sweep ramps 0 -> 1 oldest -> newest.
 * @details A trail with its first recorded point mirrors tween(Orientation)
 * for a lone snapshot.
 */
inline void test_tween_vectortrail_single_sample_reaches_one() {
  Animation::VectorTrail<8> trail;
  trail.record(math::Vector(1, 0, 0));

  std::vector<float> ts;
  tween(trail, [&](const math::Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_SIZE_OR_RETURN(ts, 1);
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);

  trail.record(math::Vector(0, 1, 0));
  trail.record(math::Vector(0, 0, 1));
  ts.clear();
  tween(trail, [&](const math::Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_SIZE_OR_RETURN(ts, 3);
  HS_EXPECT_NEAR(ts.front(), 0.0f, 1e-6f); // oldest = tail
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);  // newest = head
}

/**
 * @brief Verifies QuantizedVectorTrail round-trips the sampled unit vectors
 * within SNORM3_COMPONENT_BOUND per component, preserves 0 and ±1 exactly, clamps
 * out-of-domain components, and keeps Trail's oldest-first ring semantics.
 */
inline void test_quantized_vector_trail_roundtrip_and_ring() {
  constexpr float QUANT_ERR = SNORM3_COMPONENT_BOUND;

  Animation::QuantizedVectorTrail<8> trail;
  trail.record(math::Vector(1, 0, 0));
  trail.record(math::Vector(0, -1, 0));
  HS_EXPECT_EQ(trail.get(0).x, 1.0f);
  HS_EXPECT_EQ(trail.get(0).y, 0.0f);
  HS_EXPECT_EQ(trail.get(1).y, -1.0f);

  for (int i = 0; i < 32; ++i) {
    math::Vector v =
        math::Vector::from_spherical(0.37f + 0.19f * i, 0.11f + 0.09f * i);
    trail.record(v);
    math::Vector d = trail.get(trail.length() - 1);
    HS_EXPECT_NEAR(d.x, v.x, QUANT_ERR);
    HS_EXPECT_NEAR(d.y, v.y, QUANT_ERR);
    HS_EXPECT_NEAR(d.z, v.z, QUANT_ERR);
  }

  trail.clear();
  trail.record(math::Vector(1.5f, -2.0f, 0.25f));
  HS_EXPECT_EQ(trail.get(0).x, 1.0f);
  HS_EXPECT_EQ(trail.get(0).y, -1.0f);
  HS_EXPECT_NEAR(trail.get(0).z, 0.25f, QUANT_ERR);

  Animation::QuantizedVectorTrail<4> ring;
  for (int i = 0; i < 6; ++i)
    ring.record(math::Vector(0, 0, 0.1f * i));
  HS_EXPECT_EQ(ring.length(), static_cast<size_t>(4));
  HS_EXPECT_NEAR(ring.get(0).z, 0.2f, QUANT_ERR); // oldest retained = 3rd
  HS_EXPECT_NEAR(ring.get(3).z, 0.5f, QUANT_ERR); // newest = last recorded
  ring.expire();
  HS_EXPECT_EQ(ring.length(), static_cast<size_t>(3));
  HS_EXPECT_NEAR(ring.get(0).z, 0.3f, QUANT_ERR);

  std::vector<float> ts;
  tween(ring, [&](const math::Vector &, float t) { ts.push_back(t); });
  HS_EXPECT_SIZE_OR_RETURN(ts, 3);
  HS_EXPECT_NEAR(ts.front(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ts[1], 0.5f, 1e-6f);
  HS_EXPECT_NEAR(ts.back(), 1.0f, 1e-6f);
}
