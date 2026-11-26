/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Fields test cases. */

#include "tests/transformers/fixture.h"

namespace hs_test::transformers_tests {

/**
 * @brief Verifies a field manager with nothing spawned evaluates to 0
 *        everywhere and reports a 0 bound.
 */
inline void test_field_transformer_no_entities_is_zero() {
  Timeline tl;
  TestFieldTransformer ft(tl);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ft.field_bound(), 0.0f, 1e-6f);
  HS_EXPECT_EQ(ft.active_count(), 0);
}

/**
 * @brief Verifies spawned entities superpose by summation, scaled by the
 *        sample point, and field_bound sums the per-entity bounds.
 */
inline void test_field_transformer_sums_and_bounds() {
  Timeline tl;
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), 2.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.field(math::Vector(-0.5f, 0, 0)), -1.0f, 1e-5f);

  HS_EXPECT_TRUE(ft.spawn(0, -3.0f, 100) != nullptr);
  HS_EXPECT_EQ(ft.active_count(), 2);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), -1.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.operator()(math::Vector(1, 0, 0)), -1.0f, 1e-5f);
  HS_EXPECT_NEAR(ft.field_bound(), 5.0f, 1e-5f);
}

/**
 * @brief Verifies active_params exposes spawned entities in spawn order.
 */
inline void test_field_transformer_active_params_spawn_order() {
  Timeline tl;
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, -3.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.active_params(0).value, 2.0f, 1e-6f);
  HS_EXPECT_NEAR(ft.active_params(1).value, -3.0f, 1e-6f);
}

/**
 * @brief Verifies a completed field entity's pool slot is reclaimed and a
 *        fresh spawn reuses it.
 * @details With both CAPACITY 2 slots filled, the post-completion spawn
 *          succeeds only if the then() reclaim fired.
 */
inline void test_field_transformer_slot_reclaimed() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);

  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 2) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, 5.0f, 100) != nullptr);
  HS_EXPECT_TRUE(ft.spawn(0, 9.0f, 100) == nullptr);

  for (int i = 0; i < 6; ++i)
    tl.step(cv);

  HS_EXPECT_EQ(ft.active_count(), 1);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), 5.0f, 1e-5f);
  HS_EXPECT_TRUE(ft.spawn(0, 1.0f, 100) != nullptr);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), 6.0f, 1e-5f);
}

/**
 * @brief Verifies reclaim_storage() carries a live entity across an arena
 *        rewind: the re-claimed slots keep their addresses and contents.
 * @details Rewinds via set_offset to the pre-init_storage mark, then replays
 *          the pool allocation.
 */
inline void test_field_transformer_storage_survives_arena_rewind() {
  Timeline tl;
  const size_t mark = persistent_arena.get_offset();
  TestFieldTransformer ft(tl);
  ft.init_storage(persistent_arena);
  HS_EXPECT_TRUE(ft.spawn(0, 2.0f, 100) != nullptr);

  persistent_arena.set_offset(mark);
  ft.reclaim_storage(persistent_arena);
  HS_EXPECT_EQ(ft.active_count(), 1);
  HS_EXPECT_NEAR(ft.field(math::Vector(1, 0, 0)), 2.0f, 1e-5f);
}

// ============================================================================
// DominantFieldAccumulator — pass-through, no stacking, smooth sign crossing
// ============================================================================

/**
 * @brief Verifies the documented blend identities: nothing added reports 0, a
 *        lone field passes through exactly, and repeats of one value neither
 *        stack nor shift.
 */
inline void test_dominant_field_identities() {
  HS_EXPECT_EQ(DominantFieldAccumulator().value(), 0.0f);

  const float SINGLES[] = {0.7f, -0.7f, 4.5f, -0.013f, 1.0f};
  for (float s : SINGLES) {
    DominantFieldAccumulator acc;
    acc.add(s);
    HS_EXPECT_NEAR(acc.value(), s, 1e-6f);
  }

  for (float s : {0.7f, -2.5f}) {
    DominantFieldAccumulator acc;
    for (int n = 1; n <= 5; ++n) {
      acc.add(s);
      HS_EXPECT_NEAR(acc.value(), s, 1e-6f);
    }
  }

  DominantFieldAccumulator tiny;
  tiny.add(1e-5f);
  HS_EXPECT_EQ(tiny.value(), 0.0f);
  DominantFieldAccumulator above_floor;
  above_floor.add(1e-3f);
  HS_EXPECT_NEAR(above_floor.value(), 1e-3f, 1e-9f);
}

/**
 * @brief Verifies the strongest contribution dominates a mixed overlap without
 *        the two stacking.
 */
inline void test_dominant_field_strongest_wins_without_stacking() {
  DominantFieldAccumulator acc;
  acc.add(3.0f);
  acc.add(1.0f);
  const float blended = acc.value();
  HS_EXPECT_NEAR(blended, 2.8f, 1e-5f);
  HS_EXPECT_LT(blended, 3.0f); // not the plain max
  HS_EXPECT_GT(blended, 2.0f); // but far nearer the strong field

  DominantFieldAccumulator lopsided;
  lopsided.add(10.0f);
  lopsided.add(0.1f);
  HS_EXPECT_NEAR(lopsided.value(), 10.0f, 0.05f);
}

/**
 * @brief Verifies the blend crosses an opposite-signed cancellation smoothly
 *        rather than jumping the way a hard max by magnitude would.
 * @details A fixed +1 field is swept against a partner from -1.5 to -0.5,
 *          across which a hard max would step by 2. Every consecutive step
 *          must stay far below that, the output must change sign, and it must
 *          move substantially overall.
 */
inline void test_dominant_field_sign_crossing_is_continuous() {
  const int STEPS = 1000;
  const float LO = -1.5f, HI = -0.5f;

  auto blend = [](float b) {
    DominantFieldAccumulator acc;
    acc.add(1.0f);
    acc.add(b);
    return acc.value();
  };

  float prev = blend(LO);
  const float first = prev;
  float worst_step = 0.0f;
  for (int i = 1; i <= STEPS; ++i) {
    const float b = LO + (HI - LO) * (float)i / (float)STEPS;
    const float cur = blend(b);
    worst_step = fold_worst(worst_step, std::fabs(cur - prev));
    prev = cur;
  }
  const float last = prev;

  HS_EXPECT_LT(worst_step, 0.01f);
  HS_EXPECT_LT(first, 0.0f);
  HS_EXPECT_GT(last, 0.0f);
  HS_EXPECT_GT(last - first, 0.5f);
  HS_EXPECT_NEAR(blend(-1.0f), 0.0f, 1e-6f);
}

// ============================================================================
// bump_field — cap falloff, fast reject, envelope gating
// ============================================================================

/** @brief Verifies the bump threshold tracks the effective footprint. */
inline void test_bump_field_threshold_sync() {
  Animation::BumpParams p;
  p.radius = 0.8f;
  p.envelope = 0.25f;
  p.sync();
  HS_EXPECT_NEAR(p.cos_radius, std::cos(0.2f), 1e-6f);

  p.envelope = 0.75f;
  p.sync();
  HS_EXPECT_NEAR(p.cos_radius, std::cos(0.6f), 1e-6f);
}

/**
 * @brief Verifies the drape push along the axis meridian: zero for the ring
 *        through the center, peak bow between center and edge, zero at the
 *        footprint edge, exactly 0 outside the fast-reject cap.
 * @details On the meridian a point at offset y pushes by
 *          sign(y) * (R - |y|) * sin(pi |y| / R).
 */
inline void test_bump_field_drapes_over_ball() {
  const float C_LAT = math::PI_F / 3.0f;
  Animation::BumpParams p;
  p.center = math::Vector(std::sin(C_LAT), std::cos(C_LAT), 0.0f);
  p.axis = math::Vector(0, 1, 0);
  p.radius = 0.5f;
  p.amplitude = 1.0f;
  p.envelope = 1.0f;
  p.sync();

  // Sample the meridian at colatitude offset d from the center.
  auto at = [&](float d) {
    return bump_field(
        math::Vector(std::sin(C_LAT + d), std::cos(C_LAT + d), 0.0f), p);
  };

  HS_EXPECT_NEAR(at(0.0f), 0.0f, 2e-3f);
  HS_EXPECT_NEAR(at(p.radius * 0.5f), p.radius * 0.5f, 2e-3f);
  HS_EXPECT_NEAR(at(-p.radius * 0.5f), -p.radius * 0.5f, 2e-3f);
  HS_EXPECT_NEAR(at(p.radius * 0.98f), 0.0f, 5e-3f);
  HS_EXPECT_NEAR(at(p.radius + 0.05f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(at(-p.radius - 0.05f), 0.0f, 1e-7f);

  // Half gain halves the drape; a huge gain saturates at the full boundary
  // clearance (R - |y|); zero gain switches the bump off entirely.
  p.amplitude = 0.5f;
  HS_EXPECT_NEAR(at(p.radius * 0.5f), p.radius * 0.25f, 2e-3f);
  p.amplitude = 100.0f;
  const float y_low = p.radius * 0.1f;
  HS_EXPECT_NEAR(at(y_low), p.radius - y_low, 2e-3f);
  p.amplitude = 0.0f;
  HS_EXPECT_NEAR(at(p.radius * 0.5f), 0.0f, 1e-7f);
}

/**
 * @brief Verifies the bulge profile along one ring is the boundary arc scaled
 *        by the ring's drape weight, so the bow is round rather than pointed.
 */
inline void test_bump_field_round_bulge_along_ring() {
  const float C_LAT = math::PI_F / 2.0f;
  Animation::BumpParams p;
  p.center = math::Vector(std::sin(C_LAT), std::cos(C_LAT), 0.0f);
  p.axis = math::Vector(0, 1, 0);
  p.radius = 0.5f;
  p.envelope = 1.0f;
  p.sync();

  // Ring half a footprint below the center (drape weight sin(pi/2) = 1),
  // sampled at azimuth offset 0.3 rad about the axis. There x ~ 0.3 * sin of
  // the ring's colatitude, and the push is sqrt(R^2 - x^2) - |y|.
  const float y = p.radius * 0.5f;
  const float lat = C_LAT + y;
  const float AZ = 0.3f;
  const math::Vector v(std::sin(lat) * std::cos(AZ), std::cos(lat),
                       std::sin(lat) * std::sin(AZ));
  const float x = AZ * std::sin(lat);
  const float expected = std::sqrt(p.radius * p.radius - x * x) - y;
  HS_EXPECT_NEAR(bump_field(v, p), expected, 5e-3f);
}

/** @brief Verifies the cached-ring evaluator matches the generic bump field. */
inline void test_bump_field_precomputed_y_parity() {
  const float C_LAT = 1.1f;
  Animation::BumpParams p;
  p.center = math::Vector(std::sin(C_LAT), std::cos(C_LAT), 0.0f);
  p.axis = math::Vector(0, 1, 0);
  p.radius = 0.55f;
  p.amplitude = 0.8f;
  p.envelope = 0.7f;
  p.sync();

  const float center_colat =
      math::fast_acos(hs::clamp(math::dot(p.axis, p.center), -1.0f, 1.0f));
  for (float lat = 0.4f; lat <= 1.8f; lat += 0.07f) {
    for (float az = -0.6f; az <= 0.6f; az += 0.09f) {
      math::Vector v(std::sin(lat) * std::cos(az), std::cos(lat),
                     std::sin(lat) * std::sin(az));
      const float y = lat - center_colat;
      HS_EXPECT_NEAR(bump_field_with_y(v, p, y), bump_field(v, p), 2e-3f);
    }
  }
}

/** @brief A bump with precomputed y stays finite and bounded near its center. */
inline void test_bump_field_with_y_near_center_stays_bounded() {
  Animation::BumpParams p;
  p.center = p.axis = math::Y_AXIS;
  p.radius = 0.55f;
  p.amplitude = 0.8f;
  p.envelope = 0.7f;
  p.sync();

  for (float angle : {0.0001f, 0.0005f, 0.001f}) {
    const math::Vector point(sinf(angle), cosf(angle), 0.0f);
    const float field = bump_field_with_y(point, p, angle);
    HS_EXPECT_TRUE(std::isfinite(field));
    HS_EXPECT_LE(std::abs(field), p.field_bound());
  }
}

/**
 * @brief Verifies the lifecycle envelope scales the footprint (a half-envelope
 *        cap clears to its shrunken boundary) and gates the field to exactly 0
 *        when fully faded.
 */
inline void test_bump_field_envelope_gates() {
  const float C_LAT = math::PI_F / 3.0f;
  Animation::BumpParams p;
  p.center = math::Vector(std::sin(C_LAT), std::cos(C_LAT), 0.0f);
  p.axis = math::Vector(0, 1, 0);
  p.radius = 0.5f;
  p.sync();

  // Half the shrunken footprint below the center on the axis meridian.
  const float D = C_LAT + 0.125f;
  const math::Vector v(std::sin(D), std::cos(D), 0.0f);
  p.envelope = 0.5f; // effective radius 0.25
  p.sync();
  HS_EXPECT_NEAR(bump_field(v, p), 0.125f, 2e-3f);
  p.envelope = 0.0f;
  p.sync();
  HS_EXPECT_NEAR(bump_field(v, p), 0.0f, 1e-7f);
}

/**
 * @brief Sweeps footprint, envelope, gain, cap centre and stack axis to verify
 *        |bump_field| never exceeds field_bound().
 */
inline void test_bump_field_bound_is_conservative() {
  const hs_test::ScopedRng saved_rng(20260803);
  float worst_ratio = 0.0f;
  int nonfinite = 0;
  int over_bound = 0;
  for (int trial = 0; trial < 96; ++trial) {
    Animation::BumpParams p;
    p.center = rand_unit();
    p.axis = rand_unit();
    p.radius = hs::rand_f(0.0f, 1.2f);
    p.envelope = hs::rand_f(0.0f, 1.0f);
    constexpr float MAX_AMPLITUDE = 6.4f;
    p.amplitude = hs::rand_f(0.0f, MAX_AMPLITUDE);
    p.sync();
    const float bound = p.field_bound();

    for (int i = 0; i < 512; ++i) {
      const float f = bump_field(rand_unit(), p);
      nonfinite += !std::isfinite(f);
      over_bound += std::fabs(f) > bound;
      if (bound > 0.001f)
        worst_ratio = std::max(worst_ratio, std::fabs(f) / bound);
    }
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_EQ(over_bound, 0);
  // The sweep must approach the bound, otherwise it would pass on a field that
  // never leaves zero.
  HS_EXPECT_GT(worst_ratio, 0.5f);
}

// ============================================================================
// noise_product_field — amplitude linearity and separability; short-circuit;
// field_bound stays conservative across the parameter space
// ============================================================================

/**
 * @brief Verifies amplitude linearity, rank-1 separability and the near-zero
 *        amplitude short-circuit.
 */
inline void test_noise_product_field_is_linear_and_separable() {
  Animation::NoiseProductParams p;
  p.amplitude = 0.25f;
  p.scale1 = 1.5f;
  p.scale2 = 3.0f;
  p.time = 2.0f;
  p.noise.SetSeed(1234);

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                                  math::Vector(0.4f, 0.6f, 0.7f).normalized()};
  constexpr size_t NUM_SAMPLES = sizeof(samples) / sizeof(samples[0]);
  float total = 0.0f;
  float base[NUM_SAMPLES] = {};
  for (size_t i = 0; i < NUM_SAMPLES; ++i) {
    const math::Vector &v = samples[i];
    base[i] = noise_product_field(v, p);
    total += std::fabs(base[i]);
  }
  HS_EXPECT_GT(total, 1e-5f);

  // Independent of the octave formula: amplitude scales the field linearly and
  // its sign mirrors it, so amplitude cannot leak into where the noise is
  // sampled.
  p.amplitude = -0.5f;
  for (size_t i = 0; i < NUM_SAMPLES; ++i)
    HS_EXPECT_NEAR(noise_product_field(samples[i], p), -2.0f * base[i], 1e-6f);

  float matrix[2][2]{};
  for (int row = 0; row < 2; ++row)
    for (int col = 0; col < 2; ++col) {
      p.scale1 = row == 0 ? 1.5f : 2.75f;
      p.scale2 = col == 0 ? 3.0f : 4.5f;
      matrix[row][col] = noise_product_field(samples[2], p);
    }
  HS_EXPECT_GT(std::fabs(matrix[0][0] - matrix[1][0]), 1e-5f);
  HS_EXPECT_GT(std::fabs(matrix[0][0] - matrix[0][1]), 1e-5f);
  HS_EXPECT_NEAR(matrix[0][0] * matrix[1][1], matrix[0][1] * matrix[1][0],
                 1e-7f);

  p.amplitude = .0005f;
  HS_EXPECT_EQ(noise_product_field(samples[0], p), 0.f);
  HS_EXPECT_EQ(p.field_bound(), .0005f);
}

/**
 * @brief Sweeps seeds, octave scales, field time, amplitude, and sample points
 *        to verify |noise_product_field| never exceeds field_bound().
 * @details Outputs within [-1, 1] on both octaves guarantee
 * field_bound() = |amplitude|; the sweep checks both octave ranges and the
 * product bound. Near-zero point draws are skipped.
 */
inline void test_noise_product_field_bound_is_conservative() {
  const hs_test::ScopedRng saved_rng(20260801);
  float worst_ratio = 0.0f;
  float worst_octave = 0.0f;
  int nonfinite = 0;
  int over_bound = 0;
  int octave_out_of_range = 0;
  for (int trial = 0; trial < 96; ++trial) {
    const float amplitude = hs::rand_f(-2.0f, 2.0f);
    const float scale1 = hs::rand_f(0.25f, 12.0f);
    const float scale2 = hs::rand_f(0.25f, 12.0f);
    const float time = hs::rand_f(-200.0f, 200.0f);
    const int seed = hs::rand_int(0, 100000);

    Animation::NoiseProductParams p;
    p.amplitude = amplitude;
    p.scale1 = scale1;
    p.scale2 = scale2;
    p.time = time;
    p.noise.SetSeed(seed);
    const float bound = p.field_bound();

    for (int i = 0; i < 512; ++i) {
      const float rx = hs::rand_f(-1.0f, 1.0f);
      const float ry = hs::rand_f(-1.0f, 1.0f);
      const float rz = hs::rand_f(-1.0f, 1.0f);
      const math::Vector r(rx, ry, rz);
      if (r.length() < 0.1f)
        continue;
      const math::Vector v = r.normalized();

      const float f = noise_product_field(v, p);
      nonfinite += !std::isfinite(f);
      over_bound += std::fabs(f) > bound;

      const float n1 = p.octave1(v, p.time);
      const float n2 = p.octave2(v, p.time);
      octave_out_of_range += std::fabs(n1) > 1.0f;
      octave_out_of_range += std::fabs(n2) > 1.0f;
      worst_octave = std::max(worst_octave, std::fabs(n1));
      worst_octave = std::max(worst_octave, std::fabs(n2));
      if (bound > 0.001f)
        worst_ratio = std::max(worst_ratio, std::fabs(f) / bound);
    }
  }
  HS_EXPECT_EQ(nonfinite, 0);
  HS_EXPECT_EQ(over_bound, 0);
  HS_EXPECT_EQ(octave_out_of_range, 0);
  // The sweep must reach a substantial fraction of the bound, otherwise it
  // would pass on a field that never leaves zero.
  HS_EXPECT_GT(worst_octave, 0.5f);
  HS_EXPECT_GT(worst_ratio, 0.1f);
}

// ============================================================================
// Animation::BallDrop — pole-to-pole traversal, envelope, slot reclaim
// ============================================================================

/**
 * @brief Verifies a spawned BallDrop traverses its bump from the frame's pole
 *        to the opposite pole along the requested meridian, with an envelope
 *        that rises monotonically to full strength at mid-fall and falls back
 *        to zero, and frees its pool slot on completion.
 */
inline void test_ball_drop_traverses_and_reclaims() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  math::Orientation<> ori;
  const math::Vector POLE(0, 1, 0);
  const int DURATION = 40;

  BallDropTransformer<1> balls(tl);
  balls.init_storage(persistent_arena);
  balls.template_params.radius = 0.5f;
  HS_EXPECT_TRUE(balls.spawn(0, ori, POLE, 0.0f, DURATION) != nullptr);

  tl.step(cv);
  balls.prepare_frame();
  const float early_bound = balls.field_bound();
  HS_EXPECT_GT(early_bound, 0.0f);
  // The full-envelope plateau reads with rounding jitter.
  constexpr float ULP_SLACK = 1e-5f;
  float previous = early_bound;
  for (int i = 1; i < DURATION / 2; ++i) {
    tl.step(cv);
    balls.prepare_frame();
    HS_EXPECT_GE(balls.field_bound(), previous - ULP_SLACK);
    previous = balls.field_bound();
  }
  // Mid-fall: the bump sits at (1, 0, 0) at full envelope. Half a footprint
  // below the center the push is R - R/2; the mirrored point reads its negation.
  HS_EXPECT_NEAR(balls.field(POLE), 0.0f, 1e-5f);
  const float HALF_R = 0.5f * 0.5f;
  const math::Vector below(std::cos(HALF_R), -std::sin(HALF_R), 0.0f);
  const math::Vector above(std::cos(HALF_R), std::sin(HALF_R), 0.0f);
  HS_EXPECT_NEAR(balls.field(below), HALF_R, 3e-3f);
  HS_EXPECT_NEAR(balls.field(above), -HALF_R, 3e-3f);
  const float mid_bound = balls.field_bound();
  HS_EXPECT_GE(mid_bound, balls.field(below));
  HS_EXPECT_GT(mid_bound, early_bound);

  for (int i = DURATION / 2; i < 3 * DURATION / 4; ++i) {
    tl.step(cv);
    balls.prepare_frame();
    HS_EXPECT_LE(balls.field_bound(), previous + ULP_SLACK);
    previous = balls.field_bound();
  }
  // Three quarters down, the bump has left the mid-fall footprint and pushes
  // half a footprint below its new center.
  const float LATE_PHI = 3.0f * math::PI_F / 4.0f + HALF_R;
  const math::Vector late(std::sin(LATE_PHI), std::cos(LATE_PHI), 0.0f);
  HS_EXPECT_NEAR(balls.field(below), 0.0f, 1e-5f);
  HS_EXPECT_GT(balls.field(late), 0.0f);
  for (int i = 3 * DURATION / 4; i < DURATION - 1; ++i) {
    tl.step(cv);
    balls.prepare_frame();
    HS_EXPECT_LE(balls.field_bound(), previous + ULP_SLACK);
    previous = balls.field_bound();
  }
  HS_EXPECT_LT(previous, mid_bound);
  for (int i = DURATION - 1; i <= DURATION + 2; ++i)
    tl.step(cv);
  balls.prepare_frame();
  HS_EXPECT_EQ(balls.field_bound(), 0.0f);
  HS_EXPECT_EQ(balls.active_count(), 0);
  HS_EXPECT_TRUE(balls.spawn(0, ori, POLE, 1.0f, DURATION) != nullptr);
}

// ============================================================================
// Animation::NoiseProduct — initial linear trajectory
// ============================================================================

/**
 * @brief NoiseProduct publishes its initial linear trajectory.
 */
inline void test_noise_product_integrates_time() {
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);
  Animation::NoiseProductParams p;
  p.speed = 0.03f;
  Animation::NoiseProduct anim(p);
  for (int i = 0; i < 5; ++i)
    anim.step(cv);
  HS_EXPECT_NEAR(p.time, 5 * 0.03f, 1e-5f);
}

} // namespace hs_test::transformers_tests
