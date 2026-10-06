/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// Motion repetition, reanchoring and co-driving.

/** @brief Internal-angle allowance for float drift after 600 cycles, radians. */
constexpr float MOTION_WARP_TOL = 1e-4f;

/**
 * @brief Verifies a repeating Motion does not drift across many cycles.
 * @details A repeating Motion advances its Orientation by relative deltas taken
 * between consecutive path frames, where each frame is a pure function of the
 * path parameter (point + tangent). Because the frame depends only on the phase,
 * the per-cycle product of deltas telescopes — there is no accumulating
 * quaternion chain to warp the traced curve. The decisive, precession-immune
 * signature is the set of rotation-INVARIANT internal angles between heads
 * sampled at fixed phases within a cycle: a rigid drift (holonomy) leaves them
 * unchanged, so any growth is genuine warp. A late cycle is compared against the
 * ideal Lissajous internal angles within accumulated float drift.
 */
inline void test_motion_repeating_does_not_drift() {
  using Ori = math::Orientation<16>;
  constexpr int duration = 40;
  ProceduralPath path;
  // lissajous(.,.,.,0) == +Y, so the identity-start orientation places the head
  // on the path at phase 0.
  path.f = [](float t) {
    return math::lissajous(1.06f, 1.06f, 0.0f, t * 5.909f);
  };

  Ori o; // identity, single frame; orient(+Y) starts on the path
  const math::Vector node_v = math::Y_AXIS;

  Timeline tl;
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));

  math::Vector late_heads[duration + 1]; // indexed by phase 1..duration
  const int cycles = 600;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      if (c == cycles - 1)
        late_heads[fr] = o.orient(node_v);
    }
  }

  // Each phase's internal angle to the phase-1 anchor must match the ideal
  // Lissajous internal angle; a rigid precession leaves these untouched, so any
  // growth is genuine warp. Interior phases only (the boundary frame rewinds).
  const int anchor = 1;
  const math::Vector ideal_anchor = path.f((float)anchor / duration);
  float max_error = 0.0f;
  for (int fr = 2; fr < duration; ++fr) {
    const math::Vector ideal_fr = path.f((float)fr / duration);
    const float ERROR =
        fabsf(math::angle_between(late_heads[anchor], late_heads[fr]) -
              math::angle_between(ideal_anchor, ideal_fr));
    max_error = hs_test::fold_worst(max_error, ERROR);
  }
  std::printf("  late-cycle motion angle drift: %.9g rad\n", max_error);
  HS_EXPECT_LT(max_error, MOTION_WARP_TOL);
}

/** @brief Verifies reanchor preserves orientation continuity when the live path changes. */
inline void test_motion_reanchor_after_path_swap() {
  const auto step_after_swap = [](bool reanchor) {
    ProceduralPath path;
    path.f = [](float t) {
      const float ANGLE = math::TWO_PI_F * t;
      return math::Vector(cosf(ANGLE), sinf(ANGLE), 0);
    };
    math::Orientation<16> orientation;
    Animation::Motion<288, 16> motion(orientation, path, 1000, true);
    for (int frame = 0; frame < 5; ++frame)
      motion.step(fake_canvas());
    const math::Vector BEFORE = orientation.orient(math::X_AXIS);
    path.f = [](float t) {
      const float ANGLE = math::TWO_PI_F * t + math::PI_F * 0.5f;
      return math::Vector(cosf(ANGLE), sinf(ANGLE), 0);
    };
    if (reanchor)
      motion.reanchor();
    motion.step(fake_canvas());
    return math::angle_between(BEFORE, orientation.orient(math::X_AXIS));
  };
  HS_EXPECT_LT(step_after_swap(true), 0.02f);
  HS_EXPECT_GT(step_after_swap(false), 1.0f);
}

/**
 * @brief A co-driver sharing a repeating Motion's Orientation survives the
 * repeat seam.
 * @details Motion re-seats via a relative delta; the co-driver's accumulated
 * rotation persists across the seam. With a CLOSED path Motion's per-cycle
 * contribution telescopes to identity, so the only thing that should move the
 * shared orientation at a seam is the co-driver's own small step — never a
 * large snap-back. Assert the probe's per-frame angular step stays bounded
 * across many seams while its cumulative travel is large (so the co-driver is
 * provably active, not a no-op).
 */
inline void test_motion_codriven_survives_repeat_seam() {
  using Ori = math::Orientation<16>;
  const int duration = 30;
  ProceduralPath path;
  // A closed great circle: path(0) == path(1) with matching tangent, so Motion's
  // per-cycle delta product is identity.
  path.f = [](float t) {
    float a = 2.0f * math::PI_F * t;
    return math::Vector(std::cos(a), std::sin(a), 0.0f);
  };

  Ori o; // identity
  Timeline tl;
  // Repeating Motion + a repeating co-driver rotation about Y, both driving `o`.
  tl.add(0, Animation::Motion<288, 16>(o, path, duration, /*repeat=*/true));
  tl.add(0, Animation::Rotation<288, 16>(o, math::Y_AXIS, 2.0f * math::PI_F,
                                         duration, math::ease_linear,
                                         /*repeat=*/true));

  const math::Vector probe = math::Z_AXIS;
  math::Vector prev = o.orient(probe);
  float max_step = 0.0f;
  float total_travel = 0.0f;
  const int cycles = 8;
  for (int c = 0; c < cycles; ++c) {
    for (int fr = 1; fr <= duration; ++fr) {
      tl.step(fake_canvas());
      math::Vector cur = o.orient(probe);
      float step = math::angle_between(prev, cur);
      max_step = std::max(max_step, step);
      total_travel += step;
      prev = cur;
    }
  }

  // No single frame (seams included) snaps the orientation; 0.8 clears the
  // largest legitimate one-frame step but is far below a multi-radian snap.
  HS_EXPECT_LT(max_step, 0.8f);
  HS_EXPECT_GT(total_travel, 5.0f);
}
