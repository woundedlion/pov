/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Deforming-lens displacement and fold invariance.

template <typename Policy, bool Track = true>
inline PB::SphereSample lens_path_stage(const PB::SphereSample &input,
                                        const SphereOpMirrorFrame &frame) {
  using Bound =
      typename PB::Stage::Lens<Policy,
                               Track>::template Bind<SphereOpMirrorBinding>;
  return Bound::run(input, frame, Bound::prepare(frame));
}

/** @brief Carrier displacement adds the normalized shortest angular distance. */
inline void test_lens_displacement_metric() {
  const PB::SphereSample input{1.0003f * math::X_AXIS, 0.375f};
  const auto same = PB::Kernel::lens(input, input.dir, true);
  const auto quarter = PB::Kernel::lens(input, 0.9997f * math::Y_AXIS, true);
  const auto opposite = PB::Kernel::lens(input, -input.dir, true);
  HS_EXPECT_NEAR(same.path_length, 0.375f, 1e-4f);
  HS_EXPECT_NEAR(quarter.path_length, 0.375 + acos(-1.0) / 2, 1e-4f);
  HS_EXPECT_NEAR(opposite.path_length, 0.375 + acos(-1.0), 1e-4f);
  const auto oblique =
      PB::Kernel::lens(input, 1.0004f * math::Vector(0.6f, 0.8f, 0), true);
  HS_EXPECT_NEAR(oblique.path_length, 0.375 + acos(0.6), 1e-4f);
  HS_EXPECT_EQ(PB::Kernel::lens(input, math::Y_AXIS).path_length, 0.375f);
}

/** @brief Static and animated twists contribute their geometric displacement. */
inline void test_twist_displacement_paths() {
  using Animated = PB::Lens::AnimatedTwist<TwistMirrorProvider>;
  const PB::SphereSample input{math::Vector(sqrtf(0.75f), 0.5f, 0), 0.125f};
  const auto ctx = shared_resources().context();
  SphereOpMirrorFrame frame;
  frame.twist.twist_angle = 4;
  frame.twist.distribution = 0;
  const double expected = 0.125 + acos(0.25 + 0.75 * cos(1.0));
  const auto typed = lens_path_stage<Animated>(input, frame);
  const auto prepared = In::Op::LensTwist::prepare(ctx, frame.twist, {});
  const auto native = In::Op::LensTwist::run(input, ctx, frame.twist, prepared);
  HS_EXPECT_NEAR(typed.path_length, expected, 2e-3f);
  HS_EXPECT_NEAR(native.path_length, expected, 2e-3f);
  const double a = lenses::TWIST_DISTRIBUTION * 0.999;
  const double angle =
      0.5 * lenses::TWIST_ANGLE * (a == 0 ? 0.5 : atanh(a * 0.5) / atanh(a));
  const auto fixed = lens_path_stage<PB::Lens::Twist>(input, frame);
  HS_EXPECT_NEAR(fixed.path_length, 0.125 + acos(0.25 + 0.75 * cos(angle)),
                 3e-3f);
  HS_EXPECT_EQ((lens_path_stage<Animated, false>(input, frame).path_length),
               input.path_length);
  const auto pole = lens_path_stage<Animated>({math::Y_AXIS, 0.25f}, frame);
  HS_EXPECT_NEAR(pole.path_length, 0.25f, 1e-4f);
  frame.twist.twist_angle = 0;
  HS_EXPECT_NEAR(lens_path_stage<Animated>(input, frame).path_length,
                 input.path_length, 1e-3f);
}

/** @brief Mobius spin and push contribute independently known angular distances. */
inline void test_mobius_displacement_paths() {
  using Policy = PB::Lens::MobiusGenerator<MobiusMirrorProvider>;
  const PB::SphereSample input{math::X_AXIS, 0.125f};
  const auto ctx = shared_resources().context();
  SphereOpMirrorFrame frame;
  const auto native = [&](const PB::SphereSample &p) {
    return In::Op::LensMobius::run(
        p, ctx, frame.mobius,
        In::Op::LensMobius::prepare(ctx, frame.mobius, {}));
  };
  for (const float spin :
       {0.0f, 0.5f * math::PI_F, math::PI_F, math::TWO_PI_F}) {
    frame.mobius.spin_y = spin;
    const double expected = 0.125 + acos(cos(double(spin)));
    HS_EXPECT_NEAR(native(input).path_length, expected, 1e-3f);
    HS_EXPECT_NEAR(lens_path_stage<Policy>(input, frame).path_length, expected,
                   1e-3f);
  }
  frame.mobius.spin_y = 0;
  frame.mobius.push_y = 1;
  const double expected = 0.125 + atan(sinh(1.0));
  HS_EXPECT_NEAR(native(input).path_length, expected, 1e-4f);
  HS_EXPECT_NEAR(lens_path_stage<Policy>(input, frame).path_length, expected,
                 1e-4f);
  HS_EXPECT_NEAR(native({math::Y_AXIS, 0.25f}).path_length, 0.25f, 1e-4f);
  HS_EXPECT_EQ((lens_path_stage<Policy, false>(input, frame).path_length),
               input.path_length);
}

/** @brief Fold policies, native folds and the camera preserve existing displacement. */
inline void test_fold_displacement_invariance() {
  const PB::SphereSample input{math::Vector(0.6f, 0.5f, sqrtf(0.39f)), 0.375f};
  const auto ctx = shared_resources().context();
  const SphereOpMirrorFrame frame;
  const auto check = [&]<typename Policy>() {
    HS_EXPECT_EQ(lens_path_stage<Policy>(input, frame).path_length,
                 input.path_length);
    static_assert(
        !PB::ComposedDetail::StagePathTracked<PB::Stage::Lens<Policy>>::value);
  };
  check.template operator()<PB::Lens::Kaleidoscope>();
  check.template operator()<PB::Lens::TetrahedralKaleidoscope>();
  check.template operator()<PB::Lens::OctahedralKaleidoscope>();
  check.template operator()<PB::Lens::DodecahedralKaleidoscope>();
  check.template operator()<PB::Lens::TriangularPrismKaleidoscope>();
  check.template operator()<PB::Lens::SquarePrismKaleidoscope>();
  check.template operator()<PB::Lens::PentagonalPrismKaleidoscope>();
  check.template operator()<PB::Lens::HexagonalPrismKaleidoscope>();
  check.template operator()<PB::Lens::OctagonalPrismKaleidoscope>();
  check.template operator()<PB::Lens::Glitch>();
  for (uint8_t symmetry = 0; symmetry < 9; ++symmetry) {
    In::Op::KaleidoscopeChainParams params;
    params.symmetry = symmetry;
    const auto output = In::Op::LensKaleidoscope::run(input, ctx, params, {});
    HS_EXPECT_EQ(output.path_length, input.path_length);
  }
  HS_EXPECT_EQ(In::Op::LensGlitch::run(input, ctx, {}, {}).path_length,
               input.path_length);
  const auto camera =
      PB::Kernel::rotate_dir(input, math::make_rotation(math::X_AXIS, 0.7f));
  HS_EXPECT_EQ(camera.path_length, input.path_length);
}

/** @brief Each deforming lens adds its displacement across intervening folds. */
inline void test_lens_displacement_composition() {
  using Policy = PB::Lens::MobiusGenerator<MobiusMirrorProvider>;
  SphereOpMirrorFrame frame;
  frame.mobius.spin_y = 0.4f;
  const auto first = lens_path_stage<Policy>({math::X_AXIS, 0.25f}, frame);
  const auto folded = lens_path_stage<PB::Lens::Kaleidoscope>(first, frame);
  const auto second = lens_path_stage<Policy>(folded, frame);
  HS_EXPECT_NEAR(first.path_length, 0.65f, 1e-4f);
  HS_EXPECT_EQ(folded.path_length, first.path_length);
  HS_EXPECT_NEAR(second.path_length, 1.05f, 2e-4f);
}

/** @brief Lens displacement drives visible downstream opacity without surface steps. */
inline void test_lens_displacement_modulation() {
  const auto ctx = shared_resources().context();
  In::Op::ModulationParams params;
  params.modulation.source = Pattern::SourceKind::DISPLACEMENT;
  params.modulation.evaluation = Pattern::Evaluation::DIRECT;
  params.modulation.mechanic = {ColorModulation::Mechanic::OPACITY, 1};
  const auto prepared = In::Op::ModulateColorV1::prepare(ctx, params, {});
  const auto alpha = [&](const PB::SphereSample &sphere) {
    const PB::ColorSample color{Color4(Pixel(30000, 20000, 10000), 1),
                                PB::Kernel::sample(sphere, 0)};
    return In::Op::ModulateColorV1::run(color, ctx, params, prepared)
        .color.alpha;
  };
  SphereOpMirrorFrame frame;
  const PB::SphereSample input{math::Vector(sqrtf(0.75f), 0.5f, 0), 0};
  HS_EXPECT_EQ(alpha(input), 0.5f);
  frame.twist.twist_angle = 4;
  frame.twist.distribution = 0;
  const auto twisted =
      In::Op::LensTwist::run(input, ctx, frame.twist,
                             In::Op::LensTwist::prepare(ctx, frame.twist, {}));
  const double twist_arc = acos(0.25 + 0.75 * cos(1.0));
  HS_EXPECT_NEAR(alpha(twisted), 0.5 * (1 + twist_arc), 1e-3f);
  HS_EXPECT_GT(alpha(twisted), alpha(input) + 0.2f);
  frame.mobius.spin_y = 0.5f;
  const PB::SphereSample equator{math::X_AXIS, 0};
  const auto mobius = In::Op::LensMobius::run(
      equator, ctx, frame.mobius,
      In::Op::LensMobius::prepare(ctx, frame.mobius, {}));
  HS_EXPECT_NEAR(alpha(mobius), 0.75f, 1e-4f);
  for (uint8_t symmetry = 0; symmetry < 9; ++symmetry) {
    In::Op::KaleidoscopeChainParams fold;
    fold.symmetry = symmetry;
    HS_EXPECT_EQ(alpha(In::Op::LensKaleidoscope::run(input, ctx, fold, {})),
                 alpha(input));
  }
  HS_EXPECT_EQ(alpha(In::Op::LensGlitch::run(input, ctx, {}, {})),
               alpha(input));
}

/** @brief Typed metadata recognizes enabled lens tracking and explicit opt-outs. */
inline void test_lens_displacement_metadata() {
  using Policy = PB::Lens::MobiusGenerator<MobiusMirrorProvider>;
  static_assert(
      PB::ComposedDetail::StagePathTracked<PB::Stage::Lens<Policy>>::value);
  static_assert(!PB::ComposedDetail::StagePathTracked<
                PB::Stage::Lens<Policy, false>>::value);
  static_assert(
      PB::ComposedDetail::IsLensStage<PB::Stage::Lens<Policy, false>>::value);
  static_assert(!PB::ComposedDetail::StagePathTracked<
                PB::Stage::Lens<PB::Lens::Kaleidoscope>>::value);
  HS_EXPECT_TRUE(PB::ComposedDetail::MetadataFor<MobiusGridSpec>::PATH_TRACKED);
  HS_EXPECT_FALSE(PB::ComposedDetail::MetadataFor<GridWormsSpec>::PATH_TRACKED);
}
