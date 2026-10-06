/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/render/pullback/composed_effect.h"

struct EdgeFadeSpec : Pullback::Spec {
#ifdef HS_EDGE_FADE_NO_DISTANCE
  static constexpr auto PROJECTION =
      Pullback::ProjectionKind::FOLDED_SINUSOIDAL;
#endif
#ifdef HS_EDGE_FADE_DISABLED
  static constexpr auto COVERAGE = Pullback::ProjectionCoverageMode::NONE;
#else
  static constexpr auto COVERAGE = Pullback::ProjectionCoverageMode::EDGE_FADE;
#endif
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Project<
          typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              COVERAGE, B, Pullback::EdgeValueParams>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
struct EdgeFadeEffect
    : Pullback::ComposedEffect<32, 16, EdgeFadeEffect, EdgeFadeSpec> {};
static_assert(sizeof(EdgeFadeEffect) > 0);
