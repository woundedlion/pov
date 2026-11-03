/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/render/pullback/composed_effect.h"

template <typename B>
using EdgeFadeProjection =
#ifdef HS_EDGE_FADE_AUTHORED
    Pullback::AuthoredProjection<
#ifdef HS_EDGE_FADE_NO_DISTANCE
        Pullback::ProjectionKind::FOLDED_SINUSOIDAL,
#else
        Pullback::ProjectionKind::STEREOGRAPHIC,
#endif
        Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>;
#elif defined(HS_EDGE_FADE_NO_DISTANCE)
    Pullback::Projection::FoldedSinusoidal<Pullback::ProjectionProvider<B>>;
#else
    Pullback::Projection::Stereographic<Pullback::ProjectionProvider<B>>;
#endif

struct EdgeFadeSpec : Pullback::Spec {
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Project<EdgeFadeProjection<B>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::GridSourceParams, B>::Type,
                              Pullback::Weight::Projection,
#ifdef HS_EDGE_FADE_DISABLED
                              Pullback::ProjectionCoverage::None>,
#else
                              Pullback::ProjectionCoverage::EdgeFade<
                                  Pullback::ValueProvider<
                                      B, Pullback::EdgeValueParams>>>,
#endif
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
struct EdgeFadeEffect
    : Pullback::ComposedEffect<32, 16, EdgeFadeEffect, EdgeFadeSpec> {
public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<32, 16, EdgeFadeEffect,
                                    EdgeFadeSpec>::choreography_hooks()
        .template extend<EdgeFadeEffect>();
  }
};
static_assert(sizeof(EdgeFadeEffect) > 0);
