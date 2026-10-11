/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/render/pullback/composed_effect.h"

struct UnknownCoverage {};
struct MetadataSpec : Pullback::Spec {
#ifdef HS_METADATA_REDUNDANT_SCALAR
  static constexpr auto HUE = Pullback::HueMode::NOISE;
#endif
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, void,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::Glitch>,
          Pullback::Stage::Placed<
              Pullback::CodeEmission::OUT_OF_LINE_FLASH,
              Pullback::Stage::Displace<
                  Pullback::Surface::PeriodicRipple<Pullback::SurfaceProvider<
                      B, Pullback::PeriodicRippleParams, false>>>,
#ifdef HS_METADATA_MISSING_PROJECT
              void
#else
              Pullback::Stage::Project<Pullback::Projection::Equirectangular<
                  Pullback::ProjectionProvider<B>>>
#endif
              >>,
      Pullback::Stage::Sample<Pullback::Source::Grid<Pullback::SourceProvider<
                                  B, Pullback::GridSourceParams>>,
                              Pullback::Weight::Projection,
#ifdef HS_METADATA_UNKNOWN_COVERAGE
                              UnknownCoverage>,
#else
                              Pullback::ProjectionCoverage::WeightSquared>,
#endif
      Pullback::Stage::Transfer<Pullback::Transfer::IsoContour<
          Pullback::ValueProvider<B, Pullback::IsoValueParams>>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B,
#ifdef HS_METADATA_UNTRACKED_HUE
                                  Pullback::HueMode::PATH_LENGTH,
#else
                                  Pullback::HueMode::NOISE,
#endif
                                  Pullback::Color::BrightnessEnvelope::BELL>>>>;
};
using Metadata = Pullback::MetadataFor<MetadataSpec>;
static_assert(sizeof(Metadata) > 0);
#if !defined(HS_METADATA_UNKNOWN_COVERAGE) &&                                  \
    !defined(HS_METADATA_MISSING_PROJECT)
static_assert(Metadata::PROJECTION ==
              Pullback::ProjectionKind::EQUIRECTANGULAR);
static_assert(Metadata::COVERAGE ==
              Pullback::ProjectionCoverageMode::WEIGHT_SQUARED);
static_assert(Metadata::TRANSFER == Pullback::TransferKind::ISO_CONTOUR);
static_assert(Metadata::FIELD_COVERAGE == Pullback::FieldCoverageKind::NONE);
static_assert(Metadata::SURFACE_PLACEMENT ==
              Pullback::SurfacePlacement::AFTER_LENS);
static_assert(std::is_same_v<Metadata::LensPolicy, Pullback::Lens::Glitch>);
#ifndef HS_METADATA_UNTRACKED_HUE
static_assert(Metadata::HUE == Pullback::HueMode::NOISE);
#endif
static_assert(Metadata::BRIGHTNESS ==
              Pullback::Color::BrightnessEnvelope::BELL);
static_assert(!Metadata::PATH_TRACKED);
struct MetadataEffect
    : Pullback::ComposedEffect<8, 8, MetadataEffect, MetadataSpec,
                               Control::HookOwnership::INHERITED> {};
static_assert(sizeof(MetadataEffect) > 0);
#endif
