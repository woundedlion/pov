/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>

#include "core/render/ray/march.h"
#include "core/render/scan/shader.h"
#include "core/render/sdf/periodic_surface.h"

#ifndef HS_SPHERICAL_PROFILE_SAMPLES
#define HS_SPHERICAL_PROFILE_SAMPLES 1
#endif
#ifndef HS_SPHERICAL_PROFILE_GRADIENT
#define HS_SPHERICAL_PROFILE_GRADIENT 0
#endif

/** @brief Experimental periodic fields under the segmented profile driver. */
template <int W, int H> class SphericalExperiment : public Effect {
public:
  SphericalExperiment() : Effect(W, H) {}

  void init() override {
    hs::log("spherical experiment: samples=%u gradient=%u positions=36 cases=8",
            HS_SPHERICAL_PROFILE_SAMPLES, HS_SPHERICAL_PROFILE_GRADIENT);
  }

  void draw_frame() override {
    const unsigned CASE = (frame / 36) % 8;
    const unsigned POSITION = frame % 36;
    if (POSITION == 0)
      hs::log("Preset: %u/8", CASE + 1);
    constexpr float PERIODS[] = {0.7f, 1.5f, 3.0f};
    constexpr float ISOS[] = {-0.5f, 0.0f, 0.5f};
    const math::Vector CENTERS[] = {{0.0f, 0.0f, 0.0f},
                                    {0.23f, 0.41f, -0.17f},
                                    {0.5f, 0.5f, 0.5f},
                                    {-0.37f, 0.19f, 0.73f}};
    const float PERIOD = PERIODS[POSITION / 12];
    const float ISO = ISOS[(POSITION / 4) % 3];
    const math::Vector CENTER = CENTERS[POSITION % 4] * PERIOD;
    constexpr int BUDGETS[] = {8, 12, 16, 24};
    Raycast::TraceLimits limits;
    limits.max_queries = BUDGETS[CASE % 4];
    limits.max_steps = 1024;
    limits.max_refinements = 1024;
    limits.position_tolerance = PERIOD * 1e-4f;
    Canvas canvas(*this);
    samples = queries = peak_queries = surfaces = unresolved = exhausted =
        invalid = 0;
    {
      HS_PROFILE(spherical_draw);
      if (CASE < 4)
        draw(canvas, SDF::CosineSurface{PERIOD, ISO, {}}, CENTER, limits);
      else
        draw(canvas, SDF::GyroidSurface{PERIOD, ISO, {}}, CENTER, limits);
    }
    hs::log(
        "spherical counts: frame=%u case=%u position=%u samples=%lu "
        "queries=%lu peak=%lu surface=%lu unresolved=%lu exhausted=%lu invalid=%lu",
        frame + 1, CASE, POSITION, samples, queries, peak_queries, surfaces,
        unresolved, exhausted, invalid);
    ++frame;
  }

private:
  unsigned frame = 0;
  unsigned long samples = 0, queries = 0, peak_queries = 0, surfaces = 0;
  unsigned long unresolved = 0, exhausted = 0, invalid = 0;

  template <typename Surface>
  void draw(Canvas &canvas, const Surface &surface, const math::Vector &center,
            const Raycast::TraceLimits &limits) {
    Scan::Shader::draw_cached<W, H, HS_SPHERICAL_PROFILE_SAMPLES>(
        canvas, [&](const math::Vector &direction) HS_HOT_FLASH_MEMBER {
          const Raycast::Ray RAY{
              center, direction, {0.0f, 4.0f * surface.period}};
          const auto RESULT = Raycast::surface_search(surface, RAY, {}, limits);
          ++samples;
          queries += RESULT.counters.queries;
          peak_queries =
              std::max(peak_queries,
                       static_cast<unsigned long>(RESULT.counters.queries));
          surfaces += RESULT.has_surface;
          unresolved += RESULT.status == Raycast::TraceStatus::UNRESOLVED;
          exhausted += RESULT.status == Raycast::TraceStatus::BUDGET_EXHAUSTED;
          invalid += RESULT.status == Raycast::TraceStatus::INVALID_QUERY;
          if (!RESULT.has_surface)
            return Color4{};
          float value = 1.0f - RESULT.contribution.t / RAY.interval.far;
          if constexpr (HS_SPHERICAL_PROFILE_GRADIENT != 0) {
            const auto NORMAL = surface.normal(RAY.at(RESULT.contribution.t));
            value *= 0.2f + 0.8f * fabsf(math::dot(NORMAL, direction));
          }
          const uint16_t LEVEL =
              static_cast<uint16_t>(hs::clamp(value, 0.0f, 1.0f) * 65535.0f);
          return Color4(Pixel(LEVEL, LEVEL, LEVEL), 1.0f);
        });
  }
};
