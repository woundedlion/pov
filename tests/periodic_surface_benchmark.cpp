/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include <array>
#include <chrono>
#include <cstdio>
#include "core/math/display_geometry.h"
#include "core/render/ray/march.h"
#include "core/render/sdf/periodic_surface.h"

namespace {

constexpr int WIDTH = 288;
constexpr int HEIGHT = 144;
constexpr int LIVE_COLUMNS = WIDTH / 2 + 2;
constexpr int LIVE_ROWS = HEIGHT / 2 + 1;
constexpr int SAMPLES = LIVE_COLUMNS * LIVE_ROWS;
constexpr std::array<int, 5> BUDGETS = {8, 12, 16, 24, 1024};

struct Metrics {
  uint64_t rays = 0;
  uint64_t queries = 0;
  uint64_t unresolved = 0;
  uint64_t reference_unresolved = 0;
  uint64_t hit_disagreement = 0;
  uint64_t inside_starts = 0;
  int peak_queries = 0;
  double image_error = 0.0;
  double total_ms = 0.0;
  double peak_ms = 0.0;
};

float depth_pixel(const Raycast::TraceResult &result, float far) {
  return result.has_surface ? 1.0f - result.contribution.t / far : 0.0f;
}

bool unresolved(const Raycast::TraceResult &result) {
  return result.status != Raycast::TraceStatus::SURFACE &&
         result.status != Raycast::TraceStatus::RANGE_COMPLETE;
}

template <typename Surface> void measure(const char *name) {
  std::array<math::Vector, SAMPLES> directions;
  std::array<Raycast::TraceResult, SAMPLES> reference;
  std::array<Metrics, BUDGETS.size()> metrics{};
  constexpr std::array<float, 3> PERIODS = {0.7f, 1.5f, 3.0f};
  constexpr std::array<float, 3> ISOVALUES = {-0.5f, 0.0f, 0.5f};
  constexpr std::array<math::Vector, 4> CENTERS = {
      math::Vector{0.0f, 0.0f, 0.0f}, math::Vector{0.23f, 0.41f, -0.17f},
      math::Vector{0.5f, 0.5f, 0.5f}, math::Vector{-0.37f, 0.19f, 0.73f}};
  for (int side = 0; side < 2; ++side) {
    for (int y = 0; y < LIVE_ROWS; ++y)
      for (int x = 0; x < LIVE_COLUMNS; ++x)
        directions[y * LIVE_COLUMNS + x] = math::Vector::from_spherical(
            math::TWO_PI_F * ((x - 1 + side * WIDTH / 2 + WIDTH) % WIDTH) /
                WIDTH,
            math::DisplayGeometry<HEIGHT>::row_to_phi(static_cast<float>(y)));
    for (float period : PERIODS) {
      for (float iso : ISOVALUES) {
        for (const math::Vector &center : CENTERS) {
          Surface surface;
          surface.period = period;
          surface.iso = iso;
          const math::Vector CENTER = center * period;
          const float FAR = 4.0f * period;
          Raycast::TraceLimits limits;
          limits.max_queries = 1024;
          limits.max_steps = 1024;
          limits.max_refinements = 1024;
          limits.position_tolerance = period * 1e-4f;
          for (int i = 0; i < SAMPLES; ++i)
            reference[i] = Raycast::surface_search(
                surface, {CENTER, directions[i], {0.0f, FAR}}, {}, limits);
          for (size_t budget_index = 0; budget_index < BUDGETS.size();
               ++budget_index) {
            limits.max_queries = BUDGETS[budget_index];
            Metrics &m = metrics[budget_index];
            const auto START = std::chrono::steady_clock::now();
            for (int i = 0; i < SAMPLES; ++i) {
              const auto RESULT = Raycast::surface_search(
                  surface, {CENTER, directions[i], {0.0f, FAR}}, {}, limits);
              ++m.rays;
              m.queries += RESULT.counters.queries;
              m.peak_queries =
                  std::max(m.peak_queries, RESULT.counters.queries);
              m.unresolved += unresolved(RESULT);
              m.reference_unresolved += unresolved(reference[i]);
              m.inside_starts += surface.field(CENTER) < 0.0f;
              m.hit_disagreement +=
                  RESULT.has_surface != reference[i].has_surface;
              m.image_error += fabsf(depth_pixel(RESULT, FAR) -
                                     depth_pixel(reference[i], FAR));
            }
            const double ELAPSED = std::chrono::duration<double, std::milli>(
                                       std::chrono::steady_clock::now() - START)
                                       .count();
            m.total_ms += ELAPSED;
            m.peak_ms = std::max(m.peak_ms, ELAPSED);
          }
        }
      }
    }
  }
  for (size_t i = 0; i < BUDGETS.size(); ++i) {
    const Metrics &M = metrics[i];
    std::printf("%s,%d,%llu,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.6f\n", name,
                BUDGETS[i], static_cast<unsigned long long>(M.rays),
                static_cast<double>(M.queries) / M.rays, M.peak_queries,
                static_cast<double>(M.unresolved) / M.rays,
                static_cast<double>(M.reference_unresolved) / M.rays,
                static_cast<double>(M.hit_disagreement) / M.rays,
                M.image_error / M.rays, M.total_ms / 72.0, M.peak_ms,
                static_cast<double>(M.inside_starts) / M.rays);
  }
}

} // namespace

int main() {
  std::puts("surface,budget,rays,mean_queries,peak_queries,unresolved_fraction,"
            "reference_unresolved_fraction,hit_disagreement_fraction,"
            "depth_image_mae,mean_frame_ms,peak_frame_ms,inside_fraction");
  measure<SDF::CosineSurface>("cosine");
  measure<SDF::GyroidSurface>("gyroid");
}
