/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file octet_trace.h
 * @brief Prepared octet lattice tracing and compositing. */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

#include "core/render/ray/camera.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"

/** @brief Front-to-back ray traces of the 3D and 4D octet trusses. */
namespace SDF::OctetTrace {

/** @brief One ray's covered plane crossings; the caller owns it. */
struct CrossingStorage {
  static constexpr int CAPACITY = 64;    ///< Maximum crossings per ray.
  std::array<float, CAPACITY> distances; ///< Ray distance of each crossing.
  std::array<float, CAPACITY> coverages; ///< Coverage of each crossing.
};

/** @brief Premultiplied ray color and how its traversal ended. */
using Sample = Raycast::TraceSample;

/**
 * @brief Composites an octet adapter's plane crossings front to back.
 * @details Matches Raycast::trace_events over the adapter, whose single merge
 * group keeps the first distance and largest coverage of coincident crossings.
 * @tparam Events Octet event adapter (OctetEvents or OctetEvents4).
 * @param events Initialized event streams for the ray.
 * @param interval Ray-parameter range to trace.
 * @param limits Candidate and layer budgets.
 * @param appearance Depth appearance applied to each layer.
 * @return Premultiplied color and trace status.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace_events(const Events &events, Raycast::Interval interval,
             const Raycast::TraceLimits &limits,
             const Raycast::Appearance &appearance) {
  constexpr size_t SLOTS = Events::OWNER_CAPACITY;
  std::array<float, SLOTS> next;
  std::array<float, SLOTS> step;
  std::array<uint8_t, SLOTS> stream;
  size_t count = 0;
  for (uint8_t i = 0; i < Events::STREAM_COUNT; ++i)
    if (events.active(i)) {
      next[count] = events.next[i];
      step[count] = events.step[i];
      stream[count++] = i;
    }
  for (; count < SLOTS; ++count)
    next[count] = INFINITY;
  Sample result;
  LayerComposite composite;
  int candidates = 0;
  int layers = 0;
  bool pending = false;
  float pending_t = 0.0f;
  float pending_coverage = 0.0f;
  float group_end = 0.0f;
  const auto flush = [&]() __attribute__((always_inline)) {
    if (layers >= limits.max_layers) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      return false;
    }
    ++layers;
    HS_PROFILE_DEEP(hl_layer_composite);
    appearance.composite(composite, pending_t, pending_coverage);
    pending = false;
    if (composite.saturated()) {
      result.status = Raycast::TraceStatus::SATURATED;
      return false;
    }
    return true;
  };
  while (true) {
    HS_PROFILE_DEEP(hl_event_step);
    size_t first = 0;
    for (size_t k = 1; k < SLOTS; ++k)
      if (next[k] < next[first])
        first = k;
    const float t = next[first];
    if (!(t <= interval.far))
      break;
    if (pending && t > group_end && !flush()) {
      result.color = composite.premultiplied();
      return result;
    }
    if (candidates >= limits.max_candidates) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++candidates;
    if (t >= interval.near) {
      uint32_t feature;
      const float coverage = events.coverage(stream[first], t, feature);
      if (coverage <= 0.0f) {
        HS_PROFILE_DEEP(hl_event_miss);
      } else if (!pending) {
        pending = true;
        pending_t = t;
        pending_coverage = coverage;
        group_end =
            Raycast::merge_window_end(t, Raycast::MERGE_RELATIVE_TOLERANCE);
      } else if (coverage > pending_coverage) {
        pending_coverage = coverage;
      }
    }
    next[first] = t + step[first];
  }
  if (pending)
    flush();
  result.color = composite.premultiplied();
  return result;
}

/** @brief Covered crossing view over the octet trace storage. */
using CoveredCrossings =
    Raycast::CoveredCrossings<CrossingStorage::CAPACITY, float *>;

namespace detail::reference {

/**
 * @brief trace_events() by per-stream walks over CoveredCrossings.
 * @details Each stream walks its own crossings with the same accumulated
 * distances trace_events() pops. A ray with more crossings than the candidate
 * budget defers to trace_events(), which truncates them in distance order.
 * @tparam Events Octet event adapter (OctetEvents or OctetEvents4).
 * @param events Initialized event streams for the ray.
 * @param interval Ray-parameter range to trace.
 * @param limits Candidate and layer budgets.
 * @param appearance Depth appearance applied to each layer.
 * @param storage Per-ray crossing scratch.
 * @return Premultiplied color and trace status.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace_sorted(const Events &events, Raycast::Interval interval,
             const Raycast::TraceLimits &limits,
             const Raycast::Appearance &appearance, CrossingStorage &storage) {
  CoveredCrossings covered(storage.distances, storage.coverages);
  int crossings = 0;
  const int budget =
      std::min(limits.max_candidates, CoveredCrossings::CAPACITY);
  for (uint8_t stream = 0; stream < Events::STREAM_COUNT; ++stream) {
    if (!events.active(stream))
      continue;
    const float step = events.step[stream];
    for (float t = events.next[stream]; t <= interval.far; t += step) {
      if (++crossings > budget)
        return trace_events(events, interval, limits, appearance);
      if (t < interval.near)
        continue;
      uint32_t feature;
      const float coverage = events.coverage(stream, t, feature);
      if (coverage > 0.0f)
        covered.insert(t, coverage);
    }
  }
  return covered.composite<true, Sample>(appearance, limits);
}

} // namespace detail::reference

/**
 * @brief The 3D octet trace with every owner's strut pairs in registers.
 * @details Ownership follows SDF::OctetEvents (a plane-family pair goes to the
 * family the ray crosses faster, the lower on ties), so the fastest owner holds
 * three pairs, the next two and the third one. Owners walk with OctetEvents'
 * arithmetic, so the covered crossings are OctetEvents'.
 * @param direction Unit view direction.
 * @param camera Camera the projection was prepared for.
 * @param projection The octet planes projected onto view directions.
 * @param footprint Pixel footprint widening the struts with distance.
 * @param limits Candidate and layer budgets.
 * @param appearance Fog, palette and gain of each composited layer.
 * @param storage Scratch the covered crossings are sorted in.
 * @return The premultiplied composite and how the traversal ended.
 */
__attribute__((always_inline)) inline Sample
trace_3d(const math::Vector &direction, const Raycast::PreparedCamera &camera,
         const OctetEvents::PreparedProjection &projection,
         const Raycast::Footprint &footprint,
         const Raycast::TraceLimits &limits,
         const Raycast::Appearance &appearance, CrossingStorage &storage) {
  const float near = camera.interval.near;
  const float far = camera.interval.far;
  const float wire_radius = projection.wire_radius;
  std::array<float, 4> speeds;
  std::array<float, 4> positions;
  for (int family = 0; family < 4; ++family) {
    speeds[family] = math::dot(direction, projection.normals[family]);
    positions[family] =
        projection.offsets[family] + camera.radial_start * speeds[family];
  }
  std::array<uint8_t, 4> order;
  for (uint8_t family = 0; family < 4; ++family) {
    uint8_t rank = 0;
    for (uint8_t other = 0; other < 4; ++other)
      rank += other < family ? fabsf(speeds[other]) >= fabsf(speeds[family])
                             : fabsf(speeds[other]) > fabsf(speeds[family]);
    order[rank] = family;
  }
  CoveredCrossings covered(storage.distances, storage.coverages);
  int crossings = 0;
  const int budget =
      std::min(limits.max_candidates, CoveredCrossings::CAPACITY);
  const auto walk =
      [&]<size_t PAIRS>(uint8_t owner, const std::array<uint8_t, PAIRS> &others)
          __attribute__((always_inline)) {
            const float speed = speeds[owner];
            if (speed == 0.0f)
              return true;
            const float numerator_value = speed * speed * projection.spacing2;
            std::array<float, PAIRS> denominators, bases, rates;
            for (size_t pair = 0; pair < PAIRS; ++pair) {
              const uint8_t other = others[pair];
              const float a = speeds[std::min(owner, other)];
              const float b = speeds[std::max(owner, other)];
              denominators[pair] = a * a + b * b + (2.0f / 3.0f) * a * b;
              bases[pair] = positions[other];
              rates[pair] = speeds[other];
            }
            const float position = positions[owner] + near * speed;
            const float plane =
                speed > 0.0f ? ceilf(position) : floorf(position);
            const float inverse = 1.0f / speed;
            const float step = fabsf(inverse);
            const float first = near + (plane - position) * inverse;
            if (!Raycast::finite(first) || !Raycast::finite(step) ||
                !(step > 0.0f))
              return true;
            for (float t = first; t <= far; t += step) {
              if (++crossings > budget)
                return false;
              if (t < near)
                continue;
              float numerator = INFINITY;
              float denominator = 1.0f;
              for (size_t pair = 0; pair < PAIRS; ++pair) {
                const float u = bases[pair] + t * rates[pair];
                const float residual = u - roundf(u);
                const float n = residual * residual * numerator_value;
                if (n * denominator < numerator * denominators[pair]) {
                  numerator = n;
                  denominator = denominators[pair];
                }
              }
              const float width = footprint.at(t);
              const float support = wire_radius + .5f * width;
              if (numerator > support * support * denominator)
                continue;
              const float field = sqrtf(numerator / denominator) - wire_radius;
              const float coverage = footprint_coverage(width, field);
              if (coverage > 0.0f)
                covered.insert(t, coverage);
            }
            return true;
          };
  constexpr uint8_t OTHERS[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  const uint8_t fastest = order[0];
  const std::array<uint8_t, 3> fastest_pairs{
      OTHERS[fastest][0], OTHERS[fastest][1], OTHERS[fastest][2]};
  const std::array<uint8_t, 2> second_pairs{std::min(order[2], order[3]),
                                            std::max(order[2], order[3])};
  const std::array<uint8_t, 1> third_pair{order[3]};
  if (!walk.template operator()<3>(fastest, fastest_pairs) ||
      !walk.template operator()<2>(order[1], second_pairs) ||
      !walk.template operator()<1>(order[2], third_pair)) {
    const SDF::OctetEvents events(projection, direction, camera.radial_start,
                                  near, footprint);
    return trace_events(events, camera.interval, limits, appearance);
  }
  return covered.composite<true, Sample>(appearance, limits);
}

/** @brief A coordinate pair and its complementary free coordinates. */
struct OctetPair4 {
  uint8_t i; ///< First constrained coordinate.
  uint8_t j; ///< Second constrained coordinate.
  uint8_t k; ///< First free coordinate.
  uint8_t l; ///< Second free coordinate.
};

/** @brief Ranks ambient coordinate magnitudes with stable index ties.
 * @param ambient Ambient ray direction.
 * @return Coordinate indices in ascending magnitude. */
__attribute__((always_inline)) inline std::array<uint8_t, 4>
canonical_octet_order(const std::array<float, 4> &ambient) {
  // Magnitudes rank by their IEEE bits, which order non-negative floats.
  std::array<uint32_t, 4> magnitudes;
  for (int axis = 0; axis < 4; ++axis) {
    std::memcpy(&magnitudes[axis], &ambient[axis], sizeof(float));
    magnitudes[axis] &= 0x7FFFFFFFu;
  }
  std::array<uint8_t, 4> order;
  for (uint8_t axis = 0; axis < 4; ++axis) {
    uint8_t rank = 0;
    for (uint8_t other = 0; other < 4; ++other)
      rank += other < axis ? magnitudes[other] <= magnitudes[axis]
                           : magnitudes[other] < magnitudes[axis];
    order[rank] = axis;
  }
  return order;
}

/** @brief Prepares canonical plane positions, speeds and inverse speeds.
 * @param origins Canonical ray origin.
 * @param rates Canonical coordinate rates.
 * @param near_value Near ray bound.
 * @param positions Receives family positions.
 * @param speeds Receives family speeds.
 * @param inverses Receives reciprocal speeds. */
__attribute__((always_inline)) inline void prepare_octet_planes(
    const std::array<float, 4> &origins, const std::array<float, 4> &rates,
    float near_value, std::array<float, 4> &positions,
    std::array<float, 4> &speeds, std::array<float, 4> &inverses) {
  float position = 0.0f, speed = 0.0f;
  std::array<float, 4> near;
  for (int m = 0; m < 4; ++m) {
    near[m] = origins[m] + near_value * rates[m];
    position += near[m];
    speed += rates[m];
  }
  positions[0] = 0.5f * position;
  speeds[0] = 0.5f * speed;
  for (int m = 0; m < 3; ++m) {
    positions[m + 1] = positions[0] - near[m];
    speeds[m + 1] = speeds[0] - rates[m];
  }
  const float product = speeds[0] * speeds[1] * speeds[2] * speeds[3];
  if (product > 1e-30f) {
    const float inverse = 1.0f / product;
    const float first_pair = speeds[0] * speeds[1];
    const float second_pair = speeds[2] * speeds[3];
    inverses = {
        inverse * speeds[1] * second_pair, inverse * speeds[0] * second_pair,
        inverse * first_pair * speeds[3], inverse * first_pair * speeds[2]};
  } else {
    for (int f = 0; f < 4; ++f)
      inverses[f] = 1.0f / speeds[f];
  }
}

/** @brief Constants and retained crossings for the canonical octet walk. */
struct OctetWalk4 {
  const std::array<float, 4> &components; ///< Canonical unit direction.
  const std::array<float, 4> &rates;      ///< Scaled coordinate rates.
  const std::array<float, 4> &origins;    ///< Scaled canonical origin.
  const std::array<float, 4> &positions; ///< Plane positions at the near bound.
  const std::array<float, 4> &speeds;    ///< Plane crossing rates.
  const std::array<float, 4> &inverses;  ///< Reciprocal plane rates.
  float near_value;                      ///< Near ray bound.
  float far;                             ///< Far ray bound.
  float inverse_scale;                   ///< Reciprocal framework scale.
  float scale;                           ///< Framework world scale.
  float wire_radius;                     ///< Strut radius in world units.
  float support_base;                    ///< Support at the ray origin.
  float support_rate;                    ///< Support growth along the ray.
  const Raycast::Footprint &footprint;   ///< Pixel footprint.
  int &crossings;                        ///< Counted candidate crossings.
  int budget;                            ///< Maximum candidate crossings.
  CoveredCrossings &covered;             ///< Retained crossings.
};

/** @brief Evaluates the exact owned-class distance of a surviving crossing.
 * @param state Canonical walk constants.
 * @param sign Constrained coordinate sign.
 * @param pairs Owned line classes.
 * @param bounds Integer distance bounds.
 * @param threshold Support threshold in fixed-point units.
 * @param support2 Squared filtered support.
 * @param t Ray distance. */
template <size_t CLASSES>
__attribute__((always_inline)) inline void
retain_octet_crossing(const OctetWalk4 &state, float sign,
                      const std::array<OctetPair4, CLASSES> &pairs,
                      const std::array<uint32_t, CLASSES> &bounds,
                      uint32_t threshold, float support2, float t) {
  std::array<float, 4> residual;
  std::array<float, 4> rounded;
  for (int m = 0; m < 4; ++m) {
    const float q = state.origins[m] + state.rates[m] * t;
    rounded[m] = roundf(q);
    residual[m] = q - rounded[m];
  }
  int parity = 0;
  for (int m = 0; m < 4; ++m)
    parity += static_cast<int>(rounded[m]);
  const bool odd = parity & 1;
  float numerator = INFINITY;
  float denominator = 1.0f;
  for (size_t c = 0; c < CLASSES; ++c) {
    const auto &pair = pairs[c];
    float rk = residual[pair.k];
    float rl = residual[pair.l];
    // A class the plane bound clears cannot reach the support, and
    // neither can a farther class it would otherwise have beaten.
    if (bounds[c] > threshold)
      continue;
    const float along =
        state.components[pair.i] + sign * state.components[pair.j];
    const float denominator_value = 1.0f - 0.5f * along * along;
    float across = residual[pair.i] - sign * residual[pair.j];
    // |across| <= 1, and ties round to even: -1, 0 or 1 as the
    // strict half-way comparisons give.
    const float shift = rintf(across);
    across -= shift;
    if ((shift != 0.0f) != odd) {
      const float cost = 0.5f - fabsf(across);
      const float cost_k = 1.0f - 2.0f * fabsf(rk);
      const float cost_l = 1.0f - 2.0f * fabsf(rl);
      if (cost <= cost_k && cost <= cost_l)
        across -= copysignf(1.0f, across);
      else if (cost_k <= cost_l)
        rk -= copysignf(1.0f, rk);
      else
        rl -= copysignf(1.0f, rl);
    }
    const float offset2 = 0.5f * across * across + rk * rk + rl * rl;
    const float transverse =
        0.5f * (state.components[pair.i] - sign * state.components[pair.j]);
    const float dot = across * transverse + rk * state.components[pair.k] +
                      rl * state.components[pair.l];
    const float n = fmaxf(0.0f, offset2 * denominator_value - dot * dot);
    if (n * denominator < numerator * denominator_value) {
      numerator = n;
      denominator = denominator_value;
    }
  }
  if (numerator > support2 * denominator)
    return;
  const float field =
      state.scale * sqrtf(numerator / denominator) - state.wire_radius;
  const float width = state.footprint.at(t);
  const float coverage = footprint_coverage(width, field);
  // The counted walk can round its last crossing past the far bound.
  if (coverage > 0.0f && t <= state.far)
    state.covered.insert(t, coverage);
}

/** @brief Walks one canonical octet plane family through its candidate budget.
 * @param state Canonical walk constants and retained crossings.
 * @param flipped Flipped coordinate, or minus one for all positive.
 * @param sign Constrained coordinate sign.
 * @param pairs Owned line classes.
 * @return Whether the candidate budget was sufficient. */
template <size_t CLASSES>
__attribute__((always_inline)) inline bool
walk_octet_family(const OctetWalk4 &state, int flipped, float sign,
                  const std::array<OctetPair4, CLASSES> &pairs) {
  // Only the largest-coordinate sum can be parallel in the canonical unit frame.
  if constexpr (CLASSES == 1) {
    const auto &pair = pairs[0];
    const float along =
        state.components[pair.i] + sign * state.components[pair.j];
    const float denominator = 1.0f - 0.5f * along * along;
    uint32_t bits;
    std::memcpy(&bits, &denominator, sizeof(bits));
    const bool owned = bits - 1 < 0x7F7FFFFFu;
    if (!owned)
      return true;
  }
  const float position = state.positions[flipped + 1];
  // Canonical owner speeds are non-negative.
  const float speed = state.speeds[flipped + 1];
  if (speed == 0.0f)
    return true;
  const float plane = ceilf(position);
  const float inverse = state.inverses[flipped + 1];
  const float step = inverse;
  const float first = state.near_value + (plane - position) * inverse;
  if (!Raycast::finite(first) || !Raycast::finite(step) || !(step > 0.0f))
    return true;
  // The first plane lies at or past the near bound. A crossing within rounding of the far bound
  // carries no fog-weighted opacity either way.
  const int count = first <= state.far
                        ? static_cast<int>((state.far - first) * speed) + 1
                        : 0;
  state.crossings += count;
  if (state.crossings > state.budget)
    return false;
  if (count == 0)
    return true;
  // The plane bound runs in 16-bit fixed point: the low 16 bits of Q * 2^16
  // wrap exactly as Q mod 1, so their signed value is Q's nearest-integer
  // residual. FIXED_ERROR bounds quantization and t drift over the crossing budget.
  constexpr float FIXED_ONE = 65536.0f;
  constexpr float FIXED_ERROR = 2.0e-3f;
  const float inverse_plane = inverse * state.inverse_scale;
  const float threshold_scale =
      (FIXED_ONE * FIXED_ONE / 0.9999f) * inverse_plane * inverse_plane;
  const float THRESHOLD_MARGIN = 3 * FIXED_ERROR * FIXED_ONE * FIXED_ONE;
  const float advance_step = count > 1 ? step : 0.0f;
  std::array<uint32_t, 4> starts, advances;
  for (int m = 0; m < 4; ++m) {
    starts[m] = static_cast<uint32_t>(static_cast<int32_t>(
        (state.origins[m] + state.rates[m] * first) * FIXED_ONE));
    advances[m] = static_cast<uint32_t>(
        static_cast<int32_t>(state.rates[m] * advance_step * FIXED_ONE));
  }
  float t = first - step;
  for (uint32_t index = 0; index < static_cast<uint32_t>(count); ++index) {
    t += step;
    const float support = state.support_base + state.support_rate * t;
    const float support2 = support * support;
    std::array<uint32_t, 4> fixed;
    std::array<uint32_t, 4> squares;
    for (int m = 0; m < 4; ++m) {
      fixed[m] = starts[m] + index * advances[m];
      const int32_t residual_value =
          static_cast<int16_t>(static_cast<uint16_t>(fixed[m]));
      squares[m] = static_cast<uint32_t>(residual_value * residual_value);
    }
    // Each class's line lies no nearer than half its squared across
    // residual plus its free coordinates' squared residuals.
    const uint32_t threshold = static_cast<uint32_t>(
        fminf(support2 * threshold_scale + THRESHOLD_MARGIN, 4294967040.0f));
    uint32_t free_nearest = UINT32_MAX;
    if constexpr (CLASSES == 6) {
      free_nearest = std::min(
          std::min(squares[0] + squares[1], squares[2] + squares[3]),
          std::min(squares[0], squares[1]) + std::min(squares[2], squares[3]));
    } else if constexpr (CLASSES == 2) {
      free_nearest = std::min(squares[pairs[0].k], squares[pairs[1].k]) +
                     std::min(squares[pairs[0].l], squares[pairs[1].l]);
    } else {
      for (const auto &pair : pairs)
        free_nearest =
            std::min(free_nearest, squares[pair.k] + squares[pair.l]);
    }
    if (free_nearest > threshold) [[likely]]
      continue;
    std::array<uint32_t, CLASSES> bounds;
    uint32_t nearest = UINT32_MAX;
    for (size_t c = 0; c < CLASSES; ++c) {
      const auto &pair = pairs[c];
      const int32_t across = static_cast<int16_t>(
          static_cast<uint16_t>(sign > 0.0f ? fixed[pair.i] - fixed[pair.j]
                                            : fixed[pair.i] + fixed[pair.j]));
      bounds[c] = (static_cast<uint32_t>(across * across) >> 1) +
                  squares[pair.k] + squares[pair.l];
      nearest = std::min(nearest, bounds[c]);
    }
    if (nearest > threshold)
      continue;
    retain_octet_crossing(state, sign, pairs, bounds, threshold, support2, t);
  }
  return true;
}

/**
 * @brief The 4D octet trace in the ray's canonical frame.
 * @details Reflecting the ray to non-negative components and ordering them by
 * magnitude (exact D4 symmetries) fixes SDF::OctetEvents4's ownership, so each
 * crossing evaluates a fixed class list with OctetEvents4's arithmetic. A
 * crossing whose every owned class fails the |n . v| line-distance bound
 * against the support is skipped before the class search.
 * @param direction Unit view direction.
 * @param camera Camera the projection was prepared for.
 * @param framework The D4 framework's cell size and strut radius.
 * @param projection The framework's slice projected onto view directions.
 * @param footprint Pixel footprint widening the struts with distance.
 * @param limits Candidate and layer budgets.
 * @param appearance Fog, palette and gain of each composited layer.
 * @param storage Scratch the covered crossings are sorted in.
 * @return The premultiplied composite and how the traversal ended.
 */
__attribute__((always_inline)) inline Sample
trace_4d(const math::Vector &direction, const Raycast::PreparedCamera &camera,
         const OctetFramework4 &framework,
         const OctetEvents4::PreparedProjection &projection,
         const Raycast::Footprint &footprint,
         const Raycast::TraceLimits &limits,
         const Raycast::Appearance &appearance, CrossingStorage &storage) {
  const float near_value = camera.interval.near;
  const float far = camera.interval.far;
  const float inverse_scale = projection.inverse_scale;
  const float scale = SDF::OctetFramework4::HALF_CUBE * framework.cell_size;
  const float wire_radius = framework.wire_radius;
  std::array<float, 4> ambient;
  for (int axis = 0; axis < 4; ++axis)
    ambient[axis] = math::dot(direction, projection.embedding[axis]);
  const auto order = canonical_octet_order(ambient);
  std::array<float, 4> components, rates, origins;
  for (int m = 0; m < 4; ++m) {
    const uint8_t axis = order[m];
    const float local = ambient[axis] * inverse_scale;
    const float start = projection.origin[axis] + camera.radial_start * local;
    components[m] = fabsf(ambient[axis]);
    rates[m] = fabsf(local);
    origins[m] = ambient[axis] < 0.0f ? -start : start;
  }
  std::array<float, 4> positions, speeds, inverses;
  prepare_octet_planes(origins, rates, near_value, positions, speeds, inverses);
  const float support_rate = .5f * footprint.angular_radius * inverse_scale;
  const float support_base =
      wire_radius * inverse_scale + support_rate * footprint.radial_start;
  CoveredCrossings covered(storage.distances, storage.coverages);
  int crossings = 0;
  const int budget =
      std::min(limits.max_candidates, CoveredCrossings::CAPACITY);
  const OctetWalk4 walk{components,    rates,     origins,     positions,
                        speeds,        inverses,  near_value,  far,
                        inverse_scale, scale,     wire_radius, support_base,
                        support_rate,  footprint, crossings,   budget,
                        covered};
  constexpr std::array<OctetPair4, 6> DIFFERENCES{{{0, 1, 2, 3},
                                                   {0, 2, 1, 3},
                                                   {0, 3, 1, 2},
                                                   {1, 2, 0, 3},
                                                   {1, 3, 0, 2},
                                                   {2, 3, 0, 1}}};
  constexpr std::array<OctetPair4, 3> FIRST_SUMS{
      {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}}};
  constexpr std::array<OctetPair4, 2> SECOND_SUMS{{{1, 2, 0, 3}, {1, 3, 0, 2}}};
  constexpr std::array<OctetPair4, 1> THIRD_SUM{{{2, 3, 0, 1}}};
  if (!walk_octet_family(walk, -1, -1.0f, DIFFERENCES) ||
      !walk_octet_family(walk, 0, 1.0f, FIRST_SUMS) ||
      !walk_octet_family(walk, 1, 1.0f, SECOND_SUMS) ||
      !walk_octet_family(walk, 2, 1.0f, THIRD_SUM)) {
    const SDF::OctetEvents4 events(framework, projection, direction,
                                   camera.radial_start, near_value, footprint);
    return trace_events(events, camera.interval, limits, appearance);
  }
  return covered.composite<true, Sample>(appearance, limits);
}

} // namespace SDF::OctetTrace
