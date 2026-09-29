/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cstring>

#include "core/render/ray/camera.h"
#include "core/render/ray/shade.h"
#include "core/render/sdf/framework.h"
#include "core/render/sdf/cellular_wire.h"
#include "core/render/sdf/affine_lattice.h"
#include "core/render/sdf/periodic_shells.h"

namespace HyperLatticeExperimental {

enum class Geometry : uint8_t {
  OCTET,
  DIAMOND,
  HEXAGONAL,
  RHOMBIC,
  AFFINE_CUBIC,
  SHELLS
};

/** @brief Frame settings; camera distances and near fading use world units. */
struct Settings {
  Raycast::SamplingDomain domain = Raycast::SamplingDomain::SPATIAL_3D;
  float cell_size = 1.0f;
  float wire_radius = 0.055f;
  float radial_start = 0.0f;
  float far_distance = 7.0f;
  float near_fade = 0.5f;
  float aa_strength = 1.0f;
  math::Vec4 center{};
  math::Mat4 embedding = math::Mat4::identity();
  float pixel_half_angle = 0.0f;
  const BakedPalette *palette = nullptr;
  Geometry geometry = Geometry::OCTET;
  float shear = .55f;
  float stretch = 1.4f;
  float shell_radius = .30f;
};

struct Prepared {
  Raycast::PreparedCamera camera;
  Raycast::Footprint footprint;
  Raycast::Appearance appearance;
  Raycast::TraceLimits limits;
  SDF::OctetFramework octet;
  SDF::OctetEvents::PreparedProjection octet_projection{};
  SDF::OctetFramework4 octet4;
  SDF::OctetEvents4::PreparedProjection octet4_projection{};
  const SDF::CellularWire::Geometry *cellular = nullptr;
  SDF::PreparedPeriodicShells periodic_shells;
  bool valid = false;
  Geometry geometry = Geometry::OCTET;
  float shear = .55f;
  float stretch = 1.4f;
  float shell_radius = .30f;
};

/** @brief Premultiplied ray color and how its traversal ended. */
struct Sample {
  Pixel color;
  Raycast::TraceStatus status = Raycast::TraceStatus::RANGE_COMPLETE;
};

HS_FLASH_INLINE inline Prepared prepare(const Settings &settings) {
  Prepared result;
  if (!Raycast::finite(settings.cell_size) || settings.cell_size <= 0.0f ||
      !Raycast::finite(settings.near_fade) || settings.near_fade <= 0.0f ||
      !Raycast::finite(settings.aa_strength) || settings.aa_strength < 0.0f ||
      !Raycast::finite(settings.pixel_half_angle) ||
      settings.pixel_half_angle < 0.0f || !settings.palette)
    return result;

  result.camera.domain = settings.domain;
  result.camera.center = settings.center;
  result.camera.embedding = settings.embedding;
  result.camera.radial_start = settings.radial_start;
  result.camera.interval = {0.0f, settings.far_distance};
  result.footprint = {settings.pixel_half_angle * settings.aa_strength,
                      settings.radial_start};
  result.appearance = {1.0f / settings.far_distance, 0.0f,
                       1.0f / settings.near_fade, settings.palette};
  result.limits.max_candidates = 64;
  result.limits.max_layers = 32;
  result.octet.cell_size = settings.cell_size;
  result.octet.wire_radius = settings.wire_radius;
  result.octet4.cell_size = settings.cell_size;
  result.octet4.wire_radius = settings.wire_radius;
  result.geometry = settings.geometry;
  result.shear = settings.shear;
  result.stretch = settings.stretch;
  result.shell_radius = settings.shell_radius;
  result.valid = result.camera.valid() &&
                 Raycast::finite(result.footprint.angular_radius) &&
                 Raycast::finite(result.appearance.inv_far) &&
                 Raycast::finite(result.appearance.near_inv_span) &&
                 (settings.domain == Raycast::SamplingDomain::SLICE_4D
                      ? result.octet4.valid()
                      : result.octet.valid());
  if (!result.valid)
    return result;
  if (settings.geometry != Geometry::OCTET) {
    const bool CELLULAR = settings.geometry == Geometry::DIAMOND ||
                          settings.geometry == Geometry::HEXAGONAL ||
                          settings.geometry == Geometry::RHOMBIC;
    result.valid =
        (CELLULAR && settings.domain == Raycast::SamplingDomain::SPATIAL_3D) ||
        settings.geometry == Geometry::AFFINE_CUBIC ||
        settings.geometry == Geometry::SHELLS;
    if (CELLULAR)
      result.cellular = &SDF::CellularWire::geometry(
          settings.geometry == Geometry::DIAMOND
              ? SDF::CellularWire::Kind::DIAMOND
          : settings.geometry == Geometry::HEXAGONAL
              ? SDF::CellularWire::Kind::HEXAGONAL
              : SDF::CellularWire::Kind::RHOMBIC);
    if (settings.geometry == Geometry::SHELLS) {
      result.periodic_shells =
          SDF::prepare_periodic_shells(result.camera, settings.cell_size,
                                       settings.shell_radius, result.footprint);
      result.valid = result.periodic_shells.valid;
    }
    return result;
  }
  const auto &E = settings.embedding.m;
  if (settings.domain == Raycast::SamplingDomain::SPATIAL_3D) {
    const auto FAMILIES = result.octet.plane_families();
    auto &projection = result.octet_projection;
    const float SPACING = FAMILIES[0].spacing;
    projection.spacing2 = SPACING * SPACING;
    projection.wire_radius = result.octet.wire_radius;
    const math::Vector ORIGIN(settings.center[0], settings.center[1],
                              settings.center[2]);
    for (size_t i = 0; i < FAMILIES.size(); ++i) {
      const math::Vector NORMAL = FAMILIES[i].normal / SPACING;
      projection.offsets[i] = math::dot(ORIGIN - result.octet.origin, NORMAL);
      projection.normals[i] = {
          NORMAL.x * E[0][0] + NORMAL.y * E[1][0] + NORMAL.z * E[2][0],
          NORMAL.x * E[0][1] + NORMAL.y * E[1][1] + NORMAL.z * E[2][1],
          NORMAL.x * E[0][2] + NORMAL.y * E[1][2] + NORMAL.z * E[2][2]};
    }
  } else {
    auto &projection = result.octet4_projection;
    projection.inverse_scale =
        1.0f / (SDF::OctetFramework4::HALF_CUBE * settings.cell_size);
    for (int i = 0; i < 4; ++i) {
      projection.embedding[i] = math::Vector(E[i][0], E[i][1], E[i][2]);
      projection.origin[i] = (settings.center[i] - result.octet4.origin[i]) *
                             projection.inverse_scale;
    }
  }
  return result;
}

/**
 * @brief Composites an octet adapter's plane crossings front to back.
 * @details Matches Raycast::trace_events over the adapter, whose single merge
 * group keeps the first distance and largest coverage of coincident crossings.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace(const Events &events, Raycast::Interval interval,
      const Raycast::TraceLimits &limits,
      const Raycast::Appearance &appearance) {
  constexpr size_t SLOTS = Events::OWNER_CAPACITY;
  constexpr float RELATIVE_TOLERANCE = 1.0e-4f;
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
    const float T = next[first];
    if (!(T <= interval.far))
      break;
    if (pending && T > group_end && !flush()) {
      result.color = composite.premultiplied();
      return result;
    }
    if (candidates >= limits.max_candidates) {
      result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++candidates;
    if (T >= interval.near) {
      uint32_t feature;
      const float COVERAGE = events.coverage(stream[first], T, feature);
      if (COVERAGE <= 0.0f) {
        HS_PROFILE_DEEP(hl_event_miss);
      } else if (!pending) {
        pending = true;
        pending_t = T;
        pending_coverage = COVERAGE;
        group_end = T + RELATIVE_TOLERANCE * fmaxf(1.0f, T);
      } else if (COVERAGE > pending_coverage) {
        pending_coverage = COVERAGE;
      }
    }
    next[first] = T + step[first];
  }
  if (pending)
    flush();
  result.color = composite.premultiplied();
  return result;
}

/**
 * @brief Covered crossings in distance order, composited as trace() groups.
 * @details An uncovered crossing never closes a merge group, so the groups are
 * runs of covered crossings within the relative tolerance of their first
 * distance, each one layer at that distance with the run's largest coverage.
 */
struct CoveredCrossings {
  static constexpr int CAPACITY = 64;
  std::array<float, CAPACITY> distances;
  std::array<float, CAPACITY> coverages;
  int count = 0;

  __attribute__((always_inline)) void insert(float t, float coverage) {
    int slot = count++;
    for (; slot > 0 && distances[slot - 1] > t; --slot) {
      distances[slot] = distances[slot - 1];
      coverages[slot] = coverages[slot - 1];
    }
    distances[slot] = t;
    coverages[slot] = coverage;
  }

  __attribute__((always_inline)) Sample
  composite(const Raycast::TraceLimits &limits,
            const Raycast::Appearance &appearance) const {
    constexpr float RELATIVE_TOLERANCE = 1.0e-4f;
    Sample result;
    LayerComposite layers;
    int composited = 0;
    for (int index = 0; index < count;) {
      if (composited >= limits.max_layers) {
        result.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
        break;
      }
      ++composited;
      const float T = distances[index];
      const float GROUP_END = T + RELATIVE_TOLERANCE * fmaxf(1.0f, T);
      float coverage = coverages[index];
      for (++index; index < count && distances[index] <= GROUP_END; ++index)
        coverage = fmaxf(coverage, coverages[index]);
      appearance.composite(layers, T, coverage);
      if (layers.saturated()) {
        result.status = Raycast::TraceStatus::SATURATED;
        break;
      }
    }
    result.color = layers.premultiplied();
    return result;
  }
};

/**
 * @brief trace() by per-stream walks over CoveredCrossings.
 * @details Each stream walks its own crossings with the same accumulated
 * distances trace() pops. A ray with more crossings than the candidate budget
 * defers to trace(), which truncates them in distance order.
 */
template <typename Events>
__attribute__((always_inline)) inline Sample
trace_sorted(const Events &events, Raycast::Interval interval,
             const Raycast::TraceLimits &limits,
             const Raycast::Appearance &appearance) {
  CoveredCrossings covered;
  int crossings = 0;
  const int BUDGET =
      std::min(limits.max_candidates, CoveredCrossings::CAPACITY);
  for (uint8_t stream = 0; stream < Events::STREAM_COUNT; ++stream) {
    if (!events.active(stream))
      continue;
    const float STEP = events.step[stream];
    for (float t = events.next[stream]; t <= interval.far; t += STEP) {
      if (++crossings > BUDGET)
        return trace(events, interval, limits, appearance);
      if (t < interval.near)
        continue;
      uint32_t feature;
      const float COVERAGE = events.coverage(stream, t, feature);
      if (COVERAGE > 0.0f)
        covered.insert(t, COVERAGE);
    }
  }
  return covered.composite(limits, appearance);
}

/**
 * @brief The 3D octet trace with every owner's strut pairs in registers.
 * @details SDF::OctetEvents gives a pair of plane families to the family the
 * ray crosses faster, the lower family on ties. Ranking the families that way
 * fixes each owner's share: the fastest owns three pairs, the next two and the
 * third one. Each owner walks its crossings over its pairs in family order,
 * with OctetEvents' arithmetic, so the covered crossings are OctetEvents'.
 */
__attribute__((always_inline)) inline Sample
trace_octet_3d(const math::Vector &direction, const Prepared &prepared) {
  const auto &camera = prepared.camera;
  const auto &projection = prepared.octet_projection;
  const auto &footprint = prepared.footprint;
  const float NEAR = camera.interval.near;
  const float FAR = camera.interval.far;
  const float WIRE_RADIUS = projection.wire_radius;
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
  CoveredCrossings covered;
  int crossings = 0;
  const int BUDGET =
      std::min(prepared.limits.max_candidates, CoveredCrossings::CAPACITY);
  const auto walk =
      [&]<size_t PAIRS>(uint8_t owner, const std::array<uint8_t, PAIRS> &others)
          __attribute__((always_inline)) {
            const float SPEED = speeds[owner];
            if (SPEED == 0.0f)
              return true;
            std::array<float, PAIRS> numerators, denominators, bases, rates;
            for (size_t pair = 0; pair < PAIRS; ++pair) {
              const uint8_t OTHER = others[pair];
              const float A = speeds[std::min(owner, OTHER)];
              const float B = speeds[std::max(owner, OTHER)];
              numerators[pair] = SPEED * SPEED * projection.spacing2;
              denominators[pair] = A * A + B * B + (2.0f / 3.0f) * A * B;
              bases[pair] = positions[OTHER];
              rates[pair] = speeds[OTHER];
            }
            const float POSITION = positions[owner] + NEAR * SPEED;
            const float PLANE =
                SPEED > 0.0f ? ceilf(POSITION) : floorf(POSITION);
            const float INVERSE = 1.0f / SPEED;
            const float STEP = fabsf(INVERSE);
            const float FIRST = NEAR + (PLANE - POSITION) * INVERSE;
            if (!Raycast::finite(FIRST) || !Raycast::finite(STEP) ||
                !(STEP > 0.0f))
              return true;
            for (float t = FIRST; t <= FAR; t += STEP) {
              if (++crossings > BUDGET)
                return false;
              if (t < NEAR)
                continue;
              float numerator = INFINITY;
              float denominator = 1.0f;
              for (size_t pair = 0; pair < PAIRS; ++pair) {
                const float U = bases[pair] + t * rates[pair];
                const float RESIDUAL = U - roundf(U);
                const float N = RESIDUAL * RESIDUAL * numerators[pair];
                if (N * denominator < numerator * denominators[pair]) {
                  numerator = N;
                  denominator = denominators[pair];
                }
              }
              const float SUPPORT = WIRE_RADIUS + .5f * footprint.at(t);
              if (numerator > SUPPORT * SUPPORT * denominator)
                continue;
              const float FIELD = sqrtf(numerator / denominator) - WIRE_RADIUS;
              const float WIDTH = footprint.at(t);
              const float COVERAGE =
                  WIDTH > 0.0f ? hs::clamp(0.5f - FIELD / WIDTH, 0.0f, 1.0f)
                               : (FIELD <= 0.0f ? 1.0f : 0.0f);
              if (COVERAGE > 0.0f)
                covered.insert(t, COVERAGE);
            }
            return true;
          };
  constexpr uint8_t OTHERS[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  const uint8_t FASTEST = order[0];
  const std::array<uint8_t, 3> FASTEST_PAIRS{
      OTHERS[FASTEST][0], OTHERS[FASTEST][1], OTHERS[FASTEST][2]};
  const std::array<uint8_t, 2> SECOND_PAIRS{std::min(order[2], order[3]),
                                            std::max(order[2], order[3])};
  const std::array<uint8_t, 1> THIRD_PAIR{order[3]};
  if (!walk.template operator()<3>(FASTEST, FASTEST_PAIRS) ||
      !walk.template operator()<2>(order[1], SECOND_PAIRS) ||
      !walk.template operator()<1>(order[2], THIRD_PAIR)) {
    const SDF::OctetEvents events(projection, direction, camera.radial_start,
                                  NEAR, footprint);
    return trace(events, camera.interval, prepared.limits, prepared.appearance);
  }
  return covered.composite(prepared.limits, prepared.appearance);
}

/**
 * @brief The 4D octet trace in the ray's canonical frame.
 * @details The D4 lattice is invariant under coordinate permutations and sign
 * changes, both exact in floating point. Reflecting the ray so each direction
 * component is non-negative and ordering the components by magnitude fixes
 * SDF::OctetEvents4's ownership: the all-positive family owns the six
 * difference classes and the families flipping coordinates 0, 1 and 2 own the
 * sum classes whose smaller coordinate they flip. Each crossing then evaluates
 * a fixed class list with OctetEvents4's arithmetic.
 *
 * A class's lines lie in its owner's planes. A ray meets a line in the plane
 * of its crossing Q no closer than |n . v| times the line's distance from Q,
 * and that distance is at least the residual length of the class's two free
 * coordinates, so a crossing whose every owned class fails that bound against
 * the support is skipped before the class search.
 */
__attribute__((always_inline)) inline Sample
trace_octet_4d(const math::Vector &direction, const Prepared &prepared) {
  struct Pair {
    uint8_t i, j, k, l;
  };
  const auto &camera = prepared.camera;
  const auto &projection = prepared.octet4_projection;
  const auto &footprint = prepared.footprint;
  const float NEAR = camera.interval.near;
  const float FAR = camera.interval.far;
  const float INVERSE_SCALE = projection.inverse_scale;
  const float SCALE =
      SDF::OctetFramework4::HALF_CUBE * prepared.octet4.cell_size;
  const float WIRE_RADIUS = prepared.octet4.wire_radius;
  std::array<float, 4> ambient;
  for (int axis = 0; axis < 4; ++axis)
    ambient[axis] = math::dot(direction, projection.embedding[axis]);
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
  std::array<float, 4> components, rates, origins;
  for (int m = 0; m < 4; ++m) {
    const uint8_t AXIS = order[m];
    const float LOCAL = ambient[AXIS] * INVERSE_SCALE;
    const float START = projection.origin[AXIS] + camera.radial_start * LOCAL;
    components[m] = fabsf(ambient[AXIS]);
    rates[m] = fabsf(LOCAL);
    origins[m] = ambient[AXIS] < 0.0f ? -START : START;
  }
  // Plane coordinates and speeds of the all-positive family (index 0) and
  // the families flipping coordinates 0, 1 and 2, with one division for the
  // four inverse speeds when none is zero.
  std::array<float, 4> positions, speeds, inverses;
  {
    float position = 0.0f, speed = 0.0f;
    std::array<float, 4> near;
    for (int m = 0; m < 4; ++m) {
      near[m] = origins[m] + NEAR * rates[m];
      position += near[m];
      speed += rates[m];
    }
    positions[0] = 0.5f * position;
    speeds[0] = 0.5f * speed;
    for (int m = 0; m < 3; ++m) {
      positions[m + 1] = positions[0] - near[m];
      speeds[m + 1] = speeds[0] - rates[m];
    }
    const float PRODUCT = speeds[0] * speeds[1] * speeds[2] * speeds[3];
    if (PRODUCT > 1e-30f) {
      const float INVERSE = 1.0f / PRODUCT;
      const float FIRST_PAIR = speeds[0] * speeds[1];
      const float SECOND_PAIR = speeds[2] * speeds[3];
      inverses = {
          INVERSE * speeds[1] * SECOND_PAIR, INVERSE * speeds[0] * SECOND_PAIR,
          INVERSE * FIRST_PAIR * speeds[3], INVERSE * FIRST_PAIR * speeds[2]};
    } else {
      for (int f = 0; f < 4; ++f)
        inverses[f] = 1.0f / speeds[f];
    }
  }
  const float SUPPORT_RATE = .5f * footprint.angular_radius * INVERSE_SCALE;
  const float SUPPORT_BASE =
      WIRE_RADIUS * INVERSE_SCALE + SUPPORT_RATE * footprint.radial_start;
  CoveredCrossings covered;
  int crossings = 0;
  const int BUDGET =
      std::min(prepared.limits.max_candidates, CoveredCrossings::CAPACITY);
  const auto walk = [&]<size_t CLASSES>(
                        int flipped, float sign,
                        const std::array<Pair, CLASSES>
                            &pairs) __attribute__((always_inline)) {
    std::array<float, CLASSES> transverses, denominators, dks, dls;
    bool owned = false;
    for (size_t c = 0; c < CLASSES; ++c) {
      const auto &PAIR = pairs[c];
      const float ALONG = components[PAIR.i] + sign * components[PAIR.j];
      denominators[c] = 1.0f - 0.5f * ALONG * ALONG;
      transverses[c] = 0.5f * (components[PAIR.i] - sign * components[PAIR.j]);
      dks[c] = components[PAIR.k];
      dls[c] = components[PAIR.l];
      owned |= denominators[c] > 0.0f;
    }
    if (!owned)
      return true;
    const float position = positions[flipped + 1];
    const float speed = speeds[flipped + 1];
    if (speed == 0.0f)
      return true;
    const float PLANE = speed > 0.0f ? ceilf(position) : floorf(position);
    const float INVERSE = inverses[flipped + 1];
    const float STEP = fabsf(INVERSE);
    const float FIRST = NEAR + (PLANE - position) * INVERSE;
    if (!Raycast::finite(FIRST) || !Raycast::finite(STEP) || !(STEP > 0.0f))
      return true;
    // Squared |n . v| over the unit ray, with a margin for rounding.
    const float PLANE_SHARE = 0.9999f * (speed * SCALE) * (speed * SCALE);
    // The plane bound runs in 16-bit fixed point: the low 16 bits of Q * 2^16
    // wrap exactly as Q mod 1, so their signed value is Q's residual to the
    // nearest integer, and a wrapped sum or difference of two is a class's
    // across residual. Quantizing the start and step, and float drift in t,
    // move a residual by under FIXED_ERROR over the candidate budget, which
    // loosens each squared term by at most that much.
    constexpr float FIXED_ONE = 65536.0f;
    constexpr float FIXED_ERROR = 1.5e-3f;
    const float THRESHOLD_SCALE = FIXED_ONE * FIXED_ONE / PLANE_SHARE;
    const float THRESHOLD_MARGIN = 3 * FIXED_ERROR * FIXED_ONE * FIXED_ONE;
    std::array<uint32_t, 4> starts, advances;
    for (int m = 0; m < 4; ++m) {
      starts[m] = static_cast<uint32_t>(static_cast<int32_t>(
          roundf((origins[m] + rates[m] * FIRST) * FIXED_ONE)));
      advances[m] = static_cast<uint32_t>(
          static_cast<int32_t>(roundf(rates[m] * STEP * FIXED_ONE)));
    }
    // The first plane lies at or past NEAR. Counting the crossings up front
    // keeps float compares out of the walk; one within rounding of FAR
    // carries no fog-weighted opacity either way.
    const int COUNT =
        FIRST <= FAR ? static_cast<int>((FAR - FIRST) * fabsf(speed)) + 1 : 0;
    crossings += COUNT;
    if (crossings > BUDGET)
      return false;
    float t = FIRST - STEP;
    for (uint32_t index = 0; index < static_cast<uint32_t>(COUNT); ++index) {
      t += STEP;
      const float SUPPORT = SUPPORT_BASE + SUPPORT_RATE * t;
      const float SUPPORT2 = SUPPORT * SUPPORT;
      std::array<uint32_t, 4> fixed;
      std::array<uint32_t, 4> squares;
      for (int m = 0; m < 4; ++m) {
        fixed[m] = starts[m] + index * advances[m];
        const int32_t RESIDUAL =
            static_cast<int16_t>(static_cast<uint16_t>(fixed[m]));
        squares[m] = static_cast<uint32_t>(RESIDUAL * RESIDUAL);
      }
      // Each class's line lies no nearer than half its squared across
      // residual plus its free coordinates' squared residuals.
      std::array<uint32_t, CLASSES> bounds;
      uint32_t nearest = UINT32_MAX;
      for (size_t c = 0; c < CLASSES; ++c) {
        const auto &PAIR = pairs[c];
        const int32_t ACROSS = static_cast<int16_t>(
            static_cast<uint16_t>(sign > 0.0f ? fixed[PAIR.i] - fixed[PAIR.j]
                                              : fixed[PAIR.i] + fixed[PAIR.j]));
        bounds[c] = (static_cast<uint32_t>(ACROSS * ACROSS) >> 1) +
                    squares[PAIR.k] + squares[PAIR.l];
        nearest = std::min(nearest, bounds[c]);
      }
      // Clamped below 2^32; any threshold at or past 2^31 rejects nothing.
      const uint32_t THRESHOLD = static_cast<uint32_t>(
          fminf(SUPPORT2 * THRESHOLD_SCALE + THRESHOLD_MARGIN, 4294967040.0f));
      if (nearest > THRESHOLD)
        continue;
      std::array<float, 4> residual;
      std::array<float, 4> rounded;
      for (int m = 0; m < 4; ++m) {
        const float Q = origins[m] + rates[m] * t;
        rounded[m] = roundf(Q);
        residual[m] = Q - rounded[m];
      }
      int parity = 0;
      for (int m = 0; m < 4; ++m)
        parity += static_cast<int>(rounded[m]);
      const bool ODD = parity & 1;
      float numerator = INFINITY;
      float denominator = 1.0f;
      for (size_t c = 0; c < CLASSES; ++c) {
        const auto &PAIR = pairs[c];
        float rk = residual[PAIR.k];
        float rl = residual[PAIR.l];
        // A class the plane bound clears cannot reach the support, and
        // neither can a farther class it would otherwise have beaten.
        if (!(denominators[c] > 0.0f) || bounds[c] > THRESHOLD)
          continue;
        float across = residual[PAIR.i] - sign * residual[PAIR.j];
        const float SHIFT = across > 0.5f    ? 1.0f
                            : across < -0.5f ? -1.0f
                                             : 0.0f;
        across -= SHIFT;
        if ((SHIFT != 0.0f) != ODD) {
          const float COST = 0.5f - fabsf(across);
          const float COST_K = 1.0f - 2.0f * fabsf(rk);
          const float COST_L = 1.0f - 2.0f * fabsf(rl);
          if (COST <= COST_K && COST <= COST_L)
            across -= copysignf(1.0f, across);
          else if (COST_K <= COST_L)
            rk -= copysignf(1.0f, rk);
          else
            rl -= copysignf(1.0f, rl);
        }
        const float OFFSET2 = 0.5f * across * across + rk * rk + rl * rl;
        const float DOT = across * transverses[c] + rk * dks[c] + rl * dls[c];
        const float N = fmaxf(0.0f, OFFSET2 * denominators[c] - DOT * DOT);
        if (N * denominator < numerator * denominators[c]) {
          numerator = N;
          denominator = denominators[c];
        }
      }
      if (numerator > SUPPORT2 * denominator)
        continue;
      const float FIELD = SCALE * sqrtf(numerator / denominator) - WIRE_RADIUS;
      const float WIDTH = footprint.at(t);
      const float COVERAGE = WIDTH > 0.0f
                                 ? hs::clamp(0.5f - FIELD / WIDTH, 0.0f, 1.0f)
                                 : (FIELD <= 0.0f ? 1.0f : 0.0f);
      if (COVERAGE > 0.0f)
        covered.insert(t, COVERAGE);
    }
    return true;
  };
  constexpr std::array<Pair, 6> DIFFERENCES{{{0, 1, 2, 3},
                                             {0, 2, 1, 3},
                                             {0, 3, 1, 2},
                                             {1, 2, 0, 3},
                                             {1, 3, 0, 2},
                                             {2, 3, 0, 1}}};
  constexpr std::array<Pair, 3> FIRST_SUMS{
      {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}}};
  constexpr std::array<Pair, 2> SECOND_SUMS{{{1, 2, 0, 3}, {1, 3, 0, 2}}};
  constexpr std::array<Pair, 1> THIRD_SUM{{{2, 3, 0, 1}}};
  if (!walk(-1, -1.0f, DIFFERENCES) || !walk(0, 1.0f, FIRST_SUMS) ||
      !walk(1, 1.0f, SECOND_SUMS) || !walk(2, 1.0f, THIRD_SUM)) {
    const SDF::OctetEvents4 events(prepared.octet4, projection, direction,
                                   camera.radial_start, NEAR, footprint);
    return trace(events, camera.interval, prepared.limits, prepared.appearance);
  }
  return covered.composite(prepared.limits, prepared.appearance);
}

/** @brief shade() for a valid octet frame of the matching domain. */
template <bool SLICE_4D>
__attribute__((always_inline)) inline Sample
shade_octet(const math::Vector &direction, const Prepared &prepared) {
  if constexpr (SLICE_4D) {
    return trace_octet_4d(direction, prepared);
  } else {
    return trace_octet_3d(direction, prepared);
  }
}

template <bool SLICE_4D>
HS_HOT_FLASH_MEMBER Sample shade(const math::Vector &direction,
                                 const Prepared &prepared) {
  const auto &camera = prepared.camera;
  if (!prepared.valid ||
      (camera.domain == Raycast::SamplingDomain::SLICE_4D) != SLICE_4D)
    return {{}, Raycast::TraceStatus::INVALID_QUERY};
  if (prepared.geometry != Geometry::OCTET) {
    if (prepared.geometry == Geometry::SHELLS) {
      if constexpr (SLICE_4D) {
        if (prepared.periodic_shells.slice_march) {
          const auto SHELLS = SDF::trace_periodic_shells_4d(
              prepared.periodic_shells, camera, direction, prepared.limits,
              prepared.appearance);
          return {SHELLS.color, SHELLS.status};
        }
      } else if (prepared.periodic_shells.single_owner) {
        const auto SHELLS = SDF::trace_periodic_shells_3d(
            prepared.periodic_shells, camera, direction, prepared.limits,
            prepared.appearance);
        return {SHELLS.color, SHELLS.status};
      }
    }
    Raycast::ShadedTrace sample;
    if (prepared.geometry == Geometry::AFFINE_CUBIC)
      sample = SDF::shade_affine_lattice(
          camera, direction, prepared.octet.cell_size,
          prepared.octet.wire_radius, prepared.footprint, prepared.limits,
          prepared.appearance, prepared.shear, prepared.stretch);
    else if (prepared.geometry == Geometry::SHELLS)
      sample = SDF::shade_periodic_shells(prepared.periodic_shells, camera,
                                          direction, prepared.limits,
                                          prepared.appearance);
    else {
      sample = SDF::CellularWire::shade(
          *prepared.cellular, prepared.octet.cell_size,
          prepared.octet.wire_radius, camera, prepared.footprint,
          prepared.limits, prepared.appearance, direction);
    }
    return {sample.color.color * sample.color.alpha, sample.trace.status};
  }
  return shade_octet<SLICE_4D>(direction, prepared);
}

} // namespace HyperLatticeExperimental
