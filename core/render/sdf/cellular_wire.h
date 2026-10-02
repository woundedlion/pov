/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file cellular_wire.h
 * @brief Cellular wire geometry and ray shading. */

#include "render/ray/camera.h"
#include "render/ray/shade.h"

namespace SDF::CellularWire {

enum class Kind { DIAMOND, HEXAGONAL, RHOMBIC };

constexpr math::Vector periods(Kind kind, float cell_size) {
  return kind == Kind::HEXAGONAL
             ? math::Vector(3 * cell_size, 1.7320508075688772f * cell_size,
                            cell_size)
             : math::Vector(cell_size, cell_size, cell_size);
}

struct Edge {
  math::Vector a;
  math::Vector b;
};

struct Hit {
  float t;
  float coverage;
  uint32_t feature;
};

/** @brief Caller-owned crossing scratch for one cellular ray. */
struct HitStorage {
  std::array<Hit, 32> hits;
};

/** @brief Unique edges owned by their midpoint's rectangular translation cell. */
struct Geometry {
  std::array<Edge, 32> edges{};
  int count = 0;
  math::Vector period{};
  math::Vector lower{};
  math::Vector upper{};

  Geometry() = default;

  constexpr explicit Geometry(Kind kind) : period(periods(kind, 1)) {
    const auto add = [&](math::Vector a, math::Vector b) {
      const math::Vector MID = (a + b) * .5f;
      const auto tile = [](float value) {
        const int WHOLE = static_cast<int>(value);
        return static_cast<float>(WHOLE - (value < WHOLE ? 1 : 0));
      };
      const math::Vector OFFSET(tile(MID.x / period.x) * period.x,
                                tile(MID.y / period.y) * period.y,
                                tile(MID.z / period.z) * period.z);
      a = a - OFFSET;
      b = b - OFFSET;
      for (int i = 0; i < count; ++i) {
        const auto equal = [](const math::Vector &u, const math::Vector &v) {
          return math::dot(u - v, u - v) < 1e-10f;
        };
        if ((equal(a, edges[i].a) && equal(b, edges[i].b)) ||
            (equal(a, edges[i].b) && equal(b, edges[i].a)))
          return;
      }
      HS_CHECK(count < static_cast<int>(edges.size()),
               "cellular wire edge capacity");
      edges[count++] = {a, b};
    };
    const math::Vector FCC[] = {
        {0, 0, 0}, {0, .5f, .5f}, {.5f, 0, .5f}, {.5f, .5f, 0}};
    if (kind == Kind::DIAMOND) {
      const math::Vector BONDS[] = {{.25f, .25f, .25f},
                                    {.25f, -.25f, -.25f},
                                    {-.25f, .25f, -.25f},
                                    {-.25f, -.25f, .25f}};
      for (const auto &site : FCC)
        for (const auto &bond : BONDS)
          add(site, site + bond);
    } else if (kind == Kind::RHOMBIC) {
      for (const auto &site : FCC)
        for (int x = -1; x <= 1; x += 2)
          for (int y = -1; y <= 1; y += 2)
            for (int z = -1; z <= 1; z += 2) {
              const math::Vector CORNER = site + math::Vector(x, y, z) * .25f;
              add(CORNER, site + math::Vector(x * .5f, 0, 0));
              add(CORNER, site + math::Vector(0, y * .5f, 0));
              add(CORNER, site + math::Vector(0, 0, z * .5f));
            }
    } else {
      constexpr float H = .8660254037844386f;
      const math::Vector VERTICES[] = {{1, 0, 0},  {.5f, H, 0},   {-.5f, H, 0},
                                       {-1, 0, 0}, {-.5f, -H, 0}, {.5f, -H, 0}};
      const math::Vector CENTERS[] = {{0, 0, 0}, {1.5f, H, 0}};
      for (const auto &center : CENTERS)
        for (int i = 0; i < 6; ++i) {
          const auto A = center + VERTICES[i];
          add(A, center + VERTICES[(i + 1) % 6]);
          add(A, A + math::Vector(0, 0, 1));
        }
    }
    for (int i = 0; i < count; ++i)
      for (const auto &p : {edges[i].a, edges[i].b}) {
        lower = {std::min(lower.x, p.x), std::min(lower.y, p.y),
                 std::min(lower.z, p.z)};
        upper = {std::max(upper.x, p.x), std::max(upper.y, p.y),
                 std::max(upper.z, p.z)};
      }
  }
};

inline constexpr Geometry DIAMOND_GEOMETRY HS_PROGMEM_UNIQUE(CELLULAR_DIAMOND){
    Kind::DIAMOND};
inline constexpr Geometry
    HEXAGONAL_GEOMETRY HS_PROGMEM_UNIQUE(CELLULAR_HEXAGONAL){Kind::HEXAGONAL};
inline constexpr Geometry RHOMBIC_GEOMETRY HS_PROGMEM_UNIQUE(CELLULAR_RHOMBIC){
    Kind::RHOMBIC};

inline const Geometry &geometry(Kind kind) {
  switch (kind) {
  case Kind::DIAMOND:
    return DIAMOND_GEOMETRY;
  case Kind::HEXAGONAL:
    return HEXAGONAL_GEOMETRY;
  case Kind::RHOMBIC:
    return RHOMBIC_GEOMETRY;
  }
  return DIAMOND_GEOMETRY;
}

HS_O3_BEGIN

/** @brief Closest points of a finite ray interval and a finite round strut. */
inline float closest(const Raycast::Ray &ray, const Edge &edge, float &t) {
  const auto V = edge.b - edge.a;
  const auto W = ray.origin - edge.a;
  const float LENGTH2 = math::dot(V, V);
  const float DV = math::dot(ray.direction, V);
  const float DW = math::dot(ray.direction, W);
  const float VW = math::dot(V, W);
  const float DENOM = LENGTH2 - DV * DV;
  float s = DENOM > 1e-8f * LENGTH2
                ? hs::clamp((VW - DV * DW) / DENOM, 0.0f, 1.0f)
                : 0.0f;
  t = hs::clamp(DV * s - DW, ray.interval.near, ray.interval.far);
  s = hs::clamp((VW + t * DV) / LENGTH2, 0.0f, 1.0f);
  t = hs::clamp(DV * s - DW, ray.interval.near, ray.interval.far);
  const auto DELTA = ray.at(t) - (edge.a + V * s);
  return sqrtf(fmaxf(0.0f, math::dot(DELTA, DELTA)));
}

struct BoxRay {
  float origin[3];
  float inverse[3];

  explicit BoxRay(const Raycast::Ray &ray)
      : origin{ray.origin.x, ray.origin.y, ray.origin.z} {
    const float D[] = {ray.direction.x, ray.direction.y, ray.direction.z};
    for (int axis = 0; axis < 3; ++axis)
      inverse[axis] = fabsf(D[axis]) > 1e-12f ? 1.0f / D[axis] : 0.0f;
  }
};

inline bool box_overlap(const BoxRay &ray, const math::Vector &lo,
                        const math::Vector &hi, float near, float far) {
  const float L[] = {lo.x, lo.y, lo.z};
  const float H[] = {hi.x, hi.y, hi.z};
  for (int axis = 0; axis < 3; ++axis) {
    const float INVERSE = ray.inverse[axis];
    if (INVERSE == 0.0f) {
      if (ray.origin[axis] < L[axis] || ray.origin[axis] > H[axis])
        return false;
    } else {
      float a = (L[axis] - ray.origin[axis]) * INVERSE;
      float b = (H[axis] - ray.origin[axis]) * INVERSE;
      if (a > b)
        std::swap(a, b);
      near = fmaxf(near, a);
      far = fminf(far, b);
      if (near > far)
        return false;
    }
  }
  return true;
}

/** @brief Bounded analytic strut traversal; geometry is prepared once per frame. */
HS_HOT_FLASH_MEMBER inline Raycast::ShadedTrace
shade(const Geometry &geometry, float cell_size, float wire_radius,
      const Raycast::PreparedCamera &camera,
      const Raycast::Footprint &footprint, const Raycast::TraceLimits &limits,
      const Raycast::Appearance &appearance, HitStorage &storage,
      const math::Vector &direction) {
  Raycast::ShadedTrace result{};
  auto &trace = result.trace;
  const auto WORLD_ORIGIN = camera.point3(direction * camera.radial_start);
  const auto &M = camera.embedding.m;
  const math::Vector WORLD_DIRECTION(
      M[0][0] * direction.x + M[0][1] * direction.y + M[0][2] * direction.z,
      M[1][0] * direction.x + M[1][1] * direction.y + M[1][2] * direction.z,
      M[2][0] * direction.x + M[2][1] * direction.y + M[2][2] * direction.z);
  Raycast::Ray ray{WORLD_ORIGIN, WORLD_DIRECTION, camera.interval};
  const float REQUESTED_SUPPORT =
      wire_radius + .5f * footprint.at(camera.interval.far);
  if (camera.domain != Raycast::SamplingDomain::SPATIAL_3D || !camera.valid() ||
      !ray.valid() || !Raycast::finite(cell_size) || cell_size <= 0 ||
      !Raycast::finite(wire_radius) || wire_radius <= 0 ||
      !Raycast::finite(REQUESTED_SUPPORT) || geometry.count <= 0 ||
      !Raycast::finite(footprint.angular_radius) ||
      footprint.angular_radius < 0 ||
      !Raycast::finite(footprint.radial_start) || footprint.radial_start < 0 ||
      !appearance.palette) {
    trace.status = Raycast::TraceStatus::INVALID_QUERY;
    return result;
  }
  const bool TRUNCATED = REQUESTED_SUPPORT > .49f * cell_size;
  if (TRUNCATED) {
    ray.interval.far =
        footprint.angular_radius > 0
            ? 2 * (.49f * cell_size - wire_radius) / footprint.angular_radius -
                  footprint.radial_start
            : ray.interval.near;
    if (ray.interval.far <= ray.interval.near) {
      trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      return result;
    }
  }
  const float SUPPORT = wire_radius + .5f * footprint.at(ray.interval.far);
  const auto PERIOD = geometry.period * cell_size;
  const float P[] = {PERIOD.x, PERIOD.y, PERIOD.z};
  const float O[] = {ray.origin.x, ray.origin.y, ray.origin.z};
  const float D[] = {ray.direction.x, ray.direction.y, ray.direction.z};
  float cell[3];
  float next[3];
  float delta[3];
  for (int axis = 0; axis < 3; ++axis) {
    cell[axis] = floorf((O[axis] + D[axis] * ray.interval.near) / P[axis]);
    const float BOUND = (cell[axis] + (D[axis] >= 0 ? 1 : 0)) * P[axis];
    next[axis] =
        fabsf(D[axis]) > 1e-12f ? (BOUND - O[axis]) / D[axis] : INFINITY;
    delta[axis] = fabsf(D[axis]) > 1e-12f ? P[axis] / fabsf(D[axis]) : INFINITY;
  }
  const BoxRay BOX_RAY(ray);
  LayerComposite composite;
  float start = ray.interval.near;
  float last_t = -INFINITY;
  while (start < ray.interval.far) {
    if (trace.counters.steps >= limits.max_steps) {
      trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++trace.counters.steps;
    const float END =
        fminf(fminf(next[0], next[1]), fminf(next[2], ray.interval.far));
    auto &hits = storage.hits;
    int count = 0;
    for (int x = -1; x <= 1; ++x)
      for (int y = -1; y <= 1; ++y)
        for (int z = -1; z <= 1; ++z) {
          const math::Vector OFFSET((cell[0] + x) * P[0], (cell[1] + y) * P[1],
                                    (cell[2] + z) * P[2]);
          const math::Vector PAD(SUPPORT, SUPPORT, SUPPORT);
          if (!box_overlap(BOX_RAY, OFFSET + geometry.lower * cell_size - PAD,
                           OFFSET + geometry.upper * cell_size + PAD, start,
                           END))
            continue;
          for (int i = 0; i < geometry.count; ++i) {
            const Edge EDGE{geometry.edges[i].a * cell_size + OFFSET,
                            geometry.edges[i].b * cell_size + OFFSET};
            const math::Vector LO(fminf(EDGE.a.x, EDGE.b.x),
                                  fminf(EDGE.a.y, EDGE.b.y),
                                  fminf(EDGE.a.z, EDGE.b.z));
            const math::Vector HI(fmaxf(EDGE.a.x, EDGE.b.x),
                                  fmaxf(EDGE.a.y, EDGE.b.y),
                                  fmaxf(EDGE.a.z, EDGE.b.z));
            if (!box_overlap(BOX_RAY, LO - PAD, HI + PAD, start, END))
              continue;
            if (trace.counters.candidates >= limits.max_candidates) {
              trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
              goto shade_hits;
            }
            ++trace.counters.candidates;
            float t;
            const float DISTANCE = closest(ray, EDGE, t);
            if (t < start || t >= END)
              continue;
            const float AA = footprint.at(t);
            const float COVERAGE =
                AA > 0
                    ? hs::clamp(.5f - (DISTANCE - wire_radius) / AA, 0.0f, 1.0f)
                    : (DISTANCE <= wire_radius ? 1.0f : 0.0f);
            if (COVERAGE <= 0)
              continue;
            if (count == static_cast<int>(hits.size())) {
              trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
              goto shade_hits;
            }
            hits[count++] = {t, COVERAGE, static_cast<uint32_t>(i)};
          }
        }
  shade_hits:
    std::sort(hits.begin(), hits.begin() + count,
              [](const auto &a, const auto &b) { return a.t < b.t; });
    for (int i = 0; i < count; ++i) {
      auto hit = hits[i];
      while (i + 1 < count &&
             hits[i + 1].t - hit.t <= limits.position_tolerance)
        hit.coverage = fmaxf(hit.coverage, hits[++i].coverage);
      if (hit.t - last_t <= limits.position_tolerance)
        continue;
      if (trace.counters.layers >= limits.max_layers) {
        trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
        goto finished;
      }
      ++trace.counters.layers;
      trace.contribution = {hit.t, hit.coverage, 0,     hit.feature,
                            0,     false,        false, {}};
      last_t = hit.t;
      appearance.composite(composite, hit.t, hit.coverage);
      if (composite.saturated()) {
        trace.status = Raycast::TraceStatus::SATURATED;
        goto finished;
      }
    }
    if (trace.status == Raycast::TraceStatus::BUDGET_EXHAUSTED)
      break;
    if (END >= ray.interval.far)
      break;
    for (int axis = 0; axis < 3; ++axis)
      if (next[axis] <= END) {
        cell[axis] += D[axis] >= 0 ? 1 : -1;
        next[axis] += delta[axis];
      }
    start = END;
  }
finished:
  if (TRUNCATED && trace.status == Raycast::TraceStatus::RANGE_COMPLETE)
    trace.status = Raycast::TraceStatus::BUDGET_EXHAUSTED;
  result.color = composite.finish();
  return result;
}

HS_O3_END

} // namespace SDF::CellularWire
