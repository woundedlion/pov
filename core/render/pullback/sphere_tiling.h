/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_tiling.h @brief Radially projected solid faces and edges. */

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <new>
#include <numbers>
#include <span>
#include <utility>
#include "core/memory.h"
#include "core/mesh/solids.h"
#include "core/pattern/sphere_lut.h"
#include "core/platform/cycle_counting.h"
#include "core/render/pullback/sphere_tiling_edges.h"
#include "core/render/sdf/face.h"

namespace Pullback::Source {

/**
 * @brief Cubemap-indexed spherical faces and edges of a registry solid, built
 *        into a caller-owned block.
 * @details build() generates the mesh in the block's free space, keeps its
 * compact copy and every other transient in two scratch arenas, then writes
 * the vertices, edges, coarse and narrow edge indices and, optionally, the
 * face records and face index into the block. The geometry views the block
 * and is valid while the block is.
 */
class SphereTilingGeometry {
public:
  static constexpr uint16_t NO_FACE = UINT16_MAX; ///< Uncovered direction.

  /// Records a build writes.
  enum class Content : uint8_t {
    EDGES,          ///< Edge arcs and their indices only.
    FACES_AND_EDGES ///< Edges plus face polygons and the face index.
  };

  /// Edge indices a build writes.
  enum class Indices : uint8_t {
    COARSE, ///< The coarse index, which serves every band.
    NARROW, ///< The narrow index alone, or the coarse one without it.
    BOTH    ///< Both indices.
  };

  /** @brief Inward angular distances from a face boundary. */
  struct FaceDepth {
    float sample;      ///< Distance at the sampled direction, in radians.
    float center;      ///< Distance at the face center, in radians.
    float across_seam; ///< Center depth on the other side of a crossed seam.
  };

  /** @brief Coarse edge-index resolution for a solid.
   * @param edge_count Number of distinct mesh edges.
   * @return Cells along one cube-face axis. */
  static constexpr int edge_steps(size_t edge_count) {
    return edge_count <= 180 ? 16 : edge_count <= 540 ? 20 : 24;
  }

  /** @brief Whether a solid also gets the thin-band narrow index.
   * @param edge_count Number of distinct mesh edges.
   * @return True above the edge count where the narrow index pays off. */
  static constexpr bool narrow_indexed(size_t edge_count) {
    return edge_count > 180;
  }

  /** @brief Builds a registry solid's geometry into a block.
   * @tparam CONTENT Records to build.
   * @param registry_index Index into the solids registry.
   * @param block Destination storage; also hosts mesh generation.
   * @param scratch_a Transient storage for the compact mesh and index rows.
   * @param scratch_b Transient storage for the topology and edge pairs.
   * @param indices Edge indices to build. A solid without a narrow index
   *        gets the coarse one instead.
   * @return Geometry viewing @p block. Traps when the block is too small. */
  template <Content CONTENT = Content::EDGES>
  static HS_FLASH_MEMBER SphereTilingGeometry
  build(size_t registry_index, std::span<uint8_t> block, Arena &scratch_a,
        Arena &scratch_b, Indices indices = Indices::BOTH) {
    ScratchScope keep_a(scratch_a);
    ScratchScope keep_b(scratch_b);
    Block output(block);
    SphereTilingGeometry geometry;
    const size_t topology_mark = scratch_b.get_offset();
    MeshCopy mesh;
    {
      HS_PROFILE(tiling_mesh);
      mesh = copy_mesh(Solids::get_entry(registry_index), output, scratch_a,
                       scratch_b);
    }
    std::span<const EdgeKey> edge_keys;
    {
      HS_PROFILE(tiling_edges);
      geometry.build_vertices(mesh, output);
      if constexpr (CONTENT == Content::EDGES)
        scratch_a.set_offset(mesh.positions_mark);
      edge_keys =
          unique_edges(mesh, CONTENT == Content::EDGES ? scratch_a : scratch_b);
      geometry.build_edges(edge_keys, output);
    }
    if constexpr (CONTENT == Content::FACES_AND_EDGES) {
      geometry.build_faces(mesh, edge_keys, output);
      geometry.build_face_index(output);
    }
    const float reach = geometry.face_center_reach(mesh, edge_keys);
    scratch_a.set_offset(mesh.positions_mark);
    scratch_b.set_offset(topology_mark);
    const auto coarse = [&] {
      HS_PROFILE(tiling_coarse);
      geometry.build_edge_index(reach, output, scratch_a, scratch_b);
    };
    const auto narrow = [&] {
      if (narrow_indexed(geometry.edges.size())) {
        HS_PROFILE(tiling_narrow);
        geometry.build_narrow_index(output, scratch_a, scratch_b);
      }
    };
    if (indices == Indices::NARROW) {
      narrow();
      if (geometry.narrow_cells.empty())
        coarse();
    } else {
      coarse();
      if (indices == Indices::BOTH)
        narrow();
    }
    geometry.used = output.used();
    return geometry;
  }

  /** @brief Block bytes the geometry occupies.
   * @return Offset of the last record's end from the block start. */
  size_t bytes() const { return used; }

  /** @brief This geometry viewing a copy of its block.
   * @param from Start of the block the geometry was built in.
   * @param to Start of a byte-identical copy of that block.
   * @return Geometry viewing @p to. */
  SphereTilingGeometry rebased(const uint8_t *from, uint8_t *to) const {
    SphereTilingGeometry copy = *this;
    const auto move = [&]<typename T>(std::span<const T> &view) {
      if (!view.empty())
        view = {
            reinterpret_cast<const T *>(
                to + (reinterpret_cast<const uint8_t *>(view.data()) - from)),
            view.size()};
    };
    move(copy.vertices);
    move(copy.edges);
    move(copy.edge_offsets);
    move(copy.edge_ids);
    move(copy.narrow_cells);
    move(copy.narrow_candidates);
    move(copy.faces);
    move(copy.corners);
    move(copy.edge_faces);
    move(copy.face_offsets);
    move(copy.face_ids);
    return copy;
  }

  /** @brief Finds the nearest projected face containing a unit direction.
   * @param unit_direction Unit direction in the solid's frame.
   * @return Face ID, or NO_FACE when uncovered or built without faces. */
  uint16_t face_at(const math::Vector &unit_direction) const {
    if (face_offsets.empty())
      return NO_FACE;
    const int cell = tiling_edges_detail::cell_at(unit_direction, FACE_STEPS);
    uint16_t owner = NO_FACE;
    float nearest = 0.0f;
    for (uint32_t i = face_offsets[cell]; i < face_offsets[cell + 1]; ++i) {
      const uint16_t id = face_ids[i];
      const Face &face = faces[id];
      const float score = math::dot(face.selector, unit_direction);
      if (score > nearest && contains(face, unit_direction)) {
        owner = id;
        nearest = score;
      }
    }
    return owner;
  }

  /** @brief Measures inward depth within a projected face.
   * @param face Face ID returned by face_at(), excluding NO_FACE.
   * @param unit_direction Unit direction in that face.
   * @param filter_radius Angular radius across the nearest seam, in radians.
   * @return Sample depth and each covered face's center depth. */
  FaceDepth face_depth(uint16_t face, const math::Vector &unit_direction,
                       float filter_radius) const {
    const Face &polygon = faces[face];
    const auto [distance, seam] = boundary_at(polygon, unit_direction);
    float across = polygon.center_depth;
    if (distance < filter_radius) {
      const auto [first, second] = edge_faces[seam];
      const uint16_t neighbor = first == face ? second : first;
      if (neighbor != NO_FACE)
        across = faces[neighbor].center_depth;
    }
    return {distance, polygon.center_depth, across};
  }

  /** @brief Views the edge geometry with its indices.
   * @return Indexed projected edges owned by the block. */
  SphereTilingEdgeView edge_view() const {
    return {vertices, edges, edge_offsets,
            edge_ids, steps, {narrow_cells, narrow_candidates}};
  }

private:
  static constexpr int FACE_STEPS = 12;
  static constexpr int FACE_CELLS = 6 * FACE_STEPS * FACE_STEPS;
  static constexpr float SWEEP_MARGIN = 1.0e-3f;
  static constexpr float SPAN_SLACK = 1.0e-3f;
  static constexpr float PRUNE_SCALE = 1.0001f;
  static constexpr float PRUNE_FLOOR = 1.0e-7f;
  static constexpr float SERIES_TOLERANCE = 1.0e-5f;
  static constexpr float SERIES_CHORD = 0.25f;
  static constexpr float SERIES_ANGLE = 1.0f;
  static constexpr float CHART_EPSILON = 1.0e-6f;
  static constexpr float CAP_ANGLE_PAD = 1.0e-4f;
  static constexpr float NARROW_PAD = 1.0e-4f;

  struct Point {
    float x;
    float y;
    uint16_t edge;
  };

  struct Face {
    math::Vector axis;
    math::Vector u;
    math::Vector v;
    math::Vector selector;
    float radius;
    float center_depth;
    size_t first;
    uint16_t count;
  };

  using EdgeIds = std::pair<uint16_t, uint16_t>;
  // Lower vertex in the high half, upper in the low half.
  using EdgeKey = uint32_t;

  static constexpr EdgeKey edge_key(uint16_t a, uint16_t b) {
    return static_cast<EdgeKey>(std::min(a, b)) << 16 | std::max(a, b);
  }

  // Bump allocation over the caller's block; overflow traps.
  class Block {
    uint8_t *start;
    uint8_t *next;
    uint8_t *limit;

  public:
    HS_HOT_INLINE explicit Block(std::span<uint8_t> bytes)
        : start(bytes.data()), next(bytes.data()),
          limit(bytes.data() + bytes.size()) {}

    template <typename T> HS_HOT_INLINE std::span<T> free() const {
      const auto address = reinterpret_cast<uintptr_t>(next);
      const auto aligned = (address + alignof(T) - 1) & ~(alignof(T) - 1);
      const auto end = reinterpret_cast<uintptr_t>(limit);
      if (aligned >= end)
        return {};
      return {reinterpret_cast<T *>(aligned), (end - aligned) / sizeof(T)};
    }

    template <typename T> HS_HOT_INLINE T *take(size_t count) {
      const auto room = free<T>();
      if (count == 0)
        return room.data();
      // no death case: the budget test builds every admitted solid in place.
      HS_CHECK(count <= room.size(),
               "sphere tiling: geometry exceeds its block (%lu > %lu)",
               static_cast<unsigned long>(count),
               static_cast<unsigned long>(room.size()));
      next = reinterpret_cast<uint8_t *>(room.data() + count);
      return room.data();
    }

    uint8_t *end() const { return limit; }
    uint8_t *position() const { return next; }
    void rewind(uint8_t *mark) { next = mark; }
    size_t used() const { return static_cast<size_t>(next - start); }
  };

  // Variable-length record list growing into a block's free space.
  template <typename T> class Tail {
    Block &block;
    std::span<T> room;
    size_t count = 0;

  public:
    HS_HOT_INLINE explicit Tail(Block &target)
        : block(target), room(target.free<T>()) {}

    HS_HOT_INLINE void push(T value) {
      // no death case: the budget test builds every admitted solid in place.
      HS_CHECK(count < room.size(),
               "sphere tiling: geometry exceeds its block (%lu > %lu)",
               static_cast<unsigned long>(count + 1),
               static_cast<unsigned long>(room.size()));
      room[count++] = value;
    }

    T &back() { return room[count - 1]; }
    void pop() { --count; }
    size_t size() const { return count; }
    std::span<T> written() const { return room.first(count); }

    HS_HOT_INLINE std::span<const T> commit() {
      if (count == 0)
        return {};
      return {block.take<T>(count), count};
    }
  };

  struct MeshCopy {
    std::span<const math::Vector> positions;
    std::span<const uint16_t> indices;
    std::span<const uint8_t> counts;
    size_t positions_mark;
  };

  std::span<const math::Vector> vertices;
  std::span<const SphereTilingEdge> edges;
  std::span<const uint32_t> edge_offsets;
  std::span<const uint16_t> edge_ids;
  int steps = 0;
  std::span<const uint16_t> narrow_cells;
  std::span<const uint16_t> narrow_candidates;
  std::span<const Face> faces;
  std::span<const Point> corners;
  std::span<const EdgeIds> edge_faces;
  std::span<const uint32_t> face_offsets;
  std::span<const uint16_t> face_ids;
  size_t used = 0;

  // Generates in the block's free space, which later records overwrite.
  static HS_FLASH_MEMBER MeshCopy copy_mesh(const Solids::Entry &entry,
                                            Block &block, Arena &scratch_a,
                                            Arena &scratch_b) {
    const auto room = block.free<uint8_t>();
    constexpr size_t ALIGN = alignof(std::max_align_t);
    const auto base = reinterpret_cast<uintptr_t>(room.data());
    const uintptr_t top =
        (reinterpret_cast<uintptr_t>(block.end()) - 2 * sizeof(Arena)) &
        ~(alignof(Arena) - 1);
    // no death case: the budget test builds every admitted solid in place.
    HS_CHECK(!room.empty() && top > base + 2 * ALIGN,
             "sphere tiling: no block room for mesh generation");
    const size_t half = ((top - base) / 2) & ~(ALIGN - 1);
    auto *arenas = reinterpret_cast<Arena *>(top);
    Arena *a = ::new (static_cast<void *>(arenas)) Arena(room.data(), half);
    Arena *b = ::new (static_cast<void *>(arenas + 1))
        Arena(room.data() + half, top - base - half);
    MeshCopy copy;
    {
      const PolyMesh mesh = entry.generate(*a, *b);
      // no death case: the tiling test builds every solid in the fixed registry.
      HS_CHECK(mesh.num_faces() < NO_FACE, "sphere tiling: too many faces");
      copy.positions_mark = scratch_a.get_offset();
      auto *positions =
          scratch_a.allocate_n<math::Vector>(mesh.vertices.size());
      std::copy(mesh.vertices.begin(), mesh.vertices.end(), positions);
      auto *indices = scratch_b.allocate_n<uint16_t>(mesh.faces.size());
      std::copy(mesh.faces.begin(), mesh.faces.end(), indices);
      auto *counts = scratch_b.allocate_n<uint8_t>(mesh.face_counts.size());
      std::copy(mesh.face_counts.begin(), mesh.face_counts.end(), counts);
      copy.positions = {positions, mesh.vertices.size()};
      copy.indices = {indices, mesh.faces.size()};
      copy.counts = {counts, mesh.face_counts.size()};
    }
    std::destroy_at(b);
    std::destroy_at(a);
    return copy;
  }

  HS_FLASH_MEMBER __attribute__((flatten)) void
  build_vertices(const MeshCopy &mesh, Block &block) {
    auto *unit = block.take<math::Vector>(mesh.positions.size());
    for (size_t i = 0; i < mesh.positions.size(); ++i)
      unit[i] = mesh.positions[i].normalized();
    vertices = {unit, mesh.positions.size()};
  }

  // Sorted distinct edge keys, in `scratch`.
  static HS_FLASH_MEMBER __attribute__((flatten)) std::span<const EdgeKey>
  unique_edges(const MeshCopy &mesh, Arena &scratch) {
    const size_t start = scratch.get_offset();
    auto *keys = scratch.allocate_n<EdgeKey>(mesh.indices.size());
    size_t offset = 0;
    size_t count = 0;
    for (uint8_t sides : mesh.counts) {
      const auto ids = mesh.indices.subspan(offset, sides);
      uint16_t previous = ids.back();
      for (uint16_t id : ids) {
        keys[count++] = edge_key(previous, id);
        previous = id;
      }
      offset += sides;
    }
    sort_keys(keys, count);
    size_t distinct = 0;
    for (size_t i = 0; i < count; ++i)
      if (distinct == 0 || keys[i] != keys[distinct - 1])
        keys[distinct++] = keys[i];
    // no death case: the tiling test indexes edges for every fixed registry solid.
    HS_CHECK(distinct <= SphereTilingNarrowIndex::INTERIOR,
             "sphere tiling: too many edges");
    scratch.set_offset(start);
    scratch.allocate_n<EdgeKey>(distinct);
    return {keys, distinct};
  }

  // Heapsort: std::sort would instantiate library code outside flash.
  static HS_HOT_INLINE void sort_keys(EdgeKey *keys, size_t count) {
    const auto sift = [keys](size_t root, size_t end) HS_HOT_INLINE {
      for (size_t child = 2 * root + 1; child < end; child = 2 * root + 1) {
        if (child + 1 < end && keys[child] < keys[child + 1])
          ++child;
        if (keys[child] <= keys[root])
          return;
        std::swap(keys[root], keys[child]);
        root = child;
      }
    };
    for (size_t root = count / 2; root-- > 0;)
      sift(root, count);
    for (size_t end = count; end-- > 1;) {
      std::swap(keys[0], keys[end]);
      sift(0, end);
    }
  }

  HS_FLASH_MEMBER __attribute__((flatten)) void
  build_edges(std::span<const EdgeKey> edge_keys, Block &block) {
    auto *arcs = block.take<SphereTilingEdge>(edge_keys.size());
    for (size_t i = 0; i < edge_keys.size(); ++i) {
      const auto first = static_cast<uint16_t>(edge_keys[i] >> 16);
      const auto last = static_cast<uint16_t>(edge_keys[i] & UINT16_MAX);
      const auto &from = vertices[first];
      const auto &to = vertices[last];
      const auto cross = math::cross(from, to);
      const float length = cross.length();
      if (length > 0.0f) {
        const auto normal = cross / length;
        const float arc = atan2f(length, math::dot(from, to));
        arcs[i] = {first,
                   last,
                   normal,
                   math::cross(normal, from),
                   math::cross(to, normal),
                   1.0f / arc};
      } else {
        arcs[i] = {first, last, {}, {}, {}, 0.0f};
      }
    }
    edges = {arcs, edge_keys.size()};
  }

  HS_FLASH_MEMBER void build_faces(const MeshCopy &mesh,
                                   std::span<const EdgeKey> edge_keys,
                                   Block &block) {
    auto *seams = block.take<EdgeIds>(edge_keys.size());
    std::fill_n(seams, edge_keys.size(), EdgeIds{NO_FACE, NO_FACE});
    auto *polygons = block.take<Face>(mesh.counts.size());
    Tail<Point> points(block);
    size_t offset = 0;
    for (size_t face = 0; face < mesh.counts.size(); ++face) {
      const auto ids = mesh.indices.subspan(offset, mesh.counts[face]);
      polygons[face] = build_face(mesh.positions, ids, edge_keys, seams,
                                  static_cast<uint16_t>(face), points);
      offset += mesh.counts[face];
    }
    corners = points.commit();
    faces = {polygons, mesh.counts.size()};
    edge_faces = {seams, edge_keys.size()};
  }

  HS_FLASH_MEMBER Face build_face(std::span<const math::Vector> positions,
                                  std::span<const uint16_t> ids,
                                  std::span<const EdgeKey> edge_keys,
                                  EdgeIds *seams, uint16_t face_id,
                                  Tail<Point> &points) {
    math::Vector center{};
    math::Vector axis{};
    math::Vector normal{};
    const math::Vector &origin = positions[ids.front()];
    math::Vector previous = positions[ids.back()] - origin;
    for (uint16_t id : ids) {
      const auto &v = positions[id];
      center += v;
      axis += vertices[id];
      const math::Vector local = v - origin;
      normal += math::cross(previous, local);
      previous = local;
    }
    center = center / static_cast<float>(ids.size());
    axis = math::normalized_or(axis, vertices[ids.front()]);
    const float normal_length = normal.length();
    normal = normal_length > 0.0f ? normal / normal_length : axis;
    const auto min_depth = [&](const math::Vector &candidate) {
      float depth = 1.0f;
      for (uint16_t id : ids)
        depth = std::min(depth, math::dot(candidate, vertices[id]));
      return depth;
    };
    float depth = min_depth(axis);
    if (depth <= CHART_EPSILON) {
      const math::Vector aligned =
          math::dot(normal, center) < 0.0f ? -normal : normal;
      axis = aligned;
      depth = min_depth(axis);
    }

    const auto count = static_cast<uint16_t>(ids.size());
    Face face{axis, {}, {}, {}, -1.0f, 0.0f, points.size(), count};
    if (depth <= CHART_EPSILON)
      return face;

    const float plane_offset = math::dot(normal, center);
    face.selector = fabsf(plane_offset) > CHART_EPSILON
                        ? normal / plane_offset
                        : axis / math::dot(axis, center);
    face.u = math::tangent_axis(axis);
    face.v = math::cross(axis, face.u);
    face.radius = acosf(std::min(1.0f, depth));
    uint16_t previous_id = ids.back();
    for (uint16_t id : ids) {
      Point point;
      SDF::gnomonic_project(vertices[id], axis, face.u, face.v, point.x,
                            point.y);
      const auto match = std::lower_bound(edge_keys.begin(), edge_keys.end(),
                                          edge_key(previous_id, id));
      point.edge = static_cast<uint16_t>(match - edge_keys.begin());
      auto &[first, second] = seams[point.edge];
      if (first == NO_FACE)
        first = face_id;
      else if (first != face_id && second == NO_FACE)
        second = face_id;
      points.push(point);
      previous_id = id;
    }
    corners = points.written();
    const auto center_direction = math::normalized_or(center, axis);
    face.center_depth = boundary_at(face, center_direction).first;
    return face;
  }

  HS_FLASH_MEMBER void build_face_index(Block &block) {
    auto *offsets = block.take<uint32_t>(FACE_CELLS + 1);
    offsets[0] = 0;
    Tail<uint16_t> ids(block);
    for_each_cell(FACE_STEPS,
                  [&](int cell, const math::Vector &center, float radius) {
                    for (size_t id = 0; id < faces.size(); ++id) {
                      const Face &face = faces[id];
                      if (face.radius >= 0.0f && math::dot(face.axis, center) >=
                                                     cosf(face.radius + radius))
                        ids.push(static_cast<uint16_t>(id));
                    }
                    offsets[cell + 1] = static_cast<uint32_t>(ids.size());
                  });
    face_ids = ids.commit();
    face_offsets = {offsets, FACE_CELLS + 1};
  }

  // Index per cell: every edge whose distance from the cell center can be
  // within the cell's diameter of the nearest edge. A first pass keeps each
  // cell's nearest chord in its offset slot; the second lists the edges
  // within the cell's limit.
  HS_FLASH_MEMBER __attribute__((flatten)) void
  build_edge_index(float nearest_reach, Block &block, Arena &scratch,
                   Arena &rows) {
    using namespace tiling_edges_detail;
    ScratchScope transient(scratch);
    steps = edge_steps(edges.size());
    const int cells = 6 * steps * steps;
    auto *offsets = block.take<uint32_t>(cells + 1);
    offsets[0] = 0;
    auto *row_reach = scratch.make_n<float>(6 * steps);
    GridRow grid(steps, scratch);
    // Nearest chords up to this were found among the reached edges.
    const float resolved = chord_squared(nearest_reach);
    float widest = 0.0f;
    sweep_cells(
        steps, nearest_reach, [&](int, int) { return nearest_reach; }, scratch,
        rows,
        [&](int face, int x, int y, std::span<const uint16_t> near) {
          grid.load(face, y);
          float nearest = nearest_chord(grid.centers[x], near);
          if (near.empty() || nearest > resolved)
            nearest = nearest_chord_full(grid.centers[x]);
          offsets[(face * steps + y) * steps + x + 1] =
              std::bit_cast<uint32_t>(nearest);
          const float angle = cell_limit(nearest, grid.radius_chords[x]).angle;
          float &reach = row_reach[face * steps + y];
          reach = std::max(reach, angle);
          widest = std::max(widest, angle);
        });
    Tail<uint16_t> ids(block);
    sweep_cells(
        steps, widest,
        [&](int face, int y) { return row_reach[face * steps + y]; }, scratch,
        rows,
        [&](int face, int x, int y, std::span<const uint16_t> near) {
          grid.load(face, y);
          const auto &center = grid.centers[x];
          uint32_t &slot = offsets[(face * steps + y) * steps + x + 1];
          const float nearest = std::bit_cast<float>(slot);
          const float radius_chord = grid.radius_chords[x];
          const CellLimit bounds = cell_limit(nearest, radius_chord);
          float limit = bounds.exact ? bounds.upper : -1.0f;
          for (uint16_t id : near) {
            const auto &edge = edges[id];
            const float sine = math::dot(edge.normal, center);
            if (beyond(sine * sine, bounds.upper))
              continue;
            const float chord = edge_chord_squared(edge, center, vertices);
            if (chord > bounds.upper)
              continue;
            if (chord > bounds.lower) {
              if (limit < 0.0f)
                limit = coarse_limit(nearest, radius_chord);
              if (chord > limit)
                continue;
            }
            ids.push(id);
          }
          slot = static_cast<uint32_t>(ids.size());
        });
    edge_ids = ids.commit();
    edge_offsets = {offsets, static_cast<size_t>(cells) + 1};
  }

  // A coarse cell's admission limit, as the full scan computes it.
  static HS_HOT_INLINE float coarse_limit(float nearest, float radius_chord) {
    using namespace tiling_edges_detail;
    const float radius = chord_angle(radius_chord) + CAP_ANGLE_PAD;
    return chord_squared(
        std::min(math::PI_F, chord_angle(nearest) + 2.0f * radius));
  }

  // Bounds on coarse_limit and an angle at least its own.
  struct CellLimit {
    float lower;
    float upper;
    float angle;
    bool exact;
  };

  static HS_HOT_INLINE CellLimit cell_limit(float nearest, float radius_chord) {
    using namespace tiling_edges_detail;
    float nearest_angle;
    float radius_angle;
    float estimate;
    if (series_angle(nearest, nearest_angle) &&
        series_angle(radius_chord, radius_angle)) {
      const float angle = nearest_angle + 2.0f * (radius_angle + CAP_ANGLE_PAD);
      if (series_chord(angle, estimate))
        return {estimate * (1.0f - SERIES_TOLERANCE),
                estimate * (1.0f + SERIES_TOLERANCE),
                angle * (1.0f + SERIES_TOLERANCE), false};
    }
    const float limit = coarse_limit(nearest, radius_chord);
    const float radius = chord_angle(radius_chord) + CAP_ANGLE_PAD;
    return {limit, limit,
            std::min(math::PI_F, chord_angle(nearest) + 2.0f * radius), true};
  }

  // Series forms of chord_angle and chord_squared over the domains they
  // accept; truncation and rounding stay far inside SERIES_TOLERANCE.
  static HS_HOT_INLINE bool series_angle(float chord, float &angle) {
    if (chord > SERIES_CHORD)
      return false;
    const float z = 0.5f * sqrtf(chord);
    const float z2 = z * z;
    angle = 2.0f * z *
            (1.0f + z2 * (1.0f / 6.0f +
                          z2 * (3.0f / 40.0f + z2 * (5.0f / 112.0f +
                                                     z2 * (35.0f / 1152.0f)))));
    return true;
  }

  static HS_HOT_INLINE bool series_chord(float angle, float &chord) {
    if (angle > SERIES_ANGLE)
      return false;
    const float a2 = angle * angle;
    chord = a2 * (1.0f - a2 * (1.0f / 12.0f -
                               a2 * (1.0f / 360.0f - a2 * (1.0f / 20160.0f -
                                                           a2 / 1814400.0f))));
    return true;
  }

  // Smallest squared chord from a direction to the listed edges.
  HS_HOT_INLINE float nearest_chord(const math::Vector &direction,
                                    std::span<const uint16_t> near) const {
    using namespace tiling_edges_detail;
    size_t first = 0;
    float least = 2.0f;
    for (size_t i = 0; i < near.size(); ++i) {
      const float sine = math::dot(edges[near[i]].normal, direction);
      if (sine * sine < least) {
        least = sine * sine;
        first = i;
      }
    }
    if (near.empty())
      return 4.0f;
    float nearest = edge_chord_squared(edges[near[first]], direction, vertices);
    for (size_t i = 0; i < near.size(); ++i) {
      const auto &edge = edges[near[i]];
      const float sine = math::dot(edge.normal, direction);
      if (i != first && !beyond(sine * sine, nearest))
        nearest =
            std::min(nearest, edge_chord_squared(edge, direction, vertices));
    }
    return nearest;
  }

  // Smallest squared chord from a direction to every edge.
  HS_FLASH_MEMBER float
  nearest_chord_full(const math::Vector &direction) const {
    using namespace tiling_edges_detail;
    float nearest = 4.0f;
    for (const auto &edge : edges)
      nearest =
          std::min(nearest, edge_chord_squared(edge, direction, vertices));
    return nearest;
  }

  // Whether a squared plane sine proves an edge's squared chord exceeds limit.
  static HS_HOT_INLINE bool beyond(float sine_squared, float limit) {
    return limit < 2.0f &&
           sine_squared >
               limit * (1.0f - 0.25f * limit) * PRUNE_SCALE + PRUNE_FLOOR;
  }

  // Largest distance from a face's vertex centroid to its nearest edge, which
  // estimates the farthest any direction lies from every edge.
  HS_FLASH_MEMBER float
  face_center_reach(const MeshCopy &mesh,
                    std::span<const EdgeKey> edge_keys) const {
    using namespace tiling_edges_detail;
    float farthest = 0.0f;
    size_t offset = 0;
    for (uint8_t sides : mesh.counts) {
      const auto ids = mesh.indices.subspan(offset, sides);
      offset += sides;
      math::Vector sum{};
      for (uint16_t id : ids)
        sum += vertices[id];
      const auto center = math::normalized_or(sum, vertices[ids.front()]);
      float nearest = 4.0f;
      uint16_t previous = ids.back();
      for (uint16_t id : ids) {
        const auto match = std::lower_bound(edge_keys.begin(), edge_keys.end(),
                                            edge_key(previous, id));
        nearest = std::min(nearest,
                           edge_chord_squared(edges[match - edge_keys.begin()],
                                              center, vertices));
        previous = id;
      }
      farthest = std::max(farthest, chord_angle(nearest));
    }
    return farthest;
  }

  // Index at SphereTilingNarrowIndex::STEPS: every edge within OUTER of a
  // cell, flagged when its arc contains the whole cell's projected feet.
  // Left empty when the candidate offsets outgrow their 15-bit field.
  HS_FLASH_MEMBER __attribute__((flatten)) void
  build_narrow_index(Block &block, Arena &scratch, Arena &rows) {
    using namespace tiling_edges_detail;
    constexpr int STEPS = SphereTilingNarrowIndex::STEPS;
    constexpr float REACH = 2.0f / STEPS * std::numbers::sqrt2_v<float> * 0.5f +
                            SphereTilingNarrowIndex::OUTER + NARROW_PAD;
    ScratchScope transient(scratch);
    const float limit_sine = sinf(SphereTilingNarrowIndex::OUTER + NARROW_PAD);
    uint8_t *const mark = block.position();
    auto *cells = block.take<uint16_t>(6 * STEPS * STEPS);
    GridRow grid(STEPS, scratch);
    Tail<uint16_t> candidates(block);
    bool overflow = false;
    sweep_cells(
        STEPS, REACH, [](int, int) { return REACH; }, scratch, rows,
        [&](int face, int x, int y, std::span<const uint16_t> near) {
          uint16_t &slot = cells[(face * STEPS + y) * STEPS + x];
          slot = 0;
          if (near.empty() || overflow)
            return;
          grid.load(face, y);
          const auto &center = grid.centers[x];
          const math::Vector corners[] = {grid.top_left[x], grid.top_right[x],
                                          grid.bottom_left[x],
                                          grid.bottom_right[x]};
          const float radius_chord = grid.radius_chords[x];
          float angle;
          float estimate;
          float limit_chord = -1.0f;
          float upper;
          float lower;
          if (series_angle(radius_chord, angle) &&
              series_chord(angle + SphereTilingNarrowIndex::OUTER + NARROW_PAD,
                           estimate)) {
            upper = estimate * (1.0f + SERIES_TOLERANCE);
            lower = estimate * (1.0f - SERIES_TOLERANCE);
          } else {
            limit_chord =
                chord_squared(chord_angle(radius_chord) +
                              SphereTilingNarrowIndex::OUTER + NARROW_PAD);
            upper = lower = limit_chord;
          }
          const size_t begin = candidates.size();
          for (uint16_t id : near) {
            const auto &edge = edges[id];
            const float sine = math::dot(edge.normal, center);
            if (beyond(sine * sine, upper))
              continue;
            const float chord = edge_chord_squared(edge, center, vertices);
            if (chord > upper)
              continue;
            if (chord > lower) {
              if (limit_chord < 0.0f)
                limit_chord =
                    chord_squared(chord_angle(radius_chord) +
                                  SphereTilingNarrowIndex::OUTER + NARROW_PAD);
              if (chord > limit_chord)
                continue;
            }
            bool front = true, back = true, interior = edge.inv_length > 0.0f;
            for (const auto &corner : corners) {
              const float corner_sine = math::dot(edge.normal, corner);
              front &= corner_sine > limit_sine;
              back &= corner_sine < -limit_sine;
              interior &= math::dot(edge.from_a, corner) > NARROW_PAD &&
                          math::dot(edge.from_b, corner) > NARROW_PAD;
            }
            if (front || back)
              continue;
            candidates.push(static_cast<uint16_t>(
                id | (interior ? SphereTilingNarrowIndex::INTERIOR : 0)));
          }
          const size_t count = candidates.size() - begin;
          if (count == 1) {
            slot = candidates.back() | SphereTilingNarrowIndex::LAST;
            candidates.pop();
          } else if (count > 1) {
            overflow = begin + 1 >= SphereTilingNarrowIndex::LAST;
            slot = static_cast<uint16_t>(begin + 1);
            candidates.back() |= SphereTilingNarrowIndex::LAST;
          }
        });
    if (overflow) {
      block.rewind(mark);
      return;
    }
    narrow_cells = {cells, static_cast<size_t>(6 * STEPS * STEPS)};
    narrow_candidates = candidates.commit();
  }

  // Centers, corners and squared corner radii of one row of a cube-map grid,
  // each corner computed from its own cell's coordinates.
  struct GridRow {
    int steps;
    float step;
    int face = -1;
    int row = -1;
    math::Vector *centers;
    math::Vector *top_left;
    math::Vector *top_right;
    math::Vector *bottom_left;
    math::Vector *bottom_right;
    float *radius_chords;

    HS_HOT_INLINE GridRow(int cells, Arena &scratch)
        : steps(cells), step(2.0f / cells),
          centers(scratch.allocate_n<math::Vector>(5 * cells)),
          top_left(centers + cells), top_right(top_left + cells),
          bottom_left(top_right + cells), bottom_right(bottom_left + cells),
          radius_chords(scratch.allocate_n<float>(cells)) {}

    // Corners shared by neighbors are reused only when bit-identical.
    HS_HOT_INLINE void load(int next_face, int next_row) {
      if (next_face == face && next_row == row)
        return;
      const float v = -1.0f + step * next_row;
      const float below = v + step;
      const bool shared = next_face == face && next_row == row + 1 &&
                          -1.0f + step * row + step == v;
      face = next_face;
      row = next_row;
      for (int x = 0; x < steps; ++x) {
        const float u = -1.0f + step * x;
        top_left[x] = shared ? bottom_left[x]
                             : Pattern::sphere_face_direction(face, u, v);
        bottom_left[x] = Pattern::sphere_face_direction(face, u, below);
      }
      for (int x = 0; x < steps; ++x) {
        const float u = -1.0f + step * x;
        const float right = u + step;
        const bool next = x + 1 < steps && -1.0f + step * (x + 1) == right;
        top_right[x] = next ? top_left[x + 1]
                            : Pattern::sphere_face_direction(face, right, v);
        bottom_right[x] =
            next ? bottom_left[x + 1]
                 : Pattern::sphere_face_direction(face, right, below);
        centers[x] = Pattern::sphere_face_direction(face, u + 0.5f * step,
                                                    v + 0.5f * step);
        radius_chords[x] = std::max(
            std::max(math::distance_squared(centers[x], top_left[x]),
                     math::distance_squared(centers[x], top_right[x])),
            std::max(math::distance_squared(centers[x], bottom_left[x]),
                     math::distance_squared(centers[x], bottom_right[x])));
      }
    }
  };

  // Cells of one face whose centers an edge's arc can come within reach of.
  struct CellRange {
    uint16_t edge;
    uint8_t x0;
    uint8_t x1;
    uint8_t y0;
    uint8_t y1;
    bool long_arc;
  };

  // Components of a vector along a cube face's normal and its u and v chart
  // axes, matching Pattern::sphere_face_direction.
  struct ChartTerms {
    float w;
    float u;
    float v;
  };

  static HS_HOT_INLINE ChartTerms chart_terms(const math::Vector &a, int face) {
    switch (face) {
    case 0:
      return {a.x, a.z, a.y};
    case 1:
      return {-a.x, -a.z, a.y};
    case 2:
      return {a.y, a.x, a.z};
    case 3:
      return {-a.y, a.x, -a.z};
    case 4:
      return {a.z, a.x, a.y};
    default:
      return {-a.z, -a.x, a.y};
    }
  }

  // Cells along one axis whose centers lie in [lo, hi].
  static HS_HOT_INLINE bool center_span(float lo, float hi, int steps,
                                        uint8_t &first, uint8_t &last) {
    const float scale = 0.5f * steps;
    const int from = static_cast<int>(
        ceilf(hs::clamp(lo, -1.0f, 1.0f) * scale + scale - 0.5f - SPAN_SLACK));
    const int to = static_cast<int>(
        floorf(hs::clamp(hi, -1.0f, 1.0f) * scale + scale - 0.5f + SPAN_SLACK));
    if (to < 0 || from >= steps || from > to)
      return false;
    first = static_cast<uint8_t>(std::max(from, 0));
    last = static_cast<uint8_t>(std::min(to, steps - 1));
    return true;
  }

  // Cells along one chart axis covering a cap whose unit center has
  // components c along the axis and w along the face normal.
  static HS_HOT_INLINE bool axis_span(float c, float w, float sine, int steps,
                                      uint8_t &first, uint8_t &last) {
    float lo = -1.0f;
    float hi = 1.0f;
    const float off_axis = c * c + w * w;
    if (w > 0.0f && sine * sine < off_axis) {
      const float spread = sine / sqrtf(off_axis - sine * sine);
      const float slope = c / w;
      if (1.0f + slope * spread > 0.0f)
        lo = std::max(lo, (slope - spread) / (1.0f + slope * spread));
      if (1.0f - slope * spread > 0.0f)
        hi = std::min(hi, (slope + spread) / (1.0f - slope * spread));
    }
    return center_span(lo, hi, steps, first, last);
  }

  // Cells of a face within reach of an edge's arc, bounded through the cap of
  // radius half-length plus reach about the arc's midpoint.
  HS_HOT_INLINE bool cell_range(const SphereTilingEdge &edge, float reach,
                                int face, int steps, CellRange &range) const {
    const auto &a = vertices[edge.a];
    const auto &b = vertices[edge.b];
    const bool arc = edge.inv_length > 0.0f;
    const auto mid = arc ? math::normalized_or(a + b, a) : a;
    const float radius = (arc ? 0.5f / edge.inv_length : 0.0f) + reach;
    const auto [w, cu, cv] = chart_terms(mid, face);
    // Face directions satisfy w >= |c|; a cap moves each by at most radius.
    if (w - fabsf(cu) < -2.0f * radius || w - fabsf(cv) < -2.0f * radius)
      return false;
    const float sine = radius < 0.5f * math::PI_F ? sinf(radius) : 1.0f;
    range.long_arc = arc && 0.5f / edge.inv_length > 0.5f * reach;
    return axis_span(cu, w, sine, steps, range.x0, range.x1) &&
           axis_span(cv, w, sine, steps, range.y0, range.y1);
  }

  // Row-line interval [lo, hi] of u where a plane's signed sine at
  // P = n + u a_u + v a_v satisfies |T.P| <= sine |P|, with T.P = g + d u.
  // Returns false when the band is unbounded along the row.
  static HS_HOT_INLINE bool plane_band(float g, float d, float v, float sine,
                                       float &lo, float &hi, bool &empty) {
    const float s2 = sine * sine;
    const float a = d * d - s2;
    if (a <= 0.0f)
      return false;
    const float b = g * d;
    const float c = g * g - s2 * (1.0f + v * v);
    const float discriminant = b * b - a * c;
    empty = discriminant < 0.0f;
    if (empty)
      return true;
    const float root = sqrtf(discriminant);
    const float inverse = 1.0f / a;
    lo = (-b - root) * inverse;
    hi = (-b + root) * inverse;
    return true;
  }

  // Narrows a row's cap span to centers within reach of the arc: inside the
  // band about its great circle and not beyond either endpoint's plane.
  static HS_HOT_INLINE bool arc_span(const SphereTilingEdge &edge, int face,
                                     bool long_arc, float v, float sine,
                                     int steps, uint8_t &first, uint8_t &last) {
    float lo = -1.0f;
    float hi = 1.0f;
    float from;
    float to;
    bool empty = false;
    const auto row_terms = [&](const math::Vector &axis, float &g, float &d) {
      const ChartTerms terms = chart_terms(axis, face);
      g = terms.w + v * terms.v;
      d = terms.u;
    };
    float g;
    float d;
    row_terms(edge.normal, g, d);
    if (plane_band(g, d, v, sine, from, to, empty)) {
      if (empty)
        return false;
      lo = std::max(lo, from);
      hi = std::min(hi, to);
    }
    // Past an endpoint the arc is no nearer than that endpoint's plane.
    for (const auto &tangent : {edge.from_a, edge.from_b}) {
      if (!long_arc)
        break;
      row_terms(tangent, g, d);
      if (plane_band(g, d, v, sine, from, to, empty) && !empty) {
        if (d > 0.0f)
          lo = std::max(lo, from);
        else
          hi = std::min(hi, to);
      }
    }
    uint8_t span_first;
    uint8_t span_last;
    if (lo > hi || !center_span(lo, hi, steps, span_first, span_last))
      return false;
    first = std::max(first, span_first);
    last = std::min(last, span_last);
    return first <= last;
  }

  // Visits every cell of a grid in face, row, column order with the ascending
  // IDs of the edges that may come within reach of its center, a superset of
  // those that do. `cap_reach` bounds every row's `row_reach(face, row)`.
  template <typename RowReach, typename Visit>
  HS_FLASH_MEMBER __attribute__((flatten)) void
  sweep_cells(int steps, float cap_reach, RowReach &&row_reach, Arena &lists,
              Arena &work, Visit &&visit) const {
    ScratchScope transient(lists);
    auto *ranges = lists.allocate_n<CellRange>(edges.size());
    std::array<uint32_t, SphereTilingNarrowIndex::STEPS + 1> begin;
    std::array<uint32_t, SphereTilingNarrowIndex::STEPS + 1> starts;
    std::array<uint32_t, SphereTilingNarrowIndex::STEPS> next;
    const float step = 2.0f / steps;
    for (int face = 0; face < 6; ++face) {
      size_t count = 0;
      for (size_t id = 0; id < edges.size(); ++id)
        if (cell_range(edges[id], cap_reach + SWEEP_MARGIN, face, steps,
                       ranges[count]))
          ranges[count++].edge = static_cast<uint16_t>(id);
      ScratchScope face_scope(work);
      const size_t slots = std::max<size_t>(count, 1);
      auto *order = work.allocate_n<uint16_t>(slots);
      auto *live = work.allocate_n<uint16_t>(slots);
      auto *merged = work.allocate_n<uint16_t>(slots);
      auto *spans = work.allocate_n<uint8_t>(2 * slots);
      const size_t room =
          (work.get_capacity() - work.get_offset()) / sizeof(uint16_t) -
          alignof(uint16_t);
      auto *row = work.allocate_n<uint16_t>(room);
      // Ranges by first row, each row's in ascending edge order.
      std::fill_n(starts.begin(), steps + 1, 0u);
      for (size_t i = 0; i < count; ++i)
        ++starts[ranges[i].y0 + 1];
      for (int y = 0; y < steps; ++y)
        starts[y + 1] += starts[y];
      std::copy_n(starts.begin(), steps, next.begin());
      for (size_t i = 0; i < count; ++i)
        order[next[ranges[i].y0]++] = static_cast<uint16_t>(i);
      size_t live_count = 0;
      for (int y = 0; y < steps; ++y) {
        size_t merged_count = 0;
        size_t fresh = starts[y];
        for (size_t i = 0; i < live_count; ++i) {
          const uint16_t kept = live[i];
          for (; fresh < starts[y + 1] && order[fresh] < kept; ++fresh)
            merged[merged_count++] = order[fresh];
          if (ranges[kept].y1 >= y)
            merged[merged_count++] = kept;
        }
        for (; fresh < starts[y + 1]; ++fresh)
          merged[merged_count++] = order[fresh];
        std::swap(live, merged);
        live_count = merged_count;
        const float v = -1.0f + step * (y + 0.5f);
        const float reach = row_reach(face, y) + SWEEP_MARGIN;
        const float sine = reach < 0.5f * math::PI_F ? sinf(reach) : 1.0f;
        std::fill_n(begin.begin(), steps + 1, 0u);
        for (size_t i = 0; i < live_count; ++i) {
          const CellRange &range = ranges[live[i]];
          uint8_t &first = spans[2 * i];
          uint8_t &last = spans[2 * i + 1];
          first = range.x0;
          last = range.x1;
          if (!arc_span(edges[range.edge], face, range.long_arc, v, sine, steps,
                        first, last)) {
            first = 1;
            last = 0;
            continue;
          }
          for (int x = first; x <= last; ++x)
            ++begin[x + 1];
        }
        for (int x = 0; x < steps; ++x)
          begin[x + 1] += begin[x];
        // no death case: the budget test builds every admitted solid.
        HS_CHECK(begin[steps] <= room,
                 "sphere tiling: sweep row exceeds scratch");
        std::copy_n(begin.begin(), steps, next.begin());
        for (size_t i = 0; i < live_count; ++i)
          for (int x = spans[2 * i]; x <= spans[2 * i + 1]; ++x)
            row[next[x]++] = ranges[live[i]].edge;
        for (int x = 0; x < steps; ++x)
          visit(face, x, y,
                std::span<const uint16_t>(row + begin[x],
                                          begin[x + 1] - begin[x]));
      }
    }
  }

  bool contains(const Face &face, const math::Vector &direction) const {
    if (face.radius < 0.0f || math::dot(face.axis, direction) <= CHART_EPSILON)
      return false;

    float px;
    float py;
    SDF::gnomonic_project(direction, face.axis, face.u, face.v, px, py);
    bool inside = false;
    Point a = corners[face.first + face.count - 1];
    for (size_t index = face.first; index < face.first + face.count; ++index) {
      const Point b = corners[index];
      const float dx = b.x - a.x;
      const float dy = b.y - a.y;
      const float cross = dx * (py - a.y) - dy * (px - a.x);
      if (fabsf(cross) <= CHART_EPSILON * (fabsf(dx) + fabsf(dy)) &&
          px >= std::min(a.x, b.x) - CHART_EPSILON &&
          px <= std::max(a.x, b.x) + CHART_EPSILON &&
          py >= std::min(a.y, b.y) - CHART_EPSILON &&
          py <= std::max(a.y, b.y) + CHART_EPSILON)
        return true;
      if ((a.y > py) != (b.y > py) && (cross > 0.0f) == (dy > 0.0f))
        inside = !inside;
      a = b;
    }
    return inside;
  }

  std::pair<float, uint16_t> boundary_at(const Face &face,
                                         const math::Vector &direction) const {
    using namespace tiling_edges_detail;
    uint16_t nearest = corners[face.first].edge;
    float best = edge_chord_squared(edges[nearest], direction, vertices);
    for (size_t i = face.first + 1; i < face.first + face.count; ++i) {
      const uint16_t edge = corners[i].edge;
      const float distance =
          edge_chord_squared(edges[edge], direction, vertices);
      if (distance < best) {
        best = distance;
        nearest = edge;
      }
    }
    return {chord_angle(best), nearest};
  }

  template <class Visitor>
  static HS_HOT_INLINE void for_each_cell(int steps, Visitor visit) {
    const float step = 2.0f / steps;
    for (int cell = 0; cell < 6 * steps * steps; ++cell) {
      const int side = cell / (steps * steps);
      const float u = -1.0f + step * (cell % steps);
      const float v = -1.0f + step * ((cell / steps) % steps);
      const auto center = Pattern::sphere_face_direction(side, u + 0.5f * step,
                                                         v + 0.5f * step);
      float largest_chord = 0.0f;
      for (float x : {u, u + step})
        for (float y : {v, v + step})
          largest_chord =
              std::max(largest_chord,
                       math::distance_squared(
                           center, Pattern::sphere_face_direction(side, x, y)));
      const float radius =
          tiling_edges_detail::chord_angle(largest_chord) + CAP_ANGLE_PAD;
      visit(cell, center, radius);
    }
  }
};

} // namespace Pullback::Source
