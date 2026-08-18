/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "render/filter/world_orient.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_vertex_replicate.h
 * @brief Filter::World::VertexReplicate: emits one copy of the geometry per
 * vertex of a solid, rotated from the first vertex onto each.
 */

namespace Filter {
namespace World {

/**
 * @brief Replicates geometry onto the vertices of a solid.
 * @details Precomputes rotation quaternions from vertex[0] to each other vertex.
 * Every copy carries the source age unchanged (replication is spatial).
 */
template <int N> class VertexReplicate : public Is3D {
public:
  /**
   * @brief Builds from a vertex array, precomputing rotations vertices[0] → each.
   * @tparam VertexArray Indexable container of N unit vectors.
   * @param vertices Vertex positions; rotations map vertices[0] onto each vertex.
   */
  template <typename VertexArray> VertexReplicate(const VertexArray &vertices) {
    set_vertices(vertices);
  }

  /**
   * @brief Rebuilds the rotations from a new vertex array.
   * @tparam VertexArray Indexable container of N unit vectors.
   * @param vertices Vertex positions; rotations map vertices[0] onto each vertex.
   */
  template <typename VertexArray>
  void set_vertices(const VertexArray &vertices) {
    for (int i = 0; i < N; ++i)
      rotations[i] = make_rotation(vertices[0], vertices[i]);
  }

  /**
   * @brief Emits one rotated copy of the point per stored vertex rotation.
   * @param v World-space point to replicate.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    for (int i = 0; i < N; ++i) {
      pass(rotate(v, rotations[i]), color, age, alpha);
    }
  }

  /**
   * @brief Re-emits a clip-cull edge under every stored vertex rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any vertex copy of the edge could intersect the clip band.
   * @details Mirrors plot(). The rotations move latitude, so culling by the
   *          un-rotated endpoints would drop copies the fan-out places inside a
   *          segment band (docs/segmented_stateful_effects_spec.md).
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    for (int i = 0; i < N; ++i) {
      if (forward_rotated_edge(a, b, pb, rotations[i], forward))
        return true;
    }
    return false;
  }

private:
  std::array<Quaternion, N>
      rotations; /**< Rotation from vertices[0] to each vertex. */
};

} // namespace World
} // namespace Filter
