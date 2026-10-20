/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "core/render/filter/world_orient.h"
#include "core/math/geometry.h"
#include "core/color/color.h"

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
  static constexpr bool rigid_world_transform = true; ///< Preserves geodesics.
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
    HS_CHECK(std::size(vertices) >= static_cast<size_t>(N),
             "VertexReplicate: vertex array is smaller than replica count");
    for (int i = 0; i < N; ++i)
      rotations[i] = math::make_rotation(vertices[0], vertices[i]);
  }

  /** @brief Number of independently walked spatial copies.
   * @return Available spatial copy count. */
  int copy_count() const { return N; }

  /** @brief Selects the spatial copy used by timed operations.
   * @param index Copy index, clamped to the available copies. */
  void select_copy(int index) { selected_copy = hs::clamp(index, 0, N - 1); }

  /** @brief Evaluates the selected spatial copy.
   * @param cv Canvas identifying the frame.
   * @param v Source point.
   * @param t Frame fraction, unchanged by spatial replication.
   * @return Point transformed for the selected copy. */
  HS_HOT_INLINE math::Vector at_time(const Canvas &cv, const math::Vector &v,
                                     float t) const {
    (void)cv;
    (void)t;
    return selected_point(v);
  }

  /** @brief Emits only the selected copy for a timed walk.
   * @param cv Canvas identifying the frame.
   * @param v Source point.
   * @param color Source color.
   * @param age Source age.
   * @param alpha Source opacity.
   * @param t Frame fraction.
   * @param pass Downstream world continuation. */
  template <typename PassFnT>
  HS_HOT_INLINE void plot_at(Canvas &cv, const math::Vector &v,
                             const ::Pixel &color, float age, float alpha,
                             float t, PassFnT &&pass) const {
    pass(at_time(cv, v, t), color, age, alpha);
  }

  /** @brief Culls only the selected spatial copy.
   * @param a,b Source edge endpoints.
   * @param pb Optional source planar basis.
   * @param forward Downstream cull continuation.
   * @return The selected copy's visibility verdict. */
  template <typename FwdFn>
  bool cull_edge_timed(const math::Vector &a, const math::Vector &b,
                       const math::Basis *pb, FwdFn &&forward) const {
    if (!pb)
      return forward(selected_point(a), selected_point(b), nullptr);
    const math::Basis transformed{selected_point(pb->u), selected_point(pb->v),
                                  selected_point(pb->w)};
    return forward(selected_point(a), selected_point(b), &transformed);
  }

  /**
   * @brief Emits one rotated copy of the point per stored vertex rotation.
   * @param v World-space point to replicate.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const math::Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    for (int i = 0; i < N; ++i) {
      pass(math::rotate(v, rotations[i]), color, age, alpha);
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
   * @details Mirrors plot().
   */
  template <typename FwdFn>
  bool cull_edge(const math::Vector &a, const math::Vector &b,
                 const math::Basis *pb, FwdFn &&forward) const {
    for (int i = 0; i < N; ++i) {
      if (forward_rotated_edge(a, b, pb, rotations[i], forward))
        return true;
    }
    return false;
  }

private:
  /** @brief Transforms a point into the selected copy.
   * @param v Source point.
   * @return Selected copy position. */
  HS_HOT_INLINE math::Vector selected_point(const math::Vector &v) const {
    return math::rotate(v, rotations[selected_copy]);
  }

  int selected_copy = 0; /**< Copy fixed for the current timed walk. */

  std::array<math::Quaternion, N>
      rotations; /**< Rotation from vertices[0] to each vertex. */
};

} // namespace World
} // namespace Filter
