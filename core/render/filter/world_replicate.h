/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/filter/pipeline.h"
#include "core/math/geometry.h"
#include "core/color/color.h"

/**
 * @file world_replicate.h
 * @brief Filter::World::Replicate: emits evenly-spaced copies of the geometry
 * rotated around the Y axis.
 */

namespace Filter {
namespace World {

/**
 * @brief Selected-copy timed interface for spatial replication stages.
 * @tparam Derived Stage supplying copy_count() and selected_point().
 */
template <typename Derived> class SelectedCopy : public Is3D {
public:
  /** @brief Selects the spatial copy used by timed operations.
   * @param index Copy index, clamped to the available copies. */
  void select_copy(int index) {
    selected_copy = hs::clamp(index, 0, derived().copy_count() - 1);
  }

  /** @brief Evaluates the selected spatial copy.
   * @param cv Canvas identifying the frame.
   * @param v Source point.
   * @param t Frame fraction, unchanged by spatial replication.
   * @return Point transformed for the selected copy. */
  HS_HOT_INLINE math::Vector at_time(const Canvas &cv, const math::Vector &v,
                                     float t) const {
    (void)cv;
    (void)t;
    return derived().selected_point(v);
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
      return forward(derived().selected_point(a), derived().selected_point(b),
                     nullptr);
    const math::Basis transformed{derived().selected_point(pb->u),
                                  derived().selected_point(pb->v),
                                  derived().selected_point(pb->w)};
    return forward(derived().selected_point(a), derived().selected_point(b),
                   &transformed);
  }

protected:
  int selected_copy = 0; /**< Copy fixed for the current timed walk. */

private:
  /** @brief Accesses the spatial replication stage.
   * @return Concrete stage. */
  const Derived &derived() const { return static_cast<const Derived &>(*this); }
};

/**
 * @brief Replicates geometry by rotating it around the Y-axis.
 * @details Every copy shares the source `age` (replication is spatial, not
 * temporal).
 */
template <int W> class Replicate : public SelectedCopy<Replicate<W>> {
public:
  static constexpr bool rigid_world_transform = true; ///< Preserves geodesics.
  /**
   * @brief Builds a replicator emitting @p count evenly-spaced Y-axis copies.
   * @param count Desired copy count; clamped to [1, W].
   */
  Replicate(int count) { set_count(count); }

  /**
   * @brief Changes the number of evenly spaced copies.
   * @param new_count Desired copy count; clamped to [1, W].
   */
  void set_count(int new_count) {
    count = hs::clamp(new_count, 1, W);
    this->selected_copy = std::min(this->selected_copy, count - 1);
    step = math::make_rotation(math::Y_AXIS, 2 * math::PI_F / count);
  }
  /** @brief Number of independently walked spatial copies.
   * @return Available spatial copy count. */
  int copy_count() const { return count; }

  /**
   * @brief Emits the point plus count-1 rotated copies around the Y axis.
   * @param v World-space point to replicate; must be near unit length.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const math::Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    math::Vector r = v;
    pass(r, color, age, alpha);
    for (int i = 1; i < count; i++) {
      r = math::rotate(r, step);
      // First-order renormalization; valid only while r stays near unit length.
      r = r * (1.5f - 0.5f * math::dot(r, r));
      pass(r, color, age, alpha);
    }
  }

  /**
   * @brief Re-emits a clip-cull edge under every copy's Y-axis rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any copy of the edge could intersect the clip band.
   * @details Mirrors plot()'s rotation loop so the cull sees the longitudes the
   *          copies are drawn at, not the source geometry's.
   */
  template <typename FwdFn>
  bool cull_edge(const math::Vector &a, const math::Vector &b,
                 const math::Basis *pb, FwdFn &&forward) const {
    if (forward(a, b, pb))
      return true;
    math::Vector ra = a, rb = b;
    math::Basis rp;
    if (pb)
      rp = *pb;
    for (int i = 1; i < count; i++) {
      ra = math::rotate(ra, step).normalized();
      rb = math::rotate(rb, step).normalized();
      if (pb) {
        rp = math::rotate(rp, step);
        renormalize(rp);
        if (forward(ra, rb, &rp))
          return true;
      } else if (forward(ra, rb, nullptr)) {
        return true;
      }
    }
    return false;
  }

private:
  friend class SelectedCopy<Replicate<W>>;

  /** @brief Transforms a point into the selected copy.
   * @param v Source point.
   * @return Selected copy position. */
  HS_HOT_INLINE math::Vector selected_point(const math::Vector &v) const {
    math::Vector r = v;
    for (int i = 0; i < this->selected_copy; ++i) {
      r = math::rotate(r, step);
      r = r * (1.5f - 0.5f * math::dot(r, r));
    }
    return r;
  }

  /**
   * @brief Renormalizes a chained-rotation basis.
   * @param b Basis to restore to unit axes in place.
   * @details cull_edge rotates the basis once per copy off the previous copy's,
   * so up to W - 1 chained rotations drift the axes off unit length.
   */
  static void renormalize(math::Basis &b) {
    b.u = b.u.normalized();
    b.v = b.v.normalized();
    b.w = b.w.normalized();
  }

  int count;             /**< Number of copies emitted, in [1, W]. */
  math::Quaternion step; /**< Per-copy Y-axis rotation (2*pi / count). */
};

} // namespace World
} // namespace Filter
