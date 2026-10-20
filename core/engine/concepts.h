/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#pragma once

/**
 * @file concepts.h
 * @brief Callable wrappers, render callback aliases, and pipeline concepts.
 */

#include <concepts>
#include <cstddef>     // std::nullptr_t
#include <memory>      // std::addressof
#include <type_traits> // std::remove_cvref_t
#include <utility>     // std::forward
#include "core/math/3dmath.h"
#include "core/color/pixel.h"       // Pixel
#include "core/platform/platform.h" // Fn
#include "core/render/temporal.h"

namespace math {
struct Basis;
}
class Canvas;

// FunctionRef<Sig>: non-owning borrow for call-scoped callbacks; must not
// outlive the callable. Fn<Sig, Cap>: owning, heap-free inline storage for
// callbacks kept past the creating scope.

struct Fragment;
struct FragmentRegisters;
template <typename Signature> class FunctionRef;

namespace hs {
/** @brief Diverges when an empty FunctionRef is invoked. */
[[noreturn]] void function_ref_empty_call();
} // namespace hs

/**
 * @brief Non-owning, type-erased reference to any callable matching Ret(Args...).
 * @tparam Ret Return type of the wrapped callable.
 * @tparam Args Argument types of the wrapped callable.
 * @details Borrows the callable via a void* context plus a thunk pointer; zero
 * heap allocation. The referenced callable must outlive the FunctionRef.
 */
template <typename Ret, typename... Args> class FunctionRef<Ret(Args...)> {
  // Empty state's thunk: an empty ref diverges instead of calling through null.
  [[noreturn]] static Ret empty_thunk(void *, Args &&...) {
    ::hs::function_ref_empty_call();
  }

  void *ctx = nullptr;
  Ret (*thunk)(void *, Args &&...) = &empty_thunk;

public:
  /**
   * @brief Constructs an empty FunctionRef that refers to no callable.
   */
  FunctionRef() = default;

  /**
   * @brief Constructs an empty FunctionRef from a null pointer literal.
   */
  FunctionRef(std::nullptr_t) {}

  // Explicit copy/move — prevent the generic Callable template from matching.

  /**
   * @brief Copy-constructs a reference to the same callable.
   * @param other Source FunctionRef to copy.
   */
  FunctionRef(const FunctionRef &other) noexcept = default;

  /**
   * @brief Move-constructs a reference to the same callable.
   * @param other Source FunctionRef to move from.
   */
  FunctionRef(FunctionRef &&other) noexcept = default;

  /**
   * @brief Copy-assigns to refer to the same callable as other.
   * @param other Source FunctionRef to copy.
   * @return Reference to this FunctionRef.
   */
  FunctionRef &operator=(const FunctionRef &other) noexcept = default;

  /**
   * @brief Move-assigns to refer to the same callable as other.
   * @param other Source FunctionRef to move from.
   * @return Reference to this FunctionRef.
   */
  FunctionRef &operator=(FunctionRef &&other) noexcept = default;

  /**
   * @brief Wraps a plain function pointer.
   * @param func Function pointer with signature Ret(Args...); stored in ctx.
   * @details The function-pointer <-> void* round-trip is
   * conditionally-supported ([expr.reinterpret.cast]). A null `func` produces
   * an empty ref.
   */
  FunctionRef(Ret (*func)(Args...)) noexcept
      : ctx(reinterpret_cast<void *>(func)) {
    static_assert(
        sizeof(func) == sizeof(void *),
        "FunctionRef requires function and object pointers to share a "
        "width (true on all supported targets)");
    if (func == nullptr)
      return;
    thunk = [](void *ptr, Args &&...args) -> Ret {
      return (reinterpret_cast<Ret (*)(Args...)>(ptr))(
          std::forward<Args>(args)...);
    };
  }

  /**
   * @brief Stores a non-throwing function pointer by value.
   * @param func Function to call; null yields an empty reference.
   */
  FunctionRef(Ret (*func)(Args...) noexcept) noexcept
      : FunctionRef(static_cast<Ret (*)(Args...)>(func)) {}

  /**
   * @brief Wraps a non-const lvalue callable (functor or lambda).
   * @tparam Callable Type of the callable; must be invocable as Ret(Args...) and
   * not itself a FunctionRef.
   * @param callable Callable whose address is stored; must outlive this ref.
   */
  template <typename Callable>
    requires std::is_invocable_r_v<Ret, Callable &, Args...> &&
             (!std::is_base_of_v<FunctionRef, std::decay_t<Callable>>)
  FunctionRef(Callable &callable) noexcept : ctx(std::addressof(callable)) {
    thunk = [](void *ptr, Args &&...args) -> Ret {
      if constexpr (std::is_void_v<Ret>) {
        (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
      } else {
        return (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
      }
    };
  }

  /**
   * @brief Wraps a const lvalue callable (functor or lambda).
   * @tparam Callable Type of the callable; must be const-invocable with Args...
   * and not itself a FunctionRef.
   * @param callable Const callable whose address is stored; must outlive this
   * ref. The const is cast away into ctx and restored in the thunk.
   * @details Also binds temporaries for immediate-use borrows. A `mutable`
   * temporary lambda binds to neither overload; drop the `mutable` or pass a
   * named lvalue.
   */
  template <typename Callable>
    requires std::is_invocable_r_v<Ret, const Callable &, Args...> &&
             (!std::is_base_of_v<FunctionRef, std::decay_t<Callable>>)
  FunctionRef(const Callable &callable) noexcept
      : ctx(const_cast<void *>(
            static_cast<const void *>(std::addressof(callable)))) {
    thunk = [](void *ptr, Args &&...args) -> Ret {
      if constexpr (std::is_void_v<Ret>) {
        (*static_cast<const Callable *>(ptr))(std::forward<Args>(args)...);
      } else {
        return (*static_cast<const Callable *>(ptr))(
            std::forward<Args>(args)...);
      }
    };
  }

  /**
   * @brief Invokes the wrapped callable.
   * @param args Arguments forwarded to the callable.
   * @return Result of the wrapped callable.
   */
  inline Ret operator()(Args... args) const {
    return thunk(ctx, std::forward<Args>(args)...);
  }

  /**
   * @brief Tests whether this FunctionRef refers to a callable.
   * @return True if a callable is bound, false if empty.
   */
  [[nodiscard]] explicit operator bool() const { return thunk != &empty_thunk; }
};

/**
 * @brief A FunctionRef meant to be STORED past the call that builds it (e.g. a
 * class member invoked across many frames), not just borrowed for one call.
 * @tparam Signature The callable signature `Ret(Args...)`.
 * @details Identical to FunctionRef except it refuses to bind an rvalue
 * temporary, which would dangle. Adds no data members.
 */
template <typename Signature> class StoredFunctionRef;

/**
 * @brief StoredFunctionRef for the signature `Ret(Args...)`.
 * @tparam Ret Return type.
 * @tparam Args Argument types.
 */
template <typename Ret, typename... Args>
class StoredFunctionRef<Ret(Args...)> : public FunctionRef<Ret(Args...)> {
public:
  using FunctionRef<Ret(Args...)>::FunctionRef;

  StoredFunctionRef() = default;
  /** @brief Constructs an empty reference. */
  StoredFunctionRef(std::nullptr_t) noexcept {}

  // Reject rvalue temporaries the base would accept; the guards keep lvalue
  // callables and copy/move on the inherited ctors.
  template <typename Callable,
            typename = std::enable_if_t<
                !std::is_lvalue_reference_v<Callable> &&
                !std::is_same_v<std::decay_t<Callable>, std::nullptr_t> &&
                !(std::is_pointer_v<std::decay_t<Callable>> &&
                  std::is_convertible_v<Callable, Ret (*)(Args...)>) &&
                !std::is_same_v<std::decay_t<Callable>, StoredFunctionRef>>>
  StoredFunctionRef(Callable &&) = delete;
};

// Borrow-only aliases; stored callbacks use StoredFunctionRef.
/// Screen trail colour: (x, y, t) with t the point's age in [0, 1].
using ScreenTrailFn = FunctionRef<Color4(float, float, float)>;
/// World trail colour: (position, t) with t the point's age in [0, 1].
using WorldTrailFn = FunctionRef<Color4(const math::Vector &, float)>;
/// Fragment shader: shades the fragment at a sphere position.
using FragmentShaderFn = FunctionRef<void(const math::Vector &, Fragment &)>;
/// Time-aware fragment shader: (position, fragment, sweep fraction of the
/// shaded column in [0, 1]).
using TimedShaderFn =
    FunctionRef<void(const math::Vector &, Fragment &, float)>;
/// Vertex shader: rewrites a fragment before projection.
using VertexShaderRef = FunctionRef<void(Fragment &)>;
// Deferred per-control-point shader: receives the (position-shaded) fragment's
// shading registers and its original pre-shader position. It runs after
// projection, so it cannot move the fragment.
/// Deferred shader: (shading registers, pre-shader position).
using DeferredShaderRef =
    FunctionRef<void(FragmentRegisters, const math::Vector &)>;
/// Orientation history visitor: (orientation, t) with t in [0, 1].
using TweenFn = FunctionRef<void(const math::Quaternion &, float)>;
/// Position history visitor: (position, t) with t in [0, 1].
using VectorTweenFn = FunctionRef<void(const math::Vector &, float)>;
// Rasterizer's clip-cull predicate: does the (world-transformed) edge a-b, with
// optional planar basis, intersect the clip band?
/// Clip-cull predicate: (a, b, planar basis or null) intersects the clip band.
using CullEdgePredRef = FunctionRef<bool(
    const math::Vector &, const math::Vector &, const math::Basis *)>;

/**
 * @brief Deterministic ownership mask for dissolve transitions.
 * @details Ownership is hashed from an integer key pair, not a pixel
 * coordinate. Two masks with the same threshold/salt and opposite `invert`
 * partition key pairs exactly. The salt must derive from frame counters/seeds,
 * never wall time (sim/device parity).
 */
struct DissolveMask {
  uint32_t threshold; /**< Owned fraction in [0, 65536]. */
  uint32_t salt;      /**< Per-frame/per-transition hash salt. */
  bool invert; /**< True owns the complement (elements at/above threshold). */

  /**
   * @brief Ownership of the element identified by the key pair (a, b).
   * @param a,b Key components; order-sensitive.
   * @return True when this mask owns the element.
   */
  bool owns(int a, int b) const {
    uint32_t h = static_cast<uint32_t>(a) * 0x9E3779B1u ^
                 static_cast<uint32_t>(b) * 0x85EBCA77u ^ salt;
    h *= 0x27D4EB2Fu;
    h ^= h >> 15;
    return ((h & 0xFFFFu) < threshold) != invert;
  }
};

/**
 * @brief Whether P writes the framebuffer through a cached base.
 * @tparam P Pipeline type; one without the member is not a direct-raster path.
 * @return `P::direct_raster_path`, or false when P lacks it.
 */
template <typename P> consteval bool pipeline_direct_raster_path() {
  if constexpr (requires { P::direct_raster_path; })
    return P::direct_raster_path;
  else
    return false;
}

/**
 * @brief Whether P exposes the three plot() overloads PipelineRef erases.
 * @tparam P Candidate pipeline type.
 * @details Spelled with the same argument types the erasing thunks pass, so it
 * accepts exactly what they compile against.
 */
template <typename P>
concept Plottable =
    requires(P &p, Canvas &cv, const math::Vector &v, const Pixel &c) {
      p.plot(cv, 0.0f, 0.0f, c, 0.0f, 0.0f);
      p.plot(cv, 0, 0, c, 0.0f, 0.0f);
      p.plot(cv, v, c, 0.0f, 0.0f);
    };

/**
 * @brief Dispatches a clip query through world stages or a bare plot provider.
 * @param pipeline Pipeline queried.
 * @param a Edge start (unit sphere point, pre-world-transform).
 * @param b Edge end (unit sphere point, pre-world-transform).
 * @param pb Optional planar basis for the edge (null = geodesic).
 * @param pred Rasterizer's screen-row-span vs clip-band test.
 * @return Whether the edge could intersect the clip band.
 */
template <typename PipelineT, typename Pred>
inline bool pipeline_could_intersect_clip(PipelineT &pipeline,
                                          const math::Vector &a,
                                          const math::Vector &b,
                                          const math::Basis *pb, Pred &&pred) {
  if constexpr (requires {
                  pipeline.could_intersect_clip(a, b, pb,
                                                std::forward<Pred>(pred));
                }) {
    return pipeline.could_intersect_clip(a, b, pb, std::forward<Pred>(pred));
  } else {
    static_assert(
        !requires { PipelineT::any_crosses_segments; },
        "pipeline exposes any_crosses_segments but not "
        "could_intersect_clip (signature drift)");
    return pred(a, b, pb);
  }
}

/**
 * @brief Non-owning, type-erased handle to a rasterizer pipeline.
 * @details Forwards 2D screen-space and 3D world-space plot() calls to the
 * wrapped object; borrows the target and must not outlive it. A direct-raster
 * sink converts only through the Canvas-taking constructor, which checks
 * prepared_for().
 */
class PipelineRef {
  /** @brief Tag selecting the type-erasing constructor. */
  template <bool Timed> struct Erase {};

  void *ctx;
  void (*plot2d)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot2d_int)(void *, Canvas &, int, int, const Pixel &, float, float);
  void (*plot3d)(void *, Canvas &, const math::Vector &, const Pixel &, float,
                 float);
  bool (*cull)(void *, const math::Vector &, const math::Vector &,
               const math::Basis *, CullEdgePredRef);
  /** @brief Shared entries used only by geometry-callback draws. */
  struct TimedEntries {
    /** @brief Evaluates world stages without projection. */
    math::Vector (*world)(void *, Canvas &, const math::Vector &, float);
    /** @brief Projects an explicitly timed point. */
    math::PixelCoords (*probe)(void *, Canvas &, const math::Vector &, float);
    /** @brief Plots an explicitly timed point. */
    float (*plot)(void *, Canvas &, const math::Vector &, const Pixel &, float,
                  float, float, const Temporal::ColumnPrediction *);
    /** @brief Plots a sample after the world stages. */
    float (*transformed_plot)(void *, Canvas &, const math::Vector &,
                              const Pixel &, float, float, float,
                              const Temporal::ColumnPrediction *);
    /** @brief Plots an already-projected timed sample. */
    float (*projected_plot)(void *, Canvas &, math::PixelCoords, const Pixel &,
                            float, float, float,
                            const Temporal::ColumnPrediction *);
    bool rigid;             ///< World stages can run at controls.
    bool dense_screen_safe; ///< Screen stages retain a one-pixel footprint.
    /** @brief Returns the angular world-stage motion bound. */
    float (*bound)(void *);
    /** @brief Returns the spatial copy count. */
    int (*copies)(void *);
    /** @brief Selects one spatial copy. */
    void (*select)(void *, int);
  };
  /** @brief Shared timed entries, absent on ordinary refs. */
  const TimedEntries *timed_entries = nullptr;

  /** @brief Binds ordinary plot and cull entries. */
  template <typename T>
  PipelineRef(T &t, Erase<false>)
      : ctx(std::addressof(t)), world_transform_is_identity([]() HS_HOT_INLINE {
          if constexpr (requires { T::world_transform_is_identity; })
            return T::world_transform_is_identity;
          else
            return true;
        }()),
        has_temporal_motion([]() HS_HOT_INLINE {
          if constexpr (requires { T::has_temporal_motion; })
            return T::has_temporal_motion;
          else
            return true;
        }()) {
    plot2d = [](void *pipeline, Canvas &cv, float x, float y, const Pixel &c,
                float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot2d_int = [](void *pipeline, Canvas &cv, int x, int y, const Pixel &c,
                    float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot3d = [](void *pipeline, Canvas &cv, const math::Vector &v,
                const Pixel &c, float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, v, c, age, alpha);
    };
    cull = still_cull_entry<T>;
  }

  /** @brief Applies an erased pipeline's ordinary edge cull.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param a,b Edge endpoints.
   * @param pb Optional planar basis.
   * @param pred Downstream visibility predicate.
   * @return Whether the edge may intersect the clip. */
  template <typename T>
  HS_O3_FN HS_NOINLINE_NOCLONE static bool
  still_cull_entry(void *pipeline, const math::Vector &a, const math::Vector &b,
                   const math::Basis *pb, CullEdgePredRef pred) {
    return pipeline_could_intersect_clip(*static_cast<T *>(pipeline), a, b, pb,
                                         pred);
  }

  /** @brief Evaluates an erased pipeline's world stages.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param time Frame fraction.
   * @return World-space point. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static math::Vector
  timed_world_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                    float time) {
    if constexpr (requires(T &p) { p.world_at(cv, v, time); })
      return static_cast<T *>(pipeline)->world_at(cv, v, time);
    else
      return v;
  }

  /** @brief Projects an erased pipeline without writing.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param time Frame fraction.
   * @return Projected column and row. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static math::PixelCoords
  timed_probe_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                    float time) {
    return static_cast<T *>(pipeline)->probe(cv, v, time);
  }

  /** @brief Emits an erased timed point through its acceptance gate.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Opacity.
   * @param time Frame fraction.
   * @param prediction Optional column acceptance window.
   * @return Unwrapped landed column. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static float
  timed_plot_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                   const Pixel &c, float age, float alpha, float time,
                   const Temporal::ColumnPrediction *prediction) {
    auto &sink = *static_cast<T *>(pipeline);
    if constexpr (requires {
                    sink.plot_at(cv, v, c, age, alpha, time, prediction);
                  }) {
      return sink.plot_at(cv, v, c, age, alpha, time, prediction);
    } else {
      const auto p = sink.probe(cv, v, time);
      const auto unwrap = [&](auto &frame) HS_HOT_INLINE {
        return prediction ? prediction->unwrap(p.x, frame.width()) : p.x;
      };
      const float landed = unwrap(cv);
      if (!prediction ||
          (prediction->record_sample && prediction->accepts(landed)))
        sink.plot_at(cv, v, c, age, alpha, time);
      return landed;
    }
  }

  /** @brief Dispatches an already-transformed sample.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param v Transformed point.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @param prediction Optional part acceptance and boundary evaluator.
   * @return Actual unwrapped column. */
  template <typename T>
  HS_HOT_INLINE static float
  transformed_plot_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                         const Pixel &c, float age, float alpha, float t,
                         const Temporal::ColumnPrediction *prediction) {
    if constexpr (requires(T &p) {
                    p.plot_transformed_at(cv, v, c, age, alpha, t, prediction);
                  })
      return static_cast<T *>(pipeline)->plot_transformed_at(
          cv, v, c, age, alpha, t, prediction);
    else
      return timed_plot_entry<T>(pipeline, cv, v, c, age, alpha, t, prediction);
  }

  /** @brief Dispatches an already-projected timed sample.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @param prediction Optional landing window and boundary evaluator.
   * @return Actual unwrapped column. */
  template <typename T>
  HS_HOT_INLINE static float
  projected_plot_entry(void *pipeline, Canvas &cv, math::PixelCoords p,
                       const Pixel &c, float age, float alpha, float t,
                       const Temporal::ColumnPrediction *prediction) {
    if constexpr (requires(T &v) {
                    v.plot_projected_at(cv, p, c, age, alpha, t, prediction);
                  })
      return static_cast<T *>(pipeline)->plot_projected_at(cv, p, c, age, alpha,
                                                           t, prediction);
    else {
      static_cast<T *>(pipeline)->plot(cv, p.x, p.y, c, age, alpha);
      return p.x;
    }
  }

  /** @brief Applies an erased pipeline's timed cull chain.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param a,b Source edge endpoints.
   * @param pb Optional planar basis.
   * @param pred Downstream visibility predicate.
   * @return Whether the edge can intersect the clip. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static bool
  timed_cull_entry(void *pipeline, const math::Vector &a, const math::Vector &b,
                   const math::Basis *pb, CullEdgePredRef pred) {
    if constexpr (requires(T &p) {
                    p.could_intersect_timed_clip(a, b, pb, pred);
                  })
      return static_cast<T *>(pipeline)->could_intersect_timed_clip(a, b, pb,
                                                                    pred);
    return pipeline_could_intersect_clip(*static_cast<T *>(pipeline), a, b, pb,
                                         pred);
  }

  /** @brief Adds timed entries to the ordinary plot interface. */
  template <typename T>
  HS_HOT_FLASH_MEMBER PipelineRef(T &t, Erase<true>)
      : PipelineRef(t, Erase<false>{}) {
    cull = timed_cull_entry<T>;
    static constexpr TimedEntries entries{
        timed_world_entry<T>,
        timed_probe_entry<T>,
        timed_plot_entry<T>,
        transformed_plot_entry<T>,
        projected_plot_entry<T>,
        []() HS_HOT_INLINE {
          if constexpr (requires { T::rigid_world_stages; })
            return T::rigid_world_stages;
          else
            return false;
        }(),
        []() HS_HOT_INLINE {
          if constexpr (requires {
                          T::total_segment_margin;
                          T::any_2d_history;
                          T::any_reads_outside_band;
                        })
            return T::total_segment_margin <= 1 && !T::any_2d_history &&
                   !T::any_reads_outside_band;
          else
            return false;
        }(),
        [](void *pipeline) HS_HOT_INLINE {
          if constexpr (requires(T &p) { p.cull_bound(); })
            return static_cast<T *>(pipeline)->cull_bound();
          else
            return 0.0f;
        },
        [](void *pipeline) HS_HOT_INLINE {
          if constexpr (requires(T &p) { p.copy_count(); })
            return static_cast<T *>(pipeline)->copy_count();
          else
            return 1;
        },
        [](void *pipeline, int copy) HS_HOT_INLINE {
          if constexpr (requires(T &p) { p.select_copy(copy); })
            static_cast<T *>(pipeline)->select_copy(copy);
        }};
    timed_entries = &entries;
  }

public:
  /** @brief Moving-stage angular bound of a timed pipeline.
   * @return Bound in radians, or zero for an untimed ref. */
  float cull_bound() const {
    return timed_entries ? timed_entries->bound(ctx) : 0.0f;
  }
  /** @brief Number of independently walked spatial copies.
   * @return Copy count, or one for an untimed ref. */
  int copy_count() const {
    return timed_entries ? timed_entries->copies(ctx) : 1;
  }
  /** @brief Selects one fixed spatial copy.
   * @param copy Flattened copy index. */
  void select_copy(int copy) const {
    if (timed_entries)
      timed_entries->select(ctx, copy);
  }
  /** @brief Whether world transformations can run at controls.
   * @return True for a fully rigid, history-free world chain. */
  bool rigid_world_stages() const {
    return timed_entries && timed_entries->rigid;
  }

  /** @brief Whether screen stages permit a one-pixel dense-run cull pad.
   * @return True for bounded, history-free screen chains. */
  bool dense_screen_safe() const {
    return timed_entries && timed_entries->dense_screen_safe;
  }

  /** @brief Plots after already-evaluated world transformations.
   * @param cv Frame context.
   * @param v Transformed point.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @param prediction Optional part acceptance and boundary evaluator.
   * @return Actual unwrapped column. */
  HS_HOT_INLINE float
  plot_transformed_at(Canvas &cv, const math::Vector &v, const Pixel &c,
                      float age, float alpha, float t,
                      const Temporal::ColumnPrediction *prediction) const {
    return timed_entries->transformed_plot(ctx, cv, v, c, age, alpha, t,
                                           prediction);
  }

  /** @brief Emits a projected sample through screen stages.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @param prediction Optional landing window and boundary evaluator.
   * @return Actual unwrapped column. */
  HS_HOT_INLINE float
  plot_projected_at(Canvas &cv, math::PixelCoords p, const Pixel &c, float age,
                    float alpha, float t,
                    const Temporal::ColumnPrediction *prediction) const {
    return timed_entries->projected_plot(ctx, cv, p, c, age, alpha, t,
                                         prediction);
  }

  /** @brief Whether the erased pipeline leaves world positions unchanged. */
  bool world_transform_is_identity;
  /** @brief Whether a stage evaluates world motion over a frame. */
  bool has_temporal_motion;

  /**
   * @brief Wraps any object exposing 2D and 3D plot() methods.
   * @tparam T Pipeline type; must satisfy Plottable and be neither PipelineRef
   *         nor a direct-raster sink.
   * @param t Pipeline object whose address is stored; must outlive this ref.
   */
  template <typename T>
    requires(Plottable<T> && !std::same_as<std::decay_t<T>, PipelineRef> &&
             !pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t) : PipelineRef(t, Erase<false>{}) {}

  /** @brief Borrows a pipeline with explicit-time projection and plotting.
   * @param t Pipeline exposing probe and plot_at.
   * @details The pipeline must outlive this ref. */
  template <typename T>
    requires(!pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Temporal::Timed) : PipelineRef(t, Erase<true>{}) {}

  /**
   * @brief Erases a direct-raster sink for a draw into @p cv.
   * @param t Sink whose address is stored; must outlive this ref.
   * @param cv Canvas the draw writes into.
   * @details Checks prepared_for(@p cv): the sink writes through a framebuffer
   * base cached by prepare(), and a stale base is the buffer being scanned out.
   */
  template <typename T>
    requires(Plottable<T> && pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Canvas &cv) : PipelineRef(t, Erase<false>{}) {
    HS_CHECK(t.prepared_for(cv),
             "direct raster pipeline not prepared for this canvas");
  }

  /** @brief Borrows a prepared direct sink with explicit-time entry points.
   * @param t Direct sink.
   * @param cv Prepared canvas. */
  template <typename T>
    requires(pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, [[maybe_unused]] Canvas &cv, Temporal::Timed)
      : PipelineRef(t, Erase<true>{}) {
    assert(t.prepared_for(cv));
  }

  /** @brief Evaluates world stages without projection or history recording.
   * @param cv Frame context.
   * @param v Source point.
   * @param t Frame fraction.
   * @return Transformed world point. */
  math::Vector world_at(Canvas &cv, const math::Vector &v, float t) const {
    assert(timed_entries);
    return timed_entries->world(ctx, cv, v, t);
  }

  /** @brief Projects without plotting or recording.
   * @param cv Frame context.
   * @param v World point.
   * @param t Frame fraction.
   * @return Column and row through the world stages. */
  math::PixelCoords probe(Canvas &cv, const math::Vector &v, float t) const {
    assert(timed_entries);
    return timed_entries->probe(ctx, cv, v, t);
  }

  /** @brief Plots a world point at its supplied frame fraction.
   * @param cv Frame context.
   * @param v World point.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Source opacity.
   * @param t Frame fraction.
   * @param prediction Optional pre-write landing window.
   * @return Landed column, unwrapped near the prediction. */
  float plot_at(Canvas &cv, const math::Vector &v, const Pixel &c, float age,
                float alpha, float t,
                const Temporal::ColumnPrediction *prediction = nullptr) const {
    assert(timed_entries);
    return timed_entries->plot(ctx, cv, v, c, age, alpha, t, prediction);
  }

  /** @brief Solves and plots one geometry callback at each arm pass.
   * @tparam W,H Canvas dimensions.
   * @param cv Frame context.
   * @param geometry Pure position callback.
   * @param bound Combined geometry and world-stage motion bound.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Source opacity.
   * @param linear_pair The callback interpolates two lattice positions. */
  template <int W, int H, typename Geometry>
  void plot_timed(Canvas &cv, Geometry &&geometry, Temporal::MotionBound bound,
                  const Pixel &c, float age, float alpha,
                  bool linear_pair = false) const {
    auto position = [&](float t) HS_HOT_INLINE {
      Temporal::count_geometry<
          &Temporal::GeometryCounters::callback_evaluations>();
      return geometry(t);
    };
    if (!Temporal::temporal_sampling_enabled()) {
      plot_at(cv, position(1.0f), c, age, alpha, 1.0f);
      return;
    }
    auto project = [&](float t)
                       HS_HOT_INLINE { return probe(cv, position(t), t); };
    Temporal::count_geometry<&Temporal::GeometryCounters::end_probes>();
    const auto end = project(1.0f);
    const auto passes = Temporal::solve_capture<W>(
        project, end, bound.columns<W, H>(end.y), linear_pair);
    auto evaluate_boundary = [&](float t, math::PixelCoords &p, Pixel &color,
                                 float &sample_age, float &opacity)
                                 HS_HOT_INLINE {
                                   p = project(t);
                                   color = c;
                                   sample_age = age;
                                   opacity = alpha;
                                   return true;
                                 };
    Temporal::BoundaryEvaluator boundary(evaluate_boundary);
    for (int i = 0; i < passes.count; ++i) {
      const float t = passes.fraction[i];
      Temporal::ColumnPrediction prediction{
          end.x, std::numeric_limits<float>::infinity()};
      prediction.boundary = &boundary;
      plot_at(cv, position(t), c, age, alpha, t, &prediction);
    }
  }

  /**
   * @brief Plots a pixel at floating-point screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) const {
    plot2d(ctx, cv, x, y, c, age, alpha);
  }
  /**
   * @brief Plots a pixel at integer screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   * @details Preserves integer coordinates when forwarding to the wrapped
   * pipeline.
   */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age,
            float alpha) const {
    plot2d_int(ctx, cv, x, y, c, age, alpha);
  }
  /**
   * @brief Plots a pixel at a 3D world-space position.
   * @param cv Target canvas to draw into.
   * @param v World-space position to project and plot.
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   */
  void plot(Canvas &cv, const math::Vector &v, const Pixel &c, float age,
            float alpha) const {
    plot3d(ctx, cv, v, c, age, alpha);
  }

  /**
   * @brief Clip-cull query forwarded to the wrapped pipeline.
   * @param a Edge start (unit sphere point, pre-world-transform).
   * @param b Edge end (unit sphere point, pre-world-transform).
   * @param pb Optional planar basis for the edge (null = geodesic).
   * @param pred Rasterizer's screen-row-span vs clip-band test.
   * @return Whether any world-transformed copy of the edge could intersect the
   *         clip band (see Pipeline::could_intersect_clip).
   */
  bool could_intersect_clip(const math::Vector &a, const math::Vector &b,
                            const math::Basis *pb, CullEdgePredRef pred) const {
    return cull(ctx, a, b, pb, pred);
  }
};

/// Owned path function mapping a parameter to a sphere point.
using PlotFn = Fn<math::Vector(float), 16>;
// 16 B: holds two pointers on the 64-bit host, four on device.
/// Owned sprite draw callback: (canvas, opacity in [0, 1]).
using SpriteFn = Fn<void(Canvas &, float), 16>;
/// Owned timer callback, invoked with the frame's canvas.
using TimerFn = Fn<void(Canvas &), 16>;
/// Owned scalar function of one float.
using ScalarFn = Fn<float(float), 32>;
/// Easing curve mapping progress in [0, 1] to eased progress.
using EasingFn = float (*)(float);

/**
 * @brief Concept for a two-level orientation history: a sequence of frames,
 * each carrying its own sub-frame quaternion history.
 * @tparam T Candidate type; must expose length() and get(index) yielding a
 * frame.
 * @details Each frame needs a static CAPACITY plus its own length()/get(),
 * since deep_tween flattens both levels; a bare Orientation is rejected (use
 * tween()). The trail length is unsigned, the frame length signed.
 */
template <typename T>
concept Tweenable = requires(const T &t, size_t i) {
  { t.length() } -> std::unsigned_integral;
  { t.get(i).length() } -> std::signed_integral;
  { t.get(i).get(0) } -> std::convertible_to<math::Quaternion>;
  {
    std::remove_cvref_t<decltype(t.get(i))>::CAPACITY
  } -> std::convertible_to<size_t>;
};
