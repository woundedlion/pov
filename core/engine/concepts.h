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
#include "math/3dmath.h"
#include "color/pixel.h"       // Pixel
#include "platform/platform.h" // Fn

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

  /** @brief Stores a non-throwing function pointer by value. */
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

template <typename Ret, typename... Args>
class StoredFunctionRef<Ret(Args...)> : public FunctionRef<Ret(Args...)> {
public:
  using FunctionRef<Ret(Args...)>::FunctionRef;

  StoredFunctionRef() = default;
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
using ScreenTrailFn = FunctionRef<Color4(float, float, float)>;
using WorldTrailFn = FunctionRef<Color4(const math::Vector &, float)>;
using FragmentShaderFn = FunctionRef<void(const math::Vector &, Fragment &)>;
using VertexShaderRef = FunctionRef<void(Fragment &)>;
// Deferred per-control-point shader: receives the (position-shaded) fragment's
// shading registers and its original pre-shader position. It runs after
// projection, so it cannot move the fragment.
using DeferredShaderRef =
    FunctionRef<void(FragmentRegisters, const math::Vector &)>;
using TweenFn = FunctionRef<void(const math::Quaternion &, float)>;
using VectorTweenFn = FunctionRef<void(const math::Vector &, float)>;
// Rasterizer's clip-cull predicate: does the (world-transformed) edge a-b, with
// optional planar basis, intersect the clip band? Routed through the pipeline so
// world stages transform the edge before it is tested.
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

/** @brief Dispatches a clip query through world stages or a bare plot provider. */
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
  struct Erase {};

  void *ctx;
  void (*plot2d)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot2d_int)(void *, Canvas &, int, int, const Pixel &, float, float);
  void (*plot3d)(void *, Canvas &, const math::Vector &, const Pixel &, float,
                 float);
  bool (*cull)(void *, const math::Vector &, const math::Vector &,
               const math::Basis *, CullEdgePredRef);

  template <typename T>
  PipelineRef(T &t, Erase)
      : ctx(std::addressof(t)), world_transform_is_identity([] {
          if constexpr (requires { T::world_transform_is_identity; })
            return T::world_transform_is_identity;
          else
            return true;
        }()) {
    plot2d = [](void *pipeline, Canvas &cv, float x, float y, const Pixel &c,
                float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot2d_int = [](void *pipeline, Canvas &cv, int x, int y, const Pixel &c,
                    float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot3d = [](void *pipeline, Canvas &cv, const math::Vector &v,
                const Pixel &c, float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, v, c, age, alpha);
    };
    cull = [](void *pipeline, const math::Vector &a, const math::Vector &b,
              const math::Basis *pb, CullEdgePredRef pred) -> bool {
      return pipeline_could_intersect_clip(*static_cast<T *>(pipeline), a, b,
                                           pb, pred);
    };
  }

public:
  /** @brief Whether the erased pipeline leaves world positions unchanged. */
  bool world_transform_is_identity;

  /**
   * @brief Wraps any object exposing 2D and 3D plot() methods.
   * @tparam T Pipeline type; must satisfy Plottable and be neither PipelineRef
   *         nor a direct-raster sink.
   * @param t Pipeline object whose address is stored; must outlive this ref.
   */
  template <typename T>
    requires(Plottable<T> && !std::same_as<std::decay_t<T>, PipelineRef> &&
             !pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t) : PipelineRef(t, Erase{}) {}

  /**
   * @brief Erases a direct-raster sink for a draw into @p cv.
   * @param t Sink whose address is stored; must outlive this ref.
   * @param cv Canvas the draw writes into.
   * @details Checks prepared_for(@p cv): the sink writes through a framebuffer
   * base cached by prepare(), and a stale base is the buffer being scanned out.
   */
  template <typename T>
    requires(Plottable<T> && pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Canvas &cv) : PipelineRef(t, Erase{}) {
    HS_CHECK(t.prepared_for(cv),
             "direct raster pipeline not prepared for this canvas");
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

using PlotFn = Fn<math::Vector(float), 16>;
// 16 B: holds two pointers on the 64-bit host, four on device.
using SpriteFn = Fn<void(Canvas &, float), 16>;
using TimerFn = Fn<void(Canvas &), 16>;
using ScalarFn = Fn<float(float), 32>;
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
