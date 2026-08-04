/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#pragma once

/**
 * @file concepts.h
 * @brief Callable wrappers (FunctionRef, StoredFunctionRef) and the shader,
 *        trail and sprite callback aliases the render pipeline is written
 *        against.
 */

#include <concepts>
#include <cstddef>     // std::nullptr_t
#include <memory>      // std::addressof
#include <type_traits> // std::remove_cvref_t
#include <utility>     // std::forward
#include "math/3dmath.h"
#include "color/color.h"     // Pixel
#include "engine/platform.h" // Fn
#include <cassert>

struct Basis; // core/math/geometry.h; used only as const Basis* below
class Canvas; // core/render/canvas.h; used only behind references below

// ---------------------------------------------------------------------------
// Callable wrappers — two complementary types:
//
//   FunctionRef<Sig>   Non-owning, borrows the callable. Zero overhead.
//                      Use for parameters that are only invoked during the
//                      call (e.g. pipeline pass callbacks, shader functors).
//                      Must NOT outlive the referenced callable.
//
//   Fn<Sig, Cap>       Owning, inline storage (teensy::inplace_function on
//                      Teensy, hs::inplace_function on host/WASM — both
//                      heap-free). Use for stored
//                      callbacks that must persist beyond the creating scope
//                      (e.g. registered timers, sprite functions).
// ---------------------------------------------------------------------------

struct Fragment;
template <typename Signature> class FunctionRef;

/**
 * @brief Non-owning, type-erased reference to any callable matching Ret(Args...).
 * @tparam Ret Return type of the wrapped callable.
 * @tparam Args Argument types of the wrapped callable.
 * @details Borrows the callable via a void* context plus a thunk pointer; zero
 * heap allocation. The referenced callable must outlive the FunctionRef.
 */
template <typename Ret, typename... Args> class FunctionRef<Ret(Args...)> {
  void *ctx = nullptr;
  Ret (*thunk)(void *, Args...) = nullptr;

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
   * @details The function-pointer <-> void* round-trip is only
   * *conditionally-supported* by the standard ([expr.reinterpret.cast]) but holds
   * on every target this engine builds for (ARM Cortex-M7, x86-64, wasm32 all
   * share pointer width); the static_assert turns any future target where that
   * stops being true into a compile error. A null `func` produces an *empty* ref
   * (thunk stays null), matching a default-constructed FunctionRef, so a null
   * func cannot install a thunk that dereferences null on the first call.
   */
  FunctionRef(Ret (*func)(Args...)) noexcept
      : ctx(reinterpret_cast<void *>(func)) {
    static_assert(
        sizeof(func) == sizeof(void *),
        "FunctionRef requires function and object pointers to share a "
        "width (true on all supported targets)");
    if (func == nullptr)
      return;
    thunk = [](void *ptr, Args... args) -> Ret {
      return (reinterpret_cast<Ret (*)(Args...)>(ptr))(
          std::forward<Args>(args)...);
    };
  }

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
    thunk = [](void *ptr, Args... args) -> Ret {
      if constexpr (std::is_void_v<Ret>) {
        (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
      } else {
        return (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
      }
    };
  }

  /**
   * @brief Wraps a const lvalue callable (functor or lambda).
   * @tparam Callable Type of the callable; must be *const*-invocable with Args...
   * (the thunk invokes it through a `const Callable*`, so a mutable-only callable
   * is rejected here rather than failing inside the thunk) and not itself a
   * FunctionRef.
   * @param callable Const callable whose address is stored; must outlive this
   * ref. The const is cast away into ctx and restored in the thunk.
   * @details This overload also binds rvalues (temporary lambdas), so the
   * immediate-use borrow `take_callback([](...){ ... })` keeps working — the whole
   * purpose of a function_ref-style type. StoredFunctionRef refuses temporaries
   * for callables kept past the call.
   */
  template <typename Callable>
    requires std::is_invocable_r_v<Ret, const Callable &, Args...> &&
             (!std::is_base_of_v<FunctionRef, std::decay_t<Callable>>)
  FunctionRef(const Callable &callable) noexcept
      : ctx(const_cast<void *>(
            static_cast<const void *>(std::addressof(callable)))) {
    thunk = [](void *ptr, Args... args) -> Ret {
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
   * @details Per-pixel hot path (trail/transform callbacks): the null check is a
   * debug-only assert, deliberately NOT HS_CHECK to avoid an always-on per-call
   * branch on-device.
   */
  inline Ret operator()(Args... args) const {
    assert(thunk != nullptr &&
           "FunctionRef called on null/default-constructed ref");
    return thunk(ctx, std::forward<Args>(args)...);
  }

  /**
   * @brief Tests whether this FunctionRef refers to a callable.
   * @return True if a callable is bound, false if empty.
   */
  [[nodiscard]] explicit operator bool() const { return thunk != nullptr; }
};

/**
 * @brief A FunctionRef meant to be STORED past the call that builds it (e.g. a
 * class member invoked across many frames), not just borrowed for one call.
 * @tparam Signature The callable signature `Ret(Args...)`.
 * @details Identical to FunctionRef except it refuses to bind an rvalue temporary:
 * binding a temporary into something kept alive past the full expression dangles.
 * Storing sites use this type so the lifetime contract is enforced by the type
 * instead of a hand-rolled `= delete` at each site; plain FunctionRef stays right
 * for call-scoped parameters. Adds no data members, so it remains the same
 * two-pointer trivially-copyable payload as FunctionRef.
 */
template <typename Signature> class StoredFunctionRef;

template <typename Ret, typename... Args>
class StoredFunctionRef<Ret(Args...)> : public FunctionRef<Ret(Args...)> {
public:
  using FunctionRef<Ret(Args...)>::FunctionRef;

  // Reject rvalue temporaries the base would accept; the guards keep lvalue
  // callables and copy/move on the inherited ctors.
  template <typename Callable,
            typename = std::enable_if_t<
                !std::is_lvalue_reference_v<Callable> &&
                !std::is_same_v<std::decay_t<Callable>, StoredFunctionRef>>>
  StoredFunctionRef(Callable &&) = delete;
};

// These aliases are plain (borrow-only) FunctionRefs: they accept a temporary so
// call-scoped parameters stay ergonomic, and are used only as by-value/by-const-ref
// borrowed parameters. A class member that keeps a callable alive across frames
// must instead be typed StoredFunctionRef<Signature>, which `= delete`s the rvalue
// overload so a dangling bind to a temporary fails to compile (see the Timeline's
// tween members in animation.h for the canonical storing site).
using ScreenTrailFn = FunctionRef<Color4(float, float, float)>;
using WorldTrailFn = FunctionRef<Color4(const Vector &, float)>;
using FragmentShaderFn = FunctionRef<void(const Vector &, Fragment &)>;
using VertexShaderRef = FunctionRef<void(Fragment &)>;
// Deferred per-control-point shader: receives the (position-shaded) fragment
// and its original pre-shader position. Plot::ParticleSystem::draw runs it only
// for trails the segment cull keeps, skipping trails that render nothing.
using DeferredShaderRef = FunctionRef<void(Fragment &, const Vector &)>;
using TweenFn = FunctionRef<void(const Quaternion &, float)>;
using VectorTweenFn = FunctionRef<void(const Vector &, float)>;
// Rasterizer's clip-cull predicate: does the (world-transformed) edge a-b, with
// optional planar basis, intersect the clip band? Routed through the pipeline so
// world stages transform the edge before it is tested.
using CullEdgePredRef =
    FunctionRef<bool(const Vector &, const Vector &, const Basis *)>;

/**
 * @brief Deterministic ownership mask for dissolve transitions.
 * @details Ownership is hashed from an integer key pair, not from a pixel
 * coordinate. The wireframe path (Plot::Mesh::draw's edge-list overload) keys
 * owns() on an edge's endpoint indices and skips an unowned edge before any of
 * its geometry work; no solid-mesh scan consumes a mask. Two draws with the same threshold/salt
 * and opposite `invert` partition the domain exactly (every element owned by
 * one of them), which is what caps a two-mesh transition at one mesh's
 * rasterize cost per frame. The salt must derive from frame counters/seeds,
 * never wall time: the mask is part of the sim/device parity surface.
 */
struct DissolveMask {
  uint32_t threshold; /**< Owned fraction in [0, 65536]. */
  uint32_t salt;      /**< Per-frame/per-transition hash salt. */
  bool invert; /**< True owns the complement (elements at/above threshold). */

  /**
   * @brief Ownership of the element identified by the key pair (a, b).
   * @param a,b Key components; order-sensitive. The only consumer passes an
   *        edge's endpoint indices, emitting each edge once with one
   *        orientation.
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
 * @brief Non-owning, type-erased handle to a rasterizer pipeline.
 * @details Forwards plot() calls (2D screen-space or 3D world-space) to the
 * wrapped object's plot() methods. Like FunctionRef, it borrows the target and
 * must not outlive it. The erasure is a code-size choice: one indirect plot() call
 * per pixel, but the scanline machine instantiates once per <W,H> instead of once
 * per (shape x shader-lambda x filter-stack), which would explode device flash /
 * wasm size across the effects' distinct shader closures and filter-stack types.
 */
class PipelineRef {
  void *ctx;
  void (*plot2d)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot2d_int)(void *, Canvas &, int, int, const Pixel &, float, float);
  void (*plot3d)(void *, Canvas &, const Vector &, const Pixel &, float, float);
  bool (*cull)(void *, const Vector &, const Vector &, const Basis *,
               CullEdgePredRef);

public:
  /**
   * @brief Wraps any object exposing 2D and 3D plot() methods.
   * @tparam T Pipeline type; excluded from being PipelineRef itself.
   * @param t Pipeline object whose address is stored; must outlive this ref.
   * @details Excludes PipelineRef itself (mirroring FunctionRef): without this,
   * copying from a non-const lvalue binds this template (exact T& match) instead
   * of the implicit copy ctor, wrapping a ref-to-a-ref — an extra plot()
   * indirection per pixel and a dangling ctx if the source dies first.
   */
  template <typename T>
    requires(!std::same_as<std::decay_t<T>, PipelineRef>)
  PipelineRef(T &t) : ctx(std::addressof(t)) {
    plot2d = [](void *pipeline, Canvas &cv, float x, float y, const Pixel &c,
                float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot2d_int = [](void *pipeline, Canvas &cv, int x, int y, const Pixel &c,
                    float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot3d = [](void *pipeline, Canvas &cv, const Vector &v, const Pixel &c,
                float age, float alpha) {
      static_cast<T *>(pipeline)->plot(cv, v, c, age, alpha);
    };
    cull = [](void *pipeline, const Vector &a, const Vector &b, const Basis *pb,
              CullEdgePredRef pred) -> bool {
      // Real pipelines route the edge through their world stages; a bare
      // plot-provider (test stub) has no world transform, so test it directly.
      if constexpr (requires {
                      static_cast<T *>(pipeline)->could_intersect_clip(a, b, pb,
                                                                       pred);
                    })
        return static_cast<T *>(pipeline)->could_intersect_clip(a, b, pb, pred);
      else {
        // A filter pipeline (has any_crosses_segments) must answer the clip
        // query; only a bare plot stub legitimately falls through to pred. Catch
        // a could_intersect_clip signature drift that would else silently
        // degrade world-aware culling to raw-geometry culling.
        static_assert(
            !requires { T::any_crosses_segments; },
            "pipeline exposes any_crosses_segments but not "
            "could_intersect_clip (signature drift)");
        return pred(a, b, pb);
      }
    };
  }

  /**
   * @brief Plots a pixel at floating-point screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Normalized trail age in [0, 1].
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
   * @param age Normalized trail age in [0, 1].
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
   * @param age Normalized trail age in [0, 1].
   * @param alpha Coverage/opacity in [0, 1].
   */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
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
  bool could_intersect_clip(const Vector &a, const Vector &b, const Basis *pb,
                            CullEdgePredRef pred) const {
    return cull(ctx, a, b, pb, pred);
  }
};

using PlotFn = Fn<Vector(float), 16>;
using SpriteFn = Fn<void(Canvas &, float), 16>;
using TimerFn = Fn<void(Canvas &), 16>;
// 32: ScalarFn holds the wave/shift builders' captures, larger than 16 B.
using ScalarFn = Fn<float(float), 32>;
using EasingFn = float (*)(float);

/**
 * @brief Concept for a two-level orientation history: a sequence of frames,
 * each carrying its own sub-frame quaternion history.
 * @tparam T Candidate type; must expose length() and get(index) yielding a
 * frame.
 * @details Matches OrientationTrail (get -> Orientation). The element must
 * itself be frame-structured — a static CAPACITY plus its own length()/get() —
 * because deep_tween flattens both levels. A bare Orientation, whose get()
 * yields a Quaternion, is rejected here rather than deep in instantiation; use
 * tween() for that. The trail's own length() is consumed as an unsigned count —
 * a signed type could wrap negative into a huge loop bound; the frame's is
 * signed, matching Orientation's int index APIs.
 */
template <typename T>
concept Tweenable = requires(const T &t, size_t i) {
  { t.length() } -> std::unsigned_integral;
  { t.get(i).length() } -> std::signed_integral;
  { t.get(i).get(0) } -> std::convertible_to<Quaternion>;
  {
    std::remove_cvref_t<decltype(t.get(i))>::CAPACITY
  } -> std::convertible_to<size_t>;
};
