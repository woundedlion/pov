/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#pragma once
#include <concepts>
#include <functional>
#include "3dmath.h"
#include "canvas.h"

// ---------------------------------------------------------------------------
// Callable wrappers — two complementary types:
//
//   FunctionRef<Sig>   Non-owning, borrows the callable. Zero overhead.
//                      Use for parameters that are only invoked during the
//                      call (e.g. pipeline pass callbacks, shader functors).
//                      Must NOT outlive the referenced callable.
//
//   Fn<Sig, Cap>       Owning, inline storage (teensy::inplace_function on
//                      Teensy, std::function on WASM). Use for stored
//                      callbacks that must persist beyond the creating scope
//                      (e.g. registered timers, sprite functions).
// ---------------------------------------------------------------------------

struct Fragment;
struct Color4;
template <typename Signature> class FunctionRef;

template <typename Ret, typename... Args> class FunctionRef<Ret(Args...)> {
  void *ctx_ = nullptr;
  Ret (*thunk_)(void *, Args...) = nullptr;

public:
  FunctionRef() = default;
  FunctionRef(std::nullptr_t) {}

  // Explicit copy/move — prevent the generic Callable template from matching
  FunctionRef(const FunctionRef &other) noexcept = default;
  FunctionRef(FunctionRef &&other) noexcept = default;
  FunctionRef &operator=(const FunctionRef &other) noexcept = default;
  FunctionRef &operator=(FunctionRef &&other) noexcept = default;

  FunctionRef(Ret (*func)(Args...)) noexcept
      : ctx_(reinterpret_cast<void *>(func)) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (reinterpret_cast<Ret (*)(Args...)>(ptr))(
          std::forward<Args>(args)...);
    };
  }

  template <typename Callable>
    requires std::invocable<Callable, Args...> &&
             (!std::same_as<std::decay_t<Callable>, FunctionRef>)
  FunctionRef(Callable &callable) noexcept : ctx_(std::addressof(callable)) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
    };
  }

  template <typename Callable>
    requires std::invocable<Callable, Args...> &&
             (!std::same_as<std::decay_t<Callable>, FunctionRef>)
  FunctionRef(const Callable &callable) noexcept
      : ctx_(const_cast<void *>(
            static_cast<const void *>(std::addressof(callable)))) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (*static_cast<const Callable *>(ptr))(std::forward<Args>(args)...);
    };
  }

  inline Ret operator()(Args... args) const {
    // Per-pixel hot path (trail/transform callbacks): debug-only guard,
    // deliberately NOT HS_CHECK to avoid an always-on per-call branch on-device.
    assert(thunk_ != nullptr && "FunctionRef called on null/default-constructed ref");
    return thunk_(ctx_, std::forward<Args>(args)...);
  }

  [[nodiscard]] explicit operator bool() const { return thunk_ != nullptr; }
};

using ScreenTrailFn = FunctionRef<Color4(float, float, float)>;
using WorldTrailFn = FunctionRef<Color4(const Vector &, float)>;
using TransformFn = FunctionRef<Vector(const Vector &)>;
using SpaceTransformRef = FunctionRef<Vector(const Vector &)>;
using ColorTransformRef = FunctionRef<Pixel(const Pixel &, float)>;
using FragmentShaderFn = FunctionRef<void(const Vector &, Fragment &)>;
using VertexShaderRef = FunctionRef<void(Fragment &)>;
using BlendFn = FunctionRef<Pixel(const Pixel &, const Pixel &)>;
using TweenFn = FunctionRef<void(const Quaternion &, float)>;
using VectorTweenFn = FunctionRef<void(const Vector &, float)>;

// Non-owning, type-erased handle to a rasterizer pipeline. Forwards plot()
// calls (2D screen-space or 3D world-space) to the wrapped object's plot()
// methods. Like FunctionRef, it borrows the target and must not outlive it.
// Used by Plot::*/Scan::*::draw() so call sites can pass any pipeline without
// templating. This erasure is a deliberate code-size choice: erasing the
// pipeline here (and the shader, via FragmentShaderFn) at the draw() boundary
// costs one indirect plot() call per pixel, but instantiates the whole scanline
// machine once per <W,H> instead of once per (shape x shader-lambda x
// filter-stack) — the latter would explode device flash / wasm size across the
// 28 effects' distinct shader closures and variadic filter-stack types.
// Per-pixel inlining is intentionally traded for a bounded binary.
class PipelineRef {
  void *ctx_;
  void (*plot2d_)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot3d_)(void *, Canvas &, const Vector &, const Pixel &, float, float);

public:
  template <typename T> PipelineRef(T &t) : ctx_(std::addressof(t)) {
    plot2d_ = [](void *ctx, Canvas &cv, float x, float y, const Pixel &c,
                 float age, float alpha) {
      static_cast<T *>(ctx)->plot(cv, x, y, c, age, alpha);
    };
    plot3d_ = [](void *ctx, Canvas &cv, const Vector &v, const Pixel &c,
                 float age, float alpha) {
      static_cast<T *>(ctx)->plot(cv, v, c, age, alpha);
    };
  }

  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) const {
    plot2d_(ctx_, cv, x, y, c, age, alpha);
  }
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age,
            float alpha) const {
    plot2d_(ctx_, cv, static_cast<float>(x), static_cast<float>(y), c, age,
            alpha);
  }
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
            float alpha) const {
    plot3d_(ctx_, cv, v, c, age, alpha);
  };
};

using PlotFn = Fn<Vector(float), 16>;
using SpriteFn = Fn<void(Canvas &, float), 8>;
using TimerFn = Fn<void(Canvas &), 16>;
using ScalarFn = Fn<float(float), 32>;
using EasingFn = float (*)(float);

/**
 * @brief Concept for any object that maintains a history or sequence accessible
 * by index. Matches Orientation (get -> Quaternion) and OrientationTrail (get
 * -> Orientation).
 */
template <typename T>
concept Tweenable = requires(const T &t, size_t i) {
  { t.length() } -> std::convertible_to<size_t>;
  { t.get(i) }; // Return type is deduced (Quaternion or Orientation)
};
