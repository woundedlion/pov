/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#pragma once
#include <concepts>
#include <functional>
#include "3dmath.h"
#include "canvas.h"

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
    return thunk_(ctx_, std::forward<Args>(args)...);
  }

  explicit operator bool() const { return thunk_ != nullptr; }
};

using ScreenTrailFn = FunctionRef<Color4(float, float, float)>;
using WorldTrailFn = FunctionRef<Color4(const Vector &, float)>;
using TransformFn = FunctionRef<Vector(const Vector &)>;
using FragmentShaderFn = FunctionRef<void(const Vector &, Fragment &)>;
using VertexShaderRef = FunctionRef<void(Fragment &)>;
using BlendFn = FunctionRef<Pixel(const Pixel &, const Pixel &)>;
using TweenFn = FunctionRef<void(const Quaternion &, float)>;
using VectorTweenFn = FunctionRef<void(const Vector &, float)>;

class PipelineRef {
  void *ctx_;
  void (*plot2d_)(void *, Canvas &, float, float, const Pixel &, float, float,
                  uint8_t);
  void (*plot3d_)(void *, Canvas &, const Vector &, const Pixel &, float, float,
                  uint8_t);

public:
  template <typename T> PipelineRef(T &t) : ctx_(std::addressof(t)) {
    plot2d_ = [](void *ctx, Canvas &cv, float x, float y, const Pixel &c,
                 float age, float alpha, uint8_t tag) {
      static_cast<T *>(ctx)->plot(cv, x, y, c, age, alpha, tag);
    };
    plot3d_ = [](void *ctx, Canvas &cv, const Vector &v, const Pixel &c,
                 float age, float alpha, uint8_t tag) {
      static_cast<T *>(ctx)->plot(cv, v, c, age, alpha, tag);
    };
  }

  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha, uint8_t tag = 0) const {
    plot2d_(ctx_, cv, x, y, c, age, alpha, tag);
  }
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) const {
    plot2d_(ctx_, cv, static_cast<float>(x), static_cast<float>(y), c, age,
            alpha, tag);
  }
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) const {
    plot3d_(ctx_, cv, v, c, age, alpha, tag);
  };
};

using PlotFn = Fn<Vector(float), 16>;
using SpriteFn = Fn<void(Canvas &, float), 48>;
using TimerFn = Fn<void(Canvas &), 16>;
using ScalarFn = Fn<float(float), 32>;
using HarmonicWaveFn = Fn<float(int, int, float, float), 8>;

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
