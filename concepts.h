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

  FunctionRef(Ret (*func)(Args...)) noexcept
      : ctx_(reinterpret_cast<void *>(func)) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (reinterpret_cast<Ret (*)(Args...)>(ptr))(
          std::forward<Args>(args)...);
    };
  }

  template <typename Callable>
    requires std::invocable<Callable, Args...>
  FunctionRef(Callable &callable) noexcept : ctx_(std::addressof(callable)) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
    };
  }

  template <typename Callable>
    requires std::invocable<Callable, Args...>
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

using PlotFn = std::function<Vector(float)>;
using SpriteFn = std::function<void(Canvas &, float)>;
using TimerFn = std::function<void(Canvas &)>;
using ScalarFn = std::function<float(float)>;
using HarmonicWaveFn = std::function<float(int, int, float, float)>;

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
