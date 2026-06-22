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
//                      Teensy, hs::inplace_function on host/WASM — both
//                      heap-free). Use for stored
//                      callbacks that must persist beyond the creating scope
//                      (e.g. registered timers, sprite functions).
// ---------------------------------------------------------------------------

struct Fragment;
struct Color4;
template <typename Signature> class FunctionRef;

/**
 * @brief Non-owning, type-erased reference to any callable matching Ret(Args...).
 * @tparam Ret Return type of the wrapped callable.
 * @tparam Args Argument types of the wrapped callable.
 * @details Borrows the callable via a void* context plus a thunk pointer; zero
 * heap allocation. The referenced callable must outlive the FunctionRef.
 */
template <typename Ret, typename... Args> class FunctionRef<Ret(Args...)> {
  void *ctx_ = nullptr;
  Ret (*thunk_)(void *, Args...) = nullptr;

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
   * @param func Function pointer with signature Ret(Args...); stored in ctx_.
   * @details The function-pointer <-> void* round-trip is only
   * *conditionally-supported* by the standard ([expr.reinterpret.cast]), but it
   * holds on every target this engine builds for — ARM Cortex-M7, x86-64
   * (native tests), and wasm32 all use a single pointer width. The static_assert
   * turns any future target where that stops being true into a compile error
   * instead of silent corruption.
   */
  FunctionRef(Ret (*func)(Args...)) noexcept
      : ctx_(reinterpret_cast<void *>(func)) {
    static_assert(sizeof(func) == sizeof(void *),
                  "FunctionRef requires function and object pointers to share a "
                  "width (true on all supported targets)");
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (reinterpret_cast<Ret (*)(Args...)>(ptr))(
          std::forward<Args>(args)...);
    };
  }

  /**
   * @brief Wraps a non-const lvalue callable (functor or lambda).
   * @tparam Callable Type of the callable; must be invocable with Args... and
   * not itself a FunctionRef.
   * @param callable Callable whose address is stored; must outlive this ref.
   */
  template <typename Callable>
    requires std::invocable<Callable, Args...> &&
             (!std::same_as<std::decay_t<Callable>, FunctionRef>)
  FunctionRef(Callable &callable) noexcept : ctx_(std::addressof(callable)) {
    thunk_ = [](void *ptr, Args... args) -> Ret {
      return (*static_cast<Callable *>(ptr))(std::forward<Args>(args)...);
    };
  }

  /**
   * @brief Wraps a const lvalue callable (functor or lambda).
   * @tparam Callable Type of the callable; must be invocable with Args... and
   * not itself a FunctionRef.
   * @param callable Const callable whose address is stored; must outlive this
   * ref. The const is cast away into ctx_ and restored in the thunk.
   * @details This overload deliberately binds rvalues (temporary lambdas) too,
   * so the idiomatic immediate-use borrow — `take_callback([](...){ ... })` for
   * a parameter invoked only during the call — keeps working. That is the whole
   * purpose of a function_ref-style type; an `= delete`d rvalue overload would
   * outlaw the safe common case to catch the rarer store-past-its-lifetime
   * misuse, which the class-level lifetime contract already documents.
   */
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

  /**
   * @brief Invokes the wrapped callable.
   * @param args Arguments forwarded to the callable.
   * @return Result of the wrapped callable.
   * @details Per-pixel hot path (trail/transform callbacks): the null check is a
   * debug-only assert, deliberately NOT HS_CHECK to avoid an always-on per-call
   * branch on-device.
   */
  inline Ret operator()(Args... args) const {
    assert(thunk_ != nullptr && "FunctionRef called on null/default-constructed ref");
    return thunk_(ctx_, std::forward<Args>(args)...);
  }

  /**
   * @brief Tests whether this FunctionRef refers to a callable.
   * @return True if a callable is bound, false if empty.
   */
  [[nodiscard]] explicit operator bool() const { return thunk_ != nullptr; }
};

/**
 * @brief A FunctionRef meant to be STORED past the call that builds it (e.g. a
 * class member invoked across many frames), not just borrowed for one call.
 * @tparam Signature The callable signature `Ret(Args...)`.
 * @details Identical to FunctionRef except it refuses to bind an rvalue
 * temporary. FunctionRef's const-lvalue ctor deliberately accepts temporaries so
 * call-scoped parameter borrows stay ergonomic (see FunctionRef), but binding a
 * temporary into something kept alive past the full expression is a dangling
 * reference. Storing sites use this type so the lifetime contract is enforced by
 * the type instead of a hand-rolled `= delete` at each site; plain FunctionRef
 * stays the right choice for call-scoped parameters. Adds no data members, so it
 * remains the same two-pointer trivially-copyable payload as FunctionRef.
 */
template <typename Signature> class StoredFunctionRef;

template <typename Ret, typename... Args>
class StoredFunctionRef<Ret(Args...)> : public FunctionRef<Ret(Args...)> {
public:
  using FunctionRef<Ret(Args...)>::FunctionRef;

  // Reject rvalue temporaries the base would otherwise accept: a stored ref must
  // outlive the call that builds it. The lvalue-reference and self-type guards
  // keep lvalue callables and StoredFunctionRef copy/move binding to the
  // inherited ctors instead of this deleted overload.
  template <typename Callable,
            typename = std::enable_if_t<
                !std::is_lvalue_reference_v<Callable> &&
                !std::is_same_v<std::decay_t<Callable>, StoredFunctionRef>>>
  StoredFunctionRef(Callable &&) = delete;
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

/**
 * @brief Non-owning, type-erased handle to a rasterizer pipeline.
 * @details Forwards plot() calls (2D screen-space or 3D world-space) to the
 * wrapped object's plot() methods. Like FunctionRef, it borrows the target and
 * must not outlive it. Used by the Plot and Scan draw() entry points so call
 * sites can pass any pipeline without templating. This erasure is a deliberate code-size
 * choice: erasing the pipeline here (and the shader, via FragmentShaderFn) at
 * the draw() boundary costs one indirect plot() call per pixel, but instantiates
 * the whole scanline machine once per <W,H> instead of once per (shape x
 * shader-lambda x filter-stack) — the latter would explode device flash / wasm
 * size across the 27 effects' distinct shader closures and variadic filter-stack
 * types. Per-pixel inlining is intentionally traded for a bounded binary.
 */
class PipelineRef {
  void *ctx_;
  void (*plot2d_)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot3d_)(void *, Canvas &, const Vector &, const Pixel &, float, float);

public:
  /**
   * @brief Wraps any object exposing 2D and 3D plot() methods.
   * @tparam T Pipeline type; excluded from being PipelineRef itself.
   * @param t Pipeline object whose address is stored; must outlive this ref.
   * @details Excludes PipelineRef itself (mirroring FunctionRef): without this,
   * copying from a non-const lvalue binds this template (exact T& match) instead
   * of the implicit copy ctor, wrapping a ref-to-a-ref — an extra plot()
   * indirection per pixel and a dangling ctx_ if the source dies first.
   */
  template <typename T>
    requires(!std::same_as<std::decay_t<T>, PipelineRef>)
  PipelineRef(T &t) : ctx_(std::addressof(t)) {
    plot2d_ = [](void *ctx, Canvas &cv, float x, float y, const Pixel &c,
                 float age, float alpha) {
      static_cast<T *>(ctx)->plot(cv, x, y, c, age, alpha);
    };
    plot3d_ = [](void *ctx, Canvas &cv, const Vector &v, const Pixel &c,
                 float age, float alpha) {
      static_cast<T *>(ctx)->plot(cv, v, c, age, alpha);
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
    plot2d_(ctx_, cv, x, y, c, age, alpha);
  }
  /**
   * @brief Plots a pixel at integer screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Normalized trail age in [0, 1].
   * @param alpha Coverage/opacity in [0, 1].
   * @details Promotes the integer coordinates to float and forwards to the 2D
   * plot thunk.
   */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age,
            float alpha) const {
    plot2d_(ctx_, cv, static_cast<float>(x), static_cast<float>(y), c, age,
            alpha);
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
    plot3d_(ctx_, cv, v, c, age, alpha);
  };
};

using PlotFn = Fn<Vector(float), 16>;
// One Cap for every platform — no per-arch split. The sprite draw closures all
// capture [this, idx, slot]; on a 64-bit host that is this(8) + two ints = 16 B,
// and pointer alignment rounds any [this, data] capture up to 16 regardless, so
// 16 is the floor there. 32-bit targets (device + WASM) pack the same capture
// into 12 B, so 16 covers them with room to spare. It is NOT 8: every target —
// the device included — builds these effects (the .ino factory runs init() ->
// spawn_sprite()), and 12 B never fit 8; the historical "8" predates them.
using SpriteFn = Fn<void(Canvas &, float), 16>;
using TimerFn = Fn<void(Canvas &), 16>;
using ScalarFn = Fn<float(float), 32>;
using EasingFn = float (*)(float);

/**
 * @brief Concept for any object that maintains a history or sequence accessible
 * by index.
 * @tparam T Candidate type; must expose length() and get(index).
 * @details Matches Orientation (get -> Quaternion) and OrientationTrail (get ->
 * Orientation). The return type of get() is left deduced.
 */
template <typename T>
concept Tweenable = requires(const T &t, size_t i) {
  { t.length() } -> std::convertible_to<size_t>;
  { t.get(i) }; // Return type is deduced (Quaternion or Orientation)
};
