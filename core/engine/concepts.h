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
#include "core/render/spacetime.h"

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
using ScreenTrailRef = FunctionRef<Color4(float, float, float)>;
/// World trail colour: (position, t) with t the point's age in [0, 1].
using WorldTrailRef = FunctionRef<Color4(const math::Vector &, float)>;
/// Fragment shader: shades the fragment at a sphere position.
using FragmentShaderRef = FunctionRef<void(const math::Vector &, Fragment &)>;
/// Time-aware fragment shader: (position, fragment, sweep fraction of the
/// shaded column in [0, 1]).
using TimedShaderRef =
    FunctionRef<void(const math::Vector &, Fragment &, float)>;
/// Vertex shader: rewrites a fragment before projection.
using VertexShaderRef = FunctionRef<void(Fragment &)>;
/// Deferred per-control-point shader: shading registers and original pre-shader
/// position, invoked after projection so it cannot move the fragment.
using DeferredShaderRef =
    FunctionRef<void(FragmentRegisters, const math::Vector &)>;
/// Orientation history visitor: (orientation, t) with t in [0, 1].
using TweenRef = FunctionRef<void(const math::Quaternion &, float)>;
/// Position history visitor: (position, t) with t in [0, 1].
using VectorTweenRef = FunctionRef<void(const math::Vector &, float)>;
/// Clip-cull predicate: whether world-transformed edge (a, b), with an optional
/// planar basis, intersects the clip band.
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

#include "core/render/pipeline_ref.h"

// 16 B aliases hold two pointers on the 64-bit host, four on device.
/// Owned path function mapping a parameter to a sphere point.
using PlotFn = Fn<math::Vector(float), 16>;
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
