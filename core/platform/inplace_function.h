/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#pragma once

/**
 * @file inplace_function.h
 * @brief hs::inplace_function — heap-free, inline-storage callable for the
 *        host/WASM build, behind Fn<Sig,Cap> (see platform.h); modeled on SG14
 *        stdext::inplace_function.
 *
 * The buffer is fixed at Capacity bytes: a closure that overflows it is a hard
 * *compile error*, not a silent heap allocation. Because Capacity counts bytes,
 * a pointer-capturing closure is wider on the 64-bit host than on the 32-bit
 * device; callsites pick a fixed Cap with headroom for the wider host closure
 * (see SpriteFn in concepts.h).
 *
 * Included only from platform.h's non-ARDUINO branch, after hs::check_fail is
 * declared.
 */

#include <cstddef>
#include <new>
#include <type_traits>
#include <utility>

namespace hs {

/**
 * @brief Diverges when an empty inplace_function is invoked.
 * @details One out-of-line routine every instantiation's empty vtable entry
 * calls, so the trapping empty state costs a call per signature rather than a
 * formatted call site.
 */
[[noreturn]] void inplace_function_empty_call();

// Alignment defaults to a pointer, not max_align_t: the captures here are
// pointers/ints/floats/small PODs (max align == a pointer), so pointer alignment
// keeps the object to one pointer of overhead instead of rounding every Fn up to
// 16 B and inflating Fn-bearing animation types past TimelineEvent::MAX_ANIM_SIZE.
// A rare over-aligned capture trips the alignof(D) <= Alignment static_assert
// below — and the threshold is target-dependent, mirroring the Capacity skew:
// alignof(void *) is 8 on the 64-bit host but 4 on wasm32, so a capture wanting
// 8-byte alignment (double, int64_t) compiles natively and fails only in WASM.
template <typename Signature, size_t Capacity = 16,
          size_t Alignment = alignof(void *)>
class inplace_function; // primary template intentionally undefined

namespace detail {

template <typename T> struct is_inplace_function : std::false_type {};

template <typename Signature, size_t Capacity, size_t Alignment>
struct is_inplace_function<inplace_function<Signature, Capacity, Alignment>>
    : std::true_type {};

// Type-erased operation table, one shared instance per captured callable type.
template <typename R, typename... Args> struct ipf_vtable {
  using invoke_ptr_t = R (*)(void *, Args &&...);
  using copy_ptr_t = void (*)(void *, const void *);
  using move_ptr_t = void (*)(void *, void *);

  invoke_ptr_t invoke;
  copy_ptr_t copy;
  move_ptr_t move;
};

// Concrete operations for a captured callable C placed in the inline buffer.
template <typename C, typename R, typename... Args> struct ipf_ops {
  // The C was placement-new'd into an unsigned char buffer, which is not
  // pointer-interconvertible with it; std::launder makes each access refer to
  // the created object rather than the byte array.
  static R invoke(void *storage, Args &&...args) {
    C &callable = *std::launder(static_cast<C *>(storage));
    if constexpr (std::is_void_v<R>) {
      callable(std::forward<Args>(args)...);
    } else {
      return callable(std::forward<Args>(args)...);
    }
  }
  static void copy(void *dst, const void *src) {
    ::new (dst) C(*std::launder(static_cast<const C *>(src)));
  }
  static void move(void *dst, void *src) {
    ::new (dst) C(std::move(*std::launder(static_cast<C *>(src))));
  }

  // C++17: static constexpr data members are implicitly inline, so this needs no
  // out-of-line definition and yields one shared vtable address per (C,R,Args).
  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move};
};

// Empty-state operations: calling an empty inplace_function is a fail-fast trap,
// matching the project's no-silent-fallback habit (std::function would instead
// throw bad_function_call). copy/move are no-ops on the empty buffer.
template <typename R, typename... Args> struct ipf_empty_ops {
  static R invoke(void *, Args &&...) {
    // Unconditional [[noreturn]] call: no trailing return needed even for R!=void.
    ::hs::inplace_function_empty_call();
  }
  static void copy(void *, const void *) {}
  static void move(void *, void *) {}

  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move};
};

} // namespace detail

/**
 * @brief Owning, heap-free type-erased callable with a fixed inline buffer.
 * @tparam R Return type of the call signature.
 * @tparam Args Argument types of the call signature.
 * @tparam Capacity Inline storage budget in bytes for the captured callable.
 * @tparam Alignment Inline storage alignment.
 * @details Construction static_asserts that the decayed callable fits within
 *          Capacity/Alignment, so an oversized closure is a compile error rather
 *          than a heap allocation. The stored callable must be trivially
 *          destructible, which makes inplace_function trivially destructible in
 *          turn, so an ArenaVector may hold one. operator() is const-qualified
 *          (the buffer is mutable) to match std::function, the fallback this
 *          replaces.
 */
template <typename R, typename... Args, size_t Capacity, size_t Alignment>
class inplace_function<R(Args...), Capacity, Alignment> {
  using vtable_t = detail::ipf_vtable<R, Args...>;

  static const vtable_t *empty_vtable() noexcept {
    return &detail::ipf_empty_ops<R, Args...>::value;
  }

  const vtable_t *vtable = empty_vtable();
  alignas(Alignment) mutable unsigned char storage[Capacity];

public:
  // No ctor initializes storage: the empty vtable's copy/move are no-ops and its
  // invoke traps, so no operation ever reads an empty function's buffer; the
  // value-carrying ctor placement-news over it.

  /** @brief Constructs an empty function; calling it traps. */
  inplace_function() noexcept {}
  /** @brief Constructs an empty function (nullptr overload). */
  inplace_function(std::nullptr_t) noexcept {}

  /**
   * @brief Constructs from any compatible callable, stored inline.
   * @param c Callable invocable as R(Args...); copied/moved into the buffer.
   */
  template <typename C,
            typename = std::enable_if_t<
                !detail::is_inplace_function<std::decay_t<C>>::value &&
                std::is_invocable_r_v<R, std::decay_t<C> &, Args...>>>
  inplace_function(C &&c) {
    using D = std::decay_t<C>;
    static_assert(
        sizeof(D) <= Capacity,
        "callable too large for inplace_function Capacity — raise the "
        "Fn<Sig,Cap> capacity or shrink the capture");
    static_assert(alignof(D) <= Alignment,
                  "callable over-aligned for inplace_function storage");
    static_assert(
        std::is_trivially_destructible_v<D>,
        "inplace_function never runs the stored callable's "
        "destructor, which is also what keeps inplace_function itself "
        "trivially destructible for ArenaVector's destructor-skipping "
        "contract — store only trivially destructible callables.");
    static_assert(std::is_copy_constructible_v<D>,
                  "inplace_function requires a copy-constructible callable");
    static_assert(
        std::is_nothrow_move_constructible_v<D>,
        "inplace_function's move ctor/assign are noexcept and forward "
        "to the stored type's move — a throwing move would "
        "std::terminate. Store only nothrow-movable callables.");
    static_assert(std::is_nothrow_copy_constructible_v<D>,
                  "inplace_function's copy assignment placement-news over the "
                  "old object, so a throwing copy would leave the buffer unset "
                  "while vtable points at the new type — UB on the next call. "
                  "Store only nothrow-copyable callables — any qualifying type "
                  "is accepted (in practice lambdas capturing PODs/pointers, "
                  "which are trivially so).");
    ::new (storage) D(std::forward<C>(c));
    vtable = &detail::ipf_ops<D, R, Args...>::value;
  }

  inplace_function(const inplace_function &o) noexcept : vtable(o.vtable) {
    vtable->copy(storage, o.storage);
  }
  inplace_function(inplace_function &&o) noexcept : vtable(o.vtable) {
    vtable->move(storage, o.storage);
    o.vtable = empty_vtable();
  }

  inplace_function &operator=(const inplace_function &o) noexcept {
    if (this != &o) {
      // Trivially destructible: the incoming callable is placement-new'd over
      // the old one, which needs no destructor call.
      vtable = o.vtable;
      vtable->copy(storage, o.storage);
    }
    return *this;
  }
  inplace_function &operator=(inplace_function &&o) noexcept {
    if (this != &o) {
      vtable = o.vtable;
      vtable->move(storage, o.storage);
      o.vtable = empty_vtable();
    }
    return *this;
  }
  template <typename C,
            typename = std::enable_if_t<
                !detail::is_inplace_function<std::decay_t<C>>::value &&
                std::is_invocable_r_v<R, std::decay_t<C> &, Args...>>>
  inplace_function &operator=(C &&c) {
    // Build a temporary so the capacity/copyability static_asserts live in one
    // place (the converting constructor), then move it in.
    *this = inplace_function(std::forward<C>(c));
    return *this;
  }
  inplace_function &operator=(std::nullptr_t) noexcept {
    vtable = empty_vtable();
    return *this;
  }

  /** @brief Invokes the stored callable; traps if empty. */
  R operator()(Args... args) const {
    return vtable->invoke(storage, std::forward<Args>(args)...);
  }

  /** @brief True iff a callable is stored. */
  explicit operator bool() const noexcept { return vtable != empty_vtable(); }
};

} // namespace hs
