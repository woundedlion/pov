/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#pragma once

/**
 * @file inplace_function.h
 * @brief hs::inplace_function — heap-free, inline-storage callable for the
 *        host/WASM build, behind Fn<Sig,Cap>; modeled on SG14
 *        stdext::inplace_function.
 *
 * A closure that overflows the Capacity-byte buffer is a compile error.
 * Capacity counts bytes, so a pointer-capturing closure is wider on the 64-bit
 * host than on the 32-bit device.
 */

#include <cstddef>
#include <new>
#include <type_traits>
#include <utility>

namespace hs {

/**
 * @brief Diverges when an empty inplace_function is invoked.
 */
[[noreturn]] void inplace_function_empty_call();

// Alignment defaults to a pointer, not max_align_t. alignof(void *) is 8 on the
// 64-bit host but 4 on wasm32, so an 8-byte-aligned capture (double, int64_t)
// fails to compile only in WASM.
template <typename Signature, size_t Capacity = 16,
          size_t Alignment = alignof(void *)>
class inplace_function; // primary template intentionally undefined

namespace detail {

template <typename T> struct is_inplace_function : std::false_type {};

template <typename Signature, size_t Capacity, size_t Alignment>
struct is_inplace_function<inplace_function<Signature, Capacity, Alignment>>
    : std::true_type {};

/**
 * @brief Type-erased operation table for an inplace_function's stored
 *        callable; one shared instance per captured callable type.
 * @tparam R Call return type.
 * @tparam Args Call argument types.
 */
template <typename R, typename... Args> struct ipf_vtable {
  /// Calls the callable in a storage buffer.
  using invoke_ptr_t = R (*)(void *, Args &&...);
  /// Copy-constructs the callable from the second buffer into the first.
  using copy_ptr_t = void (*)(void *, const void *);
  /// Move-constructs the callable from the second buffer into the first.
  using move_ptr_t = void (*)(void *, void *);

  invoke_ptr_t invoke; ///< Invokes the stored callable.
  copy_ptr_t copy;     ///< Copies the stored callable into another buffer.
  move_ptr_t move;     ///< Moves the stored callable into another buffer.
};

/**
 * @brief ipf_vtable entries for a captured callable of type C placed in the
 *        inline buffer.
 * @tparam C Stored callable type.
 * @tparam R Call return type.
 * @tparam Args Call argument types.
 */
template <typename C, typename R, typename... Args> struct ipf_ops {
  // The byte buffer is not pointer-interconvertible with C; launder each access.
  /**
   * @brief Calls the C stored in `storage`.
   * @param storage Buffer holding a live C.
   * @param args Call arguments, forwarded.
   * @return The callable's result.
   */
  static R invoke(void *storage, Args &&...args) {
    C &callable = *std::launder(static_cast<C *>(storage));
    if constexpr (std::is_void_v<R>) {
      callable(std::forward<Args>(args)...);
    } else {
      return callable(std::forward<Args>(args)...);
    }
  }
  /**
   * @brief Copy-constructs a C into `dst`.
   * @param dst Uninitialized destination buffer.
   * @param src Buffer holding a live C.
   */
  static void copy(void *dst, const void *src) {
    ::new (dst) C(*std::launder(static_cast<const C *>(src)));
  }
  /**
   * @brief Move-constructs a C into `dst`; `src` stays live but moved-from.
   * @param dst Uninitialized destination buffer.
   * @param src Buffer holding a live C.
   */
  static void move(void *dst, void *src) {
    ::new (dst) C(std::move(*std::launder(static_cast<C *>(src))));
  }

  /// The shared vtable for C.
  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move};
};

/**
 * @brief ipf_vtable entries for an empty inplace_function: invoke traps;
 *        copy/move are no-ops.
 * @tparam R Call return type.
 * @tparam Args Call argument types.
 */
template <typename R, typename... Args> struct ipf_empty_ops {
  /**
   * @brief Traps via inplace_function_empty_call().
   * @return Never returns.
   */
  static R invoke(void *, Args &&...) { ::hs::inplace_function_empty_call(); }
  /** @brief No-op; an empty function holds no callable. */
  static void copy(void *, const void *) {}
  /** @brief No-op; an empty function holds no callable. */
  static void move(void *, void *) {}

  /// The shared empty-state vtable.
  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move};
};

} // namespace detail

/**
 * @brief Owning, heap-free type-erased callable with a fixed inline buffer.
 * @tparam R Return type of the call signature.
 * @tparam Args Argument types of the call signature.
 * @tparam Capacity Inline storage budget in bytes for the captured callable.
 * @tparam Alignment Inline storage alignment.
 * @details The stored callable must be trivially destructible, so
 *          inplace_function is too. operator() is const-qualified (the buffer
 *          is mutable).
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

  /**
   * @brief Copies the callable of `o` into this buffer.
   * @param o Source function.
   */
  inplace_function(const inplace_function &o) noexcept : vtable(o.vtable) {
    vtable->copy(storage, o.storage);
  }
  /**
   * @brief Moves the callable of `o` into this buffer and leaves `o` empty.
   * @param o Source function.
   */
  inplace_function(inplace_function &&o) noexcept : vtable(o.vtable) {
    vtable->move(storage, o.storage);
    o.vtable = empty_vtable();
  }

  inplace_function &operator=(const inplace_function &o) noexcept {
    if (this != &o) {
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
    *this = inplace_function(std::forward<C>(c));
    return *this;
  }
  inplace_function &operator=(std::nullptr_t) noexcept {
    vtable = empty_vtable();
    return *this;
  }

  /**
   * @brief Invokes the stored callable; traps if empty.
   * @param args Call arguments, forwarded.
   * @return The callable's result.
   */
  R operator()(Args... args) const {
    return vtable->invoke(storage, std::forward<Args>(args)...);
  }

  /** @brief True iff a callable is stored. */
  explicit operator bool() const noexcept { return vtable != empty_vtable(); }
};

} // namespace hs
