#pragma once
// ---------------------------------------------------------------------------
// hs::inplace_function — heap-free, inline-storage callable for the host/WASM
// build. A drop-in replacement for the std::function fallback behind Fn<Sig,Cap>
// (see platform.h), modeled on SG14 stdext::inplace_function.
//
// The buffer is fixed at Capacity bytes: a closure that overflows it is a hard
// *compile error*, not a silent heap allocation. Because Capacity counts bytes,
// a pointer-capturing closure is wider on the 64-bit host than on the 32-bit
// device; callsites tune Cap accordingly (see SpriteFn in concepts.h).
//
// Included only from platform.h's non-ARDUINO branch, after hs::check_fail is
// declared.
// ---------------------------------------------------------------------------

#include <cstddef>
#include <new>
#include <type_traits>
#include <utility>

namespace hs {

// Alignment defaults to a pointer, not max_align_t: the captures here are
// pointers/ints/floats/small PODs (max align == a pointer), so pointer alignment
// keeps the object to one pointer of overhead (e.g. a 16 B capture is 24 B total,
// matching teensy's footprint) instead of rounding every Fn up to 16 B and
// inflating Fn-bearing animation types past TimelineEvent::MAX_ANIM_SIZE. A
// rare over-aligned capture trips the alignof(D) <= Alignment static_assert below.
template <typename Signature, size_t Capacity = 16,
          size_t Alignment = alignof(void *)>
class inplace_function; // primary template intentionally undefined

namespace detail {

// Type-erased operation table, one shared instance per captured callable type.
template <typename R, typename... Args> struct ipf_vtable {
  using invoke_ptr_t = R (*)(void *, Args &&...);
  using copy_ptr_t = void (*)(void *, const void *);
  using move_ptr_t = void (*)(void *, void *);
  using destroy_ptr_t = void (*)(void *);

  invoke_ptr_t invoke;
  copy_ptr_t copy;
  move_ptr_t move;
  destroy_ptr_t destroy;
};

// Concrete operations for a captured callable C placed in the inline buffer.
template <typename C, typename R, typename... Args> struct ipf_ops {
  static R invoke(void *storage, Args &&...args) {
    return (*static_cast<C *>(storage))(std::forward<Args>(args)...);
  }
  static void copy(void *dst, const void *src) {
    ::new (dst) C(*static_cast<const C *>(src));
  }
  static void move(void *dst, void *src) {
    ::new (dst) C(std::move(*static_cast<C *>(src)));
    static_cast<C *>(src)->~C();
  }
  static void destroy(void *storage) { static_cast<C *>(storage)->~C(); }

  // C++17: static constexpr data members are implicitly inline, so this needs no
  // out-of-line definition and yields one shared vtable address per (C,R,Args).
  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move,
                                                &destroy};
};

// Empty-state operations: calling an empty inplace_function is a fail-fast trap,
// matching the project's no-silent-fallback habit (std::function would instead
// throw bad_function_call). copy/move/destroy are no-ops on the empty buffer.
template <typename R, typename... Args> struct ipf_empty_ops {
  static R invoke(void *, Args &&...) {
    // Unconditional [[noreturn]] call: no trailing return needed even for R!=void.
    ::hs::check_fail(__FILE__, __LINE__, "vtable != empty",
                     "empty hs::inplace_function called");
  }
  static void copy(void *, const void *) {}
  static void move(void *, void *) {}
  static void destroy(void *) {}

  static constexpr ipf_vtable<R, Args...> value{&invoke, &copy, &move,
                                                &destroy};
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
 *          than a heap allocation. operator() is const-qualified (the buffer is
 *          mutable) to match std::function, the fallback this replaces.
 */
template <typename R, typename... Args, size_t Capacity, size_t Alignment>
class inplace_function<R(Args...), Capacity, Alignment> {
  using vtable_t = detail::ipf_vtable<R, Args...>;

  static const vtable_t *empty_vtable() noexcept {
    return &detail::ipf_empty_ops<R, Args...>::value;
  }

  const vtable_t *vtable_ = empty_vtable();
  alignas(Alignment) mutable unsigned char storage_[Capacity];

  template <typename C>
  using is_self = std::is_same<std::decay_t<C>, inplace_function>;

public:
  /** @brief Constructs an empty function; calling it traps. */
  inplace_function() noexcept = default;
  /** @brief Constructs an empty function (nullptr overload). */
  inplace_function(std::nullptr_t) noexcept {}

  /**
   * @brief Constructs from any compatible callable, stored inline.
   * @param c Callable invocable as R(Args...); copied/moved into the buffer.
   */
  template <typename C,
            typename = std::enable_if_t<
                !is_self<C>::value &&
                std::is_invocable_r_v<R, std::decay_t<C> &, Args...>>>
  inplace_function(C &&c) {
    using D = std::decay_t<C>;
    static_assert(sizeof(D) <= Capacity,
                  "callable too large for inplace_function Capacity — raise the "
                  "Fn<Sig,Cap> capacity or shrink the capture");
    static_assert(alignof(D) <= Alignment,
                  "callable over-aligned for inplace_function storage");
    static_assert(std::is_copy_constructible_v<D>,
                  "inplace_function requires a copy-constructible callable");
    static_assert(std::is_nothrow_move_constructible_v<D>,
                  "inplace_function's move ctor/assign are noexcept and forward "
                  "to the stored type's move — a throwing move would "
                  "std::terminate. Store only nothrow-movable callables.");
    static_assert(std::is_nothrow_copy_constructible_v<D>,
                  "inplace_function's copy assignment destroys the old object "
                  "before copy-constructing the new one, so a throwing copy "
                  "would leave the buffer empty while vtable_ points at the new "
                  "type — UB on the next destroy. Store only nothrow-copyable "
                  "callables — any qualifying type is accepted (in practice "
                  "lambdas capturing PODs/pointers, which are trivially so).");
    ::new (storage_) D(std::forward<C>(c));
    vtable_ = &detail::ipf_ops<D, R, Args...>::value;
  }

  inplace_function(const inplace_function &o) : vtable_(o.vtable_) {
    vtable_->copy(storage_, o.storage_);
  }
  inplace_function(inplace_function &&o) noexcept : vtable_(o.vtable_) {
    vtable_->move(storage_, o.storage_);
    o.vtable_ = empty_vtable();
  }

  inplace_function &operator=(const inplace_function &o) {
    if (this != &o) {
      // Destroy-then-copy is safe only because the converting constructor's
      // is_nothrow_copy_constructible_v static_assert guarantees copy() cannot
      // throw — otherwise this would leave the buffer half-constructed.
      vtable_->destroy(storage_);
      vtable_ = o.vtable_;
      vtable_->copy(storage_, o.storage_);
    }
    return *this;
  }
  inplace_function &operator=(inplace_function &&o) noexcept {
    if (this != &o) {
      vtable_->destroy(storage_);
      vtable_ = o.vtable_;
      vtable_->move(storage_, o.storage_);
      o.vtable_ = empty_vtable();
    }
    return *this;
  }
  template <typename C, typename = std::enable_if_t<
                            !is_self<C>::value &&
                            std::is_invocable_r_v<R, std::decay_t<C> &, Args...>>>
  inplace_function &operator=(C &&c) {
    // Build a temporary so the capacity/copyability static_asserts live in one
    // place (the converting constructor), then move it in.
    *this = inplace_function(std::forward<C>(c));
    return *this;
  }
  inplace_function &operator=(std::nullptr_t) noexcept {
    vtable_->destroy(storage_);
    vtable_ = empty_vtable();
    return *this;
  }

  ~inplace_function() { vtable_->destroy(storage_); }

  /** @brief Invokes the stored callable; traps if empty. */
  R operator()(Args... args) const {
    return vtable_->invoke(storage_, std::forward<Args>(args)...);
  }

  /** @brief True iff a callable is stored. */
  explicit operator bool() const noexcept { return vtable_ != empty_vtable(); }
};

} // namespace hs
