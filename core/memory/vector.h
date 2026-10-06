/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/memory.h.

// ============================================================================
// Arena Structures
// ============================================================================

/**
 * @brief Logs an ArenaVector grow that abandons its previous block.
 * @param bytes Size of the abandoned block.
 * @param old_capacity Element capacity before the grow.
 * @param new_capacity Element capacity after the grow.
 * @details Also feeds note_arena_vector_abandon(), so the OOM trap's running
 * total covers grows as well as move-assignments. Out-of-line and non-template
 * so the device image carries one copy for every element type. The leak is
 * permanent until the arena is reset, so the line ships in release: without it a
 * persistent-arena grow surfaces only as a later, innocent-looking allocation
 * trapping on OOM.
 */
FLASHMEM void log_arena_vector_grow(size_t bytes, size_t old_capacity,
                                    size_t new_capacity);

/**
 * @brief Whether T is a sanctioned inline callable safe to store in an
 * ArenaVector despite not being trivially destructible.
 * @details Only the device's teensy::inplace_function needs the exemption; the
 * host/WASM hs::inplace_function is trivially destructible and passes on its own.
 * Stored captures remain subject to ArenaVector's element destructor contract.
 */
template <typename T> struct is_arena_inplace_fn : std::false_type {};
#ifdef ARDUINO
template <typename R, typename... Args, size_t Cap, size_t Align>
struct is_arena_inplace_fn<teensy::inplace_function<R(Args...), Cap, Align>>
    : std::true_type {};
#endif

/**
 * @brief Arena-backed vector with a capacity fixed between bind() calls.
 *        Move-only.
 * @tparam T Element type; must satisfy the element destructor contract below.
 * @details CAPACITY CONTRACT: appending never grows the block — push_back()
 * traps once element_count reaches capacity(). Only bind() changes capacity,
 * and a grow there allocates a fresh block and abandons the old one until the
 * arena is reset (log_arena_vector_grow() reports the leaked bytes in release).
 * Element access follows std::vector: operator[] and back() require a valid
 * index/non-empty vector and check that precondition only in debug builds.
 *
 * ELEMENT DESTRUCTOR CONTRACT: ArenaVector does NOT run element
 * destructors — clear(), move, move-assign and going out of scope all leave
 * stored elements un-destructed. Storage is owned and reclaimed by the arena
 * (reset/compaction), not by this handle, and an arena can be reset out from
 * under a still-live ArenaVector (see bind()'s stale-binding handling), so a
 * handle-driven destructor pass could run on already-reclaimed memory. Only store
 * types whose destructor need not run for correctness: trivially-destructible
 * PODs, or Fn<> whose stored captures are themselves trivial. A type owning
 * heap/handles outside the arena must not be stored here — notably a raw
 * std::function, which would leak when this handle skips its destructor.
 */
template <typename T> class ArenaVector {
  // ArenaSpan borrows our backing data and (in debug builds) our arena
  // generation stamp for its own use-after-free check.
  template <typename U> friend class ArenaSpan;

private:
  T *elements;             /**< Pointer to the arena-allocated backing block. */
  size_t element_count;    /**< Number of constructed elements. */
  size_t element_capacity; /**< Maximum element count the block can hold. */
  bool bound = false; /**< Whether the vector has been bound to an arena. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the block was allocated. */
  /**
   * @brief Per-vector counter bumped on every bind(), including storage reuse.
   * @details A bind resets the element count and can replace the backing block
   * without changing the arena generation. This counter invalidates spans
   * captured before either path.
   */
  uint32_t rebind_generation = 0;

  /**
   * @brief Debug-only use-after-free check against the source arena.
   * @details Asserts if the source arena was reset out from under this vector,
   * or rewound (set_offset/ScratchScope) below the backing block — a rewind
   * reclaims the bytes without bumping the generation.
   */
  void check_alive() const {
    const size_t bytes = element_capacity * sizeof(T);
    assert(!stamp.arena_reset() && "ArenaVector use-after-free!");
    assert(!stamp.block_uncovered(elements, bytes) &&
           "ArenaVector use-after-free (arena rewound below block)!");
    assert(!stamp.block_reissued(elements, bytes) &&
           "ArenaVector use-after-free (block reclaimed by a rewind and "
           "reissued)!");
  }
#else
  /**
   * @brief No-op use-after-free check in release builds.
   */
  void check_alive() const {}
#endif

  /**
   * @brief Asserts that the vector has been bound to an arena.
   */
  void check_bound() const {
    assert(bound && "Attempted to access unbound ArenaVector!");
  }

  /**
   * @brief Transfers @p other's storage and bookkeeping into this vector,
   *        leaving @p other in a pristine unbound state.
   */
  void steal_from(ArenaVector &other) noexcept {
    elements = other.elements;
    element_count = other.element_count;
    element_capacity = other.element_capacity;
    bound = other.bound;
#ifndef NDEBUG
    stamp = other.stamp;
    rebind_generation = other.rebind_generation;
#endif
    other.elements = nullptr;
    other.element_count = 0;
    other.element_capacity = 0;
    other.bound = false;
#ifndef NDEBUG
    other.stamp.clear();
    // rebind_generation stays: spans snapshotted it and still view live data.
#endif
  }

public:
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T &;
  using const_reference = const T &;
  using pointer = T *;
  using const_pointer = const T *;
  using iterator = T *;
  using const_iterator = const T *;

  /**
   * @brief Default-constructs an unbound vector.
   * @details Must call bind() before use.
   */
  ArenaVector() : elements(nullptr), element_count(0), element_capacity(0) {}

  /**
   * @brief Deleted copy constructor.
   * @details Implicit shallow copying is disabled to prevent memory aliasing.
   */
  ArenaVector(const ArenaVector &) = delete;
  /**
   * @brief Deleted copy assignment.
   * @return Reference to this (never invoked).
   * @details Implicit shallow copying is disabled to prevent memory aliasing.
   */
  ArenaVector &operator=(const ArenaVector &) = delete;

  /**
   * @brief Move constructor.
   * @param other Source vector; left in a pristine unbound state.
   */
  ArenaVector(ArenaVector &&other) noexcept { steal_from(other); }

  /**
   * @brief Move assignment.
   * @param other Source vector; left in a pristine unbound state.
   * @return Reference to this.
   */
  ArenaVector &operator=(ArenaVector &&other) noexcept {
    if (this != &other) {
      // Overwriting a bound handle drops its block; the arena reclaims it only
      // on the next reset.
      if (bound && element_capacity > 0)
        note_arena_vector_abandon(element_capacity * sizeof(T));
      steal_from(other);
    }
    return *this;
  }

#ifndef NDEBUG
  /** @brief Allocation/reuse generation for debug lifetime diagnostics. */
  uint32_t debug_binding_generation() const { return rebind_generation; }
#endif

  /**
   * @brief Constructs and binds the vector with an exact capacity.
   * @param arena Arena to allocate the backing block from.
   * @param exact_capacity Element count to allocate; appending never grows it.
   */
  ArenaVector(Arena &arena, size_t exact_capacity)
      : elements(nullptr), element_count(0), element_capacity(0) {
    bind(arena, exact_capacity);
  }

  /**
   * @brief Binds the vector to an arena, allocating its backing block.
   * @param arena Arena to allocate from.
   * @param min_capacity Minimum element count to reserve.
   * @details If already bound with at least that capacity, resets size for reuse
   * and keeps the larger prior capacity, so capacity() and the push_back
   * overflow guard report the block actually held rather than this request.
   * A grow against the same arena/generation reallocates a fresh block and
   * abandons the old one until the next reset; a stale binding (arena reset or
   * a different arena) trips a debug-only contract assert.
   */
  void bind(Arena &arena, size_t min_capacity) {
    static_assert(
        std::is_trivially_destructible_v<T> || is_arena_inplace_fn<T>::value,
        "ArenaVector never runs element destructors, so T must own no "
        "state outside the arena buffer: store a trivially-destructible "
        "type or a sanctioned Fn<> (no std::function/std::string).");
#ifndef NDEBUG
    // Rebinding a still-bound vector after its source arena was reset, or to a
    // different arena, is a contract violation (the old block is already dead). A
    // same-arena/same-generation grow is not stale and reallocates below.
    assert((!bound || element_capacity == 0 ||
            (stamp.source_arena == &arena &&
             stamp.birth_generation == arena.get_generation())) &&
           "ArenaVector::bind() on a stale binding: clear the handle before "
           "resetting or changing its arena");
#endif
    // The generation assert above misses a rewind, which leaves the reuse path
    // below handing back bytes the arena has already reissued.
    check_alive();
    // Same arena, still live, and big enough → reuse the block in place.
    if (bound && element_capacity >= min_capacity) {
      element_count = 0;
#ifndef NDEBUG
      // Reuse dangles any span snapshotted before this point; bump so its
      // check_alive() trips (the arena generation alone won't).
      rebind_generation++;
#endif
      return;
    }
    // Otherwise (unbound, or a grow that abandons the old block) → allocate
    // fresh. A grow leaks the old block until the next arena reset/compaction; a
    // zero-capacity binding owns no block, so growing out of one leaks nothing.
    if (bound && element_capacity > 0)
      log_arena_vector_grow(element_capacity * sizeof(T), element_capacity,
                            min_capacity);
    if (min_capacity > 0) {
      elements = arena.allocate_n<T>(min_capacity);
    } else {
      elements = nullptr;
    }
    element_count = 0;
    element_capacity = min_capacity;
    bound = true;
#ifndef NDEBUG
    if (min_capacity > 0)
      stamp.record(arena);
    rebind_generation++;
#endif
  }

  /**
   * @brief Reports whether the vector is bound to an arena.
   * @return True iff currently bound to an arena.
   */
  bool is_bound() const { return bound; }

  /**
   * @brief Appends a copy of an element.
   * @param value Element to copy-construct at the end.
   */
  void push_back(const T &value) {
    check_alive();
    check_bound();
    HS_CHECK(
        element_count < element_capacity,
        "ArenaVector push_back exact capacity exceeded! count=%lu capacity=%lu",
        static_cast<unsigned long>(element_count),
        static_cast<unsigned long>(element_capacity));
    new (&elements[element_count]) T(value);
    element_count++;
  }

  /**
   * @brief Bulk-append from a contiguous source.
   * @param src Pointer to the first source element.
   * @param count Number of elements to copy.
   * @details T must be trivially copyable (memcpy'd). An empty append is a
   * no-op that skips memcpy to avoid null-pointer UB.
   */
  void append_bulk(const T *src, size_t count) {
    static_assert(
        std::is_trivially_copyable_v<T>,
        "append_bulk memcpy's the source; T must be trivially copyable");
    check_alive();
    check_bound();
    // Subtractive, wrap-proof form: `element_count + count` could wrap for a
    // colossal count.
    HS_CHECK(
        count <= element_capacity - element_count,
        "ArenaVector bulk append exceeds capacity! count=%lu append=%lu capacity=%lu",
        static_cast<unsigned long>(element_count),
        static_cast<unsigned long>(count),
        static_cast<unsigned long>(element_capacity));
    // Skip memcpy on an empty append: a null src with count 0 is formal UB.
    if (count == 0)
      return;
    memcpy(static_cast<void *>(elements + element_count), src,
           count * sizeof(T));
    element_count += count;
  }

  /**
   * @brief Constructs an element in place at the end.
   * @tparam Args Constructor argument types forwarded to T.
   * @param args Arguments forwarded to T's constructor.
   * @return Reference to the newly constructed element.
   */
  template <typename... Args> T &emplace_back(Args &&...args) {
    check_alive();
    check_bound();
    HS_CHECK(
        element_count < element_capacity,
        "ArenaVector emplace_back exact capacity exceeded! count=%lu capacity=%lu",
        static_cast<unsigned long>(element_count),
        static_cast<unsigned long>(element_capacity));
    T *ptr = new (&elements[element_count]) T(std::forward<Args>(args)...);
    element_count++;
    return *ptr;
  }

  /**
   * @brief Element access by index.
   * @param i Index in [0, size()).
   * @return Mutable reference to the element at index i.
   */
  T &operator[](size_t i) {
    check_alive();
    check_bound();
    assert(i < element_count);
    return elements[i];
  }
  /**
   * @brief Element access by index (const).
   * @param i Index in [0, size()).
   * @return Const reference to the element at index i.
   */
  const T &operator[](size_t i) const {
    check_alive();
    check_bound();
    assert(i < element_count);
    return elements[i];
  }

  /**
   * @brief Returns the number of stored elements.
   * @return Current element count.
   */
  size_t size() const { return element_count; }
  /**
   * @brief Returns the maximum element count.
   * @return Capacity in elements.
   */
  size_t capacity() const { return element_capacity; }
  /**
   * @brief Reports whether the vector is empty.
   * @return True iff size() == 0.
   */
  bool is_empty() const { return element_count == 0; }
  /**
   * @brief Reports whether the vector is empty.
   * @return True iff size() == 0.
   * @details Container-requirement spelling of is_empty().
   */
  bool empty() const { return element_count == 0; }

  /**
   * @brief Accesses the last element.
   * @return Mutable reference to the element at index size() - 1.
   */
  T &back() {
    check_alive();
    check_bound();
    assert(element_count > 0);
    return elements[element_count - 1];
  }
  /**
   * @brief Accesses the last element (const).
   * @return Const reference to the element at index size() - 1.
   */
  const T &back() const {
    check_alive();
    check_bound();
    assert(element_count > 0);
    return elements[element_count - 1];
  }

  /**
   * @brief Resets the vector to empty without destroying elements.
   * @details No check_bound() here on purpose: clear() only resets
   * element_count (it neither frees nor touches elements), so it is a defined
   * no-op on an unbound vector. MeshState::clear() relies on this to reset
   * members that may never have been bound.
   */
  void clear() {
    check_alive();
#ifndef NDEBUG
    rebind_generation++;
#endif
    element_count = 0;
  }

  /**
   * @brief Returns a pointer to the backing storage.
   * @return Mutable pointer to the first element, or nullptr if unbound or
   * moved-from.
   * @details No check_bound() on data()/begin()/end() on purpose: an unbound (or
   * moved-from) vector is elements==nullptr with element_count==0, a
   * well-defined EMPTY range that callers rely on (std::span(vertices.data(),
   * size()), std::sort(data(), data()+count) on a size-0 vector). The
   * use-after-free guard (check_alive) still applies.
   */
  T *data() {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a pointer to the backing storage (const).
   * @return Const pointer to the first element, or nullptr if unbound or
   * moved-from.
   */
  const T *data() const {
    check_alive();
    return elements;
  }

  /**
   * @brief Returns an iterator to the first element.
   * @return Mutable pointer to the first element.
   */
  T *begin() {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns an iterator past the last element.
   * @return Mutable pointer one past the last element.
   */
  T *end() {
    check_alive();
    return elements ? elements + element_count : nullptr;
  }
  /**
   * @brief Returns a const iterator to the first element.
   * @return Const pointer to the first element.
   */
  const T *begin() const {
    check_alive();
    return elements;
  }
  /**
   * @brief Returns a const iterator past the last element.
   * @return Const pointer one past the last element.
   */
  const T *end() const {
    check_alive();
    return elements ? elements + element_count : nullptr;
  }
};
