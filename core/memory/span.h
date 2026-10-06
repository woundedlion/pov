/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/memory.h.

// ============================================================================
// Non-Owning Span (Explicit Borrow)
// ============================================================================

/**
 * @brief A read-only, non-owning view into arena-allocated data.
 * @tparam T Element type viewed by the span.
 * @details Makes owned (ArenaVector) vs borrowed data visible at the type level.
 *
 * LIFETIME CONTRACT: a span snapshots its source vector's elements pointer at
 * construction. Debug checks fault on an arena RESET, a source-vector RE-GROW or CLEAR,
 * a rewind below the borrowed block, or reclamation and reissue of that block.
 * The arena generation tracks reset, the per-vector rebind counter tracks grow and clear,
 * and the block stamp tracks rewind and reissue. A MOVE of the source vector is not
 * tracked: the span keeps its snapshotted elements (runtime-safe) but its debug
 * stamps reference the moved-from husk, so re-take the span after growing or
 * moving its source. Outliving the source VECTOR OBJECT (not its arena block —
 * a stack-local vector going out of scope) is worse than untracked: in debug the
 * staleness check reads source_vec, so the span must not outlive it.
 */
template <typename T> class ArenaSpan {
  const T *elements;    /**< Snapshotted pointer to the borrowed data. */
  size_t element_count; /**< Number of viewed elements. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the block was allocated. */
  const ArenaVector<T> *source_vec =
      nullptr; /**< Source vector for re-grow check. */
  uint32_t source_rebind_generation =
      0; /**< Vector rebind counter at construction. */

  /**
   * @brief Debug-only stale-span check against arena and vector stamps.
   * @details Asserts on arena reset, source-vector re-grow, a rewind below the
   * borrowed block, or reclamation and reissue of that block.
   */
  void check_alive() const {
    const size_t bytes = element_count * sizeof(T);
    assert(!stamp.arena_reset() && "ArenaSpan use-after-free!");
    assert((!source_vec ||
            source_vec->rebind_generation == source_rebind_generation) &&
           "ArenaSpan source vector re-grown out from under span!");
    assert(!stamp.block_uncovered(elements, bytes) &&
           "ArenaSpan use-after-free (arena rewound below span)!");
    assert(!stamp.block_reissued(elements, bytes) &&
           "ArenaSpan use-after-free (borrowed block reclaimed by a rewind and "
           "reissued)!");
  }
#else
  /**
   * @brief No-op stale-span check in release builds.
   */
  void check_alive() const {}
#endif

public:
  // Read-only view: the mutable spellings alias the const ones.
  using value_type = T;
  using size_type = size_t;
  using difference_type = std::ptrdiff_t;
  using reference = const T &;
  using const_reference = const T &;
  using pointer = const T *;
  using const_pointer = const T *;
  using iterator = const T *;
  using const_iterator = const T *;

#ifndef NDEBUG
  /** @brief Source binding generation snapshotted by this debug view. */
  uint32_t debug_binding_generation() const { return source_rebind_generation; }
#endif

  /**
   * @brief Default-constructs an empty span.
   */
  ArenaSpan() : elements(nullptr), element_count(0) {}

  /**
   * @brief Constructs a span borrowing from an ArenaVector (explicit borrow).
   * @param source Vector to borrow data and (in debug) lifetime stamps from.
   * @details In debug builds the span inherits the vector's arena-generation
   * stamp, so accessing it after the arena is reset/compacted trips the same
   * use-after-free check the vector itself has.
   */
  explicit ArenaSpan(const ArenaVector<T> &source)
      : elements(source.data()), element_count(source.size())
#ifndef NDEBUG
        ,
        stamp(source.stamp), source_vec(&source),
        source_rebind_generation(source.rebind_generation)
#endif
  {
  }

  /**
   * @brief Deleted constructor from a temporary ArenaVector.
   * @details In debug builds source_vec would dangle when the temporary dies,
   * even while its arena storage remains live.
   */
  explicit ArenaSpan(const ArenaVector<T> &&) = delete;

  /**
   * @brief Copy and copy-assignment duplicate the borrow verbatim.
   * @details A span copy carries the same data pointer, size, and (in debug) the
   * source's lifetime stamps, so the copy trips the same staleness check as the
   * original. Only construction from a temporary ArenaVector (above) is forbidden.
   */
  ArenaSpan(const ArenaSpan &) = default;
  ArenaSpan &operator=(const ArenaSpan &) = default;

  /**
   * @brief Element access by index.
   * @param i Index in [0, size()).
   * @return Const reference to the element at index i.
   */
  const T &operator[](size_t i) const {
    check_alive();
    assert(i < element_count);
    return elements[i];
  }
  /**
   * @brief Returns the number of viewed elements.
   * @return Element count.
   */
  size_t size() const {
    check_alive();
    return element_count;
  }
  /**
   * @brief Reports whether the span is empty.
   * @return True iff size() == 0.
   */
  bool is_empty() const {
    check_alive();
    return element_count == 0;
  }
  /**
   * @brief Reports whether the span is empty.
   * @return True iff size() == 0.
   * @details Container-requirement spelling of is_empty().
   */
  bool empty() const {
    check_alive();
    return element_count == 0;
  }
  /**
   * @brief Returns a pointer to the borrowed storage.
   * @return Const pointer to the first element.
   */
  const T *data() const {
    check_alive();
    return elements;
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
