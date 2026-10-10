/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <cstdint>
#include "math/spherical_field.h"
#include "memory.h"
#include "render/filter/feedback_cap_plane.h"
#include "render/filter/feedback_style.h"

/**
 * @file feedback_warp_cache.h
 * @brief Filter::Pixel::FeedbackWarpCache: the feedback filter's persistent
 * coarse warp field, reused while its inputs are unchanged.
 */

namespace Filter {

namespace Pixel {

/**
 * @brief Arena-owned warp-field buffers with the key they were filled for.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details acquire() is the only way to mark the buffers valid: it stores the
 * key only after its fill returns. The cache also holds the projected origin
 * of every lattice sample of its layout, filled once by init_storage().
 */
template <int W, int H> class FeedbackWarpCache {
  using SphereField = hs::SphericalFieldLayout<W, H>;
  using CapOffset = typename FeedbackCapPlane<W, H>::CapOffset;

public:
  /// Field coordinates of a lattice origin.
  using Coordinates = typename SphereField::Coordinates;

  /** @brief Inputs the coarse warp field is a pure function of (stock
   *  transforms only); equal keys make the cached field reusable. */
  struct Key {
    ::Feedback::SpaceFn space_fn;        ///< Spatial warp transform.
    const Animation::NoiseParams *noise; ///< Bound noise generator, or null.
    uint32_t noise_config; ///< Generator configuration key (its seed).
    float amplitude;       ///< Noise amplitude.
    float frequency;       ///< Noise spatial frequency.
    float speed;           ///< Noise temporal speed.
    float scale;           ///< Noise spatial scale.
    float time;            ///< Noise sampling time; 0 when unbound.
    int field_y_begin;     ///< First field ring of the band.
    int field_y_end;       ///< Last field ring of the band, inclusive.
    bool operator==(const Key &) const = default;
  };

  /** @brief The warp field's per-cell offsets. */
  struct Buffers {
    int16_t *x_offsets;   /**< Column offsets, one per cell. */
    int16_t *y_offsets;   /**< Row offsets, one per cell. */
    CapOffset *cell_caps; /**< Cap-plane offsets, one per polar-ring cell. */
  };

  /**
   * @brief Allocates the cache from the persistent arena and fills its
   * lattice origins.
   * @tparam CELL_COUNT Cell count of the warp field.
   * @param arena Persistent arena.
   * @param field Spherical layout the cache serves.
   * @param cap_cells Cell count of its polar rings; 0 allocates no cap cells.
   * @param origin Maps a lattice position and its exact field coordinates to
   *        the coordinates offsets are measured from.
   * @details Call again after any arena reset. Leaves the cache invalid.
   */
  template <int CELL_COUNT, typename OriginFnT>
  HS_COLD_MEMBER void init_storage(Arena &arena, const SphereField &field,
                                   int cap_cells, OriginFnT &&origin) {
#ifndef NDEBUG
    HS_CHECK(
        !column_offsets ||
            !stamp.block_alive(column_offsets, CELL_COUNT * sizeof(int16_t)),
        "feedback filter: storage already initialized");
    stamped_cells = CELL_COUNT;
    stamped_cap_cells = cap_cells;
#endif
    column_offsets = arena.allocate_n<int16_t>(CELL_COUNT);
    row_offsets = arena.allocate_n<int16_t>(CELL_COUNT);
    origin_points = arena.allocate_n<Coordinates>(CELL_COUNT);
    cap_field =
        cap_cells > 0 ? arena.allocate_n<CapOffset>(cap_cells) : nullptr;
    valid = false;
#ifndef NDEBUG
    stamp.record(arena);
#endif
    hs::SphericalField<Coordinates, W, H> origins(origin_points, field);
    origins.populate(0, field.ring_count() - 1, origin);
  }

  /** @brief Whether init_storage() has run. */
  bool ready() const { return column_offsets != nullptr; }

  /** @brief Whether the buffers hold the field filled for @p key. */
  bool holds(const Key &key) const { return valid && key == stored_key; }

  /** @brief Projected origin of every lattice sample of the cache's layout,
   *  or nullptr before init_storage(). */
  const Coordinates *origins() const { return origin_points; }

  /**
   * @brief Returns buffers holding the field for @p key, filling them through
   * @p fill when the cache holds another key's.
   * @param key Inputs of the field the caller needs; nullptr when the field
   *        is not cacheable.
   * @param scratch Caller-owned buffers filled and returned when @p key is
   *        nullptr.
   * @param fill Called as fill(buffers); must write every cell the field's
   *        band covers.
   * @return The cache's buffers, or @p scratch when @p key is nullptr.
   * @details The cache is invalid while @p fill runs on its buffers; the key
   * is stored, and the cache marked valid, only after @p fill returns. A
   * nullptr @p key leaves the cache untouched.
   */
  template <typename FillFnT>
  __attribute__((always_inline)) Buffers acquire(const Key *key,
                                                 const Buffers &scratch,
                                                 FillFnT &&fill) {
    if (key && holds(*key))
      return {column_offsets, row_offsets, cap_field};
    const Buffers buffers =
        key ? Buffers{column_offsets, row_offsets, cap_field} : scratch;
    if (key)
      valid = false;
    fill(buffers);
    if (key) {
      stored_key = *key;
      valid = true;
    }
    return buffers;
  }

  /**
   * @brief Debug-only use-after-free check on the arena-owned buffers.
   */
  void check_storage_alive() const {
#ifndef NDEBUG
    HS_ASSERT_BLOCK_ALIVE(stamp, column_offsets,
                          stamped_cells * sizeof(int16_t),
                          "Pixel::Feedback warp cache");
    HS_ASSERT_BLOCK_ALIVE(stamp, row_offsets, stamped_cells * sizeof(int16_t),
                          "Pixel::Feedback warp cache");
    HS_ASSERT_BLOCK_ALIVE(stamp, origin_points,
                          stamped_cells * sizeof(Coordinates),
                          "Pixel::Feedback warp cache");
    if (cap_field) {
      HS_ASSERT_BLOCK_ALIVE(stamp, cap_field,
                            stamped_cap_cells * sizeof(CapOffset),
                            "Pixel::Feedback warp cache");
    }
#endif
  }

private:
  Key stored_key{};   /**< Key the buffers were last filled for. */
  bool valid = false; /**< True once a fill has completed. */
  int16_t *column_offsets = nullptr;    /**< Arena-owned column offsets. */
  int16_t *row_offsets = nullptr;       /**< Arena-owned row offsets. */
  Coordinates *origin_points = nullptr; /**< Arena-owned lattice origins. */
  CapOffset *cap_field = nullptr; /**< Arena-owned polar-ring cap offsets. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when the buffers were allocated. */
  int stamped_cells = 0; /**< Cells in each per-cell buffer. */
  int stamped_cap_cells = 0; /**< Cells in the cap-offset buffer. */
#endif
};

} // namespace Pixel

} // namespace Filter
