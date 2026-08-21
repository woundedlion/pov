/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/pipeline.h"
#include "color/color.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"

/**
 * @file screen_trails.h
 * @brief Filter::Screen::Trails: an arena-backed set of fading screen-space
 * trail points, re-emitted and aged once per frame.
 */

namespace Filter {
namespace Screen {

/**
 * @brief Manages 2D screen-space trails.
 */
template <int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  // Trail points are seeded from and re-emitted into the same band, so they
  // never sample a neighbor segment.
  static constexpr bool crosses_segments = false;
  static constexpr bool reads_outside_band = false;

  /**
   * @brief Constructs a screen trail buffer with the given fade lifetime.
   * @param lifetime Per-frame fade divisor in frames; must be positive.
   */
  Trails(int lifetime) : lifetime(lifetime) {
    HS_CHECK(lifetime > 0, "Screen::Trails: lifetime %d must be positive",
             lifetime);
  }

  /**
   * @brief Allocates the decay-pixel storage from the persistent arena.
   * @param arena Persistent arena supplying MAX_PIXELS DecayPixel slots.
   */
  void init_storage(Arena &arena) {
    points = arena.allocate_n<DecayPixel>(MAX_PIXELS);
    num_pixels = 0;
#ifndef NDEBUG
    stamp.record(arena);
#endif
  }

  /**
   * @brief Forwards the current sample and seeds a fading screen trail point.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1]; samples at or below MIN_TRAIL_ALPHA
   * seed no trail point but are still forwarded downstream. World::Trails
   * deliberately differs, seeding on every sample.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   * @details At MAX_PIXELS the slot at index 0 is evicted; decay()'s unordered
   * compaction leaves that point's age arbitrary.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(x, y, color, age, alpha);

    if (alpha <= MIN_TRAIL_ALPHA)
      return;

    float ttl = static_cast<float>(lifetime) - age;
    if (ttl > 0.0f && points) {
      check_storage_alive();
      if (num_pixels == MAX_PIXELS) {
        num_pixels--;
        points[0] = points[num_pixels];
      }
      points[num_pixels++] = {x, y, ttl};
    }
  }

  /**
   * @brief Re-emits each buffered trail point colored by @p trailFn.
   * @param trailFn Callback producing trail color/alpha from (x, y, t).
   * @param alpha Global blend alpha in [0, 1].
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   * @details The unused Canvas parameter satisfies the 2D flush signature; ages
   * all points one frame via decay() after emission.
   */
  template <typename PassFnT>
  void flush(Canvas &, const ScreenTrailFn &trailFn, float alpha,
             PassFnT &&pass) {
    HS_CHECK(points, "Screen::Trails needs init_storage() from effect init()");
    check_storage_alive();
    for (int i = 0; i < num_pixels; ++i) {
      float t = hs::clamp(1.0f - (points[i].ttl / lifetime), 0.0f, 1.0f);
      Color4 color = trailFn(points[i].x, points[i].y, t);
      if (color.alpha > MIN_TRAIL_ALPHA) {
        float age = lifetime - points[i].ttl;
        if (age < 0.0f)
          age = 0.0f;
        pass(points[i].x, points[i].y, color.color, age, alpha * color.alpha);
      }
    }
    decay();
  }

  /**
   * @brief Ages every point one frame and swap-removes dead slots.
   * @details Unordered compaction: a dead slot is overwritten by the last live
   * point. ttl decrements by 1 per frame (whole-frame model), so a point
   * survives ceil(ttl) frames.
   */
  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--points[i].ttl <= 0.0f) {
        points[i] = points[--num_pixels];
        i--;
      }
    }
  }

private:
  /**
   * @brief Trail alpha below which a sample seeds nothing and a buffered point
   * emits nothing.
   * @details Same value as World::Trails::MIN_TRAIL_ALPHA, but gates seeding as
   * well as emission.
   */
  static constexpr float MIN_TRAIL_ALPHA = 0.001f;

  /** @brief One screen trail point: position plus remaining lifetime. */
  struct DecayPixel {
    float x, y, ttl; /**< Pixel position and remaining lifetime in frames. */
  };
  int lifetime;                 /**< Per-frame fade divisor in frames. */
  DecayPixel *points = nullptr; /**< Arena-owned array of live trail points. */
  int num_pixels = 0;           /**< Number of live points in points. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when points was allocated. */
#endif

  /**
   * @brief Debug-only use-after-free check on the arena-owned point array.
   * @details A compaction that resets or rewinds the persistent arena without a
   * fresh init_storage() leaves points dangling; every plot()/flush() then reads
   * and writes MAX_PIXELS DecayPixels through it.
   */
  void check_storage_alive() const {
#ifndef NDEBUG
    constexpr size_t BYTES = MAX_PIXELS * sizeof(DecayPixel);
    assert(!stamp.arena_reset() && "Screen::Trails use-after-free!");
    assert(!stamp.block_uncovered(points, BYTES) &&
           "Screen::Trails use-after-free (arena rewound below block)!");
    assert(!stamp.block_reissued(points, BYTES) &&
           "Screen::Trails use-after-free (block reclaimed by a rewind and "
           "reissued)!");
#endif
  }
};

} // namespace Screen
} // namespace Filter
