/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include "render/filter/pipeline.h"
#include "math/geometry.h"
#include "color/color.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"

/**
 * @file world_trails.h
 * @brief Filter::World::Trails: an arena-backed ring of quantized world-space
 * trail samples, re-emitted and aged once per frame.
 */

namespace Filter {
namespace World {

/**
 * @brief Manages 3D world-space trails.
 */
template <int Capacity> class Trails : public Is3DWithHistory {
public:
  static_assert(Capacity > 0, "World::Trails capacity must be positive");
  static constexpr bool emits_nonunit_world = true;
  static constexpr bool reads_outside_band = false;

  /** @brief One quantized trail sample: unit vector plus remaining lifetime. */
  struct Item {
    int16_t x, y, z; /**< Quantized unit vector components (6 bytes). */
    uint8_t ttl;     /**< Remaining lifetime in frames (1 byte). */
    uint8_t pad;     /**< Padding for 8-byte alignment (1 byte). */
  };
  static_assert(sizeof(Item) == 8, "World::Trails::Item must be 8 bytes");

  /**
   * @brief Constructs a world trail buffer with the given fade lifetime.
   * @param lifetime Per-frame fade divisor in frames; must be in [1, 255].
   * @details Upper bound is structural: ttl is a uint8_t, so lifetime > 255
   * would wrap the trail length.
   */
  Trails(int lifetime) : lifetime(lifetime) {
    HS_CHECK(lifetime > 0 && lifetime <= 255,
             "World::Trails: lifetime %d outside [1, 255]", lifetime);
  }

  /**
   * @brief Retunes the trail length at runtime (e.g. from a "Trail Len" slider).
   * @param new_lifetime New fade divisor in frames; must be in [1, 255].
   * @details Same bounds as the constructor; buffered points keep their ttl and
   * age out under the new length within a few frames.
   */
  void set_lifetime(int new_lifetime) {
    HS_CHECK(new_lifetime > 0 && new_lifetime <= 255,
             "World::Trails: lifetime %d outside [1, 255]", new_lifetime);
    lifetime = new_lifetime;
  }

  /**
   * @brief Allocates ring-buffer storage from the persistent arena.
   * @param arena Persistent arena supplying Capacity Item slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet).
   */
  void init_storage(Arena &arena) {
    items = arena.allocate_n<Item>(Capacity);
    head = tail = count = 0;
#ifndef NDEBUG
    stamp.record(arena);
#endif
  }

  /**
   * @brief Forwards the current sample and seeds a fading trail point.
   * @param v World-space point on the unit sphere.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames), non-negative; ttl = lifetime - age,
   * seeded only if positive.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged and NOT gated: a
   * transparent sample still consumes a ring slot. Screen::Trails deliberately
   * differs, dropping samples at its own MIN_TRAIL_ALPHA.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   * @details A point seeded here is still live for this frame's flush(), which
   * re-emits it at t = 0. Fresh samples therefore composite twice per frame:
   * once at the caller's color/alpha, once at the trailFn's.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    assert(age >= 0.0f &&
           "World::Trails: a negative age narrows ttl past the uint8_t range");
    pass(v, color, age, alpha);

    // round, not truncate (ttl is an integer byte)
    int ttl = lifetime - static_cast<int>(age + 0.5f);
    if (ttl > 0 && items) {
      check_storage_alive();
      push_back(encode(v, static_cast<uint8_t>(ttl)));
    }
  }

  /**
   * @brief Re-emits each buffered point colored by @p trailFn, then ages every
   * point one frame and culls the dead.
   * @param trailFn Callback producing trail color/alpha from (point, t).
   * @param alpha Global blend alpha in [0, 1].
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   * @details Emits before aging, so a point still renders on the frame its ttl
   * reaches 1 rather than being culled unseen.
   */
  template <typename PassFnT>
  void flush(const WorldTrailFn &trailFn, float alpha, PassFnT &&pass) {
    HS_CHECK(items, "World::Trails needs init_storage() from effect init()");
    check_storage_alive();
    for (size_t i = 0; i < count; ++i) {
      const auto &item = at(i);
      Vector v = decode(item);
      float t = hs::clamp(
          1.0f - (static_cast<float>(item.ttl) / static_cast<float>(lifetime)),
          0.0f, 1.0f);
      Color4 c = trailFn(v, t);

      if (c.alpha > MIN_TRAIL_ALPHA) {
        int age = lifetime - static_cast<int>(item.ttl);
        if (age < 0)
          age = 0;
        pass(v, c.color, static_cast<float>(age), c.alpha * alpha);
      }
    }

    for (size_t i = 0; i < count;) {
      Item &item = at(i);
      if (item.ttl > 0)
        item.ttl--;
      if (item.ttl == 0) {
        // swap-remove the logical-last live item (index count-1) into the dead
        // slot; only tail retreats, head stays put.
        item = at(count - 1);
        tail = (tail + Capacity - 1) % Capacity;
        count--;
      } else {
        ++i;
      }
    }
  }

  /**
   * @brief Returns the number of live trail points currently buffered.
   * @return Count of buffered Item entries.
   */
  size_t size() const { return count; }

private:
  /**
   * @brief Trail alpha below which flush() emits nothing for a buffered point.
   * @details Gates emission only; plot() seeds every sample whatever its
   * alpha.
   */
  static constexpr float MIN_TRAIL_ALPHA = 0.001f;

  Item *items = nullptr; /**< Ring-buffer storage (arena-owned). */
  size_t head = 0, tail = 0,
         count = 0; /**< Ring-buffer head, tail, and live count. */
  int lifetime;     /**< Per-frame fade divisor in frames. */
#ifndef NDEBUG
  ArenaBlockStamp stamp; /**< Arena state when items was allocated. */
#endif

  /**
   * @brief Debug-only use-after-free check on the arena-owned ring buffer.
   * @details A compaction that resets or rewinds the persistent arena without a
   * fresh init_storage() leaves items dangling; every plot()/flush() then reads
   * and writes Capacity Items through it.
   */
  void check_storage_alive() const {
#ifndef NDEBUG
    constexpr size_t BYTES = Capacity * sizeof(Item);
    assert(stamp.block_alive(items, BYTES) && "World::Trails use-after-free!");
#endif
  }

  static constexpr float Q =
      32767.0f; /**< Quantization scale for unit-vector components. */
  /**
   * @brief Encodes a unit vector and ttl into a quantized Item.
   * @param v World-space point; each component is clamped to [-1, 1] before scaling.
   * @param ttl Remaining lifetime in frames.
   * @return Packed Item with quantized coordinates.
   */
  static Item encode(const Vector &v, uint8_t ttl) {
    // clamp before scaling: an unclamped component past 1 overflows int16
    return {static_cast<int16_t>(hs::clamp(v.x, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.y, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.z, -1.0f, 1.0f) * Q), ttl, 0};
  }
  /**
   * @brief Decodes a quantized Item back into a near-unit vector.
   * @param item Packed trail sample.
   * @return Reconstructed world-space point.
   * @note Only *near* unit length (int16 quantization), and not renormalized;
   * a World::Trails must not precede a unit-assuming World filter.
   */
  static Vector decode(const Item &item) {
    constexpr float INV_Q = 1.0f / Q;
    return Vector(item.x * INV_Q, item.y * INV_Q, item.z * INV_Q);
  }

  /**
   * @brief Returns the i-th logical live item.
   * @param i Index into the live range [0, count).
   * @return Mutable reference to the buffered Item.
   */
  Item &at(size_t i) { return items[physical_index(i)]; }
  /**
   * @brief Returns the i-th logical live item.
   * @param i Index into the live range [0, count).
   * @return Const reference to the buffered Item.
   */
  const Item &at(size_t i) const { return items[physical_index(i)]; }

  /** @brief Maps a live logical index onto its physical ring slot. */
  size_t physical_index(size_t i) const {
    const size_t index = head + i;
    return index >= Capacity ? index - Capacity : index;
  }

  /**
   * @brief Appends an item, evicting a live item of arbitrary age at capacity.
   * @param item Encoded trail sample to push.
   */
  void push_back(const Item &item) {
    if (count == Capacity) {
      pop_front();
    }
    items[tail] = item;
    tail = (tail + 1) % Capacity;
    count++;
  }

  /** @brief Drops the logical head, whose age is arbitrary after compaction. */
  void pop_front() {
    head = (head + 1) % Capacity;
    count--;
  }
};

} // namespace World
} // namespace Filter
