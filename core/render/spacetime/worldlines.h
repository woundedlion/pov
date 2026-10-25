/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file worldlines.h
 * @brief Per-vertex frame captures carried from one frame's end to the next
 * frame's start, and the interpolant each vertex's path takes between them,
 * chosen once per frame from an analytic error bound.
 */

#include "core/render/spacetime/motion.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace Spacetime {

namespace detail {
/**
 * @brief Packs a capture for storage.
 * @param v Position.
 * @return @p v unchanged.
 */
HS_HOT_INLINE inline math::Vector pack_capture(const math::Vector &v,
                                               math::Vector *) {
  return v;
}

/**
 * @brief Packs a unit-range capture as snorm16.
 * @param v Position with components in [-1, 1].
 * @return The snorm16 triple.
 */
HS_HOT_INLINE inline math::Snorm3 pack_capture(const math::Vector &v,
                                               math::Snorm3 *) {
  return math::Snorm3::encode(v);
}

/**
 * @brief Unpacks a stored capture.
 * @param v Stored position.
 * @return @p v.
 */
HS_HOT_INLINE inline math::Vector unpack_capture(const math::Vector &v) {
  return v;
}

/**
 * @brief Unpacks a snorm16 capture.
 * @param v Stored position.
 * @return The decoded position, within 2^-16 per component.
 */
HS_HOT_INLINE inline math::Vector unpack_capture(const math::Snorm3 &v) {
  return v.decode();
}

} // namespace detail

/**
 * @brief Per-vertex frame captures and interpolants of one geometry.
 * @tparam Stored Capture storage: math::Snorm3 for unit-range positions, or
 *         math::Vector to keep captures exact.
 * @details Storage is arena-backed and sized once for the geometry's most
 * vertices. Each frame, begin_frame() turns the previous frame's end captures
 * into this frame's start captures, so a vertex costs one new evaluation, its
 * end; choose() picks its interpolant from an analytic bound, and a vertex
 * choosing QUADRATIC adds one more, at the frame's middle.
 */
template <typename Stored = math::Snorm3> class Worldlines {
public:
  /**
   * @brief Arena bytes bind() takes.
   * @param capacity Most vertices.
   * @return The storage's size, alignment slack included.
   */
  static constexpr size_t storage_bytes(size_t capacity) {
    return 3 * (capacity * sizeof(Stored) + alignof(Stored)) + capacity;
  }

  /**
   * @brief Allocates the captures.
   * @param arena Arena that outlives every frame the motion spans.
   * @param capacity Most vertices.
   */
  HS_COLD_MEMBER void bind(Arena &arena, size_t capacity) {
    slots = capacity;
    captures[0] = arena.allocate_n<Stored>(capacity);
    captures[1] = arena.allocate_n<Stored>(capacity);
    middles = arena.allocate_n<Stored>(capacity);
    kinds = arena.allocate_n<Interpolant>(capacity);
    count = 0;
    carried = false;
  }

  /**
   * @brief Starts a frame.
   * @param vertices Vertices this frame, at most the bound capacity.
   * @param carry Whether vertex `i` is the same vertex as in the previous
   *        frame, so its end capture is this frame's start; otherwise each
   *        vertex starts where its end capture puts it, without motion.
   */
  HS_HOT_INLINE void begin_frame(size_t vertices, bool carry) {
    HS_CHECK(vertices <= slots, "worldline frame exceeds bound capacity");
    count = vertices;
    carried = carry;
    current ^= 1;
  }

  /**
   * @brief Records a vertex's position at the frame's end.
   * @param i Vertex index.
   * @param p Position at t = 1.
   */
  HS_HOT_INLINE void capture_end(size_t i, const math::Vector &p) {
    const Stored packed =
        detail::pack_capture(p, static_cast<Stored *>(nullptr));
    captures[current][i] = packed;
    if (!carried)
      captures[current ^ 1][i] = packed;
    kinds[i] = Interpolant::CHORD;
  }

  /**
   * @brief Holds a vertex at its end capture for this frame, as after a
   *        discrete change.
   * @param i Vertex index, its end already captured.
   */
  HS_HOT_INLINE void hold(size_t i) {
    captures[current ^ 1][i] = captures[current][i];
    kinds[i] = Interpolant::CHORD;
  }

  /**
   * @brief Records a vertex's position at the frame's middle.
   * @param i Vertex index.
   * @param p Position at t = 1/2.
   */
  HS_HOT_INLINE void capture_middle(size_t i, const math::Vector &p) {
    middles[i] = detail::pack_capture(p, static_cast<Stored *>(nullptr));
  }

  /**
   * @brief Picks and records a vertex's interpolant for this frame.
   * @param i Vertex index, its end captured.
   * @param bound Derivative bounds of its path over the frame.
   * @param tolerance Largest acceptable distance from the path.
   * @param on_arc Whether the path is a circle path() can name.
   * @return The interpolant; QUADRATIC asks for capture_middle().
   */
  HS_HOT_INLINE Interpolant choose(size_t i, const PathBound &bound,
                                   float tolerance, bool on_arc) {
    kinds[i] = choose_interpolant(bound, tolerance, on_arc);
    return kinds[i];
  }

  /** @brief Vertices this frame.
   * @return The count begin_frame() set. */
  size_t size() const { return count; }

  /** @brief Most vertices the storage holds.
   * @return The capacity bind() took. */
  size_t capacity() const { return slots; }

  /**
   * @brief A vertex's position at the frame's start.
   * @param i Vertex index.
   * @return The previous frame's end capture, or this frame's without carry.
   */
  HS_HOT_INLINE math::Vector start(size_t i) const {
    return detail::unpack_capture(captures[current ^ 1][i]);
  }

  /**
   * @brief A vertex's position at the frame's end.
   * @param i Vertex index.
   * @return Its end capture.
   */
  HS_HOT_INLINE math::Vector end(size_t i) const {
    return detail::unpack_capture(captures[current][i]);
  }

  /**
   * @brief A vertex's position at the frame's middle.
   * @param i Vertex index, its middle captured this frame.
   * @return Its middle capture.
   */
  HS_HOT_INLINE math::Vector middle(size_t i) const {
    return detail::unpack_capture(middles[i]);
  }

  /**
   * @brief A vertex's interpolant this frame.
   * @param i Vertex index.
   * @return What choose() recorded, or CHORD.
   */
  HS_HOT_INLINE Interpolant interpolant(size_t i) const { return kinds[i]; }

  /**
   * @brief A vertex's path this frame.
   * @param i Vertex index.
   * @param axis Unit axis of an ARC path, or null for the great circle
   *        through its captures.
   * @return The path its interpolant describes.
   */
  HS_HOT_INLINE Worldline worldline(size_t i,
                                    const math::Vector *axis = nullptr) const {
    switch (kinds[i]) {
    case Interpolant::CHORD:
      return Worldline::chord(start(i), end(i));
    case Interpolant::QUADRATIC:
      return Worldline::quadratic(start(i), middle(i), end(i));
    case Interpolant::ARC:
      break;
    }
    return Worldline::arc(start(i), end(i), axis);
  }

private:
  Stored *captures[2] = {nullptr, nullptr}; ///< End captures, by frame parity.
  Stored *middles = nullptr;                ///< Middle captures.
  Interpolant *kinds = nullptr;             ///< Interpolant per vertex.
  size_t slots = 0;                         ///< Capacity.
  size_t count = 0;                         ///< Vertices this frame.
  int current = 0;      ///< Index of this frame's end captures.
  bool carried = false; ///< Start captures carry from the previous frame.
};

} // namespace Spacetime
