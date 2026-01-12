/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cassert>
#include <array>
#include <memory>
#include <FastLED.h>
#include <variant>
#include "3dmath.h"
#include "FastNoiseLite.h"
#include "static_circular_buffer.h"

class Canvas;
class Orientation;

typedef CRGB Pixel;
typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;

#include "color.h" // Must be included before concepts using Color4

/**
 * @brief Concept for a function that determines color based on position and/or time.
 * Signature: Color4 f(const Vector& v, float t)
 */
template<typename F>
concept ColorFn = requires(F f, const Vector & v, float t) {
  { f(v, t) } -> std::convertible_to<Color4>;
};



/**
 * @brief Concept for a function that generates a trail color.
 * Signature: Color4 f(float x, float y, float t)
 */
template<typename F>
concept TrailFn = requires(F f, float x, float y, float t) {
  { f(x, y, t) } -> std::convertible_to<Color4>;
} || requires(F f, const Vector& v, float t) {
  { f(v, t) } -> std::convertible_to<Color4>;
};

/**
 * @brief Concept for a function that plots a path or curve.
 * Signature: Vector f(float t)
 */
template<typename F>
concept PlotFn = requires(F f, float t) {
  { f(t) } -> std::convertible_to<Vector>;
};

/**
 * @brief Concept for a scalar function (e.g. easing or wave).
 * Signature: float f(float t)
 */
template<typename F>
concept ScalarFn = requires(F f, float t) {
  { f(t) } -> std::convertible_to<float>;
};

/**
 * @brief Concept for a function that transforms a vector.
 * Signature: Vector f(const Vector& v)
 */
template<typename F>
concept TransformFn = requires(F f, Vector v) {
  { f(v) } -> std::convertible_to<Vector>;
};

/**
 * @brief Concept for a function that draws a sprite to the canvas.
 * Signature: void f(Canvas& canvas, float opacity)
 */
template<typename F>
concept SpriteFn = std::invocable<F, Canvas&, float>;

/**
 * @brief Concept for a timer callback function.
 * Signature: void f(Canvas& canvas)
 */
template<typename F>
concept TimerFn = std::invocable<F, Canvas&>;

/**
 * @brief Concept for a function called during tweening.
 * @details Can be called with EITHER (const Quaternion&, float) OR (const Orientation&, float).
 * This allows the same concept to support tweening an Orientation (yielding Quaternions)
 * or an OrientationTrail (yielding Orientations).
 */
template<typename F>
concept TweenFn = std::invocable<F, const Quaternion&, float> ||
std::invocable<F, const Orientation&, float>;

/**
 * @brief Concept for any object that maintains a history or sequence accessible by index.
 * Matches Orientation (get -> Quaternion) and OrientationTrail (get -> Orientation).
 */
template<typename T>
concept Tweenable = requires(const T & t, size_t i) {
  { t.length() } -> std::convertible_to<size_t>;
  { t.get(i) }; // Return type is deduced (Quaternion or Orientation)
};

#include "geometry.h"
#include "filter.h"
#include "draw.h"
#include "animation.h"
#include "FieldSampler.h"
