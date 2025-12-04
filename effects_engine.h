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
#include <array>
#include "3dmath.h"
#include "FastNoiseLite.h"
#include "static_circular_buffer.h"

typedef CRGB Pixel;
typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;

template<typename F>
concept ColorFn = std::invocable<F, const Vector&, float>
&& std::convertible_to<std::invoke_result_t<F, const Vector&, float>, Pixel>;

template<typename F>
concept TrailFn = std::invocable<F, float, float, float>
&& std::convertible_to<std::invoke_result_t<F, float, float, float>, Pixel>;

template<typename F>
concept PlotFn = std::invocable<F, float>
&& std::convertible_to<std::invoke_result_t<F, float>, Vector>;

template<typename F>
concept ScalarFn = std::invocable<F, float>
&& std::convertible_to<std::invoke_result_t<F, float>, float>;

template<typename F>
concept SpriteFn = std::invocable<F, Canvas&, float>;

template<typename F>
concept TimerFn = std::invocable<F, Canvas&>;

template<typename F>
concept TweenFn = std::invocable<F, const Quaternion&, float>;

#include "geometry.h"
#include "color.h"
#include "filter.h"
#include "draw.h"
#include "animation.h"