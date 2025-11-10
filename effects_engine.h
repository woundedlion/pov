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
typedef std::function<Pixel(const Vector&, float_t)> ColorFn;
typedef std::function<Pixel(float_t x, float_t y, float_t t)> TrailFn;
typedef std::function<float_t(float_t)> EasingFn;
typedef std::function<Vector(float_t)> PlotFn;
typedef std::function<float_t(float_t)> ShiftFn;
typedef std::function<void(Canvas&, float_t)> SpriteFn;
typedef std::function<void(Canvas&)> TimerFn;
typedef std::function<float_t(float_t)> MutateFn;
typedef std::function<float_t(float_t)> WaveFn;

#include "geometry.h"
#include "color.h"
#include "filter.h"
#include "draw.h"
#include "animation.h"