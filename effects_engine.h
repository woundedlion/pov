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

typedef std::function<Pixel(const Vector&, float)> ColorFn;
typedef std::function<Pixel(float x, float y, float t)> TrailFn;

typedef std::function<Vector(float)> PlotFn;

typedef std::function<float(float)> EasingFn;
typedef std::function<float(float)> ShiftFn;
typedef std::function<float(float)> MutateFn;
typedef std::function<float(float)> WaveFn;

typedef std::function<void(Canvas&, float)> SpriteFn;
typedef std::function<void(Canvas&)> TimerFn;

#include "geometry.h"
#include "color.h"
#include "filter.h"
#include "draw.h"
#include "animation.h"