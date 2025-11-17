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

typedef std::function<Pixel(const Vector&, double)> ColorFn;
typedef std::function<Pixel(double x, double y, double t)> TrailFn;

typedef std::function<Vector(double)> PlotFn;

typedef std::function<double(double)> EasingFn;
typedef std::function<double(double)> ShiftFn;
typedef std::function<double(double)> MutateFn;
typedef std::function<double(double)> WaveFn;

typedef std::function<void(Canvas&, double)> SpriteFn;
typedef std::function<void(Canvas&)> TimerFn;

#include "geometry.h"
#include "color.h"
#include "filter.h"
#include "draw.h"
#include "animation.h"