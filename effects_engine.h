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
#include "platform.h"
#include <variant>
#include "3dmath.h"
#include "FastNoiseLite.h"
#include "static_circular_buffer.h"
#include "canvas.h"

template <int W, int HISTORY = W + 1> class Orientation;

#include "geometry.h" // Provides Fragment, ShaderResult, Vector
#include "color.h" // Must be included before concepts using Color4
#include "solids.h"
#include "concepts.h"
#include "filter.h"
#include "scan.h"
#include "plot.h"
#include "animation.h"
