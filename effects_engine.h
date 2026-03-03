/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <list>
#include <map> // for map
#include <memory>
#include <string>  // for string
#include <variant> // for variant
#include <vector>

#include "3dmath.h"
#include "FastNoiseLite.h"
#include "canvas.h"
#include "static_circular_buffer.h"

template <int W> class Orientation;

#include "geometry.h" // Provides Fragment, ShaderResult, Vector
#include "concepts.h" // Concepts needs fragment
#include "color.h"
#include "animation.h"
#include "transformers.h"
#include "generators.h"

#include "filter.h"
#include "platform.h"
#include "plot.h"
#include "scan.h"
#include "solids.h"
#include "palettes.h"

#include "presets.h"
#include "mesh.h"
#include "waves.h"
