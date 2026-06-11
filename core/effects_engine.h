/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Umbrella header for the effects engine: pulls in the full public API
// (geometry, color, animation, plotting, palettes, presets, ...) that every
// effect in effects/ includes via this one file.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <string_view>

#include "3dmath.h"
#include "FastNoiseLite.h"
#include "canvas.h"
#include "static_circular_buffer.h"

#include "geometry.h" // Provides Fragment and defines Orientation<CAP = 4>
#include "reaction_graph.h"
#include "concepts.h" // Concepts needs fragment
#include "color.h"
#include "animation.h"
#include "transformers.h"

#include "filter.h"
#include "platform.h"
#include "plot.h"
#include "scan.h"
#include "mesh.h"
#include "hankin.h"
#include "conway.h"
#include "solids.h"
#include "palettes.h"

#include "presets.h"
#include "waves.h"
#include "styles.h"
