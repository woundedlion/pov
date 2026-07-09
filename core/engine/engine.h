/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Engine API umbrella: pulls in the full public API (geometry, color,
// animation, plotting, palettes, presets, ...) that every effect in effects/
// includes via this one file. Distinct from core/engine/effects.h, which is the effect
// *roster* (it pulls in all 27 effect headers); include this from an effect,
// never that.

// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro — otherwise assert-stripping would depend
// on a prior TU having pulled in platform.h (see canvas.h's identical note).
#include "engine/platform.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <memory>
#include <string_view>

#include "math/3dmath.h"
#include "vendor/FastNoiseLite.h"
#include "render/canvas.h"
#include "engine/static_circular_buffer.h"

#include "math/geometry.h"
#include "render/shading.h" // Fragment + mesh-topology shading
#include "engine/reaction_graph.h"
#include "engine/concepts.h" // Concepts needs fragment
#include "color/color.h"
#include "animation/animation.h"
#include "engine/transformers.h"

#include "render/filter.h"
#include "render/plot.h"
#include "render/scan.h"
#include "mesh/mesh.h"
#include "mesh/hankin.h"
#include "mesh/conway.h"
#include "mesh/solids.h"
#include "color/palettes.h"

#include "engine/presets.h"
#include "math/waves.h"
#include "engine/styles.h"
