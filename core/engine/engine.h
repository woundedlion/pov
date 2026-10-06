/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file engine.h
 * @brief Engine API umbrella header for handwritten effects.
 */

// platform.h supplies the NDEBUG fallback and must precede <cassert>.
#include "platform/platform.h"

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
#include "containers/static_circular_buffer.h"

#include "math/geometry.h"
#include "render/shading.h" // Fragment + mesh-topology shading
#include "spatial/reaction_graph.h"
#include "engine/concepts.h"
#include "memory.h"
#include "color/color.h"
#include "color/palette_cycler.h"
#include "animation/animation.h"
#include "animation/transformer.h"
#include "control/params.h"

#include "render/filter.h"
#include "render/plot.h"
#include "render/scan.h"
#include "mesh/mesh.h"
#include "mesh/hankin.h"
#include "mesh/conway.h"
#include "mesh/solids.h"
#include "color/palettes.h"

#include "control/presets.h"
#include "control/choreography.h"
#include "math/waves.h"
#include "control/registry.h"
