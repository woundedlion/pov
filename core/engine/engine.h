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
#include "core/platform/platform.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <memory>
#include <string_view>

#include "core/math/3dmath.h"
#include "core/vendor/FastNoiseLite.h"
#include "core/render/canvas.h"
#include "core/containers/static_circular_buffer.h"

#include "core/math/geometry.h"
#include "core/render/shading.h" // Fragment + mesh-topology shading
#include "core/spatial/reaction_graph.h"
#include "core/engine/concepts.h"
#include "core/memory.h"
#include "core/color/color.h"
#include "core/color/palette_cycler.h"
#include "core/animation/animation.h"
#include "core/animation/transformer.h"
#include "core/control/params.h"

#include "core/render/filter.h"
#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/mesh/mesh.h"
#include "core/mesh/hankin.h"
#include "core/mesh/conway.h"
#include "core/mesh/solids.h"
#include "core/color/palettes.h"

#include "core/control/presets.h"
#include "core/control/choreography.h"
#include "core/math/waves.h"
#include "core/control/registry.h"
