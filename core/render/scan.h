/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file scan.h
 * @brief The scanline rasterizer: Scan::rasterize, the SDF-backed draw
 *        primitives, the mesh path, the full-screen shader and the raymarcher:
 *        umbrella over render/scan/.
 */

#include "render/scan/raster.h"
#include "render/scan/shapes.h"
#include "render/scan/mesh.h"
#include "render/scan/shader.h"
#include "render/scan/volume.h"
