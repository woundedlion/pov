/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file filter.h
 * @brief The filter pipeline: Pipeline composition plus the World, Screen and
 * Pixel filter families.
 * @details The stage roster is a library surface, deliberately wider than the
 * set of stages the shipping effects instantiate; a stage with no current user
 * is composable inventory, not dead code.
 */

#include "render/filter/pipeline.h"
#include "render/filter/world_orient.h"
#include "render/filter/world_orient_slice.h"
#include "render/filter/world_hole.h"
#include "render/filter/world_replicate.h"
#include "render/filter/world_vertex_replicate.h"
#include "render/filter/world_mobius.h"
#include "render/filter/world_trails.h"
#include "render/filter/screen_splat.h"
#include "render/filter/screen_anti_alias.h"
#include "render/filter/screen_direct_aa_sink.h"
#include "render/filter/screen_trails.h"
#include "render/filter/screen_blur.h"
#include "render/filter/pixel_chromatic_shift.h"
