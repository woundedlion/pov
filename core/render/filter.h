/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file filter.h
 * @brief The filter pipeline: Pipeline composition plus the World, Screen and
 * Pixel filter families.
 */

#include "core/render/filter/pipeline.h"
#include "core/render/filter/splat.h"
#include "core/render/filter/world_orient.h"
#include "core/render/filter/world_orient_slice.h"
#include "core/render/filter/world_hole.h"
#include "core/render/filter/world_replicate.h"
#include "core/render/filter/world_vertex_replicate.h"
#include "core/render/filter/world_mobius.h"
#include "core/render/filter/world_trails.h"
#include "core/render/filter/screen_anti_alias.h"
#include "core/render/filter/screen_direct_aa_sink.h"
#include "core/render/filter/screen_trails.h"
#include "core/render/filter/screen_blur.h"
#include "core/render/filter/pixel_chromatic_shift.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/feedback_style.h"
