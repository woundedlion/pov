/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file composed_effect.h
 * @brief Shared machinery for the composed-effect family: the parameter
 *        providers and present-only instance storage derived from a ranked
 *        stage Spec, with the engine's preset choreography and palette lifecycle.
 */

#include "animation/orientation.h"
#include "math/mobius.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>

#include "color/effect_palette_recipes.h"
#include "control/choreography.h"
#include "color/palette_cycler.h"
#include "control/registry.h"
#include "memory.h"
#include "render/scan.h"
#include "math/noise_field.h"
#include "render/pullback.h"
#include "render/pullback/runtime_seeds.h"
#include "render/pullback/composed_resources.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct ComposedFrameWhiteBox;
}
#endif

namespace Pullback {

#include "render/pullback/composed_providers.h"
#include "render/pullback/composed_parameters.h"
#include "render/pullback/composed_policies.h"
#include "render/pullback/composed_descriptors.h"
#include "render/pullback/composed_runtime.h"
} // namespace Pullback
