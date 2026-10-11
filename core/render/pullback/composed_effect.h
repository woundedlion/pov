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

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>

#include "core/control/choreography.h"
#include "core/color/palette_cycler.h"
#include "core/control/registry.h"
#include "core/memory.h"
#include "core/render/scan.h"
#include "core/math/noise_field.h"
#include "core/render/pullback.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/pullback/composed_resources.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct ComposedFrameWhiteBox;
}
#endif

namespace Pullback {

#include "core/render/pullback/composed_providers.h"
#include "core/render/pullback/composed_parameters.h"
#include "core/render/pullback/composed_policies.h"
#include "core/render/pullback/composed_descriptors.h"
#include "core/render/pullback/composed_runtime.h"
} // namespace Pullback
