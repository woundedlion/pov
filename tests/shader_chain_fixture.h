/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shader-chain test fixtures.
 */
/** @file
 * @brief Shared shader-chain program and effect fixtures.
 */
#pragma once

#include "core/math/mobius.h"
#include <bit>
#include <cstdio>
#include <cstring>
#include <memory>
#include <limits>
#include <numbers>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/interpreter.h"
#include "tests/composed_frame_fixture.h"
#include "tests/json_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/vec_test_util.h"
#include "tests/test_harness.h"
#include "workbench/shader/chain_host.h"

namespace hs_test::shader_chain_tests {

namespace PB = Pullback;
namespace In = Pullback::Interp;

#include "tests/shader_chain/fixtures.h"
#include "tests/shader_chain/effect_fixture.h"

} // namespace hs_test::shader_chain_tests
