/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <string_view>

#include "effects/ShaderBall.h"

#define HS_SHADER_DERIVED_EFFECT_LIST(X) X(Shader)

template <int W, int H> class ShaderWorkbench : public ShaderBall<W, H> {
public:
  static constexpr std::string_view EFFECT_ID = "shader";
};

template <int W, int H> using Shader = ShaderWorkbench<W, H>;

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Shader)
