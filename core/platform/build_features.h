/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file build_features.h
 * @brief Compile-time feature switches, canvas dimensions and size budgets the
 *        rest of the engine derives from, defaulted per target and validated
 *        here.
 */

/** @brief Canvas width in pixels (override via -DCANVAS_W). */
#ifndef CANVAS_W
#define CANVAS_W 288
#endif
/** @brief Canvas height in pixels (override via -DCANVAS_H). */
#ifndef CANVAS_H
#define CANVAS_H 144
#endif

// Size of the real device arena block. Deliberately not overridable: host
// harnesses widen HS_GLOBAL_ARENA_BYTES, so this is the figure memory.h ties
// DEVICE_GLOBAL_ARENA_SIZE to.
#define HS_DEVICE_ARENA_BYTES 305152

#ifndef HS_GLOBAL_ARENA_BYTES
#define HS_GLOBAL_ARENA_BYTES HS_DEVICE_ARENA_BYTES
#endif

#ifndef HS_TIMELINE_MAX_ANIM_BYTES
#define HS_TIMELINE_MAX_ANIM_BYTES 112
#endif

#ifndef HS_ENABLE_EFFECT_REGISTRY
#if defined(__EMSCRIPTEN__)
#define HS_ENABLE_EFFECT_REGISTRY 1
#else
#define HS_ENABLE_EFFECT_REGISTRY 0
#endif
#endif

#ifndef HS_ENABLE_PARAM_GUI_BRIDGE
#if defined(__EMSCRIPTEN__)
#define HS_ENABLE_PARAM_GUI_BRIDGE 1
#else
#define HS_ENABLE_PARAM_GUI_BRIDGE 0
#endif
#endif

// Arena-backed ParamDef storage for effects that outgrow ParamList's inline
// array (Effect::use_parameter_storage). Changes ParamList/Effect layout, so it
// must hold for every TU in an image.
#ifndef HS_EXTERNAL_PARAM_STORAGE
#if defined(ARDUINO)
#define HS_EXTERNAL_PARAM_STORAGE 1
#else
#define HS_EXTERNAL_PARAM_STORAGE 0
#endif
#endif

#ifndef HS_ENABLE_TEST_HOOKS
#define HS_ENABLE_TEST_HOOKS 0
#endif

#ifndef HS_ENABLE_TEST_ORACLES
#define HS_ENABLE_TEST_ORACLES 0
#endif

#ifndef HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND
#if defined(__EMSCRIPTEN__)
#define HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND 1
#else
#define HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND 0
#endif
#endif

#ifndef HS_ENABLE_SHADER_WORKBENCH
#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_ORACLES
#define HS_ENABLE_SHADER_WORKBENCH 1
#else
#define HS_ENABLE_SHADER_WORKBENCH 0
#endif
#endif

#ifndef HS_ENABLE_CHAIN_INTERPRETER
#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_ORACLES
#define HS_ENABLE_CHAIN_INTERPRETER 1
#else
#define HS_ENABLE_CHAIN_INTERPRETER 0
#endif
#endif

#ifndef HS_ENABLE_STRUCTURAL_AUDITS
#define HS_ENABLE_STRUCTURAL_AUDITS 0
#endif

#ifndef HS_ENABLE_EFFECT_CONTROL_API
#if defined(HS_PROFILE_ENABLE)
#define HS_ENABLE_EFFECT_CONTROL_API 1
#else
#define HS_ENABLE_EFFECT_CONTROL_API 0
#endif
#endif

#if (HS_ENABLE_EFFECT_REGISTRY != 0) && (HS_ENABLE_EFFECT_REGISTRY != 1)
#error "HS_ENABLE_EFFECT_REGISTRY must be 0 or 1"
#endif
#if (HS_ENABLE_PARAM_GUI_BRIDGE != 0) && (HS_ENABLE_PARAM_GUI_BRIDGE != 1)
#error "HS_ENABLE_PARAM_GUI_BRIDGE must be 0 or 1"
#endif
#if (HS_EXTERNAL_PARAM_STORAGE != 0) && (HS_EXTERNAL_PARAM_STORAGE != 1)
#error "HS_EXTERNAL_PARAM_STORAGE must be 0 or 1"
#endif
#if (HS_ENABLE_TEST_HOOKS != 0) && (HS_ENABLE_TEST_HOOKS != 1)
#error "HS_ENABLE_TEST_HOOKS must be 0 or 1"
#endif
#if (HS_ENABLE_TEST_ORACLES != 0) && (HS_ENABLE_TEST_ORACLES != 1)
#error "HS_ENABLE_TEST_ORACLES must be 0 or 1"
#endif
#if (HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND != 0) &&                       \
    (HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND != 1)
#error "HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND must be 0 or 1"
#endif
#if (HS_ENABLE_SHADER_WORKBENCH != 0) && (HS_ENABLE_SHADER_WORKBENCH != 1)
#error "HS_ENABLE_SHADER_WORKBENCH must be 0 or 1"
#endif
#if defined(ARDUINO) && HS_ENABLE_SHADER_WORKBENCH
#error "Shader workbench is simulator-only"
#endif
#if defined(ARDUINO) && HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND
#error "HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND is simulator-only"
#endif
#if (HS_ENABLE_CHAIN_INTERPRETER != 0) && (HS_ENABLE_CHAIN_INTERPRETER != 1)
#error "HS_ENABLE_CHAIN_INTERPRETER must be 0 or 1"
#endif
#if defined(ARDUINO) && HS_ENABLE_CHAIN_INTERPRETER
#error "Chain interpreter is simulator-only"
#endif
#if (HS_ENABLE_STRUCTURAL_AUDITS != 0) && (HS_ENABLE_STRUCTURAL_AUDITS != 1)
#error "HS_ENABLE_STRUCTURAL_AUDITS must be 0 or 1"
#endif
#if (HS_ENABLE_EFFECT_CONTROL_API != 0) && (HS_ENABLE_EFFECT_CONTROL_API != 1)
#error "HS_ENABLE_EFFECT_CONTROL_API must be 0 or 1"
#endif
