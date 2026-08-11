/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_GLOBAL_ARENA_BYTES
#define HS_GLOBAL_ARENA_BYTES 305152
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

#ifndef HS_ENABLE_TEST_HOOKS
#define HS_ENABLE_TEST_HOOKS 0
#endif

#ifndef HS_ENABLE_TEST_ORACLES
#define HS_ENABLE_TEST_ORACLES 0
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
#if (HS_ENABLE_TEST_HOOKS != 0) && (HS_ENABLE_TEST_HOOKS != 1)
#error "HS_ENABLE_TEST_HOOKS must be 0 or 1"
#endif
#if (HS_ENABLE_TEST_ORACLES != 0) && (HS_ENABLE_TEST_ORACLES != 1)
#error "HS_ENABLE_TEST_ORACLES must be 0 or 1"
#endif
#if (HS_ENABLE_STRUCTURAL_AUDITS != 0) && (HS_ENABLE_STRUCTURAL_AUDITS != 1)
#error "HS_ENABLE_STRUCTURAL_AUDITS must be 0 or 1"
#endif
#if (HS_ENABLE_EFFECT_CONTROL_API != 0) && (HS_ENABLE_EFFECT_CONTROL_API != 1)
#error "HS_ENABLE_EFFECT_CONTROL_API must be 0 or 1"
#endif
