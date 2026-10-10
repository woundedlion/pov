/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file memory.h
 * @brief Global arena budget and the arena memory facilities.
 */

// platform.h supplies the Arduino IDE/VMicro NDEBUG fallback before <cassert>.
#include "platform/platform.h"
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <new>
#include <cassert>
#include <utility>
#include <concepts>

// DEVICE_GLOBAL_ARENA_SIZE is the device budget on every build target.
/// Device global arena budget, bytes.
constexpr size_t DEVICE_GLOBAL_ARENA_SIZE = HS_DEVICE_ARENA_BYTES;
/// This build's global arena budget, bytes.
constexpr size_t GLOBAL_ARENA_SIZE = HS_GLOBAL_ARENA_BYTES;

/// Default scratch arena A capacity, bytes.
constexpr size_t DEFAULT_SCRATCH_A_SIZE = 16 * 1024;
/// Default scratch arena B capacity, bytes.
constexpr size_t DEFAULT_SCRATCH_B_SIZE = 16 * 1024;
static_assert(GLOBAL_ARENA_SIZE >
                  DEFAULT_SCRATCH_A_SIZE + DEFAULT_SCRATCH_B_SIZE,
              "HS_GLOBAL_ARENA_BYTES must exceed the default scratch split");
/// Persistent arena capacity on this build's default split, bytes.
constexpr size_t DEFAULT_PERSISTENT_SIZE =
    GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;
// Persistent budget on the device's default split, on every build target.
/// Device persistent arena budget on the default split, bytes.
constexpr size_t DEVICE_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// Browser budget for footprint checks, independent of the native harness arena.
// Mirror of CMakeLists.txt's WASM arena configuration.
/// Browser (WASM) global arena budget, bytes.
constexpr size_t WASM_GLOBAL_ARENA_SIZE = 512 * 1024;
#if defined(__EMSCRIPTEN__)
static_assert(WASM_GLOBAL_ARENA_SIZE == GLOBAL_ARENA_SIZE,
              "the module's HS_GLOBAL_ARENA_BYTES moved; update "
              "WASM_GLOBAL_ARENA_SIZE");
#endif
/// Browser persistent arena budget on the default split, bytes.
constexpr size_t WASM_PERSISTENT_BUDGET =
    WASM_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

#include "memory/arena.h"
#include "memory/vector.h"
#include "memory/span.h"
#include "memory/scratch.h"
#include "memory/persist.h"
#include "memory/generate.h"
