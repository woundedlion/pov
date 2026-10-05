/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file memory.h
 * @brief Arena allocator, the engine's global arena budget, the containers
 *        built on top of it, and the scratch-scoped generate() wrapper.
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

// Arena budgets come from platform/platform.h and the build definitions.
// Native 64-bit pointer-containing structs can exceed device footprints.
// DEVICE_GLOBAL_ARENA_SIZE uses the device budget across build targets.
constexpr size_t DEVICE_GLOBAL_ARENA_SIZE = HS_DEVICE_ARENA_BYTES;
constexpr size_t GLOBAL_ARENA_SIZE = HS_GLOBAL_ARENA_BYTES;

constexpr size_t DEFAULT_SCRATCH_A_SIZE = 16 * 1024;
constexpr size_t DEFAULT_SCRATCH_B_SIZE = 16 * 1024;
// An HS_GLOBAL_ARENA_BYTES override at or below the scratch split would wrap the
// unsigned subtraction below, giving the persistent arena a capacity far larger
// than global_arena_block; its two-argument constructor passes size as its own
// extent, so nothing traps and the first allocation writes past the block.
static_assert(GLOBAL_ARENA_SIZE >
                  DEFAULT_SCRATCH_A_SIZE + DEFAULT_SCRATCH_B_SIZE,
              "HS_GLOBAL_ARENA_BYTES must exceed the default scratch split");
constexpr size_t DEFAULT_PERSISTENT_SIZE =
    GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;
// Persistent budget on the real device split (from DEVICE_GLOBAL_ARENA_SIZE, not
// the host-inflated GLOBAL_ARENA_SIZE) so an effect's default-split footprint
// static_assert checks the true device figure even in the host suite.
constexpr size_t DEVICE_PERSISTENT_BUDGET =
    DEVICE_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

// Browser budget for footprint checks, independent of the native harness arena.
// Mirror of CMakeLists.txt's WASM arena configuration.
constexpr size_t WASM_GLOBAL_ARENA_SIZE = 512 * 1024;
#if defined(__EMSCRIPTEN__)
static_assert(WASM_GLOBAL_ARENA_SIZE == GLOBAL_ARENA_SIZE,
              "the module's HS_GLOBAL_ARENA_BYTES moved; update "
              "WASM_GLOBAL_ARENA_SIZE");
#endif
constexpr size_t WASM_PERSISTENT_BUDGET =
    WASM_GLOBAL_ARENA_SIZE - DEFAULT_SCRATCH_A_SIZE - DEFAULT_SCRATCH_B_SIZE;

#include "memory/arena.h"
#include "memory/vector.h"
#include "memory/span.h"
#include "memory/scratch.h"
#include "memory/persist.h"
#include "memory/generate.h"
