/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "platform/platform.h"
#include "animation/animation.h"
#include "render/canvas.h"

#include "engine/concepts.h"
#include "platform/inplace_function.h"

// Engine timeline and framebuffer storage.

/** @brief NOLOAD DMAMEM storage; Timeline clears every slot at runtime. */
DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
/** @brief Single live-Timeline guard for the shared event array. */
bool global_timeline_live = false;
/** @brief Shared singleton playhead cursor into global_timeline_events. */
uint32_t global_timeline_t = 0;
/** @brief Shared singleton event count for global_timeline_events. */
int global_timeline_num_events = 0;
/**
 * @brief Modulo-2^32 count of animations rejected because the timeline was full.
 * @details Never reset, including across Timeline instances.
 */
uint32_t global_timeline_dropped = 0;
/**
 * @brief Whether the current saturation episode has already logged a drop.
 * @details Cleared whenever the event table empties.
 */
bool global_timeline_drop_logged = false;
/** @brief First rotating pixel buffer for the effect framebuffer. */
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
/** @brief Second rotating pixel buffer for the effect framebuffer. */
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];
/** @brief Single-live-Effect guard for the shared buffer_a/buffer_b. */
bool Effect::s_alive = false;

#ifdef ARDUINO
#include <exception>
#endif

// Device handlers suppress the C++ demangler; host builds use toolchain handlers.
#ifdef ARDUINO
/**
 * @brief Fail-fast handler for a pure-virtual call on the device.
 * @details Flushes the log then traps.
 */
extern "C" void __cxa_pure_virtual() {
  hs::flush_log();
  __builtin_trap();
}
#if TEENSYDUINO < 162
namespace __gnu_cxx {
/**
 * @brief Fail-fast terminate handler for the device.
 * @details Flushes the log then traps; avoids linking the C++ demangler.
 */
void __verbose_terminate_handler() {
  hs::flush_log();
  __builtin_trap();
}
} // namespace __gnu_cxx
#else
// Teensyduino >= 1.62 supplies a strong handler; override it via set_terminate
// before setup().
namespace {
const std::terminate_handler REPLACED_TERMINATE_HANDLER =
    std::set_terminate([] {
      hs::flush_log();
      __builtin_trap();
    });
} // namespace
#endif
#endif

namespace hs {
[[noreturn]] HS_COLD void function_ref_empty_call() {
  check_fail(HS_CHECK_SITE("thunk != empty_thunk"), "empty FunctionRef called");
}
[[noreturn]] HS_COLD void inplace_function_empty_call() {
  check_fail(HS_CHECK_SITE("vtable != empty"),
             "empty hs::inplace_function called");
}
} // namespace hs
