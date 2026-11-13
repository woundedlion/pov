/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "core/platform/platform.h"
#include "core/animation/animation.h"
#include "core/render/canvas.h"

#include "core/engine/concepts.h"
#include "core/math/lenses.h"
#include "core/platform/inplace_function.h"

/**
 * @file static_storage.cpp
 * @brief Engine storage and fail-fast handlers.
 */

namespace lenses::DodecahedralTable {
namespace {
/** @brief Run-encoded classifier maps. */
constinit const auto FACE_RUNS_0 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs0) =
    make_face_runs<0>();
constinit const auto FACE_RUNS_1 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs1) =
    make_face_runs<1>();
constinit const auto FACE_RUNS_2 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs2) =
    make_face_runs<2>();
constinit const auto FACE_RUNS_3 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs3) =
    make_face_runs<3>();
constinit const auto FACE_RUNS_4 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs4) =
    make_face_runs<4>();
constinit const auto FACE_RUNS_5 HS_PROGMEM_UNIQUE(dodecahedral_fold_runs5) =
    make_face_runs<5>();
} // namespace
constinit const std::array<std::span<const ClassifierRun>, 6>
    CLASSIFIER_RUNS HS_PROGMEM_UNIQUE(dodecahedral_fold_runs) = {
        FACE_RUNS_0, FACE_RUNS_1, FACE_RUNS_2,
        FACE_RUNS_3, FACE_RUNS_4, FACE_RUNS_5};

#ifndef ARDUINO
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_0 HS_PROGMEM_UNIQUE(dodecahedral_fold_face0) = make_face<0>();
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_1 HS_PROGMEM_UNIQUE(dodecahedral_fold_face1) = make_face<1>();
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_2 HS_PROGMEM_UNIQUE(dodecahedral_fold_face2) = make_face<2>();
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_3 HS_PROGMEM_UNIQUE(dodecahedral_fold_face3) = make_face<3>();
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_4 HS_PROGMEM_UNIQUE(dodecahedral_fold_face4) = make_face<4>();
constinit const std::array<uint8_t, FACE_STEPS * FACE_STEPS>
    FACE_5 HS_PROGMEM_UNIQUE(dodecahedral_fold_face5) = make_face<5>();
constinit const std::array<const uint8_t *, 6> FACES HS_PROGMEM_UNIQUE(
    dodecahedral_fold_faces) = {FACE_0.data(), FACE_1.data(), FACE_2.data(),
                                FACE_3.data(), FACE_4.data(), FACE_5.data()};
#endif
} // namespace lenses::DodecahedralTable

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
  hs::check_fail(HS_CHECK_SITE("device runtime termination"),
                 "pure virtual call or termination");
}
#if TEENSYDUINO < 162
namespace __gnu_cxx {
/**
 * @brief Fail-fast terminate handler for the device.
 * @details Flushes the log then traps; avoids linking the C++ demangler.
 */
void __verbose_terminate_handler() {
  hs::check_fail(HS_CHECK_SITE("device runtime termination"),
                 "pure virtual call or termination");
}
} // namespace __gnu_cxx
#else
// Teensyduino >= 1.62 supplies a strong handler; override it via set_terminate
// before setup().
namespace {
const std::terminate_handler REPLACED_TERMINATE_HANDLER =
    std::set_terminate([] {
      hs::check_fail(HS_CHECK_SITE("device runtime termination"),
                     "C++ termination");
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
