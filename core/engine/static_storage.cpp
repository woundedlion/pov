/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "platform/platform.h"
#include "animation/animation.h"
#include "render/canvas.h"

// The large static buffers below are defined here, not next to their
// declarations: this TU is linked into every target, so gathering them keeps
// every DMAMEM/large-static placement decision in one file the linker map points
// at. Look here, not in timeline.h / canvas.h, for where the storage actually
// lands. The arena block is the one exception -- it is file-local to memory.cpp,
// which partitions it.

/** @brief Shared event array backing every Timeline instance. */
DMAMEM TimelineEvent global_timeline_events[TIMELINE_MAX_EVENTS];
/**
 * @brief Single live-Timeline guard.
 * @details The event array above is shared by every Timeline instance, so the
 * "only one alive" invariant is one global flag.
 */
bool global_timeline_live = false;
/**
 * @brief Shared singleton playhead cursor into global_timeline_events.
 * @details Free global so every Timeline instance reads/writes the same cursor
 * (see timeline.h).
 */
uint32_t global_timeline_t = 0;
/**
 * @brief Shared singleton event count for global_timeline_events.
 * @details Free global so every Timeline instance reads/writes the same count
 * (see timeline.h).
 */
int global_timeline_num_events = 0;
/**
 * @brief Monotonic count of animations rejected because the timeline was full.
 * @details Never reset, including across Timeline instances (see timeline.h).
 */
uint32_t global_timeline_dropped = 0;
/** @brief Front pixel buffer for the double-buffered effect framebuffer. */
DMAMEM Pixel Effect::buffer_a[MAX_W * MAX_H];
/** @brief Back pixel buffer for the double-buffered effect framebuffer. */
DMAMEM Pixel Effect::buffer_b[MAX_W * MAX_H];
/** @brief Single-live-Effect guard for the shared buffer_a/buffer_b (see Effect). */
bool Effect::s_alive = false;
