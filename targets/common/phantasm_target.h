/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Build-time budgets shared by the Phantasm-class target sketches
 * (targets/Phantasm, targets/Profile).
 */
#pragma once

#include <stddef.h>

/** Per-effect heap-object budget for the Phantasm playlist, in bytes. */
static constexpr size_t HS_PHANTASM_EFFECT_HEAP_BYTES = 3584;
