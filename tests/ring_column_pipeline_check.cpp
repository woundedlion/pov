/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/engine/engine.h"
#include "core/render/filter.h"
#include "core/render/plot/shapes.h"

#if defined(HS_RING_COLUMN_WORLD)
using RingPipeline =
    Pipeline<8, 8, Filter::World::Trails<16>, Filter::Screen::AntiAlias<8, 8>>;
#elif defined(HS_RING_COLUMN_DIRECT)
using RingPipeline = Filter::Screen::DirectAntiAliasSink<8, 8>;
#else
using RingPipeline = Pipeline<8, 8, Filter::Screen::AntiAlias<8, 8>>;
#endif

void draw_ring(RingPipeline &pipeline, Canvas &canvas,
               const Plot::GeometryCallback &geometry) {
  Plot::draw_ring_columns<8, 8>(pipeline, canvas, geometry,
                                [](const math::Vector &, Fragment &) {});
}
