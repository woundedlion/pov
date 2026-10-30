/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Chained edges share the last slot; unchained edges reserve both slots. */
inline void test_column_chunk_capacity_boundary() {
  constexpr int W = 72, H = 36;
  const ScopedTemporalSampling sampling(false);
  for (bool chained : {true, false}) {
    StubEffect effect(W, H);
    Canvas canvas(effect);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    sink.prepare(canvas);
    Plot::PointRoute<W, H> route;
    route.prepare(canvas, 0.0f);
    ScratchScope storage(scratch_arena_a);
    Plot::Mesh::ColumnChunk<4, float> chunk(scratch_arena_a);
    size_t flushes = 0, stores = 0;
    bool overflow = false;
    const auto at = [&](size_t slot, float) { return chunk.position(slot); };
    const Plot::GeometryCallback geometry{at, {}, true};
    const auto shade = [](float) { return Color4(Pixel(65535, 0, 0), 1.0f); };
    const auto flush = [&] {
      ++flushes;
      if (!overflow)
        chunk.flush<W, H>(
            pipeline, canvas,
            [&](size_t slot) { return chunk.attribute<0>(slot); },
            [](size_t, float, float) {
              return Color4(Pixel(65535, 0, 0), 1.0f);
            },
            geometry, 0.0f);
    };
    const auto store = [&](size_t slot, size_t source) {
      HS_EXPECT_LT(slot, size_t{4});
      if (slot >= 4) {
        overflow = true;
        return;
      }
      if (!chained && source == 3)
        HS_EXPECT_EQ(flushes, size_t{1});
      ++stores;
      chunk.position(slot) =
          math::pixel_to_vector<W, H>(10.0f + float(source), 18.0f);
      chunk.attribute<0>(slot) = float(source);
    };
    chunk.begin_chain(0, flush);
    for (uint16_t edge = 0; edge < 3; ++edge) {
      const uint16_t head = !chained && edge == 2 ? 3 : edge;
      const Plot::Mesh::Edge source{head, uint16_t(head + 1)};
      const math::Vector ends[2] = {
          math::pixel_to_vector<W, H>(10.0f + head, 18.0f),
          math::pixel_to_vector<W, H>(11.0f + head, 18.0f)};
      chunk.plot_edge<W, H>(route, sink, canvas, source, ends, ends, 0, shade,
                            store, flush);
    }
    HS_EXPECT_FALSE(overflow);
    HS_EXPECT_EQ(flushes, chained ? size_t{0} : size_t{1});
    HS_EXPECT_EQ(stores, chained ? size_t{4} : size_t{5});
    flush();
    HS_EXPECT_EQ(flushes, chained ? size_t{1} : size_t{2});
  }
}
