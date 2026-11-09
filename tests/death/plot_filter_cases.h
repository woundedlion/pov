/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Plot filter death cases.

/** @brief Traps when a pose table is built without storage. */
inline void case_pose_table_build_without_storage() {
  Spacetime::PoseTable table;
  const math::Quaternion history[] = {math::Quaternion{}};
  table.build(history, 1);
}

/**
 * @brief Death case: a face vertex index past the edge-dedup bitset must trap.
 * @details The face-walk draw() overload traps on a vertex index beyond the
 *          TriangularBitset<128> capacity.
 */
inline void case_plot_mesh_vertex_over_capacity() {
  constexpr int W = 32, H = 16;
  OverCapacityMockMesh mesh;
  DeathEffect fx;
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::Mesh::draw<W, H>(pipe, c, mesh, [](const math::Vector &, Fragment &) {
  }); // index 130 -> trap
}

/** @brief Timed mesh edges reject a start index outside the vertex pool. */
inline void case_plot_column_edge_start_out_of_range() {
  plot_column_edge_index_fixture(true);
}

/** @brief Timed mesh edges reject an end index outside the vertex pool. */
inline void case_plot_column_edge_end_out_of_range() {
  plot_column_edge_index_fixture(false);
}

/** @brief An open face dual cannot supply a complete wireframe. */
inline void case_plot_four_regular_open_mesh() {
  configure_arenas_default();
  MeshState mesh;
  build_meshstate_solid<Solids::Octahedron>(mesh, persistent_arena,
                                            opaque(true));
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, mesh.faces.size());
  Plot::Mesh::extract_four_regular_edges(mesh, edges, scratch_arena_a);
}

/** @brief A tetrahedron's dual has odd cycles and cannot be two-colored. */
inline void case_plot_four_regular_non_bipartite() {
  configure_arenas_default();
  MeshState mesh;
  build_meshstate_solid<Solids::Tetrahedron>(mesh, persistent_arena);
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, mesh.faces.size());
  Plot::Mesh::extract_four_regular_edges(mesh, edges, scratch_arena_a);
}

/** @brief Face-edge lookup rejects an absent pair. */
inline void case_plot_find_missing_edge() {
  configure_arenas_default();
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, 1);
  edges.push_back({opaque<uint16_t>(0), opaque<uint16_t>(1)});
  const auto index = Plot::Mesh::find_edge_index(edges, opaque<uint16_t>(1),
                                                 opaque<uint16_t>(2));
  (void)index;
}

/** @brief Death case: extract_edges with an over-capacity vertex index must trap. */
inline void case_plot_extract_edges_vertex_over_capacity() {
  OverCapacityMockMesh mesh;
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(scratch_arena_a, 8);
  Plot::Mesh::extract_edges(mesh, edges); // index 130 -> trap
}

/**
 * @brief Death case: a plot rejects a canvas that is not its <W, H>.
 * @details The sink strides the framebuffer by its own W, so a canvas of a
 *          different size writes past the row it means to.
 */
inline void case_plot_canvas_dim_mismatch() {
  constexpr int W = 32, H = 16;
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 2);
  Fragment f;
  f.pos = math::Vector(1, 0, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 1, 0);
  points.push_back(f);

  DeathEffect fx(W, opaque(H + 1));
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(pipe, c, points,
                        [](const math::Vector &, Fragment &) {});
}

/** @brief Default raster sampling rejects a balanced-sampling request. */
inline void case_raster_default_policy_balanced_sampling() {
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 1);
  DeathEffect effect;
  Canvas canvas(effect);
  Pipeline<32, 16> pipeline;
  Plot::rasterize<32, 16>(pipeline, canvas, points,
                          [](const math::Vector &, Fragment &) {},
                          {.balanced_sampling = opaque(true)});
}

inline void case_raster_empty_path_null_shader() {
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 1);
  DeathEffect fx;
  Canvas canvas(fx);
  Pipeline<32, 16> pipe;
  Plot::rasterize<32, 16>(pipe, canvas, points, FragmentShaderRef{});
}

inline void case_raster_single_nonfinite_point() {
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 1);
  Fragment point;
  point.pos.x = opaque(std::numeric_limits<float>::quiet_NaN());
  points.push_back(point);
  DeathEffect fx;
  Canvas canvas(fx);
  Pipeline<32, 16> pipe;
  Plot::rasterize<32, 16>(pipe, canvas, points,
                          [](const math::Vector &, Fragment &) {});
}

/**
 * @brief Death case: a plot window over a multi-segment polyline must trap.
 * @details The window narrows one segment's arc fraction.
 */
inline void case_plot_window_multi_segment() {
  constexpr int W = 32, H = 16;
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 4);
  Fragment f;
  f.pos = math::Vector(1, 0, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 1, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 0, 1);
  points.push_back(f); // 2 segments under a window -> HS_CHECK

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(
      pipe, c, points, [](const math::Vector &, Fragment &) {},
      {.plot_t_start = opaque(0.25f), .plot_t_end = opaque(0.75f)});
}

/** @brief Scratch sizing rejects zero downsample before dividing. */
inline void case_feedback_scratch_zero() {
  (void)Filter::Pixel::Feedback<32, 16>::uncached_scratch_bytes(opaque(0));
}

/**
 * @brief Death case: a feedback downsample that doesn't divide the resolution must trap.
 * @details The trap fires before any scratch allocation, so no buffers are
 *          needed.
 */
inline void case_feedback_downsample_indivisible() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  Canvas c(fx);
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.downsample = opaque(5); // 32 % 5 != 0 -> HS_CHECK
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  pipe.begin_frame(c, 1.0f);
}

inline void case_feedback_uncached_scratch_budget() {
  DeathEffect fx;
  {
    Canvas first(fx);
    first(0, 0) = Pixel(65535, 65535, 65535);
  }
  fx.advance_display();
  Canvas canvas(fx);
  static uint8_t storage[16];
  scratch_arena_a.rebind(storage, sizeof(storage));
  ::Feedback::Style style{};
  Pipeline<32, 16, Filter::Pixel::Feedback<32, 16>> pipe{
      Filter::Pixel::Feedback<32, 16>(style)};
  pipe.begin_frame(canvas, 1.0f);
}

/** @brief Space transform reading each direction's source a column west.
 *  @param v Target direction.
 *  @return Source direction. */
inline math::Vector death_feedback_turn(const math::Vector &v,
                                        const ::Feedback::Style &) {
  return math::rotate(
      v, math::make_rotation(math::Y_AXIS, -2.0f * math::PI_F / 32.0f));
}

/** @brief A lit 32x16 feedback pipeline whose frames split at column time. */
struct DeathSplitFeedback {
  using Terminal = Filter::Pixel::Feedback<32, 16>;
  using Frames = Pipeline<32, 16, Terminal>;
  ::Feedback::Style style = [] {
    ::Feedback::Style s{};
    s.space_fn = &death_feedback_turn;
    return s;
  }();
  DeathEffect fx{32, 16, pipeline_config<Frames>()};
  Frames pipe{Terminal(style)};

  /** @brief Turns on column time and lights two frames of history. */
  DeathSplitFeedback() {
    Spacetime::set_flat_surface(false);
    for (int frame = 0; frame < 2; ++frame) {
      {
        Canvas c(fx);
        for (int y = 0; y < 16; ++y)
          for (int x = 0; x < 32; ++x)
            c(x, y) = Pixel(30000, 30000, 30000);
      }
      fx.advance_display();
    }
  }
};

/** @brief Death case: a frame flushed while a split frame awaits its finish
 *  must trap. */
inline void case_feedback_split_frame_unfinished() {
  DeathSplitFeedback split;
  {
    Canvas c(split.fx);
    (void)split.pipe.begin_split_frame(c, 1.0f);
  }
  split.fx.advance_display();
  Canvas c(split.fx);
  (void)split.pipe.begin_frame(c, 1.0f);
}

/**
 * @brief Death case: retuning a screen trail to a non-positive lifetime must trap.
 * @details flush()'s fade progress divides by the lifetime.
 */
inline void case_screen_trails_set_lifetime_nonpositive() {
  Filter::Screen::Trails<8> trails(4);
  trails.set_lifetime(opaque(0));
}

/**
 * @brief Death case: a world trail lifetime past the ttl byte must trap.
 * @details World::Trails packs the remaining lifetime into a uint8_t ttl.
 */
inline void case_world_trails_lifetime_over_max() {
  Filter::World::Trails<8> trails(opaque(256)); // lifetime > 255 -> HS_CHECK
  (void)trails;
}

/** @brief Death case: seeding a screen trail before init_storage() must trap. */
inline void case_screen_trails_plot_without_storage() {
  Filter::Screen::Trails<8> trails(4); // no init_storage() -> HS_CHECK
  trails.plot(1.0f, 1.0f, Pixel(1, 1, 1), 0.0f, 1.0f,
              [](float, float, const Pixel &, float, float) {});
}

/** @brief Death case: seeding a world trail before init_storage() must trap. */
inline void case_world_trails_plot_without_storage() {
  Filter::World::Trails<8> trails(4); // no init_storage() -> HS_CHECK
  trails.plot(math::Vector(0, 1, 0), Pixel(1, 1, 1), 0.0f, 1.0f,
              [](const math::Vector &, const Pixel &, float, float) {});
}

inline void case_plot_open_loop_seam() {
  constexpr int W = 32, H = 16;
  Fragment seam;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipeline;
  Plot::draw_fragments<W, H>(
      pipeline, canvas, nullptr, [](const math::Vector &, Fragment &) {},
      {.capacity = 2, .loop_seam = &seam}, [](Fragments &) {});
}

inline void case_raster_point_projection_pair_mismatch() {
  float rows[3] = {};
  float cols[3] = {};
  (void)Plot::PointProjections::paired({rows, opaque<size_t>(3)},
                                       {cols, opaque<size_t>(2)});
}

inline void case_raster_nonfinite_point() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope scope(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  Fragment point;
  point.pos = math::Vector(1, 0, 0);
  points.push_back(point);
  point.pos.x = opaque(std::numeric_limits<float>::quiet_NaN());
  points.push_back(point);
  DeathEffect effect(W, H);
  Canvas canvas(effect);
  Pipeline<W, H> pipeline;
  Plot::rasterize<W, H>(pipeline, canvas, points,
                        [](const math::Vector &, Fragment &) {});
}

inline void case_raster_edge_flags_short() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope scope(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 3);
  for (int i = 0; i < 3; ++i) {
    Fragment f;
    f.pos = math::Vector(1, 0, 0);
    points.push_back(f);
  }
  uint8_t flags[1] = {Plot::RasterOptions::EDGE_VISIBLE};
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipeline;
  Plot::rasterize<W, H>(
      pipeline, canvas, points, [](const math::Vector &, Fragment &) {},
      {.projection = Plot::RasterProjection::geodesic(flags)});
}

/**
 * @brief Death case: hoisted point projections shorter than the polyline must
 *        trap.
 * @details rasterize indexes point_rows/point_cols by point index, not edge
 *          index.
 */
inline void case_raster_point_projections_short() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope sc(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 3);
  for (int i = 0; i < 3; ++i) {
    Fragment f;
    f.pos = math::Vector(1, 0, 0);
    points.push_back(f);
  }
  float rows[3] = {0.0f, 0.0f, 0.0f};
  float cols[3] = {0.0f, 0.0f, 0.0f};
  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(
      pipe, c, points, [](const math::Vector &, Fragment &) {},
      {.point_projections = Plot::PointProjections::paired(
           {rows, opaque<size_t>(2)}, {cols, opaque<size_t>(2)})});
}

/**
 * @brief Death case: a negative feedback fade must trap in sync_hue.
 * @details logf of a negative fade yields NaN, which whitens every feedback
 *          pixel. The guard also catches a NaN fade.
 */
inline void case_feedback_negative_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(-0.25f); // logf(negative) -> HS_CHECK
  style.sync_hue();
}

/**
 * @brief Death case: an infinite feedback fade must trap in sync_hue.
 * @details -logf(+INFINITY) is -inf, whose cos/sin reduction is NaN.
 */
inline void case_feedback_infinite_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(std::numeric_limits<float>::infinity());
  style.sync_hue();
}

/** @brief Death case: planar chord storage for no vertices traps. */
inline void case_planar_chords_empty_storage() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, opaque(0));
}

/** @brief Death case: a polyline past the bound chord storage traps. */
inline void case_planar_chords_over_capacity() {
  static uint8_t arena_buf[Plot::PlanarChords<16, 8>::storage_bytes(2)];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  chords.prepare(canvas.clip());
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, opaque(3), math::Basis{}, Color4{},
                     shader);
}

/** @brief Draws a timed star with optional motion preparation.
 * @param mode Zero for clip preparation, one for motion, two for motion then clip.
 */
inline void draw_timed_planar_chords(int mode) {
  constexpr int W = 32, H = 16, VERTICES = 6;
  static uint8_t bytes[Plot::PlanarChords<W, H>::storage_bytes(VERTICES) +
                       sizeof(Fragment) * (VERTICES + 1) + 64];
  Arena arena(bytes, sizeof(bytes));
  Plot::PlanarChords<W, H> chords;
  chords.init_storage(arena, VERTICES);
  Fragments points;
  points.bind(arena, VERTICES + 1);
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  math::Basis projection_basis;
  const auto &plane =
      *Plot::PlanarProjection::edge_basis(basis, 0.2f, projection_basis);
  using Star = Plot::Star<Plot::PlanarProjection>;
  Star::sample_chart_positions(points, chords.chart_x(), chords.chart_y(),
                               basis, 0.2f, 3, 0.0f, Star::radius_trig(0.2f),
                               Star::step_trig(3), plane);
  hs_test::StubEffect effect(W, H);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  sink.prepare(canvas);
  if (mode != 0)
    chords.prepare(canvas, 0.0f);
  if (mode != 1)
    chords.prepare(canvas.clip());
  chords.draw_closed(sink, canvas, points, VERTICES, plane,
                     Color4(Pixel(65535, 65535, 65535), 0.5f),
                     {.start = math::RotationMatrix{},
                      .middle = math::RotationMatrix{},
                      .max_intervals = 2,
                      .alternating_radii = true,
                      .reference_chords = false,
                      .contour_speed = 0.0f});
}

/** @brief A timed chord draw rejects clip-only preparation. */
inline void case_planar_chords_motion_unprepared() {
  draw_timed_planar_chords(0);
}

/** @brief Clip preparation retires earlier motion preparation. */
inline void case_planar_chords_motion_reset() { draw_timed_planar_chords(2); }

/** @brief A chord draw rejects a missing closing vertex. */
inline void case_planar_chords_missing_closing_vertex() {
  static uint8_t arena_buf[1024];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  chords.prepare(canvas.clip());
  Fragments points;
  points.bind(arena, 2);
  points.push_back(Fragment{math::Y_AXIS});
  points.push_back(Fragment{math::Y_AXIS});
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

/** @brief Death case: a chord draw rejects an unprepared direct sink. */
inline void case_planar_chords_sink_unprepared() {
  static uint8_t arena_buf[1024];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  chords.prepare(canvas.clip());
  Fragments points;
  points.bind(arena, 3);
  points.push_back(Fragment{math::Y_AXIS});
  points.push_back(Fragment{math::Y_AXIS});
  for (int i = 0; i < 2; ++i) {
    chords.chart_x()[i] = 0.0f;
    chords.chart_y()[i] = 0.0f;
  }
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

/** @brief Death case: a chord draw before prepare traps. */
inline void case_planar_chords_unprepared() {
  static uint8_t arena_buf[Plot::PlanarChords<16, 8>::storage_bytes(2)];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

inline void case_planar_chords_stale_clip() {
  static uint8_t arena_buf[Plot::PlanarChords<16, 8>::storage_bytes(2)];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  ClipRegion clip = canvas.clip();
  clip.y_start = opaque(1);
  chords.prepare(clip);
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

/** @brief Death case: band-split storage for under two points traps. */
inline void case_planar_band_split_empty_storage() {
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, opaque(1));
}

/** @brief Death case: a band split cannot append to a populated destination. */
inline void case_planar_band_split_nonempty_output() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  out.bind(arena, 5);
  out.push_back(Fragment{});
  Fragments ring;
  band_split.split(out, ring, opaque(1), 1, math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief Death case: a band split past its bound storage traps. */
inline void case_planar_band_split_over_capacity() {
  static uint8_t arena_buf[1024];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  Fragments ring;
  ring.bind(arena, 3);
  for (int i = 0; i < 3; ++i)
    ring.push_back(Fragment{});
  band_split.split(out, ring, opaque(2), 4, math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief A band split rejects a point product beyond the signed-int range. */
inline void case_planar_band_split_large_product() {
  static uint8_t arena_buf[1024];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out, ring;
  ring.bind(arena, 2);
  ring.push_back(Fragment{});
  ring.push_back(Fragment{});
  band_split.split(out, ring, 1, opaque(INT_MAX), math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief Death case: a band split with no pieces per edge traps. */
inline void case_planar_band_split_zero_pieces() {
  static uint8_t arena_buf[8192];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  out.bind(arena, 5);
  Fragments ring;
  ring.bind(arena, 9);
  for (int i = 0; i < 9; ++i)
    ring.push_back(Fragment{});
  band_split.split(out, ring, 8, opaque(0), math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief Death case: a band split over a ring shorter than edges + 1 traps. */
inline void case_planar_band_split_short_ring() {
  static uint8_t arena_buf[8192];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  out.bind(arena, 5);
  Fragments ring;
  ring.bind(arena, 2);
  ring.push_back(Fragment{});
  ring.push_back(Fragment{});
  band_split.split(out, ring, opaque(2), 1, math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/**
 * @brief Death case: a face palette shader crossfade without a palette.
 */
inline void case_face_palette_shader_blend_without_palette() {
  static BakedPalette to;
  FacePaletteShader shader;
  shader.set_blend(opaque<const BakedPalette *>(nullptr), &to, 0.5f, 0.25f);
}

/**
 * @brief Death case: erasing a direct-raster sink prepared for no canvas.
 * @details PipelineRef drops prepared_for(), so a draw taking the erased handle
 *          cannot check the cached base; the erasure checks it instead.
 */
inline void case_pipeline_ref_erase_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  PipelineRef erased(sink, c);
  (void)erased;
}

inline void case_star_mismatched_chart() {
  using Star = Plot::Star<Plot::PlanarProjection>;
  Fragments points;
  float x[10], y[10];
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  const auto wrong = math::make_basis(math::Quaternion(), math::Y_AXIS);
  Star::sample_chart_positions(points, x, y, basis, 0.5f, 5, 0.0f,
                               Star::radius_trig(0.5f), Star::step_trig(5),
                               wrong);
}

inline void case_star_mismatched_radius_cache() {
  using Star = Plot::Star<Plot::PlanarProjection>;
  Fragments points;
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  Star::sample_positions(points, basis, 0.5f, 5, 0.0f, Star::radius_trig(0.75f),
                         Star::step_trig(5));
}

inline void case_star_mismatched_step_cache() {
  using Star = Plot::Star<Plot::PlanarProjection>;
  Fragments points;
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  Star::sample_positions(points, basis, 0.5f, 5, 0.0f, Star::radius_trig(0.5f),
                         Star::step_trig(6));
}

inline void case_screen_storage_twice() {
  static uint8_t storage[8192];
  Arena arena(storage, sizeof(storage));
  Filter::Screen::Trails<4> stage(10);
  stage.init_storage(arena);
  stage.init_storage(arena);
}

inline void case_world_storage_twice() {
  static uint8_t storage[8192];
  Arena arena(storage, sizeof(storage));
  Filter::World::Trails<4> stage(10);
  stage.init_storage(arena);
  stage.init_storage(arena);
}

inline void case_feedback_storage_twice() {
  static uint8_t storage[8192];
  Arena arena(storage, sizeof(storage));
  ::Feedback::Style style{};
  Filter::Pixel::Feedback<16, 8> stage(style);
  stage.init_storage(arena);
  stage.init_storage(arena);
}

inline void case_direct_sink_wrong_dimensions() {
  hs_test::StubEffect effect(16, 8);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<32, 8> sink;
  sink.prepare(canvas);
}

inline void case_direct_sink_unprepared_plot() {
  hs_test::StubEffect effect(16, 8);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.plot(canvas, 2.0f, 2.0f, Pixel(65535, 0, 0), 0, 1);
}

inline void case_direct_sink_unprepared_integer_plot() {
  hs_test::StubEffect effect(16, 8);
  Canvas canvas(effect);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.plot(canvas, 2, 2, Pixel(65535, 0, 0), 0, 1);
}

inline void case_vertex_replicate_short_input() {
  std::array<math::Vector, 2> vertices{};
  Filter::World::VertexReplicate<3> replicate(vertices);
}

/** @brief Worldline frames reject counts exceeding bound storage. */
inline void case_worldline_frame_capacity() {
  uint8_t bytes[1024];
  Arena arena(bytes, sizeof(bytes));
  Spacetime::Worldlines<> lines;
  lines.bind(arena, 1);
  lines.begin_frame(1, false);
  lines.begin_frame(opaque<size_t>(2), false);
}

/** @brief Point routing rejects a plot before preparation. */
inline void case_point_route_without_prepare() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipe;
  Plot::PointRoute<W, H> valid;
  valid.prepare(canvas, 0.0f);
  const math::Vector point(1, 0, 0);
  const auto start = [&] { return point; };
  const auto shade = [] { return Color4(Pixel(1, 0, 0), 1.0f); };
  valid.plot(pipe, canvas, point, start, shade);
  Plot::PointRoute<W, H> unprepared;
  unprepared.plot(pipe, canvas, point, start, shade);
}

/** @brief Timed pipeline erasure rejects an unprepared direct sink. */
inline void case_pipeline_ref_timed_erase_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<W, H> ready;
  ready.prepare(canvas);
  PipelineRef valid(ready, canvas, Spacetime::Timed{});
  (void)valid;
  Filter::Screen::DirectAntiAliasSink<W, H> unprepared;
  PipelineRef invalid(unprepared, canvas, Spacetime::Timed{});
  (void)invalid;
}
