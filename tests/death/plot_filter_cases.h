/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Plot filter death fixtures and guard cases.

inline void case_scan_ring_stack_too_many_slots() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(b, 1.0f, 0.05f, knots, 4, 0.0f, nullptr);
  const int8_t slot_by_ring[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, slot_by_ring, opaque(INT8_MAX + 1), table,
      [](int, const math::Vector &, Fragment &) {});
}

inline void case_scan_ring_stack_too_many_rings() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(b, 1.0f, 0.05f, knots, 4, 0.0f, nullptr);
  int8_t slot_by_ring[256];
  for (int8_t &s : slot_by_ring)
    s = -1;
  slot_by_ring[0] = 0;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, opaque(256), &ring, slot_by_ring, 1, table,
      [](int, const math::Vector &, Fragment &) {}); // 256 > 255 -> trap
}

inline void case_scan_ring_stack_callback_ring() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::DistortedRing ring(
      b, 1.0f, 0.05f, [](float) { return 0.0f; }, opaque(0.01f), 0.0f);
  const int8_t slot_by_ring[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, slot_by_ring, 1, table,
      [](int, const math::Vector &, Fragment &) {}); // no knots -> trap
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
  if (index == opaque<uint16_t>(42))
    std::printf("x");
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
  static uint8_t arena_buf[256];
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

/** @brief Death case: a chord draw before prepare traps. */
inline void case_planar_chords_unprepared() {
  static uint8_t arena_buf[256];
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
  static uint8_t arena_buf[256];
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
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  Fragments ring;
  band_split.split(out, ring, opaque(2), 4, math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief Death case: the coherent shader rejects a zero block edge. */
inline void case_scan_block_coherent_zero_block() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  const math::Vector position = math::UP;
  Scan::Shader::draw_block_coherent<W, H, 1>(
      canvas, opaque(0), &position, scratch_arena_a,
      [](const math::Vector &) { return Scan::Shader::BlockCell<1>{0}; },
      [](const math::Vector &, const Scan::Shader::BlockCandidates<1> &) {
        return Color4{};
      });
}

/**
 * @brief Death case: the shader rejects a clip beyond its LUT width.
 * @details Direct construction exercises the downstream guard independently
 *          of the Effect setters.
 */
inline void case_scan_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = opaque(W + 1);
  cr.y_end = H;
  cr.margin = 0;
  cr.w = W;
  cr.h = H;
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: the shader rejects a render band past the canvas rows.
 * @details The rows the guard admits subscript the canvas, so a band reaching
 *          past H must be rejected even though the phi LUT holds those rows.
 */
inline void case_scan_clip_rows_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = W;
  cr.y_end = opaque(H + 1);
  cr.margin = 0;
  cr.w = W;
  cr.h = opaque(H + 1);
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later Face claimed.
 * @details The second build retargets the first face's spans, so the first no
 *          longer describes its own geometry.
 */
inline void case_face_scratch_retargeted() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[6];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 6; ++i) {
    float a = (2.0f * math::PI_F * i) / 6.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
  }
  for (int i = 0; i < 3; ++i) {
    idx_a[i] = static_cast<uint16_t>(i);
    idx_b[i] = static_cast<uint16_t>(i + 3);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face first(std::span<const math::Vector>(verts, 6),
                  std::span<const uint16_t>(idx_a, 3), scratch, HV, H);
  SDF::Face second(std::span<const math::Vector>(verts, 6),
                   std::span<const uint16_t>(idx_b, 3), scratch, HV, H);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/** @brief Death case: a virtual-height Face must use the active display grid. */
inline void case_face_virtual_height_mismatched_geometry() {
  constexpr int H = 16;
  const math::Vector verts[] = {math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const uint16_t indices[] = {0, 1, 2};
  static SDF::FaceScratchBuffer scratch;
  SDF::Face face(verts, indices, scratch, opaque(H + hs::H_OFFSET + 1), H);
  (void)face;
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later, culled Face
 *        claimed.
 * @details The second build fills the scratch buffer and only then culls on
 *          collapsed area, so the first face's spans are retargeted even though
 *          nothing will be drawn for the second.
 */
inline void case_face_scratch_retargeted_by_culled_face() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[4];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * math::PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx_a[i] = static_cast<uint16_t>(i);
    // Coincident vertices enclose no area: culled after the scratch fill.
    idx_b[i] = 3;
  }
  verts[3] = basis.v;
  static SDF::FaceScratchBuffer scratch;
  SDF::Face first(std::span<const math::Vector>(verts, 4),
                  std::span<const uint16_t>(idx_a, 3), scratch, HV, H);
  SDF::Face second(std::span<const math::Vector>(verts, 4),
                   std::span<const uint16_t>(idx_b, 3), scratch, HV, H);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/**
 * @brief Death case: a scan rejects a canvas that is not its <W, H>.
 * @details Direct construction exercises the guard independently of the draw
 *          primitives that call it.
 */
inline void case_scan_canvas_dim_mismatch() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, opaque(H + 1));
  Canvas c(fx);
  Scan::check_canvas_dims<W, H>(c);
}

/**
 * @brief Death case: a scan rejects a direct-raster sink prepared for no
 *        canvas.
 * @details The sink writes through its cached framebuffer base, which the
 *          double-buffered canvas leaves pointing at the buffer being scanned
 *          out until prepare() runs.
 */
inline void case_scan_pipeline_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  Scan::check_pipeline_prepared(sink, c);
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

/** @brief Death case: absent face offsets must trap. */
inline void case_scan_mesh_missing_offsets() {
  scan_mesh_invalid_fixture(true);
}

/** @brief Death case: a face index past the vertex pool must trap. */
inline void case_scan_mesh_face_index_out_of_range() {
  scan_mesh_invalid_fixture(false);
}

/**
 * @brief Death case: a class bake naming a class the class table does not hold.
 * @details ArenaVector::operator[] only asserts; Scan::Mesh bounds the id
 *          per face.
 */
inline void case_scan_mesh_class_id_out_of_range() {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[4096];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  mesh.vertices.bind(geom, 3);
  mesh.vertices.push_back(math::Vector(1.0f, 0.0f, 0.0f));
  mesh.vertices.push_back(math::Vector(0.0f, 1.0f, 0.0f));
  mesh.vertices.push_back(math::Vector(0.0f, 0.0f, 1.0f));
  mesh.face_counts.bind(geom, 1);
  mesh.face_counts.push_back(static_cast<uint8_t>(3));
  mesh.face_offsets.bind(geom, 1);
  mesh.face_offsets.push_back(static_cast<uint16_t>(0));
  mesh.faces.bind(geom, 3);
  for (uint16_t v = 0; v < 3; ++v)
    mesh.faces.push_back(v);

  MeshOps::MeshClassBake bake;
  bake.classes.bind(geom, 1);
  bake.face_recs.bind(geom, 1);
  bake.face_recs.push_back({opaque<uint8_t>(0), 0, 0}); // no class founded

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Scan::Mesh::draw<W, H>(
      pipe, c, mesh, [](const math::Vector &, Fragment &) {}, scratch, &bake);
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
  sink.plot(canvas, 2, 2, Pixel(65535, 0, 0), 0, 1);
}

inline void case_vertex_replicate_short_input() {
  std::array<math::Vector, 2> vertices{};
  Filter::World::VertexReplicate<3> replicate(vertices);
}
