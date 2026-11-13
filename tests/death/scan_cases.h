/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Scan and radial-surface death cases.

/**
 * @brief Builds an invalid scan mesh with either missing offsets or a bad index.
 * @param omit_offsets Selects the missing-offset guard before index validation.
 */
inline void scan_mesh_invalid_fixture(bool omit_offsets) {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[1024];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  mesh.vertices.bind(geom, 3);
  for (int i = 0; i < 3; ++i)
    mesh.vertices.push_back(math::Vector(0.0f, 1.0f, 0.0f));
  mesh.face_counts.bind(geom, 1);
  mesh.face_counts.push_back(opaque<uint8_t>(3));
  if (!omit_offsets) {
    mesh.face_offsets.bind(geom, 1);
    mesh.face_offsets.push_back(opaque<uint16_t>(0));
  }
  mesh.faces.bind(geom, 3);
  mesh.faces.push_back(opaque<uint16_t>(0));
  mesh.faces.push_back(opaque<uint16_t>(1));
  mesh.faces.push_back(opaque<uint16_t>(3));

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Scan::Mesh::draw<W, H>(
      pipe, c, mesh, [](const math::Vector &, Fragment &) {}, scratch);
}

/** @brief Binds a single axis triangle in the supplied arena.
 * @param mesh Receives the triangle.
 * @param arena Owns its vertex and face storage. */
inline void build_death_triangle(MeshState &mesh, Arena &arena) {
  mesh.vertices.bind(arena, 3);
  mesh.vertices.push_back(math::X_AXIS);
  mesh.vertices.push_back(math::Y_AXIS);
  mesh.vertices.push_back(math::Z_AXIS);
  mesh.face_counts.bind(arena, 1);
  mesh.face_counts.push_back(3);
  mesh.face_offsets.bind(arena, 1);
  mesh.face_offsets.push_back(0);
  mesh.faces.bind(arena, 3);
  for (uint16_t index = 0; index < 3; ++index)
    mesh.faces.push_back(index);
}

/** @brief Controls the scan selector and malformed face data. */
struct ScanSelectorIndexOptions {
  bool visibility = false; ///< Selects visibility instead of moving faces.
  uint16_t bad_index = 0;  ///< Invalid face vertex index.
  bool omit_offsets =
      false; ///< Removes face offsets instead of corrupting an index.
};

/** @brief Selects a valid mesh before rejecting malformed face data.
 * @param options Selector and malformed face data controls. */
inline void
scan_selector_index_fixture(const ScanSelectorIndexOptions &options) {
  const auto &[visibility, bad_index, omit_offsets] = options;
  constexpr int W = 32, H = 16;
  uint8_t geometry_storage[1024], frame_storage[4096];
  Arena geometry(geometry_storage, sizeof(geometry_storage));
  Arena frame(frame_storage, sizeof(frame_storage));
  MeshState mesh;
  build_death_triangle(mesh, geometry);
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Scan::ColumnBasis basis;
  uint8_t moving[3] = {1, 1, 1};
  uint16_t columns[2] = {};
  const auto select = [&] {
    if (visibility)
      Scan::select_mesh_visibility<W, H>(
          canvas.clip(), mesh, [](const math::Vector &p) { return p; }, 0.0f,
          frame, basis);
    else
      Scan::pick_moving_faces<W>(canvas, mesh, moving, 0.0f, columns);
  };
  select();
  death_admitted();
  if (omit_offsets)
    mesh.face_offsets.clear();
  else
    mesh.faces[2] = opaque<uint16_t>(bad_index);
  select();
}

/** @brief Radial surface family shared by the radial surface death cases. */
using DeathRadialSurface = Scan::RadialSurface<4, 8>;

inline void case_scan_ring_stack_too_many_slots() {
  constexpr int W = 32, H = 16;
  const float KNOTS[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, 1.0f, 0.05f, KNOTS, 4, 0.0f,
                          nullptr);
  const int8_t SLOT_BY_RING[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, SLOT_BY_RING, opaque(INT8_MAX + 1), table,
      [](int, const math::Vector &, Fragment &) {});
}

inline void case_scan_ring_stack_too_many_rings() {
  constexpr int W = 32, H = 16;
  const float KNOTS[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(DEATH_IDENTITY_BASIS, 1.0f, 0.05f, KNOTS, 4, 0.0f,
                          nullptr);
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
  SDF::DistortedRing ring(
      DEATH_IDENTITY_BASIS, 1.0f, 0.05f, [](float) { return 0.0f; },
      opaque(0.01f), 0.0f);
  const int8_t SLOT_BY_RING[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, SLOT_BY_RING, 1, table,
      [](int, const math::Vector &, Fragment &) {}); // no knots -> trap
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
 * @brief Death case: a scan rejects an empty time-aware shader.
 */
inline void case_scan_timed_shader_empty() {
  Scan::check_fragment_shader(TimedShaderRef(nullptr));
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
 * @brief Death case: vertex deformation without the frame's column-time
 *        state.
 */
inline void case_scan_deformation_without_column_time() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  MeshState mesh;
  const uint16_t COLUMNS[2] = {0, 0};
  Scan::deform_mesh<W, H>(
      c, mesh, true, COLUMNS, nullptr, nullptr, scratch_arena_a, 0.1f,
      [](size_t, float) { return math::Vector(1, 0, 0); }, [] {});
}

/** @brief Death case: absent face offsets must trap. */
inline void case_scan_mesh_missing_offsets() {
  scan_mesh_invalid_fixture(true);
}

/** @brief Visibility selection rejects a mesh without face offsets. */
inline void case_scan_visibility_missing_offsets() {
  scan_selector_index_fixture(
      {.visibility = true, .bad_index = 2, .omit_offsets = true});
}

/** @brief Moving-face selection rejects a mesh without face offsets. */
inline void case_scan_moving_missing_offsets() {
  scan_selector_index_fixture(
      {.visibility = false, .bad_index = 2, .omit_offsets = true});
}

/** @brief Death case: a face index past the vertex pool must trap. */
inline void case_scan_mesh_face_index_out_of_range() {
  scan_mesh_invalid_fixture(false);
}

/** @brief Visibility selection rejects an index just past the vertex pool. */
inline void case_scan_visibility_face_index_out_of_range() {
  scan_selector_index_fixture({.visibility = true, .bad_index = 3});
}

/** @brief Visibility selection rejects the largest face index. */
inline void case_scan_visibility_face_index_max() {
  scan_selector_index_fixture({.visibility = true, .bad_index = UINT16_MAX});
}

/** @brief Moving-face selection rejects an index just past the vertex pool. */
inline void case_scan_moving_face_index_out_of_range() {
  scan_selector_index_fixture({.visibility = false, .bad_index = 3});
}

/** @brief Moving-face selection rejects the largest face index. */
inline void case_scan_moving_face_index_max() {
  scan_selector_index_fixture({.visibility = false, .bad_index = UINT16_MAX});
}

/**
 * @brief Death case: a class bake naming a class the class table does not hold.
 */
inline void case_scan_mesh_class_id_out_of_range() {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[4096];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  build_death_triangle(mesh, geom);

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

/** @brief Traps when visibility bits belong to a differently sized mesh. */
inline void case_scan_mesh_visibility_dimensions() {
  uint8_t bits = 1;
  Scan::detail::MeshPass pass;
  pass.visible_faces = &bits;
  pass.visible_vertices = &bits;
  pass.visibility_face_count = 1;
  pass.visibility_vertex_count = 3;
  Scan::check_mesh_pass(pass, 1, 3);
  death_admitted();
  Scan::check_mesh_pass(pass, 2, 3);
}

/** @brief Traps when a deformation pass belongs to a differently sized mesh. */
inline void case_scan_mesh_selection_dimensions() {
  uint8_t selected = 0;
  Scan::detail::MeshPass pass;
  pass.selected_faces = &selected;
  pass.selection_face_count = 1;
  pass.selection_vertex_count = 3;
  Scan::check_mesh_pass(pass, 1, 3);
  death_admitted();
  Scan::check_mesh_pass(pass, 1, 4);
}

/** @brief Active deformation rejects an absent azimuth table. */
inline void case_scan_deformation_without_azimuth() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  MeshState mesh;
  const uint16_t COLUMNS[2] = {};
  Scan::ColumnBasis basis;
  Scan::ColumnBasisScope scope(canvas, basis);
  const auto vertex = [](size_t, float) { return math::Vector(1, 0, 0); };
  Scan::deform_mesh<W, H>(canvas, mesh, false, COLUMNS, nullptr, nullptr,
                          scratch_arena_a, 0.1f, vertex, [] {});
  death_admitted();
  Scan::deform_mesh<W, H>(canvas, mesh, true, COLUMNS, nullptr, nullptr,
                          scratch_arena_a, 0.1f, vertex, [] {});
}

/** @brief Moving mesh planning rejects an absent column azimuth table. */
inline void case_scan_motion_without_azimuth() {
  constexpr int W = 32, H = 16;
  uint8_t geometry_storage[1024], frame_storage[32768], start_storage[4096];
  Arena geometry(geometry_storage, sizeof(geometry_storage));
  Arena frame(frame_storage, sizeof(frame_storage));
  Arena start(start_storage, sizeof(start_storage));
  MeshState mesh;
  build_death_triangle(mesh, geometry);
  const auto rotation = math::make_rotation(math::Z_AXIS, 0.2f);
  math::Vector previous[3];
  for (size_t index = 0; index < 3; ++index)
    previous[index] = math::rotate(mesh.vertices[index], rotation);
  const Scan::ModelPath path{mesh.vertices.data(), previous};
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  math::Orientation<> orientation;
  auto basis = Scan::prepare_column_basis<W, H>(frame, orientation, true);
  Scan::ColumnBasisScope scope(canvas, basis);
  const Scan::NoRipples ripples;
  const auto camera = Scan::EndpointCamera::none();
  {
    ScratchScope frame_scope(frame), start_scope(start);
    const auto plan = Scan::plan_motion<W>(canvas, mesh, path, 0.2f, ripples,
                                           camera, frame, start);
    if (plan.active)
      death_admitted();
  }
  basis.azimuth = nullptr;
  Scan::plan_motion<W>(canvas, mesh, path, 0.2f, ripples, camera, frame, start);
}

/** @brief Solid rasterization rejects an installed column-time basis. */
inline void case_scan_solid_with_column_time() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipeline;
  const auto basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 0, -1));
  const SDF::PlanarPolygon shape(basis, 0.2f, 5, 0.0f);
  const Color4 color(Pixel(1, 0, 0), 1.0f);
  Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color);
  Scan::ColumnBasis time;
  Scan::ColumnBasisScope scope(canvas, time);
  Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color);
}

/** @brief Traps on a radial surface with more sides than its capacity. */
inline void case_radial_surface_side_capacity() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  surface.bind(scratch_arena_a, opaque(9), 0.1f, 1);
}

/** @brief Traps on a radial surface with no sweep intervals. */
inline void case_radial_surface_interval_capacity() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  surface.bind(scratch_arena_a, 4, 0.1f, opaque(0));
}

/** @brief Traps on an inner ratio past the convex limit cos(pi/sides). */
inline void case_radial_surface_inner_ratio() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  surface.bind(scratch_arena_a, 4, opaque(0.9f), 1);
}

/** @brief Traps on a replacement sweep whose knot count mismatches. */
inline void case_radial_surface_knot_count() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  surface.bind(scratch_arena_a, 4, 0.1f, 1);
  const DeathRadialSurface::Knot KNOTS[3] = {{1, 0, 0}, {1, 0, 0}, {1, 0, 0}};
  surface.set_knots({KNOTS, size_t(opaque(3))});
}

/** @brief Traps on a zero multiplier after the first knot. */
inline void case_radial_surface_noninitial_zero_multiplier() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  const DeathRadialSurface::Knot knots[2] = {{1, 0, 0}, {opaque(0.0f), 0, 0}};
  surface.prepare(scratch_arena_a, 4, 0.1f, knots);
}

/** @brief Traps on adjacent knot phases a full edge sector apart. */
inline void case_radial_surface_interval_crosses_sectors() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  const float phase = opaque(math::PI_F / 4);
  const DeathRadialSurface::Knot knots[2] = {{1, 0, 0},
                                             {cosf(phase), sinf(phase), phase}};
  surface.prepare(scratch_arena_a, 4, 0.1f, knots);
}

/** @brief Traps when an unbound candidate mask is prepared. */
inline void case_radial_surface_mask_unbound() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  const DeathRadialSurface::Knot KNOTS[2] = {{1, 0, 0}, {1, 0, 0}};
  surface.prepare(scratch_arena_a, 4, 0.1f, KNOTS);
  Scan::RadialSurfaceCandidateMask<> mask;
  (void)mask.prepare(surface, 0.0f);
}

/** @brief Traps on a candidate mask with a negative footprint radius. */
inline void case_radial_surface_mask_negative_footprint() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  const DeathRadialSurface::Knot KNOTS[2] = {{1, 0, 0}, {1, 0, 0}};
  surface.prepare(scratch_arena_a, 4, 0.1f, KNOTS);
  Scan::RadialSurfaceCandidateMask<> mask;
  mask.bind(scratch_arena_a);
  (void)mask.prepare(surface, opaque(-0.01f));
}

/** @brief Traps on a candidate mask with a negative complex motion bound. */
inline void case_radial_surface_mask_negative_motion() {
  ScratchScope scope(scratch_arena_a);
  DeathRadialSurface surface;
  const DeathRadialSurface::Knot KNOTS[2] = {{1, 0, 0}, {1, 0, 0}};
  surface.prepare(scratch_arena_a, 4, 0.1f, KNOTS);
  Scan::RadialSurfaceCandidateMask<> mask;
  mask.bind(scratch_arena_a);
  (void)mask.prepare(surface, 0.0f, opaque(-0.01f));
}

/** @brief Constructs one flat knot ring in a pool.
 * @param pool Pool receiving the ring. */
template <typename Pool> inline void commit_stack_ring(Pool &pool) {
  auto row = pool.next();
  for (int k = 0; k <= 16; ++k)
    row.shift_row[k] = 0.0f;
  const auto basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  pool.commit(0, 1.0f, 16, basis, 0.5f, 0.03f, row.shift_row, 16, 0.0f,
              nullptr);
}
/** @brief A ring pool cannot begin over live rings. */
inline void case_scan_ring_pool_live_begin() {
  using Pool = Scan::RingStackPool<1, 33>;
  alignas(std::max_align_t) uint8_t buf[Pool::BYTES];
  Arena arena(buf, sizeof(buf));
  Pool pool;
  pool.init_storage(arena);
  pool.begin_frame();
  commit_stack_ring(pool);
  pool.begin_frame();
}
/** @brief A full ring pool has no pending slot. */
inline void case_scan_ring_pool_full() {
  using Pool = Scan::RingStackPool<1, 33>;
  alignas(std::max_align_t) uint8_t buf[Pool::BYTES];
  Arena arena(buf, sizeof(buf));
  Pool pool;
  pool.init_storage(arena);
  pool.begin_frame();
  commit_stack_ring(pool);
  static_cast<void>(pool.next());
}
/** @brief Alternate storage must fit a complete knot row. */
inline void case_scan_alternates_small_capacity() {
  alignas(std::max_align_t) uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Scan::KnotAlternates<32> alternates;
  alternates.allocate(arena, arena, 2, 32);
  alternates.allocate(arena, arena, 2, opaque(31));
}
/** @brief Bound alternates must cover every pool slot. */
inline void case_scan_ring_pool_few_alternate_slots() {
  using Pool = Scan::RingStackPool<2, 33>;
  alignas(std::max_align_t) uint8_t buf[Pool::BYTES + 4096];
  Arena arena(buf, sizeof(buf));
  typename Pool::Alternates alternates;
  alternates.allocate(arena, arena, 2, 32);
  Pool pool;
  pool.init_storage(arena);
  pool.begin_frame(&alternates);
  alternates.allocate(arena, arena, opaque(1), 32);
  pool.begin_frame(&alternates);
}

/** @brief Packed offsets must fit the prefix representation. */
inline void case_scan_alternates_prefix_overflow() {
  alignas(std::max_align_t) uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  Scan::KnotAlternates<32> alternates;
  alternates.allocate(arena, arena, 1, 32);
  alternates.allocate(arena, arena, 1, opaque(65536));
}
