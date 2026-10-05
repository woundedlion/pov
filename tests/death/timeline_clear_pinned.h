/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/**
 * @brief Death case: clear()ing a pinned event must trap.
 * @details Animation surface — the third teardown path, alongside
 *          case_timeline_pinned_relocation (move_into) and
 *          case_timeline_pinned_completion (step's destroy branch). The public
 *          clear() would otherwise free an event whose animation pointer the
 *          caller still holds. ~Timeline reaches the same events through the
 *          unguarded reset_storage(), which is safe because no retained handle
 *          spans the instance boundary.
 */
inline void case_timeline_clear_pinned() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear));
  global_timeline_events[0].pinned = opaque(true);
  tl.clear(); // HS_CHECK(!pinned) -> trap
}

/**
 * @brief Death case: clear()ing from a completion callback must trap.
 * @details Animation surface — step() runs post_callback() and only afterwards
 *          destroys the event, so a clear() inside that callback would free the
 *          callable whose frame is still executing. The trap sits at the top of
 *          clear(), ahead of destroy_events().
 */
inline void case_timeline_clear_during_step() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear).then([&tl]() {
    tl.clear();
  }));
  tl.step(canvas); // t=1: completes -> callback -> clear() while stepping
}

/** @brief Death case: finite parameter animations reject the -1 sentinel. */
inline void case_finite_param_perpetual_duration() {
  float value = 0.0f;
  Animation::Transition transition(value, 1.0f, opaque(-1), math::ease_linear);
  if (transition.done())
    std::printf("x");
}

/** @brief Death case: a Transition target must be finite. */
inline void case_transition_nonfinite_target() {
  float value = 0.0f;
  Animation::Transition transition(
      value, opaque(std::numeric_limits<float>::quiet_NaN()), 1,
      math::ease_linear);
  if (transition.done())
    std::printf("x");
}

/** @brief Adds an event from a clear hook, violating the hook contract. */
inline void add_event_from_clear_hook(void *ctx) {
  static float value = 0.0f;
  static_cast<Timeline *>(ctx)->add(
      0, Animation::Transition(value, 1.0f, 1, math::ease_linear));
}

/** @brief Death case: clear hooks must not mutate timeline event storage. */
inline void case_timeline_clear_hook_adds_event() {
  Timeline tl;
  tl.add_clear_hook(&tl, add_event_from_clear_hook);
  tl.clear();
}

/**
 * @brief Death case: scheduling a segue sprite with no free timeline slot must
 *        trap.
 * @details Animation surface — every segue policy's schedule() returns the next
 *          transition's delay whether or not its sprite landed, so a dropped add
 *          leaves the sphere dark for a whole transition while the effect
 *          advances on schedule. The budget guard traps at the schedule.
 */
inline void case_segue_sprite_no_slot() {
  Timeline tl;
  float sink = 0.0f;
  while (Timeline::remaining() > 0)
    tl.add(0, Animation::Transition(sink, 1.0f, 1000, math::ease_linear));
  Segue::schedule_faded_sprite(tl, [](Canvas &, float) {}, 4, 1);
}

/** @brief Death case: a segue must target the already-flipped front slot. */
inline void case_mesh_carousel_unflipped_slot() {
  Timeline tl;
  MeshCarousel<> carousel;
  carousel.schedule_segue(tl, 1, [](Canvas &, float) {}, 4, 1);
}

/**
 * @brief Death case: a second simultaneously-live Timeline must trap.
 * @details Animation surface — every Timeline shares the single global event
 *          array, so a second live instance would silently stomp the first's
 *          events; the construction guard traps instead. The real app holds
 *          exactly one (the old effect is destroyed before the next is built).
 */
inline void case_timeline_double_construct() {
  Timeline a;
  Timeline b; // second live ctor -> HS_CHECK(!global_timeline_live) -> trap
  if (global_timeline_num_events == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: narrowing an index past the int16 topology range must trap.
 * @details Mesh-topology surface — both conway.h and hankin.h route every output
 *          vertex/face-index narrowing through this shared MeshOps guard, so a
 *          future MeshLimits::MAX_VERTICES bump traps at the bench instead of silently wrapping
 *          an index and corrupting topology.
 */
inline void case_mesh_narrow_index() {
  size_t over = static_cast<size_t>(INT16_MAX) + 1;
  uint16_t i = MeshOps::narrow_index(opaque(over)); // > INT16_MAX -> HS_CHECK
  if (i == 0xBEEF)
    std::printf("x");
}

/** @brief Death case: medial rejects an input used as both outputs. */
inline void case_medial_aliases_input() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::medial(mesh, mesh, mesh.vertices, target, temp);
}

/**
 * @brief Death case: needle rejects one arena passed as both target and temp.
 * @details needle applies dual then kis, swapping the arenas between legs, so a
 *          single arena has the second leg reading the block the first is
 *          overwriting.
 */
inline void case_needle_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::needle(mesh, arena, arena);
}

/** @brief Death case: zip rejects one arena passed as both target and temp. */
inline void case_zip_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::zip(mesh, arena, arena);
}

/** @brief Death case: gyro rejects one arena passed as both target and temp. */
inline void case_gyro_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::gyro(mesh, arena, arena);
}

/** @brief Death case: bevel rejects one arena passed as both target and temp. */
inline void case_bevel_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::bevel(mesh, arena, arena);
}

/** @brief Death case: ambo rejects one arena passed as both target and temp. */
inline void case_ambo_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::ambo(mesh, arena, arena);
}

/**
 * @brief Death case: MeshOps::compile rejects one arena passed as both the
 *        geometry arena and the scratch arena.
 */
inline void case_mesh_compile_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh src;
  MeshState dst;
  MeshOps::compile(src, dst, arena, arena);
}

/** @brief Death case: MeshOps::clone rejects a self-aliased destination. */
inline void case_mesh_ops_clone_aliases_dst() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::clone(mesh, mesh, arena);
}

/** @brief Death case: MeshState::clone rejects a self-aliased destination. */
inline void case_mesh_state_clone_aliases_dst() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  MeshState mesh;
  MeshState::clone(mesh, mesh, arena);
}

/**
 * @brief Death case: MeshOps::transform rejects a self-aliased destination.
 * @details set_borrowed() drops the source's owned topology before it is read,
 *          so a self-transform would report F faces over reclaimed arena
 *          bytes.
 */
inline void case_mesh_transform_aliases_source() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  MeshState mesh;
  MeshOps::transform(mesh, mesh, arena);
}

/** @brief Death case: chamfer rejects the face-collapse endpoint. */
inline void case_chamfer_collapsed_endpoint() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::chamfer(mesh, target, temp, opaque(1.0f));
}

/** @brief Death case: snub rejects the face-collapse endpoint. */
inline void case_snub_collapsed_endpoint() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::snub(mesh, target, temp, opaque(1.0f));
}

/** @brief Death case: a Conway morph operator rejects an empty mesh. */
inline void case_conway_empty_mesh() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Death case: a Conway morph operator rejects an open one-face mesh. */
inline void case_conway_degenerate_mesh() {
  static uint8_t source_buf[1024];
  static uint8_t target_buf[4096];
  static uint8_t temp_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  mesh.vertices.bind(source, 3);
  mesh.vertices.push_back(math::Vector(1, 0, 0));
  mesh.vertices.push_back(math::Vector(0, 1, 0));
  mesh.vertices.push_back(math::Vector(0, 0, 1));
  mesh.face_counts.bind(source, 1);
  mesh.face_counts.push_back(3);
  mesh.faces.bind(source, 3);
  mesh.faces.push_back(0);
  mesh.faces.push_back(1);
  mesh.faces.push_back(2);
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Death case: only the target arena is small enough to exhaust. */
inline void case_conway_target_exhausted() {
  static uint8_t source_buf[65536];
  static uint8_t target_buf[16];
  static uint8_t temp_buf[65536];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  build_solid<Solids::Tetrahedron>(mesh, source);
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Vertex-bit payload words a tetrahedron bake needs (3 per vertex). */
inline constexpr size_t TETRAHEDRON_BAKE_WORDS = 3 * 4;

/**
 * @brief Builds a tetrahedron plus a relax bake that matches it exactly.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param bits Payload storage, TETRAHEDRON_BAKE_WORDS words, filled with the
 *        mesh's own vertex bits; must outlive the relax_baked() call.
 * @return Bake relax_baked() accepts until the caller perturbs one field.
 * @details The source mesh's own vertices make both source and payload hashes
 *          match.
 */
inline MeshOps::RelaxBake
build_matching_relax_bake(PolyMesh &mesh, Arena &arena, uint32_t *bits) {
  build_solid<Solids::Tetrahedron>(mesh, arena);
  uint32_t output_hash = MeshOps::FNV1A_BASIS;
  for (size_t i = 0; i < mesh.vertices.size(); ++i) {
    const math::Vector &v = mesh.vertices[i];
    bits[3 * i] = std::bit_cast<uint32_t>(v.x);
    bits[3 * i + 1] = std::bit_cast<uint32_t>(v.y);
    bits[3 * i + 2] = std::bit_cast<uint32_t>(v.z);
    for (size_t k = 0; k < 3; ++k)
      output_hash = MeshOps::fnv1a_step(output_hash, bits[3 * i + k]);
  }
  MeshOps::RelaxBake bake{};
  bake.name = "death_tetrahedron";
  bake.vertex_bits = bits;
  bake.vertex_count = static_cast<uint16_t>(mesh.vertices.size());
  bake.face_count = static_cast<uint16_t>(mesh.get_face_counts_size());
  bake.index_count = static_cast<uint16_t>(mesh.get_faces_size());
  bake.iterations = 0;
  bake.source_hash = MeshOps::relax_source_hash(mesh);
  bake.topology_hash = MeshOps::relax_topology_hash(mesh);
  bake.output_hash = output_hash;
  return bake;
}

/**
 * @brief Death case: relax_baked rejects a bake whose vertex count differs from
 *        the source mesh.
 * @details Baked-payload surface — the dimension check is what stops a payload
 *          baked against different geometry from being read past its end.
 */
inline void case_relax_baked_dimension_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.vertex_count = opaque<uint16_t>(bake.vertex_count + 1);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: relax_baked rejects a bake whose topology hash differs.
 * @details Baked-payload surface — dimensions alone do not pin connectivity, so
 *          this check is what stops a bake replaying onto a mesh of the same
 *          size but different face wiring.
 */
inline void case_relax_baked_topology_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.topology_hash = opaque(bake.topology_hash ^ 1u);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: relax_baked rejects a bake from different source vertices.
 * @details Topology and dimensions do not identify parameterized operators whose
 *          connectivity is invariant across parameter values.
 */
inline void case_relax_baked_source_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  mesh.vertices[0].x += opaque(4.0f / MeshOps::RELAX_SOURCE_SCALE);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: relax_baked rejects a payload whose re-hash differs from
 *        the bake's output hash.
 * @details Baked-payload surface — the only check covering the vertex words
 *          themselves, so a corrupt or truncated flash payload stops here
 *          rather than shipping as geometry.
 */
inline void case_relax_baked_output_hash_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.output_hash = opaque(bake.output_hash ^ 1u);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Builds a one-face PolyMesh with independently sized count and index
 * data.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param side_count Value stored in the face-count array.
 * @param num_indices Number of entries stored in the flat face-index array.
 */
inline void build_mismatched_polymesh(PolyMesh &mesh, Arena &arena,
                                      uint8_t side_count, size_t num_indices) {
  mesh.vertices.bind(arena, 4);
  for (size_t i = 0; i < 4; ++i)
    mesh.vertices.push_back(math::Vector{});
  mesh.face_counts.bind(arena, 1);
  mesh.face_counts.push_back(opaque(side_count));
  mesh.faces.bind(arena, num_indices);
  for (size_t i = 0; i < num_indices; ++i)
    mesh.faces.push_back(static_cast<uint16_t>(i));
}

/**
 * @brief Death case: HalfEdgeMesh rejects counts shorter than its flat index
 * list.
 */
inline void case_half_edge_face_counts_short() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  build_mismatched_polymesh(mesh, arena, 3, opaque<size_t>(4));
  HalfEdgeMesh half_edges(arena, mesh);
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: HalfEdgeMesh rejects counts longer than its flat index
 * list.
 */
inline void case_half_edge_face_counts_long() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  build_mismatched_polymesh(mesh, arena, 4, opaque<size_t>(3));
  HalfEdgeMesh half_edges(arena, mesh);
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: MeshOps::compile rejects counts shorter than its index
 * list.
 */
inline void case_mesh_compile_face_counts_short() {
  static uint8_t src_buf[1024];
  static uint8_t dst_buf[1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh mesh;
  MeshState compiled;
  build_mismatched_polymesh(mesh, src_arena, 3, opaque<size_t>(4));
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: MeshOps::compile rejects counts longer than its index
 * list.
 */
inline void case_mesh_compile_face_counts_long() {
  static uint8_t src_buf[1024];
  static uint8_t dst_buf[1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh mesh;
  MeshState compiled;
  build_mismatched_polymesh(mesh, src_arena, 4, opaque<size_t>(3));
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: MeshOps::compile rejects a face whose span ends past the
 *        16-bit index range even though its start offset still fits.
 */
inline void case_mesh_compile_face_span_over_16bit() {
  static uint8_t src_buf[256 * 1024];
  static uint8_t dst_buf[256 * 1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  // The last triangle starts at 65535 and ends at 65538.
  constexpr size_t FACES = 21846;
  PolyMesh mesh;
  mesh.vertices.bind(src_arena, 3);
  for (int i = 0; i < 3; ++i)
    mesh.vertices.push_back(math::Vector{});
  mesh.face_counts.bind(src_arena, FACES);
  mesh.faces.bind(src_arena, FACES * 3);
  for (size_t f = 0; f < FACES; ++f) {
    mesh.face_counts.push_back(opaque<uint8_t>(3));
    for (uint16_t k = 0; k < 3; ++k)
      mesh.faces.push_back(opaque(k));
  }
  MeshState compiled;
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: update_hankin rejects a retained topology that no
 *        classification of this pattern's output ever produced.
 * @details The topology array survives an angle re-solve on purpose, so a mesh
 *          pointed at a new pattern would otherwise carry class ids that no
 *          longer match the faces written into it.
 */
inline void case_update_hankin_stale_topology() {
  static uint8_t geom_buf[64 * 1024];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, geom, scratch);

  MeshState mesh;
  mesh.topology.bind(geom, compiled.face_counts.size());
  for (size_t i = 0; i < compiled.face_counts.size(); ++i)
    mesh.topology.push_back(opaque<uint16_t>(0));
  MeshOps::update_hankin(compiled, mesh, geom, opaque(0.0f));
  if (mesh.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: update_hankin rejects retained topology from a pattern
 *        whose census matches but whose connectivity does not.
 * @details Cube- and octahedron-seeded patterns agree on face, index and vertex
 *          counts, so only the topology key separates them.
 */
inline void case_update_hankin_dual_seed_topology() {
  static uint8_t geom_buf[192 * 1024];
  static uint8_t scratch_buf[128 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  CompiledHankin cube_pattern;
  MeshOps::compile_hankin(cube, cube_pattern, geom, scratch);

  MeshState mesh;
  MeshOps::update_hankin(cube_pattern, mesh, geom, opaque(0.0f));
  MeshOps::classify_faces_by_topology(mesh, scratch, scratch, geom);

  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, geom);
  CompiledHankin octa_pattern;
  MeshOps::compile_hankin(octa, octa_pattern, geom, scratch);
  MeshOps::update_hankin(octa_pattern, mesh, geom, opaque(0.0f));
  if (mesh.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: update_hankin rejects a borrowed-mode topology from a
 *        different compiled pattern.
 * @details A borrowed MeshState reports its topology through the view, not the
 *          owned array, and update_hankin drops that view on entry — so the
 *          reuse check has to sample the size before the drop or a borrowed
 *          mesh walks past it.
 */
inline void case_update_hankin_borrowed_stale_topology() {
  static uint8_t geom_buf[192 * 1024];
  static uint8_t scratch_buf[128 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  CompiledHankin cube_pattern;
  MeshOps::compile_hankin(cube, cube_pattern, geom, scratch);

  MeshState source;
  MeshOps::update_hankin(cube_pattern, source, geom, opaque(0.0f));
  MeshOps::classify_faces_by_topology(source, scratch, scratch, geom);

  MeshState mesh;
  mesh.set_borrowed(ArenaSpan<uint8_t>(source.face_counts),
                    ArenaSpan<uint16_t>(source.faces),
                    ArenaSpan<uint16_t>(source.face_offsets),
                    ArenaSpan<uint16_t>(source.topology), source.topology_key);

  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, geom);
  CompiledHankin octa_pattern;
  MeshOps::compile_hankin(octa, octa_pattern, geom, scratch);
  MeshOps::update_hankin(octa_pattern, mesh, geom, opaque(0.0f));
  if (mesh.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: update_hankin rejects a non-finite contact angle.
 * @details Hankin surface — the half-angle sine and cosine carry a NaN into
 *          every star point, and normalized_or's dot(v, v) < EPS guard is
 *          false for NaN, so the whole pattern mesh reaches the rasterizer.
 */
inline void case_update_hankin_nonfinite_angle() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  CompiledHankin compiled;
  MeshState mesh;
  MeshOps::update_hankin(compiled, mesh, arena,
                         opaque(std::numeric_limits<float>::quiet_NaN()));
}

/**
 * @brief Death case: reading a ParamDef with an unknown target type must trap.
 * @details An unknown tag has no supported value representation and traps
 *          before the descriptor reads the target.
 */
inline void case_param_def_unknown_get_target_type() {
  float storage = 0.5f;
  ParamDef def;
  def.target = &storage;
  def.target_type = static_cast<ParamDef::TargetType>(opaque<uint8_t>(9));
  if (def.get_from(&storage) == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: writing a ParamDef with an unknown target type must trap.
 * @details An unknown tag has no supported value representation and traps
 *          before the descriptor writes the target.
 */
inline void case_param_def_unknown_set_target_type() {
  float storage = 0.5f;
  ParamDef def;
  def.target = &storage;
  def.target_type = static_cast<ParamDef::TargetType>(opaque<uint8_t>(9));
  struct InternalWriter : ParamHost {
    using ParamHost::write_parameter_unchecked;
  };
  InternalWriter::write_parameter_unchecked(def, opaque(1.0f));
  if (storage == opaque(42.0f))
    std::printf("x");
}
