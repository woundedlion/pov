/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Mesh death cases.

/**
 * @brief Death case: CompiledHankin::clone rejects a self-aliased destination.
 * @details Each vector is rebound from the arena before the copy, so a
 *          self-clone memcpy's a block onto itself from a stale source pointer.
 */
inline void case_hankin_clone_aliases_dst() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  CompiledHankin compiled;
  CompiledHankin::clone(compiled, compiled, arena);
}

/**
 * @brief Death case: a face-offsets span with the wrong length must trap.
 * @details The accessors index offsets by face, so an offsets array that is
 *          not one entry per face would read past its end.
 */
inline void case_mesh_state_set_borrowed_offsets_count_mismatch() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 1);
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 3);
  for (uint16_t i = 0; i < 3; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 2);
  offsets.push_back(opaque<uint16_t>(0));
  offsets.push_back(opaque<uint16_t>(3));
  MeshState m;
  // 2 offsets for 1 face -> HS_CHECK
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: face offsets that do not span the flat faces list must
 *        trap.
 * @details Mesh-borrow surface — the last offset plus that face's count must
 *          reach the end of the flat list, or a walk of the final face reads
 *          short of the data the view claims to cover.
 */
inline void case_mesh_state_set_borrowed_offsets_short_span() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 1);
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 4);
  for (uint16_t i = 0; i < 4; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 1);
  offsets.push_back(opaque<uint16_t>(0));
  MeshState m;
  // 0 + 3 != 4 -> HS_CHECK
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: an empty topology span carrying a non-zero key must trap.
 * @details The key names the connectivity a topology was classified for.
 */
inline void case_mesh_state_set_borrowed_keyed_empty_topology() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 1);
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 3);
  for (uint16_t i = 0; i < 3; ++i)
    faces.push_back(opaque(i));
  MeshState m;
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces), {}, {},
                 opaque<uint32_t>(0x1234u));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: face offsets that are not the counts' prefix sum must trap.
 * @details Mesh-borrow surface — the count and span checks pass on the endpoints
 *          alone, so an interior offset off the prefix sum would walk one face
 *          over another's indices. The audit walk catches it.
 */
inline void case_mesh_state_set_borrowed_offsets_not_prefix_sum() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 2);
  counts.push_back(opaque<uint8_t>(3));
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 6);
  for (uint16_t i = 0; i < 6; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 2);
  offsets.push_back(opaque<uint16_t>(1)); // prefix sum starts at 0
  offsets.push_back(opaque<uint16_t>(3));
  MeshState m;
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a zero-side face must trap while building half-edges.
 * @details Mesh-topology surface — a zero-count face emits no half-edges yet
 *          still claims a face slot, whose half_edge entry would then point at
 *          the next face's loop. The trailing triangle keeps the flat index
 *          list non-empty so the pairing scratch is a real allocation.
 */
inline void case_half_edge_zero_side_face() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {0, 3};
  const uint16_t indices[] = {0, 1, 2};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 3, counts, 2, indices, 3);
  HalfEdgeMesh half_edges(arena, mesh); // face 0 has zero sides -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: >2 half-edges on one undirected edge must trap.
 * @details Mesh-topology surface — three faces share edge (0,1), so the pairing
 *          pass sees a run of three where a 2-manifold allows at most two.
 */
inline void case_half_edge_non_manifold_edge() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 1, 3, 0, 1, 4};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 5, counts, 3, indices, 9);
  HalfEdgeMesh half_edges(arena, mesh); // 3 half-edges on (0,1) -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: two faces wound the same way around a shared edge must
 *        trap.
 * @details Mesh-topology surface — both triangles traverse edge (0,1) in the
 *          same direction, so the undirected pairing key matches.
 */
inline void case_half_edge_inconsistent_winding() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 1, 3};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 4, counts, 2, indices, 6);
  HalfEdgeMesh half_edges(arena, mesh); // both edges run 0->1 -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: a face side count past uint8_t must trap.
 * @details Traps instead of wrapping the uint8_t face_counts entry.
 */
inline void case_mesh_narrow_face_count() {
  uint8_t c = MeshOps::narrow_face_count(opaque(UINT8_MAX + 1)); // -> HS_CHECK
  if (c == 0xEE)
    std::printf("x");
}

/** @brief Death case: an open mesh must trap the closed-manifold requirement. */
inline void case_mesh_require_closed_manifold() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3};
  const uint16_t indices[] = {0, 1, 2};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 3, counts, 1, indices, 3);
  HalfEdgeMesh half_edges(arena, mesh);
  // unpaired -> HS_CHECK
  MeshOps::require_closed_manifold(half_edges, arena, "death");
}

/**
 * @brief Death case: a bowtie vertex must trap the closed-manifold requirement.
 * @details Mesh-topology surface — two tetrahedra joined at vertex 0 are closed
 *          and edge-manifold, so only the fan pass catches them.
 */
inline void case_mesh_require_vertex_manifold() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3, 3, 3, 3, 3, 3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2,
                              0, 4, 5, 0, 5, 6, 0, 6, 4, 4, 6, 5};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 7, counts, 8, indices, 24);
  HalfEdgeMesh half_edges(arena, mesh);
  // split fan at vertex 0 -> HS_CHECK
  MeshOps::require_closed_manifold(half_edges, arena, "death");
}

/**
 * @brief Death case: a half-edge mesh whose faces have different side counts
 *        must trap even when the census matches.
 * @details Mesh-topology surface — the reuse overloads walk face loops from the
 *          half-edge mesh while sizing and indexing from the source mesh, so a
 *          {3,4} pairing against a {4,3} source shares (V,F,I) yet emits every
 *          face from the wrong span.
 */
inline void case_mesh_require_matching_face_sides() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t he_counts[] = {3, 4};
  const uint8_t mesh_counts[] = {4, 3};
  const uint16_t indices[] = {0, 1, 2, 3, 4, 5, 6};
  PolyMesh he_source;
  build_polymesh(he_source, arena, 7, he_counts, 2, indices, 7);
  HalfEdgeMesh half_edges(arena, he_source);
  PolyMesh mesh;
  build_polymesh(mesh, arena, 7, mesh_counts, 2, indices, 7);
  // face 1 starts at 3 in he_source, 4 in mesh -> HS_CHECK
  MeshOps::require_matching_half_edges(half_edges, mesh, "death");
}

/**
 * @brief Death case: side counts that outrun the flat index list must trap.
 * @details Mesh-topology surface -- the entry census compares the half-edge
 *          count against the flat index length, so a source whose own side
 *          counts sum past that length clears it and then indexes past the
 *          end of every later face.
 */
inline void case_mesh_require_matching_half_edge_census() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t he_counts[] = {3, 3};
  const uint8_t mesh_counts[] = {3, 4};
  const uint16_t indices[] = {0, 1, 2, 3, 4, 5};
  PolyMesh he_source;
  build_polymesh(he_source, arena, 7, he_counts, 2, indices, 6);
  HalfEdgeMesh half_edges(arena, he_source);
  PolyMesh mesh;
  // counts sum to 7 over a 6-entry index list -> HS_CHECK
  build_polymesh(mesh, arena, 7, mesh_counts, 2, indices, 6);
  MeshOps::require_matching_half_edges(half_edges, mesh, "death");
}

/**
 * @brief Death case: loops naming a different mesh's faces must trap.
 * @details Mesh-topology surface -- census, side counts and vertex range all
 *          match between two meshes over the same vertex set, so only the loop
 *          head vertices distinguish a connectivity built from the wrong mesh.
 */
inline void case_mesh_require_matching_half_edge_loops() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3};
  const uint16_t he_indices[] = {0, 1, 2, 3, 4, 5};
  const uint16_t mesh_indices[] = {0, 2, 1, 3, 4, 5};
  PolyMesh he_source;
  build_polymesh(he_source, arena, 6, counts, 2, he_indices, 6);
  HalfEdgeMesh half_edges(arena, he_source);
  PolyMesh mesh;
  // face 0 winds 0,2,1 in mesh and 0,1,2 in he_source -> HS_CHECK
  build_polymesh(mesh, arena, 6, counts, 2, mesh_indices, 6);
  MeshOps::require_matching_half_edges(half_edges, mesh, "death");
}

/** @brief Death case: reconcile endpoints with different sizes must trap. */
inline void case_reconcile_vertices_size_mismatch() {
  static uint8_t identity_buf[64];
  static uint8_t target_buf[64];
  static uint8_t scratch_buf[64];
  Arena identity_arena(identity_buf, sizeof(identity_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh identity;
  identity.vertices.bind(identity_arena, 1);
  identity.vertices.push_back(math::Vector{});
  PolyMesh authored;
  PolyMesh out;
  MeshOps::reconcile_vertices(identity, authored, out, target, scratch);
}

/** @brief Rejects output aliasing a reconciliation input. */
inline void case_reconcile_vertices_aliased_output() {
  static uint8_t target_bytes[128], scratch_bytes[128];
  Arena target(target_bytes, sizeof(target_bytes));
  Arena scratch(scratch_bytes, sizeof(scratch_bytes));
  PolyMesh identity, authored;
  MeshOps::reconcile_vertices(identity, authored, identity, target, scratch);
}

/** @brief Death case: reconciling a vertex-less endpoint pair must trap. */
inline void case_reconcile_vertices_empty() {
  static uint8_t target_buf[64];
  static uint8_t scratch_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh identity;
  PolyMesh authored;
  PolyMesh out;
  MeshOps::reconcile_vertices(identity, authored, out, target, scratch);
}

/** @brief Death case: narrowing an index past the int16 topology range must trap. */
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

/**
 * @brief Death case: relax_baked rejects a bake whose vertex count differs from
 *        the source mesh.
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
 * @details Dimensions alone do not pin connectivity.
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
 * @details Covers the baked vertex words themselves.
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
 * @details The topology array survives an angle re-solve, so a mesh pointed at
 *          a new pattern would carry stale class ids.
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
 * @details A borrowed MeshState reports its topology through the view, which
 *          update_hankin drops on entry.
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
 * @details normalized_or's dot(v, v) < EPS guard is false for NaN.
 */
inline void case_update_hankin_nonfinite_angle() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  CompiledHankin compiled;
  MeshState mesh;
  MeshOps::update_hankin(compiled, mesh, arena,
                         opaque(std::numeric_limits<float>::quiet_NaN()));
}

/** @brief Reconciliation cannot retain output in its scratch arena. */
inline void case_reconcile_aliased_arenas() {
  static uint8_t storage[64];
  Arena arena(storage, sizeof(storage));
  PolyMesh identity, authored, out;
  MeshOps::reconcile_vertices(identity, authored, out, arena, arena);
}
