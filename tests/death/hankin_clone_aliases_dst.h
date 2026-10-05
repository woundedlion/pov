/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

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
 * @details Mesh-borrow surface — the accessors index offsets by face, so an
 *          offsets array that is not one entry per face would read past its end
 *          on the solid scan path; set_borrowed rejects it at the install site.
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
 * @details Mesh-borrow surface — the key names the connectivity a topology was
 *          classified for, so a key with no span behind it would hand a
 *          downstream reuse check a classification the mesh does not carry.
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
 * @brief Builds a PolyMesh from an explicit face-count and flat index list.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param num_verts Vertex count; positions are all zero (never read).
 * @param counts Per-face side counts.
 * @param num_faces Number of entries in @p counts.
 * @param indices Flat per-face vertex index list.
 * @param num_indices Number of entries in @p indices.
 */
inline void build_polymesh(PolyMesh &mesh, Arena &arena, size_t num_verts,
                           const uint8_t *counts, size_t num_faces,
                           const uint16_t *indices, size_t num_indices) {
  mesh.vertices.bind(arena, num_verts);
  for (size_t i = 0; i < num_verts; ++i)
    mesh.vertices.push_back(math::Vector{});
  mesh.face_counts.bind(arena, num_faces);
  for (size_t i = 0; i < num_faces; ++i)
    mesh.face_counts.push_back(opaque(counts[i]));
  mesh.faces.bind(arena, num_indices);
  for (size_t i = 0; i < num_indices; ++i)
    mesh.faces.push_back(opaque(indices[i]));
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
 *          Pairing the first two would leave the third silently unpaired.
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
 *          same direction, so the undirected pairing key matches and they would
 *          otherwise pair into a mesh that passes require_closed_manifold while
 *          every vertex_orbit walk through the pair runs backwards.
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
 * @details Mesh-topology surface — every operator narrows its output valence
 *          through this shared guard, so a high-valence orbit traps instead of
 *          wrapping the uint8_t face_counts entry.
 */
inline void case_mesh_narrow_face_count() {
  uint8_t c = MeshOps::narrow_face_count(opaque(UINT8_MAX + 1)); // -> HS_CHECK
  if (c == 0xEE)
    std::printf("x");
}

/**
 * @brief Death case: an open mesh must trap the closed-manifold requirement.
 * @details Mesh-topology surface — operators size their output pools from
 *          E = I/2, so a lone triangle's three unpaired half-edges are rejected
 *          up front instead of overrunning a pool far from the cause.
 */
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
 *          and edge-manifold, so only the fan pass catches them; the orbit
 *          scaffolding would otherwise emit one face from the first fan and
 *          silently drop the second.
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

inline void case_noise_hue_bake_invalid_scale() {
  static std::array<int8_t, HueNoiseLutView::SIZE> output{};
  FastNoiseLite noise;
  HueNoiseBakeCache cache;
  (void)cache.refresh(output, noise, opaque(0.0f), 0.0f);
}

inline void run_invalid_recipe_step(const Solids::OpStep &step) {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  (void)Solids::build_steps(opaque<uint8_t>(1), &step, 1, a, b);
}

inline void case_recipe_bake_wrong_op() {
  const MeshOps::RelaxBake bake{};
  run_invalid_recipe_step({Solids::Op::AMBO, 0.0f, 0.0f, &bake});
}

inline void case_recipe_twist_wrong_op() {
  run_invalid_recipe_step({Solids::Op::AMBO, 0.0f, opaque(0.1f)});
}

inline void case_recipe_bake_live_iterations() {
  const MeshOps::RelaxBake bake{};
  run_invalid_recipe_step({Solids::Op::RELAX, opaque(1.0f), 0.0f, &bake});
}

inline void case_pullback_mobius_degenerate() {
  Pullback::Interp::Op::MobiusChainParams params;
  params.a_re = opaque(0.0f);
  params.d_re = opaque(0.0f);
  Pullback::Interp::Op::LensMobius::State state;
  Pullback::Interp::FrameContext context{};
  (void)Pullback::Interp::Op::LensMobius::prepare(context, params, state);
}

inline void case_lattice_shells_oob() {
  SDF::Lattice::Settings settings;
  settings.shells = static_cast<SDF::Lattice::ShellCount>(3);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: zero lattice softness would divide by zero in shading. */
inline void case_lattice_zero_softness() {
  SDF::Lattice::Settings settings;
  settings.softness = opaque(0.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: zero lattice cell size has no inverse transform. */
inline void case_lattice_zero_cell_size() {
  SDF::Lattice::Settings settings;
  settings.cell_size = opaque(0.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: a negative AA strength gives invalid crossing widths. */
inline void case_lattice_negative_aa() {
  SDF::Lattice::Settings settings;
  settings.aa_strength = opaque(-1.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: a HyperLattice frame without crossing scratch traps. */
inline void case_hyperlattice_frame_without_crossings() {
  const HyperLatticeDetail::FrameState frame{};
  (void)HyperLatticeDetail::prepare_trace(frame);
}

inline void case_mindsplatter_profile_preset_oob() {
  MindSplatter<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(SIZE_MAX));
}

/**
 * @brief Death case: a HANKIN step with no contact angle must trap.
 * @details Recipe-replay surface — the zero default collapses every star point
 *          onto its corner, so an authored step that forgot its angle replays
 *          as a flat tiling instead of failing.
 */
inline void case_apply_step_hankin_no_angle() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::HANKIN, opaque(0.0f)}};
  PolyMesh mesh = Solids::build_steps(opaque<uint8_t>(1), steps, 1, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a BEVEL step with no depth must trap.
 * @details Recipe-replay surface — the composite lowers to ambo, truncate(t),
 *          so a zero default is the depthless truncate the lowered replay
 *          already traps on; the authored replay must not diverge from it.
 */
inline void case_apply_step_bevel_no_depth() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::BEVEL, opaque(0.0f)}};
  const Solids::Recipe recipe = {opaque<uint8_t>(1), steps, 1};
  PolyMesh mesh = Solids::build_recipe(recipe, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a NaN endpoint fed to slerp must trap.
 * @details Math-core surface — the NaN poisons interpolation through both
 *          branches into the final strict normalized(), which traps rather than
 *          emitting a NaN direction into geometry. Proves the non-finite input
 *          is caught at the slerp seam, not just at bare normalize().
 */
inline void case_slerp_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector dst{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  math::Vector v =
      math::slerp(bad, dst, opaque(0.5f)); // NaN -> normalized() -> HS_CHECK
  if (v.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a NaN source must trap.
 * @details A NaN component fails the unit-vector precondition before rotation
 *          arithmetic, complementing the finite non-unit input case.
 */
inline void case_make_rotation_vectors_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector from{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector to{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  math::Quaternion q = math::make_rotation(from, to);
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(axis, theta) with a NaN angle must trap.
 * @details Math-core surface — cos/sin of a NaN poison the quaternion, and its
 *          normalized() traps on the NaN magnitude.
 */
inline void case_make_rotation_angle_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector axis{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Quaternion q =
      math::make_rotation(axis, nan); // NaN quat -> normalized() -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a NaN normal must trap.
 * @details Geometry surface — rotate(normal,.).normalized() is the first strict
 *          normalize in the basis construction and traps on the NaN-poisoned
 *          vector rather than returning a garbage frame.
 */
inline void case_make_basis_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  math::Vector normal{nan, opaque(0.0f), opaque(0.0f)};
  math::Basis b = math::make_basis(math::Quaternion(),
                                   normal); // NaN -> normalized() -> HS_CHECK
  if (b.u.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an active noise_transform fed a non-finite direction must
 *        trap, not propagate NaN/Inf into the rendered geometry.
 * @details The structural audit rejects non-finite directions before noise sampling.
 */
inline void case_noise_transform_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Animation::NoiseParams p;
  p.amplitude = opaque(0.5f); // active path (skips the zero-amplitude no-op)
  p.scale = opaque(4.0f);
  p.time = opaque(1.0f);
  math::Vector v{nan, opaque(0.0f), opaque(0.0f)};
  math::Vector r = noise_transform(v, p);
  if (r.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a non-unit source must trap.
 * @details Math-core surface — the d-based parallel/antiparallel branches assume
 *          |from| = |to| = 1, so a finite but non-unit input must trap at the
 *          unit-vector guard rather than silently skewing the rotation angle.
 */
inline void case_make_rotation_nonunit() {
  math::Vector from{opaque(2.0f), opaque(0.0f), opaque(0.0f)}; // |from| = 2
  math::Vector to{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Quaternion q = math::make_rotation(from, to); // |from| != 1 -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a live-source Driver built with a null speed pointer must trap.
 * @details Animation surface — the guard traps rather than dereferencing the
 *          null pointer in the member-init list.
 */
inline void case_driver_null_speed_src() {
  static float mutant = 0.0f;
  Animation::Driver d(mutant, opaque<const float *>(nullptr),
                      1.0f); // -> HS_CHECK
  (void)d;
  if (mutant == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a second TransformerPool::init_storage() must trap.
 * @details Transformer surface — a re-init would hand the pool a second block
 *          while spawned animations still hold Params references into the first,
 *          and would silently double-charge the persistent arena.
 */
inline void case_transformer_pool_init_storage_twice() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.init_storage(persistent_arena); // entities already set -> HS_CHECK
}

/** @brief Unpinned perpetual spawns trap even when the timeline is full. */
inline void case_transformer_unpinned_full() {
  configure_arenas_default();
  Timeline timeline;
  NoiseTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    timeline.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  transformer.spawn(0);
}

/**
 * @brief Death case: spawning before init_storage() must trap.
 * @details Transformer surface — the slot scan indexes the entity block, so a
 *          spawn on an un-initialized pool would dereference null instead of
 *          reporting the missed init() wiring.
 */
inline void case_transformer_pool_spawn_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  Animation::Ripple *p =
      rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: preparing frame state before init_storage() must trap.
 * @details Transformer surface — prepare_frame() is the ordering contract's
 *          other half: an un-initialized pool has no active slots, so it would
 *          silently do nothing and leave the composition reading state that was
 *          never prepared, instead of reporting the missed init() wiring.
 */
inline void case_transformer_pool_prepare_frame_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.prepare_frame(); // -> HS_CHECK
}

/**
 * @brief Death case: a pausable spawn with no pause flag must trap.
 * @details Transformer surface -- spawn_pausable() exists only to hand the
 *          animation a gate to read every frame; a null flag would schedule an
 *          event that can never pause, which is what plain spawn() is for.
 */
inline void case_transformer_pool_spawn_pausable_null_flag() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  Animation::Ripple *p = rt.spawn_pausable(nullptr, 0, math::Vector(0, 1, 0),
                                           0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: an out-of-range active index must trap.
 * @details Transformer surface — active_params() indexes the compact active list,
 *          which is shorter than CAPACITY, so an index taken from the slot domain
 *          (or from a stale count) would read a dead slot as if it were live.
 */
inline void case_transformer_pool_active_index_oob() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  const Animation::RippleParams &p = rt.active_params(opaque(0)); // -> HS_CHECK
  if (p.amplitude == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: reclaimed storage landing at a new address must trap.
 * @details Transformer surface — spawned animations hold Params references into
 *          the slots, so the post-reset replay must re-land the blocks exactly
 *          where init_storage() put them. Here the arena is NOT reset first, so
 *          the replay appends past the originals and every live reference would
 *          be left pointing at abandoned bytes.
 */
inline void case_transformer_pool_reclaim_storage_moved() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.reclaim_storage(persistent_arena); // blocks land elsewhere -> HS_CHECK
}

/**
 * @brief Death case: spawning after the pool's arena was reclaimed must trap.
 * @details Transformer surface — init_storage() must run after
 *          configure_arenas(), which rebinds the arena and hands its bytes out
 *          again. The slot pointers stay non-null across that, so the watermark
 *          is what catches the ordering, in every build.
 */
inline void case_transformer_pool_arena_reclaimed() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  configure_arenas_default(); // rebinds under the live pool
  Animation::Ripple *p =
      rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/** @brief Rejects a margin below the filter pipeline requirement. */
inline void case_effect_margin_below_pipeline() {
  struct MarginEffect : Effect {
    MarginEffect() : Effect(32, 16, {.margin = 3, .required_margin = 3}) {}
    void draw_frame() override {}
  } effect;
  effect.set_margin(2);
}

inline void case_transformer_pinned_owner_order() {
  configure_arenas_default();
  Timeline timeline;
  using Pool = NoiseTransformer<1>;
  alignas(Pool) static uint8_t storage[sizeof(Pool)];
  Pool *first = new (storage) Pool(timeline);
  Pool second(timeline);
  first->init_storage(persistent_arena);
  second.init_storage(persistent_arena);
  first->spawn_pinned(0);
  second.spawn_pinned(0);
  first->~Pool();
}
