/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Shared death-test fixtures.

/**
 * @brief Recursive generator body: each level calls generate() once more.
 * @param target Arena forwarded to the next generate().
 * @param remaining Levels of nesting still to open; stops at 0.
 * @return Always 0.
 */
inline int nested_generate(Arena &target, Arena &, Arena &, int remaining) {
  if (remaining <= 0)
    return 0;
  return hs::generate(target, nested_generate, remaining - 1);
}

/**
 * @brief Trivial Cloneable whose clone() allocates from the destination arena,
 *        so a Persist restore measurably grows the persistent arena.
 */
struct PersistProbe {
  uint8_t *storage = nullptr; /**< Stand-in for arena-backed object state. */
  /**
   * @brief Clones by allocating fresh storage from @p arena.
   * @param src Source probe (unused beyond the Cloneable interface).
   * @param dst Destination probe receiving freshly allocated storage.
   * @param arena Arena the clone allocates from.
   */
  static void clone(const PersistProbe &src, PersistProbe &dst, Arena &arena) {
    (void)src;
    dst.storage = static_cast<uint8_t *>(arena.allocate(opaque<size_t>(8)));
  }
};

/**
 * @brief Cloneable payload that allocates nothing, so only the distinct-arena
 *        guard can trap a same-arena Persist.
 */
struct FlatProbe {
  int value = 0; /**< Whole payload; clone() copies it without an arena. */
  /**
   * @brief Clones by plain copy, leaving the arena offset untouched.
   * @param src Source probe.
   * @param dst Destination probe.
   * @param arena Unused; the payload needs no storage.
   */
  static void clone(const FlatProbe &src, FlatProbe &dst, Arena &arena) {
    (void)arena;
    dst.value = src.value;
  }
};

/** @brief Draw callback for the OpLeg construction death cases; never runs. */
inline void death_opleg_draw(Canvas &, MeshState &,
                             const Animation::OpLeg::Shading &) {}

/**
 * @brief A palette handoff complete enough to clear the OpLeg handoff guard.
 * @return A handoff naming a default bank and a one-face departed palette.
 * @details Lets a case reach a guard that sits behind the handoff check.
 */
inline Animation::OpLeg::PaletteHandoff death_opleg_handoff() {
  static const BakedPaletteBank bank;
  static const uint8_t face_palette[1] = {0};
  return {.bank = &bank,
          .prev_face_palette = face_palette,
          .prev_faces = 1,
          .prev_face_centroid = nullptr,
          .correspondence = Animation::OpLeg::FaceCorrespondence::GEOMETRIC};
}

/** @brief A well-formed non-settling graph edge for the OpLeg death cases. */
inline constexpr ConwayGraph::EdgeSpec death_opleg_edge{
    .from_node = 0,
    .to_node = 0,
    .seed_solid = 0,
    .op = ConwayGraph::MorphOp::TRUNCATE,
    .t_from = 0.0f,
    .t_to = 0.4f,
    .twist_from = 0.0f,
    .twist_to = 0.0f,
    .settle = false,
    .reseed = ConwayGraph::Reseed::NONE,
    .bridge = false};

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

inline void run_invalid_recipe_step(const Solids::OpStep &step) {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  (void)Solids::build_steps(opaque<uint8_t>(1), &step, 1, a, b);
}

inline void chain_invalid_layout(unsigned variant) {
  using namespace Pullback::Interp;
  ChainProgram program;
  OperatorDescriptor descriptor = OPERATOR_TABLE[0];
  alignas(std::max_align_t) uint8_t block_a[512], block_b[512];
  size_t capacity = sizeof(block_a);
  switch (variant) {
  case 0:
    descriptor.runtime.param.align = 0;
    break;
  case 1:
    descriptor.runtime.prepared.align = 3;
    break;
  case 2:
    descriptor.runtime.state.align = 2 * alignof(std::max_align_t);
    break;
  case 3:
    descriptor.runtime.param.size = 0;
    break;
  case 4:
    descriptor.runtime.state.size = 3;
    break;
  case 5:
    capacity = std::numeric_limits<uint32_t>::max();
    break;
  }
  program.bind_storage(block_a, block_b, capacity, {&descriptor, 1});
}

/** @brief Adds an event from a clear hook, violating the hook contract. */
inline void add_event_from_clear_hook(void *ctx) {
  static float value = 0.0f;
  static_cast<Timeline *>(ctx)->add(
      0, Animation::Transition(value, 1.0f, 1, math::ease_linear));
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
 * @brief Concrete Effect for the canvas death cases.
 * @details Defaults to 32x16 and exposes register_param via reg.
 */
struct DeathEffect : public Effect {
  using Effect::register_param;
  /**
   * @brief Constructs the effect at the requested resolution.
   * @param width Canvas width in pixels.
   * @param height Canvas height in pixels.
   */
  DeathEffect(int width = 32, int height = 16) : Effect(width, height) {}
  /**
   * @brief Draws one frame (no-op; the death cases never render).
   */
  void draw_frame() override {}
  /**
   * @brief Registers a parameter over the unit range, exposing register_param.
   * @param n Parameter name.
   * @param p Pointer to the backing float storage.
   */
  void reg(const char *n, float *p) { register_param(n, p, 0.0f, 1.0f); }
  void reg_float_options(float *p, const char *const *options, int count) {
    register_param("options", p, 0.0f, 1.0f, false, false, options, count);
  }
  /** @brief Registers a typed enum with the requested option count. */
  template <typename Enum> void reg_enum(Enum *p, int count) {
    static constexpr const char *OPTIONS[] = {"zero"};
    register_param("enum", p, OPTIONS, nullptr, count);
  }
  /**
   * @brief Registers an integer parameter, exposing register_int_param.
   * @param n Parameter name.
   * @param p Pointer to the backing integer storage.
   * @param min Minimum value, inclusive.
   * @param max Maximum value, inclusive.
   */
  template <typename Integer>
  void reg_int(const char *n, Integer *p, int min, int max) {
    register_int_param(n, p, min, max);
  }
};

/**
 * @brief Initializes a particle system with the requested maximum lifetime.
 * @param max_life Maximum lifetime passed to ParticleSystem::init.
 */
inline void init_particle_system_with_lifetime(float max_life) {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, max_life);
}

/** @brief Pipeline stub for the particle-render lifetime death case. */
struct DeathPlotPipeline {
  void plot(Canvas &, const math::Vector &, const Pixel &, float, float) {}
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

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

/**
 * @brief Minimal duck-typed mesh: one 2-gon face whose second index (130)
 *        exceeds the TriangularBitset<128> capacity.
 * @details The trap fires before any vertex or pipeline access, so the vertex
 *          store only needs to satisfy the interface.
 */
struct OverCapacityMockMesh {
  /**
   * @brief Stand-in vertex store satisfying the mesh interface.
   */
  struct Verts {
    /**
     * @brief Returns a fixed vertex for any index.
     * @return A constant Vector{0,1,0}.
     */
    math::Vector operator[](size_t) const {
      return math::Vector{0.0f, 1.0f, 0.0f};
    }
    /**
     * @brief Reports the vertex count.
     * @return Always 1.
     */
    size_t size() const { return 1; }
  } vertices;
  uint8_t fc[1];  /**< Face-counts data: a single 2-gon face. */
  uint16_t fi[2]; /**< Face-index data; second entry is over-capacity. */
  /**
   * @brief Builds the mock mesh with one over-capacity 2-gon face.
   * @details Stores the over-capacity index at runtime so the optimizer can't
   *          prove the trap at compile time and reshape the case (see opaque).
   */
  OverCapacityMockMesh() : fc{2}, fi{0, opaque<uint16_t>(130)} {}
  /**
   * @brief Returns the face-counts array.
   * @return Pointer to the face-counts data.
   */
  const uint8_t *get_face_counts_data() const { return fc; }
  /**
   * @brief Returns the number of faces.
   * @return Always 1.
   */
  size_t get_face_counts_size() const { return 1; }
  /**
   * @brief Returns the flat face-index array.
   * @return Pointer to the face-index data.
   */
  const uint16_t *get_faces_data() const { return fi; }
  /**
   * @brief Returns the flat face-index array length.
   * @return Always 2.
   */
  size_t get_faces_size() const { return 2; }
};
