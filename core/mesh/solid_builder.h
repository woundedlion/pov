/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/mesh/solid_generators.h.

#if defined(HS_RELAX_BAKE_VERIFY)
/** @brief Payloads relax_baked() has re-derived and matched this run. */
inline int relax_bakes_verified = 0;
#endif

/**
 * @brief Fluent builder for chaining Conway operators with automatic arena
 * swapping.
 * @details Each method runs `mesh = op(mesh, output_arena, scratch_arena)`,
 * swaps the two arenas, then rewinds the arena the next step writes into back
 * to its chain-start offset. The seed may sit in either arena; it sits below
 * both marks and stays for the life of the chain.
 */
class SolidBuilder {
  PolyMesh mesh; /**< Mesh being built; updated in place by each operator. */
  Arena *output_arena;  /**< Current output arena (swapped per op). */
  Arena *scratch_arena; /**< Current scratch arena (swapped per op). */
  size_t output_mark;   /**< output_arena's offset when the chain started. */
  size_t scratch_mark;  /**< scratch_arena's offset when the chain started. */

  /**
   * @brief Swaps the arena roles and reclaims the one the next step writes.
   * @details Every operator returns its output in `output_arena`, so once the
   * roles swap the live mesh is in `scratch_arena` and everything the other
   * arena holds above its start mark is a spent intermediate.
   */
  void advance() {
    std::swap(output_arena, scratch_arena);
    std::swap(output_mark, scratch_mark);
    output_arena->set_offset(output_mark);
  }

public:
  /**
   * @brief Constructs a builder seeded with an initial mesh and arena pair.
   * @param seed Starting mesh, moved into the builder.
   * @param a Initial output arena.
   * @param b Initial scratch arena.
   */
  HS_COLD_MEMBER SolidBuilder(PolyMesh seed, Arena &a, Arena &b)
      : mesh(std::move(seed)), output_arena(&a), scratch_arena(&b),
        output_mark(a.get_offset()), scratch_mark(b.get_offset()) {
    // One arena for both roles would put the live mesh in the arena advance()
    // rewinds.
    HS_CHECK(&a != &b, "SolidBuilder: output and scratch must differ");
  }

  /**
   * @brief Applies the dual operator (faces become vertices and vice versa).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &dual() {
    mesh = MeshOps::dual(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the kis operator (raise a pyramid on every face).
   * @return Reference to this builder for chaining.
   */
  HS_COLD_MEMBER SolidBuilder &kis() {
    mesh = MeshOps::kis(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the ambo operator (rectification: new vertex per edge).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &ambo() {
    mesh = MeshOps::ambo(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the truncate operator (cut corners off each vertex).
   * @param t Truncation depth in [0, 1] along each edge (the fraction at which
   *   each cut point sits). t == 0.5 short-circuits to ambo; t > 0.5 crosses
   *   the cuts into self-intersecting faces. See MeshOps::truncate.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &truncate(float t = MeshOps::TRUNCATE_DEFAULT_T) {
    mesh = MeshOps::truncate(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the expand operator (cantellation: push faces outward).
   * @param t Expansion amount in [0, 1); default places square faces at the canonical
   * gap.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &expand(float t = MeshOps::EXPAND_DEFAULT_T) {
    mesh = MeshOps::expand(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the chamfer operator (replace edges with hexagons).
   * @param t Fraction each face corner moves toward the face centroid, in
   *   [0, 1); at 1 a face collapses to its centroid, so it traps. See
   *   MeshOps::chamfer.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &chamfer(float t = MeshOps::CHAMFER_DEFAULT_T) {
    mesh = MeshOps::chamfer(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the snub operator (chiral expansion with a twist).
   * @param t Inset factor of each face toward its centroid, in [0, 1); at 1 a
   *   face collapses to its centroid, so it traps.
   * @param twist Per-face rotation about the face normal, in radians; 0
   *   disables the twist pass. See MeshOps::snub.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &snub(float t = MeshOps::SNUB_DEFAULT_T,
                     float twist = MeshOps::SNUB_DEFAULT_TWIST) {
    mesh = MeshOps::snub(mesh, *output_arena, *scratch_arena, t, twist);
    advance();
    return *this;
  }
  /**
   * @brief Applies the gyro operator (dual of snub; pentagonal chiral
   * subdivision).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &gyro() {
    mesh = MeshOps::gyro(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Relaxes vertex positions toward a regular configuration.
   * @param iterations Upper bound on the smoothing passes; relax stops early
   *   on convergence. The default is a light smoothing cap. Must be non-negative;
   *   0 is a normalize-only pass-through. See MeshOps::relax.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &relax(int iterations = MeshOps::RELAX_DEFAULT_ITERATIONS) {
    mesh = MeshOps::relax(mesh, *output_arena, *scratch_arena, iterations);
    advance();
    return *this;
  }
  /**
   * @brief Applies a host-generated relaxation payload.
   * @details The relaxed vertices load bit-identically on host and device.
   * Under HS_RELAX_BAKE_EXTRACT or HS_RELAX_BAKE_VERIFY the payload is instead
   * reproduced live with `bake.iterations` as the iteration cap: EXTRACT logs
   * the bits and measured guards, VERIFY asserts them against the committed
   * payload.
   * @param bake Payload whose guarded topology must match the current mesh.
   * @return Reference to this builder for chaining.
   */
  HS_COLD_MEMBER SolidBuilder &relax_baked(const MeshOps::RelaxBake &bake) {
#if defined(HS_RELAX_BAKE_EXTRACT) || defined(HS_RELAX_BAKE_VERIFY)
#if defined(HS_RELAX_BAKE_VERIFY)
    HS_CHECK(MeshOps::relax_topology_hash(mesh) == bake.topology_hash,
             "relax bake verify: source topology differs");
    MeshOps::check_relax_bake_source(mesh, bake);
#else
    const uint32_t topology_hash = MeshOps::relax_topology_hash(mesh);
    const uint32_t source_hash = MeshOps::relax_source_hash(mesh);
    const float source_margin = MeshOps::relax_source_quantization_margin(mesh);
#endif
    mesh = MeshOps::relax(mesh, *output_arena, *scratch_arena, bake.iterations);
    uint32_t output_hash = MeshOps::FNV1A_BASIS;
    for (const math::Vector &v : mesh.vertices) {
      output_hash = MeshOps::relax_output_hash(
          output_hash, std::bit_cast<uint32_t>(v.x),
          std::bit_cast<uint32_t>(v.y), std::bit_cast<uint32_t>(v.z));
    }
#if defined(HS_RELAX_BAKE_VERIFY)
    HS_CHECK(mesh.vertices.size() == bake.vertex_count &&
                 mesh.face_counts.size() == bake.face_count &&
                 mesh.faces.size() == bake.index_count,
             "relax bake verify: dimensions differ");
    HS_CHECK(output_hash == bake.output_hash,
             "relax bake verify: output differs");
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
      HS_CHECK(std::bit_cast<uint32_t>(mesh.vertices[i].x) ==
                       bake.vertex_bits[3 * i] &&
                   std::bit_cast<uint32_t>(mesh.vertices[i].y) ==
                       bake.vertex_bits[3 * i + 1] &&
                   std::bit_cast<uint32_t>(mesh.vertices[i].z) ==
                       bake.vertex_bits[3 * i + 2],
               "relax bake verify: vertex differs");
    }
    ++relax_bakes_verified;
#else // HS_RELAX_BAKE_EXTRACT: emit the payload for the generated header.
    hs::log("RELAX_BAKE_BEGIN %s %d %lu %lu %lu %08lx %08lx %08lx %d %08lx "
            "%08lx %08lx",
            bake.name, static_cast<int>(bake.iterations),
            static_cast<unsigned long>(mesh.vertices.size()),
            static_cast<unsigned long>(mesh.face_counts.size()),
            static_cast<unsigned long>(mesh.faces.size()),
            static_cast<unsigned long>(topology_hash),
            static_cast<unsigned long>(source_hash),
            static_cast<unsigned long>(output_hash),
            static_cast<int>(MeshOps::RELAX_SOURCE_SCALE),
            static_cast<unsigned long>(
                std::bit_cast<uint32_t>(MeshOps::RELAX_SOURCE_BIAS)),
            static_cast<unsigned long>(
                std::bit_cast<uint32_t>(MeshOps::RELAX_SOURCE_MIN_MARGIN)),
            static_cast<unsigned long>(std::bit_cast<uint32_t>(source_margin)));
    for (const math::Vector &v : mesh.vertices)
      hs::log("RELAX_BAKE_DATA %08lx %08lx %08lx",
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.x)),
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.y)),
              static_cast<unsigned long>(std::bit_cast<uint32_t>(v.z)));
    hs::log("RELAX_BAKE_END");
#endif
    advance();
    return *this;
#else
    mesh = MeshOps::relax_baked(mesh, *output_arena, bake);
    advance();
    return *this;
#endif
  }
  /**
   * @brief Applies the meta operator (kis composed with join).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &meta() {
    mesh = MeshOps::meta(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the needle operator (kis of the dual; n = kd).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &needle() {
    mesh = MeshOps::needle(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the zip operator (dual of kis).
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &zip() {
    mesh = MeshOps::zip(mesh, *output_arena, *scratch_arena);
    advance();
    return *this;
  }
  /**
   * @brief Applies the bevel operator (truncate composed with ambo).
   * @param t Truncation depth forwarded to the truncate step, in [0, 1]. At
   *   exactly 0.5 the chain is ambo(ambo). See MeshOps::bevel.
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &bevel(float t = MeshOps::BEVEL_DEFAULT_T) {
    mesh = MeshOps::bevel(mesh, *output_arena, *scratch_arena, t);
    advance();
    return *this;
  }
  /**
   * @brief Applies the Hankin star-pattern operator to each face.
   * @param angle Contact angle of the star pattern, in radians in [0, pi/2].
   * @return Reference to this builder for chaining.
   */
  SolidBuilder &hankin(float angle) {
    // *scratch_arena may hold mesh itself; hankin's ScratchScope marks above
    // the input, so compiling the topology into it leaves the input intact.
    mesh = MeshOps::hankin(mesh, *output_arena, *scratch_arena, angle);
    advance();
    return *this;
  }

  /**
   * @brief Finalizes the chain and yields the built mesh.
   * @return The accumulated PolyMesh, moved out of the builder.
   */
  HS_COLD_MEMBER PolyMesh build() { return std::move(mesh); }
};
