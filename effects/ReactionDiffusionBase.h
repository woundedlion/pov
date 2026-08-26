/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ReactionDiffusionBase.h
 * @brief CRTP base carrying the scaffolding the spherical reaction-diffusion
 *        effects share.
 */

#include "core/engine/engine.h"

/**
 * @brief CRTP base for the spherical reaction-diffusion effects (BZ, Gray-Scott).
 * @tparam Derived Concrete effect supplying the system-specific render().
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 * @details Both systems run on the same 7680-node Fibonacci-lattice K-NN graph,
 * share a Languid random-walk view orientation, build the cached node positions
 * once at init (the lattice is static), and interpolate with the same Wendland
 * C2 kernel, and seed their fields with the same saturated blobs. This base
 * captures exactly that shared scaffolding. The physics (Lotka-Volterra
 * 3-species vs Gray-Scott 2-species), params, palette, and rendering are
 * fundamentally different and stay in the derived classes.
 *
 * Dispatch is static (CRTP): draw_frame() forwards to Derived::render(); derived
 * classes befriend this base so render() can stay private.
 */
template <typename Derived, int W, int H>
class ReactionDiffusionBase : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N; /**< Lattice node count. */
  static constexpr int RD_K =
      ReactionGraph::RD_K; /**< K-NN neighbors per node. */

  /**
   * @brief Constructs the base with POV column strobing.
   * @details Leaves the framebuffer non-persistent: each frame is fully
   * repainted from lattice state.
   */
  HS_COLD_MEMBER ReactionDiffusionBase() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Advances one animation frame and dispatches to the derived renderer.
   * @details Advances the orientation timeline, then statically dispatches to
   * Derived::render() for the system-specific physics and drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(rd_timeline_step);
      timeline.step(canvas);
    }
    static_cast<Derived *>(this)->render(canvas);
  }

protected:
  /**
   * @brief Q16 full-scale factor: maps the [0, 65535] state to [0.0, 1.0].
   */
  static constexpr float Q16_SCALE = 65535.0f;
  static constexpr float Q16_INV =
      1.0f / Q16_SCALE; /**< Reciprocal of Q16_SCALE. */

  /**
   * @brief Converts a Q16 fixed-point sample to a normalized float.
   * @param v Q16 value in [0, 65535].
   * @return Concentration in [0.0, 1.0].
   */
  static inline float from_q16(uint16_t v) { return v * Q16_INV; }
  /**
   * @brief Converts a normalized concentration to a Q16 fixed-point sample.
   * @param v Concentration; clamped to [0.0, 1.0] before scaling.
   * @return Q16 value in [0, 65535], rounded to nearest.
   * @details Rounds to nearest (+0.5f); truncating would bias the dynamics down
   *          by dropping sub-LSB positive updates. clamp bounds the input so
   *          65535.5 -> 65535 with no overflow.
   */
  static __attribute__((always_inline)) inline uint16_t to_q16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 1.0f) * Q16_SCALE + 0.5f);
  }

  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG =
      ReactionGraph::D_AVG; /**< Mean inter-node spacing, sqrt(4π / RD_N). */
  static constexpr float KERNEL_R = 1.5f * D_AVG; /**< Kernel support radius. */
  static constexpr float INV_R2 =
      1.0f /
      (KERNEL_R * KERNEL_R); /**< Reciprocal of the squared support radius. */
  /**
   * @brief Total-weight floor below which a sampled kernel is treated as empty.
   * @details Sits far below the smallest legitimate total weight (the center
   * node alone contributes ~0.8 for any on-sphere query), so it only trips on a
   * degenerate/empty walk and never culls a real pixel.
   */
  static constexpr float KERNEL_MIN_TOTAL_WEIGHT = 1e-4f;

  /**
   * @brief Squared Euclidean distance between two points (no sqrt).
   * @param a First point.
   * @param b Second point.
   * @return Squared distance between `a` and `b`.
   */
  HS_O3_FN static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * @brief Invokes `on_weight(w)` with the Wendland C2 weight of one squared
   *        distance, skipping nodes outside the support radius.
   * @tparam OnWeight Callable accepting the kernel weight.
   * @param d2 Squared distance from the query to the node.
   * @param on_weight Callable invoked as `on_weight(weight)` only when the node
   * lies inside the support radius.
   * @details Passes the weight to a callback rather than returning it so the
   * out-of-support test stays a branch around the caller's accumulation.
   */
  template <typename OnWeight>
  static __attribute__((always_inline)) void
  with_wendland_weight(float d2, OnWeight &&on_weight) {
    float u = 1.0f - d2 * INV_R2;
    if (u > 0)
      on_weight(u * u);
  }

  /**
   * @brief Refines a cubemap-LUT seed by argmin over the seed and its one-ring.
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Node positions in the same frame as `rv`, indexed by node id.
   * @param seed Seed node id from the cubemap LUT.
   * @return The id of the closest node among the seed and its direct neighbors.
   * That is the true nearest node whenever the seed is that node or adjacent to
   * it, which is all the cubemap LUT promises; the two-hop equatorial seeds it
   * tolerates (tests/test_reaction_graph.h) stay unrecovered.
   * @details Off the render path: both systems center their stencils with
   * refine_render_center. This unconditional walk is the independent oracle
   * tests/test_effects.h measures that certified early-out against.
   */
  HS_O3_FN static int refine_center(const Vector &rv, const Vector *nodes,
                                    int seed) {
    float best_d2 = dist2(rv, nodes[seed]);
    int best = seed;
    for_each_neighbor(seed, [&](int ni) {
      float d = dist2(rv, nodes[ni]);
      if (d < best_d2) {
        best_d2 = d;
        best = ni;
      }
    });
    return best;
  }

  /** @brief Nodes at each pole packed tighter than the bulk minimum spacing. */
  static constexpr int POLE_BAND = 4;
  /**
   * @brief Minimum nearest-neighbor spacing outside the polar band, as a
   *        fraction of D_AVG; the lattice measures 0.9038.
   */
  static constexpr float BULK_MIN_SPACING_FRAC = 0.903f;
  /** @brief Half the bulk minimum spacing, squared. */
  static constexpr float BULK_CERTIFIED_D2 = 0.000333f;
  /** @brief Half the polar-band minimum spacing, squared. */
  static constexpr float POLE_CERTIFIED_D2 = 0.00012f;
  static_assert(BULK_CERTIFIED_D2 <= 0.25f * (BULK_MIN_SPACING_FRAC * D_AVG) *
                                         (BULK_MIN_SPACING_FRAC * D_AVG),
                "bulk certificate radius exceeds half the lattice's minimum "
                "bulk node spacing");
  // The two nodes closest to each pole sit 2/sqrt(RD_N - 1) apart, so a quarter
  // of that squared is the polar half-spacing bound.
  static_assert(POLE_CERTIFIED_D2 <= 1.0f / (RD_N - 1),
                "polar certificate radius exceeds half the pole-pair spacing");

  /**
   * @brief Refines a cubemap-LUT seed to the nearest node, skipping the
   *        six-neighbor walk when the seed is provably already nearest.
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Node positions in the same frame as `rv`, indexed by node id.
   * @param seed Seed node id from the cubemap LUT.
   * @return The id of the nearest node among the seed and its neighbors.
   * @details The certificates are properties of the shared lattice, so either
   * system may center its render stencil with this. A query lying within half
   * the distance from the seed to its closest lattice neighbor cannot be nearer
   * to any other node, so the seed is the argmin by triangle inequality. The
   * polar band carries its own, smaller certificate because the Fibonacci
   * lattice packs those nodes far tighter than D_AVG.
   */
  HS_O3_FN static int refine_render_center(const Vector &rv,
                                           const Vector *nodes, int seed) {
    float best_d2 = dist2(rv, nodes[seed]);
    float safe_d2 = seed >= POLE_BAND && seed < RD_N - POLE_BAND
                        ? BULK_CERTIFIED_D2
                        : POLE_CERTIFIED_D2;
    if (best_d2 < safe_d2)
      return seed;
    int best = seed;
    for_each_neighbor(seed, [&](int ni) {
      float d = dist2(rv, nodes[ni]);
      if (d < best_d2) {
        best_d2 = d;
        best = ni;
      }
    });
    return best;
  }

  /**
   * @brief Refines a cubemap-LUT seed to the nearest node and runs the
   *        Wendland C2 kernel walk in one stencil pass.
   * @tparam OnWeight Callable accepting (node_index, weight).
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Node positions in the same frame as `rv`, indexed by node id.
   * @param seed Seed node id from the cubemap LUT.
   * @param on_weight Callable invoked as `on_weight(node_index, weight)` for
   * every node inside the support radius.
   * @details Off the render path: its only caller is GSReactionDiffusion's
   * interpolate_b, itself a tests/test_effects.h oracle. It stands as the
   * per-sample reference the shared-stencil shaders are bounded against. The
   * seed stencil's squared distances are computed once while tracking the
   * argmin: when the seed is already nearest they feed the kernel weights
   * directly; otherwise the kernel re-walks the refined center's stencil.
   */
  template <typename OnWeight>
  HS_O3_FN static void refine_and_accumulate(const Vector &rv,
                                             const Vector *nodes, int seed,
                                             OnWeight &&on_weight) {
    float d2s[RD_K + 1];
    int ids[RD_K + 1];
    d2s[0] = dist2(rv, nodes[seed]);
    ids[0] = seed;
    int n = 1;
    int best = 0;
    for_each_neighbor(seed, [&](int ni) {
      float d = dist2(rv, nodes[ni]);
      ids[n] = ni;
      d2s[n] = d;
      if (d < d2s[best])
        best = n;
      ++n;
    });
    if (best == 0) {
      for (int i = 0; i < n; ++i)
        with_wendland_weight(d2s[i], [&](float w) { on_weight(ids[i], w); });
      return;
    }
    kernel_accumulate(rv, nodes, ids[best], on_weight);
  }

  /**
   * @brief Gathers a render stencil: the center node and its one-ring.
   * @tparam OnSlot Callable accepting (slot, node_index).
   * @param nodes Node positions in the query's frame, indexed by node id.
   * @param center Center node id, already refined.
   * @param positions Receives the RD_K + 1 stencil positions.
   * @param on_slot Callable invoked as `on_slot(slot, node_index)` once per
   * stencil slot, filling whatever per-node payload the caller interpolates.
   * @details Hoisted out of the sub-sample loop: every sub-sample of a pixel
   * re-weights this one stencil.
   * @note always_inline rather than HS_O3_FN: the mismatched optimize
   * attribute would keep GCC from inlining this into the -O3 shader body.
   */
  template <typename OnSlot>
  static __attribute__((always_inline)) void
  gather_stencil(const Vector *nodes, int center, Vector *positions,
                 OnSlot &&on_slot) {
    positions[0] = nodes[center];
    on_slot(0, center);
    int k = 1;
    for_each_neighbor(center, [&](int ni) {
      positions[k] = nodes[ni];
      on_slot(k, ni);
      ++k;
    });
  }

  /**
   * @brief Runs the Wendland C2 kernel over a gathered stencil.
   * @tparam OnWeight Callable accepting (slot, weight).
   * @param rv Sub-sample direction (unit vector on the sphere).
   * @param positions Stencil positions from gather_stencil.
   * @param on_weight Callable invoked as `on_weight(slot, weight)` for every
   * stencil slot inside the support radius.
   * @note always_inline is load-bearing on both this walk and the weight
   * callback it wraps: out of line, GCC spends extra ITCM in both shaders.
   */
  template <typename OnWeight>
  static __attribute__((always_inline)) void
  accumulate_stencil(const Vector &rv, const Vector *positions,
                     OnWeight &&on_weight) {
    for (int j = 0; j < RD_K + 1; ++j)
      with_wendland_weight(
          dist2(rv, positions[j]),
          [&](float w) __attribute__((always_inline)) { on_weight(j, w); });
  }

  /**
   * @brief Vertex-shader seed: tags a fragment with its cubemap-LUT node id.
   * @param frag Fragment whose pos seeds frag.v0 with the nearest-node id.
   * @details Shared by both systems' render() vertex shaders; the fragment
   * shader refines this face-quantized seed to the true nearest node (see
   * refine_render_center).
   */
  void seed_face_lut(Fragment &frag) {
    Vector rv = inverse_orientation.apply(frag.pos);
    frag.v0 = static_cast<float>(cube_lut.lookup(rv));
  }

private:
  /**
   * @brief Installs a Languid random-walk animation of the view orientation.
   * @details Shared by both systems; seeds OpenSimplex2 noise and adds the
   * random-walk animation to the timeline.
   */
  void init_orientation_animation() {
    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /**
   * @brief Builds the cached lattice node positions (RD_N entries) into `nodes`.
   * @param nodes Output array of at least RD_N positions to fill.
   * @details Called once at init; the Fibonacci lattice is static, not per-frame.
   */
  static void build_nodes(Vector *nodes) {
    for (int i = 0; i < RD_N; ++i)
      nodes[i] = ReactionGraph::node(i);
  }

  /**
   * @brief Rotates the lattice nodes into world space by one quaternion.
   * @param nodes Lattice-space node positions, count entries.
   * @param world Output world-space positions, count entries.
   * @param count Node count.
   * @param q Lattice-to-world rotation (the current view orientation).
   * @details Rotating the lattice once per frame lets the per-sub-sample kernel
   * walks compare world-space queries directly instead of un-orienting every
   * query.
   */
  HS_O3_FN static void orient_nodes(const Vector *nodes, Vector *world,
                                    int count, const Quaternion &q) {
    for (int i = 0; i < count; ++i)
      world[i] = rotate(nodes[i], q);
  }

protected:
  /**
   * @brief Splits the global arena for a reaction-diffusion effect, checking
   *        both halves against the device budget at compile time.
   * @tparam StateT Persistent per-node species sample type.
   * @tparam NSPECIES Persistent per-node species arrays.
   * @tparam PERSISTENT_BYTES Persistent arena size to carve.
   * @tparam EXTRA_PERSISTENT_BYTES Derived-only persistent tenants beyond the
   * cubemap LUT, the species state and the node array.
   * @tparam SCRATCH_PEAK_BYTES Largest of render()'s disjoint scratch phases.
   * @tparam SCRATCH_PEAK_TENANTS Separate allocations that phase makes.
   * @details The node-array term also bounds the equal-size transient lattice
   * cube_lut.build() carves and rewinds before init_lattice() allocates the
   * resident one, so the build() peak is covered by the same assert.
   */
  template <typename StateT, size_t NSPECIES, size_t PERSISTENT_BYTES,
            size_t EXTRA_PERSISTENT_BYTES, size_t SCRATCH_PEAK_BYTES,
            size_t SCRATCH_PEAK_TENANTS>
  HS_COLD_MEMBER static void configure_rd_arenas() {
    constexpr size_t CUBE_LUT_BYTES = 6u * ReactionGraph::CubemapLUT::RES *
                                      ReactionGraph::CubemapLUT::RES *
                                      sizeof(uint16_t);
    constexpr size_t STATE_BYTES = NSPECIES * RD_N * sizeof(StateT);
    constexpr size_t NODE_BYTES = RD_N * sizeof(Vector);
    // allocate() may skip up to alignof - 1 bytes ahead of each block; one pad
    // per persistent tenant (the species arrays, the LUT, the node array).
    constexpr size_t ALIGN_SLACK_BYTES =
        NSPECIES * alignof(StateT) + alignof(uint16_t) + alignof(Vector);
    static_assert(CUBE_LUT_BYTES + STATE_BYTES + NODE_BYTES +
                          ALIGN_SLACK_BYTES + EXTRA_PERSISTENT_BYTES <=
                      PERSISTENT_BYTES,
                  "RD persistent arena too small for LUT + state + build peak");
    // Same per-block pad as the persistent side, one per scratch tenant, at the
    // widest alignment allocate() can be asked for.
    constexpr size_t SCRATCH_SLACK_BYTES =
        SCRATCH_PEAK_TENANTS * (alignof(std::max_align_t) - 1);
    static_assert(SCRATCH_PEAK_BYTES + SCRATCH_SLACK_BYTES <=
                      DEVICE_GLOBAL_ARENA_SIZE - PERSISTENT_BYTES,
                  "RD scratch arena too small for render()'s phase peak");
    configure_arenas(PERSISTENT_BYTES, GLOBAL_ARENA_SIZE - PERSISTENT_BYTES, 0);
  }

  /**
   * @brief One frame's oriented lattice, bundled with the scratch scope it is
   *        carved from.
   * @details Ties the scope to the array the way init_lattice() ties the node
   * reserve to the fill: a caller cannot take the world nodes without taking
   * the reclaim. Nests LIFO like any other ScratchScope, so allocations made
   * from the same arena after the handle is constructed must not outlive it.
   */
  class OrientedLattice {
  public:
    /**
     * @brief Opens the scope, carves RD_N nodes, and fills them.
     * @param arena Scratch arena to carve from.
     * @param lattice Source lattice positions (RD_N entries).
     * @param q Lattice-to-world rotation.
     */
    OrientedLattice(Arena &arena, const Vector *lattice, const Quaternion &q)
        : scope(arena), world(static_cast<Vector *>(arena.allocate(
                            RD_N * sizeof(Vector), alignof(Vector)))) {
      orient_nodes(lattice, world, RD_N, q);
    }

    /**
     * @brief World-space node positions.
     * @return RD_N entries, live until this handle is destroyed.
     */
    Vector *get() const { return world; }

  private:
    ScratchScope scope; /**< Reclaims `world` on destruction. */
    Vector *world;      /**< Oriented positions carved from `scope`. */
  };

  /**
   * @brief Carves the per-frame oriented lattice out of scratch and fills it.
   * @return Handle owning both the scratch scope and the world-space nodes.
   */
  [[nodiscard]] OrientedLattice orient_lattice() {
    const Quaternion &current = orientation.get();
    inverse_orientation = RotationMatrix(current.conjugate());
    return OrientedLattice(scratch_arena_a, nodes, current);
  }

  /**
   * @brief Reserves and fills the shared lattice, then arms the view animation.
   * @details Bundles the three steps every derived `init()` must perform
   * together — reserve the RD_N node array, fill it with the static Fibonacci
   * lattice, and install the orientation random-walk — so a derived class cannot
   * perform a subset and silently ship a frozen or empty view. MUST be called
   * after configure_arenas() and after the derived class's own persistent
   * allocations, since the node array shares the persistent arena.
   */
  HS_COLD_MEMBER void init_lattice() {
    // configure_arenas() must have sized the persistent arena to hold the shared
    // node array before this runs; trap an unconfigured or under-sized arena by
    // contract here rather than as a later generic allocation OOM.
    HS_CHECK(persistent_arena.get_capacity() - persistent_arena.get_offset() >=
             RD_N * sizeof(Vector));
    // for_each_neighbor and the RD_K-degree Laplacian read every neighbor slot
    // unguarded; trap a deficient lattice here instead of on the hot path.
    for (int i = 0; i < RD_N; ++i)
      for (int k = 0; k < RD_K; ++k)
        HS_CHECK(ReactionGraph::neighbors[i][k] >= 0 &&
                 ReactionGraph::neighbors[i][k] < RD_N);
    nodes = static_cast<Vector *>(
        persistent_arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);
    init_orientation_animation();
  }

  /**
   * @brief Invokes `fn(ni)` for each of the RD_K neighbor indices of `node`.
   * @tparam Fn Callable accepting a neighbor node id.
   * @param node Center node id whose neighbors are visited.
   * @param fn Callable invoked once per neighbor index.
   * @details Reads all RD_K slots unguarded: the lattice is full-degree, with
   * every slot a valid node index (verified at init_lattice).
   */
  template <typename Fn> static void for_each_neighbor(int node, Fn &&fn) {
    for (int k = 0; k < RD_K; ++k)
      fn(ReactionGraph::neighbors[node][k]);
  }

  /**
   * @brief Saturates `clusters` blobs — a random node plus its neighbors — in
   *        one Q16 species field.
   * @param field Per-node Q16 samples of one species, RD_N entries.
   * @param clusters Blob count.
   * @note Draws exactly one hs::rand_int(0, RD_N) per cluster from the global
   *       deterministic RNG, so the blob count and the order fields are seeded
   *       in are both part of the global-determinism contract.
   */
  HS_COLD_MEMBER static void seed_blobs(uint16_t *field, int clusters) {
    for (int i = 0; i < clusters; ++i) {
      int center = hs::rand_int(0, RD_N);
      field[center] = 65535;
      for_each_neighbor(center, [&](int nb) { field[nb] = 65535; });
    }
  }

private:
  /**
   * @brief Wendland C2 kernel walk over `center` and its neighbors.
   * @tparam OnWeight Callable accepting (node_index, weight).
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Node positions in the same frame as `rv`, indexed by node id.
   * @param center Center node id of the kernel.
   * @param on_weight Callable invoked as `on_weight(node_index, weight)` for
   * every node inside the support radius.
   * @details The caller accumulates whatever fields it needs and applies the
   * total-weight guard, so this stays agnostic to species count and fixed-point.
   */
  template <typename OnWeight>
  HS_O3_FN static void kernel_accumulate(const Vector &rv, const Vector *nodes,
                                         int center, OnWeight &&on_weight) {
    auto visit = [&](int i) {
      with_wendland_weight(dist2(rv, nodes[i]),
                           [&](float w) { on_weight(i, w); });
    };
    visit(center);
    for_each_neighbor(center, visit);
  }

protected:
  Orientation<> orientation; /**< Current view orientation on the sphere. */
  RotationMatrix inverse_orientation{Quaternion()};
  FastNoiseLite noise; /**< Noise source driving the orientation walk. */
  ReactionGraph::CubemapLUT
      cube_lut;      /**< Cubemap LUT for fast nearest-node seeding. */
  Timeline timeline; /**< Animation timeline advancing the orientation. */
  Vector *nodes =
      nullptr; /**< Fixed Fibonacci-lattice node positions (RD_N), built once by init_lattice() and shared by both systems. */
};
