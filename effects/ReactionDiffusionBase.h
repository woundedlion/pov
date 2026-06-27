/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

#include <array>
#include <cstring>
#include <utility>

/**
 * @brief CRTP base for the spherical reaction-diffusion effects (BZ, Gray-Scott).
 * @tparam Derived Concrete effect supplying the system-specific render().
 * @tparam W Framebuffer width in pixels.
 * @tparam H Framebuffer height in pixels.
 * @details Both systems run on the same 7680-node Fibonacci-lattice K-NN graph,
 * share a Languid random-walk view orientation, build the cached node positions
 * once at init (the lattice is static), and interpolate with the same Wendland
 * C2 kernel. This base captures exactly that shared scaffolding. The physics
 * (Lotka-Volterra 3-species Q8 vs Gray-Scott 2-species Q16), state
 * representation, seeding, params, palette, and rendering are fundamentally
 * different and stay in the derived classes.
 *
 * Dispatch is static (CRTP): draw_frame() forwards to Derived::render() with no
 * virtual call, so the abstraction is zero-overhead. Derived classes befriend
 * this base so render() can stay private.
 */
template <typename Derived, int W, int H>
class ReactionDiffusionBase : public Effect {
public:
  static constexpr int RD_N = ReactionGraph::RD_N; /**< Lattice node count. */
  static constexpr int RD_K = ReactionGraph::RD_K; /**< K-NN neighbors per node. */

  /**
   * @brief Constructs the base, disabling pixel persistence.
   * @details Each frame is fully repainted from lattice state, so the
   * framebuffer need not persist between frames.
   */
  ReactionDiffusionBase() : Effect(W, H) { persist_pixels = false; }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Advances one animation frame and dispatches to the derived renderer.
   * @details Advances the orientation timeline, then statically dispatches to
   * Derived::render() for the system-specific physics and drawing.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    static_cast<Derived *>(this)->render(canvas);
  }

protected:
  // Wendland C2 compact kernel: w(d) = max(0, 1 - d²/R²)²
  static constexpr float D_AVG = ReactionGraph::D_AVG; /**< Mean inter-node spacing, sqrt(4π / RD_N). */
  static constexpr float KERNEL_R = 1.5f * D_AVG;      /**< Kernel support radius. */
  static constexpr float INV_R2 = 1.0f / (KERNEL_R * KERNEL_R); /**< Reciprocal of the squared support radius. */
  /**
   * @brief Total-weight floor below which a sampled kernel is treated as empty.
   * @details The center node alone contributes w = (1 - d²·INV_R2)² where d is
   * the query-to-nearest-node distance. Since the nearest lattice node is always
   * well within KERNEL_R (spacing ≈ D_AVG, radius = 1.5·D_AVG), that center
   * weight stays near ~0.8 for every on-sphere query, so the real total weight is
   * O(1) and never legitimately approaches this floor. The guard therefore only
   * trips on a degenerate/empty walk; it sits far below the smallest legitimate
   * total weight, not between two live values, so it cannot cull a real pixel.
   */
  static constexpr float KERNEL_MIN_TOTAL_WEIGHT = 1e-4f;

  /**
   * @brief Squared Euclidean distance between two points (no sqrt).
   * @param a First point.
   * @param b Second point.
   * @return Squared distance between `a` and `b`.
   */
  static float dist2(const Vector &a, const Vector &b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
  }

  /**
   * @brief Refines a cubemap-LUT seed to its closest lattice node.
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Cached lattice node positions, indexed by node id.
   * @param center_node Seed node id from the cubemap LUT.
   * @return Node id of the genuine nearest node among the seed and its neighbors.
   * @details The CubemapLUT quantizes the query direction to a face cell, so its
   * seed can be a not-quite-nearest node; centering the kernel on the genuine
   * nearest node removes that quantization bias. Shared by both systems; without
   * it the interpolation inherits the cubemap grid artifact.
   */
  static int refine_nearest_node(const Vector &rv, const Vector *nodes,
                                 int center_node) {
    float best_d = dist2(rv, nodes[center_node]);
    int best_node = center_node;
    for_each_neighbor(center_node, [&](int ni) {
      float d = dist2(rv, nodes[ni]);
      if (d < best_d) {
        best_d = d;
        best_node = ni;
      }
    });
    return best_node;
  }

  /**
   * @brief Vertex-shader seed: tags a fragment with its cubemap-LUT node id.
   * @param frag Fragment whose pos seeds frag.v0 with the nearest-node id.
   * @details Shared by both systems' render() vertex shaders; the fragment
   * shader refines this face-quantized seed to the true nearest node per
   * sub-sample (see refine_nearest_node).
   */
  void seed_face_lut(Fragment &frag) {
    Vector rv = orientation.unorient(frag.pos);
    frag.v0 = static_cast<float>(cube_lut.lookup(rv));
  }

  /**
   * @brief Lands the latest ping-pong generation back into persistent state.
   * @tparam T Fixed-point species sample type.
   * @tparam N Species count.
   * @param persistent Destination persistent buffers, one per species.
   * @param latest Latest-generation buffers after the substep loop.
   * @param count Samples per buffer (RD_N).
   * @details An even substep count already leaves the latest generation in the
   * persistent buffers; an odd count leaves it in scratch, so copy each species
   * buffer back before the caller's ScratchScope pops. Agnostic to species
   * count and fixed-point width, so both systems share one land-back.
   */
  template <typename T, size_t N>
  static void land_back(const std::array<T *, N> &persistent,
                        const std::array<T *, N> &latest, size_t count) {
    if (latest[0] == persistent[0])
      return;
    for (size_t i = 0; i < N; ++i)
      std::memcpy(persistent[i], latest[i], count * sizeof(T));
  }

  /**
   * @brief Runs `steps` ping-ponged physics substeps and lands the final
   *        generation back in persistent state.
   * @tparam T Fixed-point species sample type.
   * @tparam N Species count.
   * @tparam StepFn Callable invoked as step(cur, nxt) to advance one substep.
   * @param steps Number of physics substeps to advance.
   * @param persistent Persistent state buffers, one per species.
   * @param scratch Scratch ping-pong buffers, one per species.
   * @param step Physics kernel reading the N `cur` buffers and writing the N
   * `nxt` buffers for one substep.
   * @details Ping-pongs between the persistent state and the scratch buffers. An
   * odd `steps` count leaves the final generation in scratch, so land_back
   * copies it back before the caller's ScratchScope pops. Agnostic to species
   * count and fixed-point width, so both systems share one substep driver.
   */
  template <typename T, size_t N, typename StepFn>
  void advance_substeps(int steps, const std::array<T *, N> &persistent,
                        const std::array<T *, N> &scratch, StepFn &&step) {
    std::array<T *, N> cur = persistent;
    std::array<T *, N> nxt = scratch;
    for (int k = 0; k < steps; ++k) {
      step(cur, nxt);
      for (size_t i = 0; i < N; ++i)
        std::swap(cur[i], nxt[i]);
    }
    land_back(persistent, cur, RD_N);
  }

  /**
   * @brief Installs a Languid random-walk animation of the view orientation.
   * @details Shared by both systems; seeds OpenSimplex2 noise and adds the
   * random-walk animation to the timeline.
   */
  void init_orientation_animation() {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
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
   * @brief Reserves and fills the shared lattice, then arms the view animation.
   * @details Bundles the three steps every derived `init()` must perform
   * together — reserve the RD_N node array, fill it with the static Fibonacci
   * lattice, and install the orientation random-walk — so a derived class cannot
   * perform a subset and silently ship a frozen or empty view. MUST be called
   * after configure_arenas() and after the derived class's own persistent
   * allocations, since the node array shares the persistent arena.
   */
  void init_lattice() {
    // configure_arenas() must have sized the persistent arena to hold the shared
    // node array before this runs; trap an unconfigured or under-sized arena by
    // contract here rather than as a later generic allocation OOM.
    HS_CHECK(persistent_arena.get_capacity() >= RD_N * sizeof(Vector));
    nodes = static_cast<Vector *>(
        persistent_arena.allocate(RD_N * sizeof(Vector), alignof(Vector)));
    build_nodes(nodes);
    init_orientation_animation();
  }

  /**
   * @brief Invokes `fn(ni)` for each valid K-NN neighbor index of `node`.
   * @tparam Fn Callable accepting a neighbor node id.
   * @param node Center node id whose neighbors are visited.
   * @param fn Callable invoked once per valid neighbor index.
   * @details The lambda inlines, so this is identical codegen to a hand-written
   * neighbor loop. Systems needing several fields per node fuse them into one
   * walk here rather than re-reading the neighbor list once per species.
   *
   * Negative entries are the K-NN sentinel for an unfilled slot and are skipped,
   * so the effective Laplacian degree is the count of valid neighbors (≤ RD_K).
   * A deficient node diffuses slower, but is strictly on the stable side of the
   * explicit-Euler bound (which assumes full degree RD_K). The shipped Fibonacci
   * lattice is full-degree (pinned by test_indices_in_range), so this is a
   * robustness path, not an active inhomogeneity today.
   */
  template <typename Fn> static void for_each_neighbor(int node, Fn &&fn) {
    for (int k = 0; k < RD_K; ++k) {
      int ni = ReactionGraph::neighbors[node][k];
      if (ni >= 0)
        fn(ni);
    }
  }

  /**
   * @brief Wendland C2 kernel walk over `center` and its neighbors.
   * @tparam OnWeight Callable accepting (node_index, weight).
   * @param rv Query direction (unit vector on the sphere).
   * @param nodes Cached lattice node positions, indexed by node id.
   * @param center Center node id of the kernel.
   * @param on_weight Callable invoked as `on_weight(node_index, weight)` for
   * every node inside the support radius.
   * @details The caller accumulates whatever fields it needs and applies the
   * total-weight guard, so this stays agnostic to species count and fixed-point.
   */
  template <typename OnWeight>
  static void kernel_accumulate(const Vector &rv, const Vector *nodes,
                                int center, OnWeight &&on_weight) {
    auto visit = [&](int i) {
      float u = 1.0f - dist2(rv, nodes[i]) * INV_R2;
      if (u > 0)
        on_weight(i, u * u);
    };
    visit(center);
    for_each_neighbor(center, visit);
  }

  Orientation<> orientation;        /**< Current view orientation on the sphere. */
  FastNoiseLite noise;              /**< Noise source driving the orientation walk. */
  ReactionGraph::CubemapLUT cube_lut; /**< Cubemap LUT for fast nearest-node seeding. */
  Timeline timeline;                /**< Animation timeline advancing the orientation. */
  Vector *nodes = nullptr; /**< Fixed Fibonacci-lattice node positions (RD_N), built once by init_lattice() and shared by both systems. */
};
