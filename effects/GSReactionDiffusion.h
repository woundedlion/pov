/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GSReactionDiffusion.h
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 */

#include <array>
#include <cmath>
#include <utility>
#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "effects/ReactionDiffusionBase.h"

// Unit-test accessor (tests/test_effects.h) reaching the private Q16
// conversions and one Gray-Scott substep that the smoke/determinism harness
// cannot pin.
namespace hs_test {
namespace effects_tests {
struct GSWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Gray-Scott reaction-diffusion on a Fibonacci lattice sphere.
 *
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 *
 * @details
 * Two species (A, B) evolve via Gray-Scott dynamics (A·B² autocatalysis with
 * feed/kill) on the shared 7680-node lattice, producing spots/stripes/mazes.
 * Persistent state is Q16 (uint16_t) for the cubic reaction-term precision;
 * substeps integrate in float and quantize back once per frame. Shared
 * lattice/orientation/kernel scaffolding lives in ReactionDiffusionBase.
 *
 * A reaction runs until its field has all but stopped moving, then dissolves off
 * the sphere and reseeds at new sites, so each cycle grows a different form from
 * the same constants. Editing the constants dissolves the current field too.
 *
 * Memory budget (persistent arena, configured 174 KB):
 *   - Cubemap LUT:                  6 × 64² × 2B = 49,152 B
 *   - State:   2 arrays × 7680 × 2B (Q16)        = 30,720 B
 *   - Node XYZ: 7680 × 12B                       = 92,160 B  (fixed lattice, built once)
 *   - Palette LUT: 256 × 8B + 4B                 =  2,052 B  (the extra tenant vs BZ)
 *   - Total:                                       174,084 B (170 KB)
 *
 * Scratch arena (per frame, disjoint phases):
 *   - Physics: float ping-pong 4 × 7680 × 4B                             = 122,880 B
 *   - Raster:  oriented lattice 7680 × 12B + cull flags 2 × 7680 × 1B    = 107,520 B
 */
template <int W, int H>
class GSReactionDiffusion
    : public ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H> {
  using Base = ReactionDiffusionBase<GSReactionDiffusion<W, H>, W, H>;
  friend Base; // draw_frame() forwards to render()

  // Bring dependent-base names into scope (template base requires this).
  using Base::cube_lut;
  using Base::dist2;
  using Base::for_each_neighbor;
  using Base::from_q16;
  using Base::init_lattice;
  using Base::orient_lattice;
  using Base::Q16_INV;
  using Base::RD_K;
  using Base::RD_N;
  using Base::refine_and_accumulate;
  using Base::refine_center;
  using Base::refine_render_center;
  using Base::register_param;
  using Base::seed_blobs;
  using Base::seed_face_lut;
  using Base::to_q16;
  using Base::with_wendland_weight;

public:
  /**
   * @brief Default-constructs the effect; all setup is deferred to init().
   */
  GSReactionDiffusion() = default;

  /**
   * @brief One-time setup: arenas, GUI params, A/B state, cubemap LUT, lattice.
   * @details Carves the persistent arena, registers the GUI params, allocates
   * and seeds the A/B state, and builds the cubemap LUT and lattice nodes once.
   */
  void init() override {
    constexpr size_t PALETTE_BYTES =
        BakedPalette::required_arena_bytes(); // palette.bake
    constexpr size_t PERSISTENT_BYTES = 174 * 1024;
    // render()'s scratch peaks at the larger of the physics phase (4 float
    // ping-pong buffers) and the raster phase (the oriented lattice + cull
    // flags); the two run under disjoint scopes.
    constexpr size_t PHYSICS_SCRATCH_BYTES = 4u * RD_N * sizeof(float);
    constexpr size_t RASTER_SCRATCH_BYTES =
        RD_N * sizeof(Vector) + 2u * RD_N * sizeof(uint8_t);
    constexpr size_t SCRATCH_BYTES =
        PHYSICS_SCRATCH_BYTES > RASTER_SCRATCH_BYTES ? PHYSICS_SCRATCH_BYTES
                                                     : RASTER_SCRATCH_BYTES;
    Base::template configure_rd_arenas<uint16_t, 2, PERSISTENT_BYTES,
                                       PALETTE_BYTES, SCRATCH_BYTES>();

    register_param("Feed", &params.feed, 0.0f, 0.1f);
    register_param("Kill", &params.k, 0.0f, 0.1f);
    // The six-step temporal block scales the Euler timestep by 5/3, which at
    // the joint Speed/diffusion maximum exceeds the linear diffusion bound;
    // step_physics' per-substep [0,1] clamp bounds that corner.
    register_param("dA", &params.d_a, 0.0f, 0.05f);
    register_param("dB", &params.d_b, 0.0f, 0.05f);
    register_param("Speed", &params.dt, 0.1f, 3.0f);

    state.A = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    state.B = static_cast<uint16_t *>(
        persistent_arena.allocate(RD_N * sizeof(uint16_t), alignof(uint16_t)));
    for (int i = 0; i < RD_N; i++) {
      state.A[i] = 65535;
      state.B[i] = 0;
    }

    palette.bake(persistent_arena, make_palette());

    cube_lut.build(persistent_arena);
    init_lattice();
    seed_blobs(state.B, NUM_SEED_CLUSTERS);
    reaction_edited(); // latch the defaults; frame 1 is not an edit
  }

private:
  // Test seam: lets unit tests reach the Q16 helpers, step_physics, and params
  // without exposing them to production callers.
  friend struct ::hs_test::effects_tests::GSWhiteBox;

  /**
   * @brief Fully-saturated B blobs seeded per reaction.
   * @details Nucleation sites for the GS instability to grow from; without them
   * the uniform A=1/B=0 field never moves.
   */
  static constexpr int NUM_SEED_CLUSTERS = 30;
  /** @brief Physics work performed per frame at the original 8 fps cadence. */
  static constexpr int BASELINE_STEPS_PER_FRAME = 16;
  /** @brief Rendered frames the dissolve takes to convert every node back to
   * rest; 6.4 s at the 16 fps cadence. */
  static constexpr int DISSOLVE_FRAMES = 103;
  /**
   * @brief Rendered frames a fresh reaction runs before the stabilization
   * detector arms; 9.6 s at the 16 fps cadence.
   * @details A young field has only NUM_SEED_CLUSTERS active sites, so its mean
   * |dB| sits under MEAN_DB_STABLE and would read as stalled at birth.
   */
  static constexpr int MIN_GROW_FRAMES = 154;
  /** @brief Consecutive sub-floor frames that count as stabilized. */
  static constexpr int STABLE_HOLD_FRAMES = 39;
  /** @brief Fraction of the dissolve spent fading a node before conversion. */
  static constexpr float DISSOLVE_FADE_FRACTION = 1.0f / 16.0f;
  /**
   * @brief Mean per-node |dB| per frame below which the field counts as
   * settled, at DEFAULT_DT and BASELINE_STEPS_PER_FRAME; the detector rescales
   * it by params.dt / DEFAULT_DT and by EVOLUTION_STEPS_PER_FRAME /
   * BASELINE_STEPS_PER_FRAME, so neither the 30x Speed range nor the frame's
   * physics budget moves the stabilization point.
   * @details Loose relative to the 1.1e-6..4.0e-6 Q16 chatter a converged field
   * floors at; fires at ~222 baseline frames.
   */
  static constexpr float MEAN_DB_STABLE = 2.0e-4f;
  /** @brief Speed the stabilization floor is calibrated at. */
  static constexpr float DEFAULT_DT = 2.5f;
  /** @brief Original-size substep equivalents advanced per rendered frame. */
  static constexpr int EVOLUTION_STEPS_PER_FRAME = 10;
  /**
   * @brief Euler integrations performed per rendered frame.
   * @details Six 5/3-sized integrations cover the same simulated interval as
   * ten original-size integrations.
   */
  static constexpr int STEPS_PER_FRAME = 6;
  static constexpr float STEP_DT_SCALE =
      static_cast<float>(EVOLUTION_STEPS_PER_FRAME) / STEPS_PER_FRAME;
  static_assert(STEPS_PER_FRAME % 2 == 0,
                "GS ping-pong state must land in its input buffers");

  /**
   * @brief Advances the frame-scratch species state through all substeps.
   * @param state Input buffers receiving the final even-numbered generation.
   * @param scratch Alternate frame-scratch ping-pong buffers.
   * @param step Physics kernel invoked for each substep.
   */
  template <typename StepFn>
  static void advance_substeps(const std::array<float *, 2> &state,
                               const std::array<float *, 2> &scratch,
                               StepFn &&step) {
    std::array<float *, 2> cur = state;
    std::array<float *, 2> nxt = scratch;
    for (int k = 0; k < STEPS_PER_FRAME; ++k) {
      step(cur, nxt);
      std::swap(cur[0], nxt[0]);
      std::swap(cur[1], nxt[1]);
    }
  }
  /**
   * @brief Lower bound of the B render band: below this, pixels are transparent;
   * [B_COLOR_FLOOR, B_COLOR_FLOOR + 1/B_COLOR_SCALE] maps to the full palette
   * range [0,1].
   */
  static constexpr float B_COLOR_FLOOR = 0.1f;
  static constexpr float B_COLOR_SCALE = 4.0f; /**< Slope mapping B above the
                                                    floor into palette t. */
  /**
   * @brief Cull threshold; coincides with the color floor so there is no band
   * between the two: a pixel is either fully transparent or on the gradient.
   * @details A cull below the floor would map b in [cull, floor) to t==0,
   * rendering an opaque flat plateau of the lowest palette color.
   */
  static constexpr float B_CULL_THRESHOLD = B_COLOR_FLOOR;

  HS_COLD_MEMBER static GenerativePalette make_palette() {
    return GenerativePalette{EffectPaletteRecipes::gs_reaction_diffusion(
        EffectPaletteRecipes::random_base_turns())};
  }

  /**
   * @brief Fades nodes approaching the dissolve frontier, then holds them at
   * rest.
   * @param phase Dissolve progress in [0, 1]; the cleared fraction.
   * @details The band ahead of the frontier eases A/B toward rest over
   * DISSOLVE_FADE_FRACTION of the window. The swept set is re-cleared every
   * frame because B is autocatalytic and otherwise refills from its neighbors.
   */
  void convert_below(float phase) {
    for (int i = 0; i < RD_N; i++) {
      float h = hash01(static_cast<uint32_t>(i), transition.dissolve_seed);
      if (h < phase) {
        state.A[i] = 65535;
        state.B[i] = 0;
      } else if (h < phase + DISSOLVE_FADE_FRACTION) {
        float keep = (h - phase) * (1.0f / DISSOLVE_FADE_FRACTION);
        state.A[i] = static_cast<uint16_t>(
            65535.0f - (65535.0f - state.A[i]) * keep + 0.5f);
        state.B[i] =
            static_cast<uint16_t>(static_cast<float>(state.B[i]) * keep + 0.5f);
      }
    }
  }

  /**
   * @brief Ends a dissolve by seeding the next reaction at fresh cluster sites.
   * @details The field is already at rest (every node converted), so seeding
   * alone reproduces init()'s starting condition at new random sites;
   * feed/k are the user's and are left alone.
   */
  HS_COLD_MEMBER void start_reaction() {
    transition.dissolve_frames = -1;
    transition.grow_frames = 0;
    transition.stable_frames = 0;
    {
      HS_PROFILE(grd_palette_rebake);
      palette.rebake(make_palette());
    }
    seed_blobs(state.B, NUM_SEED_CLUSTERS);
  }

  /**
   * @brief Starts a dissolve at a fresh node ordering.
   */
  void begin_dissolve() {
    transition.dissolve_frames = 0;
    transition.dissolve_seed = static_cast<uint32_t>(hs::random()());
  }

  /**
   * @brief Reports whether the reaction constants moved since the last frame.
   * @return True on the first frame after any of feed/k/dA/dB changes.
   * @details Latches the current values, so a slider drag reports once per
   * distinct value. Speed is excluded: it sets the integration rate, not the
   * reaction.
   */
  bool reaction_edited() {
    bool changed =
        params.feed != transition.last_feed || params.k != transition.last_k ||
        params.d_a != transition.last_d_a || params.d_b != transition.last_d_b;
    transition.last_feed = params.feed;
    transition.last_k = params.k;
    transition.last_d_a = params.d_a;
    transition.last_d_b = params.d_b;
    return changed;
  }

  /**
   * @brief Runs the reaction lifecycle: edit and stabilization detection, then
   *        the dissolve.
   * @param mean_db Mean per-node |dB| between the frame's start and end states.
   *        Substep motion that cancels over the frame reads as zero.
   * @details Dissolves when the user edits the reaction, or once the field has
   * stalled for STABLE_HOLD_FRAMES. Edits mid-dissolve are absorbed by the
   * in-flight one, which reseeds into whatever the constants read at its end.
   */
  void advance_transition(float mean_db) {
    if (transition.dissolve_frames >= 0) {
      transition.dissolve_frames++;
      convert_below(static_cast<float>(transition.dissolve_frames) /
                    DISSOLVE_FRAMES);
      // Latch edits made mid-dissolve; this dissolve already covers them.
      reaction_edited();
      if (transition.dissolve_frames >= DISSOLVE_FRAMES)
        start_reaction();
      return;
    }
    if (reaction_edited()) {
      begin_dissolve();
      return;
    }
    transition.grow_frames++;
    if (transition.grow_frames < MIN_GROW_FRAMES)
      return;
    // Per-frame |dB| scales with the timestep, so the floor tracks Speed.
    const float floor_db = MEAN_DB_STABLE * (params.dt * (1.0f / DEFAULT_DT)) *
                           (static_cast<float>(EVOLUTION_STEPS_PER_FRAME) /
                            BASELINE_STEPS_PER_FRAME);
    transition.stable_frames =
        mean_db < floor_db ? transition.stable_frames + 1 : 0;
    if (transition.stable_frames >= STABLE_HOLD_FRAMES)
      begin_dissolve();
  }

  /**
   * @brief Advances one Gray-Scott substep into the next buffers (Jacobi).
   * @param c_a Current A field (read-only), float in [0, 1] per node.
   * @param c_b Current B field (read-only), float in [0, 1] per node.
   * @param n_a Next A field (write target), float in [0, 1] per node.
   * @param n_b Next B field (write target), float in [0, 1] per node.
   * @details Gray-Scott: dA/dt = dA·∇²A - A·B² + feed·(1-A);
   * dB/dt = dB·∇²B + A·B² - (k+feed)·B. Double-buffered Jacobi: reads current
   * buffers, writes next; the caller owns the ping-pong (see render()). The
   * [0, 1] clamp saturates explicit-Euler overshoot past the stability bound
   * (see "Speed"); substeps stay in float so the Q16 state quantizes once per
   * frame, not once per substep.
   */
  HS_O3_FN void step_physics(const float *__restrict c_a,
                             const float *__restrict c_b, float *__restrict n_a,
                             float *__restrict n_b) {
    const float feed = params.feed;
    const float k = params.k;
    const float d_a = params.d_a;
    const float d_b = params.d_b;
    const float dt = params.dt * STEP_DT_SCALE;
    for (int i = 0; i < RD_N; i++) {
      float a = c_a[i];
      float b = c_b[i];

      float sum_a = 0, sum_b = 0;
      for_each_neighbor(i, [&](int ni) {
        sum_a += c_a[ni];
        sum_b += c_b[ni];
      });
      float l_a = sum_a - RD_K * a;
      float l_b = sum_b - RD_K * b;

      float abb = a * b * b;
      n_a[i] =
          hs::clamp(a + (d_a * l_a - abb + feed * (1.0f - a)) * dt, 0.0f, 1.0f);
      n_b[i] =
          hs::clamp(b + (d_b * l_b + abb - (k + feed) * b) * dt, 0.0f, 1.0f);
    }
  }

  /**
   * @brief Kernel-weighted sample of the B concentration at a point.
   * @param p Query point on the sphere.
   * @param seed Seed node id from the cubemap LUT, refined to the true nearest
   * inside the fused stencil walk.
   * @param nodes Node positions in the same frame as `p`.
   * @return Support-radius weighted average of B in [0, 1]; 0 if no node is
   * within the support radius.
   * @details Off the render path: shade_pixel gathers the stencil once per
   * pixel and re-weights it inline. This one-sample form is the oracle
   * tests/test_effects.h bounds that shared stencil against.
   */
  float interpolate_b(const Vector &p, int seed, const Vector *nodes) const {
    float tw = 0, wb = 0;
    refine_and_accumulate(p, nodes, seed, [&](int i, float w) {
      wb += from_q16(state.B[i]) * w;
      tw += w;
    });
    // Guard the division: returning 0 (not a 0/0 NaN) keeps the value cullable
    // by render()'s `b < B_CULL_THRESHOLD` test downstream — a NaN compares false
    // there, slips through, and poisons palette.get().
    if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
      return 0.0f;
    return wb / tw;
  }

  /**
   * @brief Shades one pixel's four sub-samples through an inlinable typed path.
   * @tparam Grid Scan::Shader::SsaaGrid type supplying the sub-pixel offsets.
   * @param seed Cubemap-LUT seed node id, or -1 for a culled pixel.
   * @param center_rv World-space direction at the pixel center.
   * @param world_nodes Oriented lattice node positions.
   * @param grid Row's SSAA sub-pixel grid.
   * @param x Pixel column.
   * @return The finished alpha-premultiplied pixel.
   * @details Accepts seeds inside a proven nearest-node radius immediately;
   * boundary pixels check all six neighbors. The center stencil is shared
   * across the four sub-pixel samples.
   */
  template <typename Grid>
  HS_O3_FN Pixel shade_pixel(int seed, const Vector &center_rv,
                             const Vector *world_nodes, const Grid &grid,
                             int x) const {
    if (seed < 0)
      return Pixel(0, 0, 0);

    int center = refine_render_center(center_rv, world_nodes, seed);
    Vector spos[RD_K + 1];
    uint16_t sb[RD_K + 1];
    spos[0] = world_nodes[center];
    sb[0] = state.B[center];
    int k = 1;
    for_each_neighbor(center, [&](int ni) {
      spos[k] = world_nodes[ni];
      sb[k] = state.B[ni];
      ++k;
    });

    uint32_t accum_r = 0, accum_g = 0, accum_b = 0;
    for (int i = 0; i < 4; ++i) {
      Vector v = grid.at(x, i);
      float tw = 0.0f, wb = 0.0f;
      for (int j = 0; j < RD_K + 1; ++j)
        with_wendland_weight(dist2(v, spos[j]),
                             [&](float w) __attribute__((always_inline)) {
                               wb += sb[j] * w;
                               tw += w;
                             });
      if (tw <= Base::KERNEL_MIN_TOTAL_WEIGHT)
        continue;
      float b = wb * (Q16_INV / tw);
      if (b < B_CULL_THRESHOLD)
        continue;

      float t = hs::clamp((b - B_COLOR_FLOOR) * B_COLOR_SCALE, 0.0f, 1.0f);
      Pixel sample = palette.get_color_unit(t);
      accum_r += (static_cast<uint32_t>(sample.r) + 2u) >> 2;
      accum_g += (static_cast<uint32_t>(sample.g) + 2u) >> 2;
      accum_b += (static_cast<uint32_t>(sample.b) + 2u) >> 2;
    }
    return Pixel(static_cast<uint16_t>(accum_r > 65535u ? 65535u : accum_r),
                 static_cast<uint16_t>(accum_g > 65535u ? 65535u : accum_g),
                 static_cast<uint16_t>(accum_b > 65535u ? 65535u : accum_b));
  }

  /**
   * @brief Builds per-node two-ring "renderable" flags for the B field.
   * @param b Per-node B concentrations, Q16.
   * @param hot1 Scratch: per-node flag, set when any of {node, neighbors}
   *        reaches the threshold.
   * @param hot2 Output: per-node flag, set when any node within two hops
   *        reaches the threshold.
   * @param count Node count.
   * @param threshold Q16 render floor.
   * @details A kernel sample is a convex average over the refined stencil and
   * the refined center is at most one hop from the seed, so a seed whose
   * two-ring sits entirely below the floor cannot produce a renderable sample —
   * culling on !hot2[seed] is exact, not approximate.
   */
  HS_O3_FN static void fill_hot_flags(const uint16_t *b, uint8_t *hot1,
                                      uint8_t *hot2, int count,
                                      uint16_t threshold) {
    for (int i = 0; i < count; ++i) {
      bool hot = b[i] >= threshold;
      for (int k = 0; k < RD_K && !hot; ++k)
        hot = b[ReactionGraph::neighbors[i][k]] >= threshold;
      hot1[i] = hot;
    }
    for (int i = 0; i < count; ++i) {
      bool hot = hot1[i];
      for (int k = 0; k < RD_K && !hot; ++k)
        hot = hot1[ReactionGraph::neighbors[i][k]];
      hot2[i] = hot;
    }
  }

  /**
   * @brief Advances the sim STEPS_PER_FRAME substeps and rasterizes the B field.
   * @param canvas Destination canvas to draw the sphere into.
   * @details Rasterizes the B field onto the sphere via the orientation-aware
   * SSAA shader pipeline after advancing the simulation.
   */
  void render(Canvas &canvas) {
    HS_PROFILE(grd_render);
    ScratchScope frame_guard(scratch_arena_a);
    float mean_db = 0.0f;
    {
      // The substeps ping-pong in float; the Q16 state converts in once and
      // quantizes back once per frame, not once per substep.
      HS_PROFILE(grd_simulate);
      ScratchScope physics_guard(scratch_arena_a);
      float *cur_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *cur_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *nxt_a = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));
      float *nxt_b = static_cast<float *>(
          scratch_arena_a.allocate(RD_N * sizeof(float), alignof(float)));

      for (int i = 0; i < RD_N; i++) {
        cur_a[i] = from_q16(state.A[i]);
        cur_b[i] = from_q16(state.B[i]);
      }
      advance_substeps(std::array<float *, 2>{cur_a, cur_b},
                       std::array<float *, 2>{nxt_a, nxt_b},
                       [&](auto &cur, auto &nxt) {
                         step_physics(cur[0], cur[1], nxt[0], nxt[1]);
                       });
      uint32_t db_sum_q16 = 0;
      for (int i = 0; i < RD_N; i++) {
        state.A[i] = to_q16(cur_a[i]);
        uint16_t next_b = to_q16(cur_b[i]);
        int db = static_cast<int>(next_b) - state.B[i];
        db_sum_q16 += static_cast<uint32_t>(db < 0 ? -db : db);
        state.B[i] = next_b;
      }
      mean_db = static_cast<float>(db_sum_q16) * (Q16_INV / RD_N);
    }
    advance_transition(mean_db);

    // Physics scratch is popped; the raster phase reuses the arena for the
    // oriented lattice so the kernel walks stay in world space, plus the
    // two-ring cull flags.
    HS_PROFILE(grd_rasterize);
    Vector *world_nodes = orient_lattice();
    uint8_t *hot1 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    uint8_t *hot2 = static_cast<uint8_t *>(scratch_arena_a.allocate(RD_N, 1));
    fill_hot_flags(state.B, hot1, hot2, RD_N, to_q16(B_CULL_THRESHOLD));

    // Seed the cubemap lookup once per pixel center; a seed whose two-ring
    // sits below the render floor is culled for the whole pixel (v0 = -1).
    auto vertex_shader = [&](Fragment &frag) {
      seed_face_lut(frag);
      if (!hot2[static_cast<int>(frag.v0)])
        frag.v0 = -1.0f;
    };

    auto pixel_shader = [&](Fragment &frag, const auto &grid, int x) -> Pixel {
      return shade_pixel(static_cast<int>(frag.v0), frag.pos, world_nodes, grid,
                         x);
    };

    {
      HS_PROFILE(grd_shader_draw);
      Scan::Shader::draw_grid<W, H>(canvas, vertex_shader, pixel_shader);
    }
  }

  /**
   * @brief Persistent Q16 state buffers for the two species.
   */
  struct {
    uint16_t *A = nullptr,
             *B = nullptr; /**< Per-node A/B concentrations, Q16. */
  } state;

  /**
   * @brief Auto-transition state: current reaction's lifetime and dissolve
   *        progress.
   */
  struct {
    int grow_frames = 0;        /**< Frames since this reaction was seeded. */
    int stable_frames = 0;      /**< Consecutive sub-floor frames. */
    int dissolve_frames = -1;   /**< Dissolve progress; -1 when inactive. */
    uint32_t dissolve_seed = 0; /**< Per-transition node-order hash seed. */
    /** Reaction constants as of the last frame; reaction_edited() latches them
     *  to spot a user edit, seeded from params in init() so frame 1 is
     *  clean. */
    float last_feed = 0.0f, last_k = 0.0f, last_d_a = 0.0f, last_d_b = 0.0f;
  } transition;

  /** @brief 16-bit LUT baked from the generative palette mapping B to RGB. */
  BakedPalette palette;

  /**
   * @brief GUI-tunable Gray-Scott parameters.
   */
  struct Params {
    float feed = 0.04f; /**< Feed rate of A. */
    float k = 0.06f;    /**< Kill rate of B. */
    float d_a = 0.02f;  /**< Diffusion coefficient of A. */
    float d_b = 0.01f;  /**< Diffusion coefficient of B. */
    float dt = 2.5f;    /**< Integration timestep (Speed slider). */
  } params;
  static_assert(Params{}.dt == DEFAULT_DT,
                "the stabilization floor is calibrated at the Speed default");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(GSReactionDiffusion)
