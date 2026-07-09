/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

#include <algorithm>

/**
 * @brief Renders Hankin interlace patterns over Platonic/Archimedean solids.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Sweeps the interlace angle continuously and periodically morphs
 * between solids. Faces are colored by topology class via shuffled mesh
 * palettes.
 */
template <int W, int H> class HankinSolids : public Effect {
public:
  /**
   * @brief Constructs the effect with a W x H canvas and empty filter pipeline.
   */
  FLASHMEM HankinSolids() : Effect(W, H, {.strobe = true}), filters() {}

  /**
   * @brief Sizes arenas, registers params, loads the first solid, and starts
   * the interlace sweep/morph cycle.
   */
  void init() override {
    // scratch_a (24 KB) covers its largest non-overlapping peak: the render path
    // (draw_mesh transforms the front mesh into scratch_a, then Scan::Mesh::draw
    // stacks an SDF::FaceScratchBuffer on top), which for the heaviest hankin
    // mesh exceeds the generation/classify peak. scratch_b (32 KB) covers the
    // larger of generation intermediates or the morph-cycle Persist set. The
    // device H=144 high-water isn't in CI; confirm the peak on hardware.
    configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                     32 * 1024);
    register_param("Intensity", &params.intensity, 0.0f, 5.0f);
    register_animated_param("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);
    register_param("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    solid_idx = 0;

    palette_bank_.bake_all(persistent_arena);

    load_shape(carousel.current(), compiled_hankin,
               palettes_slots[carousel.front_index()], solid_idx,
               params.hankin_angle);

    start_hankin_cycle();
  }

  /**
   * @brief Advances the timeline by one frame and renders into the canvas.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;

  /**
   * @brief Builds and finalizes the idx-th base solid from the simple-solids
   * collection.
   * @param idx Index into the simple-solids collection.
   * @param a Scratch arena used during generation.
   * @param b Secondary scratch arena used during generation.
   * @return The finalized polygon mesh.
   */
  PolyMesh generate_base_solid(int idx, Arena &a, Arena &b) {
    auto solids = Solids::Collections::get_simple_solids();
    hs::log("Loading shape: '%s'", solids[idx].name);
    return Solids::finalize_solid(solids[idx].generate(a, b), a);
  }

  /**
   * @brief Classifies mesh faces into topology groups.
   * @param mesh Mesh whose faces are classified; the result is stored in
   * persistent_arena.
   * @details Saves/restores the scratch high-water marks rather than hard-
   * resetting to base, so a caller's prior allocations in the shared scratch
   * arenas survive.
   */
  void classify_mesh_topology(MeshState &mesh) {
    // ScratchScope frees only this call's own allocations, preserving prior
    // caller allocations in these shared arenas that a bare reset() would drop.
    ScratchScope a_guard(scratch_arena_a);
    ScratchScope b_guard(scratch_arena_b);
    MeshOps::classify_faces_by_topology(mesh, scratch_arena_a, scratch_arena_b,
                                        persistent_arena);
  }

  /**
   * @brief Runs the full load pipeline: generate, compile, evaluate, classify,
   * and assign palette indices.
   * @param out_mesh Receives the evaluated and classified mesh.
   * @param out_hankin Receives the compiled Hankin pattern.
   * @param out_palette_idx Receives the shuffled palette indices.
   * @param idx Index of the solid to load.
   * @param angle Interlace angle in radians.
   */
  void load_shape(MeshState &out_mesh, CompiledHankin &out_hankin,
                  std::array<int, NUM_PALETTES> &out_palette_idx, int idx,
                  float angle) {
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = generate_base_solid(idx, a, b);

      out_hankin = CompiledHankin();
      MeshOps::compile_hankin(base, out_hankin, target, a);

      out_mesh.clear();
      MeshOps::update_hankin(out_hankin, out_mesh, target, angle);
    });

    classify_mesh_topology(out_mesh);
    MeshPaletteBank::shuffle_indices(out_palette_idx);
  }

  /**
   * @brief Camera-rotates and rasterizes one mesh, coloring each face by its
   * topology class and shading edges by distance.
   * @param canvas Target canvas to draw into.
   * @param mesh Source mesh in model space.
   * @param topology Per-face topology-class indices.
   * @param palette_idx Maps topology class to a palette in the mesh bank.
   * @param opacity Output alpha in [0, 1].
   */
  void draw_mesh(Canvas &canvas, const MeshState &mesh,
                 const ArenaVector<int> &topology,
                 const std::array<int, NUM_PALETTES> &palette_idx,
                 float opacity) {
    if (mesh.vertices.is_empty() || opacity < 0.01f)
      return;

    ScratchScope scratch_a_guard(scratch_arena_a);
    MeshState rotated_mesh;
    OrientTransformer camera(orientation);
    MeshOps::transform(mesh, rotated_mesh, scratch_arena_a, camera);

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      f.color = shade_mesh_topology(f, topology.data(),
                                    static_cast<int>(topology.size()),
                                    palette_bank_, palette_idx, params.intensity,
                                    opacity);
    };

    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  /**
   * @brief Makes the staged (morph-target) compiled mesh active and clears
   * staging.
   */
  void promote_staged_hankin() {
    compiled_hankin = std::move(compiled_hankin_staging);
    compiled_hankin_staging = CompiledHankin();
  }

  /**
   * @brief Schedules one interlace-angle sweep plus the sprite that re-evaluates
   * and draws the front mesh each frame.
   * @details Chains into a morph cycle when the sweep ends. Both the sweep and
   * the render sprite are gated on the same pause flag so grabbing the slider
   * holds the frame instead of blanking it.
   */
  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_linear, false, &anims_paused_)
                        .then([this]() {
                          this->start_morph_cycle();
                        }));

    int front = carousel.front_index();
    // Snapshot the angle-independent counts for the per-frame HS_CHECK below.
    hankin_vertex_count_ = compiled_hankin.static_vertices.size() +
                           compiled_hankin.dynamic_vertices.size();
    hankin_face_count_ = compiled_hankin.face_counts.size();
    timeline.add(
        0, Animation::Sprite(
               [this, front](Canvas &c, float opacity) {
                 // update_hankin re-binds the slot's vectors against
                 // persistent_arena every frame; the angle never changes the
                 // vertex/face counts, so bind reuses the blocks in place.
                 MeshOps::update_hankin(compiled_hankin, carousel.slot(front),
                                        persistent_arena, params.hankin_angle);
                 // Always-on guard: grown counts would leak persistent_arena
                 // every frame on a permanent install.
                 const MeshState &s = carousel.slot(front);
                 HS_CHECK(s.vertices.size() == hankin_vertex_count_ &&
                              s.face_counts.size() == hankin_face_count_,
                          "HankinSolids: per-frame mesh counts changed; the "
                          "persistent re-bind would grow the arena");
                 draw_mesh(c, s, s.topology, palettes_slots[front], opacity);
               },
               DURATION, 0, ease_linear, 0, ease_linear, &anims_paused_));
  }

  /**
   * @brief Builds the next solid into the back slot and schedules a morph to it.
   * @details On completion swaps slots, compacts arenas, and restarts the
   * hankin cycle.
   */
  void start_morph_cycle() {
    constexpr int MORPH_FRAMES = 16;
    auto solids = Solids::Collections::get_simple_solids();
    const int n = static_cast<int>(solids.size());
    // Random next solid, never the current one: offset in [1, n) keeps it distinct.
    int next_idx = (solid_idx + 1 + hs::rand_int(0, n - 1)) % n;

    int old_front = carousel.front_index();
    int new_slot = 1 - old_front;

    load_shape(carousel.slot(new_slot), compiled_hankin_staging,
               palettes_slots[new_slot], next_idx, params.hankin_angle);

    // Set before creating MeshMorph: the draw callbacks reference these.
    morph_old_slot_ = old_front;
    morph_new_slot_ = new_slot;

    timeline.add(
        0, Animation::MeshMorph(
               carousel.slot(old_front), carousel.slot(new_slot),
               persistent_arena, draw_morph_outgoing_fn_,
               draw_morph_incoming_fn_, MORPH_FRAMES, ease_in_out_sin)
               .then([this, next_idx, new_slot]() {
                 solid_idx = next_idx;
                 carousel.set_front(new_slot);
                 promote_staged_hankin();
                 // Manual compaction; only the front slot survives. Clear the
                 // back slot to a fresh MeshState first so its now-dangling
                 // ArenaVectors aren't restored across the reset. Survivors go to
                 // separate scratch arenas so neither overflows.
                 carousel.slot(1 - new_slot) = MeshState();
                 {
                   Persist<CompiledHankin> ph(compiled_hankin, scratch_arena_b,
                                              persistent_arena);
                   Persist<MeshState> pf(carousel.current(), scratch_arena_a,
                                         persistent_arena);
                   Persist<MeshPaletteBank> pp(palette_bank_, scratch_arena_b,
                                               persistent_arena);
                   persistent_arena.reset();
                   hs::log("morph_cycle_then: finished arena compaction");
                 }

                 MeshOps::update_hankin(compiled_hankin, carousel.current(),
                                        persistent_arena, params.hankin_angle);
                 start_hankin_cycle();
               }));
  }

  MeshCarousel<> carousel; /**< Double-slot mesh store for front/back solids;
                                transitions are MeshMorph-driven, so the segue
                                is unused. */
  CompiledHankin compiled_hankin;         /**< Active during the hankin cycle. */
  CompiledHankin compiled_hankin_staging; /**< Built during the morph cycle. */
  std::array<int, NUM_PALETTES> palettes_slots[2] = {}; /**< Per-slot palette indices; value-init so a missed shuffle reads 0, not garbage. */

  int morph_old_slot_ = 0; /**< Outgoing slot index for morph draw callbacks. */
  int morph_new_slot_ = 1; /**< Incoming slot index for morph draw callbacks. */

  size_t hankin_vertex_count_ = 0; /**< Expected per-frame vertex count for the active cycle. */
  size_t hankin_face_count_ = 0;   /**< Expected per-frame face count for the active cycle. */

  /**
   * @brief Draw callback for the outgoing mesh during a morph.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_outgoing_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_old_slot_).topology,
                  palettes_slots[morph_old_slot_], o);
      }};
  /**
   * @brief Draw callback for the incoming mesh during a morph.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_incoming_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_new_slot_).topology,
                  palettes_slots[morph_new_slot_], o);
      }};

  Orientation<> orientation; /**< Current camera orientation. */
  FastNoiseLite noise;       /**< Noise source driving the orientation walk. */
  Timeline timeline;         /**< Schedules sweeps, sprites, and morphs. */
  Pipeline<W, H> filters;    /**< Per-pixel filter pipeline applied on draw. */
  int solid_idx = 0;         /**< Index of the currently displayed solid. */

  /**
   * @brief User-adjustable rendering parameters.
   */
  struct Params {
    float intensity = 1.2f;        /**< Edge-distance shading gain. */
    float hankin_angle = PI_F / 4.0f; /**< Interlace angle in radians. */
    bool debug_bb = false;            /**< Draw face bounding boxes when true. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(HankinSolids)
