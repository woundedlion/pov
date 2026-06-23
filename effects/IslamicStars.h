/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"
#include <algorithm>

/**
 * @brief Effect that displays a sequence of Islamic-geometry polyhedra,
 *        cross-fading one shape into the next while ripples distort the mesh.
 * @tparam W Target canvas width in pixels.
 * @tparam H Target canvas height in pixels.
 */
template <int W, int H> class IslamicStars : public Effect {

public:
  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  FLASHMEM IslamicStars() : Effect(W, H), filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk, the recurring ripple bursts, and the
   *        first shape.
   */
  void init() override {
    // scratch arenas + room for BakedPaletteBank (~15KB) in persistent
    configure_arenas(GLOBAL_ARENA_SIZE - (120 + 120) * 1024, 120 * 1024,
                     120 * 1024);

    palette_bank_.bake_all(persistent_arena);

    // Set ripple defaults BEFORE registering their pointers: registerParam
    // snaps *ptr as the slider default, so these must already hold the intended
    // runtime values.
    ripple_gen.template_params.amplitude = 0.4f;
    ripple_gen.template_params.thickness = 0.7f;
    ripple_gen.template_params.decay = 0.1f;

    registerParam("Duration", &params.duration, 48.0f, 192.0f);
    registerParam("Fade", &params.fade, 0.0f, 96.0f);
    // Burst and Ripp Dur ranges are clamped to the 8-slot ripple pool: with
    // Burst <= 4 and Ripp Dur <= 144, the current burst (4) plus the still-live
    // previous burst (4) peaks at exactly 8, so no spawns are dropped.
    registerParam("Burst", &params.burst_size, 1.0f, 4.0f);
    registerParam("Ripp Amp", &ripple_gen.template_params.amplitude, 0.0f, 1.0f);
    registerParam("Ripp Width", &ripple_gen.template_params.thickness, 0.1f, 1.0f);
    registerParam("Ripp Decay", &ripple_gen.template_params.decay, 0.0f, 5.0f);
    registerParam("Ripp Dur", &ripple_duration, 30.0f, 144.0f);
    registerParam("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    // One immediate ripple burst, then a recurring one every 96 frames.
    timeline.add(0, Animation::PeriodicTimer(
                        0, [this](Canvas &canvas) { ripple(canvas); }, false));
    timeline.add(0, Animation::PeriodicTimer(
                        96, [this](Canvas &canvas) { ripple(canvas); }, true));
    spawn_shape();
  }

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return true; the strip blanks to black immediately after each column is
   *         shown, so every column reads as a sharp slice with dark gaps
   *         between columns rather than persisting across the sweep.
   */
  bool strobe_columns() const override { return true; }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    ripple_gen.prepare_frame();
    timeline.step(canvas);
  }

private:
  Orientation<> orientation;
  Timeline timeline;
  Pipeline<W, H> filters;
  RippleTransformer<8> ripple_gen;
  FastNoiseLite noise;
  float ripple_duration = 80.0f;
  int solid_idx = -1;
  MeshCarousel carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;
  // Per-slot palette indices; value-init so a missed shuffle reads palette 0,
  // not garbage. The shuffle-before-draw rule is a convention (each spawn_shape
  // shuffles the back slot), enforced the same way as HankinSolids' identical
  // palettes_slots — not by a runtime guard.
  std::array<int, NUM_PALETTES> palettes_slots[2] = {};

  /**
   * @brief Spawns one burst of burst_size ripples from a random origin,
   *        staggered 16 frames apart, each expanding over ripple_duration
   *        frames.
   * @param canvas Unused render target for the timer callback signature.
   */
  void ripple(Canvas &) {
    Vector origin = random_vector();
    for (int i = 0; i < (int)params.burst_size; i++) {
      ripple_gen.spawn(i * 16, origin, PI_F / ripple_duration,
                       static_cast<int>(ripple_duration));
    }
  }

  /**
   * @brief Orients and ripple-distorts base_state, then rasterizes it with a
   *        per-face palette lookup.
   * @param canvas Render target receiving the rasterized mesh.
   * @param opacity Sprite's current fade alpha in [0, 1].
   * @param base_state Undistorted source mesh to transform and draw.
   * @param faceIndices Maps each face to its topology class.
   * @param palette_idx Assigns a palette per topology class.
   */
  void draw_shape(Canvas &canvas, float opacity, const MeshState &base_state,
                  const ArenaVector<int> &faceIndices,
                  const std::array<int, NUM_PALETTES> &palette_idx) {
    if (opacity <= 0.005f)
      return;
    ScratchScope _a(scratch_arena_a);
    MeshState transformed_state;
    OrientTransformer camera(orientation);
    MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                       ripple_gen, camera);

    const int *raw_indices = faceIndices.data();

    const int num_faces = static_cast<int>(faceIndices.size());
    // Color each fragment by its face's topology class, shaded by edge distance
    // (gain 1.0). Shared with HankinSolids via shade_mesh_topology.
    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      frag.color = shade_mesh_topology(frag, raw_indices, num_faces,
                                       palette_bank_, palette_idx, 1.0f, opacity);
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  /**
   * @brief Advances to the next solid, generates it into the carousel's back
   *        slot with a freshly shuffled palette, makes it the front, and
   *        schedules a cross-fading sprite plus the next spawn_shape call.
   */
  void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    // The back slot is the render/generate target for this cycle. Captured once
    // so the palette shuffle, the draw_fn closure, and the generate target all
    // reference the same slot.
    int back = 1 - carousel.front_index();
    MeshPaletteBank::shuffle_indices(palettes_slots[back]);

    int idx = solid_idx; // capture for generate lambda

    auto draw_fn = [this, back](Canvas &canvas, float opacity) {
      const MeshState &mesh = carousel.slot(back);
      this->draw_shape(canvas, opacity, mesh, mesh.topology,
                       palettes_slots[back]);
    };

    // 1. Free the back slot and compact, rebaking palettes into the fresh
    // arena instead of tracking them through the evacuation (saves
    // fragmentation/OOM).
    carousel.compact_keep_front(
        [this](Arena &arena) { palette_bank_.bake_all(arena); });

    // 2. Generate new shape
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh mesh = solids[idx].generate(a, b);
      carousel.slot(back).clear();
      MeshOps::compile(mesh, carousel.slot(back), target);
    });

    // classify_faces_by_topology uses both scratch arenas as transient
    // workspace (its result lands in persistent_arena), so bracket it with
    // ScratchScope to free only its own allocations on exit — the same idiom
    // HankinSolids::classify_mesh_topology uses, rather than a bare reset() that
    // would discard any allocation a caller already holds in these arenas.
    {
      ScratchScope _a(scratch_arena_a);
      ScratchScope _b(scratch_arena_b);
      MeshOps::classify_faces_by_topology(carousel.slot(back), scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // 3. Flip front eagerly for the overlapping sprite
    carousel.set_front(back);

    // Live-read the Duration slider per shape cycle. Clamp fade to dur/2 so
    // the fade-in/out windows never overlap and the respawn delay (dur - fade)
    // stays >= dur/2 — a larger fade piles up sprites and per-frame respawns.
    int dur = static_cast<int>(params.duration);
    int fade = std::min(static_cast<int>(params.fade), dur / 2);
    // Make the comment's invariant checkable: fade never exceeds half the
    // duration, so the respawn delay (dur - fade) stays >= dur/2 and sprites
    // never pile up. Cold setup-time seam (once per shape cycle), so trap on
    // device too with HS_CHECK rather than a debug-stripped assert.
    HS_CHECK(fade <= dur / 2 && dur - fade >= dur / 2);
    timeline.add(0, Animation::Sprite(draw_fn, dur, fade, ease_linear, fade,
                                      ease_linear));

    // Topology indices are left un-reduced: draw_shape applies
    // `topoIdx % NUM_PALETTES` at render time, so reducing the stored topology
    // here would be a redundant double-modulo and an unwanted mutation of the
    // carousel's mesh state.

    const auto &entry = solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)carousel.current().vertices.size(),
            (int)carousel.current().faces.size());

    // Schedule the next shape to begin one fade-out before this one ends, so
    // the two sprites overlap during the cross-fade.
    int next_delay = dur - fade;
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  /**
   * @brief Slider-backed runtime parameters for the effect.
   */
  struct Params {
    float duration = 160.0f; /**< Shape display period, in frames. */
    float fade = 32.0f; /**< Cross-fade window length, in frames. */
    float burst_size = 4.0f; /**< Ripples per burst; float-backed for registerParam. */
    bool debug_bb = false; /**< Whether to draw mesh bounding boxes. */
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
