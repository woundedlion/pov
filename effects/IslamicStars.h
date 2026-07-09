/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/animation/animation.h"
#include "core/engine/engine.h"

/**
 * @brief Effect that displays a sequence of Islamic-geometry polyhedra,
 *        transitioning one shape into the next while ripples distort the mesh.
 * @tparam W Target canvas width in pixels.
 * @tparam H Target canvas height in pixels.
 */
template <int W, int H>
class IslamicStars : public Effect {

public:
  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  FLASHMEM IslamicStars() : Effect(W, H, {.strobe = true}), filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  void init() override {
    // Asymmetric scratch split sized to the measured recipe peaks (worst
    // generation high-water 112.3 KB in a / 76.5 KB in b; compact_keep_front
    // evacuates the front slot, up to 63.7 KB, through b). The remainder is
    // persistent: carousel slots + BakedPaletteBank (~15 KB). Budgets enforced
    // by test_solids.h's high-water sweeps.
    configure_arenas(GLOBAL_ARENA_SIZE - (114 + 80) * 1024, 114 * 1024,
                     80 * 1024);

    palette_bank_.bake_all(persistent_arena);

    // Set BEFORE registering: register_param snaps *ptr as the slider default.
    // Amplitude is held below the self-fold onset: above amplitude/thickness ~0.2
    // the rippled mesh folds over itself, so faces stop tiling the sphere and
    // stack along the view ray, multiplying rasterizer overdraw (self-occlusion).
    ripple_gen.template_params.amplitude = RIPPLE_AMP_MAX;
    ripple_gen.template_params.thickness = RIPPLE_THICKNESS;
    ripple_gen.template_params.decay = 0.1f;

    register_param("Fade", &params.fade, 0.0f, 96.0f);
    // Per-face fade length range (frames): each face draws a random fade from
    // [lo, hi] as the terminator reaches it, fraying the sweep front.
    register_param("Face Fade Lo", &carousel.segue().fade_frames_min, 0.0f, 32.0f);
    register_param("Face Fade Hi", &carousel.segue().fade_frames_max, 0.0f, 32.0f);
    // Burst/Ripp Dur ranges are clamped to the ripple pool capacity invariant
    // (see the RIPPLE* constants below).
    register_param("Burst", &params.burst_size, 1.0f, (float)BURST_MAX);
    // Amplitude slider capped at the fold-free ceiling; thickness is fixed (not a
    // slider) so amplitude/thickness can never cross the self-fold onset.
    register_param("Ripp Amp", &ripple_gen.template_params.amplitude, 0.0f,
                  RIPPLE_AMP_MAX);
    register_param("Ripp Decay", &ripple_gen.template_params.decay, 0.0f, 5.0f);
    register_param("Ripp Dur", &ripple_duration, 30.0f, (float)RIPPLE_DURATION_MAX);
    register_param("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    spawn_shape();
  }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    ripple_gen.prepare_frame();
    timeline.step(canvas);
  }

private:
  // Ripple-pool sizing: a slot is held from spawn() until the staggered ripple
  // completes. Bursts are scheduled one per shape and each shape's duration is
  // sized so its burst ends a full still window before the next one can start,
  // so only one burst is normally live; the pool holds two so a Ripp Dur/Burst
  // slider change mid-burst still cannot drop a spawn.
  static constexpr int RIPPLE_POOL_SIZE = 8;
  static constexpr int RIPPLE_STAGGER_FRAMES = 16;
  static constexpr int RIPPLE_DURATION_MAX = 143;
  static constexpr int BURST_MAX = 4;
  static constexpr int STILL_FRAMES = 16; /**< 1 s hold (16 fps) between fade and ripple stages. */
  static constexpr float RIPPLE_THICKNESS = 0.7f; /**< Fixed ripple wavelet width (radians). */
  static constexpr float RIPPLE_AMP_MAX = 0.15f;   /**< Fold-free amplitude ceiling at RIPPLE_THICKNESS (amp/thickness < ~0.2 self-fold onset). */
  static_assert(2 * BURST_MAX <= RIPPLE_POOL_SIZE,
                "IslamicStars: ripple pool must hold two overlapping bursts");

  Orientation<> orientation;
  Timeline timeline;
  Pipeline<W, H> filters;
  RippleTransformer<RIPPLE_POOL_SIZE> ripple_gen;
  FastNoiseLite noise;
  float ripple_duration = 80.0f;
  int solid_idx = -1;
  using SegueT = Segue::TerminatorSweep;
  MeshCarousel<SegueT> carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;
  // Per-slot palette indices; value-init so a missed shuffle reads palette 0,
  // not garbage.
  std::array<int, NUM_PALETTES> palettes_slots[2] = {};

  /**
   * @brief Spawns one burst of burst_size ripples from a random origin,
   *        staggered RIPPLE_STAGGER_FRAMES apart, each expanding over ripple_duration
   *        frames.
   * @param canvas Unused render target for the timer callback signature.
   */
  void ripple(Canvas &) {
    Vector origin = random_vector();
    for (int i = 0; i < (int)params.burst_size; i++) {
      if (!ripple_gen.spawn(i * RIPPLE_STAGGER_FRAMES, origin,
                            PI_F / ripple_duration,
                            static_cast<int>(ripple_duration)))
        hs::log("IslamicStars: ripple pool full, dropping spawn");
    }
  }

  /**
   * @brief Orients, ripple-distorts, and segue-shapes base_state, then
   *        rasterizes it with a per-face palette lookup.
   * @param canvas Render target receiving the rasterized mesh.
   * @param phase Segue phase in [0, 1] from the sprite envelope: rises over
   *        the incoming window, holds 1, falls over the outgoing window.
   * @param base_state Undistorted source mesh to transform and draw.
   * @param face_indices Maps each face to its topology class.
   * @param palette_idx Assigns a palette per topology class.
   * @note Draws on the exact SDF path, not the congruence-class LUT
   * (mesh_classes.h): ripple and the warping segues deform the mesh most
   * frames, and a canonical LUT under deformation either mis-shades the
   * interior gradient or pops when a face switches to the exact path. The
   * facility is for effects whose meshes hold still.
   */
  void draw_shape(Canvas &canvas, float phase, const MeshState &base_state,
                  const ArenaVector<int> &face_indices,
                  const std::array<int, NUM_PALETTES> &palette_idx) {
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase))
      return;
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state;
    OrientTransformer camera(orientation);
    MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                       ripple_gen, camera);

    const int *raw_indices = face_indices.data();
    const int num_faces = static_cast<int>(face_indices.size());

    // Per-face segues order faces by their center, recomputed per frame: from
    // world space by default (the front stays fixed in the room while the
    // mesh rotates through it), or from the untransformed mesh for segues
    // declaring LOCAL_SWEEP (the front rides the mesh). The third argument is
    // the face's palette-slot class, mapped exactly as the fragment shader
    // maps it; class-agnostic sweeps ignore it.
    constexpr bool PER_FACE =
        requires(const Vector &c) { seg.face_offset(c, 0, 0); };
    ArenaVector<float> face_offsets;
    ArenaVector<float> face_fades;
    if constexpr (PER_FACE) {
      constexpr bool LOCAL_SWEEP = requires { requires SegueT::LOCAL_SWEEP; };
      const MeshState &sweep_state =
          LOCAL_SWEEP ? base_state : transformed_state;
      const size_t faces = sweep_state.num_faces();
      face_offsets.bind(scratch_arena_a, faces);
      face_fades.bind(scratch_arena_a, faces);
      const uint16_t *fidx = sweep_state.get_faces_data();
      const uint16_t *foff = sweep_state.get_face_offsets_data();
      const uint8_t *fcnt = sweep_state.get_face_counts_data();
      for (size_t f = 0; f < faces; ++f) {
        Vector c(0.0f, 0.0f, 0.0f);
        for (int k = 0; k < fcnt[f]; ++k)
          c = c + sweep_state.vertices[fidx[foff[f] + k]];
        int cls = (f < static_cast<size_t>(num_faces))
                      ? wrap(raw_indices[f], NUM_PALETTES)
                      : 0;
        face_offsets.push_back(
            seg.face_offset(normalized_or(c, UP), static_cast<int>(f), cls));
        face_fades.push_back(seg.face_fade_frac(static_cast<int>(f)));
      }
    }

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      float p = phase;
      if constexpr (PER_FACE) {
        int fi = static_cast<int>(frag.v2);
        if (fi >= 0 && fi < static_cast<int>(face_offsets.size()))
          p = seg.face_phase(phase, face_offsets[fi], face_fades[fi]);
      }
      frag.color = shade_mesh_topology(frag, raw_indices, num_faces,
                                       palette_bank_, palette_idx, 1.0f, seg, p);
    };

    Scan::Mesh::draw<W, H>(filters, canvas, transformed_state, fragment_shader,
                           scratch_arena_a, params.debug_bb);
  }

  /**
   * @brief Advances to the next solid, generates it into the carousel's back
   *        slot with a freshly shuffled palette, makes it the front, schedules
   *        the segue and the shape's mid-display ripple burst, and queues the
   *        next spawn_shape call.
   */
  void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int back = 1 - carousel.front_index();
    MeshPaletteBank::shuffle_indices(palettes_slots[back]);

    int idx = solid_idx;

    auto draw_fn = [this, back](Canvas &canvas, float phase) {
      const MeshState &mesh = carousel.slot(back);
      this->draw_shape(canvas, phase, mesh, mesh.topology,
                       palettes_slots[back]);
    };

    // Compact the back slot, rebaking palettes into the fresh arena rather than
    // tracking them through the evacuation.
    carousel.compact_keep_front(
        [this](Arena &arena) { palette_bank_.bake_all(arena); });

    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh mesh = solids[idx].generate(a, b);
      carousel.slot(back).clear();
      MeshOps::compile(mesh, carousel.slot(back), target);
    });

    // ScratchScope frees only this call's own allocations, preserving prior
    // caller allocations in these shared arenas that a bare reset() would drop.
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(carousel.slot(back), scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Flip front eagerly for the overlapping sprite.
    carousel.set_front(back);

    // Segues with a spatial anchor (sweep axis, wave origin, spin axis) get a
    // fresh random one per transition. Safe mid-carousel: those segues are
    // sequential, so the previous sprite has already finished.
    if constexpr (requires(SegueT &s, const Vector &v) { s.retarget(v); })
      carousel.segue().retarget(random_vector());

    // Per-shape choreography: segue in, hold still one second, ripple, settle
    // one second, segue out. Duration is derived from the stage lengths so the
    // stages never overlap; the segue warps are identity on the phase-1
    // plateau, so the mesh only moves during its own stage.
    int fade = static_cast<int>(params.fade);
    int burst_span =
        (static_cast<int>(params.burst_size) - 1) * RIPPLE_STAGGER_FRAMES +
        static_cast<int>(ripple_duration);
    int duration = fade + STILL_FRAMES + burst_span + STILL_FRAMES + fade;

    int next_delay = carousel.schedule_segue(timeline, draw_fn, duration, fade);

    timeline.add(fade + STILL_FRAMES,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &canvas) { ripple(canvas); }, false));

    const auto &entry = solids[solid_idx];
    hs::log("Spawning Shape: %s (V=%d, F=%d)", entry.name,
            (int)carousel.current().vertices.size(),
            (int)carousel.current().faces.size());

    // The segue decides when the next shape starts relative to this one.
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  /**
   * @brief Slider-backed runtime parameters for the effect.
   */
  struct Params {
    float fade = 72.0f; /**< Segue window length, in frames: a 64-frame (4 s) sweep crossing plus one per-face fade tail. */
    float burst_size = 4.0f; /**< Ripples per burst; float-backed for register_param. */
    bool debug_bb = false; /**< Whether to draw mesh bounding boxes. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
