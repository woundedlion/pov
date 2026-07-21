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
  HS_COLD_MEMBER IslamicStars() : Effect(W, H, {.strobe = true}), filters(), ripple_gen(timeline) {}

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

    ripple_gen.init_storage(persistent_arena);
    palette_bank_.bake_all(persistent_arena);

    // Set BEFORE registering: register_param snaps *ptr as the slider default.
    // Amplitude starts at the fold-free ceiling (see RIPPLE_AMP_MAX).
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
    register_param("Trans Speed", &params.trans_speed, 1.0f, 8.0f);
    register_param("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    spawn_shape();
  }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(is_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(is_ripple_prepare);
      ripple_gen.prepare_frame();
    }
    {
      HS_PROFILE(is_timeline_step);
      timeline.step(canvas);
    }
  }

private:
  // Ripple-pool sizing: a slot is held from spawn() until the staggered ripple
  // completes. Only one burst is normally live; the pool holds two so a Ripp
  // Dur/Burst slider change mid-burst cannot drop a spawn.
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
  // Effective per-shape stage lengths after the Trans Speed divisor, written by
  // spawn_shape and read by the deferred ripple() callback.
  int ripple_dur_eff_ = 80;
  int ripple_stagger_eff_ = RIPPLE_STAGGER_FRAMES;
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
      if (!ripple_gen.spawn(i * ripple_stagger_eff_, origin,
                            PI_F / ripple_dur_eff_, ripple_dur_eff_))
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
   * (mesh_classes.h): ripple/segue deformation makes a canonical LUT mis-shade
   * or pop. The facility is for effects whose meshes hold still.
   */
  HS_O3_FN void draw_shape(Canvas &canvas, float phase,
                           const MeshState &base_state,
                           const ArenaVector<int> &face_indices,
                           const std::array<int, NUM_PALETTES> &palette_idx) {
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase))
      return;
    HS_PROFILE(is_draw_shape);
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state;
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
      MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                         ripple_gen, camera);
    }

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
    // phase is fixed for the whole draw call, so the segue's per-face phase
    // resolves here rather than per fragment; so does the face's palette, whose
    // slot is the same class that face_offset already needs.
    ArenaVector<float> face_phases;
    ArenaVector<const BakedPalette *> face_palettes;
    if constexpr (PER_FACE) {
      HS_PROFILE(is_face_offsets);
      constexpr bool LOCAL_SWEEP = requires { requires SegueT::LOCAL_SWEEP; };
      const MeshState &sweep_state =
          LOCAL_SWEEP ? base_state : transformed_state;
      const size_t faces = sweep_state.num_faces();
      face_phases.bind(scratch_arena_a, faces);
      face_palettes.bind(scratch_arena_a, faces);
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
        float off =
            seg.face_offset(normalized_or(c, UP), static_cast<int>(f), cls);
        float fade = seg.face_fade_frac(static_cast<int>(f));
        face_phases.push_back(seg.face_phase(phase, off, fade));
        face_palettes.push_back(&palette_bank_[palette_idx[cls]]);
      }
    }

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      if constexpr (PER_FACE) {
        int fi = static_cast<int>(frag.v2);
        if (fi >= 0 && fi < static_cast<int>(face_phases.size())) {
          frag.color = shade_mesh_topology(frag, *face_palettes[fi], 1.0f, seg,
                                           face_phases[fi]);
          return;
        }
      }
      frag.color = shade_mesh_topology(frag, raw_indices, num_faces,
                                       palette_bank_, palette_idx, 1.0f, seg,
                                       phase);
    };

    {
      HS_PROFILE(is_mesh_scan);
      Scan::Mesh::draw<W, H>(filters, canvas, transformed_state,
                             fragment_shader, scratch_arena_a, params.debug_bb);
    }
  }

  /**
   * @brief Advances to the next solid, generates it into the carousel's back
   *        slot with a freshly shuffled palette, makes it the front, schedules
   *        the segue and the shape's mid-display ripple burst, and queues the
   *        next spawn_shape call.
   */
  HS_COLD_MEMBER void spawn_shape() {
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
    // tracking them through the evacuation. The ripple pool re-claims first so
    // its slots re-land at their init_storage() addresses.
    carousel.compact_keep_front([this](Arena &arena) {
      ripple_gen.reclaim_storage(arena);
      palette_bank_.bake_all(arena);
    });

    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh mesh = solids[idx].generate(a, b);
      carousel.slot(back).clear();
      MeshOps::compile(mesh, carousel.slot(back), target, scratch_arena_a);
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
    // Trans Speed divides every stage length so the carousel can be sped up
    // (e.g. for profiling) without touching the shape geometry. Each stage keeps
    // a >=1-frame floor. The effective ripple duration/stagger are cached for the
    // deferred ripple() callback, which fires before the next shape spawns.
    const float sp = std::max(1.0f, params.trans_speed);
    int fade = std::max(1, static_cast<int>(params.fade / sp));
    int still = std::max(1, static_cast<int>(STILL_FRAMES / sp));
    ripple_dur_eff_ = std::max(8, static_cast<int>(ripple_duration / sp));
    ripple_stagger_eff_ =
        std::max(1, static_cast<int>(RIPPLE_STAGGER_FRAMES / sp));
    int burst_span =
        (static_cast<int>(params.burst_size) - 1) * ripple_stagger_eff_ +
        ripple_dur_eff_;
    int duration = fade + still + burst_span + still + fade;

    int next_delay = carousel.schedule_segue(timeline, draw_fn, duration, fade);

    timeline.add(fade + still,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &canvas) { ripple(canvas); }, false));

    // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
    const auto &entry = solids[solid_idx];
    const MeshState &spawned = carousel.current();
    hs::log("Spawning Shape: %s (V=%d, E=%d, F=%d, I=%d)", entry.name,
            (int)spawned.vertices.size(), (int)(spawned.faces.size() / 2),
            (int)spawned.face_counts.size(), (int)spawned.faces.size());

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
    float trans_speed = 1.0f; /**< Divides every per-shape stage length (fade, still holds, ripple span): 1 = shipping cadence, higher cycles shapes faster. */
    bool debug_bb = false; /**< Whether to draw mesh bounding boxes. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
