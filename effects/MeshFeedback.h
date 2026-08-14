/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file MeshFeedback.h
 * @brief Feedback filter over a selectable polyhedral wireframe.
 */

#include "core/engine/engine.h"
#include "core/render/filter_feedback.h"

// Unit-test accessor reaching the private style/noise/preset bookkeeping.
namespace hs_test {
namespace effects_tests {
struct MeshFeedbackWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Feedback effect over a selectable polyhedron.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Draws the solid's wireframe under an orientation random-walk while
 * the feedback filter warps and fades the accumulated frame. Style presets
 * switch immediately at fixed intervals while mesh emission continues. The
 * base mesh changes with the preset or live selector.
 */
template <int W, int H> class MeshFeedback : public Effect {
public:
  using Style = Feedback::Style;
  using BaseMesh = Solids::BaseMesh;

  struct Params {
    BaseMesh base_mesh = BaseMesh::ICOSAHEDRON;
    Style style = Style::ArcingLightning();
  };

  static constexpr int PRESET_FRAMES = 241;

  // Gamut boundary bracket grid bought from the persistent arena (131,072 B),
  // at the flash master's own resolution, so the copy is verbatim.
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  static constexpr float FADE_MIN = 0.0f, FADE_MAX = 0.99f;
  static constexpr float AMP_MIN = 0.0f, AMP_MAX = 30.0f;
  static constexpr float FREQ_MIN = 0.01f, FREQ_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float SCALE_MIN = 0.1f, SCALE_MAX = 50.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 0.5f;

  /** @brief True iff every preset-driven field of @p s lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Params &p) {
    const Style &s = p.style;
    return static_cast<size_t>(p.base_mesh) < Solids::BASE_MESH_COUNT &&
           s.fade >= FADE_MIN && s.fade <= FADE_MAX && s.amplitude >= AMP_MIN &&
           s.amplitude <= AMP_MAX && s.frequency >= FREQ_MIN &&
           s.frequency <= FREQ_MAX && s.speed >= SPEED_MIN &&
           s.speed <= SPEED_MAX && s.scale >= SCALE_MIN &&
           s.scale <= SCALE_MAX && s.hue_shift >= HUE_SHIFT_MIN &&
           s.hue_shift <= HUE_SHIFT_MAX;
  }
  static constexpr size_t PRESET_COUNT = 12;
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{BaseMesh::ICOSAHEDRON, Style::ArcingLightning()}},
      {{BaseMesh::DODECAHEDRON, Style::SlowFire()}},
      {{BaseMesh::TRUNCATED_ICOSAHEDRON, Style::EnergeticFire()}},
      {{BaseMesh::CUBE, Style::Smoke()}},
      {{BaseMesh::RHOMBICUBOCTAHEDRON, Style::SlowDust()}},
      {{BaseMesh::ICOSIDODECAHEDRON, Style::WavyTrails()}},
      {{BaseMesh::OCTAHEDRON, Style::MeltingHi()}},
      {{BaseMesh::TRIAKIS_OCTAHEDRON, Style::MeltingLo()}},
      {{BaseMesh::RHOMBIC_TRIACONTAHEDRON, Style::Miasma()}},
      {{BaseMesh::TRUNCATED_CUBOCTAHEDRON, Style::LooseWormhole()}},
      {{BaseMesh::SNUB_CUBE, Style::TightWormhole()}},
      {{BaseMesh::PENTAGONAL_HEXECONTAHEDRON, Style::WigglingWormhole()}},
  }};

  static_assert(
      all_presets_in_ranges(PRESETS, preset_in_ranges),
      "a MeshFeedback preset drives a style field outside its "
      "registered slider range; widen the range to accommodate the "
      "preset (the range exposes the presets, it does not clamp them)");

  /**
   * @brief Wires up noise, orientation, and the filter pipeline.
   * @details The Feedback filter binds `params.style` by reference; keep
   * `params` declared before `filters` so it is constructed first.
   */
  HS_COLD_MEMBER MeshFeedback()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        noise_params(), orientation(), timeline(),
        filters(Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H>(params.style)) {}

  /**
   * @brief One-time effect setup.
   * @details Binds shared noise into presets, builds the selected mesh, samples
   * the mesh shade, registers tunable params, and schedules the
   * noise/walk/preset timers.
   */
  void init() override {
    configure_presets(PRESETS.size());
    for (auto &e : presets.get_entries()) {
      e.params.style.noise = &noise_params;
    }

    // Configure the noise type before apply_params(): it calls sync_noise(),
    // which would otherwise propagate the default noise type on the first frame.
    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.noise.SetSeed(hs::rand_int(0, 65536));
    noise_params.sync();

    params = presets.get();

    mesh_shade = Palettes::PEACH_POP.get(0.0f);

    register_animated_param(
        "Base Mesh", &params.base_mesh, Solids::BASE_MESH_OPTIONS,
        Solids::BASE_MESH_EXPORT_OPTIONS, Solids::BASE_MESH_COUNT);
    register_animated_param("Fade", &params.style.fade, FADE_MIN, FADE_MAX);
    register_animated_param("Distort Amp", &params.style.amplitude, AMP_MIN,
                            AMP_MAX);
    register_animated_param("Distort Freq", &params.style.frequency, FREQ_MIN,
                            FREQ_MAX);
    register_animated_param("Distort Speed", &params.style.speed, SPEED_MIN,
                            SPEED_MAX);
    register_animated_param("Noise Scale", &params.style.scale, SCALE_MIN,
                            SCALE_MAX);
    register_animated_param("Hue Shift", &params.style.hue_shift, HUE_SHIFT_MIN,
                            HUE_SHIFT_MAX);
    register_param("Feedback", &feedback_enabled);
    mark_global("Feedback");

    filters.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
        persistent_arena);
    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    mesh_storage_mark = persistent_arena.get_offset();
    apply_params();

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS, walk_noise));
  }

  /**
   * @brief Renders one frame.
   * @details Advances the preset selection, applies params, steps the timeline,
   * runs the feedback decay flush, then draws the mesh. The preset switch leads
   * apply_params() so the noise scalars and the fade/hue the flush reads come
   * from the same preset, and the flush leads the mesh draw so this frame's
   * wireframe is not decayed by its own flush.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    advance_preset();

    {
      HS_PROFILE(mf_apply_params);
      apply_params();
    }

    {
      HS_PROFILE(mf_timeline_step);
      timeline.step(canvas);
    }

    {
      // Feedback-buffer warp/tap + decay flush.
      HS_PROFILE(mf_feedback_flush);
      filters.flush(canvas, 1.0f);
    }

    {
      HS_PROFILE(mf_mesh_draw);
      const Color4 shade = mesh_shade;
      Plot::Mesh::draw<W, H>(
          filters, canvas, mesh, edges,
          [&](const Vector &, Fragment &f) { f.color = shade; });
    }
  }

private:
  HS_FLASH_MEMBER bool apply_preset(const PresetChange &change) override {
    if (!presets.select(change.to))
      return false;
    presets.apply(params);
    preset_frames = 0;
    return true;
  }

  friend struct ::hs_test::effects_tests::MeshFeedbackWhiteBox;

  static_assert(Solids::MAX_SOLID_VERTICES <=
                    static_cast<size_t>(Plot::Mesh::DEDUP_CAPACITY),
                "a bound solid must fit the wireframe edge-dedup bitset");

  void rebuild_mesh(BaseMesh base_mesh) {
    mesh = MeshState();
    edges = ArenaVector<Plot::Mesh::Edge>();
    persistent_arena.set_offset(mesh_storage_mark);

    const Solids::Entry &entry =
        Solids::get_entry(static_cast<size_t>(base_mesh));
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh poly = entry.generate(a, b);
      HS_CHECK(poly.vertices.size() <= Solids::MAX_SOLID_VERTICES &&
                   poly.faces.size() <= Solids::MAX_SOLID_FACE_SLOTS &&
                   poly.face_counts.size() <= Solids::MAX_SOLID_FACES,
               "MeshFeedback selectable solid exceeds bounds");
      MeshOps::compile(poly, mesh, target, a);
    });

    edges.bind(persistent_arena, Solids::MAX_SOLID_EDGES);
    Plot::Mesh::extract_edges(mesh, edges);
    active_base_mesh = base_mesh;
    mesh_ready = true;
  }

  /**
   * @brief Pushes UI-tunable state into the live style/filters each frame.
   * @details Refreshes the noise binding and toggles the feedback filter from
   * `feedback_enabled`.
   */
  void apply_params() {
    if (!mesh_ready || params.base_mesh != active_base_mesh)
      rebuild_mesh(params.base_mesh);
    params.style.sync_noise();
    filters.template get<Filter::Pixel::Feedback<W, H>>().set_enabled(
        feedback_enabled);
  }

  /**
   * @brief Advances to the next preset every PRESET_FRAMES.
   * @details Frozen while animations are paused.
   */
  void advance_preset() {
    if (animations_paused())
      return;
    if (++preset_frames < PRESET_FRAMES)
      return;
    preset_frames = 0;
    HS_CHECK(advancePreset());
  }

  Params params;
  Presets<Params, PRESET_COUNT> presets{PRESETS};
  bool feedback_enabled = true;
  int preset_frames = 0;
  Animation::NoiseParams noise_params;
  // Dedicated walk generator: RandomWalk's ctor takes exclusive control of its
  // frequency and seed, and noise_params.sync() runs every frame off the
  // "Distort Freq" slider.
  FastNoiseLite walk_noise;

  Orientation<> orientation;
  Timeline timeline;

  Color4 mesh_shade; /**< Wireframe shade; sampled once in init(). */

  MeshState mesh;
  ArenaVector<Plot::Mesh::Edge>
      edges; /**< Unique edge list (topology is static). */
  size_t mesh_storage_mark = 0;
  BaseMesh active_base_mesh = BaseMesh::ICOSAHEDRON;
  bool mesh_ready = false;

  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H>>
      filters;

  static constexpr size_t MESH_STORAGE_BYTES =
      Solids::MAX_SOLID_VERTICES * sizeof(Vector) +
      Solids::MAX_SOLID_FACES * (sizeof(uint8_t) + sizeof(uint16_t)) +
      Solids::MAX_SOLID_FACE_SLOTS * sizeof(uint16_t) +
      Solids::MAX_SOLID_EDGES * sizeof(Plot::Mesh::Edge);
  static_assert(Filter::Pixel::Feedback<W, H>::STORAGE_BYTES +
                        gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
                        MESH_STORAGE_BYTES <=
                    DEVICE_PERSISTENT_BUDGET,
                "MeshFeedback persistent storage exceeds the "
                "default persistent partition; retune the feedback downsample, "
                "coarsen the gamut grid, or carve arenas");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
