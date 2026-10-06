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

#include "core/animation/orientation.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/render/filter/pixel_feedback.h"

namespace hs_test {
namespace effects_tests {
struct MeshFeedbackWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/** @brief MeshFeedback's parameter set: the wireframe solid and the feedback
 *  style rendering it. */
struct MeshFeedbackParams {
  Solids::BaseMesh base_mesh = Solids::BaseMesh::ICOSAHEDRON;
  Feedback::Style style = Feedback::Style::ArcingLightning();
};

/**
 * @brief Feedback effect over a selectable polyhedron.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Draws the solid's wireframe under an orientation random-walk while
 * the feedback filter warps and fades the accumulated frame. Style presets
 * switch immediately at fixed intervals while mesh emission continues. The
 * base mesh changes with the preset or live selector.
 */
template <int W, int H>
class MeshFeedback
    : public ChoreographedEffect<MeshFeedback<W, H>, MeshFeedbackParams> {
public:
  static constexpr const char *EFFECT_ID = "MeshFeedback";

  using Choreography =
      ChoreographedEffect<MeshFeedback<W, H>, MeshFeedbackParams>;
  using Params = MeshFeedbackParams;
  using Style = Feedback::Style;
  using BaseMesh = Solids::BaseMesh;

  /** Snap: Style embeds a noise binding and base_mesh rewinds the mesh arena,
      so parameters cannot blend. */
  static constexpr Segue::Preset::Snap DEPARTURE{};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 241;

  // Persistent gamut boundary bracket grid at the flash master's resolution.
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  static constexpr float FADE_MIN = 0.0f, FADE_MAX = 0.99f;
  static constexpr float AMP_MIN = 0.0f, AMP_MAX = 30.0f;
  static constexpr float FREQ_MIN = 0.01f, FREQ_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float SCALE_MIN = 0.1f, SCALE_MAX = 50.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 0.5f;
  static constexpr float POLE_RES_MIN = 0.0f, POLE_RES_MAX = 2.0f;

  /** @brief Shared registration, validation and interpolation descriptions. */
  static constexpr auto parameter_fields() {
    return std::tuple{
        Control::Field<Params, BaseMesh>{
            .id = "base_mesh",
            .member = &Params::base_mesh,
            .name = "Base Mesh",
            .spec = {.min = 0,
                     .max = static_cast<int64_t>(Solids::BASE_MESH_COUNT) - 1,
                     .animated = true,
                     .options = Solids::BASE_MESH_OPTIONS,
                     .export_options = Solids::BASE_MESH_EXPORT_OPTIONS,
                     .option_count = Solids::BASE_MESH_COUNT}},
        Control::FieldGroup{
            &Params::style,
            std::tuple{
                Control::Field<Style, float>{.id = "fade",
                                             .member = &Style::fade,
                                             .name = "Fade",
                                             .spec = {.min = FADE_MIN,
                                                      .max = FADE_MAX,
                                                      .animated = true}},
                Control::Field<Style, float>{
                    .id = "amplitude",
                    .member = &Style::amplitude,
                    .name = "Distort Amp",
                    .spec = {.min = AMP_MIN, .max = AMP_MAX, .animated = true}},
                Control::Field<Style, float>{.id = "frequency",
                                             .member = &Style::frequency,
                                             .name = "Distort Freq",
                                             .spec = {.min = FREQ_MIN,
                                                      .max = FREQ_MAX,
                                                      .animated = true}},
                Control::Field<Style, float>{.id = "speed",
                                             .member = &Style::speed,
                                             .name = "Distort Speed",
                                             .spec = {.min = SPEED_MIN,
                                                      .max = SPEED_MAX,
                                                      .animated = true}},
                Control::Field<Style, float>{.id = "scale",
                                             .member = &Style::scale,
                                             .name = "Noise Scale",
                                             .spec = {.min = SCALE_MIN,
                                                      .max = SCALE_MAX,
                                                      .animated = true}},
                Control::Field<Style, float>{.id = "hue_shift",
                                             .member = &Style::hue_shift,
                                             .name = "Hue Shift",
                                             .spec = {.min = HUE_SHIFT_MIN,
                                                      .max = HUE_SHIFT_MAX,
                                                      .animated = true}}}}};
  }

  static constexpr bool preset_in_ranges(const Params &p) {
    return Control::valid_fields(p, parameter_fields());
  }
  static constexpr size_t PRESET_COUNT = 12;
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{BaseMesh::ICOSAHEDRON, Style::ArcingLightning()}, DEPARTURE},
      {{BaseMesh::DODECAHEDRON, Style::SlowFire()}, DEPARTURE},
      {{BaseMesh::TRUNCATED_ICOSAHEDRON, Style::EnergeticFire()}, DEPARTURE},
      {{BaseMesh::CUBE, Style::Smoke()}, DEPARTURE},
      {{BaseMesh::RHOMBICUBOCTAHEDRON, Style::SlowDust()}, DEPARTURE},
      {{BaseMesh::ICOSIDODECAHEDRON, Style::WavyTrails()}, DEPARTURE},
      {{BaseMesh::OCTAHEDRON, Style::MeltingHi()}, DEPARTURE},
      {{BaseMesh::TRIAKIS_OCTAHEDRON, Style::MeltingLo()}, DEPARTURE},
      {{BaseMesh::RHOMBIC_TRIACONTAHEDRON, Style::Miasma()}, DEPARTURE},
      {{BaseMesh::TRUNCATED_CUBOCTAHEDRON, Style::LooseWormhole()}, DEPARTURE},
      {{BaseMesh::SNUB_CUBE, Style::TightWormhole()}, DEPARTURE},
      {{BaseMesh::PENTAGONAL_HEXECONTAHEDRON, Style::WigglingWormhole()},
       DEPARTURE},
  }};

  static_assert(
      all_presets_in_ranges(PRESETS, preset_in_ranges),
      "a MeshFeedback preset drives a style field outside its "
      "registered slider range; widen the range to accommodate the "
      "preset (the range exposes the presets, it does not clamp them)");

  /** @brief Startup parameters: preset 0 with the noise binding still null;
   *  init() binds the effect-owned NoiseParams. */
  static Params initial_params() { return PRESETS[0].params; }

  /**
   * @brief Wires up noise, orientation, and the filter pipeline.
   * @details The Feedback filter binds `params.style` by reference; the
   * Choreography base constructs `params` first.
   */
  HS_COLD_MEMBER MeshFeedback()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        noise_params(), orientation(),
        filters(Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H>(params.style)) {}

  /**
   * @brief One-time effect setup.
   */
  HS_COLD_MEMBER void init() override {
    begin_choreography();

    // Configure the noise type before apply_params() calls sync_noise().
    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.set_seed(hs::rand_int(0, 65536));
    noise_params.sync();

    // initial_params() cannot bind the noise pointer; adopt_params() does.
    adopt_params(params);

    mesh_shade = Palettes::PEACH_POP.get(0.0f);

    this->register_described_params();

    register_param("Feedback", &feedback_enabled);
    mark_global("Feedback");
    register_param("Pole Half-Res", &pole_half_res, POLE_RES_MIN, POLE_RES_MAX);
    mark_global("Pole Half-Res");

    filters.init_storage(persistent_arena);
    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    mesh_storage_mark = persistent_arena.get_offset();
    apply_params();

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, math::Y_AXIS, walk_noise));
  }

  /**
   * @brief Renders one frame.
   * @details The preset step precedes apply_params() so the flush reads one
   * preset's style; the flush precedes the mesh draw so it does not decay this
   * frame's wireframe.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    step_choreography();

    {
      HS_PROFILE(mf_apply_params);
      apply_params();
    }

    {
      HS_PROFILE(mf_timeline_step);
      timeline.step(canvas);
    }

    auto &frame_filters = [&]() -> decltype(auto) {
      HS_PROFILE(mf_feedback_flush);
      return filters.begin_frame(canvas, 1.0f);
    }();

    {
      HS_PROFILE(mf_mesh_draw);
      const Color4 shade = mesh_shade;
      Plot::Mesh::draw<W, H>(
          frame_filters, canvas, mesh, edges,
          [&](const math::Vector &, Fragment &f) { f.color = shade; });
    }
  }

private:
  friend Choreography;
  friend struct ::hs_test::effects_tests::MeshFeedbackWhiteBox;

  using Choreography::begin_choreography;
  using Choreography::mark_global;
  using Choreography::params;
  using Choreography::register_param;
  using Choreography::step_choreography;
  using Choreography::timeline;

  /** @brief Params for the preset at @p index, with the effect-owned noise
   *  bound into the style. */
  Params preset_params(size_t index) {
    Params p = PRESETS[index].params;
    p.style.noise = &noise_params;
    return p;
  }

  /** @brief Adopts a snap or snapshot target, re-binding the noise pointer a
   *  foreign snapshot would otherwise carry stale. */
  void adopt_params(const Params &target) {
    params = target;
    params.style.noise = &noise_params;
  }

  static_assert(Solids::MAX_SOLID_VERTICES <=
                    static_cast<size_t>(Plot::Mesh::DEDUP_CAPACITY),
                "a bound solid must fit the wireframe edge-dedup bitset");

  HS_FLASH_MEMBER void rebuild_mesh(BaseMesh base_mesh) {
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
   * @brief Pushes UI-tunable state into the mesh, style and filters.
   */
  void apply_params() {
    if (!mesh_ready || params.base_mesh != active_base_mesh)
      rebuild_mesh(params.base_mesh);
    params.style.sync_noise();
    params.style.pole_half_res = pole_half_res;
    filters.template get<Filter::Pixel::Feedback<W, H>>().set_enabled(
        feedback_enabled);
  }

  bool feedback_enabled = true;
  /** Feedback::Style::pole_half_res; global, so preset snaps keep it. */
  float pole_half_res = 1.0f;
  Animation::NoiseParams noise_params;
  // Separate walk generator: RandomWalk owns its frequency and seed, and
  // noise_params.sync() rewrites the shared one every frame.
  FastNoiseLite walk_noise;

  math::Orientation<> orientation;

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

  static constexpr size_t SCRATCH_A_PEAK_BYTES =
      Plot::Mesh::EDGE_MAX_POINTS * sizeof(Fragment) +
      Plot::rasterize_scratch_a_bytes<W>();
  static_assert(
      SCRATCH_A_PEAK_BYTES <= DEFAULT_SCRATCH_A_SIZE,
      "MeshFeedback wireframe draw exceeds the default scratch_a budget");
  static_assert(
      sizeof(TriangularBitset<Plot::Mesh::DEDUP_CAPACITY>) <=
          DEFAULT_SCRATCH_B_SIZE,
      "MeshFeedback edge extraction exceeds the default scratch_b budget");

  static constexpr size_t MESH_STORAGE_BYTES =
      Solids::MAX_SOLID_VERTICES * sizeof(math::Vector) +
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
