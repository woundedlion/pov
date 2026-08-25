/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Raymarch.h
 * @brief Twisted-torus SDFs ray-marched at every vertex of a selectable solid.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/color/noise_hue_palette.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/sdf/volume.h"

// Unit-test accessor reaching the private torus proportions, so a test can pin
// the cull sphere against the warped SDF it has to contain.
namespace hs_test {
namespace effects_tests {
struct RaymarchWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/** Solids whose vertex counts fit Raymarch's copy and timeline capacity. */
enum class RaymarchPlacementSolid : uint8_t {
  TETRAHEDRON,
  CUBE,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  TRUNCATED_TETRAHEDRON,
  CUBOCTAHEDRON,
  TRUNCATED_CUBE,
  TRUNCATED_OCTAHEDRON,
  RHOMBICUBOCTAHEDRON,
  SNUB_CUBE,
  ICOSIDODECAHEDRON,
  TRIAKIS_TETRAHEDRON,
  RHOMBIC_DODECAHEDRON,
  TRIAKIS_OCTAHEDRON,
  TETRAKIS_HEXAHEDRON,
  DELTOIDAL_ICOSITETRAHEDRON,
  DISDYAKIS_DODECAHEDRON,
  RHOMBIC_TRIACONTAHEDRON,
  TRIAKIS_ICOSAHEDRON,
  PENTAKIS_DODECAHEDRON,
  COUNT
};

/** @brief Raymarch preset and live-control state. */
struct RaymarchParams {
  RaymarchPlacementSolid base_solid =
      RaymarchPlacementSolid::DISDYAKIS_DODECAHEDRON;
  float pulse_speed = 5.0f;
  float fill = 0.75f;
  uint8_t max_steps = 18;
  float diffuse = 0.4f;
  float specular = 1.2f;
  float fresnel = 0.2f;
  float twist = 2.0f;
  float aa_mult = 0.5f;
  float hue_shift = 0.76f;
  float hue_noise_scale = 0.3f;
  float hue_noise_speed = 0.0002f;
};

/**
 * @brief Ray-marches a twisted torus SDF at each vertex of a selectable solid,
 *        shading each with a metallic headlight model and a baked OKLCH palette
 *        under a seamless torus-UV noise hue field. Each torus has an independent
 *        random-walk tumble and is auto-sized to its own nearest-neighbour gap
 *        (scaled by the live Fill param), so they pack without overlap.
 * @tparam W Effect render width in pixels.
 * @tparam H Effect render height in pixels.
 */
template <int W, int H>
class Raymarch : public ChoreographedEffect<Raymarch<W, H>, RaymarchParams> {
public:
  using Choreography = ChoreographedEffect<Raymarch<W, H>, RaymarchParams>;
  using Params = RaymarchParams;
  using PlacementSolid = RaymarchPlacementSolid;

  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "uv-surface-noise"};
  static constexpr Segue::Snap PRESET_SEGUE{};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  static constexpr int MAX_POINTS = 32;
  static constexpr size_t PLACEMENT_SOLID_COUNT =
      static_cast<size_t>(PlacementSolid::COUNT);

  static constexpr const char *PLACEMENT_SOLID_OPTIONS[] = {
      "Tetrahedron",
      "Cube",
      "Octahedron",
      "Dodecahedron",
      "Icosahedron",
      "Truncated Tetrahedron",
      "Cuboctahedron",
      "Truncated Cube",
      "Truncated Octahedron",
      "Rhombicuboctahedron",
      "Snub Cube",
      "Icosidodecahedron",
      "Triakis Tetrahedron",
      "Rhombic Dodecahedron",
      "Triakis Octahedron",
      "Tetrakis Hexahedron",
      "Deltoidal Icositetrahedron",
      "Disdyakis Dodecahedron",
      "Rhombic Triacontahedron",
      "Triakis Icosahedron",
      "Pentakis Dodecahedron"};

  static constexpr const char *PLACEMENT_SOLID_EXPORT_OPTIONS[] = {
      "RaymarchPlacementSolid::TETRAHEDRON",
      "RaymarchPlacementSolid::CUBE",
      "RaymarchPlacementSolid::OCTAHEDRON",
      "RaymarchPlacementSolid::DODECAHEDRON",
      "RaymarchPlacementSolid::ICOSAHEDRON",
      "RaymarchPlacementSolid::TRUNCATED_TETRAHEDRON",
      "RaymarchPlacementSolid::CUBOCTAHEDRON",
      "RaymarchPlacementSolid::TRUNCATED_CUBE",
      "RaymarchPlacementSolid::TRUNCATED_OCTAHEDRON",
      "RaymarchPlacementSolid::RHOMBICUBOCTAHEDRON",
      "RaymarchPlacementSolid::SNUB_CUBE",
      "RaymarchPlacementSolid::ICOSIDODECAHEDRON",
      "RaymarchPlacementSolid::TRIAKIS_TETRAHEDRON",
      "RaymarchPlacementSolid::RHOMBIC_DODECAHEDRON",
      "RaymarchPlacementSolid::TRIAKIS_OCTAHEDRON",
      "RaymarchPlacementSolid::TETRAKIS_HEXAHEDRON",
      "RaymarchPlacementSolid::DELTOIDAL_ICOSITETRAHEDRON",
      "RaymarchPlacementSolid::DISDYAKIS_DODECAHEDRON",
      "RaymarchPlacementSolid::RHOMBIC_TRIACONTAHEDRON",
      "RaymarchPlacementSolid::TRIAKIS_ICOSAHEDRON",
      "RaymarchPlacementSolid::PENTAKIS_DODECAHEDRON"};

  static constexpr std::array<Solids::BaseMesh, PLACEMENT_SOLID_COUNT>
      PLACEMENT_SOLIDS{Solids::BaseMesh::TETRAHEDRON,
                       Solids::BaseMesh::CUBE,
                       Solids::BaseMesh::OCTAHEDRON,
                       Solids::BaseMesh::DODECAHEDRON,
                       Solids::BaseMesh::ICOSAHEDRON,
                       Solids::BaseMesh::TRUNCATED_TETRAHEDRON,
                       Solids::BaseMesh::CUBOCTAHEDRON,
                       Solids::BaseMesh::TRUNCATED_CUBE,
                       Solids::BaseMesh::TRUNCATED_OCTAHEDRON,
                       Solids::BaseMesh::RHOMBICUBOCTAHEDRON,
                       Solids::BaseMesh::SNUB_CUBE,
                       Solids::BaseMesh::ICOSIDODECAHEDRON,
                       Solids::BaseMesh::TRIAKIS_TETRAHEDRON,
                       Solids::BaseMesh::RHOMBIC_DODECAHEDRON,
                       Solids::BaseMesh::TRIAKIS_OCTAHEDRON,
                       Solids::BaseMesh::TETRAKIS_HEXAHEDRON,
                       Solids::BaseMesh::DELTOIDAL_ICOSITETRAHEDRON,
                       Solids::BaseMesh::DISDYAKIS_DODECAHEDRON,
                       Solids::BaseMesh::RHOMBIC_TRIACONTAHEDRON,
                       Solids::BaseMesh::TRIAKIS_ICOSAHEDRON,
                       Solids::BaseMesh::PENTAKIS_DODECAHEDRON};

  static_assert(std::size(PLACEMENT_SOLID_OPTIONS) == PLACEMENT_SOLID_COUNT);
  static_assert(std::size(PLACEMENT_SOLID_EXPORT_OPTIONS) ==
                PLACEMENT_SOLID_COUNT);

  static constexpr Params initial_params() { return {}; }

  static constexpr bool valid_params(const Params &value) {
    return static_cast<size_t>(value.base_solid) < PLACEMENT_SOLID_COUNT &&
           value.pulse_speed >= 0.0f && value.pulse_speed <= 10.0f &&
           value.fill >= 0.3f && value.fill <= 1.3f && value.max_steps >= 4 &&
           value.max_steps <= 30 && value.diffuse >= 0.0f &&
           value.diffuse <= 1.0f && value.specular >= 0.0f &&
           value.specular <= 1.5f && value.fresnel >= 0.0f &&
           value.fresnel <= 1.0f && value.twist >= 0.0f &&
           value.twist <= 8.0f && value.aa_mult >= 0.1f &&
           value.aa_mult <= 1.5f && value.hue_shift >= -4.0f &&
           value.hue_shift <= 4.0f && value.hue_noise_scale >= 1.0f / 64.0f &&
           value.hue_noise_scale <= 8.0f && value.hue_noise_speed >= -0.001f &&
           value.hue_noise_speed <= 0.001f;
  }

  /**
   * @brief Constructs the effect at the templated render dimensions.
   */
  HS_COLD_MEMBER Raymarch()
      : Choreography(W, H,
                     pipeline_config<decltype(pipeline)>({.strobe = true})) {}

  /**
   * @brief Registers tunable params, builds the placement-solid vertices,
   *        bakes the palette LUT, and installs the camera walk and phase drivers
   *        on the timeline.
   */
  HS_COLD_MEMBER void init() override {
    begin_choreography();
    register_animated_param(
        "Base Solid", &params.base_solid, PLACEMENT_SOLID_OPTIONS,
        PLACEMENT_SOLID_EXPORT_OPTIONS, PLACEMENT_SOLID_COUNT);
    register_param("Pulse Speed", &params.pulse_speed, 0.0f, 10.0f);
    // Fraction of the half nearest-neighbour gap the ring's outer edge reaches:
    // < 1 leaves a gap, 1 makes neighbours touch, > 1 overlaps them deliberately.
    register_param("Fill", &params.fill, 0.3f, 1.3f);
    register_int_param("Max Steps", &params.max_steps, 4, 30);
    register_param("Diffuse", &params.diffuse, 0.0f, 1.0f);
    register_param("Specular", &params.specular, 0.0f, 1.5f);
    register_param("Fresnel", &params.fresnel, 0.0f, 1.0f);
    register_param("Twist", &params.twist, TWIST_OPTIONS, NUM_TWISTS);
    register_param("AA Width", &params.aa_mult, 0.1f, 1.5f);
    register_param("Hue Shift", &params.hue_shift, -4.0f, 4.0f);
    register_param("Hue Noise Scale", &params.hue_noise_scale, 1.0f / 64.0f,
                   8.0f);
    register_param("Hue Noise Speed", &params.hue_noise_speed, -0.001f, 0.001f);

    build_points();

    baked_palette.bake(persistent_arena, palette);
    volume_spins = persistent_arena.make_n<VolumeSpin>(MAX_POINTS);
    palette_state = persistent_arena.make<PaletteState>();
    palette_state->noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    palette_state->noise.SetSeed(Pullback::HUE_NOISE_SEED);
    palette_state->noise.SetFrequency(1.0f);
    prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE>(
                                 palette_state->hue_rotation_lut),
                             baked_palette);
    noise_palette.bind(&baked_palette, palette_state->hue_rotation_lut.data(),
                       palette_state->hue_noise_lut.data());
    refresh_hue_noise();

    timeline.add(0, Animation::RandomWalk<W>(camera, Y_AXIS, camera_noise));
    for (int i = 0; i < MAX_POINTS; ++i) {
      VolumeSpin &spin = volume_spins[i];
      timeline.add(0, Animation::RandomWalk<W>(
                          spin.orientation, random_vector(), spin.noise,
                          Animation::RandomWalk<W>::Options::Energetic()));
    }

    // spin_phase / palette_phase are effect-owned accumulators wrapped to [0,1)
    // each step, so the trig argument never grows. spin_phase is scaled to
    // radians by *2pi where consumed.
    timeline.add(0, Animation::Driver(spin_phase, &params.pulse_speed,
                                      1.5f / (60.0f * TWO_PI_F), true));
    timeline.add(0, Animation::Driver(palette_phase, &params.pulse_speed,
                                      0.05f / 60.0f, true));
  }

  /**
   * @brief Steps the timeline against a fresh canvas, then ray-marches the
   *        frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    step_choreography();
    refresh_points();
    {
      HS_PROFILE(rm_timeline_step);
      timeline.step(canvas);
    }
    hue_noise_phase = wrap_t(hue_noise_phase + params.hue_noise_speed);
    refresh_hue_noise();
    draw_fn(canvas);
  }

private:
  friend Choreography;
  friend struct ::hs_test::effects_tests::RaymarchWhiteBox;

  using Choreography::begin_choreography;
  using Choreography::step_choreography;
  using Choreography::params;
  using Choreography::register_animated_param;
  using Choreography::register_int_param;
  using Choreography::register_param;
  using Choreography::timeline;

  static_assert(MAX_POINTS + 3 <= Timeline::MAX_EVENTS,
                "Raymarch animations exceed the timeline capacity");

  /** Twist-count labels; the option index IS the twist count draw_fn reads. */
  static constexpr const char *TWIST_OPTIONS[] = {"0", "1", "2", "3", "4",
                                                  "5", "6", "7", "8"};
  static constexpr int NUM_TWISTS = static_cast<int>(std::size(TWIST_OPTIONS));

  // Torus proportions at scale 1: VIS_K is the visible outer ring radius,
  // UNIT_BOUNDS the bounding-sphere radius (bigger, may overlap a neighbour —
  // a few wasted ray steps, no visual overlap).
  static constexpr float MAJOR_K = 0.45f, MINOR_K = 0.14f, TWIST_K = 0.35f;
  static constexpr float VIS_K = MAJOR_K + MINOR_K;

  // Farthest point of the MINOR_K tube about the twisted centerline.
  static constexpr float UNIT_BOUNDS =
      constexpr_sqrt(MAJOR_K * MAJOR_K + TWIST_K * TWIST_K) + MINOR_K;

  /**
   * @brief Builds the vertex directions, per-vertex orientation quaternions,
   *        and nearest-neighbour gaps of the selected base solid.
   * @details Generates the solid from the registry into the scratch arenas
   *          (reclaimed on scope exit) via the shared HS_COLD builder, so draw_fn
   *          can size every torus to its own local gap.
   */
  HS_COLD_MEMBER void build_points() {
    ScratchScope a_guard(scratch_arena_a);
    ScratchScope b_guard(scratch_arena_b);
    const size_t placement_index = static_cast<size_t>(params.base_solid);
    HS_CHECK(placement_index < PLACEMENT_SOLID_COUNT,
             "Raymarch placement solid is out of range");
    const Solids::Entry &entry = Solids::get_entry(
        static_cast<size_t>(PLACEMENT_SOLIDS[placement_index]));
    active_count = Solids::build_vertex_directions(
        scratch_arena_a, scratch_arena_b, entry.name, MAX_POINTS, points.data(),
        raw_quats.data(), nn_angle.data());
    active_base_solid = params.base_solid;
  }

  HS_COLD_MEMBER void refresh_points() {
    if (params.base_solid != active_base_solid)
      build_points();
  }

  HS_FLASH_MEMBER void refresh_hue_noise() {
    if (palette_state->hue_noise_scale == params.hue_noise_scale &&
        palette_state->hue_noise_phase == hue_noise_phase)
      return;
    prepare_hue_noise_lut(
        std::span<int8_t, HueNoiseLutView::SIZE>(palette_state->hue_noise_lut),
        palette_state->noise, params.hue_noise_scale, hue_noise_phase);
    palette_state->hue_noise_scale = params.hue_noise_scale;
    palette_state->hue_noise_phase = hue_noise_phase;
  }

  struct SurfaceFrame {
    Vector normal;
    float cos_u;
    float sin_u;
    float cos_v;
    float sin_v;
  };

  static SurfaceFrame
  surface_frame(const SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> &torus,
                const Vector &loc) {
    const float radial = sqrtf(loc.x * loc.x + loc.z * loc.z);
    const float inverse_radial = 1.0f / radial;
    const auto harmonic = torus.warp.sincos_ntheta_inv(loc, inverse_radial);
    const Vector warped(loc.x, loc.y - torus.warp.amplitude * harmonic.sin_n,
                        loc.z);
    const Vector raw_normal = torus.base.normal_raw(warped, inverse_radial);
    const float inverse_tube = fast_rsqrt(dot(raw_normal, raw_normal));
    return {torus.warp.correct_normal_inv(loc, raw_normal, inverse_radial,
                                          harmonic.cos_n),
            loc.x * inverse_radial, loc.z * inverse_radial,
            (radial - torus.base.R) * inverse_tube, warped.y * inverse_tube};
  }

  /**
   * @brief Ray-marches and shades the twisted torus at every vertex for the
   *        current frame.
   * @param canvas Render target for this frame's fragments.
   */
  void draw_fn(Canvas &canvas) {
    HS_PROFILE(rm_shader_draw);

    int twist_n = static_cast<int>(params.twist + 0.5f);
    int max_steps = params.max_steps;

    // spin_phase rides in [0,1); scale to radians for make_rotation.
    float spin_angle = spin_phase * TWO_PI_F;
    Quaternion spin_q = make_rotation(X_AXIS, spin_angle);

    for (int i = 0; i < active_count; ++i) {
      // Per-vertex auto-size: fit the ring's outer edge to `fill` of this
      // vertex's half nearest-neighbour gap, so open regions get large tori and
      // tight ones stay small; at fill 1 mutual neighbours just touch.
      float outer_r = sinf(0.5f * nn_angle[i] * params.fill);
      float scale = outer_r / VIS_K;
      float major_r = scale * MAJOR_K;
      float minor_r = scale * MINOR_K;
      float twist_amp = scale * TWIST_K;
      float aa_width = minor_r * params.aa_mult;
      // Coverage runs aa_width past the surface, so the cull sphere must too.
      float bounds_radius = scale * UNIT_BOUNDS + aa_width;
      SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
          {major_r, minor_r}, {twist_n, twist_amp, major_r}};
      // Volume::draw's widest band test is min_behind < 2*aa_width, so the
      // cheap bound only has to be accurate above that.
      torus.precision = 2.0f * aa_width;

      Vector center = camera.orient(points[i]);
      Vector ray_dir = -center;

      Quaternion world_q = camera.get() * raw_quats[i] *
                           volume_spins[i].orientation.get() * spin_q;
      Vector tangent = rotate(Vector(1, 0, 0), world_q);
      // The half-vector is fixed for this torus, so it is hoisted out of the
      // per-pixel shader.
      Vector half_w = blinn_phong_half(center, tangent);

      float palette_offset =
          palette_phase + static_cast<float>(i) / active_count;

      auto frag_fn = [&](const Vector &loc, Fragment &frag) {
        SurfaceFrame surface = surface_frame(torus, loc);
        Vector n_world = rotate(surface.normal, world_q);
        float shade = shade_blinn_phong(n_world, center, half_w, params.diffuse,
                                        params.specular, params.fresnel);

        float surface_noise = noise_palette.noise_uv(
            surface.cos_u, surface.sin_u, surface.cos_v, surface.sin_v);
        float palette_t =
            wrap_t(0.5f * (surface_noise + 1.0f) + palette_offset);
        float hue_shift = params.hue_shift * surface_noise;
        Color4 c = params.hue_shift == 0.0f
                       ? baked_palette.get(palette_t)
                       : noise_palette.get(palette_t, hue_shift);
        frag.color = Color4(c.color * shade, 1.0f);
      };

      Scan::TransformedVolume vol(torus, center, world_q);
      Scan::Volume::draw<W, H>(pipeline, canvas, center, bounds_radius, ray_dir,
                               vol, frag_fn, max_steps, aa_width);
    }
  }

  struct VolumeSpin {
    Orientation<> orientation;
    FastNoiseLite noise;
  };

  struct PaletteState {
    std::array<Pixel, HueRotationLutView::SIZE> hue_rotation_lut;
    std::array<int8_t, HueNoiseLutView::SIZE> hue_noise_lut;
    FastNoiseLite noise;
    float hue_noise_scale = -1.0f;
    float hue_noise_phase = -1.0f;
  };

  FastNoiseLite camera_noise;
  Orientation<> camera;
  VolumeSpin *volume_spins = nullptr;
  PaletteState *palette_state = nullptr;
  float spin_phase = 0.0f;    // torus tumble phase, [0,1) -> [0,2pi) radians
  float palette_phase = 0.0f; // baked-palette scroll offset, [0,1) cycles
  float hue_noise_phase = 0.0f;
  int active_count = 0;
  PlacementSolid active_base_solid = PlacementSolid::COUNT;
  std::array<Vector, MAX_POINTS> points;
  std::array<Quaternion, MAX_POINTS> raw_quats;
  /** Per-vertex nearest-neighbour gaps in radians. */
  std::array<float, MAX_POINTS> nn_angle;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{EffectPaletteRecipes::raymarch()};
  BakedPalette baked_palette;
  NoiseHuePalette<BakedPalette> noise_palette;

  // init() bakes one palette LUT into the persistent arena. Effect keeps the
  // default arena split, so the total must fit the device persistent partition.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes() + MAX_POINTS * sizeof(VolumeSpin) +
      alignof(VolumeSpin) + sizeof(PaletteState) + alignof(PaletteState);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Raymarch persistent footprint exceeds the default partition");
};

#include "core/control/registry.h"
REGISTER_EFFECT(Raymarch)
