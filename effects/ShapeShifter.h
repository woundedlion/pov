/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ShapeShifter.h
 * @brief Phase-modulated concentric shapes drawn across the sphere.
 */

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/render/plot/contour_layout.h"
#include "core/render/plot/flat_contour_stack.h"
#include "core/render/plot/timed_contour_stack.h"
#include "core/render/scan/radial_star_surface.h"

namespace hs_test {
namespace shapeshifter_oracle_tests {
struct ShapeShifterWhiteBox;
} // namespace shapeshifter_oracle_tests
} // namespace hs_test

/** @brief Tunable rendering state stored by each ShapeShifter preset. */
struct ShapeShifterParams {
  /** @brief Plot primitives available through the Shape slider. */
  using ShapeType = Plot::ContourShape;

  /** @brief Waveforms available through the Function slider. */
  enum class PhaseFunction : uint8_t { SINE, TRIANGLE, SAWTOOTH, SQUARE };

  /** @brief Per-shape alpha functions available to presets. */
  enum class AlphaFalloff : uint8_t { CONSTANT_HALF, TOWARD_EQUATOR };

  /** @brief Radial distributions available to presets. */
  enum class RadiusSpacing : uint8_t { UNIFORM, SCREEN_BALANCED };

  ShapeType shape{}; ///< Plot primitive drawn for every contour.
  /// Concentric contour count; truncated and clamped to [1, DRAW_LIMIT].
  float count{};
  float sides{};            ///< Polygon or star side count, truncated.
  PhaseFunction function{}; ///< Waveform sampled across the contour radii.
  float amplitude{};        ///< Peak contour rotation, in radians.
  /// Sweep rate; phase advances speed / amplitude turns per frame, capped at
  /// half a turn.
  float speed{};
  bool opposite{}; ///< Toggles the rotation sense past the equator.
  AlphaFalloff alpha_falloff{}; ///< Per-contour alpha function.
  RadiusSpacing spacing{};      ///< Radial distribution of the contours.
};

/**
 * @brief Draws phase-modulated concentric shapes across the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The default P1 renders a continuous radial surface. Other presets
 * draw concentric Plot primitives using the selected Spacing law and sample a
 * selectable waveform at successive radii, producing an animated radial twist. Presets
 * cycle through a Segue::Preset::Fade choreography: the whole effect fades
 * through zero opacity and the parameters snap inside the dark frame.
 */
template <int W, int H>
class ShapeShifter
    : public ChoreographedEffect<ShapeShifter<W, H>, ShapeShifterParams> {
  using Choreography =
      ChoreographedEffect<ShapeShifter<W, H>, ShapeShifterParams>;
  friend Choreography;

public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::ShapeShifter::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  using Params = ShapeShifterParams;           ///< Live parameter struct.
  using ShapeType = Params::ShapeType;         ///< Shape enum.
  using PhaseFunction = Params::PhaseFunction; ///< Waveform enum.
  using AlphaFalloff = Params::AlphaFalloff;   ///< Alpha-function enum.
  using RadiusSpacing = Params::RadiusSpacing; ///< Spacing enum.

  /// Number of `ShapeType` values.
  static constexpr int NUM_SHAPES =
      static_cast<int>(ShapeType::SPHERICAL_STAR) + 1;
  /// Number of `PhaseFunction` values.
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(PhaseFunction::SQUARE) + 1;
  /// Number of `AlphaFalloff` values.
  static constexpr int NUM_ALPHA_FALLOFFS =
      static_cast<int>(AlphaFalloff::TOWARD_EQUATOR) + 1;
  /// Number of `RadiusSpacing` values.
  static constexpr int NUM_RADIUS_SPACINGS =
      static_cast<int>(RadiusSpacing::SCREEN_BALANCED) + 1;
  static constexpr int MAX_SHAPES = 288; ///< Contour storage capacity.
  /** @brief Rendered contours and Count slider are capped at two per row. */
  static constexpr int DRAW_LIMIT = std::min(MAX_SHAPES, 2 * H);
  static constexpr int SURFACE_SIDES = 7; ///< Continuous star lobe count.

  /** @brief Constructs the effect on a WxH canvas. */
  HS_FLASH_MEMBER ShapeShifter()
      : Choreography(
            W, H, pipeline_config<decltype(plot_filters)>({.strobe = true})) {}

  /** @brief Registers the GUI sliders and starts the preset choreography. */
  HS_FLASH_MEMBER void init() override {
    ArenaSplit{SURFACE_SCRATCH_A_BYTES, SURFACE_SCRATCH_B_BYTES}.configure();
    params.count = std::min(params.count, static_cast<float>(DRAW_LIMIT));
    this->template register_described_params<&ShapeShifter::global_fields>(
        *this);
    this->register_described_params();

    layout.init_storage(persistent_arena);
    star_chords.init_storage(persistent_arena, Plot::CONTOUR_MAX_VERTICES);
    timed_stack.init_storage(persistent_arena);
    flat_stack.init_storage(persistent_arena);
    prepare_count(hs::clamp(static_cast<int>(params.count), 1, DRAW_LIMIT));
    timeline.add(0, Animation::RandomWalk<W>(orientation, math::X_AXIS, noise,
                                             {}, hs::rand_int(0, 65536)));
    begin_choreography();
  }

  /** @brief Advances the waveform and draws the full radial shape stack. */
  HS_HOT_FLASH_MEMBER void draw_frame() override {
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(ss_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(ss_timeline_step);
      timeline.step(canvas);
    }
    step_choreography();
    advance_phase();
    plot_filters.prepare(canvas);
    draw_all(canvas);
  }

#if HS_ENABLE_EFFECT_CONTROL_API && defined(HS_PROFILE_SHAPESHIFTER_COUNT)
  void after_profile_preset() {
    static_assert(HS_PROFILE_SHAPESHIFTER_COUNT >= 1 &&
                  HS_PROFILE_SHAPESHIFTER_COUNT <= DRAW_LIMIT);
    params.count = static_cast<float>(HS_PROFILE_SHAPESHIFTER_COUNT);
  }
#endif

  /** @brief Shared registration, validation and interpolation descriptions.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto parameter_fields() {
    return std::tuple{
        Control::Field<Params, ShapeType>{
            .id = "shape",
            .member = &Params::shape,
            .name = "Shape",
            .spec = {.min = 0,
                     .max = static_cast<int64_t>(NUM_SHAPES) - 1,
                     .animated = true,
                     .options = SHAPE_OPTIONS,
                     .export_options = SHAPE_EXPORT_OPTIONS,
                     .option_count = NUM_SHAPES}},
        Control::Field<Params, float>{
            .id = "count",
            .member = &Params::count,
            .name = "Count",
            .spec = {.min = 1.0f,
                     .max = static_cast<float>(DRAW_LIMIT),
                     .animated = true},
            .validation_max = static_cast<float>(MAX_SHAPES)},
        Control::Field<Params, float>{
            .id = "sides",
            .member = &Params::sides,
            .name = "Sides",
            .spec = {.min = SIDES_MIN, .max = SIDES_MAX, .animated = true}},
        Control::Field<Params, PhaseFunction>{
            .id = "function",
            .member = &Params::function,
            .name = "Function",
            .spec = {.min = 0,
                     .max = static_cast<int64_t>(NUM_FUNCTIONS) - 1,
                     .animated = true,
                     .options = FUNCTION_OPTIONS,
                     .export_options = FUNCTION_EXPORT_OPTIONS,
                     .option_count = NUM_FUNCTIONS}},
        Control::Field<Params, float>{.id = "amplitude",
                                      .member = &Params::amplitude,
                                      .name = "Amplitude",
                                      .spec = {.min = AMPLITUDE_MIN,
                                               .max = AMPLITUDE_MAX,
                                               .animated = true}},
        Control::Field<Params, float>{
            .id = "speed",
            .member = &Params::speed,
            .name = "Speed",
            .spec = {.min = SPEED_MIN, .max = SPEED_MAX, .animated = true}},
        Control::Field<Params, bool>{
            .id = "opposite",
            .member = &Params::opposite,
            .name = "Opposite",
            .spec = {.min = 0, .max = 1, .animated = true}},
        Control::Field<Params, AlphaFalloff>{
            .id = "alpha_falloff",
            .member = &Params::alpha_falloff,
            .name = "Alpha Falloff",
            .spec = {.min = 0,
                     .max = static_cast<int64_t>(NUM_ALPHA_FALLOFFS) - 1,
                     .animated = true,
                     .options = ALPHA_FALLOFF_OPTIONS,
                     .export_options = ALPHA_FALLOFF_EXPORT_OPTIONS,
                     .option_count = NUM_ALPHA_FALLOFFS}},
        Control::Field<Params, RadiusSpacing>{
            .id = "spacing",
            .member = &Params::spacing,
            .name = "Spacing",
            .spec = {.min = 0,
                     .max = static_cast<int64_t>(NUM_RADIUS_SPACINGS) - 1,
                     .animated = true,
                     .options = SPACING_OPTIONS,
                     .export_options = SPACING_EXPORT_OPTIONS,
                     .option_count = NUM_RADIUS_SPACINGS}}};
  }

private:
  friend struct ::hs_test::shapeshifter_oracle_tests::ShapeShifterWhiteBox;

#if HS_ENABLE_TEST_HOOKS
  bool test_force_contours = false;
#endif

  HS_HOT_INLINE bool surface_route() const {
#if HS_ENABLE_TEST_HOOKS
    if (test_force_contours)
      return false;
#endif
    return surface_params(params);
  }

  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::step_choreography;
  using Choreography::timeline;

  /** @brief Adopts a snap target; the radial sweep restarts at phase zero. */
  void adopt_params(const Params &target) {
    params = target;
    params.count = std::min(params.count, static_cast<float>(DRAW_LIMIT));
    phase = 0.0f;
    twist_ready = false;
  }

  /** @brief Receives a fading departure's opacity each frame. */
  void set_preset_opacity(float value) { preset_opacity = value; }

  static constexpr float ALPHA_MIN = 0.0f;
  static constexpr float ALPHA_MAX = 1.0f;
  static constexpr float SIDES_MIN = 3.0f;
  static constexpr float SIDES_MAX =
      static_cast<float>(Plot::CONTOUR_MAX_SIDES);
  static constexpr float AMPLITUDE_MIN = 0.1f;
  static constexpr float AMPLITUDE_MAX = 10.0f;
  static constexpr float SPEED_MIN = 0.0f;
  static constexpr float SPEED_MAX = 0.16f;
  static constexpr int PRESET_FRAMES = 480;
  /** Departure fades through black, so the two parameter sets never render on
      the same frame. */
  static constexpr Segue::Preset::Fade DEPARTURE{16};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES =
      PRESET_FRAMES - DEPARTURE.frames;
  static constexpr const char *SHAPE_OPTIONS[] = {
      "Planar Polygon", "Spherical Polygon", "Flower", "Planar Star",
      "Spherical Star"};
  static constexpr const char *SHAPE_EXPORT_OPTIONS[] = {
      "ShapeType::PLANAR_POLYGON", "ShapeType::SPHERICAL_POLYGON",
      "ShapeType::FLOWER", "ShapeType::PLANAR_STAR",
      "ShapeType::SPHERICAL_STAR"};
  static constexpr const char *FUNCTION_OPTIONS[] = {"Sine", "Triangle",
                                                     "Sawtooth", "Square"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "PhaseFunction::SINE", "PhaseFunction::TRIANGLE",
      "PhaseFunction::SAWTOOTH", "PhaseFunction::SQUARE"};
  static constexpr const char *ALPHA_FALLOFF_OPTIONS[] = {"Constant 0.5",
                                                          "Toward Equator"};
  static constexpr const char *ALPHA_FALLOFF_EXPORT_OPTIONS[] = {
      "AlphaFalloff::CONSTANT_HALF", "AlphaFalloff::TOWARD_EQUATOR"};
  static constexpr const char *SPACING_OPTIONS[] = {"Uniform",
                                                    "Screen Balanced"};
  static constexpr const char *SPACING_EXPORT_OPTIONS[] = {
      "RadiusSpacing::UNIFORM", "RadiusSpacing::SCREEN_BALANCED"};
  static_assert(std::size(SHAPE_OPTIONS) == NUM_SHAPES &&
                std::size(SHAPE_EXPORT_OPTIONS) == NUM_SHAPES);
  static_assert(std::size(FUNCTION_OPTIONS) == NUM_FUNCTIONS &&
                std::size(FUNCTION_EXPORT_OPTIONS) == NUM_FUNCTIONS);
  static_assert(std::size(ALPHA_FALLOFF_OPTIONS) == NUM_ALPHA_FALLOFFS &&
                std::size(ALPHA_FALLOFF_EXPORT_OPTIONS) == NUM_ALPHA_FALLOFFS);
  static_assert(std::size(SPACING_OPTIONS) == NUM_RADIUS_SPACINGS &&
                std::size(SPACING_EXPORT_OPTIONS) == NUM_RADIUS_SPACINGS);

  /** @brief Half a turn per frame: the waveform's per-frame Nyquist limit. */
  static constexpr float NYQUIST_PHASE_STEP = 0.5f;

  void advance_phase() {
    // Dividing by amplitude holds the contour's sweep velocity constant across
    // the Amplitude slider.
    const float step =
        fminf(params.speed / params.amplitude, NYQUIST_PHASE_STEP);
    phase = math::wrap_t(phase + step);
    twist_step = twist_ready ? step : 0.0f;
    twist_ready = true;
  }

  static float alpha_falloff_at(AlphaFalloff falloff, float radius_t,
                                int count) {
    if (falloff == AlphaFalloff::CONSTANT_HALF)
      return 0.5f;
    const int steps_to_equator = (count - 1) / 2;
    if (steps_to_equator == 0)
      return 1.0f;
    const float distance_from_pole = hs::clamp(
        Plot::folded_radius_t(radius_t) * static_cast<float>(count) - 0.5f,
        0.0f, static_cast<float>(steps_to_equator));
    const float equator_alpha = 2.0f / static_cast<float>(count);
    return 1.0f - (1.0f - equator_alpha) * distance_from_pole /
                      static_cast<float>(steps_to_equator);
  }

  struct AlphaFalloffModifier {
    AlphaFalloff falloff;
    int count = 1;

    Color4 shade(Color4 color, float radius_t) const {
      color.alpha *= alpha_falloff_at(falloff, radius_t, count);
      return color;
    }
  };

  HS_FLASH_MEMBER void prepare_count(int count) {
    layout.prepare(count, params.spacing == RadiusSpacing::SCREEN_BALANCED);

    MirrorModifier mirror;
    AlphaFalloffModifier constant{AlphaFalloff::CONSTANT_HALF};
    AlphaFalloffModifier toward_equator{AlphaFalloff::TOWARD_EQUATOR, count};
    StaticPalette<ProceduralPalette, Coords<MirrorModifier>,
                  Colors<AlphaFalloffModifier>, false>
        constant_source;
    StaticPalette<ProceduralPalette, Coords<MirrorModifier>,
                  Colors<AlphaFalloffModifier>, false>
        toward_equator_source;
    constant_source.bind(&Palettes::RICH_SUNSET, &mirror, &constant);
    toward_equator_source.bind(&Palettes::RICH_SUNSET, &mirror,
                               &toward_equator);
    if (prepared_count == 0) {
      baked_constant.bake(persistent_arena, constant_source);
      baked_toward_equator.bake(persistent_arena, toward_equator_source);
    } else {
      baked_toward_equator.rebake(toward_equator_source);
    }
    prepared_count = count;
    prepared_spacing = params.spacing;
  }

  void prepare_waveform(PhaseFunction function, int count) {
    layout.prepare_waveform(
        phase, count, function == PhaseFunction::SINE,
        [function](float t) { return evaluate(function, t); });
  }

  const BakedPalette &selected_palette() const {
    return params.alpha_falloff == AlphaFalloff::TOWARD_EQUATOR
               ? baked_toward_equator
               : baked_constant;
  }

  HS_FLASH_MEMBER void
  draw_planar_star_pole_cap(Canvas &canvas, const math::Basis &basis,
                            float geometry_radius_t, float palette_radius_t,
                            int sides, const BakedPalette &palette) {
    const float radius = 2.0f * geometry_radius_t;
    const auto cap = math::get_antipode(basis, radius);
    constexpr float MIN_CAP_RADIUS = 8.0f / W;
    constexpr float CAP_EDGE_OVERLAP = 8.0f / W;
    const float cap_radius = fmaxf(
        MIN_CAP_RADIUS, cap.second * Plot::STAR_INNER_RATIO + CAP_EDGE_OVERLAP);

    Color4 color = palette.get(palette_radius_t);
    color.alpha *=
        fminf(1.0f, alpha * static_cast<float>(sides)) * preset_opacity;
    auto shader = [&](const math::Vector &, Fragment &fragment) {
      fragment.color = color;
    };
    Scan::Circle::draw<W, H>(Scan::ScanRef(plot_filters, canvas), canvas,
                             cap.first, cap_radius, shader);
  }

  HS_FLASH_MEMBER void draw_planar_star_pole_caps(Canvas &canvas,
                                                  const math::Basis &basis,
                                                  int count, int sides,
                                                  const BakedPalette &palette) {
    const float radius_t = 0.5f / static_cast<float>(count);
    draw_planar_star_pole_cap(canvas, basis, layout.spaced_radius_t(0),
                              radius_t, sides, palette);
    // At count 1 both caps resolve to identical arguments; drawing the second
    // would composite the same pixels twice.
    if (count == 1)
      return;
    draw_planar_star_pole_cap(canvas, basis, layout.spaced_radius_t(count - 1),
                              1.0f - radius_t, sides, palette);
  }

  /**
   * @brief Draws the default P1 surface or the selected midpoint-sampled contours.
   * @param canvas Target canvas.
   */
  void draw_all(Canvas &canvas) {
    HS_PROFILE(ss_draw_all);
    const int count = hs::clamp(static_cast<int>(params.count), 1, DRAW_LIMIT);
    if (count != prepared_count || params.spacing != prepared_spacing)
      prepare_count(count);
    const BakedPalette &palette = selected_palette();
    const int sides =
        hs::clamp(static_cast<int>(params.sides), static_cast<int>(SIDES_MIN),
                  static_cast<int>(SIDES_MAX));
    const ShapeType shape = selected_shape();
    const bool planar_star = shape == ShapeType::PLANAR_STAR;
    if (planar_star)
      layout.prepare_star_sides(sides);
    if (surface_route()) {
      Surface::draw(canvas, surface_frame(palette), !Spacetime::flat_surface());
      return;
    }
    const PhaseFunction function = selected_function();
    prepare_waveform(function, count);
    const Plot::ContourFrame frame{palette,
                                   orientation,
                                   shape,
                                   count,
                                   sides,
                                   params.amplitude,
                                   params.opposite,
                                   function == PhaseFunction::SINE,
                                   phase,
                                   twist_step,
                                   alpha * preset_opacity};
    if (!Spacetime::flat_surface()) {
      timed_stack.draw(canvas, plot_filters, star_chords, layout, frame);
      if (planar_star)
        draw_pole_caps_at_columns(canvas, count, sides, palette);
      return;
    }
    flat_stack.draw(canvas, plot_filters, star_chords, layout, frame);
    if (planar_star)
      draw_planar_star_pole_caps(
          canvas, math::make_basis(orientation.get(), math::X_AXIS), count,
          sides, palette);
  }

  /** @brief The continuous star surface's inputs for this frame.
   * @param palette Bound contour palette.
   * @return Phase, twist, opacity, palette and orientation. */
  Scan::SurfaceFrame surface_frame(const BakedPalette &palette) const {
    return {phase,
            twist_step,
            params.amplitude,
            params.opposite,
            alpha * preset_opacity,
            palette,
            orientation};
  }

  /** @brief Draws planar-star pole caps at their column orientations.
   * @param canvas Destination frame.
   * @param count Contour count.
   * @param sides Star tip count.
   * @param palette Bound contour palette. */
  HS_FLASH_MEMBER void draw_pole_caps_at_columns(Canvas &canvas, int count,
                                                 int sides,
                                                 const BakedPalette &palette) {
    ScratchScope column_guard(scratch_arena_a);
    const Scan::ColumnBasis columns =
        Scan::prepare_column_basis<W, H>(scratch_arena_a, orientation);
    Scan::ColumnBasisScope column_scope(canvas, columns);
    draw_planar_star_pole_caps(
        canvas, math::make_basis(orientation.get(), math::X_AXIS), count, sides,
        palette);
  }

  /** @brief Returns the selected Plot primitive. */
  ShapeType selected_shape() const { return params.shape; }

  /** @brief Returns the selected phase waveform. */
  PhaseFunction selected_function() const { return params.function; }

  /**
   * @brief Samples a phase waveform.
   * @param function Waveform to sample.
   * @param t Phase in turns.
   * @return Waveform value in [-1, 1].
   */
  static float evaluate(PhaseFunction function, float t) {
    const float wrapped = t - floorf(t);
    switch (function) {
    case PhaseFunction::SINE:
      return sinf(2.0f * math::PI_F * wrapped);
    case PhaseFunction::TRIANGLE:
      return 1.0f - 4.0f * fabsf(wrapped - 0.5f);
    case PhaseFunction::SAWTOOTH:
      return 2.0f * wrapped - 1.0f;
    case PhaseFunction::SQUARE:
      return wrapped < 0.5f ? 1.0f : -1.0f;
    }
    return 0.0f;
  }

  static constexpr size_t PRESET_COUNT = 9;
  /** @brief Whether clamped parameters select the continuous star surface.
   * @param candidate Rendering parameters.
   * @return True for the surface geometry and sampling policy. */
  static constexpr bool surface_params(const Params &candidate) {
    return candidate.shape == ShapeType::PLANAR_STAR &&
           hs::clamp(static_cast<int>(candidate.count), 1, DRAW_LIMIT) ==
               DRAW_LIMIT &&
           hs::clamp(static_cast<int>(candidate.sides),
                     static_cast<int>(SIDES_MIN),
                     static_cast<int>(SIDES_MAX)) == SURFACE_SIDES &&
           candidate.function == PhaseFunction::SINE &&
           candidate.spacing == RadiusSpacing::SCREEN_BALANCED;
  }
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{.shape = ShapeType::PLANAR_STAR,
        .count = 288.0f,
        .sides = 7.745f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.0f,
        .speed = 0.016f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::TOWARD_EQUATOR,
        .spacing = RadiusSpacing::SCREEN_BALANCED},
       DEPARTURE},
      {{.shape = ShapeType::SPHERICAL_POLYGON,
        .count = 74.644997f,
        .sides = 3.0f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.0f,
        .speed = 0.0318f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::PLANAR_STAR,
        .count = 43.327999f,
        .sides = 6.562f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.0f,
        .speed = 0.0142f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::TOWARD_EQUATOR,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::FLOWER,
        .count = 70.0f,
        .sides = 3.0f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.0f,
        .speed = 0.0186f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::PLANAR_STAR,
        .count = 72.0f,
        .sides = 4.417f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.0f,
        .speed = 0.0077f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::TOWARD_EQUATOR,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::SPHERICAL_POLYGON,
        .count = 128.0f,
        .sides = 5.561f,
        .function = PhaseFunction::SINE,
        .amplitude = 4.0f,
        .speed = 0.0405f,
        .opposite = true,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::SPHERICAL_POLYGON,
        .count = 144.0f,
        .sides = 4.001f,
        .function = PhaseFunction::SINE,
        .amplitude = 2.377f,
        .speed = 0.027086f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::SPHERICAL_POLYGON,
        .count = 144.0f,
        .sides = 3.195f,
        .function = PhaseFunction::SINE,
        .amplitude = 7.0696f,
        .speed = 0.0113f,
        .opposite = false,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{.shape = ShapeType::FLOWER,
        .count = 72.0f,
        .sides = 3.0f,
        .function = PhaseFunction::SINE,
        .amplitude = 1.8721f,
        .speed = 0.00752f,
        .opposite = true,
        .alpha_falloff = AlphaFalloff::CONSTANT_HALF,
        .spacing = RadiusSpacing::UNIFORM},
       DEPARTURE},
  }};
  static_assert(surface_params(PRESETS[0].params));

  FastNoiseLite noise;
  math::Orientation<> orientation;
  Filter::Screen::DirectAntiAliasSink<W, H> plot_filters;
  /** @brief Contours drawn at their columns' times. */
  Plot::TimedContourStack<W, H> timed_stack{orientation};
  BakedPaletteStorage baked_constant;
  BakedPaletteStorage baked_toward_equator;
  Plot::ContourLayout<MAX_SHAPES> layout;
  Plot::PlanarChords<W, H> star_chords;
  /** @brief Contours drawn at the frame's end orientation. */
  Plot::FlatContourStack<W, H> flat_stack;
  int prepared_count = 0;
  RadiusSpacing prepared_spacing = RadiusSpacing::UNIFORM;
  float alpha = 1.0f;
  float preset_opacity = 1.0f;

  /** @brief Descriptions of the effect-owned parameters preset exports omit.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto global_fields() {
    return std::tuple{Control::Field<ShapeShifter, float>{
        .id = "alpha",
        .member = &ShapeShifter::alpha,
        .name = "Alpha",
        .spec = {.min = ALPHA_MIN, .max = ALPHA_MAX, .preset = false}}};
  }
  float phase = 0.0f;
  float twist_step = 0.0f;  ///< Waveform phase advance over this frame, turns.
  bool twist_ready = false; ///< Whether the phase continues the last frame.

  static constexpr size_t SURFACE_SCRATCH_A_BYTES = 32 * 1024;
  static constexpr size_t SURFACE_SCRATCH_B_BYTES = 32 * 1024;
  static_assert(Plot::TimedContourStack<W, H>::scratch_a_bytes() <=
                    SURFACE_SCRATCH_A_BYTES,
                "ShapeShifter column storage exceeds scratch_a");
  static_assert(Plot::FlatContourStack<W, H>::scratch_a_bytes() <=
                    SURFACE_SCRATCH_A_BYTES,
                "ShapeShifter nested contour buffers exceed scratch_a");
  using Surface = Scan::RadialStarSurface<W, H, SURFACE_SIDES>;
  static_assert(Surface::scratch_a_bytes() <= SURFACE_SCRATCH_A_BYTES);
  static_assert(Surface::scratch_b_bytes() <= SURFACE_SCRATCH_B_BYTES);
  static constexpr size_t FOOTPRINT_BYTES =
      Plot::ContourLayout<MAX_SHAPES>::storage_bytes() +
      2 * BakedPalette::required_arena_bytes() +
      Plot::PlanarChords<W, H>::storage_bytes(Plot::CONTOUR_MAX_VERTICES) +
      Plot::FlatContourStack<W, H>::storage_bytes() +
      Plot::TimedContourStack<W, H>::storage_bytes();
  static_assert(FOOTPRINT_BYTES <=
                    ArenaSplit{SURFACE_SCRATCH_A_BYTES, SURFACE_SCRATCH_B_BYTES}
                        .device_persistent(),
                "ShapeShifter persistent footprint exceeds its surface "
                "partition; retune MAX_SHAPES or carve arenas");

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    using Hooks =
        decltype(Control::HookSet<
                 ShapeShifter,
                 Control::HookBinding<
                     Control::ChoreographyHook::PARAMETER_FIELDS,
                     &ShapeShifter::parameter_fields>,
                 Control::HookBinding<Control::ChoreographyHook::ADOPT_PARAMS,
                                      &ShapeShifter::adopt_params>,
                 Control::HookBinding<Control::ChoreographyHook::PRESET_OPACITY,
                                      &ShapeShifter::set_preset_opacity>>{});
#if HS_ENABLE_EFFECT_CONTROL_API && defined(HS_PROFILE_SHAPESHIFTER_COUNT)
    return Hooks::template extend<
        ShapeShifter,
        Control::HookBinding<Control::ChoreographyHook::AFTER_PROFILE_PRESET,
                             &ShapeShifter::after_profile_preset>>();
#else
    return Hooks{};
#endif
  }
};
