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
#include "core/render/scan/radial_surface_renderer.h"
#include "core/render/scan/radial_surface_candidate_mask.h"

namespace hs_test {
namespace shapeshifter_oracle_tests {
struct ShapeShifterWhiteBox;
} // namespace shapeshifter_oracle_tests
} // namespace hs_test

/** @brief Tunable rendering state stored by each ShapeShifter preset. */
struct ShapeShifterParams {
  /** @brief Plot primitives available through the Shape slider. */
  enum class ShapeType : uint8_t {
    PLANAR_POLYGON,
    SPHERICAL_POLYGON,
    FLOWER,
    PLANAR_STAR,
    SPHERICAL_STAR
  };

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
  /// Sweep rate; the waveform phase advances speed / amplitude turns per frame.
  float speed{};
  bool opposite{}; ///< Toggles the rotation sense past the equator.
  AlphaFalloff alpha_falloff{}; ///< Per-contour alpha function.
  RadiusSpacing spacing{};      ///< Radial distribution of the contours.

  constexpr ShapeShifterParams() = default;
  /**
   * @brief Constructs a preset from every field.
   * @param shape Plot primitive.
   * @param count Contour count.
   * @param sides Side count.
   * @param function Phase waveform.
   * @param amplitude Peak contour rotation, in radians.
   * @param speed Sweep rate.
   * @param opposite Whether to toggle the far-side rotation sense.
   * @param alpha_falloff Per-contour alpha function.
   * @param spacing Radial distribution.
   */
  constexpr ShapeShifterParams(ShapeType shape, float count, float sides,
                               PhaseFunction function, float amplitude,
                               float speed, bool opposite,
                               AlphaFalloff alpha_falloff,
                               RadiusSpacing spacing)
      : shape(shape), count(count), sides(sides), function(function),
        amplitude(amplitude), speed(speed), opposite(opposite),
        alpha_falloff(alpha_falloff), spacing(spacing) {}
};

/**
 * @brief Draws phase-modulated concentric shapes across the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Concentric Plot primitives use the selected Spacing law and sample a selectable
 * waveform at successive radii, producing an animated radial twist. Presets
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
  /** @brief Contour count from which star edges switch to screen-step-balanced
   *  sampling: the dense planar-star path, and the balanced policy for
   *  spherical stars. */
  static constexpr float DENSE_CONTOUR_COUNT = 32.0f;

  /** @brief Constructs the effect on a WxH canvas. */
  HS_COLD_MEMBER ShapeShifter()
      : Choreography(
            W, H, pipeline_config<decltype(plot_filters)>({.strobe = true})) {}

  /** @brief Registers the GUI sliders and starts the preset choreography. */
  HS_COLD_MEMBER void init() override {
    ArenaSplit{SURFACE_SCRATCH_A_BYTES, SURFACE_SCRATCH_B_BYTES}.configure();
    params.count = std::min(params.count, static_cast<float>(DRAW_LIMIT));
    register_param("Alpha", &alpha, ALPHA_MIN, ALPHA_MAX);
    mark_global("Alpha");
    this->register_described_params();

    spaced_radius_t = persistent_arena.allocate_n<float>(MAX_SHAPES);
    phase_sin = persistent_arena.allocate_n<float>(MAX_SHAPES);
    phase_cos = persistent_arena.allocate_n<float>(MAX_SHAPES);
    waveform = persistent_arena.allocate_n<float>(MAX_SHAPES);
    folded_draw_indices = persistent_arena.allocate_n<uint16_t>(MAX_SHAPES);
    planar_star_radius_trig =
        persistent_arena
            .allocate_n<Plot::Star<Plot::PlanarProjection>::RadiusTrig>(
                MAX_SHAPES);
    star_chords.init_storage(persistent_arena, STAR_VERTICES);
    timed_filters.init_storage(persistent_arena);
    flower_split.init_storage(persistent_arena, FLOWER_MAX_POINTS);
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
    draw_all(canvas, COLUMN_TIME_PRESETS[this->getPresetIndex()]);
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

  using Choreography::begin_choreography;
  using Choreography::mark_global;
  using Choreography::params;
  using Choreography::register_param;
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
  static constexpr float SIDES_MAX = 16.0f;
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

  float phase_direction(float radius) const {
    return !params.opposite && radius > 1.0f ? -1.0f : 1.0f;
  }

  float star_phase_direction(float radius) const {
    return params.opposite && radius > 1.0f ? -1.0f : 1.0f;
  }

  static constexpr float folded_radius_t(float radius_t) {
    return radius_t <= 0.5f ? radius_t : 1.0f - radius_t;
  }

  static constexpr int folded_draw_index(int count, int ordinal) {
    if ((count & 1) != 0) {
      if (ordinal == 0)
        return count / 2;
      const int offset = (ordinal + 1) / 2;
      return (ordinal & 1) != 0 ? count / 2 + offset : count / 2 - offset;
    }
    const int offset = ordinal / 2;
    return (ordinal & 1) == 0 ? count / 2 + offset : (count - 1) / 2 - offset;
  }

  static float alpha_falloff_at(AlphaFalloff falloff, float radius_t,
                                int count) {
    if (falloff == AlphaFalloff::CONSTANT_HALF)
      return 0.5f;
    const int steps_to_equator = (count - 1) / 2;
    if (steps_to_equator == 0)
      return 1.0f;
    const float distance_from_pole =
        hs::clamp(folded_radius_t(radius_t) * static_cast<float>(count) - 0.5f,
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

  /** Maps contour quantiles to the screen-space sampling envelope. */
  HS_COLD_MEMBER static float screen_balanced_radius_t(float radius_t) {
    constexpr float DENSITY_FLOOR = 0.5f;
    constexpr float BREAK_ANGLE = math::PI_F / 6.0f;
    constexpr float DENSITY_INTEGRAL = 1.1278247916f;
    constexpr float BREAK_QUANTILE =
        DENSITY_FLOOR * BREAK_ANGLE / DENSITY_INTEGRAL;
    const bool far_side = radius_t > 0.5f;
    const float u = 2.0f * folded_radius_t(radius_t);
    const float theta = u < BREAK_QUANTILE
                            ? u * DENSITY_INTEGRAL / DENSITY_FLOOR
                            : acosf(DENSITY_INTEGRAL * (1.0f - u));
    const float folded = theta / math::PI_F;
    return far_side ? 1.0f - folded : folded;
  }

  HS_COLD_MEMBER void prepare_count(int count) {
    HS_CHECK(count >= 1 && count <= MAX_SHAPES,
             "ShapeShifter: contour count %d outside table capacity", count);
    const bool screen_balanced =
        params.spacing == RadiusSpacing::SCREEN_BALANCED;
    for (int i = 0; i < count; ++i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      phase_sin[i] = sinf(2.0f * math::PI_F * radius_t);
      phase_cos[i] = cosf(2.0f * math::PI_F * radius_t);
      spaced_radius_t[i] =
          screen_balanced ? screen_balanced_radius_t(radius_t) : radius_t;
      planar_star_radius_trig[i] =
          Plot::Star<Plot::PlanarProjection>::radius_trig(2.0f *
                                                          spaced_radius_t[i]);
      folded_draw_indices[i] =
          static_cast<uint16_t>(folded_draw_index(count, i));
    }

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

  HS_FLASH_MEMBER void prepare_waveform(PhaseFunction function, int count) {
    if (function == PhaseFunction::SINE) {
      const float phase_angle = 2.0f * math::PI_F * phase;
      const float phase_sine = sinf(phase_angle);
      const float phase_cosine = cosf(phase_angle);
      for (int i = 0; i < count; ++i)
        waveform[i] = phase_sin[i] * phase_cosine + phase_cos[i] * phase_sine;
      return;
    }

    for (int i = 0; i < count; ++i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      waveform[i] = evaluate(function, radius_t + phase);
    }
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
    draw_planar_star_pole_cap(canvas, basis, spaced_radius_t[0], radius_t,
                              sides, palette);
    // At count 1 both caps resolve to identical arguments; drawing the second
    // would composite the same pixels twice.
    if (count == 1)
      return;
    draw_planar_star_pole_cap(canvas, basis, spaced_radius_t[count - 1],
                              1.0f - radius_t, sides, palette);
  }

  /** @brief A contour's ordered slot and geometric radius. */
  struct ContourSlot {
    int i;          ///< Contour index.
    float radius_t; ///< Palette coordinate.
    float radius;   ///< Geometric radius.
  };

  /** @brief A contour's phase and shaded color. */
  struct ContourShading {
    float direction;     ///< Signed phase direction.
    float contour_phase; ///< Waveform displacement.
    Color4 shaded_color; ///< Stroke color and opacity.
  };

  /** @brief Selects a contour and refreshes its paired palette color.
   * @param ordinal Draw-order position.
   * @param count Number of contours.
   * @param continuous_star Whether contours have independent colors.
   * @param palette Palette sampled by radius.
   * @param pair_color Current paired palette color.
   * @return Contour index and radii. */
  HS_HOT_INLINE ContourSlot contour_slot(int ordinal, int count,
                                         bool continuous_star,
                                         const BakedPalette &palette,
                                         Color4 &pair_color) const {
    const int i =
        continuous_star ? count - ordinal - 1 : folded_draw_indices[ordinal];
    const float radius_t =
        (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
    const bool starts_pair = (count & 1) != 0
                                 ? ordinal == 0 || (ordinal & 1) != 0
                                 : (ordinal & 1) == 0;
    if (continuous_star || starts_pair)
      pair_color = palette.get(radius_t);
    return {i, radius_t, 2.0f * spaced_radius_t[i]};
  }

  /** @brief Tests a planar star contour's cap against the clip.
   * @param clip Drawing clip.
   * @param near_cap,far_cap Cap centers at the two poles.
   * @param radius Contour radius.
   * @param pad Angular motion padding.
   * @return Whether the contour can touch the clip. */
  static HS_HOT_INLINE bool planar_star_cap_visible(
      const ClipRegion &clip, const Plot::CapCenter &near_cap,
      const Plot::CapCenter &far_cap, float radius, float pad) {
    constexpr float AA_PAD = 2.0f * math::PI_F / W;
    const bool far_side = radius > 1.0f;
    const float cap_radius = far_side ? 2.0f - radius : radius;
    const float half_angle = cap_radius * (math::PI_F / 2.0f) + AA_PAD + pad;
    const float t2 = fminf(half_angle, math::PI_F);
    return Plot::cap_may_touch_clip<H>(clip, far_side ? far_cap : near_cap, t2,
                                       sinf(t2));
  }

  /** @brief Evaluates a contour's phase and opacity.
   * @param i Contour index.
   * @param radius Contour radius.
   * @param continuous_star Whether phase follows the spherical-star rule.
   * @param pair_color Selected palette color.
   * @param global_alpha Frame opacity.
   * @return Phase direction, displacement and shaded color. */
  HS_HOT_INLINE ContourShading contour_shading(int i, float radius,
                                               bool continuous_star,
                                               const Color4 &pair_color,
                                               float global_alpha) const {
    const float direction = continuous_star ? star_phase_direction(radius)
                                            : phase_direction(radius);
    const float contour_phase = direction * params.amplitude * waveform[i];
    Color4 shaded_color = pair_color;
    shaded_color.alpha *= global_alpha;
    return {direction, contour_phase, shaded_color};
  }

  /**
   * @brief Draws the selected number of midpoint-sampled contours.
   * @param canvas Target canvas.
   */
  void draw_all(Canvas &canvas, bool column_time = true) {
    HS_PROFILE(ss_draw_all);
    const int count = hs::clamp(static_cast<int>(params.count), 1, DRAW_LIMIT);
    if (count != prepared_count || params.spacing != prepared_spacing)
      prepare_count(count);
    const BakedPalette &palette = selected_palette();
    const int sides =
        hs::clamp(static_cast<int>(params.sides), static_cast<int>(SIDES_MIN),
                  static_cast<int>(SIDES_MAX));
    const ShapeType shape = selected_shape();
    const bool dense_contours =
        static_cast<float>(count) >= DENSE_CONTOUR_COUNT;
    const bool planar_star = shape == ShapeType::PLANAR_STAR;
    if (planar_star && sides != prepared_planar_star_sides) {
      planar_star_step_trig =
          Plot::Star<Plot::PlanarProjection>::step_trig(sides);
      prepared_planar_star_sides = sides;
    }
    const PhaseFunction function = selected_function();
    prepare_waveform(function, count);
    if (this->getPresetIndex() == 0 && shape == ShapeType::PLANAR_STAR &&
        count == std::min(288, DRAW_LIMIT) && sides == 7 &&
        params.function == PhaseFunction::SINE &&
        params.spacing == RadiusSpacing::SCREEN_BALANCED) {
      draw_surface_all(canvas, palette,
                       column_time && !Spacetime::flat_surface());
      return;
    }
    if (column_time && !Spacetime::flat_surface()) {
      draw_timed_all(canvas, count, sides, shape, palette);
      if (planar_star) {
        ScratchScope column_guard(scratch_arena_a);
        const Scan::ColumnBasis column_time =
            Scan::prepare_column_basis<W, H>(scratch_arena_a, orientation);
        Scan::ColumnBasisScope column_scope(canvas, column_time);
        draw_planar_star_pole_caps(
            canvas, math::make_basis(orientation.get(), math::X_AXIS), count,
            sides, palette);
      }
      return;
    }
    const math::Basis basis = math::make_basis(orientation.get(), math::X_AXIS);
    const ClipRegion &clip = canvas.clip();
    Plot::CapCenter near_cap{}, far_cap{};
    if (planar_star) {
      near_cap = Plot::make_cap_center(clip, basis.v);
      far_cap = Plot::make_cap_center(clip, -basis.v);
    }

    if (planar_star && dense_contours)
      star_chords.prepare(clip);
    if (shape == ShapeType::FLOWER)
      flower_band = Plot::ClipBand<W, H>::of(clip);

    Color4 pair_color;
    const float global_alpha = alpha * preset_opacity;
    const bool continuous_star = shape == ShapeType::SPHERICAL_STAR;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const ContourSlot slot =
          contour_slot(ordinal, count, continuous_star, palette, pair_color);
      const int i = slot.i;
      const float radius = slot.radius;
      if (planar_star &&
          !planar_star_cap_visible(clip, near_cap, far_cap, radius, 0.0f))
        continue;
      const ContourShading shading =
          contour_shading(i, radius, continuous_star, pair_color, global_alpha);
      const float contour_phase = shading.contour_phase;
      const Color4 shaded_color = shading.shaded_color;
      auto shader = [&](const math::Vector &, Fragment &fragment) {
        fragment.color = shaded_color;
      };
      dispatch_plot(canvas, basis, shape, radius, sides, shader, contour_phase,
                    shaded_color, i, dense_contours);
    }
    if (planar_star) {
      ScratchScope column_guard(scratch_arena_a);
      const Scan::ColumnBasis column_basis =
          Scan::prepare_column_basis<W, H>(scratch_arena_a, orientation);
      Scan::ColumnBasisScope column_scope(canvas, column_basis);
      draw_planar_star_pole_caps(canvas, basis, count, sides, palette);
    }
  }

  /** @brief Draws the continuous radial star family at each column time. */
  template <bool COUNT = false, bool MASK = true>
  HS_HOT_FLASH_MEMBER void
  draw_surface_all(Canvas &canvas, const BakedPalette &palette, bool timed,
                   Scan::SurfaceAreaStats *stats = nullptr,
                   float row_displacement = 0) {
    const int INTERVALS =
        std::clamp(int(ceilf(fabsf(params.amplitude) * 24)), 24, 240);
    ScratchScope scope_a(scratch_arena_a), scope_b(scratch_arena_b);
    const auto columns = timed ? Scan::prepare_column_basis<W, H>(
                                     scratch_arena_a, orientation, true)
                               : Scan::ColumnBasis{};
    const auto basis = math::make_basis(orientation.get(), math::X_AXIS);
    Surface near, far;
    near.bind(scratch_arena_a, int(params.sides), Plot::STAR_INNER_RATIO,
              INTERVALS);
    far.bind(scratch_arena_a, int(params.sides), Plot::STAR_INNER_RATIO,
             INTERVALS);
    Scan::RadialSurfaceEdgeAA<Surface> integrator;
    Scan::RadialSurfaceCandidateMask<> near_mask, far_mask;
    if constexpr (MASK) {
      near_mask.bind(scratch_arena_b);
      far_mask.bind(scratch_arena_b);
    }
    const float FOOTPRINT =
        .25f * math::PI_F *
            sqrtf(math::RADIANS_PER_COLUMN<W> * math::RADIANS_PER_COLUMN<W> +
                  math::RADIANS_PER_ROW<H> * math::RADIANS_PER_ROW<H>) +
        1e-5f;
    auto *knots =
        scratch_arena_b.allocate_n<typename Surface::Knot>(INTERVALS + 1);
    auto *radii = scratch_arena_b.allocate_n<float>(INTERVALS + 1);
    auto *sines = scratch_arena_b.allocate_n<float>(INTERVALS + 1);
    auto *cosines = scratch_arena_b.allocate_n<float>(INTERVALS + 1);
    auto *colors =
        scratch_arena_b.allocate_n<Scan::SurfaceColor>(INTERVALS + 1);
    auto *first = scratch_arena_b.allocate_n<uint16_t>(257);
    const float opacity = alpha * preset_opacity;
    for (int i = 0; i <= INTERVALS; ++i) {
      const float q = .5f * i / INTERVALS;
      radii[i] = math::PI_F * screen_balanced_radius_t(q);
      sines[i] = sinf(math::TWO_PI_F * q);
      cosines[i] = cosf(math::TWO_PI_F * q);
      const auto color = palette.get(1 - q);
      colors[i] = {color.color.r * color.alpha * opacity,
                   color.color.g * color.alpha * opacity,
                   color.color.b * color.alpha * opacity};
    }
    int index = 0;
    for (int i = 0; i <= 256; ++i) {
      const float radius = i * (.5f * math::PI_F / 256);
      while (index < INTERVALS - 1 && radii[index + 1] < radius)
        ++index;
      first[i] = uint16_t(index);
    }
    float sine = 0, cosine = 1, current_time = 1;
    float near_mask_time = -2, far_mask_time = -2;
    constexpr float MASK_TIME_RANGE = 16.0f / W;
    const float PHASE_MOTION =
        timed ? fabsf(params.amplitude * math::TWO_PI_F * twist_step) *
                    MASK_TIME_RANGE
              : 0;
    const float COMPLEX_MOTION = .5f * math::PI_F * PHASE_MOTION;
    const float RELATIVE_PHASE_MOTION = PHASE_MOTION * math::PI_F / INTERVALS;
    auto prepare = [&](int column) HS_HOT_INLINE {
      const float time =
          timed ? Spacetime::SamplingSurface<W>::time(column) : 1;
      current_time = time;
      const float angle = math::TWO_PI_F * (phase - twist_step * (1 - time));
      sine = sinf(angle);
      cosine = cosf(angle);
    };
    auto ensure = [&](int side) HS_HOT_FLASH_MEMBER {
      HS_PROFILE(ss_surface_knots);
      for (int i = 0; i <= INTERVALS; ++i) {
        const float wave = side ? sine * cosines[i] - cosine * sines[i]
                                : sine * cosines[i] + cosine * sines[i];
        const float turn =
            (side && !params.opposite ? -1 : 1) * params.amplitude * wave;
        float angle = turn - floorf(turn / math::TWO_PI_F) * math::TWO_PI_F;
        const float sign = angle > math::PI_F ? -1.0f : 1.0f;
        if (angle > math::PI_F)
          angle -= math::PI_F;
        float turn_sine, turn_cosine;
        math::poly_sincosf_0_pi(angle, turn_sine, turn_cosine);
        knots[i] = {radii[i] * turn_cosine * sign, radii[i] * turn_sine * sign,
                    turn};
      }
      (side ? far : near)
          .set_knots(
              std::span<const typename Surface::Knot>(knots, INTERVALS + 1));
      if constexpr (MASK) {
        float &mask_time = side ? far_mask_time : near_mask_time;
        if (fabsf(current_time - mask_time) > MASK_TIME_RANGE) {
          HS_PROFILE(ss_surface_mask);
          const auto work = (side ? far_mask : near_mask)
                                .template prepare<COUNT>(
                                    side ? far : near, FOOTPRINT,
                                    COMPLEX_MOTION, RELATIVE_PHASE_MOTION);
          mask_time = current_time;
          if constexpr (COUNT)
            if (stats) {
              ++stats->mask_prepares;
              stats->mask_candidates += work.candidates;
              stats->mask_pieces += work.pieces;
              stats->mask_row_writes += work.row_writes;
            }
          (void)work;
        }
      }
    };
    auto lower_bound = [&](float radius) HS_HOT_INLINE {
      return int(first[std::min(256, int(radius * (512 / math::PI_F)))]);
    };
    auto area = [&](const Surface &surface, auto p, auto dx, auto dy,
                    float theta, int first, auto &shade, auto *stats,
                    float radius, float extent) HS_HOT_INLINE {
      uint32_t candidate_bits = UINT32_MAX;
      int candidate_shift = 0;
      if constexpr (MASK) {
        if constexpr (COUNT)
          if (stats)
            ++stats->mask_queries;
        const auto &mask = &surface == &near ? near_mask : far_mask;
        candidate_bits = mask.interval_bits(radius, theta);
        candidate_shift = mask.interval_shift();
        if (!candidate_bits) {
          typename Surface::Stats work;
          const auto hit = surface.template sample<COUNT, true, false>(
              p, theta, lower_bound(radius), COUNT ? &work : nullptr);
          if constexpr (COUNT)
            if (stats) {
              ++stats->pixels;
              stats->intervals += work.intervals;
              stats->point_intervals += work.intervals;
              stats->max_intervals =
                  std::max(stats->max_intervals, int(work.intervals));
              stats->max_point_intervals =
                  std::max(stats->max_point_intervals, int(work.intervals));
            }
          return hit.found() ? shade(hit) : Scan::SurfaceColor{};
        }
        if constexpr (COUNT)
          if (stats)
            ++stats->mask_hits;
      }
      return integrator.template sample<COUNT>(surface, p, dx, dy, theta, first,
                                               shade, stats, radius, extent,
                                               candidate_bits, candidate_shift);
    };
    auto shade = [&](const typename Surface::Hit &hit) HS_HOT_INLINE {
      const auto a = colors[hit.interval], b = colors[hit.interval + 1];
      const float t = hit.fraction;
      return Scan::SurfaceColor{a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t,
                                a.b + (b.b - a.b) * t};
    };
    Scan::draw_radial_surface<W, H>(canvas, near, far, area, basis, columns,
                                    prepare, ensure, lower_bound, shade, stats,
                                    row_displacement);
  }

  /** @brief Most anchor intervals a flower or planar polygon edge takes. */
  static constexpr int FLOWER_ANCHOR_INTERVALS =
      Plot::PlanarChords<W, H>::MAX_RADIAL_INTERVALS;

  /** @brief Contour vertices one column walk draws. */
  static constexpr size_t COLUMN_CHUNK_VERTICES = 384;
  /** @brief Contours one column walk draws. */
  static constexpr size_t COLUMN_CHUNK_CONTOURS = 64;

  /** @brief A timed contour's twist over the frame and its shading. */
  struct TimedContour {
    math::Vector axis; /**< Axis its twist turns about. */
    float twist;       /**< Signed waveform amplitude, radians. */
    float wave;        /**< Waveform angle at the frame's end, radians. */
    float rate;        /**< Waveform angle change over the frame, radians. */
    float radius;      /**< Distance of its vertices from the axis, at most. */
    Color4 color;      /**< Stroke color and opacity. */
    math::Quaternion start_rotation;  ///< Analytic frame-start turn.
    math::Quaternion middle_rotation; ///< Analytic mid-frame turn.
    Spacetime::PathBound derivatives; ///< Cached trajectory derivative bounds.

    /** @brief Captures the contour's shared turns and derivative bounds. */
    HS_HOT_INLINE void prepare() {
      start_rotation = math::make_rotation(axis, turn_at(0.0f));
      middle_rotation = math::make_rotation(axis, turn_at(0.5f));
      derivatives = path();
    }

    /** @brief Evaluates the shared rigid turn.
     * @param position Frame-end vertex.
     * @param t Frame fraction.
     * @return Vertex at the requested fraction. */
    HS_HOT_INLINE math::Vector at(const math::Vector &position, float t) const {
      if (t == 1.0f)
        return position;
      if (t == 0.0f)
        return math::rotate(position, start_rotation);
      if (t == 0.5f)
        return math::rotate(position, middle_rotation);
      return math::rotate(position, math::make_rotation(axis, turn_at(t)));
    }

    /** @brief Twist from the frame's end back to a frame fraction.
     * @param t Frame fraction.
     * @return The angle to turn a frame-end vertex by, radians. */
    HS_HOT_INLINE float turn_at(float t) const {
      return twist * (sinf(wave - rate * (1.0f - t)) - sinf(wave));
    }

    /** @brief Bounds of a vertex's path over the frame.
     * @return The twist's derivative bounds at the contour's radius. */
    HS_HOT_INLINE Spacetime::PathBound path() const {
      const float amplitude = fabsf(twist) * rate;
      return Spacetime::PathBound::turning(amplitude, amplitude * rate,
                                           amplitude * rate * rate, radius);
    }
  };

  /** @brief A chunk of contour vertices, as Plot::Mesh reads them. */
  struct ContourMesh {
    std::span<const math::Vector> vertices; /**< Frame-end source vertices. */
  };

  /**
   * @brief Draws the contour stack at its columns' times.
   * @param canvas Target canvas.
   * @param count Contours drawn.
   * @param sides Side count.
   * @param shape Primitive.
   * @param palette Contour palette.
   * @details Every contour turns rigidly about its axis by its waveform's
   * twist, whose phase advances linearly over the frame, and the world stage
   * turns the stack by the orientation at each column's time. A contour of
   * geodesic edges is a closed chain the column walk draws, in chunks, once
   * its swept vertices can reach the walk's rows and band.
   */
  HS_HOT_FLASH_MEMBER void draw_timed_all(Canvas &canvas, int count, int sides,
                                          ShapeType shape,
                                          const BakedPalette &palette) {
    const math::Basis basis =
        math::make_basis(math::Quaternion(), math::X_AXIS);
    const bool sine = selected_function() == PhaseFunction::SINE;
    const float rate = sine ? 2.0f * math::PI_F * twist_step : 0.0f;
    const float global_alpha = alpha * preset_opacity;
    const bool continuous_star = shape == ShapeType::SPHERICAL_STAR;
    const auto cull = Plot::SweptEdgeCull<W>::template of<H>(canvas);
    const float spin = timed_filters.sweep_arc();
    const Spacetime::PoseTable &poses = *timed_filters.pose_table(canvas);
    const math::RotationMatrix end_pose = poses.at(poses.times());
    const math::Quaternion start_orientation = orientation.get(0.0f);
    const math::Quaternion middle_orientation = orientation.get(0.5f);
    const math::Quaternion end_orientation = orientation.get();
    const math::Basis world_basis =
        math::make_basis(end_orientation, math::X_AXIS);
    const bool planar = shape == ShapeType::PLANAR_STAR ||
                        shape == ShapeType::PLANAR_POLYGON ||
                        shape == ShapeType::FLOWER;
    const ClipRegion &clip = canvas.clip();
    Plot::CapCenter near_cap{}, far_cap{};
    if (shape == ShapeType::PLANAR_STAR) {
      near_cap = Plot::make_cap_center(clip, world_basis.v);
      far_cap = Plot::make_cap_center(clip, -world_basis.v);
    }
    if (planar)
      star_chords.prepare(canvas, fabsf(params.amplitude) * rate + spin +
                                      2.0f * math::PI_F / W);

    ScratchScope chunk_scope(scratch_arena_a);
    auto *vertices =
        scratch_arena_a.allocate_n<math::Vector>(COLUMN_CHUNK_VERTICES);
    auto *owners = scratch_arena_a.allocate_n<uint16_t>(COLUMN_CHUNK_VERTICES);
    auto *contours =
        scratch_arena_a.allocate_n<TimedContour>(COLUMN_CHUNK_CONTOURS);
    ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_a, COLUMN_CHUNK_VERTICES);
    auto *world = scratch_arena_a.allocate_n<math::Vector>(STAR_VERTICES + 2);
    Fragments ring;
    ring.bind(scratch_arena_a, static_cast<size_t>(STAR_VERTICES + 2));
    size_t count_vertices = 0, count_contours = 0;
    float chunk_motion = 0.0f;

    auto path = [&](size_t i, float t) HS_HOT_FLASH_MEMBER {
      const TimedContour &contour = contours[owners[i]];
      return contour.at(vertices[i], t);
    };
    auto vertex_path = [&](size_t i) HS_HOT_FLASH_MEMBER {
      return contours[owners[i]].derivatives;
    };
    auto shader = [&](const math::Vector &, Fragment &f) HS_HOT_FLASH_MEMBER {
      f.color = contours[owners[edges[static_cast<size_t>(f.v2)].u]].color;
    };
    Plot::GeometryCallback geometry;
    geometry.evaluate = path;
    geometry.vertex_path = vertex_path;
    geometry.source_positions = false;
    geometry.chained = true;
    geometry.rigid_lattice = true;
    geometry.uniform_color = true;
    geometry.capture_error = 0.001f;
    auto flush = [&]() HS_HOT_INLINE {
      if (count_vertices == 0)
        return;
      Plot::Mesh::draw_column_edges<W, H>(
          plot_filters, canvas, ContourMesh{{vertices, count_vertices}}, edges,
          shader, geometry, chunk_motion, &poses, spin);
      count_vertices = 0;
      count_contours = 0;
      chunk_motion = 0.0f;
      edges.clear();
    };

    Color4 pair_color;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const ContourSlot slot =
          contour_slot(ordinal, count, continuous_star, palette, pair_color);
      const int i = slot.i;
      const float radius_t = slot.radius_t;
      const float radius = slot.radius;
      const ContourShading shading =
          contour_shading(i, radius, continuous_star, pair_color, global_alpha);
      const float direction = shading.direction;
      const float contour_phase = shading.contour_phase;
      const Color4 shaded_color = shading.shaded_color;
      const auto frame = Plot::ring_frame(basis, radius);
      // A continuous star turns in the shared basis on both sides of the
      // equator; the other primitives are charted about the far pole past it.
      const math::Basis &turn_basis = continuous_star ? basis : frame.basis;
      TimedContour contour{math::cross(turn_basis.u, turn_basis.w).normalized(),
                           sine ? direction * params.amplitude : 0.0f,
                           2.0f * math::PI_F * (radius_t + phase),
                           rate,
                           continuous_star ? 1.0f : frame.sin_theta,
                           shaded_color,
                           {},
                           {},
                           {}};
      const float motion = fabsf(contour.twist) * rate * contour.radius;
      ring.clear();
      const bool geodesic_flower = shape == ShapeType::FLOWER &&
                                   frame.cos_theta > Plot::COS_PLANAR_ANTIPODE;
      if (planar && !geodesic_flower) {
        if (shape == ShapeType::PLANAR_STAR) {
          if (!planar_star_cap_visible(clip, near_cap, far_cap, radius,
                                       motion + spin))
            continue;
        }
        math::Basis projection_basis;
        const math::Basis *chart = &projection_basis;
        int corners = 2 * sides;
        if (shape == ShapeType::PLANAR_STAR) {
          chart = Plot::PlanarProjection::edge_basis(world_basis, radius,
                                                     projection_basis);
          Plot::Star<Plot::PlanarProjection>::sample_chart_positions(
              ring, star_chords.chart_x(), star_chords.chart_y(), world_basis,
              radius, sides, contour_phase, planar_star_radius_trig[i],
              planar_star_step_trig, *chart);
        } else {
          if (shape == ShapeType::FLOWER) {
            projection_basis = Plot::planar_chart_basis(
                math::get_antipode(world_basis, radius).first.v);
            Plot::Flower::sample(ring, world_basis, radius, sides,
                                 contour_phase);
          } else {
            chart = Plot::PlanarProjection::edge_basis(world_basis, radius,
                                                       projection_basis);
            Plot::Polygon<Plot::PlanarProjection>::sample(
                ring, world_basis, radius, sides, contour_phase);
            corners = sides;
          }
          for (int k = 0; k < corners; ++k) {
            const auto at = Plot::azimuthal_project_exact(ring[k].pos, *chart);
            star_chords.chart_x()[k] = at.first;
            star_chords.chart_y()[k] = at.second;
          }
        }
        const math::Quaternion start_turn =
            start_orientation *
            math::make_rotation(contour.axis, contour.turn_at(0.0f)) *
            end_orientation.conjugate();
        const float lever = Plot::PlanarChords<W, H>::rotation_radius(
            {star_chords.chart_x()[0], star_chords.chart_y()[0]},
            {star_chords.chart_x()[corners > 1 ? 1 : 0],
             star_chords.chart_y()[corners > 1 ? 1 : 0]});
        const float contour_speed = fabsf(contour.twist) * rate * lever + spin;
        flush();
        star_chords.draw_closed(
            plot_filters, canvas, ring, corners, *chart, shaded_color,
            math::RotationMatrix(start_turn),
            shape == ShapeType::PLANAR_STAR
                ? Plot::PlanarChords<W, H>::MAX_ANCHOR_INTERVALS
                : FLOWER_ANCHOR_INTERVALS,
            shape == ShapeType::PLANAR_STAR,
            math::RotationMatrix(
                middle_orientation *
                math::make_rotation(contour.axis, contour.turn_at(0.5f)) *
                end_orientation.conjugate()),
            shape == ShapeType::PLANAR_STAR &&
                float(count) >= DENSE_CONTOUR_COUNT,
            contour_speed);
        continue;
      }
      if (geodesic_flower)
        Plot::Flower::sample(ring, basis, radius, sides, contour_phase);
      else if (continuous_star)
        Plot::Star<Plot::GeodesicProjection>::sample_continuous_positions(
            ring, basis, radius, sides, contour_phase);
      else
        Plot::Polygon<Plot::GeodesicProjection>::sample(ring, basis, radius,
                                                        sides, contour_phase);
      const int points = static_cast<int>(ring.size()) - 1;
      for (int k = 0; k < points; ++k)
        world[k] = end_pose.apply(ring[k].pos);
      if (cull.rejects(world, points, motion + spin))
        continue;
      if (count_vertices + static_cast<size_t>(points) >
              COLUMN_CHUNK_VERTICES ||
          count_contours == COLUMN_CHUNK_CONTOURS)
        flush();
      contour.prepare();
      contours[count_contours] = contour;
      for (int k = 0; k < points; ++k) {
        vertices[count_vertices + k] = ring[k].pos;
        owners[count_vertices + k] = static_cast<uint16_t>(count_contours);
        edges.push_back(
            {static_cast<uint16_t>(count_vertices + k),
             static_cast<uint16_t>(count_vertices + (k + 1) % points)});
      }
      count_vertices += static_cast<size_t>(points);
      ++count_contours;
      chunk_motion = fmaxf(chunk_motion, motion);
    }
    flush();
  }

  /** @brief Most vertices a star contour has. */
  static constexpr int STAR_VERTICES = 2 * static_cast<int>(SIDES_MAX);

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

  static constexpr Plot::RasterConfig SAMPLED_RASTER_CONFIG{
      .single_pass = true,
      .derive_planar_arc_registers = false,
      .interpolate_registers = false,
      .sampling_policy = Plot::RasterSamplingPolicy::SELECTABLE};

  /**
   * @brief Plot-rasterizes an open polyline sampled into scratch storage.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param capacity Fragment slots to bind for the sampler.
   * @param planar_basis Azimuthal-equidistant chart for the edges, or nullptr
   * for geodesic edges.
   * @param balanced_sampling Selects the balanced raster sampling policy.
   * @param fragment_shader Per-fragment shader.
   * @param fill Callable that samples the primitive into the bound fragments.
   */
  template <typename F>
  HS_HOT_FLASH_MEMBER void
  draw_sampled(Canvas &canvas, size_t capacity, const math::Basis *planar_basis,
               bool balanced_sampling, const F &fragment_shader, auto &&fill) {
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, capacity);
    fill(points);
    Plot::rasterize<W, H, SAMPLED_RASTER_CONFIG>(
        plot_filters, canvas, points, fragment_shader,
        {.projection = planar_basis
                           ? Plot::RasterProjection::planar(*planar_basis)
                           : Plot::RasterProjection{},
         .omit_end = true,
         .balanced_sampling = balanced_sampling});
  }

  /**
   * @brief Strokes a planar star with Plot::PlanarChords.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param radius Shape radius in [0, 2].
   * @param sides Star point count.
   * @param fragment_shader Per-fragment shader, used by the pole runs.
   * @param color Contour color.
   * @param phase Star rotation in radians.
   * @param contour_index Contour slot in the baked radius-trig table.
   * @details Used from DENSE_CONTOUR_COUNT contours up.
   */
  template <typename F>
  HS_FLASH_MEMBER void
  draw_dense_planar_star(Canvas &canvas, const math::Basis &basis, float radius,
                         int sides, const F &fragment_shader,
                         const Color4 &color, float phase, int contour_index) {
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, static_cast<size_t>(sides * 2 + 2));
    math::Basis projection_basis;
    const math::Basis &planar_basis =
        *Plot::PlanarProjection::edge_basis(basis, radius, projection_basis);
    Plot::Star<Plot::PlanarProjection>::sample_chart_positions(
        points, star_chords.chart_x(), star_chords.chart_y(), basis, radius,
        sides, phase, planar_star_radius_trig[contour_index],
        planar_star_step_trig, planar_basis);
    star_chords.draw_closed(plot_filters, canvas, points, sides * 2,
                            planar_basis, color, fragment_shader);
  }

  /** @brief Each of the 2 * sides flower edges takes FLOWER_PIECE_BUDGET / sides
   *  pieces, for at most 2 * FLOWER_PIECE_BUDGET pieces per flower. */
  static constexpr int FLOWER_PIECE_BUDGET = 24;
  static constexpr int FLOWER_MAX_POINTS =
      Plot::PlanarBandSplit<W, H>::max_points(2, FLOWER_PIECE_BUDGET);
  static_assert(2 * static_cast<int>(SIDES_MAX) <= 2 * FLOWER_PIECE_BUDGET);

  /**
   * @brief Rasterizes a flower with its petal edges split against the clip
   *        band by Plot::PlanarBandSplit.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param planar_basis The flower's azimuthal-equidistant chart.
   * @param radius Shape radius in [0, 2].
   * @param sides Petal count.
   * @param fragment_shader Per-fragment shader.
   * @param phase Flower rotation in radians.
   */
  template <typename F>
  HS_FLASH_MEMBER void
  draw_banded_flower(Canvas &canvas, const math::Basis &basis,
                     const math::Basis &planar_basis, float radius, int sides,
                     const F &fragment_shader, float phase) {
    const int edges = sides * 2;
    const int pieces = std::max(1, FLOWER_PIECE_BUDGET / sides);
    ScratchScope guard(scratch_arena_a);
    Fragments ring;
    ring.bind(scratch_arena_a, static_cast<size_t>(edges + 1));
    Plot::Flower::sample(ring, basis, radius, sides, phase);
    Fragments path;
    path.bind(scratch_arena_a,
              static_cast<size_t>(
                  Plot::PlanarBandSplit<W, H>::max_points(edges, pieces)));
    const std::span<const uint8_t> flags = flower_split.split(
        path, ring, edges, pieces, planar_basis, flower_band);
    Plot::rasterize<W, H, SAMPLED_RASTER_CONFIG>(
        plot_filters, canvas, path, fragment_shader,
        {.projection = Plot::RasterProjection::planar(planar_basis, flags),
         .omit_end = true});
  }

  /**
   * @brief Samples the selected primitive and draws it.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param shape Primitive to draw.
   * @param radius Shape radius in [0, 2]; above 1 the shape wraps past the
   * pole and is charted about the antipode.
   * @param sides Polygon side, flower petal, or star point count.
   * @param fragment_shader Per-fragment shader.
   * @param shape_phase Primitive rotation in radians.
   * @param shape_color Color applied to this shape.
   * @param contour_index Index of this contour in the stack.
   * @param dense_contours Whether the stack is at or above DENSE_CONTOUR_COUNT
   * contours, selecting the screen-step-balanced star paths.
   * @details Cold (flash); the hot work is inside Plot::rasterize.
   */
  template <typename F>
  HS_FLASH_MEMBER void
  dispatch_plot(Canvas &canvas, const math::Basis &basis, ShapeType shape,
                float radius, int sides, const F &fragment_shader,
                float shape_phase, const Color4 &shape_color, int contour_index,
                bool dense_contours) {
    HS_PROFILE(ss_plot_dispatch);
    switch (shape) {
    case ShapeType::PLANAR_POLYGON: {
      math::Basis projection_basis;
      const math::Basis &planar_basis =
          *Plot::PlanarProjection::edge_basis(basis, radius, projection_basis);
      draw_sampled(canvas, static_cast<size_t>(sides + 2), &planar_basis, false,
                   fragment_shader, [&](Fragments &points) {
                     Plot::Polygon<Plot::PlanarProjection>::sample(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    }
    case ShapeType::SPHERICAL_POLYGON:
      draw_sampled(canvas, static_cast<size_t>(sides + 2), nullptr, false,
                   fragment_shader, [&](Fragments &points) {
                     Plot::Polygon<Plot::GeodesicProjection>::sample(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    case ShapeType::FLOWER: {
      math::Basis planar_basis =
          Plot::planar_chart_basis(math::get_antipode(basis, radius).first.v);
      if (flower_band.x_active) {
        draw_banded_flower(canvas, basis, planar_basis, radius, sides,
                           fragment_shader, shape_phase);
        break;
      }
      draw_sampled(canvas, static_cast<size_t>(sides * 2 + 2), &planar_basis,
                   false, fragment_shader, [&](Fragments &points) {
                     Plot::Flower::sample(points, basis, radius, sides,
                                          shape_phase);
                   });
      break;
    }
    case ShapeType::PLANAR_STAR: {
      if (dense_contours) {
        draw_dense_planar_star(canvas, basis, radius, sides, fragment_shader,
                               shape_color, shape_phase, contour_index);
        break;
      }
      math::Basis projection_basis;
      const math::Basis &planar_basis =
          *Plot::PlanarProjection::edge_basis(basis, radius, projection_basis);
      draw_sampled(canvas, static_cast<size_t>(sides * 2 + 2), &planar_basis,
                   false, fragment_shader, [&](Fragments &points) {
                     Plot::Star<Plot::PlanarProjection>::sample_positions(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    }
    case ShapeType::SPHERICAL_STAR:
      draw_sampled(
          canvas, static_cast<size_t>(sides * 2 + 2), nullptr, dense_contours,
          fragment_shader, [&](Fragments &points) {
            Plot::Star<Plot::GeodesicProjection>::sample_continuous_positions(
                points, basis, radius, sides, shape_phase);
          });
      break;
    }
  }

  static constexpr size_t PRESET_COUNT = 9;
  /** @brief Per-preset column-time sampling policy. */
  static constexpr std::array<bool, PRESET_COUNT> COLUMN_TIME_PRESETS = {
      true, true, true, true, true, true, true, true, true};
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{ShapeType::PLANAR_STAR, 288.0f, 7.745f, PhaseFunction::SINE, 1.0f,
        0.016f, false, AlphaFalloff::TOWARD_EQUATOR,
        RadiusSpacing::SCREEN_BALANCED},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 74.644997f, 3.0f, PhaseFunction::SINE,
        1.0f, 0.0318f, false, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::PLANAR_STAR, 43.327999f, 6.562f, PhaseFunction::SINE, 1.0f,
        0.0142f, false, AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::FLOWER, 70.0f, 3.0f, PhaseFunction::SINE, 1.0f, 0.0186f,
        false, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::PLANAR_STAR, 72.0f, 4.417f, PhaseFunction::SINE, 1.0f,
        0.0077f, false, AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 128.0f, 5.561f, PhaseFunction::SINE, 4.0f,
        0.0405f, true, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 4.001f, PhaseFunction::SINE,
        2.377f, 0.027086f, false, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 3.195f, PhaseFunction::SINE,
        7.0696f, 0.0113f, false, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::FLOWER, 72.0f, 3.0f, PhaseFunction::SINE, 1.8721f, 0.00752f,
        true, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
       DEPARTURE},
  }};

  /** @brief Authored bounds; Count is clamped to the canvas-specific DRAW_LIMIT. */
  static constexpr bool preset_in_ranges(const Params &preset) {
    return Control::valid_fields(preset, parameter_fields());
  }

  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "ShapeShifter preset exceeds authored parameter bounds");

  FastNoiseLite noise;
  math::Orientation<> orientation;
  Filter::Screen::DirectAntiAliasSink<W, H> plot_filters;
  /** @brief Timed contours: orientation at column times, then anti-aliasing. */
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      timed_filters{Filter::World::Orient(orientation),
                    Filter::Screen::AntiAlias<W, H>()};
  BakedPaletteStorage baked_constant;
  BakedPaletteStorage baked_toward_equator;
  float *spaced_radius_t = nullptr;
  float *phase_sin = nullptr;
  float *phase_cos = nullptr;
  float *waveform = nullptr;
  uint16_t *folded_draw_indices = nullptr;
  Plot::Star<Plot::PlanarProjection>::RadiusTrig *planar_star_radius_trig =
      nullptr;
  Plot::Star<Plot::PlanarProjection>::StepTrig planar_star_step_trig{};
  int prepared_planar_star_sides = 0;
  Plot::PlanarChords<W, H> star_chords;
  Plot::PlanarBandSplit<W, H> flower_split;
  Plot::ClipBand<W, H> flower_band;
  int prepared_count = 0;
  RadiusSpacing prepared_spacing = RadiusSpacing::UNIFORM;
  float alpha = 1.0f;
  float preset_opacity = 1.0f;
  float phase = 0.0f;
  float twist_step = 0.0f;  ///< Waveform phase advance over this frame, turns.
  bool twist_ready = false; ///< Whether the phase continues the last frame.

  static_assert(SAMPLED_RASTER_CONFIG.single_pass &&
                !SAMPLED_RASTER_CONFIG.derive_planar_arc_registers);
  /** @brief Bound on chunk storage plus the nested column walk. */
  static constexpr size_t TIMED_SCRATCH_A_BYTES =
      COLUMN_CHUNK_VERTICES *
          (sizeof(math::Vector) + sizeof(uint16_t) + sizeof(Plot::Mesh::Edge)) +
      COLUMN_CHUNK_CONTOURS * sizeof(TimedContour) +
      (STAR_VERTICES + 2) * (sizeof(math::Vector) + sizeof(Fragment)) +
      Plot::Mesh::column_scratch_a_bytes() + 128;
  static_assert(TIMED_SCRATCH_A_BYTES <= DEFAULT_SCRATCH_A_SIZE,
                "ShapeShifter column storage exceeds scratch_a");
  static constexpr size_t SCRATCH_A_PEAK_BYTES = std::max(
      (2 * static_cast<size_t>(SIDES_MAX) + 2) * sizeof(Fragment) +
          alignof(Fragment) + Plot::PlanarChords<W, H>::scratch_a_bytes(),
      (2 * static_cast<size_t>(SIDES_MAX) + 1 + FLOWER_MAX_POINTS) *
              sizeof(Fragment) +
          2 * alignof(Fragment));
  static_assert(SCRATCH_A_PEAK_BYTES <= DEFAULT_SCRATCH_A_SIZE,
                "ShapeShifter nested contour buffers exceed scratch_a");
  static constexpr size_t SURFACE_SCRATCH_A_BYTES = 32 * 1024;
  static constexpr size_t SURFACE_SCRATCH_B_BYTES = 32 * 1024;
  static constexpr int SURFACE_MAX_INTERVALS = 240;
  using Surface = Scan::RadialSurface<SURFACE_MAX_INTERVALS, 16>;
  static_assert(2 * Surface::storage_bytes(7, SURFACE_MAX_INTERVALS) +
                    W * sizeof(math::Vector) + 32 <=
                SURFACE_SCRATCH_A_BYTES);
  static_assert(2 * 64 * 32 * sizeof(uint32_t) +
                    (SURFACE_MAX_INTERVALS + 1) *
                        (sizeof(typename Surface::Knot) + 3 * sizeof(float) +
                         sizeof(Scan::SurfaceColor)) +
                    257 * sizeof(uint16_t) + 32 <=
                SURFACE_SCRATCH_B_BYTES);
  static constexpr size_t FOOTPRINT_BYTES =
      MAX_SHAPES * (4 * sizeof(float) + sizeof(uint16_t) +
                    sizeof(Plot::Star<Plot::PlanarProjection>::RadiusTrig)) +
      2 * BakedPalette::required_arena_bytes() +
      Plot::PlanarChords<W, H>::storage_bytes(STAR_VERTICES) +
      Plot::PlanarBandSplit<W, H>::storage_bytes(FLOWER_MAX_POINTS) +
      Filter::World::Orient::storage_bytes(W);
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
