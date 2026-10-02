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

#include "core/animation/orientation.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"

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

  ShapeType shape{};
  float count{};
  float sides{};
  PhaseFunction function{};
  float amplitude{};
  float speed{};
  bool opposite{};
  AlphaFalloff alpha_falloff{};
  RadiusSpacing spacing{};

  constexpr ShapeShifterParams() = default;
  constexpr ShapeShifterParams(ShapeType shape, float count, float sides,
                               PhaseFunction function, float amplitude,
                               float speed, float opposite,
                               AlphaFalloff alpha_falloff,
                               RadiusSpacing spacing)
      : shape(shape), count(count), sides(sides), function(function),
        amplitude(amplitude), speed(speed), opposite(opposite >= 0.5f),
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
  static constexpr const char *EFFECT_ID = "ShapeShifter";

  using Params = ShapeShifterParams;
  using ShapeType = Params::ShapeType;
  using PhaseFunction = Params::PhaseFunction;
  using AlphaFalloff = Params::AlphaFalloff;
  using RadiusSpacing = Params::RadiusSpacing;

  static constexpr int NUM_SHAPES =
      static_cast<int>(ShapeType::SPHERICAL_STAR) + 1;
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(PhaseFunction::SQUARE) + 1;
  static constexpr int NUM_ALPHA_FALLOFFS =
      static_cast<int>(AlphaFalloff::TOWARD_EQUATOR) + 1;
  static constexpr int NUM_RADIUS_SPACINGS =
      static_cast<int>(RadiusSpacing::SCREEN_BALANCED) + 1;
  static constexpr int MAX_SHAPES = 288;
  /** @brief Rendered contours and Count slider are capped at two per row. */
  static constexpr int DRAW_LIMIT = std::min(MAX_SHAPES, 2 * H);
  /** @brief Contour count from which star edges switch to screen-step-balanced
   *  sampling: the dense planar-star path, and the balanced policy for
   *  spherical stars. */
  static constexpr float DENSE_CONTOUR_COUNT = 32.0f;

  /** @brief Constructs the Plot-only effect on a WxH canvas. */
  HS_COLD_MEMBER ShapeShifter()
      : Choreography(
            W, H, pipeline_config<decltype(plot_filters)>({.strobe = true})) {}

  /** @brief Registers the GUI sliders and starts the preset choreography. */
  HS_COLD_MEMBER void init() override {
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
    flower_split.init_storage(persistent_arena, FLOWER_MAX_POINTS);
    prepare_count(hs::clamp(static_cast<int>(params.count), 1, DRAW_LIMIT));
    timeline.add(0, Animation::RandomWalk<W>(orientation, math::X_AXIS, noise,
                                             {}, hs::rand_int(0, 65536)));
    begin_choreography();
  }

  /** @brief Advances the waveform and draws the full radial shape stack. */
  void draw_frame() override {
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

#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    HS_CHECK(index < PRESETS.size(),
             "ShapeShifter profile preset index out of range");
    HS_CHECK(this->selectPreset(index),
             "ShapeShifter profile preset selection failed");
#ifdef HS_PROFILE_SHAPESHIFTER_COUNT
    static_assert(HS_PROFILE_SHAPESHIFTER_COUNT >= 1 &&
                  HS_PROFILE_SHAPESHIFTER_COUNT <= DRAW_LIMIT);
    params.count = static_cast<float>(HS_PROFILE_SHAPESHIFTER_COUNT);
#endif
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESETS.size()));
  }
#endif

  /** @brief Shared registration, validation and interpolation descriptions. */
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
  static constexpr int PRESET_FRAMES = 240;
  /** Every preset departs through black over 16 frames, so the two parameter
      sets never render on the same frame. */
  static constexpr Segue::Preset::Fade DEPARTURE{16};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** Dwell + departure = the 240-frame preset cadence. */
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

  /** @brief Half a turn per frame: the waveform's per-frame Nyquist limit. */
  static constexpr float NYQUIST_PHASE_STEP = 0.5f;

  void advance_phase() {
    // Dividing by amplitude holds the contour's sweep velocity constant across
    // the Amplitude slider.
    phase = math::wrap_t(
        phase + std::min(params.speed / params.amplitude, NYQUIST_PHASE_STEP));
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
    const float cap_radius = std::max(
        MIN_CAP_RADIUS, cap.second * Plot::STAR_INNER_RATIO + CAP_EDGE_OVERLAP);

    Color4 color = palette.get(palette_radius_t);
    color.alpha *=
        std::min(1.0f, alpha * static_cast<float>(sides)) * preset_opacity;
    auto shader = [&](const math::Vector &, Fragment &fragment) {
      fragment.color = color;
    };
    Scan::Circle::draw<W, H>(PipelineRef(plot_filters, canvas), canvas,
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

  /**
   * @brief Draws the selected number of midpoint-sampled contours.
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
      const int i =
          continuous_star ? count - ordinal - 1 : folded_draw_indices[ordinal];
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      const bool starts_pair = (count & 1) != 0
                                   ? ordinal == 0 || (ordinal & 1) != 0
                                   : (ordinal & 1) == 0;
      if (continuous_star || starts_pair)
        pair_color = palette.get(radius_t);
      const Color4 color = pair_color;
      const float geometry_radius_t = spaced_radius_t[i];
      const float radius = 2.0f * geometry_radius_t;
      if (planar_star) {
        constexpr float AA_PAD = 2.0f * math::PI_F / W;
        const bool far_side = radius > 1.0f;
        const float cap_radius = far_side ? 2.0f - radius : radius;
        const float half_angle = cap_radius * (math::PI_F / 2.0f) + AA_PAD;
        const float t2 = std::min(half_angle, math::PI_F);
        if (!Plot::cap_may_touch_clip<H>(clip, far_side ? far_cap : near_cap,
                                         t2, sinf(t2)))
          continue;
      }
      const float direction = continuous_star ? star_phase_direction(radius)
                                              : phase_direction(radius);
      const float contour_phase = direction * params.amplitude * waveform[i];
      Color4 shaded_color = color;
      shaded_color.alpha *= global_alpha;
      auto shader = [&](const math::Vector &, Fragment &fragment) {
        fragment.color = shaded_color;
      };
      dispatch_plot(canvas, basis, shape, radius, sides, shader, contour_phase,
                    shaded_color, i, dense_contours);
    }
    if (planar_star)
      draw_planar_star_pole_caps(canvas, basis, count, sides, palette);
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
  HS_NOINLINE_NOCLONE void
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
   * @details The path taken from DENSE_CONTOUR_COUNT contours up, where the
   * per-edge setup of Plot::rasterize's adaptive walk outweighs the edges
   * themselves.
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

  /** @brief Chart pieces a flower's edges share for the band split; each edge
   *  takes FLOWER_PIECE_BUDGET / sides of them. */
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
   * @details Cold (flash): the five-way switch instantiates a sampler lambda
   * per shape, so its body stays out of ITCM even though it runs once per
   * shape (up to DRAW_LIMIT per frame); the hot work is inside Plot::rasterize.
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
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS = {{
      {{ShapeType::PLANAR_STAR, 288.0f, 7.745f, PhaseFunction::SINE, 1.0f,
        0.016f, 0.0f, AlphaFalloff::TOWARD_EQUATOR,
        RadiusSpacing::SCREEN_BALANCED},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 74.644997f, 3.0f, PhaseFunction::SINE,
        1.0f, 0.0318f, 0.0f, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::PLANAR_STAR, 43.327999f, 6.562f, PhaseFunction::SINE, 1.0f,
        0.0142f, 0.0f, AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::FLOWER, 70.0f, 3.0f, PhaseFunction::SINE, 1.0f, 0.0186f,
        0.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::PLANAR_STAR, 72.0f, 4.417f, PhaseFunction::SINE, 1.0f,
        0.0077f, 0.0f, AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 128.0f, 5.561f, PhaseFunction::SINE, 4.0f,
        0.0405f, 1.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 4.001f, PhaseFunction::SINE,
        2.377f, 0.027086f, 0.0f, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 3.195f, PhaseFunction::SINE,
        7.0696f, 0.0113f, 0.0f, AlphaFalloff::CONSTANT_HALF,
        RadiusSpacing::UNIFORM},
       DEPARTURE},
      {{ShapeType::FLOWER, 72.0f, 3.0f, PhaseFunction::SINE, 1.8721f, 0.00752f,
        1.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM},
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

  // init() allocates the six MAX_SHAPES-sized contour tables and the planar
  // chord storage and flower band-split flags; prepare_count() bakes both palettes,
  // from the persistent arena.
  static_assert(SAMPLED_RASTER_CONFIG.single_pass &&
                !SAMPLED_RASTER_CONFIG.derive_planar_arc_registers);
  static constexpr size_t SCRATCH_A_PEAK_BYTES = std::max(
      (2 * static_cast<size_t>(SIDES_MAX) + 2) * sizeof(Fragment) +
          alignof(Fragment) + Plot::PlanarChords<W, H>::scratch_a_bytes(),
      (2 * static_cast<size_t>(SIDES_MAX) + 1 + FLOWER_MAX_POINTS) *
              sizeof(Fragment) +
          2 * alignof(Fragment));
  static_assert(SCRATCH_A_PEAK_BYTES <= DEFAULT_SCRATCH_A_SIZE,
                "ShapeShifter nested contour buffers exceed scratch_a");
  static constexpr size_t FOOTPRINT_BYTES =
      MAX_SHAPES * (4 * sizeof(float) + sizeof(uint16_t) +
                    sizeof(Plot::Star<Plot::PlanarProjection>::RadiusTrig)) +
      2 * BakedPalette::required_arena_bytes() +
      Plot::PlanarChords<W, H>::storage_bytes(STAR_VERTICES) +
      Plot::PlanarBandSplit<W, H>::storage_bytes(FLOWER_MAX_POINTS);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "ShapeShifter persistent footprint exceeds the default "
                "partition; retune MAX_SHAPES or carve arenas");
};
