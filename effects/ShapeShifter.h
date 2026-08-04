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

#include "core/engine/engine.h"

#ifdef HS_TEST_BUILD
namespace hs_test {
namespace shapeshifter_oracle_tests {
struct ShapeShifterWhiteBox;
} // namespace shapeshifter_oracle_tests
} // namespace hs_test
#endif

/**
 * @brief Draws phase-modulated concentric shapes across the sphere.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Evenly spaced Plot primitives share a center and sample a selectable
 * waveform at successive radii, producing an animated radial twist.
 */
template <int W, int H> class ShapeShifter : public Effect {
public:
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
  enum class RadiusSpacing : uint8_t { UNIFORM, ADAPTIVE };

  static constexpr int NUM_SHAPES =
      static_cast<int>(ShapeType::SPHERICAL_STAR) + 1;
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(PhaseFunction::SQUARE) + 1;
  static constexpr int NUM_ALPHA_FALLOFFS =
      static_cast<int>(AlphaFalloff::TOWARD_EQUATOR) + 1;
  static constexpr int NUM_RADIUS_SPACINGS =
      static_cast<int>(RadiusSpacing::ADAPTIVE) + 1;
  /** @brief Count slider ceiling. */
  static constexpr int MAX_SHAPES = 288;

  /** @brief Constructs the Plot-only effect on a WxH canvas. */
  HS_COLD_MEMBER ShapeShifter()
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(plot_filters)::any_crosses_segments,
                .reads_outside_band =
                    decltype(plot_filters)::any_reads_outside_band,
                .margin = decltype(plot_filters)::max_segment_margin}) {}

  /** @brief Loads the initial preset and registers its fields as GUI sliders. */
  void init() override {
    params = presets.get();

    register_animated_param("Alpha", &alpha, ALPHA_MIN, ALPHA_MAX);
    mark_global("Alpha");
    register_animated_param("Shape", &params.shape, SHAPE_OPTIONS,
                            SHAPE_EXPORT_OPTIONS, NUM_SHAPES);
    register_animated_param("Count", &params.count, 1.0f,
                            static_cast<float>(MAX_SHAPES));
    register_animated_param("Sides", &params.sides, SIDES_MIN, SIDES_MAX);
    register_animated_param("Function", &params.function, 0.0f,
                            static_cast<float>(NUM_FUNCTIONS - 1));
    register_animated_param("Amplitude", &params.amplitude, AMPLITUDE_MIN,
                            AMPLITUDE_MAX);
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("Opposite", &params.opposite);
    register_animated_param("Alpha Falloff", &params.alpha_falloff,
                            ALPHA_FALLOFF_OPTIONS, ALPHA_FALLOFF_EXPORT_OPTIONS,
                            NUM_ALPHA_FALLOFFS);
    register_animated_param("Spacing", &params.spacing, SPACING_OPTIONS,
                            SPACING_EXPORT_OPTIONS, NUM_RADIUS_SPACINGS);

    planar_radius_t = persistent_arena.allocate_n<float>(MAX_SHAPES);
    phase_sin = persistent_arena.allocate_n<float>(MAX_SHAPES);
    phase_cos = persistent_arena.allocate_n<float>(MAX_SHAPES);
    waveform = persistent_arena.allocate_n<float>(MAX_SHAPES);
    folded_draw_indices = persistent_arena.allocate_n<uint16_t>(MAX_SHAPES);
    prepare_count(hs::clamp(static_cast<int>(params.count), 1, MAX_SHAPES));
    timeline.add(0, Animation::RandomWalk<W>(orientation, X_AXIS, noise, {},
                                             hs::rand_int(0, 65536)));
    timeline.add_pausable(
        0,
        Animation::PeriodicTimer(
            PRESET_FRAMES, [this](Canvas &) { next_preset(); }, true),
        &anims_paused);
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
    advance_phase();
    plot_filters.prepare(canvas);
    draw_all(canvas);
  }

#if defined(HS_PROFILE_ENABLE) || defined(HS_TEST_BUILD)
  void profile_select_preset(size_t index) {
    HS_CHECK(index < PRESETS.size(),
             "ShapeShifter profile preset index out of range");
    setAnimationsPaused(true);
    params = PRESETS[index].params;
#ifdef HS_PROFILE_SHAPESHIFTER_COUNT
    static_assert(HS_PROFILE_SHAPESHIFTER_COUNT >= 1 &&
                  HS_PROFILE_SHAPESHIFTER_COUNT <= MAX_SHAPES);
    params.count = static_cast<float>(HS_PROFILE_SHAPESHIFTER_COUNT);
#endif
    phase = 0.0f;
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESETS.size()));
  }
#endif

private:
#ifdef HS_TEST_BUILD
  friend struct ::hs_test::shapeshifter_oracle_tests::ShapeShifterWhiteBox;
#endif

  static constexpr float ALPHA_MIN = 0.0f;
  static constexpr float ALPHA_MAX = 1.0f;
  static constexpr float SIDES_MIN = 3.0f;
  static constexpr float SIDES_MAX = 16.0f;
  static constexpr float AMPLITUDE_MIN = 0.1f;
  static constexpr float AMPLITUDE_MAX = 10.0f;
  static constexpr float SPEED_MIN = 0.0f;
  static constexpr float SPEED_MAX = 0.16f;
  static constexpr int PRESET_FRAMES = 240;

  static constexpr const char *SHAPE_OPTIONS[] = {
      "Planar Polygon", "Spherical Polygon", "Flower", "Planar Star",
      "Spherical Star"};
  static constexpr const char *SHAPE_EXPORT_OPTIONS[] = {
      "ShapeType::PLANAR_POLYGON", "ShapeType::SPHERICAL_POLYGON",
      "ShapeType::FLOWER", "ShapeType::PLANAR_STAR",
      "ShapeType::SPHERICAL_STAR"};
  static constexpr const char *ALPHA_FALLOFF_OPTIONS[] = {"Constant 0.5",
                                                          "Toward Equator"};
  static constexpr const char *ALPHA_FALLOFF_EXPORT_OPTIONS[] = {
      "AlphaFalloff::CONSTANT_HALF", "AlphaFalloff::TOWARD_EQUATOR"};
  static constexpr const char *SPACING_OPTIONS[] = {"Uniform", "Adaptive"};
  static constexpr const char *SPACING_EXPORT_OPTIONS[] = {
      "RadiusSpacing::UNIFORM", "RadiusSpacing::ADAPTIVE"};

  /** @brief Tunable rendering state stored by each preset. */
  struct Params {
    ShapeType shape;
    float count;
    float sides;
    float function;
    float amplitude;
    float speed;
    bool opposite;
    AlphaFalloff alpha_falloff;
    RadiusSpacing spacing;

    constexpr Params() = default;
    constexpr Params(ShapeType shape, float count, float sides, float function,
                     float amplitude, float speed, float opposite,
                     AlphaFalloff alpha_falloff, RadiusSpacing spacing)
        : shape(shape), count(count), sides(sides), function(function),
          amplitude(amplitude), speed(speed), opposite(opposite >= 0.5f),
          alpha_falloff(alpha_falloff), spacing(spacing) {}
  };

  void advance_phase() {
    phase = wrap_t(phase + params.speed / params.amplitude);
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

  static constexpr float shape_alpha(AlphaFalloff falloff, int index,
                                     int count) {
    if (falloff == AlphaFalloff::CONSTANT_HALF)
      return 0.5f;
    const int STEPS_TO_EQUATOR = (count - 1) / 2;
    if (STEPS_TO_EQUATOR == 0)
      return 1.0f;
    const int OPPOSITE_INDEX = count - index - 1;
    const int DISTANCE_FROM_POLE =
        index < OPPOSITE_INDEX ? index : OPPOSITE_INDEX;
    const float EQUATOR_ALPHA = 2.0f / static_cast<float>(count);
    return 1.0f - (1.0f - EQUATOR_ALPHA) *
                      static_cast<float>(DISTANCE_FROM_POLE) /
                      static_cast<float>(STEPS_TO_EQUATOR);
  }

  static constexpr float palette_position(float radius_t) {
    return 2.0f * folded_radius_t(radius_t);
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
    int count;

    Color4 shade(Color4 color, float radius_t) const {
      color.alpha *= alpha_falloff_at(falloff, radius_t, count);
      return color;
    }
  };

  /** Maps uniform contour quantiles to sine-weighted latitude density. */
  HS_COLD_MEMBER static float adaptive_planar_radius_t(float radius_t) {
    const bool far_side = radius_t > 0.5f;
    const float u = 2.0f * folded_radius_t(radius_t);
    const float theta = acosf(1.0f - u);
    const float folded = theta / PI_F;
    return far_side ? 1.0f - folded : folded;
  }

  HS_COLD_MEMBER void prepare_count(int count) {
    const bool adaptive_spacing = params.spacing == RadiusSpacing::ADAPTIVE;
    for (int i = 0; i < count; ++i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      phase_sin[i] = sinf(2.0f * PI_F * radius_t);
      phase_cos[i] = cosf(2.0f * PI_F * radius_t);
      planar_radius_t[i] =
          adaptive_spacing ? adaptive_planar_radius_t(radius_t) : radius_t;
      folded_draw_indices[i] =
          static_cast<uint16_t>(folded_draw_index(count, i));
    }

    MirrorModifier mirror;
    AlphaFalloffModifier constant{AlphaFalloff::CONSTANT_HALF, count};
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
    if (baked_palette_count == 0) {
      baked_constant.bake(persistent_arena, constant_source);
      baked_toward_equator.bake(persistent_arena, toward_equator_source);
    } else {
      baked_toward_equator.rebake(toward_equator_source);
    }
    baked_palette_count = count;
    prepared_spacing = params.spacing;
  }

  HS_FLASH_MEMBER void prepare_waveform(PhaseFunction function, int count) {
    if (function == PhaseFunction::SINE) {
      const float phase_angle = 2.0f * PI_F * phase;
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
  draw_planar_star_pole_cap(Canvas &canvas, const Basis &basis,
                            float geometry_radius_t, float palette_radius_t,
                            int sides, const BakedPalette &palette) {
    const float radius = 2.0f * geometry_radius_t;
    const auto cap = get_antipode(basis, radius);
    constexpr float MIN_CAP_RADIUS = 8.0f / W;
    if (cap.second >= MIN_CAP_RADIUS)
      return;

    Color4 color = palette.get(palette_radius_t);
    const float cap_alpha = std::min(1.0f, alpha * static_cast<float>(sides));
    color.alpha = cap_alpha;
    auto shader = [&](const Vector &, Fragment &fragment) {
      fragment.color = color;
    };
    Scan::Circle::draw<W, H>(plot_filters, canvas, cap.first, MIN_CAP_RADIUS,
                             shader);
  }

  HS_FLASH_MEMBER void draw_planar_star_pole_caps(Canvas &canvas,
                                                  const Basis &basis, int count,
                                                  int sides,
                                                  const BakedPalette &palette) {
    const float radius_t = 0.5f / static_cast<float>(count);
    draw_planar_star_pole_cap(canvas, basis, planar_radius_t[0], radius_t,
                              sides, palette);
    draw_planar_star_pole_cap(canvas, basis, planar_radius_t[count - 1],
                              1.0f - radius_t, sides, palette);
  }

  /** @brief Advances to the next preset and applies it atomically. */
  void next_preset() {
    presets.next();
    presets.apply(params);
  }

  /**
   * @brief Draws the selected number of midpoint-sampled contours.
   * @param canvas Target canvas.
   */
  void draw_all(Canvas &canvas) {
    HS_PROFILE(ss_draw_all);
    const int count = hs::clamp(static_cast<int>(params.count), 1, MAX_SHAPES);
    if (count != baked_palette_count || params.spacing != prepared_spacing)
      prepare_count(count);
    const BakedPalette &palette = selected_palette();
    const int sides =
        hs::clamp(static_cast<int>(params.sides), static_cast<int>(SIDES_MIN),
                  static_cast<int>(SIDES_MAX));
    const ShapeType shape = selected_shape();
    const bool planar_star = shape == ShapeType::PLANAR_STAR;
    const PhaseFunction function = selected_function();
    prepare_waveform(function, count);
    const Basis basis = make_basis(orientation.get(), X_AXIS);
    const ClipRegion &clip = canvas.clip();
    const bool full_width_clip = clip.x_start == 0 && clip.x_end == clip.w;
    const float near_cap_beta =
        planar_star ? acosf(hs::clamp(basis.v.y, -1.0f, 1.0f)) : 0.0f;
    const float far_cap_beta =
        planar_star ? acosf(hs::clamp(-basis.v.y, -1.0f, 1.0f)) : 0.0f;
    float cap_half_width = 0.0f;
    float near_cap_distance = 0.0f;
    float far_cap_distance = 0.0f;
    float near_cap_sin_beta = 0.0f;
    float far_cap_sin_beta = 0.0f;
    if (planar_star && !full_width_clip) {
      const float width_px = static_cast<float>(clip.x_end - clip.x_start);
      cap_half_width =
          (width_px * 0.5f + clip.margin + 1.0f) * (2.0f * PI_F) / clip.w;
      const float lambda_center =
          (clip.x_start + width_px * 0.5f) * (2.0f * PI_F) / clip.w;
      auto cap_distance = [&](const Vector &dir) {
        const float lambda = atan2f(dir.z, dir.x);
        return std::fabs(
                   wrap_t((lambda - lambda_center) / (2.0f * PI_F) + 0.5f) -
                   0.5f) *
               (2.0f * PI_F);
      };
      near_cap_distance = cap_distance(basis.v);
      far_cap_distance = cap_distance(-basis.v);
      near_cap_sin_beta = sinf(near_cap_beta);
      far_cap_sin_beta = sinf(far_cap_beta);
    }

    Color4 pair_color;
    const float global_alpha = alpha;
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
      const float geometry_radius_t =
          planar_star ? planar_radius_t[i] : radius_t;
      const float radius = 2.0f * geometry_radius_t;
      if (planar_star) {
        constexpr float AA_PAD = 2.0f * PI_F / W;
        const bool far_side = radius > 1.0f;
        const float cap_radius = far_side ? 2.0f - radius : radius;
        const float half_angle = cap_radius * (PI_F / 2.0f) + AA_PAD;
        const float t2 = std::min(half_angle, PI_F);
        const float beta = far_side ? far_cap_beta : near_cap_beta;
        const float phi_lo = std::max(beta - t2, 0.0f);
        const float phi_hi = std::min(beta + t2, PI_F);
        bool visible = phi_to_y<H>(phi_hi) >= clip.render_y_start() &&
                       phi_to_y<H>(phi_lo) < clip.render_y_end();
        if (visible && !full_width_clip && beta > t2 && PI_F - beta > t2) {
          const float sin_beta =
              far_side ? far_cap_sin_beta : near_cap_sin_beta;
          const float dlam = asinf(hs::clamp(sinf(t2) / sin_beta, 0.0f, 1.0f));
          const float distance =
              far_side ? far_cap_distance : near_cap_distance;
          visible = distance <= dlam + cap_half_width;
        }
        if (!visible)
          continue;
      }
      const float direction = continuous_star ? star_phase_direction(radius)
                                              : phase_direction(radius);
      const float contour_phase = direction * params.amplitude * waveform[i];
      auto shader = [&](const Vector &, Fragment &fragment) {
        fragment.color = color;
        fragment.color.alpha *= global_alpha;
      };
      Color4 shaded_color = color;
      shaded_color.alpha *= global_alpha;
      dispatch_plot(canvas, basis, shape, radius, sides, shader, contour_phase,
                    shaded_color);
    }
    if (planar_star)
      draw_planar_star_pole_caps(canvas, basis, count, sides, palette);
  }

#ifdef HS_TEST_BUILD
  void draw_all_reference(Canvas &canvas) {
    const int count = hs::clamp(static_cast<int>(params.count), 1, MAX_SHAPES);
    if (count != baked_palette_count || params.spacing != prepared_spacing)
      prepare_count(count);
    const BakedPalette &palette = selected_palette();
    const int sides =
        hs::clamp(static_cast<int>(params.sides), static_cast<int>(SIDES_MIN),
                  static_cast<int>(SIDES_MAX));
    const ShapeType shape = selected_shape();
    const PhaseFunction function = selected_function();
    const Basis basis = make_basis(orientation.get(), X_AXIS);
    const float global_alpha = alpha;

    const bool continuous_star = shape == ShapeType::SPHERICAL_STAR;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const int i =
          continuous_star ? count - ordinal - 1 : folded_draw_indices[ordinal];
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      const float geometry_radius_t =
          shape == ShapeType::PLANAR_STAR ? planar_radius_t[i] : radius_t;
      const float radius = 2.0f * geometry_radius_t;
      const float direction = continuous_star ? star_phase_direction(radius)
                                              : phase_direction(radius);
      const float contour_phase =
          direction * params.amplitude * evaluate(function, radius_t + phase);
      const Color4 color = palette.get(radius_t);
      auto shader = [&](const Vector &, Fragment &fragment) {
        fragment.color = color;
        fragment.color.alpha *= global_alpha;
      };
      dispatch_plot_reference(canvas, basis, shape, radius, sides, shader,
                              contour_phase);
    }
    if (shape == ShapeType::PLANAR_STAR)
      draw_planar_star_pole_caps(canvas, basis, count, sides, palette);
  }
#endif

  /** @brief Returns the selected Plot primitive. */
  ShapeType selected_shape() const { return params.shape; }

  /** @brief Returns the nearest valid Function slider selection. */
  PhaseFunction selected_function() const {
    const int selected = hs::clamp(static_cast<int>(params.function + 0.5f), 0,
                                   NUM_FUNCTIONS - 1);
    return static_cast<PhaseFunction>(selected);
  }

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
      return sinf(2.0f * PI_F * wrapped);
    case PhaseFunction::TRIANGLE:
      return 1.0f - 4.0f * fabsf(wrapped - 0.5f);
    case PhaseFunction::SAWTOOTH:
      return 2.0f * wrapped - 1.0f;
    case PhaseFunction::SQUARE:
      return wrapped < 0.5f ? 1.0f : -1.0f;
    }
    return 0.0f;
  }

  /**
   * @brief Plot-rasterizes an open polyline sampled into scratch storage.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param capacity Fragment slots to bind for the sampler.
   * @param planar_basis Azimuthal-equidistant chart for the edges, or nullptr
   * for geodesic edges.
   * @param fragment_shader Per-fragment shader.
   * @param fill Callable that samples the primitive into the bound fragments.
   */
  template <typename F>
  HS_NOINLINE_NOCLONE void
  draw_sampled(Canvas &canvas, size_t capacity, const Basis *planar_basis,
               bool balanced_sampling, const F &fragment_shader, auto &&fill) {
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, capacity);
    fill(points);
    Plot::rasterize<W, H, true, false, false, false,
                    Plot::RasterSamplingPolicy::SELECTABLE>(
        plot_filters, canvas, points, fragment_shader,
        {.planar_basis = planar_basis,
         .omit_end = true,
         .balanced_sampling = balanced_sampling});
  }

  template <typename F>
  HS_FLASH_MEMBER void
  draw_planar_star_edge(Canvas &canvas, const Vector &a, const Vector &b,
                        const Basis &planar_basis, const F &fragment_shader) {
    draw_sampled(canvas, 2, &planar_basis, true, fragment_shader,
                 [&](Fragments &points) {
                   Fragment point;
                   point.pos = a;
                   points.push_back(point);
                   point.pos = b;
                   points.push_back(point);
                 });
  }

  template <typename F>
  HS_FLASH_MEMBER void
  draw_dense_planar_star(Canvas &canvas, const Basis &basis, float radius,
                         int sides, const Color4 &color,
                         const F &fragment_shader, float phase) {
    constexpr int ANCHOR_INTERVALS = 6;
    constexpr float POLE_GUARD_ROWS = 3.0f;
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, static_cast<size_t>(sides * 2 + 2));
    Plot::Star<Plot::PlanarProjection>::sample_positions(points, basis, radius,
                                                         sides, phase);

    Basis planar_basis = basis;
    if (radius > 1.0f)
      planar_basis = Plot::planar_chart_basis(-basis.v);
    const ClipRegion &clip = canvas.clip();
    const ClipRegion::XClip x_clip = clip.x_clip();
    constexpr float TARGET_STEP = Plot::BALANCED_SCREEN_STEP_PX;
    constexpr float ALPHA_GAIN = 1.012f;
    for (int edge = 0; edge < sides * 2; ++edge) {
      const Vector &a = points[edge].pos;
      const Vector &b = points[edge + 1].pos;
      const auto p0 = Plot::azimuthal_project(a, planar_basis);
      const auto p1 = Plot::azimuthal_project(b, planar_basis);
      const float dx = p1.first - p0.first;
      const float dy = p1.second - p0.second;
      const float gap_arc = hypotf(dx, dy) / ANCHOR_INTERVALS;
      constexpr float ROWS_PER_RADIAN =
          static_cast<float>(H + hs::H_OFFSET - 1) / PI_F;

      std::array<PixelCoords, ANCHOR_INTERVALS + 1> anchors;
      float row_lo = static_cast<float>(H);
      float row_hi = 0.0f;
      for (int k = 0; k <= ANCHOR_INTERVALS; ++k) {
        Vector position;
        if (k == 0)
          position = a;
        else if (k == ANCHOR_INTERVALS)
          position = b;
        else {
          const float t = static_cast<float>(k) / ANCHOR_INTERVALS;
          position = Plot::azimuthal_unproject(
              p0.first + dx * t, p0.second + dy * t, planar_basis);
        }
        anchors[k] = vector_to_pixel<W, H>(position);
        if (k > 0) {
          const float delta = anchors[k].x - anchors[k - 1].x;
          if (delta > W * 0.5f)
            anchors[k].x -= W;
          else if (delta < -W * 0.5f)
            anchors[k].x += W;
        }
        row_lo = std::min(row_lo, anchors[k].y);
        row_hi = std::max(row_hi, anchors[k].y);
      }

      const float row_margin = gap_arc * ROWS_PER_RADIAN + 1.0f;
      if (row_lo - row_margin < POLE_GUARD_ROWS ||
          row_hi + row_margin > H - 1.0f - POLE_GUARD_ROWS) {
        draw_planar_star_edge(canvas, a, b, planar_basis, fragment_shader);
        continue;
      }

      for (int k = 0; k < ANCHOR_INTERVALS; ++k) {
        const float segment_dx = anchors[k + 1].x - anchors[k].x;
        const float segment_dy = anchors[k + 1].y - anchors[k].y;
        const float length = hypotf(segment_dx, segment_dy);
        const int samples =
            std::max(1, static_cast<int>(ceilf(length / TARGET_STEP)));
        const float step_ratio =
            std::min(TARGET_STEP, length / samples) / Plot::SCREEN_STEP_PX;
        const float sample_alpha = std::min(
            1.0f,
            ALPHA_GAIN * Plot::balanced_sample_alpha(color.alpha, step_ratio));
        for (int sample = 0; sample < samples; ++sample) {
          const float t = static_cast<float>(sample) / samples;
          float x = anchors[k].x + segment_dx * t;
          const float y = anchors[k].y + segment_dy * t;
          if (x < 0.0f)
            x += W;
          else if (x >= W)
            x -= W;
          const int y0 = static_cast<int>(floorf(y));
          if (!clip.contains_y(y0) && !clip.contains_y(y0 + 1))
            continue;
          const int x0 = static_cast<int>(floorf(x));
          const int x1 = x0 + 1 == W ? 0 : x0 + 1;
          if (x_clip.clipped(x0) && x_clip.clipped(x1))
            continue;
          plot_filters.plot(canvas, x, y, color.color, 0.0f, sample_alpha);
        }
      }
    }
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
   * @details Cold (flash): the five-way switch instantiates a sampler lambda
   * per shape, so its body stays out of ITCM even though it runs once per
   * shape (up to MAX_SHAPES per frame); the hot work is inside Plot::rasterize.
   */
  template <typename F>
  HS_FLASH_MEMBER void
  dispatch_plot(Canvas &canvas, const Basis &basis, ShapeType shape,
                float radius, int sides, const F &fragment_shader,
                float shape_phase, const Color4 &shape_color) {
    HS_PROFILE(ss_plot_dispatch);
    switch (shape) {
    case ShapeType::PLANAR_POLYGON: {
      Basis planar_basis =
          radius > 1.0f ? Plot::planar_chart_basis(-basis.v) : basis;
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
      Basis planar_basis =
          Plot::planar_chart_basis(get_antipode(basis, radius).first.v);
      draw_sampled(canvas, static_cast<size_t>(sides * 2 + 2), &planar_basis,
                   false, fragment_shader, [&](Fragments &points) {
                     Plot::Flower::sample(points, basis, radius, sides,
                                          shape_phase);
                   });
      break;
    }
    case ShapeType::PLANAR_STAR: {
      if (params.count >= 32.0f) {
        draw_dense_planar_star(canvas, basis, radius, sides, shape_color,
                               fragment_shader, shape_phase);
        break;
      }
      Basis planar_basis = basis;
      if (radius > 1.0f)
        planar_basis = Plot::planar_chart_basis(-basis.v);
      draw_sampled(canvas, static_cast<size_t>(sides * 2 + 2), &planar_basis,
                   params.count >= 32.0f, fragment_shader,
                   [&](Fragments &points) {
                     Plot::Star<Plot::PlanarProjection>::sample_positions(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    }
    case ShapeType::SPHERICAL_STAR:
      draw_sampled(
          canvas, static_cast<size_t>(sides * 2 + 2), nullptr,
          params.count >= 32.0f, fragment_shader, [&](Fragments &points) {
            Plot::Star<Plot::GeodesicProjection>::sample_continuous_positions(
                points, basis, radius, sides, shape_phase);
          });
      break;
    }
  }

#ifdef HS_TEST_BUILD
  template <typename F>
  void dispatch_plot_reference(Canvas &canvas, const Basis &basis,
                               ShapeType shape, float radius, int sides,
                               const F &fragment_shader, float shape_phase) {
    switch (shape) {
    case ShapeType::PLANAR_POLYGON:
      Plot::Polygon<Plot::PlanarProjection>::draw<W, H>(
          plot_filters, canvas, basis, radius, sides, fragment_shader,
          shape_phase);
      break;
    case ShapeType::SPHERICAL_POLYGON:
      Plot::Polygon<Plot::GeodesicProjection>::draw<W, H>(
          plot_filters, canvas, basis, radius, sides, fragment_shader,
          shape_phase);
      break;
    case ShapeType::FLOWER:
      Plot::Flower::draw<W, H>(plot_filters, canvas, basis, radius, sides,
                               fragment_shader, {}, shape_phase);
      break;
    case ShapeType::PLANAR_STAR:
      Plot::Star<Plot::PlanarProjection>::draw<W, H>(
          plot_filters, canvas, basis, radius, sides, fragment_shader,
          shape_phase);
      break;
    case ShapeType::SPHERICAL_STAR:
      Plot::Star<Plot::GeodesicProjection>::draw_continuous<W, H>(
          plot_filters, canvas, basis, radius, sides, fragment_shader,
          shape_phase);
      break;
    }
  }
#endif

  static constexpr std::array<PresetEntry<Params>, 9> PRESETS = {{
      {{ShapeType::PLANAR_STAR, 184.0f, 7.745f, 0.0f, 1.0f, 0.016f, 0.0f,
        AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::ADAPTIVE}},
      {{ShapeType::SPHERICAL_POLYGON, 74.644997f, 3.0f, 0.0f, 1.0f, 0.0318f,
        0.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
      {{ShapeType::PLANAR_STAR, 43.327999f, 6.562f, 0.0f, 1.0f, 0.0142f, 0.0f,
        AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM}},
      {{ShapeType::FLOWER, 70.0f, 3.0f, 0.0f, 1.0f, 0.0186f, 0.0f,
        AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
      {{ShapeType::PLANAR_STAR, 72.0f, 4.417f, 0.0f, 1.0f, 0.0077f, 0.0f,
        AlphaFalloff::TOWARD_EQUATOR, RadiusSpacing::UNIFORM}},
      {{ShapeType::SPHERICAL_POLYGON, 128.0f, 5.561f, 0.0f, 4.0f, 0.0405f, 1.0f,
        AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 4.001f, 0.0f, 2.377f, 0.027086f,
        0.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
      {{ShapeType::SPHERICAL_POLYGON, 144.0f, 3.195f, 0.0f, 7.0696f, 0.0113f,
        0.0f, AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
      {{ShapeType::FLOWER, 72.0f, 3.0f, 0.0f, 1.8721f, 0.00752f, 1.0f,
        AlphaFalloff::CONSTANT_HALF, RadiusSpacing::UNIFORM}},
  }};

  static constexpr bool preset_in_ranges(const Params &preset) {
    return static_cast<int>(preset.shape) >= 0 &&
           static_cast<int>(preset.shape) < NUM_SHAPES &&
           preset.count >= 1.0f &&
           preset.count <= static_cast<float>(MAX_SHAPES) &&
           preset.sides >= SIDES_MIN && preset.sides <= SIDES_MAX &&
           preset.function >= 0.0f &&
           preset.function <= static_cast<float>(NUM_FUNCTIONS - 1) &&
           preset.amplitude >= AMPLITUDE_MIN &&
           preset.amplitude <= AMPLITUDE_MAX && preset.speed >= SPEED_MIN &&
           preset.speed <= SPEED_MAX &&
           static_cast<int>(preset.alpha_falloff) >= 0 &&
           static_cast<int>(preset.alpha_falloff) < NUM_ALPHA_FALLOFFS &&
           static_cast<int>(preset.spacing) >= 0 &&
           static_cast<int>(preset.spacing) < NUM_RADIUS_SPACINGS;
  }

  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "ShapeShifter preset is outside a registered slider range");

  FastNoiseLite noise;
  Orientation<> orientation;
  Timeline timeline;
  Filter::Screen::DirectAntiAliasSink<W, H> plot_filters;
  BakedPalette baked_constant;
  BakedPalette baked_toward_equator;
  float *planar_radius_t = nullptr;
  float *phase_sin = nullptr;
  float *phase_cos = nullptr;
  float *waveform = nullptr;
  uint16_t *folded_draw_indices = nullptr;
  int baked_palette_count = 0;
  RadiusSpacing prepared_spacing = RadiusSpacing::UNIFORM;
  Presets<Params, 9> presets{PRESETS};
  Params params{};
  float alpha = 1.0f;
  float phase = 0.0f;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShapeShifter)
