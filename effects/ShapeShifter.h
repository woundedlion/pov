/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

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
  enum class ShapeType { PLANAR_POLYGON, SPHERICAL_POLYGON, FLOWER, STAR };

  /** @brief Waveforms available through the Function slider. */
  enum class PhaseFunction { SINE, TRIANGLE, SAWTOOTH, SQUARE };

  static constexpr int NUM_SHAPES = static_cast<int>(ShapeType::STAR) + 1;
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(PhaseFunction::SQUARE) + 1;
  static constexpr int MAX_SHAPES = H > 75 ? H : 75;

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

    register_animated_param("Alpha", &params.alpha, ALPHA_MIN, ALPHA_MAX);
    register_animated_param("Shape", &params.shape, 0.0f,
                            static_cast<float>(NUM_SHAPES - 1));
    register_animated_param("Count", &params.count, 1.0f,
                            static_cast<float>(MAX_SHAPES));
    register_animated_param("Sides", &params.sides, SIDES_MIN, SIDES_MAX);
    register_animated_param("Function", &params.function, 0.0f,
                            static_cast<float>(NUM_FUNCTIONS - 1));
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);

    baked_sunset.bake(persistent_arena, Palettes::RICH_SUNSET);
    timeline.add(0, Animation::RandomWalk<W>(orientation, X_AXIS, noise, {},
                                             hs::rand_int(0, 65536)));
    next_preset();
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
    phase = wrap_t(phase + params.speed);
    draw_all(canvas);
  }

private:
  static constexpr float ALPHA_MIN = 0.0f;
  static constexpr float ALPHA_MAX = 1.0f;
  static constexpr float SIDES_MIN = 3.0f;
  static constexpr float SIDES_MAX = 16.0f;
  static constexpr float SPEED_MIN = 0.0f;
  static constexpr float SPEED_MAX = 0.1f;
  static constexpr float PHASE_AMPLITUDE = 1.0f;
  static constexpr int PRESET_FRAMES = 240;

  /** @brief Tunable rendering state stored by each preset. */
  struct Params {
    float alpha;
    float shape;
    float count;
    float sides;
    float function;
    float speed;

    /** @brief Interpolates every preset field from a to b. */
    void lerp(const Params &a, const Params &b, float t) {
      alpha = hs::lerp(a.alpha, b.alpha, t);
      shape = hs::lerp(a.shape, b.shape, t);
      count = hs::lerp(a.count, b.count, t);
      sides = hs::lerp(a.sides, b.sides, t);
      function = hs::lerp(a.function, b.function, t);
      speed = hs::lerp(a.speed, b.speed, t);
    }
  };

  /** @brief Advances to the next preset and schedules a continuous blend. */
  void next_preset() {
    presets.next();
    timeline.add(0,
                 Animation::Lerp(params, presets.prev_get(), presets.get(),
                                 PRESET_FRAMES, ease_in_out_sin, &anims_paused)
                     .then([this]() { next_preset(); }));
  }

  /**
   * @brief Draws the selected shape at equally spaced midpoint radii.
   * @param canvas Target canvas.
   */
  void draw_all(Canvas &canvas) {
    HS_PROFILE(ss_draw_all);
    const int count = hs::clamp(static_cast<int>(params.count), 1, MAX_SHAPES);
    const int sides =
        hs::clamp(static_cast<int>(params.sides), static_cast<int>(SIDES_MIN),
                  static_cast<int>(SIDES_MAX));
    const ShapeType shape = selected_shape();
    const PhaseFunction function = selected_function();
    const Basis basis = make_basis(orientation.get(), X_AXIS);

    for (int i = count - 1; i >= 0; --i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      const float radius = 2.0f * radius_t;
      const float shape_phase =
          PHASE_AMPLITUDE * evaluate(function, radius_t + phase);
      const Color4 color = baked_sunset.get(radius_t);

      auto shader = [&](const Vector &, Fragment &fragment) {
        fragment.color = color;
        fragment.color.alpha *= params.alpha;
      };
      dispatch_plot(canvas, basis, shape, radius, sides, shader, shape_phase);
    }
  }

  /** @brief Returns the nearest valid Shape slider selection. */
  ShapeType selected_shape() const {
    const int selected =
        hs::clamp(static_cast<int>(params.shape + 0.5f), 0, NUM_SHAPES - 1);
    return static_cast<ShapeType>(selected);
  }

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
   * @brief Plot-rasterizes one selected primitive.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param shape Primitive to draw.
   * @param radius Shape radius in [0, 2].
   * @param sides Polygon side, flower petal, or star point count.
   * @param fragment_shader Per-fragment shader.
   * @param shape_phase Primitive rotation in radians.
   */
  template <typename F>
  __attribute__((noinline)) void
  dispatch_plot(Canvas &canvas, const Basis &basis, ShapeType shape,
                float radius, int sides, const F &fragment_shader,
                float shape_phase) {
    HS_PROFILE(ss_plot_dispatch);
    switch (shape) {
    case ShapeType::PLANAR_POLYGON:
      Plot::PlanarPolygon::draw<W, H>(plot_filters, canvas, basis, radius,
                                      sides, fragment_shader, shape_phase);
      break;
    case ShapeType::SPHERICAL_POLYGON:
      Plot::SphericalPolygon::draw<W, H>(plot_filters, canvas, basis, radius,
                                         sides, fragment_shader, shape_phase);
      break;
    case ShapeType::FLOWER:
      Plot::Flower::draw<W, H>(plot_filters, canvas, basis, radius, sides,
                               fragment_shader, {}, shape_phase);
      break;
    case ShapeType::STAR:
      Plot::Star::draw<W, H>(plot_filters, canvas, basis, radius, sides,
                             fragment_shader, shape_phase);
      break;
    }
  }

  static constexpr std::array<PresetEntry<Params>, 3> PRESETS = {{
      {{0.5f, 1.017f, 74.644997f, 3.0f, 0.0f, 0.0318f}},
      {{0.5f, 2.793f, 43.327999f, 6.562f, 0.0f, 0.0142f}},
      {{0.5f, 1.872f, 70.0f, 3.0f, 0.0f, 0.0186f}},
  }};

  static constexpr bool preset_in_ranges(const Params &preset) {
    return preset.alpha >= ALPHA_MIN && preset.alpha <= ALPHA_MAX &&
           preset.shape >= 0.0f &&
           preset.shape <= static_cast<float>(NUM_SHAPES - 1) &&
           preset.count >= 1.0f &&
           preset.count <= static_cast<float>(MAX_SHAPES) &&
           preset.sides >= SIDES_MIN && preset.sides <= SIDES_MAX &&
           preset.function >= 0.0f &&
           preset.function <= static_cast<float>(NUM_FUNCTIONS - 1) &&
           preset.speed >= SPEED_MIN && preset.speed <= SPEED_MAX;
  }

  static_assert(preset_in_ranges(PRESETS[0].params) &&
                    preset_in_ranges(PRESETS[1].params) &&
                    preset_in_ranges(PRESETS[2].params),
                "ShapeShifter preset is outside a registered slider range");

  FastNoiseLite noise;
  Orientation<> orientation;
  Timeline timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> plot_filters;
  BakedPalette baked_sunset;
  Presets<Params, 3> presets{PRESETS};
  Params params{};
  float phase = 0.0f;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShapeShifter)
