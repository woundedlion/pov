/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Feedback effect over a fixed icosahedron.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Draws the solid's wireframe under an orientation random-walk while
 * the feedback filter warps and fades the accumulated frame. Style presets
 * cycle without a hard cut: emission stops, the feedback drains to black, then
 * the next preset begins. The shape never changes and never morphs.
 */
template <int W, int H> class MeshFeedback : public Effect {
public:
  using Style = Feedback::Style;

  static constexpr int PRESET_FRAMES = 241;

  // Between presets, emission stops and the feedback fades to black over
  // DRAIN_FRAMES at DRAIN_FADE (well below any preset's fade) before the next
  // preset begins.
  static constexpr int DRAIN_FRAMES = 45;
  static constexpr float DRAIN_FADE = 0.82f;

  static constexpr float FADE_MIN = 0.0f,  FADE_MAX = 0.99f;
  static constexpr float AMP_MIN = 0.0f,   AMP_MAX = 30.0f;
  static constexpr float FREQ_MIN = 0.01f, FREQ_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float SCALE_MIN = 0.1f, SCALE_MAX = 50.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 0.1f;

  /** @brief True iff every preset-driven field of @p s lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Style &s) {
    return s.fade >= FADE_MIN && s.fade <= FADE_MAX &&
           s.amplitude >= AMP_MIN && s.amplitude <= AMP_MAX &&
           s.frequency >= FREQ_MIN && s.frequency <= FREQ_MAX &&
           s.speed >= SPEED_MIN && s.speed <= SPEED_MAX &&
           s.scale >= SCALE_MIN && s.scale <= SCALE_MAX &&
           s.hue_shift >= HUE_SHIFT_MIN && s.hue_shift <= HUE_SHIFT_MAX;
  }
  static_assert(preset_in_ranges(Style::SlowTwist()) &&
                    preset_in_ranges(Style::Churn()) &&
                    preset_in_ranges(Style::Smoke()) &&
                    preset_in_ranges(Style::Frozen()) &&
                    preset_in_ranges(Style::Shatter()) &&
                    preset_in_ranges(Style::Drift()) &&
                    preset_in_ranges(Style::Melting()) &&
                    preset_in_ranges(Style::Swirling()),
                "a MeshFeedback preset drives a style field outside its "
                "registered slider range; widen the range to accommodate the "
                "preset (the range exposes the presets, it does not clamp them)");

  /**
   * @brief Wires up palette, noise, orientation, and the filter pipeline.
   * @details The Feedback filter binds `style` by reference; keep `style`
   * declared before `filters` so it is constructed first (member-declaration
   * init order).
   */
  HS_COLD_MEMBER MeshFeedback()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}), noise_params(),
        orientation(), timeline(),
        palette(Palettes::PEACH_POP),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H>(style)) {}

  /**
   * @brief One-time effect setup.
   * @details Binds shared noise into presets, builds the icosahedron, registers
   * tunable params, and schedules the noise/walk/preset timers.
   */
  void init() override {
    for (auto &e : presets.get_entries()) {
      e.params.noise = &noise_params;
    }

    // Configure the noise type before apply_params(): it calls sync_noise(),
    // which would otherwise propagate the default noise type on the first frame.
    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.sync();

    style = presets.get();
    apply_params();

    {
      PolyMesh poly = generate(persistent_arena, [&](Arena &target, Arena &a,
                                                     Arena &b) {
        return Solids::finalize_solid(Solids::Platonic::icosahedron(a, b),
                                      target);
      });
      MeshOps::compile(poly, mesh_, persistent_arena, scratch_arena_a);
    }

    register_animated_param("Fade", &style.fade, FADE_MIN, FADE_MAX);
    register_animated_param("Distort Amp", &style.amplitude, AMP_MIN, AMP_MAX);
    register_animated_param("Distort Freq", &style.frequency, FREQ_MIN, FREQ_MAX);
    register_animated_param("Distort Speed", &style.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("Noise Scale", &style.scale, SCALE_MIN, SCALE_MAX);
    register_animated_param("Hue Shift", &style.hue_shift, HUE_SHIFT_MIN,
                          HUE_SHIFT_MAX);
    register_param("Feedback", &feedback_enabled);

    filters.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
        persistent_arena);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));
  }

  /**
   * @brief Renders one frame.
   * @details Applies params, advances the preset transition, runs the feedback
   * decay flush, draws the mesh while emitting, then advances the timeline.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(mf_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(mf_apply_params);
      apply_params();
    }

    advance_transition();

    {
      // Feedback-buffer warp/tap + decay flush.
      HS_PROFILE(mf_feedback_flush);
      filters.flush(
          canvas, [](float, float, float) { return Color4(0, 0, 0, 0); },
          1.0f);
    }

    if (emitting_) {
      HS_PROFILE(mf_mesh_draw);
      Color4 shade = palette.get(0.0f);
      Plot::Mesh::draw<W, H>(filters, canvas, mesh_,
                             [&](const Vector &, Fragment &f) { f.color = shade; });
    }

    {
      HS_PROFILE(mf_timeline_step);
      timeline.step(canvas);
    }
  }

private:
  /**
   * @brief Pushes UI-tunable state into the live style/filters each frame.
   * @details Refreshes the noise binding and toggles the feedback filter from
   * `feedback_enabled`.
   */
  void apply_params() {
    style.sync_noise();
    filters.template get<Filter::Pixel::Feedback<W, H>>().set_enabled(
        feedback_enabled);
  }

  /**
   * @brief Drives the emit/drain preset transition.
   * @details Emits for PRESET_FRAMES, then drops the fade to DRAIN_FADE and
   * stops emitting for DRAIN_FRAMES so the accumulated frame decays to black,
   * then advances to the next preset (restoring its fade) and resumes emitting.
   * Frozen while animations are paused.
   */
  void advance_transition() {
    if (animations_paused()) return;
    if (++transition_frames_ < (emitting_ ? PRESET_FRAMES : DRAIN_FRAMES))
      return;
    transition_frames_ = 0;
    emitting_ = !emitting_;
    if (emitting_) {
      presets.next();
      presets.apply(style);
    } else {
      style.fade = DRAIN_FADE;
    }
  }

  Style style;

  Presets<Style, 8> presets{std::array<PresetEntry<Style>, 8>{{{Style::SlowTwist()},
                                                               {Style::Melting()},
                                                               {Style::Swirling()},
                                                               {Style::Churn()},
                                                               {Style::Smoke()},
                                                               {Style::Frozen()},
                                                               {Style::Shatter()},
                                                               {Style::Drift()}}}};
  bool feedback_enabled = true;
  bool emitting_ = true;
  int transition_frames_ = 0;
  NoiseParams noise_params;

  Orientation<> orientation;
  Timeline timeline;
  ProceduralPalette palette;

  // The single, fixed solid; built once in init() and never recompiled.
  MeshState mesh_;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H>>
      filters;

  // init() allocates the Feedback warp-field cache from the persistent arena.
  static_assert(Filter::Pixel::Feedback<W, H>::STORAGE_BYTES <=
                    DEVICE_PERSISTENT_BUDGET,
                "MeshFeedback warp cache exceeds the default persistent "
                "partition; retune the feedback downsample or carve arenas");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
