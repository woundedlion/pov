/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <numeric>

#include "core/engine/engine.h"

/**
 * @brief Feedback effect over a morphing carousel of Platonic solids.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Draws the current mesh, cycles style presets on a hard-cut timer,
 * and cross-morphs to the next solid on a separate timer. The two cycle periods
 * are kept coprime so preset cuts and shape morphs drift out of phase rather
 * than locking.
 */
template <int W, int H> class MeshFeedback : public Effect {
public:
  using Style = Feedback::Style;

  static constexpr int MORPH_FRAMES = 160;
  static constexpr int PRESET_FRAMES = 241;
  static constexpr int SHAPE_FRAMES = 240;
  static constexpr int NO_MORPH_FRAMES = SHAPE_FRAMES - MORPH_FRAMES;

  static_assert(std::gcd(NO_MORPH_FRAMES + MORPH_FRAMES, PRESET_FRAMES) == 1,
                "the shape cycle (NO_MORPH_FRAMES + MORPH_FRAMES) and "
                "PRESET_FRAMES must stay coprime so the shape and preset cycles "
                "drift out of phase instead of locking");

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
   * @details The Feedback pixel filter reads `style` by reference. Initializers
   * run in member-DECLARATION order, and `style` is declared before `filters`,
   * so style is fully constructed before filters binds its reference. Keep style
   * declared ahead of filters.
   */
  FLASHMEM MeshFeedback()
      : Effect(W, H, decltype(filters)::any_crosses_segments), noise_params(),
        orientation(), timeline(),
        palette(Palettes::peachPop),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H>(style)) {}

  /**
   * @brief One-time effect setup.
   * @details Binds shared noise into presets, builds the first solid, registers
   * tunable params, and schedules the noise/walk/preset/shape timers.
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

    solid_idx = 0;
    {
      auto solids = Solids::Collections::get_platonic_solids();
      PolyMesh poly = generate(persistent_arena, [&](Arena &target, Arena &a,
                                                     Arena &b) {
        return Solids::finalize_solid(solids[solid_idx].generate(a, b), target);
      });
      carousel.slot(carousel.front_index()) = MeshState();
      MeshOps::compile(poly, carousel.slot(carousel.front_index()),
                       persistent_arena);
    }

    registerAnimatedParam("Fade", &style.fade, FADE_MIN, FADE_MAX);
    registerAnimatedParam("Distort Amp", &style.amplitude, AMP_MIN, AMP_MAX);
    registerAnimatedParam("Distort Freq", &style.frequency, FREQ_MIN, FREQ_MAX);
    registerAnimatedParam("Distort Speed", &style.speed, SPEED_MIN, SPEED_MAX);
    registerAnimatedParam("Noise Scale", &style.scale, SCALE_MIN, SCALE_MAX);
    registerAnimatedParam("Hue Shift", &style.hue_shift, HUE_SHIFT_MIN,
                          HUE_SHIFT_MAX);
    registerParam("Feedback", &feedback_enabled);

    filters.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
        persistent_arena);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));

    timeline.add(0, Animation::PeriodicTimer(
                        PRESET_FRAMES,
                        [this](Canvas &) {
                          if (animationsPaused())
                            return;
                          presets.next();
                          presets.apply(style);
                        },
                        true));

    arm_hold_timer();
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Renders one frame.
   * @details Applies params, runs the feedback decay flush, draws the current
   * mesh (skipped while morphing — the morph animation draws instead), then
   * advances the timeline.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    apply_params();

    filters.flush(
        canvas, [](float, float, float) { return Color4(0, 0, 0, 0); },
        1.0f);

    if (!morphing) {
      Color4 shade = palette.get(0.0f);
      Plot::Mesh::draw<W, H>(filters, canvas, carousel.current(),
                             [&](const Vector &, Fragment &f) { f.color = shade; });
    }

    timeline.step(canvas);
  }

private:
  /**
   * @brief Arms the one-shot hold timer that triggers the next shape morph.
   * @details The shape carousel is a chain of non-repeating hold timers: each
   * fires once after NO_MORPH_FRAMES, starts a morph, and the completed morph
   * re-arms the next hold. Both the initial arm (init()) and the post-morph
   * re-arm route through here, so the pause behaviour is identical on both
   * paths. While paused the callback re-arms the hold instead of starting a
   * morph: a one-shot that merely returned would leave the chain with no
   * successor and permanently stall the carousel until the effect reloads.
   * Re-arming keeps the carousel halted for the whole pause and resumes it
   * within one hold once "Pause Animation" is released.
   */
  void arm_hold_timer() {
    timeline.add(
        0, Animation::PeriodicTimer(
               NO_MORPH_FRAMES,
               [this](Canvas &) {
                 if (animationsPaused()) {
                   arm_hold_timer();
                   return;
                 }
                 start_morph();
               },
               false));
  }

  /**
   * @brief Advances to the next solid and starts the cross-morph.
   * @details Compiles the incoming solid into the carousel's back slot, runs a
   * MeshMorph from front to back over MORPH_FRAMES, and on completion swaps the
   * front, compacts, and re-arms the hold timer for the following morph.
   */
  void start_morph() {
    auto solids = Solids::Collections::get_platonic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int new_slot = 1 - carousel.front_index();

    carousel.slot(new_slot) = MeshState();
    PolyMesh poly = generate(persistent_arena, [&](Arena &target, Arena &a,
                                                   Arena &b) {
      return Solids::finalize_solid(solids[solid_idx].generate(a, b), target);
    });
    MeshOps::compile(poly, carousel.slot(new_slot), persistent_arena);

    morphing = true;
    timeline.add(
        0, Animation::MeshMorph(carousel.current(), carousel.slot(new_slot),
                                persistent_arena, draw_morph_fn_,
                                draw_morph_fn_, MORPH_FRAMES, ease_in_out_sin)
               .then([this, new_slot]() {
                 morphing = false;
                 carousel.set_front(new_slot);
                 carousel.compact();
                 // compact() reset the persistent arena; re-allocate the warp
                 // cache there (it re-populates on the next flush).
                 filters.template get<Filter::Pixel::Feedback<W, H>>()
                     .init_storage(persistent_arena);
                 arm_hold_timer();
               }));
  }

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
  NoiseParams noise_params;

  Orientation<> orientation;
  Timeline timeline;
  ProceduralPalette palette;

  // Mesh carousel + morph state (transitions are MeshMorph-driven, so the
  // carousel's segue is unused)
  MeshCarousel<> carousel;
  int solid_idx = 0;
  bool morphing = false;

  /** @brief Morph draw callback; Fn member gives FunctionRef a stable lifetime. */
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_fn_{
      [this](Canvas &c, const MeshState &m, float opacity) {
        Color4 shade = palette.get(0.0f);
        Plot::Mesh::draw<W, H>(filters, c, m,
                               [shade, opacity](const Vector &, Fragment &f) {
                                 f.color = shade;
                                 f.color.alpha *= opacity;
                               });
      }};

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H>>
      filters;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
