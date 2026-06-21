/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <numeric>

#include "core/effects_engine.h"

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
  static constexpr int PRESET_FRAMES = 241; // preset hard-cut period (prime)
  // The shape cycle is a no-morph hold followed by a morph, so its PERIOD is
  // hold + morph — and that, not the hold alone, must be coprime with
  // PRESET_FRAMES for the shape and preset cycles to stay out of phase. 241 is
  // prime, so any period that is not a multiple of it is coprime; 240
  // (= PRESET_FRAMES - 1) is coprime and maximises the beat between the cycles.
  static constexpr int SHAPE_FRAMES = 240; // shape-cycle period, coprime with PRESET_FRAMES
  static constexpr int NO_MORPH_FRAMES = SHAPE_FRAMES - MORPH_FRAMES; // 80-frame hold

  // Lock the coprimality the comment above relies on: if a future edit to either
  // period reintroduces a common factor, the two cycles re-lock and the drift is
  // lost. Trap that at compile time rather than shipping a silently re-phased
  // carousel (the smoke harness cannot observe the lost beat).
  static_assert(std::gcd(SHAPE_FRAMES, PRESET_FRAMES) == 1,
                "SHAPE_FRAMES and PRESET_FRAMES must stay coprime so the shape "
                "and preset cycles drift out of phase instead of locking");

  /**
   * @brief Wires up palette, noise, orientation, and the filter pipeline.
   * @details Constructs the World/Screen/Pixel filter stack; the Feedback pixel
   * filter reads `style` by reference. The mem-initializer list below names
   * `filters` (which captures &style) without naming `style`, which reads as if
   * filters were built first — but C++ runs initializers in member-DECLARATION
   * order, and `style` is declared well before `filters`, so style is fully
   * constructed before filters binds its reference. Keep style declared ahead of
   * filters.
   */
  FLASHMEM MeshFeedback()
      : Effect(W, H), noise_params(), orientation(), timeline(),
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
    // Bind mutable state into all presets
    for (auto &e : presets.entries) {
      e.params.noise = &noise_params;
    }
    style = presets.get();
    apply_params();

    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.sync();

    // Load first shape
    solid_idx = 0;
    {
      auto solids = Solids::Collections::get_platonic_solids();
      PolyMesh poly = generate(persistent_arena, [&](Arena &target, Arena &a,
                                                     Arena &b) {
        return Solids::finalize_solid(solids[solid_idx].generate(a, b), target);
      });
      carousel.slot(carousel.front_index()).clear();
      MeshOps::compile(poly, carousel.slot(carousel.front_index()),
                       persistent_arena);
    }

    registerParam("Fade", &style.fade, 0.0f, 0.99f);
    registerParam("Distort Amp", &style.amplitude, 0.0f, 30.0f);
    registerParam("Distort Freq", &style.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &style.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &style.scale, 0.1f, 50.0f);
    registerParam("Hue Shift", &style.hue_shift, 0.0f, 0.1f);
    registerParam("Feedback", &feedback_enabled);
    // The preset cycle drives these six style params; flag them so the standard
    // "Pause Animation" toggle gates the cycling.
    markAnimated("Fade");
    markAnimated("Distort Amp");
    markAnimated("Distort Freq");
    markAnimated("Distort Speed");
    markAnimated("Noise Scale");
    markAnimated("Hue Shift");

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));

    // Preset cycling — presets hard-cut: snap style to the new preset (no
    // crossfade). presets.apply() copies all scalar/fn-ptr/downsample fields
    // and leaves the bound noise pointer intact.
    timeline.add(0, Animation::PeriodicTimer(
                        PRESET_FRAMES,
                        [this](Canvas &) {
                          if (animationsPaused())
                            return;
                          presets.next();
                          presets.apply(style);
                        },
                        true));

    // Shape cycling
    timeline.add(
        0, Animation::PeriodicTimer(
               NO_MORPH_FRAMES, [this](Canvas &) { start_morph(); }, false));
  }

  /**
   * @brief Suppresses the engine background clear; the feedback flush owns the
   * frame.
   * @return Always false to disable the background clear.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Renders one frame.
   * @details Applies params, runs the feedback decay flush, draws the current
   * mesh (skipped while morphing — the morph animation draws instead), then
   * advances the timeline.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    apply_params();

    // feedback step
    filters.flush(
        canvas, [](float, float, float) { return Color4(0, 0, 0, 0); },
        1.0f);

    if (!morphing && carousel.current().is_bound()) {
      Plot::Mesh::draw<W, H>(filters, canvas, carousel.current(),
                             [&](const Vector &v, Fragment &f) {
                               float t_val = (v.y + 1.0f) * 0.5f;
                               f.color = palette.get(t_val);
                             });
    }

    timeline.step(canvas);
  }

private:
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

    // Generate incoming shape into back slot
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
                 // One-shot hold timer before the next morph; non-repeating so
                 // no draw_fn fires on the idle hold frames.
                 timeline.add(
                     0, Animation::PeriodicTimer(
                            NO_MORPH_FRAMES,
                            [this](Canvas &) { start_morph(); }, false));
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

  Presets<Style, 7> presets = {.entries = {{{Style::Melting()},
                                            {Style::Swirling()},
                                            {Style::Churn()},
                                            {Style::Smoke()},
                                            {Style::Frozen()},
                                            {Style::Shatter()},
                                            {Style::Drift()}}},
                               .current_idx = 0};
  bool feedback_enabled = true;
  NoiseParams noise_params;

  Orientation<> orientation;
  Timeline timeline;
  ProceduralPalette palette;

  // Mesh carousel + morph state
  MeshCarousel carousel;
  int solid_idx = 0;
  bool morphing = false;

  /** @brief Morph draw callback; Fn member gives FunctionRef a stable lifetime. */
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_fn_{
      [this](Canvas &c, const MeshState &m, float opacity) {
        Plot::Mesh::draw<W, H>(filters, c, m,
                               [this, opacity](const Vector &v, Fragment &f) {
                                 float t_val = (v.y + 1.0f) * 0.5f;
                                 f.color = palette.get(t_val);
                                 f.color.alpha *= opacity;
                               });
      }};

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H>>
      filters;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
