/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../styles.h"
#include "../solids.h"

template <int W, int H> class MeshFeedback : public Effect {
public:
  using Style = Feedback::Style;

  static constexpr int MORPH_FRAMES = 160;
  static constexpr int NO_MORPH_FRAMES =
      81; // morph period 241 is prime → coprime with preset cycle
  static constexpr int LERP_FRAMES = 0;
  static constexpr int NO_LERP_FRAMES = 241;

  FLASHMEM MeshFeedback()
      : Effect(W, H), noise_params(), orientation(), timeline(),
        palette(Palettes::peachPop),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Feedback::Filter<W, H>(style)) {}

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
    registerParam("Pause Presets", &preset_paused, preset_paused);
    registerParam("Feedback", &feedback_enabled, feedback_enabled);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));

    // Preset cycling
    timeline.add(0, Animation::PeriodicTimer(
                        NO_LERP_FRAMES + LERP_FRAMES,
                        [this](Canvas &) {
                          if (preset_paused)
                            return;
                          presets.next();
                          timeline.add(
                              0, Animation::Lerp(style, presets.prev_get(),
                                                 presets.get(), LERP_FRAMES,
                                                 ease_in_out_sin));
                        },
                        true));

    // Shape cycling
    timeline.add(
        0, Animation::PeriodicTimer(
               NO_MORPH_FRAMES, [this](Canvas &) { start_morph(); }, false));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    apply_params();

    // feedback step
    filters.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
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
                 timeline.add(0, Animation::Sprite([](Canvas &, float) {},
                                                   NO_MORPH_FRAMES)
                                     .then([this]() { start_morph(); }));
               }));
  }

  void apply_params() {
    style.sync_noise();
    filters.template get<Feedback::Filter<W, H>>().set_enabled(
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
  bool preset_paused = false;
  bool feedback_enabled = true;
  NoiseParams noise_params;

  Orientation<W> orientation;
  Timeline<W> timeline;
  ProceduralPalette palette;

  // Mesh carousel + morph state
  MeshCarousel<W> carousel;
  int solid_idx = 0;
  bool morphing = false;

  /// Morph draw callback - Fn member gives FunctionRef a stable lifetime
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
           Feedback::Filter<W, H>>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
