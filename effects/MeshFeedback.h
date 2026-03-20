/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../solids.h"

template <int W, int H> class MeshFeedback : public Effect {
public:
  struct Params {
    float fade = 0.93f;
    float amplitude = 0.9f;
    float frequency = 0.33f;
    float speed = 0.8f;
    float scale = 22.6f;
    float hue_shift = 0.03f;

    void lerp(const Params &a, const Params &b, float t) {
      fade = ::lerp(a.fade, b.fade, t);
      amplitude = ::lerp(a.amplitude, b.amplitude, t);
      frequency = ::lerp(a.frequency, b.frequency, t);
      speed = ::lerp(a.speed, b.speed, t);
      scale = ::lerp(a.scale, b.scale, t);
      hue_shift = ::lerp(a.hue_shift, b.hue_shift, t);
    }
  } params;

  static constexpr int MORPH_FRAMES = 160;
  static constexpr int NO_MORPH_FRAMES =
      81; // morph period 241 is prime → coprime with preset cycle
  static constexpr int LERP_FRAMES = 0;
  static constexpr int NO_LERP_FRAMES = 241;

  FLASHMEM MeshFeedback()
      : Effect(W, H), noise_params(), orientation(), timeline(),
        palette(Palettes::peachPop),
        filters(
            Filter::World::Orient<W>(orientation),
            Filter::Screen::AntiAlias<W, H>(),
            Filter::Pixel::Feedback<W, H, TransformerFn, HueShiftFade>(
                TransformerFn{this}, 0.95f, HueShiftFade{&params.hue_shift})) {}

  void init() override {
    params = presets.get();
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

    registerParam("Fade", &params.fade, 0.0f, 0.99f);
    registerParam("Distort Amp", &params.amplitude, 0.0f, 30.0f);
    registerParam("Distort Freq", &params.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &params.scale, 0.1f, 50.0f);
    registerParam("Hue Shift", &params.hue_shift, 0.0f, 0.1f);
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
                              0, Animation::Lerp(params, presets.prev_get(),
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
    noise_params.amplitude = params.amplitude;
    noise_params.frequency = params.frequency;
    noise_params.speed = params.speed;
    noise_params.scale = params.scale;
    noise_params.sync();
    filters.template get<FeedbackFilter>().set_fade(feedback_enabled ? params.fade : 0.0f);
  }

  Presets<Params, 5> presets = {
      .entries = {{{{0.82f, 3.15f, 1.0f, 0.02f, 1.0f, 0.035f}},
                   {{0.9f, 0.51f, 0.42f, 0.46f, 23.0f, 0.01f}},
                   {{0.58f, 2.73f, 0.07f, 0.0f, 26.0f, 0.03f}},
                   {{0.58f, 8.21f, 0.01f, 0.0f, 46.0f, 0.03f}},
                   {{0.68f, 4.98f, 0.07f, 0.2f, 5.0f, 0.03f}}}},
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

  struct TransformerFn {
    const MeshFeedback *self;
    Vector operator()(const Vector &v) const {
      return noise_transform(v, self->noise_params);
    }
  };

  struct HueShiftFade {
    const float *hue_amount;
    Pixel operator()(const Pixel &p, float fade) const {
      return hue_rotate(Color4(p * fade, 1.0f), *hue_amount).color;
    }
  };

  using FeedbackFilter = Filter::Pixel::Feedback<W, H, TransformerFn, HueShiftFade>;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           FeedbackFilter>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
