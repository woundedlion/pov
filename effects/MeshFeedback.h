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

    void lerp(const Params &a, const Params &b, float t) {
      fade = ::lerp(a.fade, b.fade, t);
      amplitude = ::lerp(a.amplitude, b.amplitude, t);
      frequency = ::lerp(a.frequency, b.frequency, t);
      speed = ::lerp(a.speed, b.speed, t);
      scale = ::lerp(a.scale, b.scale, t);
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
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H, TransformerFn>(
                    TransformerFn{this}, 0.95f)) {}

  void init() override {
    params = presets.get();
    apply_params();

    noise_params.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_params.sync();

    // Load first shape
    solid_idx = 0;
    generate_incoming_shape(carousel.front_index());

    registerParam("Fade", &params.fade, 0.0f, 0.99f);
    registerParam("Distort Amp", &params.amplitude, 0.0f, 30.0f);
    registerParam("Distort Freq", &params.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &params.scale, 0.1f, 50.0f);
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
  void generate_incoming_shape(int slot_idx) {
    auto solids = Solids::Collections::get_platonic_solids();
    PolyMesh poly = generate(persistent_arena,
        [&](Arena &target, Arena &a, Arena &b) {
          return Solids::finalize_solid(
              solids[solid_idx].generate(a, b), target);
        });
    carousel.slot(slot_idx).clear();
    MeshOps::compile(poly, carousel.slot(slot_idx), persistent_arena);
  }

  void prepare_morph_buffers(int new_slot) {
    size_t max_v = std::max(carousel.current().vertices.size(),
                            carousel.slot(new_slot).vertices.size());
    morph_buffer.preallocate(persistent_arena, max_v);
  }

  void release_morph_transients() {
    active_mesh_A = MeshState();
    active_mesh_B = MeshState();
    carousel.incoming() = MeshState();
    morph_buffer = Animation::MorphBuffer();
  }

  void compact_persistent_data() {
    Persist<MeshState> p(carousel.current(), scratch_arena_a, persistent_arena);
    persistent_arena.reset();
  }

  void start_morph() {

    auto solids = Solids::Collections::get_platonic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int new_slot = 1 - carousel.front_index();

    release_morph_transients();
    generate_incoming_shape(new_slot);
    prepare_morph_buffers(new_slot);

    morphing = true;
    timeline.add(
        0, Animation::MeshMorph(
               &active_mesh_A, &active_mesh_B, &morph_buffer, &persistent_arena,
               carousel.current(), carousel.slot(new_slot), draw_morph_fn_,
               draw_morph_fn_, MORPH_FRAMES, false, ease_in_out_sin)
               .then([this, new_slot]() {
                 carousel.set_front(new_slot);
                 morphing = false;

                 release_morph_transients();
                 compact_persistent_data();

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
    filters.next.next.set_fade(feedback_enabled ? params.fade : 0.0f);
  }

  Presets<Params, 4> presets = {
      .entries = {{{"Flames", {0.9f, 0.51f, 0.42f, 0.46f, 23.0f}},
                   {"Zebras1", {0.58f, 2.73f, 0.07f, 0.0f, 26.0f}},
                   {"Zebras2", {0.58f, 8.21f, 0.01f, 0.0f, 46.0f}},
                   {"Interdimensional", {0.68f, 4.98f, 0.07f, 0.2f, 5.0f}}}},
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
  MeshState active_mesh_A;
  MeshState active_mesh_B;
  Animation::MorphBuffer morph_buffer;
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

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H, TransformerFn>>
      filters;
};

#include "../effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
