/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include "../generators.h"
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

    // Load first Catalan solid
    auto solids = Solids::Collections::get_catalan_solids();
    solid_idx = 0;
    carousel.load([&solids, this](Arena &a, Arena &b) {
      return solids[solid_idx].generate(a, b);
    });

    registerParam("Fade", &params.fade, 0.5f, 0.99f);
    registerParam("Distort Amp", &params.amplitude, 0.0f, 30.0f);
    registerParam("Distort Freq", &params.frequency, 0.01f, 1.0f);
    registerParam("Distort Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Noise Scale", &params.scale, 0.1f, 50.0f);
    registerParam("Pause Presets", &preset_paused);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));

    // Preset cycling
    timeline.add(0, Animation::PeriodicTimer(
                        150,
                        [this](Canvas &) {
                          if (preset_paused)
                            return;
                          presets.next();
                          timeline.add(
                              0, Animation::Lerp(params, presets.prev_get(),
                                                 presets.get(), 48, ease_mid));
                        },
                        true));

    // Shape cycling
    timeline.add(0,
                 Animation::Sprite([](Canvas &, float) {}, 16).then([this]() {
                   start_morph();
                 }));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    apply_params();
    noise_params.sync();
    filters.next.next.set_fade(params.fade);

    filters.flush(
        canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); },
        1.0f);

    timeline.step(canvas);

    if (carousel.current().is_bound()) {
      Plot::Mesh::draw<W, H>(filters, canvas, carousel.current(),
                             [&](const Vector &v, Fragment &f) {
                               float t_val = (v.y + 1.0f) * 0.5f;
                               f.color = palette.get(t_val);
                             });
    }
  }

private:
  void start_morph() {
    constexpr int MORPH_FRAMES = 32;
    auto solids = Solids::Collections::get_catalan_solids();
    solid_idx = (solid_idx + 1) % solids.size();

    int new_slot = 1 - carousel.front_index();

    carousel.slot(new_slot) = MeshState();
    active_mesh_A = MeshState();
    active_mesh_B = MeshState();

    // Generate incoming shape into back slot
    {
      scratch_arena_a.reset();
      scratch_arena_b.reset();
      ScopedScratch _a(scratch_arena_a);
      ScopedScratch _b(scratch_arena_b);
      PolyMesh poly =
          solids[solid_idx].generate(scratch_arena_a, scratch_arena_b);
      carousel.slot(new_slot).clear();
      MeshOps::compile(poly, carousel.slot(new_slot), persistent_arena);
    }

    // Preallocate morph buffer
    size_t max_v = std::max(carousel.current().vertices.size(),
                            carousel.slot(new_slot).vertices.size());
    morph_buffer.preallocate(persistent_arena, max_v);

    morphing = true;
    timeline.add(
        0, Animation::MeshMorph(
               &active_mesh_A, &active_mesh_B, &morph_buffer, &persistent_arena,
               carousel.current(), carousel.slot(new_slot), draw_morph_fn_,
               draw_morph_fn_, MORPH_FRAMES, false, ease_in_out_sin)
               .then([this, new_slot]() {
                 carousel.set_front(new_slot);
                 morphing = false;
                 // Morph complete: drop temporary buffers and compact
                 active_mesh_A = MeshState();
                 active_mesh_B = MeshState();
                 carousel.incoming() = MeshState();
                 morph_buffer = Animation::MorphBuffer();
                 {
                   Persist<MeshState> p(carousel.current(), scratch_arena_a,
                                        persistent_arena);
                   persistent_arena.reset();
                 }

                 // Rest on the static shape for a moment before morphing again
                 timeline.add(0, Animation::Sprite([](Canvas &, float) {}, 16)
                                     .then([this]() { start_morph(); }));
               }));
  }

  void apply_params() {
    noise_params.amplitude = params.amplitude;
    noise_params.frequency = params.frequency;
    noise_params.speed = params.speed;
    noise_params.scale = params.scale;
  }

  Presets<Params, 4> presets = {
      .entries = {{{"Flames", {0.94f, 0.51f, 0.42f, 0.46f, 23.0f}},
                   {"Organic", {0.93f, 0.9f, 0.33f, 0.8f, 22.6f}},
                   {"Intense", {0.96f, 5.0f, 0.15f, 1.5f, 10.0f}},
                   {"Subtle", {0.90f, 0.3f, 0.60f, 0.3f, 35.0f}}}},
      .current_idx = 0};
  bool preset_paused = true;
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
  float morph_alpha = 0.0f;

  /// Draw callback for morph transitions — stored as member so FunctionRef
  /// in MeshMorph can safely reference it (effect outlives animation).
  void draw_morph_mesh(Canvas &c, const MeshState &mesh, float opacity) {
    Plot::Mesh::draw<W, H>(filters, c, mesh,
                           [this, opacity](const Vector &v, Fragment &f) {
                             float t_val = (v.y + 1.0f) * 0.5f;
                             f.color = palette.get(t_val);
                             f.color.alpha *= opacity;
                           });
  }
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_morph_mesh(c, m, o);
      }};

  struct TransformerFn {
    const MeshFeedback *self;
    Vector operator()(const Vector &v) const {
      return noise_transform(v, self->noise_params).normalize();
    }
  };

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H, TransformerFn>>
      filters;
};
