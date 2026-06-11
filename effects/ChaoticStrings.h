/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

template <int W, int H> class ChaoticStrings : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115;
  // Each trail frame is an Orientation holding up to ORIENTATION_SUBSTEPS
  // sub-frames; deep_tween() emits one Fragment per (frame, sub-frame). The
  // product is the strict worst-case vertex count per draw, and it sizes both
  // the scratch vertex buffer and the scratch arena carve below.
  static constexpr int ORIENTATION_SUBSTEPS = 16;
  static constexpr int MAX_FRAGMENTS = TRAIL_LENGTH * ORIENTATION_SUBSTEPS;

  struct Params {
    float alpha = 1.0f;
    float cycle_duration = 80.0f;
    float jitterAmp = 1.7f;
    float speed = 0.04f;
    float noiseFreq = 0.32f;
    float scaleFactor = 200.0f;
    float cycleSpeed = 0.1f;
  } params;

  struct Node {
    Orientation<ORIENTATION_SUBSTEPS> orientation;
    Animation::OrientationTrail<Orientation<ORIENTATION_SUBSTEPS>,
                                TRAIL_LENGTH>
        trail;
    Vector v;

    Node() : v(Y_AXIS) {}
  };

  FLASHMEM ChaoticStrings()
      : Effect(W, H), timeline(), filters(Filter::Screen::AntiAlias<W, H>()),
        path([](float) { return Vector(0, 1, 0); }), orientation(),
        palette_variant(), cycle_phase(0.0f), noise_xform(timeline) {}

  // Scratch A holds the per-frame Fragment buffer (MAX_FRAGMENTS) plus headroom
  // for the downstream Multiline draw. Tie the carve to that worst case so the
  // buffer can never overrun its arena, with margin above the bare requirement.
  static constexpr size_t SCRATCH_A_BYTES = 200 * 1024;
  static_assert(SCRATCH_A_BYTES >= MAX_FRAGMENTS * sizeof(Fragment),
                "scratch arena A must fit the worst-case fragment buffer");

  void init() override {

    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES, SCRATCH_A_BYTES, 0);

    node = static_cast<Node *>(
        persistent_arena.allocate(sizeof(Node), alignof(Node)));
    new (node) Node();

    // Colors
    static_palette.bind(&palette_variant, &scale_mod, &cycle_mod);
    palette_variant = Palettes::fireAndIce;

    // Parameters
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Jitter Amp", &params.jitterAmp, 0.0f, 10.0f);
    registerParam("Noise Freq", &params.noiseFreq, 0.01f, 10.0f);
    registerParam("Scale Factor", &params.scaleFactor, 1.0f, 500.0f);
    registerParam("Cycle Speed", &params.cycleSpeed, 0.0f, 1.0f);

    // Initialize noise params
    noise_xform.params.amplitude = params.jitterAmp;
    noise_xform.params.frequency = params.noiseFreq;
    noise_xform.params.speed = params.speed;
    noise_xform.params.sync();

    // Build the initial Lissajous path and start the noise transformer
    update_path();
    noise_xform.spawn(0, -1);

    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    // Retain the motion handle so the Cycle Dur slider can be applied live
    // (its duration was otherwise captured once at init).
    motion_ = timeline.add_get(
        0, Animation::Motion<W, ORIENTATION_SUBSTEPS>(
               node->orientation, path, (int)params.cycle_duration, true));

    // Bind the Driver to the live Cycle Speed slider (scale 1.0) so it reads the
    // value itself each step; no retained handle or per-frame set_speed needed
    // (the idiom DreamBalls/Liquid2D use).
    timeline.add(0, Animation::Driver(cycle_phase, &params.cycleSpeed, 1.0f));

    last_cycle_duration_ = params.cycle_duration;
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    ScratchScope _(scratch_arena_a);
    ArenaVector<Fragment> vertices(scratch_arena_a, MAX_FRAGMENTS);
    timeline.step(canvas);

    // Push live parameter changes into the active noise entities. The
    // FastNoiseLite re-sync is left to prepare_frame() (which calls
    // NoiseParams::sync on every active entity) rather than hand-rolled here.
    for (auto &e : noise_xform.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
      }
    }
    noise_xform.prepare_frame();

    // Live-apply the Cycle Dur slider to the motion
    // (guarded: set_duration reschedules from now, so calling it every frame
    // would perpetually defer the trigger).
    apply_if_changed(params.cycle_duration, last_cycle_duration_, [&](float cd) {
      if (motion_)
        motion_->set_duration((int)cd);
    });

    // Record the current orientation snapshot
    node->trail.record(node->orientation);

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      // Re-apply transforms at draw time so the trail undulates with noise
      Vector pos =
          noise_xform.transform(orientation.orient(rotate(node->v, q)));
      Fragment f;
      f.pos = normalized_or(pos, Vector(1, 0, 0));
      f.v3 = t;
      vertices.push_back(f);
    });

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      // Drive color and fade from the normalized trail parameter (v3), as the
      // sibling Comets effect does.
      float color_t = frag.v3;
      frag.color = static_palette.get(color_t);
      frag.color.alpha *= quintic_kernel(frag.v3) * params.alpha;
    };

    Plot::Multiline::draw<W, H>(filters, canvas, vertices, fragment_shader);
  }

private:
  void update_path() {
    static constexpr LissajousParams config{12.0f, 5.0f, 0, 2 * PI_F};
    path.f = [](float t) {
      return lissajous(config.m1, config.m2, config.a, t * config.domain);
    };
  }

  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ProceduralPath path;
  Orientation<> orientation;
  ScaleModifier scale_mod{200.0f, &params.scaleFactor};
  CycleModifier cycle_mod{&cycle_phase};
  ProceduralPalette palette_variant;
  StaticPalette<ProceduralPalette, Coords<ScaleModifier, CycleModifier>>
      static_palette;
  Animation::Motion<W, ORIENTATION_SUBSTEPS> *motion_ = nullptr;
  float last_cycle_duration_ = -1.0f;
  float cycle_phase = 0.0f;
  Node *node = nullptr;
  NoiseTransformer<1> noise_xform;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(ChaoticStrings)
