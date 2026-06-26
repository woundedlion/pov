/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Undulating string of orientation trails drawn as an anti-aliased
 *        multiline colored by trail position.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A single node random-walks and follows a Lissajous path while a
 *          noise transformer warps its trail each frame.
 * @note Sibling trail effects `Comets` and `RingSpin` share only the
 *       record + deep_tween skeleton (in the engine); their bodies diverge in
 *       draw primitive, transform chain, and color/fade, so each renders
 *       independently and trail-rendering fixes must be propagated by hand.
 */
template <int W, int H> class ChaoticStrings : public Effect {
public:
  static constexpr int TRAIL_LENGTH = 115;
  static constexpr int ORIENTATION_SUBSTEPS = 16;
  static constexpr int MAX_FRAGMENTS = TRAIL_LENGTH * ORIENTATION_SUBSTEPS;

  /**
   * @brief Live-tunable parameters exposed as sliders.
   */
  struct Params {
    float alpha = 1.0f;          /**< Overall trail opacity in [0, 1]. */
    float cycle_duration = 80.0f; /**< Motion cycle duration in frames. */
    float jitterAmp = 1.7f;      /**< Noise displacement amplitude. */
    float speed = 0.04f;         /**< Noise field evolution speed. */
    float noiseFreq = 0.32f;     /**< Noise spatial frequency. */
    float scaleFactor = 200.0f;  /**< Palette coordinate scale factor. */
    float cycleSpeed = 0.1f;     /**< Palette cycle phase advance per step. */
  } params;

  /**
   * @brief The single animated body: orientation, rolling trail history, and
   *        base direction vector.
   */
  struct Node {
    Orientation<ORIENTATION_SUBSTEPS> orientation; /**< Current node orientation. */
    Animation::OrientationTrail<Orientation<ORIENTATION_SUBSTEPS>,
                                TRAIL_LENGTH>
        trail; /**< Rolling history of orientations forming the drawn trail. */
    Vector v;  /**< Base direction vector, seeded to +Y. */

    /**
     * @brief Constructs a Node with its base direction set to the Y axis.
     */
    Node() : v(Y_AXIS) {}
  };

  /**
   * @brief Constructs the effect and its members.
   * @details The path seeds to +Y until update_path() installs the Lissajous
   *          curve in init().
   */
  FLASHMEM ChaoticStrings()
      : Effect(W, H), timeline(), filters(Filter::Screen::AntiAlias<W, H>()),
        path([](float) { return Vector(0, 1, 0); }), orientation(),
        palette_variant(), cycle_phase(0.0f), noise_xform(timeline) {}

  // Scratch A holds the per-frame vertices buffer and, during the draw call, the
  // Multiline fragment buffer it binds (capacity vertices.size()+2) at the same
  // time, so the worst case is both buffers live at once.
  static constexpr size_t SCRATCH_A_BYTES = 200 * 1024;
  static_assert(SCRATCH_A_BYTES >= (2 * MAX_FRAGMENTS + 2) * sizeof(Fragment),
                "scratch arena A must fit the vertices buffer and the "
                "Multiline-draw fragment buffer at once");

  /**
   * @brief Carves arenas, allocates the node, binds the palette, registers
   *        sliders, and wires up the timeline animations.
   * @details Sets up the random walk, path motion, and cycle driver animations.
   */
  void init() override {

    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES, SCRATCH_A_BYTES, 0);

    node = static_cast<Node *>(
        persistent_arena.allocate(sizeof(Node), alignof(Node)));
    new (node) Node();

    static_palette.bind(&palette_variant, &scale_mod, &cycle_mod);
    palette_variant = Palettes::fireAndIce;

    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);
    registerParam("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    registerParam("Speed", &params.speed, 0.0f, 5.0f);
    registerParam("Jitter Amp", &params.jitterAmp, 0.0f, 10.0f);
    registerParam("Noise Freq", &params.noiseFreq, 0.01f, 10.0f);
    registerParam("Scale Factor", &params.scaleFactor, 1.0f, 500.0f);
    registerParam("Cycle Speed", &params.cycleSpeed, 0.0f, 1.0f);

    noise_xform.template_params.amplitude = params.jitterAmp;
    noise_xform.template_params.frequency = params.noiseFreq;
    noise_xform.template_params.speed = params.speed;
    noise_xform.template_params.sync();

    update_path();
    noise_xform.spawn(0, -1);

    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    motion_ = timeline.add_get(
        0, Animation::Motion<W, ORIENTATION_SUBSTEPS>(
               node->orientation, path, (int)params.cycle_duration, true));

    timeline.add(0, Animation::Driver(cycle_phase, &params.cycleSpeed, 1.0f));

    last_cycle_duration_ = params.cycle_duration;
  }

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return true; the strip blanks to black immediately after each column is
   *         shown, so every column reads as a sharp slice with dark gaps
   *         between columns rather than persisting across the sweep.
   */
  bool strobe_columns() const override { return true; }

  /**
   * @brief Advances and renders one frame of the effect.
   * @details Steps the timeline, pushes live slider values into the noise
   *          entities, records the latest orientation into the trail, tweens it
   *          into per-vertex fragments warped by noise, and draws the colored
   *          multiline.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    ScratchScope scratch_a_guard(scratch_arena_a);
    ArenaVector<Fragment> vertices(scratch_arena_a, MAX_FRAGMENTS);
    timeline.step(canvas);

    // Push live slider values into the noise entities (they bind to `params`,
    // not the entities); prepare_frame() re-syncs each active entity afterward.
    for (auto &e : noise_xform.entities) {
      if (e.active) {
        e.params.frequency = params.noiseFreq;
        e.params.amplitude = params.jitterAmp;
        e.params.speed = params.speed;
      }
    }
    noise_xform.prepare_frame();

    apply_if_changed(params.cycle_duration, last_cycle_duration_, [&](float cd) {
      if (motion_)
        motion_->set_duration((int)cd);
    });

    node->trail.record(node->orientation);

    deep_tween(node->trail, [&](const Quaternion &q, float t) {
      Vector pos =
          noise_xform.transform(orientation.orient(rotate(node->v, q)));
      Fragment f;
      f.pos = normalized_or(pos, Vector(1, 0, 0));
      f.v3 = t;
      vertices.push_back(f);
    });

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      float color_t = frag.v3;
      frag.color = static_palette.get(color_t);
      frag.color.alpha *= quintic_kernel(frag.v3) * params.alpha;
    };

    Plot::Multiline::draw<W, H>(filters, canvas, vertices, fragment_shader);
  }

private:
  /**
   * @brief Installs the procedural path as a fixed Lissajous curve.
   * @details Samples the curve over its full domain.
   */
  void update_path() {
    static constexpr LissajousParams config{12.0f, 5.0f, 0, 2 * PI_F};
    path.f = [](float t) {
      return lissajous(config.m1, config.m2, config.a, t * config.domain);
    };
  }

  FastNoiseLite noise; /**< Noise source for the random walk. */
  Timeline timeline;   /**< Drives all per-frame animations. */
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters; /**< Anti-aliasing render pipeline. */
  ProceduralPath path;    /**< Lissajous path the node follows. */
  Orientation<> orientation; /**< Random-walk orientation reference frame. */
  ScaleModifier scale_mod{200.0f, &params.scaleFactor}; /**< Palette scale coordinate modifier. */
  CycleModifier cycle_mod{&cycle_phase}; /**< Palette cycle coordinate modifier. */
  ProceduralPalette palette_variant; /**< Active palette variant. */
  StaticPalette<ProceduralPalette, Coords<ScaleModifier, CycleModifier>>
      static_palette; /**< Bound palette sampled per fragment. */
  Animation::Motion<W, ORIENTATION_SUBSTEPS> *motion_ = nullptr; /**< Retained motion handle for live Cycle Dur updates. */
  float last_cycle_duration_ = -1.0f; /**< Last applied cycle duration, for change detection. */
  float cycle_phase = 0.0f; /**< Current palette cycle phase. */
  Node *node = nullptr;     /**< The single animated body, arena-allocated. */
  NoiseTransformer<1> noise_xform; /**< Warps the trail with noise each frame. */
};

#include "core/effect_registry.h"
REGISTER_EFFECT(ChaoticStrings)
