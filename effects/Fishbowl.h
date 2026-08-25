/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Fishbowl.h
 * @brief Undulating string of orientation trails drawn as an anti-aliased
 *        multiline.
 */

#include "core/engine/engine.h"
#include <array>
#include <string_view>

// Unit-test accessor reaching the private node so a test can count the vertices
// deep_tween emits against the MAX_FRAGMENTS the scratch split is sized for.
namespace hs_test {
namespace effects_tests {
struct FishbowlWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Undulating string of orientation trails drawn as an anti-aliased
 *        multiline colored by trail position.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A single node random-walks and follows a Lissajous path while a
 *          noise transformer warps its trail each frame.
 * @note Sibling trail effects `Comets` and `RingSpin` share the
 *       record + deep_tween skeleton (in the engine), and Comets shares the
 *       `Animation::TrailBody` aggregate at the same capacity and substep
 *       count. Draw primitive, transform chain and colour/fade are
 *       hand-propagated.
 */
template <int W, int H> class Fishbowl : public Effect {
public:
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"fire-trail"};
  static constexpr int TRAIL_LENGTH = 115;
  static constexpr int ORIENTATION_SUBSTEPS = 16;
  static constexpr int MAX_FRAGMENTS = 2 * TRAIL_LENGTH * ORIENTATION_SUBSTEPS;

  /** @brief Compact multiline control point carrying its palette coordinate. */
  struct TrailVertex {
    Vector pos;
    float palette_t;

    operator Fragment() const {
      Fragment fragment;
      fragment.pos = pos;
      fragment.v3 = palette_t;
      return fragment;
    }
  };

  /** @brief The single animated body: orientation, trail, base direction. */
  using Node = Animation::TrailBody<TRAIL_LENGTH, ORIENTATION_SUBSTEPS>;

  /**
   * @brief Constructs the effect and its members.
   * @details The path seeds to +Y until init() installs the Lissajous curve.
   */
  HS_COLD_MEMBER Fishbowl()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        timeline(), filters(Filter::Screen::AntiAlias<W, H>()),
        path([](float) { return Vector(0, 1, 0); }), orientation(),
        fire_palette({{0.00f, CPixel{0x000000}},
                      {0.18f, CPixel{0x260000}},
                      {0.40f, CPixel{0xD01000}},
                      {0.62f, CPixel{0xFF7800}},
                      {0.78f, CPixel{0xFFE080}},
                      {0.84f, CPixel{0xFFF8E8}},
                      {1.00f, CPixel{0x000000}}}),
        cycle_phase(0.0f), noise_xform(timeline) {}

  // Scratch A holds the per-frame vertices buffer and, during the draw call, the
  // Multiline fragment buffer it binds (capacity vertices.size()+2) plus
  // rasterize's own sub-step cache, so the worst case is all three live at once.
  static constexpr size_t SCRATCH_A_BYTES = 234 * 1024;
  static_assert(SCRATCH_A_BYTES >= MAX_FRAGMENTS * sizeof(TrailVertex) +
                                       (MAX_FRAGMENTS + 2) * sizeof(Fragment) +
                                       Plot::rasterize_scratch_a_bytes<W>(),
                "scratch arena A must fit the vertices buffer, the "
                "Multiline-draw fragment buffer and rasterize's sub-step cache "
                "at once");

  // init() allocates the noise transformer pool and the single Node
  // (Orientation + OrientationTrail) from the persistent partition;
  // static_palette binds in place with no arena storage.
  using NoiseEntity = typename NoiseTransformer<1>::Entity;
  static constexpr size_t FOOTPRINT_BYTES = sizeof(Node) + sizeof(NoiseEntity) +
                                            alignof(NoiseEntity) + sizeof(int) +
                                            alignof(int);
  // The custom split leaves the remainder of the device arena persistent.
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES,
      "Fishbowl persistent footprint exceeds its partition; "
      "retune TRAIL_LENGTH/ORIENTATION_SUBSTEPS or enlarge the split");

  /**
   * @brief Carves arenas, allocates the node, binds the palette, registers
   *        sliders, installs the Lissajous path, and wires up the timeline
   *        animations.
   * @details Sets up the random walk, path motion, and cycle driver animations.
   */
  void init() override {
    configure_presets(PRESET_IDS.size());
    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_A_BYTES, SCRATCH_A_BYTES, 0);

    noise_xform.init_storage(persistent_arena);
    node = persistent_arena.make<Node>();

    static_palette.bind(&fire_palette, &scale_mod, &cycle_mod, &duty_mod);

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Cycle Dur", &params.cycle_duration, 10.0f, 200.0f);
    register_param("Speed", &params.speed, 0.0f, 5.0f);
    register_param("Jitter Amp", &params.jitter_amp, 0.0f, 10.0f);
    register_param("Noise Scale", &params.noise_freq, 0.01f, 10.0f);
    register_param("Scale Factor", &params.scale_factor, 1.0f, 500.0f);
    register_param("Cycle Speed", &params.cycle_speed, 0.0f, 1.0f);
    register_param("Duty Cycle", &params.duty_cycle, 0.0f, 1.0f);

    noise_xform.template_params.amplitude = params.jitter_amp;
    noise_xform.template_params.frequency = params.noise_freq;
    noise_xform.template_params.speed = params.speed;

    // m2 * domain = 5 * 2*PI: path.f(1) == path.f(0), so the curve closes.
    static constexpr LissajousParams PATH_CONFIG{12.0f, 5.0f, 0, 2 * PI_F};
    path.f = [](float t) {
      return lissajous(PATH_CONFIG.m1, PATH_CONFIG.m2, PATH_CONFIG.a,
                       t * PATH_CONFIG.domain);
    };

    auto *warp = noise_xform.spawn_pinned(0, -1);
    HS_CHECK(warp, "Fishbowl: pinned noise spawn must succeed");

    timeline.add(0,
                 Animation::RandomWalk<W>(orientation, random_vector(), noise));
    motion = timeline.add_get(
        0,
        Animation::Motion<W, ORIENTATION_SUBSTEPS>(
            node->orientation, path, (int)params.cycle_duration, true),
        Timeline::Pin::PINNED);

    timeline.add(0, Animation::Driver(cycle_phase, &palette_speed, 1.0f));

    last_cycle_duration = (int)params.cycle_duration;
  }

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
    ArenaVector<TrailVertex> vertices(scratch_arena_a, MAX_FRAGMENTS);
    palette_speed = palette_phase_speed(params.cycle_speed, params.scale_factor,
                                        node->trail.length());
    {
      HS_PROFILE(fish_timeline_step);
      timeline.step(canvas);
    }

    // Push live slider values onto the noise template; prepare_frame() copies
    // them into each active entity (refresh_from) and re-syncs the generator.
    {
      HS_PROFILE(fish_noise_prepare);
      noise_xform.template_params.frequency = params.noise_freq;
      noise_xform.template_params.amplitude = params.jitter_amp;
      noise_xform.template_params.speed = params.speed;
      noise_xform.prepare_frame();
    }

    apply_if_changed((int)params.cycle_duration, last_cycle_duration,
                     [&](int cd) {
                       if (motion)
                         motion->set_duration(cd);
                     });

    node->trail.record(node->orientation);

    // A sub-LSB alpha paints nothing; skip building and rasterizing. The trail
    // is still recorded above, so motion resumes when alpha rises.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;

    const float fill_scale = palette_fill_scale(node->trail.length());

    {
      HS_PROFILE(fish_build_vertices);
      Quaternion previous_q;
      Vector previous_pos;
      float previous_t = 0.0f;
      bool have_previous = false;
      deep_tween(node->trail, [&](const Quaternion &q, float t) {
        const Vector pos = warped_position(q);
        if (have_previous) {
          const Quaternion mid_q = slerp(previous_q, q, 0.5f);
          const Vector mid_pos = warped_position(mid_q);
          if (needs_adaptive_midpoint(previous_pos, mid_pos, pos)) {
            vertices.push_back({mid_pos, 0.5f * (previous_t + t) * fill_scale});
          }
        }
        TrailVertex f;
        f.pos = pos;
        f.palette_t = t * fill_scale;
        vertices.push_back(f);
        previous_q = q;
        previous_pos = pos;
        previous_t = t;
        have_previous = true;
      });
    }

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      const float age_t = frag.v3 / fill_scale;
      frag.color = shade_trail(frag.v3, age_t);
    };

    {
      HS_PROFILE(fish_multiline_draw);
      Plot::Multiline::draw<W, H>(filters, canvas, vertices, fragment_shader);
    }
  }

private:
  HS_COLD_MEMBER bool apply_preset(const PresetChange &) override {
    params = PRESET;
    return true;
  }

  friend struct ::hs_test::effects_tests::FishbowlWhiteBox;

  /**
   * @brief Palette-coordinate modifier lighting only part of each cycle.
   * @details Compresses the whole gradient into the leading `duty_cycle`
   * fraction of every cycle; the remainder holds t = 1, the gradient's black
   * stop.
   */
  struct DutyCycleModifier {
    static constexpr bool bounded_output = true;
    static constexpr bool rebounds_input = true;

    const float *duty_cycle;

    float modify(float t) const {
      const float u = wrap_t(t);
      const float duty = hs::clamp(*duty_cycle, 0.0f, 1.0f);
      return duty > 0.0f && u < duty ? u / duty : 1.0f;
    }
  };

  /**
   * @brief Palette-coordinate scale for a partly filled trail.
   * @param trail_length Recorded trail length after this frame's record().
   * @return Recorded length over capacity, so a sample's coordinate depends on
   * its absolute slot rather than its normalized position.
   */
  static float palette_fill_scale(size_t trail_length) {
    return static_cast<float>(trail_length) / TRAIL_LENGTH;
  }

  /**
   * @brief Cycle-phase advance per frame, compensated while the trail fills.
   * @param cycle_speed Slider phase advance per frame.
   * @param scale_factor Palette coordinate scale factor.
   * @param trail_length Recorded trail length *before* this frame's record().
   * @return Phase advance to drive `cycle_phase` with this frame.
   * @details A full trail evicts its oldest sample each frame, dropping every
   * survivor one slot and shifting its palette coordinate for free; while
   * filling, slots are static and the driver must supply that shift itself.
   * Hence the pre-record length here against palette_fill_scale()'s
   * post-record one: the frame that fills the last slot still evicts nothing
   * and must stay compensated.
   */
  static float palette_phase_speed(float cycle_speed, float scale_factor,
                                   size_t trail_length) {
    if (trail_length < TRAIL_LENGTH)
      return cycle_speed - scale_factor / TRAIL_LENGTH;
    return cycle_speed;
  }

  /**
   * @brief Trail sample oriented into the view frame and warped by the noise
   *        field.
   * @param q Recorded trail orientation.
   * @return Unit position on the sphere.
   */
  Vector warped_position(const Quaternion &q) const {
    const Vector pos =
        noise_xform.transform(orientation.orient(rotate(node->v, q)));
    return normalized_or(pos, Vector(1, 0, 0));
  }

  /**
   * @brief True when a warped segment bows far enough to need a midpoint
   *        vertex.
   * @param a Segment start.
   * @param mid Warped midpoint of the segment.
   * @param b Segment end.
   * @return True once `mid` strays a quarter pixel-column from the geodesic
   * midpoint of `a` and `b`.
   */
  static bool needs_adaptive_midpoint(const Vector &a, const Vector &mid,
                                      const Vector &b) {
    constexpr float MAX_ERROR = 0.25f * 2.0f * PI_F / W;
    const Vector geodesic_mid = slerp(a, b, 0.5f);
    return angle_between(mid, geodesic_mid) > MAX_ERROR;
  }

  /**
   * @brief Samples the fire palette and applies the trail fade.
   * @param palette_t Fill-scaled palette coordinate.
   * @param age_t Normalized trail age, 0 at the tail and 1 at the head.
   * @return Palette colour faded by the age kernel and master alpha.
   * @details The gradient's black stops read as transparent, not as black
   * over the background.
   */
  Color4 shade_trail(float palette_t, float age_t) const {
    Color4 color = static_palette.get(palette_t);
    if (color.color == Pixel())
      color.alpha = 0.0f;
    color.alpha *= quintic_kernel(age_t) * params.alpha;
    return color;
  }

  FastNoiseLite noise; /**< Noise source for the random walk. */
  Timeline timeline;   /**< Drives all per-frame animations. */
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>>
      filters;               /**< Anti-aliasing render pipeline. */
  ProceduralPath path;       /**< Lissajous path the node follows. */
  Orientation<> orientation; /**< Random-walk orientation reference frame. */

  /**
   * @brief Live-tunable parameters exposed as sliders.
   */
  struct Params {
    float alpha;             /**< Overall trail opacity in [0, 1]. */
    float cycle_duration;    /**< Motion cycle duration in frames. */
    float speed;             /**< Noise field evolution speed. */
    float jitter_amp;        /**< Noise displacement amplitude. */
    float noise_freq;        /**< Noise spatial frequency. */
    float scale_factor;      /**< Palette coordinate scale factor. */
    float cycle_speed;       /**< Palette cycle phase advance per step. */
    float duty_cycle = 0.5f; /**< Lit fraction of each palette cycle. */
  };

  // All 8 Params fields are listed explicitly; a trailing omission would
  // silently fall back to the default member initializer, not the preset value.
  static constexpr Params PRESET{1.0f,     80.0f,      0.235f, 2.88f,
                                 0.26974f, 84.832001f, 0.672f, 0.5f};
  static_assert(sizeof(Params) == 8 * sizeof(float),
                "Fishbowl::Params field set changed — update PRESET's "
                "float list to match");
  Params params = PRESET;

  // Precedes scale_mod, which binds &params.scale_factor at construction.
  ScaleModifier scale_mod{
      &params.scale_factor}; /**< Palette scale coordinate modifier. */
  CycleModifier cycle_mod{
      &cycle_phase}; /**< Palette cycle coordinate modifier. */
  DutyCycleModifier duty_mod{&params.duty_cycle};
  Gradient fire_palette;
  StaticPalette<Gradient,
                Coords<ScaleModifier, CycleModifier, DutyCycleModifier>,
                Colors<>, false>
      static_palette; /**< Bound palette sampled per fragment. */
  Animation::Motion<W, ORIENTATION_SUBSTEPS> *motion =
      nullptr; /**< Retained motion handle for live Cycle Dur updates. */
  int last_cycle_duration = -1; /**< Last applied cycle duration, in frames. */
  float cycle_phase = 0.0f;     /**< Current palette cycle phase. */
  float palette_speed = 0.0f;   /**< Fill-compensated palette phase speed. */
  Node *node = nullptr; /**< The single animated body, arena-allocated. */
  NoiseTransformer<1>
      noise_xform; /**< Warps the trail with noise each frame. */
};

#include "core/control/registry.h"
REGISTER_EFFECT(Fishbowl)
