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
 * switch immediately at fixed intervals while mesh emission continues. The
 * shape never changes and never morphs.
 */
template <int W, int H> class MeshFeedback : public Effect {
public:
  using Style = Feedback::Style;

  static constexpr int PRESET_FRAMES = 241;

  // Gamut boundary bracket grid this effect buys from the persistent arena
  // (131,072 B, leaving ~124 KiB of the partition free). Resolution sets only
  // how wide the bracket starts, and the per-pixel refinement narrows it by a
  // factor of 32, so this trades arena bytes against refinement steps rather
  // than against accuracy: the next step up in either axis doubles the cost to
  // save well under one step.
  static constexpr int GAMUT_ANGLE_STEPS = 256;
  static constexpr int GAMUT_L_STEPS = 128;

  static constexpr float FADE_MIN = 0.0f, FADE_MAX = 0.99f;
  static constexpr float AMP_MIN = 0.0f, AMP_MAX = 30.0f;
  static constexpr float FREQ_MIN = 0.01f, FREQ_MAX = 1.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float SCALE_MIN = 0.1f, SCALE_MAX = 50.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 0.5f;

  /** @brief True iff every preset-driven field of @p s lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Style &s) {
    return s.fade >= FADE_MIN && s.fade <= FADE_MAX && s.amplitude >= AMP_MIN &&
           s.amplitude <= AMP_MAX && s.frequency >= FREQ_MIN &&
           s.frequency <= FREQ_MAX && s.speed >= SPEED_MIN &&
           s.speed <= SPEED_MAX && s.scale >= SCALE_MIN &&
           s.scale <= SCALE_MAX && s.hue_shift >= HUE_SHIFT_MIN &&
           s.hue_shift <= HUE_SHIFT_MAX;
  }
  static constexpr std::array<PresetEntry<Style>, 12> PRESETS = {{
      {Style::ArcingLightning()},
      {Style::SlowFire()},
      {Style::EnergeticFire()},
      {Style::Smoke()},
      {Style::SlowDust()},
      {Style::WavyTrails()},
      {Style::MeltingHi()},
      {Style::MeltingLo()},
      {Style::Miasma()},
      {Style::LooseWormhole()},
      {Style::TightWormhole()},
      {Style::WigglingWormhole()},
  }};

  /** @brief True iff every entry of PRESETS satisfies preset_in_ranges(). */
  static constexpr bool all_presets_in_ranges() {
    for (const auto &e : PRESETS)
      if (!preset_in_ranges(e.params))
        return false;
    return true;
  }
  static_assert(
      all_presets_in_ranges(),
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
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(filters)::any_crosses_segments}),
        noise_params(), orientation(), timeline(), palette(Palettes::PEACH_POP),
        filters(Filter::World::Orient<W>(orientation),
                Filter::Screen::AntiAlias<W, H>(),
                Filter::Pixel::Feedback<W, H>(style)) {}

  /**
   * @brief One-time effect setup.
   * @details Binds shared noise into presets, builds the icosahedron, samples
   * the mesh shade, registers tunable params, and schedules the
   * noise/walk/preset timers.
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
      PolyMesh poly =
          generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
            return Solids::finalize_solid(Solids::Platonic::icosahedron(a, b),
                                          target);
          });
      MeshOps::compile(poly, mesh_, persistent_arena, scratch_arena_a);
    }
    mesh_shade_ = palette.get(0.0f);

    register_animated_param("Fade", &style.fade, FADE_MIN, FADE_MAX);
    register_animated_param("Distort Amp", &style.amplitude, AMP_MIN, AMP_MAX);
    register_animated_param("Distort Freq", &style.frequency, FREQ_MIN,
                            FREQ_MAX);
    register_animated_param("Distort Speed", &style.speed, SPEED_MIN,
                            SPEED_MAX);
    register_animated_param("Noise Scale", &style.scale, SCALE_MIN, SCALE_MAX);
    register_animated_param("Hue Shift", &style.hue_shift, HUE_SHIFT_MIN,
                            HUE_SHIFT_MAX);
    register_param("Feedback", &feedback_enabled);

    filters.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
        persistent_arena);
    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);

    timeline.add(0, Animation::Noise(noise_params));
    timeline.add(
        0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise_params.noise));
  }

  /**
   * @brief Renders one frame.
   * @details Advances the preset selection, applies params, runs the feedback
   * decay flush, draws the mesh, then advances the timeline. The preset switch
   * leads apply_params() so the noise scalars and the fade/hue the flush reads
   * come from the same preset.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    advance_transition();

    {
      HS_PROFILE(mf_apply_params);
      apply_params();
    }

    {
      // Feedback-buffer warp/tap + decay flush.
      HS_PROFILE(mf_feedback_flush);
      filters.flush(canvas, 1.0f);
    }

    {
      HS_PROFILE(mf_mesh_draw);
      const Color4 shade = mesh_shade_;
      Plot::Mesh::draw<W, H>(
          filters, canvas, mesh_,
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
   * @brief Advances to the next preset every PRESET_FRAMES.
   * @details Frozen while animations are paused.
   */
  void advance_transition() {
    if (animations_paused())
      return;
    if (++transition_frames_ < PRESET_FRAMES)
      return;
    transition_frames_ = 0;
    presets.next();
    presets.apply(style);
  }

  Style style;

  Presets<Style, 12> presets{PRESETS};
  bool feedback_enabled = true;
  int transition_frames_ = 0;
  NoiseParams noise_params;

  Orientation<> orientation;
  Timeline timeline;
  ProceduralPalette palette;

  // The mesh draws one fixed color, so the palette is sampled once in init().
  Color4 mesh_shade_;

  // The single, fixed solid; built once in init() and never recompiled.
  MeshState mesh_;

  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>,
           Filter::Pixel::Feedback<W, H>>
      filters;

  // Covers the two compile-time-sized persistent tenants: the Feedback
  // warp-field cache and the gamut boundary bracket table. init() also buys the
  // icosahedron PolyMesh and its compiled MeshState from the persistent arena,
  // but those are sized at runtime and ride in the headroom left here.
  static_assert(Filter::Pixel::Feedback<W, H>::STORAGE_BYTES +
                        gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) <=
                    DEVICE_PERSISTENT_BUDGET,
                "MeshFeedback warp cache plus gamut bracket table exceed the "
                "default persistent partition; retune the feedback downsample, "
                "coarsen the gamut grid, or carve arenas");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MeshFeedback)
