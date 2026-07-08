/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Sphere-traces a morphing superquadric at each cube vertex: the
 *        exponent sweeps octahedron -> sphere -> cube and back, phase-offset
 *        by latitude so the morph travels pole-to-pole as a wave. Shaded with
 *        the metallic headlight model and a baked OKLCH palette.
 * @tparam W Effect render width in pixels.
 * @tparam H Effect render height in pixels.
 */
template <int W, int H> class MorphBlob : public Effect {
public:
  /**
   * @brief Constructs the effect at the templated render dimensions.
   */
  FLASHMEM MorphBlob() : Effect(W, H) {}

  /**
   * @brief Registers tunable params, precomputes per-vertex quaternions, bakes
   *        the palette LUT, and installs the camera walk and per-frame draw on
   *        the timeline.
   */
  void init() override {
    registerParam("Morph Speed", &params.morph_speed, 0.0f, 10.0f);
    // Upper bound keeps the cube-end circumradius (s * 3^(1/2 - 1/8) at max
    // Boxiness) below half the 70.5° vertex spacing, so no slider setting
    // overlaps a neighbour.
    registerParam("Core Size", &params.core_size, 0.1f, 0.38f);
    registerParam("Spike", &params.spike, 1.0f, 2.0f);
    registerParam("Boxiness", &params.boxiness, 2.0f, 8.0f);
    registerParam("Wave", &params.wave, 0.0f, 1.0f);
    registerParam("Max Steps", &params.max_steps, 4.0f, 30.0f);
    registerParam("Diffuse", &params.diffuse, 0.0f, 1.0f);
    registerParam("Specular", &params.specular, 0.0f, 1.5f);
    registerParam("Fresnel", &params.fresnel, 0.0f, 1.0f);
    registerParam("AA Width", &params.aa_mult, 0.1f, 1.5f);

    for (int i = 0; i < NUM_VERTS; ++i)
      raw_quats[i] = make_rotation(Y_AXIS, Solids::Cube::vertices[i]);

    baked_palette.bake(persistent_arena, palette);

    timeline.add(0, Animation::RandomWalk<W>(camera, normal, noise));

    // morph_phase / spin_phase / palette_phase are effect-owned accumulators
    // wrapped to [0,1) each step, so the trig argument never grows.
    constexpr float TWO_PI_F = 2.0f * PI_F;
    timeline.add(0, Animation::Driver(morph_phase, &params.morph_speed,
                                      1.0f / (60.0f * TWO_PI_F), true));
    timeline.add(0, Animation::Driver(spin_phase, &params.morph_speed,
                                      0.6f / (60.0f * TWO_PI_F), true));
    timeline.add(0, Animation::Driver(palette_phase, &params.morph_speed,
                                      0.05f / 60.0f, true));

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1));
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Steps the timeline against a fresh canvas, driving the per-frame ray
   *        march.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  static constexpr int NUM_VERTS = Solids::Cube::NUM_VERTS;

  /**
   * @brief Ray-marches and shades the morphing superquadric at every cube
   *        vertex for the current frame.
   * @param canvas Render target for this frame's fragments.
   * @param opacity Sprite fade alpha in [0, 1], written into each fragment's
   *                color.
   */
  void drawFn(Canvas &canvas, float opacity) {
    constexpr float TWO_PI_F = 2.0f * PI_F;

    float s = params.core_size;
    // Bounds must cover the widest shape of the sweep: the cube end's
    // circumradius at the Boxiness exponent.
    float bounds_radius =
        s * powf(3.0f, std::max(0.0f, 0.5f - 1.0f / params.boxiness));
    float aa_width = s * 0.15f * params.aa_mult;
    int max_steps = static_cast<int>(params.max_steps + 0.5f);
    float log_ratio = logf(params.boxiness / params.spike);

    // spin_phase rides in [0,1); scale to radians for make_rotation.
    Quaternion spin_q = make_rotation(X_AXIS, spin_phase * TWO_PI_F);

    for (int vi = 0; vi < NUM_VERTS; ++vi) {
      Vector home = Solids::Cube::vertices[vi];
      // Latitude phase offset: the morph sweeps pole-to-pole as a wave.
      float phase_i =
          morph_phase + params.wave * (0.5f + 0.5f * home.y);
      float t = 0.5f + 0.5f * fast_sinf(TWO_PI_F * phase_i);
      // Exponential exponent lerp: Spike at t=0, Boxiness at t=1, crossing
      // the sphere (p=2) between the two pointy extremes.
      SDF::Superquadric blob(s, params.spike * expf(t * log_ratio));

      Vector center = camera.orient(home);
      Vector ray_dir(-center.x, -center.y, -center.z);

      Quaternion world_q = camera.get() * raw_quats[vi] * spin_q;
      Vector tangent = rotate(Vector(1, 0, 0), world_q);

      auto frag_fn = [&](const Vector &loc, Fragment &frag) {
        Vector n_world = rotate(blob.normal(loc), world_q);
        // Headlight model: light coincides with the viewer, so the view vector
        // `center` serves as both light_dir and view_dir.
        float shade = shade_blinn_phong(n_world, center, center, tangent,
                                        params.diffuse, params.specular,
                                        params.fresnel);

        float axis_angle = (fast_atan2(loc.z, loc.x) + PI_F) / TWO_PI_F;
        float palette_t = fmodf(
            axis_angle + palette_phase + static_cast<float>(vi) / NUM_VERTS,
            1.0f);
        Color4 c = baked_palette.get(palette_t);
        frag.color = Color4(c.color * shade, opacity);
      };

      Scan::TransformedVolume vol(blob, center, world_q);
      Scan::Volume::draw<W, H>(pipeline, canvas, center, bounds_radius,
                               ray_dir, vol, frag_fn, max_steps, aa_width);
    }
  }

  /**
   * @brief Tunable shader and animation parameters exposed via registerParam.
   */
  struct Params {
    float morph_speed = 5.0f;
    float core_size = 0.3f;
    float spike = 1.0f;
    float boxiness = 6.0f;
    float wave = 0.35f;
    float max_steps = 18.0f;
    float diffuse = 0.4f;
    float specular = 1.2f;
    float fresnel = 0.25f;
    float aa_mult = 0.5f;
  } params;

  FastNoiseLite noise;
  Vector normal = Y_AXIS;
  Orientation<> camera;
  float morph_phase = 0.0f;   // exponent sweep phase, [0,1) cycles
  float spin_phase = 0.0f;    // blob tumble phase, [0,1) -> [0,2pi) radians
  float palette_phase = 0.0f; // baked-palette scroll offset, [0,1) cycles
  std::array<Quaternion, NUM_VERTS> raw_quats;
  Timeline timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            83};
  BakedPalette baked_palette;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MorphBlob)
