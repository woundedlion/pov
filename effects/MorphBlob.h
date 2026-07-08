/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Sphere-traces a morphing superquadric at each vertex of a rhombic
 *        triacontahedron: the exponent sweeps octahedron -> sphere -> cube and
 *        back, phase-offset by latitude so the morph travels pole-to-pole as a
 *        wave. Each blob is auto-sized to its own nearest-neighbour gap (scaled
 *        by the live Fill param), so they pack without overlap. Shaded with the
 *        metallic headlight model and a baked OKLCH palette.
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
   * @brief Registers tunable params, builds the rhombic-triacontahedron vertex
   *        set, bakes the palette LUT, and installs the camera walk and
   *        per-frame draw on the timeline.
   */
  void init() override {
    registerParam("Morph Speed", &params.morph_speed, 0.0f, 10.0f);
    // Fraction of the half nearest-neighbour gap the blob's widest sweep extent
    // (the cube end's corner circumradius) reaches: < 1 leaves a gap, 1 makes
    // neighbours touch, > 1 overlaps them deliberately.
    registerParam("Fill", &params.fill, 0.3f, 1.3f);
    registerParam("Spike", &params.spike, 1.0f, 2.0f);
    registerParam("Boxiness", &params.boxiness, 2.0f, 8.0f);
    registerParam("Wave", &params.wave, 0.0f, 1.0f);
    registerParam("Max Steps", &params.max_steps, 4.0f, 30.0f);
    registerParam("Diffuse", &params.diffuse, 0.0f, 1.0f);
    registerParam("Specular", &params.specular, 0.0f, 1.5f);
    registerParam("Fresnel", &params.fresnel, 0.0f, 1.0f);
    registerParam("AA Width", &params.aa_mult, 0.1f, 1.5f);

    build_points();

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
  /** Vertex-array capacity; the rhombic triacontahedron has 32. */
  static constexpr int MAX_POINTS = 32;

  /**
   * @brief Builds the rhombic-triacontahedron vertex directions and per-vertex
   *        orientation quaternions and nearest-neighbour gaps.
   * @details Generates the solid from the registry into the scratch arenas
   *          (reclaimed on scope exit) via the shared HS_COLD builder, so drawFn
   *          can size every blob to its own local gap.
   */
  void build_points() {
    ScratchScope a_guard(scratch_arena_a);
    ScratchScope b_guard(scratch_arena_b);
    active_count = Solids::build_vertex_directions(
        scratch_arena_a, scratch_arena_b, "rhombicTriacontahedron", MAX_POINTS,
        points.data(), raw_quats.data(), nn_angle.data());
  }

  /**
   * @brief Ray-marches and shades the morphing superquadric at every vertex for
   *        the current frame.
   * @param canvas Render target for this frame's fragments.
   * @param opacity Sprite fade alpha in [0, 1], written into each fragment's
   *                color.
   */
  void drawFn(Canvas &canvas, float opacity) {
    constexpr float TWO_PI_F = 2.0f * PI_F;

    // Widest sweep extent per unit scale: the cube end's corner circumradius.
    float diag_k = powf(3.0f, std::max(0.0f, 0.5f - 1.0f / params.boxiness));
    int max_steps = static_cast<int>(params.max_steps + 0.5f);
    float log_ratio = logf(params.boxiness / params.spike);

    // spin_phase rides in [0,1); scale to radians for make_rotation.
    Quaternion spin_q = make_rotation(X_AXIS, spin_phase * TWO_PI_F);

    for (int i = 0; i < active_count; ++i) {
      // Latitude phase offset: the morph sweeps pole-to-pole as a wave.
      float phase_i = morph_phase + params.wave * (0.5f + 0.5f * points[i].y);
      float t = 0.5f + 0.5f * fast_sinf(TWO_PI_F * phase_i);
      // Exponential exponent lerp: Spike at t=0, Boxiness at t=1, crossing
      // the sphere (p=2) between the two pointy extremes.
      float p = params.spike * expf(t * log_ratio);

      // Per-vertex auto-size: fit the cube end's circumradius to `fill` of this
      // vertex's half nearest-neighbour gap, so open regions get large blobs
      // and tight ones stay small; at fill 1 mutual neighbours just touch.
      float outer_r = sinf(0.5f * nn_angle[i] * params.fill);
      float s = outer_r / diag_k;
      SDF::Superquadric blob(s, p);
      float aa_width = s * 0.15f * params.aa_mult;
      float bounds_radius = outer_r + aa_width;

      Vector center = camera.orient(points[i]);
      Vector ray_dir(-center.x, -center.y, -center.z);

      Quaternion world_q = camera.get() * raw_quats[i] * spin_q;
      Vector tangent = rotate(Vector(1, 0, 0), world_q);

      auto frag_fn = [&](const Vector &loc, Fragment &frag) {
        Vector n_world = rotate(blob.normal(loc), world_q);
        // Headlight model: light coincides with the viewer, so the view vector
        // `center` serves as both light_dir and view_dir.
        float shade = shade_blinn_phong(n_world, center, center, tangent,
                                        params.diffuse, params.specular,
                                        params.fresnel);

        float axis_angle = (fast_atan2(loc.z, loc.x) + PI_F) / TWO_PI_F;
        float palette_t = fmodf(axis_angle + palette_phase +
                                    static_cast<float>(i) / active_count,
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
    float fill = 0.8f;
    float spike = 1.0f;
    float boxiness = 6.0f;
    float wave = 0.5f;
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
  int active_count = 0;       // vertices built (rhombic triacontahedron = 32)
  std::array<Vector, MAX_POINTS> points;
  std::array<Quaternion, MAX_POINTS> raw_quats;
  std::array<float, MAX_POINTS> nn_angle; // per-vertex nearest-neighbour gap (rad)
  Timeline timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            83};
  BakedPalette baked_palette;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MorphBlob)
