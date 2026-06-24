/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

/**
 * @brief Ray-marches a twisted torus SDF once per dodecahedron vertex, shading
 *        each with a metallic headlight model and a baked OKLCH palette.
 * @tparam W Effect render width in pixels.
 * @tparam H Effect render height in pixels.
 */
template <int W, int H> class Raymarch : public Effect {
public:
  /**
   * @brief Constructs the effect at the templated render dimensions.
   */
  FLASHMEM Raymarch() : Effect(W, H) {}

  /**
   * @brief Registers tunable params, precomputes per-vertex quaternions, bakes
   *        the palette LUT, and installs the camera walk and per-frame draw on
   *        the timeline.
   */
  void init() override {
    registerParam("Pulse Speed", &params.pulse_speed, 0.0f, 10.0f);
    registerParam("Core Size", &params.core_size, 0.1f, 0.8f);
    registerParam("Max Steps", &params.max_steps, 4.0f, 30.0f);
    registerParam("Diffuse", &params.diffuse, 0.0f, 1.0f);
    registerParam("Specular", &params.specular, 0.0f, 1.5f);
    registerParam("Fresnel", &params.fresnel, 0.0f, 1.0f);
    registerParam("Twist", &params.twist, 0.0f, 8.0f);
    registerParam("AA Width", &params.aa_mult, 0.1f, 1.5f);

    for (int i = 0; i < NUM_VERTS; ++i) {
      raw_quats[i] = make_rotation(Y_AXIS, Solids::Dodecahedron::vertices[i]);
    }

    // Bake the immutable fixed-seed palette into a 256-entry LUT so the
    // per-fragment shader does a cheap lookup instead of the full ~2000-cycle
    // dual sRGB->OKLCH conversion per fragment.
    baked_palette.bake(persistent_arena, palette);

    timeline.add(0, Animation::RandomWalk<W>(camera, normal, noise));

    // Drive the tumble and palette-scroll phases from the live Pulse Speed
    // slider rather than the global frame clock: each is an effect-owned
    // accumulator wrapped to [0,1) every step, so the trig argument never grows.
    // The scales give 1.5 rad of spin and 0.05 cycle of palette scroll per
    // animation-second at 60 fps; spin_phase is scaled to radians by *2pi where
    // consumed. Both gated on the pause flag.
    constexpr float kTwoPi = 2.0f * PI_F;
    timeline.add(0, Animation::Driver(spin_phase, &params.pulse_speed,
                                      1.5f / (60.0f * kTwoPi), true,
                                      &anims_paused_));
    timeline.add(0, Animation::Driver(palette_phase, &params.pulse_speed,
                                      0.05f / 60.0f, true, &anims_paused_));

    // The render Sprite needs no pause gate: it is perpetual and its image is
    // driven entirely by the two phase accumulators above, frozen while paused.
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1));
  }

  /**
   * @brief POV column-strobe flag — see Effect::strobe_columns.
   * @return true; the strip blanks to black immediately after each column is
   *         shown, so every column reads as a sharp slice with dark gaps
   *         between columns rather than persisting across the sweep.
   */
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
  // --- Helper functions ---

  /**
   * @brief Metallic Blinn-Phong shade factor: half-Lambert diffuse, tight
   *        specular, Fresnel rim.
   * @param normal_w Surface normal in world space (unit length).
   * @param light_dir Direction toward the light in world space (unit length).
   * @param view_dir Direction toward the viewer in world space (unit length).
   * @param tangent Surface tangent used to tilt the specular highlight off-axis.
   * @return Scalar shade factor (ambient base plus weighted diffuse/specular/
   *         Fresnel terms); the caller multiplies it by the surface color for
   *         the metallic look.
   * @details The sole caller uses a headlight model — the light and the viewer
   *          coincide, so it passes the same view vector for both light_dir and
   *          view_dir. The two parameters are kept distinct to preserve the
   *          standard Blinn-Phong signature should the light ever decouple from
   *          the camera.
   */
  float shadeBlinnPhong(const Vector &normal_w, const Vector &light_dir,
                        const Vector &view_dir, const Vector &tangent) {
    // Diffuse: half-Lambert wrap using light direction
    float ndotl = dot(normal_w, light_dir);
    float half_lam = ndotl * 0.5f + 0.5f;
    float diffuse = half_lam * half_lam;

    // Specular: light tilted off-axis along tangent for visible highlight.
    // Normalize only above TOLERANCE; a near-zero tilt stays as-is.
    Vector light = light_dir + tangent * 0.3f;
    float ll = light.length();
    if (ll > TOLERANCE)
      light /= ll;
    // Half-vector between tilted light and view direction
    Vector half = light + view_dir;
    float hl = half.length();
    if (hl > TOLERANCE)
      half /= hl;
    float ndoth = std::max(0.0f, dot(normal_w, half));
    // ndoth^32 via repeated squaring
    float spec = ndoth * ndoth; // ^2
    spec *= spec;               // ^4
    spec *= spec;               // ^8
    spec *= spec;               // ^16
    spec *= spec;               // ^32

    float fresnel = 1.0f - hs::clamp(ndotl, 0.0f, 1.0f);
    fresnel = fresnel * fresnel * fresnel;

    return 0.05f + diffuse * params.diffuse + spec * params.specular +
           fresnel * params.fresnel;
  }

  static constexpr int NUM_VERTS = Solids::Dodecahedron::NUM_VERTS;

  /**
   * @brief Ray-marches and shades the twisted torus at every dodecahedron
   *        vertex for the current frame.
   * @param canvas Render target for this frame's fragments.
   * @param opacity Sprite fade alpha in [0, 1], written into each fragment's
   *                color.
   */
  void drawFn(Canvas &canvas, float opacity) {
    constexpr float kTwoPi = 2.0f * PI_F;

    float major_r = params.core_size * 0.45f;
    float minor_r = params.core_size * 0.14f;
    int twist_n = static_cast<int>(params.twist);
    float twist_amp = minor_r * 2.5f;

    SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
        {major_r, minor_r}, {twist_n, twist_amp, major_r}};

    float bounds_radius = sqrtf((major_r + minor_r) * (major_r + minor_r) +
                                twist_amp * twist_amp) +
                          minor_r;
    float aa_width = minor_r * params.aa_mult;
    int max_steps = static_cast<int>(params.max_steps + 0.5f);

    // Tumble rotation around local X-axis (tangent). spin_phase rides in [0,1);
    // scale to radians for make_rotation (2pi-periodic, so the wrap is exact).
    float spin_angle = spin_phase * kTwoPi;
    Quaternion spin_q = make_rotation(X_AXIS, spin_angle);

    // Per-frame cost is O(NUM_VERTS * max_steps * W * H): one volumetric
    // ray-march per dodecahedron vertex, each up to max_steps distance() evals
    // per covered pixel. One of the heaviest effects (alongside Voronoi);
    // unbudgeted, acceptable on the WASM sim. For a device target, lower
    // max_steps or cap the marched pixel count against the frame budget here.
    for (int vi = 0; vi < NUM_VERTS; ++vi) {
      Vector vertex = camera.orient(Solids::Dodecahedron::vertices[vi]);
      Vector view_dir(-vertex.x, -vertex.y, -vertex.z);

      Quaternion world_q = camera.get() * raw_quats[vi] * spin_q;
      Vector tangent = rotate(Vector(1, 0, 0), world_q);

      auto frag_fn = [&](const Vector &loc, Fragment &frag) {
        Vector n_local = torus.normal(loc);
        Vector n_world = rotate(n_local, world_q);
        // Headlight model: light coincides with the viewer, so the view vector
        // `vertex` serves as both light_dir and view_dir.
        float shade = shadeBlinnPhong(n_world, vertex, vertex, tangent);

        float ring_angle = (atan2f(loc.z, loc.x) + PI_F) / (2.0f * PI_F);
        float palette_t = fmodf(
            ring_angle + palette_phase + static_cast<float>(vi) / NUM_VERTS,
            1.0f);
        Color4 c = baked_palette.get(palette_t);
        frag.color = Color4(c.color * shade, opacity);
      };

      Scan::TransformedVolume vol(torus, vertex, world_q);
      Scan::Volume::draw<W, H>(pipeline, canvas, vertex, bounds_radius,
                               view_dir, vol, frag_fn, max_steps, aa_width);
    }
  }

  /**
   * @brief Tunable shader and animation parameters exposed via registerParam.
   */
  struct Params {
    float pulse_speed = 5.0f;
    float core_size = 0.4f;
    float max_steps = 18.0f;
    float diffuse = 0.4f;
    float specular = 1.2f;
    float fresnel = 0.2f;
    float twist = 2.0f;
    float aa_mult = 0.5f;
  } params;

  FastNoiseLite noise;
  Vector normal = Y_AXIS;
  Orientation<> camera;
  float spin_phase = 0.0f;    // torus tumble phase, [0,1) -> [0,2pi) radians
  float palette_phase = 0.0f; // baked-palette scroll offset, [0,1) cycles
  std::array<Quaternion, NUM_VERTS> raw_quats;
  Timeline timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::COMPLEMENTARY,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            219};
  BakedPalette baked_palette;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Raymarch)
