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

    // Precompute frame quaternions — dodecahedron vertices are constexpr
    for (int i = 0; i < NUM_VERTS; ++i) {
      raw_quats[i] = make_rotation(Y_AXIS, Solids::Dodecahedron::vertices[i]);
    }

    // Bake the fixed-seed palette into a 256-entry LUT once. `palette` is
    // immutable (constructed with a fixed manual seed, never mutated), so the
    // per-fragment shader can do a cheap LUT lookup instead of the full
    // ~2000-cycle dual sRGB->OKLCH conversion on every ray-marched fragment.
    baked_palette.bake(persistent_arena, palette);

    timeline.add(0, Animation::RandomWalk<W>(camera, normal, noise));

    // Drive the tumble and palette-scroll phases from the live Pulse Speed
    // slider instead of reading the unbounded global frame clock: each is an
    // effect-owned accumulator wrapped to [0,1) every step, so the trig
    // argument never grows (stable precision over long runs) and the phase no
    // longer couples to the global counter at effect-switch time. The scales
    // reproduce the prior per-animation-second increments (1.5 rad of spin and
    // 0.05 cycle of palette scroll at 60 fps); spin_phase is scaled to radians
    // by *2pi where it is consumed.
    constexpr float kTwoPi = 2.0f * PI_F;
    timeline.add(0, Animation::Driver(spin_phase, &params.pulse_speed,
                                      1.5f / (60.0f * kTwoPi), true));
    timeline.add(0, Animation::Driver(palette_phase, &params.pulse_speed,
                                      0.05f / 60.0f, true));

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1));
  }

  /**
   * @brief Reports whether the engine should render its background behind this
   *        effect.
   * @return Always false; the ray march fills the frame itself.
   */
  bool show_bg() const override { return false; }

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
    // ray-march per dodecahedron vertex (NUM_VERTS = 20), each marching up to
    // max_steps (slider 4..30) shape.distance() evaluations at every covered
    // pixel. The Scan::Volume bounds_radius cull below trims pixels outside
    // each torus's screen disc, but the sphere covers a large fraction of the
    // frame, so this is the heaviest effect in the set alongside Voronoi.
    // Left unbudgeted (no per-frame cap against W*H); acceptable on the WASM
    // sim. If a device target ever runs this effect, lower max_steps or cap the
    // marched pixel count against the frame budget here.
    for (int vi = 0; vi < NUM_VERTS; ++vi) {
      Vector vertex = camera.orient(Solids::Dodecahedron::vertices[vi]);
      Vector view_dir(-vertex.x, -vertex.y, -vertex.z);

      // Compose: frame (local→world at vertex) * spin (rotate in local X)
      // Then apply camera orientation
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
  BakedPalette baked_palette; // LUT baked once from `palette` (immutable) in init()
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Raymarch)
