/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H> class Raymarch : public Effect {
public:
  FLASHMEM Raymarch() : Effect(W, H) {}

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
    for (int i = 0; i < 20; ++i) {
      raw_quats[i] = make_rotation(Y_AXIS, Solids::Dodecahedron::vertices[i]);
    }

    timeline.add(0, Animation::RandomWalk<W>(camera, normal, noise));

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->drawFn(canvas, opacity);
                        },
                        -1));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:
  // --- Helper functions ---

  /// Metallic Blinn-Phong: half-Lambert diffuse, tight specular, Fresnel rim.
  /// All components tinted by surface color (metallic reflection).
  float shadeBlinnPhong(const Vector &normal_w, const Vector &light_dir,
                        const Vector &view_dir, const Vector &tangent) {
    // Diffuse: half-Lambert wrap using light direction
    float ndotl = normal_w.x * light_dir.x + normal_w.y * light_dir.y +
                  normal_w.z * light_dir.z;
    float half_lam = ndotl * 0.5f + 0.5f;
    float diffuse = half_lam * half_lam;

    // Specular: light tilted off-axis along tangent for visible highlight
    float lx = light_dir.x + tangent.x * 0.3f;
    float ly = light_dir.y + tangent.y * 0.3f;
    float lz = light_dir.z + tangent.z * 0.3f;
    float ll = sqrtf(lx * lx + ly * ly + lz * lz);
    if (ll > TOLERANCE) {
      lx /= ll;
      ly /= ll;
      lz /= ll;
    }
    // Half-vector between tilted light and view direction
    float hx = lx + view_dir.x, hy = ly + view_dir.y, hz = lz + view_dir.z;
    float hl = sqrtf(hx * hx + hy * hy + hz * hz);
    if (hl > TOLERANCE) {
      hx /= hl;
      hy /= hl;
      hz /= hl;
    }
    float ndoth =
        std::max(0.0f, normal_w.x * hx + normal_w.y * hy + normal_w.z * hz);
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

  void drawFn(Canvas &canvas, float opacity) {
    float t = static_cast<float>(timeline.t) / 60.0f;
    float anim_t = t * params.pulse_speed;

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

    // Tumble rotation around local X-axis (tangent)
    float spin_angle = anim_t * 1.5f;
    Quaternion spin_q = make_rotation(X_AXIS, spin_angle);

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
        float shade = shadeBlinnPhong(n_world, vertex, vertex, tangent);

        float ring_angle = (atan2f(loc.z, loc.x) + PI_F) / (2.0f * PI_F);
        float palette_t = fmodf(
            ring_angle + anim_t * 0.05f + static_cast<float>(vi) / 20.0f, 1.0f);
        Color4 c = palette.get(palette_t);
        frag.color = Color4(c.color * shade, opacity);
      };

      Scan::TransformedVolume vol(torus, vertex, world_q);
      Scan::Volume::draw<W, H>(pipeline, canvas, vertex, bounds_radius,
                               view_dir, vol, frag_fn, max_steps, aa_width);
    }
  }

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
  Orientation<W> camera;
  std::array<Quaternion, 20> raw_quats;
  Timeline<W> timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::COMPLEMENTARY,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            219};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Raymarch)
