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

    // Precompute tangent frames — dodecahedron vertices are constexpr
    for (int i = 0; i < 20; ++i)
      raw_frames[i] = buildTangentFrame(Solids::Dodecahedron::vertices[i]);

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
  // --- Geometry types ---

  struct TangentFrame {
    Vector vertex;   // outward normal = torus axis
    Vector tangent;  // first tangent (spin axis)
    Vector tangent2; // second tangent (vertex × tangent)
  };

  // --- Helper functions ---

  /// Build an orthonormal tangent frame at a unit-sphere vertex.
  static TangentFrame buildTangentFrame(const Vector &vertex) {
    Vector ref = (fabsf(vertex.y) < 0.9f) ? Y_AXIS : X_AXIS;
    float d = ref.x * vertex.x + ref.y * vertex.y + ref.z * vertex.z;

    Vector tangent(ref.x - d * vertex.x, ref.y - d * vertex.y,
                   ref.z - d * vertex.z);
    float tlen = tangent.length();
    if (tlen > TOLERANCE)
      tangent = tangent * (1.0f / tlen);

    Vector tangent2(vertex.y * tangent.z - vertex.z * tangent.y,
                    vertex.z * tangent.x - vertex.x * tangent.z,
                    vertex.x * tangent.y - vertex.y * tangent.x);

    return {vertex, tangent, tangent2};
  }

  /// Transform a world-space point into the spun local frame.
  /// Orientation (spin) is applied after the tangent frame projection.
  static Vector worldToSpunLocal(const Vector &world_p, const Vector &center,
                                 const TangentFrame &frame, float spin_cos,
                                 float spin_sin) {
    float ox = world_p.x - center.x;
    float oy = world_p.y - center.y;
    float oz = world_p.z - center.z;

    float lx =
        ox * frame.tangent.x + oy * frame.tangent.y + oz * frame.tangent.z;
    float ly = ox * frame.vertex.x + oy * frame.vertex.y + oz * frame.vertex.z;
    float lz =
        ox * frame.tangent2.x + oy * frame.tangent2.y + oz * frame.tangent2.z;

    return Vector(lx, ly * spin_cos - lz * spin_sin,
                  ly * spin_sin + lz * spin_cos);
  }

  /// Transform a local-space normal back to world space (including un-spin).
  static Vector localNormalToWorld(const Vector &n_local,
                                   const TangentFrame &frame, float spin_cos,
                                   float spin_sin) {
    float ny_us = n_local.y * spin_cos + n_local.z * spin_sin;
    float nz_us = -n_local.y * spin_sin + n_local.z * spin_cos;

    return Vector(n_local.x * frame.tangent.x + ny_us * frame.vertex.x +
                      nz_us * frame.tangent2.x,
                  n_local.x * frame.tangent.y + ny_us * frame.vertex.y +
                      nz_us * frame.tangent2.y,
                  n_local.x * frame.tangent.z + ny_us * frame.vertex.z +
                      nz_us * frame.tangent2.z);
  }

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

  // Dodecahedron vertices from Solids::Dodecahedron (constexpr)
  static constexpr int NUM_VERTS = Solids::Dodecahedron::NUM_VERTS;

  // --- Main draw ---

  void drawFn(Canvas &canvas, float opacity) {
    float t = static_cast<float>(timeline.t) / 60.0f;
    float anim_t = t * params.pulse_speed;

    SDF::TwistedTorus torus{params.core_size * 0.45f, params.core_size * 0.14f,
                            static_cast<int>(params.twist),
                            params.core_size * 0.14f * 2.5f};
    float bounds_radius = sqrtf((torus.R + torus.r) * (torus.R + torus.r) +
                                torus.amplitude * torus.amplitude) +
                          torus.r;
    float aa_width = torus.r * params.aa_mult;
    int max_steps = static_cast<int>(params.max_steps + 0.5f);

    float spin_angle = anim_t * 1.5f;
    float spin_cos = cosf(spin_angle);
    float spin_sin = sinf(spin_angle);

    for (int vi = 0; vi < NUM_VERTS; ++vi) {
      const Vector &raw_vertex = Solids::Dodecahedron::vertices[vi];
      const TangentFrame &raw_frame = raw_frames[vi];
      Vector vertex = camera.orient(raw_vertex);
      TangentFrame frame = {vertex, camera.orient(raw_frame.tangent),
                            camera.orient(raw_frame.tangent2)};
      Vector view_dir(-vertex.x, -vertex.y, -vertex.z);

      // SDF is pure geometry — spin is a pre-transform applied here
      struct SpunTorus {
        const SDF::TwistedTorus &torus;
        const Vector &center;
        const TangentFrame &frame;
        float spin_cos, spin_sin;

        float distance(const Vector &p) const {
          Vector loc = worldToSpunLocal(p, center, frame, spin_cos, spin_sin);
          // Cheap bound (no trig) for the approach phase
          float bd = torus.bounding_distance(loc);
          if (bd > torus.r)
            return bd;
          // Close — full SDF with twist
          return torus.distance(loc);
        }
      } spun_torus{torus, vertex, frame, spin_cos, spin_sin};

      auto frag_fn = [&](const Vector &hit, Fragment &frag) {
        Vector loc = worldToSpunLocal(hit, vertex, frame, spin_cos, spin_sin);

        // Normal directly at march point — skip surface projection since
        // the hit is within aa_width of the surface (normal is stable).
        Vector n_local = torus.normal(loc);
        Vector n_world = localNormalToWorld(n_local, frame, spin_cos, spin_sin);
        float shade =
            shadeBlinnPhong(n_world, frame.vertex, frame.vertex, frame.tangent);

        // Ring angle from un-spun local frame (stable across spin rotation).
        // The spin only rotates local YZ; lx and lz (tangent projections)
        // are invariant, so atan2(lz, lx) gives a stable ring angle.
        float ox = hit.x - vertex.x, oy = hit.y - vertex.y,
              oz = hit.z - vertex.z;
        float lx =
            ox * frame.tangent.x + oy * frame.tangent.y + oz * frame.tangent.z;
        float lz = ox * frame.tangent2.x + oy * frame.tangent2.y +
                   oz * frame.tangent2.z;
        float ring_angle = (atan2f(lz, lx) + PI_F) / (2.0f * PI_F);

        float palette_t = fmodf(
            ring_angle + anim_t * 0.05f + static_cast<float>(vi) / 20.0f, 1.0f);
        Color4 c = palette.get(palette_t);
        frag.color = Color4(c.color * shade, opacity);
      };

      Scan::Volume::draw<W, H>(pipeline, canvas, vertex, bounds_radius,
                               view_dir, spun_torus, frag_fn, max_steps,
                               aa_width);
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
  std::array<TangentFrame, 20> raw_frames;
  Timeline<W> timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::COMPLEMENTARY,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            219};
};

#include "../effect_registry.h"
REGISTER_EFFECT(Raymarch)
