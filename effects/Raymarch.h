/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Ray-marches a twisted torus SDF at each vertex of a disdyakis
 *        dodecahedron, shading each with a metallic headlight model and a baked
 *        OKLCH palette. Each torus is auto-sized to its own nearest-neighbour
 *        gap (scaled by the live Fill param), so they pack without overlap.
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
   * @brief Registers tunable params, builds the disdyakis-dodecahedron vertex
   *        set, bakes the palette LUT, and installs the camera walk and per-frame
   *        draw on the timeline.
   */
  void init() override {
    registerParam("Pulse Speed", &params.pulse_speed, 0.0f, 10.0f);
    // Fraction of the half nearest-neighbour gap the ring's outer edge reaches:
    // < 1 leaves a gap, 1 makes neighbours touch, > 1 overlaps them deliberately.
    registerParam("Fill", &params.fill, 0.3f, 1.3f);
    registerParam("Max Steps", &params.max_steps, 4.0f, 30.0f);
    registerParam("Diffuse", &params.diffuse, 0.0f, 1.0f);
    registerParam("Specular", &params.specular, 0.0f, 1.5f);
    registerParam("Fresnel", &params.fresnel, 0.0f, 1.0f);
    registerParam("Twist", &params.twist, 0.0f, 8.0f);
    registerParam("AA Width", &params.aa_mult, 0.1f, 1.5f);

    build_points();

    baked_palette.bake(persistent_arena, palette);

    timeline.add(0, Animation::RandomWalk<W>(camera, normal, noise));

    // spin_phase / palette_phase are effect-owned accumulators wrapped to [0,1)
    // each step, so the trig argument never grows. spin_phase is scaled to
    // radians by *2pi where consumed.
    constexpr float TWO_PI_F = 2.0f * PI_F;
    timeline.add(0, Animation::Driver(spin_phase, &params.pulse_speed,
                                      1.5f / (60.0f * TWO_PI_F), true));
    timeline.add(0, Animation::Driver(palette_phase, &params.pulse_speed,
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
  /** Vertex-array capacity; the disdyakis dodecahedron has 26. */
  static constexpr int MAX_POINTS = 32;

  /**
   * @brief Builds the disdyakis-dodecahedron vertex directions and per-vertex
   *        orientation quaternions and nearest-neighbour gaps.
   * @details Generates the solid from the registry into the scratch arenas
   *          (reclaimed on scope exit) via the shared HS_COLD builder, so drawFn
   *          can size every torus to its own local gap.
   */
  void build_points() {
    ScratchScope a_guard(scratch_arena_a);
    ScratchScope b_guard(scratch_arena_b);
    active_count = Solids::build_vertex_directions(
        scratch_arena_a, scratch_arena_b, "disdyakisDodecahedron", MAX_POINTS,
        points.data(), raw_quats.data(), nn_angle.data());
  }

  /**
   * @brief Ray-marches and shades the twisted torus at every vertex for the
   *        current frame.
   * @param canvas Render target for this frame's fragments.
   * @param opacity Sprite fade alpha in [0, 1], written into each fragment's
   *                color.
   */
  void drawFn(Canvas &canvas, float opacity) {
    constexpr float TWO_PI_F = 2.0f * PI_F;

    // Torus proportions at scale 1: VIS_K is the visible outer ring radius,
    // UNIT_BOUNDS the bounding-sphere radius (bigger, may overlap a neighbour —
    // a few wasted ray steps, no visual overlap).
    constexpr float MAJOR_K = 0.45f, MINOR_K = 0.14f, TWIST_K = 0.35f;
    constexpr float VIS_K = MAJOR_K + MINOR_K; // outer ring radius at scale 1
    constexpr float UNIT_BOUNDS = 0.826003f;   // √(VIS_K²+TWIST_K²)+MINOR_K
    int twist_n = static_cast<int>(params.twist);
    int max_steps = static_cast<int>(params.max_steps + 0.5f);

    // spin_phase rides in [0,1); scale to radians for make_rotation.
    float spin_angle = spin_phase * TWO_PI_F;
    Quaternion spin_q = make_rotation(X_AXIS, spin_angle);

    for (int i = 0; i < active_count; ++i) {
      // Per-vertex auto-size: fit the ring's outer edge to `fill` of this
      // vertex's half nearest-neighbour gap, so open regions get large tori and
      // tight ones stay small; at fill 1 mutual neighbours just touch.
      float outer_r = sinf(0.5f * nn_angle[i] * params.fill);
      float scale = outer_r / VIS_K;
      float bounds_radius = scale * UNIT_BOUNDS;
      float major_r = scale * MAJOR_K;
      float minor_r = scale * MINOR_K;
      float twist_amp = scale * TWIST_K;
      SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
          {major_r, minor_r}, {twist_n, twist_amp, major_r}};
      float aa_width = minor_r * params.aa_mult;

      Vector center = camera.orient(points[i]);
      Vector ray_dir(-center.x, -center.y, -center.z);

      Quaternion world_q = camera.get() * raw_quats[i] * spin_q;
      Vector tangent = rotate(Vector(1, 0, 0), world_q);

      auto frag_fn = [&](const Vector &loc, Fragment &frag) {
        Vector n_local = torus.normal(loc);
        Vector n_world = rotate(n_local, world_q);
        // Headlight model: light coincides with the viewer, so the view vector
        // `center` serves as both light_dir and view_dir.
        float shade = shade_blinn_phong(n_world, center, center, tangent,
                                        params.diffuse, params.specular,
                                        params.fresnel);

        float ring_angle = (atan2f(loc.z, loc.x) + PI_F) / (2.0f * PI_F);
        float palette_t = fmodf(ring_angle + palette_phase +
                                    static_cast<float>(i) / active_count,
                                1.0f);
        Color4 c = baked_palette.get(palette_t);
        frag.color = Color4(c.color * shade, opacity);
      };

      Scan::TransformedVolume vol(torus, center, world_q);
      Scan::Volume::draw<W, H>(pipeline, canvas, center, bounds_radius,
                               ray_dir, vol, frag_fn, max_steps, aa_width);
    }
  }

  /**
   * @brief Tunable shader and animation parameters exposed via registerParam.
   */
  struct Params {
    float pulse_speed = 5.0f;
    float fill = 0.75f;
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
  int active_count = 0;       // vertices built (disdyakis dodecahedron = 26)
  std::array<Vector, MAX_POINTS> points;
  std::array<Quaternion, MAX_POINTS> raw_quats;
  std::array<float, MAX_POINTS> nn_angle; // per-vertex nearest-neighbour gap (rad)
  Timeline timeline;
  Pipeline<W, H> pipeline; // Empty — camera rotation applied to inputs
  GenerativePalette palette{GradientShape::STRAIGHT, HarmonyType::COMPLEMENTARY,
                            BrightnessProfile::BELL, SaturationProfile::VIBRANT,
                            219};
  BakedPalette baked_palette;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Raymarch)
