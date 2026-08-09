/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ShadierBall.h
 * @brief Slot-based sphere shader: discrete function, projection, and lens
 *        slots with continuous per-slot params, sequenced as presets.
 */

#include "core/engine/engine.h"

// Unit-test accessor for the accumulators, pipeline stages, and palette walk.
namespace hs_test {
namespace shadierball_tests {
struct ShadierBallWhiteBox;
} // namespace shadierball_tests
} // namespace hs_test

/**
 * @brief Sphere shader assembled from preset-selected slots.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each preset names a pattern Function, a sphere-to-plane
 * Projection, and a sphere Lens — discrete slot tags dispatched per pixel on
 * frame-constant copies — plus the continuous params those slots consume. A
 * random-walk camera drifts the whole view at a preset-set gain. Color comes
 * from a PaletteCycler walking the 256 prebaked triadic profiles on a
 * golden-ratio hue step.
 */
template <int W, int H> class ShadierBall : public Effect {
public:
  HS_COLD_MEMBER ShadierBall() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Registers slot and param controls and starts the palette walk.
   */
  HS_COLD_MEMBER void init() override {
    // register_param captures *ptr as the GUI default, so the live preset has
    // to be loaded before any registration below.
    active = PRESETS[0];
    // Slot tags register animated like every preset-driven dropdown: writes
    // engage the pause, and the future preset choreography owns them.
    register_animated_param("Function", &active.slots.function,
                            FUNCTION_OPTIONS, FUNCTION_EXPORT_OPTIONS,
                            NUM_FUNCTIONS);
    register_animated_param("Projection", &active.slots.projection,
                            PROJECTION_OPTIONS, PROJECTION_EXPORT_OPTIONS,
                            NUM_PROJECTIONS);
    register_animated_param("Lens", &active.slots.lens, LENS_OPTIONS,
                            LENS_EXPORT_OPTIONS, NUM_LENSES);
    Params &params = active.params;
    register_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_param("Pattern Freq", &params.pattern_freq, PATTERN_FREQ_MIN,
                   PATTERN_FREQ_MAX);
    register_param("Wave Spin", &params.wave_spin, WAVE_SPIN_MIN,
                   WAVE_SPIN_MAX);
    register_param("Pole Fade", &params.pole_fade, POLE_FADE_MIN,
                   POLE_FADE_MAX);
    register_param("Lens Mix", &params.lens_mix, LENS_MIX_MIN, LENS_MIX_MAX);
    register_param("Wander", &params.wander, WANDER_MIN, WANDER_MAX);

    // The walk runs permanently on a shadow orientation; `wander` scales its
    // per-frame motion into the camera (update_camera).
    timeline.add(0, Animation::RandomWalk<W>(walk, UP, walk_noise));

    // No pause flag: "Pause Animation" stops preset advancement (the future
    // choreography), never the palette cycle.
    palette_cycler.init_generated(persistent_arena, next_triadic_palette,
                                  &palette_hue, PALETTE_DWELL_FRAMES,
                                  PALETTE_FADE_FRAMES, ease_in_out_sin);
  }

  /**
   * @brief Advances the walk and wave phases and shades every frame's pixels.
   * @details Steps the timeline, wraps the phase accumulators, steps the
   * palette walk, integrates the wander camera, hoists the frame-constant slot
   * tags and wave state, then shades every pixel: camera, lens, projection,
   * pattern function, pole normalize, palette.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sdb_timeline_step);
      timeline.step(canvas);
    }
    const Params &P = active.params;
    // Wrap the trig phases to 2pi so fast_sinf/fast_cosf keep precise range
    // reduction.
    phase = fmodf(phase + P.speed, TWO_PI_F);
    wave_angle = fmodf(wave_angle + P.wave_spin, TWO_PI_F);

    palette_cycler.step();
    update_camera(P.wander);

    // Frame-constant hoists; the per-pixel slot switches test these copies, so
    // the branch predictor makes the dispatch nearly free.
    const Slots slots = active.slots;
    const WaveState wave{phase, wave_angle, fast_cosf(wave_angle),
                         fast_sinf(wave_angle)};
    const Quaternion view_conj = cam_conj;
    const float pattern_freq = P.pattern_freq;
    const float pole_fade = P.pole_fade;
    const float lens_mix = P.lens_mix;

    auto shader = [&](const Vector &v) -> Color4 {
      const Sample s = sample_pipeline(rotate(v, view_conj), slots, wave,
                                       pattern_freq, pole_fade, lens_mix);
      Color4 color = palette_cycler.palette().get(s.value);
      color.alpha *= s.alpha;
      return color;
    };
    {
      HS_PROFILE(sdb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
  }

private:
  // Test seam: reaches the accumulators, pipeline stages, and palette walk.
  friend struct ::hs_test::shadierball_tests::ShadierBallWhiteBox;

  /** @brief Pattern-field slot tag; one body per enumerator. */
  enum class Function : uint8_t { TWIN_WAVE, RINGS, SPIRAL, GRID };
  /** @brief Sphere-to-plane projection slot tag. */
  enum class Projection : uint8_t { EQUIRECTANGULAR, STEREOGRAPHIC, GNOMONIC };
  /** @brief Sphere-lens slot tag applied before projection. */
  enum class Lens : uint8_t { NONE, GLITCH, TWIST, KALEIDOSCOPE };

  /** @brief Slot tags one preset selects. */
  struct Slots {
    Function function;     /**< Pattern field evaluated per pixel. */
    Projection projection; /**< Sphere-to-plane map feeding the pattern. */
    Lens lens;             /**< Sphere lens ahead of the projection. */
  };

  /** @brief Continuous params the slots consume, one snapshot per preset. */
  struct Params {
    float speed;        /**< Wave phase advance, rad/frame. */
    float pattern_freq; /**< Spatial frequency of the pattern field. */
    float wave_spin;    /**< Second-wave rotation rate, rad/frame. */
    float pole_fade;    /**< Pole attenuation radius. */
    float lens_mix;     /**< Lens blend in [0, 1]; 0 and 1 skip a projection. */
    float wander;       /**< Random-walk camera gain in [0, 1]. */
  };

  /** @brief One sequencable look: slot tags plus their params. */
  struct Preset {
    Slots slots;
    Params params;
  };

  /** @brief Frame-constant travel and rotation state every function reads. */
  struct WaveState {
    float phase;     /**< Travel phase in [0, 2pi). */
    float angle;     /**< Rotation angle in [0, 2pi). */
    float angle_cos; /**< Cosine of the rotation angle. */
    float angle_sin; /**< Sine of the rotation angle. */
  };

  /** @brief Inter-stage record the function stage hands the palette stage. */
  struct Sample {
    float value; /**< Pole-normalized pattern value in [0, 1]. */
    float alpha; /**< Perceptually compensated pole opacity in [0, 1]. */
  };

  /** @brief Full per-pixel pipeline: lens, projection, function, normalize. */
  static Sample sample_pipeline(const Vector &v, const Slots &slots,
                                const WaveState &wave, float pattern_freq,
                                float pole_fade, float lens_mix) {
    const Complex z = project(v, slots.projection, slots.lens, lens_mix);
    const float r_sq = z.re * z.re + z.im * z.im;
    const float pattern = sample_function(
        slots.function, stereo_pattern_args(z, pattern_freq), wave);
    const float pole_weight = pole_attenuation(r_sq, pole_fade);
    return {(pattern * pole_weight + 1.0f) * 0.5f, pole_weight * pole_weight};
  }

  /** @brief Lens dispatch and projection to the pattern plane. */
  static Complex project(const Vector &v, Projection projection, Lens lens,
                         float lens_mix) {
    if (lens == Lens::NONE || lens_mix == 0.0f)
      return project_point(v, projection);
    const Complex lensed = project_point(apply_lens(v, lens), projection);
    if (lens_mix == 1.0f)
      return lensed;
    const Complex direct = project_point(v, projection);
    return {hs::lerp(direct.re, lensed.re, lens_mix),
            hs::lerp(direct.im, lensed.im, lens_mix)};
  }

  /** @brief Lens-slot dispatch on a frame-constant tag. */
  static Vector apply_lens(const Vector &v, Lens lens) {
    switch (lens) {
    case Lens::NONE:
      return v;
    case Lens::GLITCH:
      return glitch_lens(v);
    case Lens::TWIST:
      return twist_lens(v);
    case Lens::KALEIDOSCOPE:
      return kaleidoscope_lens(v);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Azimuthal twist: rotates (x, z) about Y by TWIST_RATE * y.
   * @param v Unit direction vector on the sphere.
   * @return Twisted unit direction; the equator is fixed, the poles turn most.
   */
  static Vector twist_lens(const Vector &v) {
    const float angle = TWIST_RATE * v.y;
    const float c = fast_cosf(angle);
    const float s = fast_sinf(angle);
    return Vector(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
  }

  /**
   * @brief Kaleidoscope: mirror-folds the azimuth into one half-sector wedge.
   * @param v Unit direction vector on the sphere.
   * @return Folded unit direction; the pattern repeats mirrored
   * 2 * KALEIDOSCOPE_SECTORS times around the Y axis.
   */
  static Vector kaleidoscope_lens(const Vector &v) {
    constexpr float SECTOR = TWO_PI_F / KALEIDOSCOPE_SECTORS;
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    float azimuth = fmodf(fast_atan2(v.z, v.x) + PI_F, SECTOR);
    if (azimuth > 0.5f * SECTOR)
      azimuth = SECTOR - azimuth;
    return Vector(radius * fast_cosf(azimuth), v.y,
                  radius * fast_sinf(azimuth));
  }

  /** @brief Projection-slot dispatch: unit direction to plane coordinates. */
  static Complex project_point(const Vector &v, Projection projection) {
    switch (projection) {
    case Projection::EQUIRECTANGULAR:
      return equirectangular(v);
    case Projection::STEREOGRAPHIC:
      return stereo(v);
    case Projection::GNOMONIC:
      return gnomonic(v);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Seam-free equirectangular-style map of a unit direction.
   * @param v Unit direction vector on the sphere.
   * @return re = |azimuth| * cos(latitude) in [0, pi], im = latitude in
   * [-pi/2, pi/2].
   * @details A raw longitude coordinate tears at the +-pi wrap (the z = 0
   * half-plane) and spins without bound at the poles, which a lens can drag
   * anywhere (the glitch lens maps the whole equator onto a pole). Folding the
   * azimuth mirrors the pattern across the wrap instead of tearing, and the
   * cos(latitude) taper — the cylindrical radius, no trig — collapses the
   * azimuth singularity so both poles stay continuous.
   */
  static Complex equirectangular(const Vector &v) {
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    return {std::fabs(fast_atan2(v.z, v.x)) * radius,
            0.5f * PI_F - fast_acos(v.y)};
  }

  /**
   * @brief Gnomonic map through the sphere center onto the y = 1 plane.
   * @param v Unit direction vector on the sphere.
   * @return Plane coordinates; antipodal directions map to the same point.
   */
  static Complex gnomonic(const Vector &v) {
    // Floor |y| so the equator ring divides by the epsilon instead of blowing
    // up to inf/NaN; the resulting ~1e3 coordinates stay inside the
    // pattern-arg clamp.
    float y = v.y;
    if (std::fabs(y) < GNOMONIC_AXIS_EPS)
      y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
    return {v.x / y, v.z / y};
  }

  /** @brief Function-slot dispatch on a frame-constant tag. */
  static float sample_function(Function function, const Complex &p,
                               const WaveState &wave) {
    switch (function) {
    case Function::TWIN_WAVE:
      return twin_wave(p, wave);
    case Function::RINGS:
      return rings(p, wave);
    case Function::SPIRAL:
      return spiral(p, wave);
    case Function::GRID:
      return grid(p, wave);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Two traveling plane waves, the second rotated by the spin angle.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param wave Frame-constant phase and rotation of the pair.
   * @return Mean of the two waves in [-1, 1].
   */
  static float twin_wave(const Complex &p, const WaveState &wave) {
    const float rotated = p.re * wave.angle_cos + p.im * wave.angle_sin;
    return 0.5f *
           (fast_sinf(p.re + wave.phase) + fast_sinf(rotated + wave.phase));
  }

  /**
   * @brief Concentric rings traveling outward from the projection origin.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param wave Frame-constant travel phase; the rotation state is unused.
   * @return Ring value in [-1, 1].
   */
  static float rings(const Complex &p, const WaveState &wave) {
    return fast_sinf(sqrtf(p.re * p.re + p.im * p.im) - wave.phase);
  }

  /**
   * @brief Rotating SPIRAL_ARMS-armed spiral about the projection origin.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param wave Frame-constant travel phase and arm rotation angle.
   * @return Spiral value in [-1, 1]; integer arm count keeps the azimuth
   * seam continuous.
   */
  static float spiral(const Complex &p, const WaveState &wave) {
    const float radius = sqrtf(p.re * p.re + p.im * p.im);
    const float azimuth = fast_atan2(p.im, p.re);
    return fast_sinf(radius - SPIRAL_ARMS * (azimuth + wave.angle) -
                     wave.phase);
  }

  /**
   * @brief Multiplicative grid lattice rotated by the spin angle.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param wave Frame-constant travel phase and lattice rotation.
   * @return Lattice value in [-1, 1].
   */
  static float grid(const Complex &p, const WaveState &wave) {
    const float a = p.re * wave.angle_cos + p.im * wave.angle_sin;
    const float b = -p.re * wave.angle_sin + p.im * wave.angle_cos;
    return fast_sinf(a + wave.phase) * fast_cosf(b - wave.phase);
  }

  /**
   * @brief Integrates a fraction of the walk's per-frame motion into the
   *        camera.
   * @param wander Walk motion gain in [0, 1]; 0 freezes the camera, 1 follows
   * the walk fully.
   * @details Integrating the incremental rotation keeps a partial gain
   * continuous through every full turn and while the gain changes.
   */
  // Per-frame, not per-pixel: the inlined slerp is too big for ITCM.
  HS_COLD_MEMBER void update_camera(float wander) {
    const Quaternion current = walk.get();
    const Quaternion delta = current * walk_prev.conjugate();
    camera =
        (slerp(Quaternion(), delta.normalized(), wander) * camera).normalized();
    walk_prev = current;
    cam_conj = camera.conjugate();
  }

  /** @brief PaletteCycler provider: the triadic profile walking the prebaked
   *  hue wheel; sequence 0 returns the starting hue unchanged. */
  HS_COLD_MEMBER static void next_triadic_palette(void *context,
                                                  uint32_t sequence,
                                                  GenerativePalette &out) {
    uint32_t &hue = *static_cast<uint32_t *>(context);
    if (sequence > 0)
      hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(hue))};
  }

  /** @brief Zero dwell: the palette shifts continuously, fade after fade. */
  static constexpr int PALETTE_DWELL_FRAMES = 0;
  /** @brief Frames each palette fade spans. */
  static constexpr int PALETTE_FADE_FRAMES = 600;
  /** @brief Hue-wheel steps per palette. Odd, so the walk visits all 256
   *  prebaked hues before repeating; 159/256 approximates the golden-ratio
   *  conjugate. */
  static constexpr uint32_t HUE_STEP = 159;

  /** @brief |v.y| floor for the gnomonic division. */
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
  /** @brief Spiral arm count; integral, so the azimuth seam stays continuous. */
  static constexpr float SPIRAL_ARMS = 3.0f;
  /** @brief Twist-lens rotation per unit of latitude height, radians. */
  static constexpr float TWIST_RATE = 3.0f;
  /** @brief Kaleidoscope sector count around the Y axis. */
  static constexpr float KALEIDOSCOPE_SECTORS = 6.0f;

  static constexpr const char *FUNCTION_OPTIONS[] = {"Twin Wave", "Rings",
                                                     "Spiral", "Grid"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE", "Function::RINGS", "Function::SPIRAL",
      "Function::GRID"};
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(std::size(FUNCTION_OPTIONS));
  static constexpr const char *PROJECTION_OPTIONS[] = {
      "Equirectangular", "Stereographic", "Gnomonic"};
  static constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
      "Projection::EQUIRECTANGULAR", "Projection::STEREOGRAPHIC",
      "Projection::GNOMONIC"};
  static constexpr int NUM_PROJECTIONS =
      static_cast<int>(std::size(PROJECTION_OPTIONS));
  static constexpr const char *LENS_OPTIONS[] = {"None", "Glitch", "Twist",
                                                 "Kaleidoscope"};
  static constexpr const char *LENS_EXPORT_OPTIONS[] = {
      "Lens::NONE", "Lens::GLITCH", "Lens::TWIST", "Lens::KALEIDOSCOPE"};
  static constexpr int NUM_LENSES = static_cast<int>(std::size(LENS_OPTIONS));

  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 0.5f;
  static constexpr float PATTERN_FREQ_MIN = 1.0f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float WAVE_SPIN_MIN = 0.0f, WAVE_SPIN_MAX = 0.05f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float LENS_MIX_MIN = 0.0f, LENS_MIX_MAX = 1.0f;
  static constexpr float WANDER_MIN = 0.0f, WANDER_MAX = 1.0f;

  /** @brief True iff every param of @p p lies within its registered range. */
  static constexpr bool preset_in_ranges(const Preset &p) {
    return p.params.speed >= SPEED_MIN && p.params.speed <= SPEED_MAX &&
           p.params.pattern_freq >= PATTERN_FREQ_MIN &&
           p.params.pattern_freq <= PATTERN_FREQ_MAX &&
           p.params.wave_spin >= WAVE_SPIN_MIN &&
           p.params.wave_spin <= WAVE_SPIN_MAX &&
           p.params.pole_fade >= POLE_FADE_MIN &&
           p.params.pole_fade <= POLE_FADE_MAX &&
           p.params.lens_mix >= LENS_MIX_MIN &&
           p.params.lens_mix <= LENS_MIX_MAX && p.params.wander >= WANDER_MIN &&
           p.params.wander <= WANDER_MAX;
  }

  static constexpr Preset PRESETS[] = {
      // Rotating twin-wave interference through the full glitch lens, drifting
      // on a languid quarter-gain walk.
      {{Function::TWIN_WAVE, Projection::STEREOGRAPHIC, Lens::GLITCH},
       {0.05f, 6.0f, 0.006f, 2.0f, 1.0f, 0.25f}},
  };
  static_assert(
      [] {
        for (const Preset &p : PRESETS)
          if (!preset_in_ranges(p))
            return false;
        return true;
      }(),
      "a ShadierBall preset drives a param outside its registered slider "
      "range; widen the range to accommodate the preset");

  /**
   * @brief Camera walk orientation.
   * @details Declared before `timeline` so it outlives the RandomWalk that
   * points here, which ~Timeline clears on teardown.
   */
  Orientation<> walk;
  Timeline timeline; /**< Drives the camera walk. */
  /** @brief OpenSimplex2 source for the walk; RandomWalk's ctor reseeds and
   *  re-frequencies it, so it is never shared. */
  FastNoiseLite walk_noise;

  Quaternion walk_prev; /**< Previous walk orientation. */
  Quaternion camera;    /**< Integrated wander camera. */
  Quaternion cam_conj;  /**< Per-frame inverse of the camera. */

  Preset active = PRESETS[0]; /**< Live slots and params; the GUI edits these
                                 directly. */
  float phase = 0.0f;         /**< Wrapped to [0, 2pi): shared wave travel. */
  float wave_angle = 0.0f;    /**< Wrapped to [0, 2pi): second-wave rotation. */
  uint32_t palette_hue = 0;   /**< Hue-wheel index of the current triadic. */
  PaletteCycler palette_cycler;

  // init() buys the palette cycler's display LUT and morph slots from the
  // persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "ShadierBall persistent footprint exceeds the default "
                "partition");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShadierBall)
