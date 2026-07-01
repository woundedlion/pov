/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "animation_core.h"

namespace Animation {

/**
 * @brief An animation that smoothly interpolates a GenerativePalette toward a
 * target palette.
 */
class ColorWipe : public AnimationBase<ColorWipe> {
public:
  /**
   * @brief Constructs a ColorWipe animation.
   * @param from_palette The GenerativePalette to animate (snapshotted on the
   * first step(), mirroring Transition).
   * @param to_palette The GenerativePalette to interpolate toward.
   * @param duration The duration in frames.
   * @param easing_fn The easing function.
   */
  ColorWipe(GenerativePalette &from_palette,
            const GenerativePalette &to_palette, int duration,
            EasingFn easing_fn)
      : AnimationBase(duration, false), cur_palette(from_palette),
        to_snap(to_palette.snapshot()), easing_fn(std::move(easing_fn)) {}

  /**
   * @brief Steps the animation, blending the palette's colors based on the time
   * factor.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    // Snapshot the start palette once, on the first step (see Transition::step).
    if (!captured) {
      from_snap = cur_palette.get().snapshot();
      captured = true;
    }
    AnimationBase::step(canvas);
    float amount = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    cur_palette.get().lerp(from_snap, to_snap, easing_fn(amount));
  }

private:
  GenerativePalette::Snapshot from_snap{}; /**< Snapshot of starting colors. */
  std::reference_wrapper<GenerativePalette>
      cur_palette;                     /**< The palette being animated. */
  GenerativePalette::Snapshot to_snap; /**< Snapshot of target colors. */
  EasingFn easing_fn;                  /**< Easing curve. */
  bool captured = false; /**< Whether from_snap was taken on the first step. */
};

/**
 * @brief Animates the Mobius parameters for a continuous loxodromic flow.
 */
class MobiusFlow : public AnimationBase<MobiusFlow> {
public:
  /**
   * @brief Constructs a MobiusFlow animation.
   * @param params Reference to the MobiusParams to animate.
   * @param num_rings Reference to the number of rings (scalar).
   * @param num_lines Reference to the number of lines (scalar).
   * @param duration Duration of the flow.
   * @param repeat Whether to repeat.
   */
  MobiusFlow(MobiusParams &params, const float &num_rings,
             const float &num_lines, int duration, bool repeat = true)
      : AnimationBase(duration, repeat), params(params), num_rings(num_rings),
        num_lines(num_lines) {
    HS_CHECK(duration >= 0, "MobiusFlow duration must be >= 0");
  }

  // Borrow contract: num_rings/num_lines are read every frame, so they must
  // outlive the Timeline; these deleted overloads reject a temporary scalar.
  MobiusFlow(MobiusParams &params, const float &&num_rings,
             const float &num_lines, int duration, bool repeat = true) = delete;
  MobiusFlow(MobiusParams &params, const float &num_rings,
             const float &&num_lines, int duration, bool repeat = true) = delete;
  MobiusFlow(MobiusParams &params, const float &&num_rings,
             const float &&num_lines, int duration, bool repeat = true) = delete;

  /**
   * @brief Steps the animation, updating params a and d.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    // Floor rings at 0 so the divisor rings + 1 stays >= 1 (num_rings == -1 would
    // make logPeriod inf and poison a/d).
    float rings = num_rings;
    if (rings < 0.0f)
      rings = 0.0f;
    float logPeriod = 5.0f / (rings + 1);
    float flowParam = progress * logPeriod;
    float scale = expf(flowParam);
    float s = sqrtf(scale);
    // Clamp lines to >= 1 before dividing (the slider bottoms out at 0 → 2π/0).
    float lines = num_lines;
    if (lines < 1.0f)
      lines = 1.0f;
    float angle = progress * (2 * PI_F / lines);

    params.get().a.re = s * cosf(angle);
    params.get().a.im = s * sinf(angle);
    params.get().d.re = (1.0f / s) * cosf(-angle);
    params.get().d.im = (1.0f / s) * sinf(-angle);
  }

private:
  std::reference_wrapper<MobiusParams> params; /**< Mobius params to animate. */
  /** @brief Live ring count driving the log period. */
  std::reference_wrapper<const float> num_rings;
  /** @brief Live line count driving the angular step. */
  std::reference_wrapper<const float> num_lines;
};

/**
 * @brief Animates the Mobius parameters for a warping effect pulling the poles
 * together.
 */
class MobiusWarp : public AnimationBase<MobiusWarp> {
public:
  /**
   * @brief Constructs a MobiusWarp animation.
   * @param params Reference to the MobiusParams to animate.
   * @param scale The magnitude of the warp effect.
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   * @param easing The easing function to use (default: ease_in_out_sin).
   */
  MobiusWarp(MobiusParams &params, float scale, int duration,
             bool repeat = true, EasingFn easing = ease_in_out_sin)
      : AnimationBase(duration, repeat), params_(params), scale_(scale),
        easing_(easing) {
    HS_CHECK(duration >= 0, "MobiusWarp duration must be >= 0");
  }

  /**
   * @brief Binds the warp magnitude to a live external float.
   * @param live_scale The external float to read each frame as the magnitude.
   * @details By default the magnitude is the `scale` captured at construction.
   * Binding makes step() read the referent every frame instead, so a GUI slider
   * wired to that float takes effect immediately rather than only at the next
   * (re)spawn. The referent must outlive the animation (e.g. an effect's param
   * member). This avoids retaining the animation pointer across frames, which
   * would dangle under timeline compaction.
   */
  void bind_scale(const float &live_scale) { scale_ref_ = &live_scale; }

  /** @brief Sets the warp magnitude (ignored while a live scale is bound). */
  void set_scale(float scale) { scale_ = scale; }

  /** @brief Sets the easing curve. */
  void set_easing(EasingFn easing) { easing_ = easing; }

  /**
   * @brief Steps the animation, updating param b.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing_(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    float s = scale_ref_ ? *scale_ref_ : scale_;
    params_.get().b.re = s * (cosf(angle) - 1.0f);
    params_.get().b.im = s * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params_; /**< Mobius params to animate. */
  float scale_;                      /**< Warp magnitude. */
  EasingFn easing_;                  /**< Easing curve. */
  const float *scale_ref_ = nullptr; /**< Optional live magnitude source. */
};

/**
 * @brief Animates the Mobius parameters for a circular warping effect.
 */
class MobiusWarpCircular : public AnimationBase<MobiusWarpCircular> {
public:
  /**
   * @brief Constructs a MobiusWarpCircular animation.
   * @param params Reference to the MobiusParams to animate.
   * @param scale The magnitude of the warp effect.
   * @param duration Duration of the warp.
   * @param repeat Whether to repeat.
   * @param easing The easing function to use (default: ease_in_out_sin).
   */
  MobiusWarpCircular(MobiusParams &params, float scale, int duration,
                     bool repeat = true, EasingFn easing = ease_in_out_sin)
      : AnimationBase(duration, repeat), params_(params), scale_(scale),
        easing_(easing) {
    HS_CHECK(duration >= 0, "MobiusWarpCircular duration must be >= 0");
  }

  /** @brief Sets the warp magnitude. */
  void set_scale(float scale) { scale_ = scale; }

  /** @brief Sets the easing curve. */
  void set_easing(EasingFn easing) { easing_ = easing; }

  /**
   * @brief Steps the animation, updating param b.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    float t_norm = static_cast<float>(t) / duration;
    float progress = easing_(hs::clamp(t_norm, 0.0f, 1.0f));
    float angle = progress * 2 * PI_F;
    params_.get().b.re = scale_ * cosf(angle);
    params_.get().b.im = -scale_ * sinf(angle);
  }

private:
  std::reference_wrapper<MobiusParams> params_; /**< Mobius params to animate. */
  float scale_;    /**< Warp magnitude. */
  EasingFn easing_; /**< Easing curve. */
};



/**
 * @brief Animates a vertex-interpolated crossfade between two meshes.
 * @details Owns its transient state (cloned meshes + SLERP buffers) via a
 * pointer to arena-allocated storage — keeps inline size small for
 * TimelineEvent. The caller provides source/dest MeshState references and an
 * Arena; MeshMorph clones both, builds nearest-vertex correspondence, and
 * interpolates each frame. Destruction runs only member destructors — there is
 * no destructor that reclaims the arena, so the transient bytes persist in the
 * arena until the caller compacts or resets it.
 *
 * Crossfade contract (intentional): only the incoming mesh (mesh_B, carrying
 * dest topology) morphs — each frame its vertices SLERP from their nearest
 * source vertex toward their dest position. The outgoing mesh (mesh_A, a clone
 * of source) holds its geometry and fades out via opacity (op_A = 1 - alpha)
 * while the incoming fades in; the opacities sum to 1 for constant total
 * brightness across the topology swap. The source is cloned (not borrowed) so
 * the animation is self-contained against the caller recycling/mutating source
 * mid-morph, consistent with the borrow contract enforced on the draw callbacks
 * below. The separate draw_outgoing/draw_incoming callbacks exist precisely to
 * shade the two halves independently (see HankinSolids).
 */
class MeshMorph : public AnimationBase<MeshMorph> {
public:
  /**
   * @brief Non-owning per-half draw callback: `void(Canvas&, const MeshState&,
   * float opacity)`. A StoredFunctionRef (not a plain FunctionRef) because it is
   * held as a member and invoked across many frames: the type rejects rvalue
   * temporaries so a dangling inline lambda is a compile error, not a silent
   * use-after-free (see the borrow contract below).
   */
  using MorphDrawFn = StoredFunctionRef<void(Canvas &, const MeshState &, float)>;

  /**
   * @brief Constructs a MeshMorph with separate shading for the two halves.
   * @param source The outgoing mesh (cloned, not borrowed).
   * @param dest The incoming mesh whose topology the morph targets.
   * @param arena Arena providing backing storage for cloned meshes and buffers.
   * @param draw_outgoing Draw callback for the fading-out source clone.
   * @param draw_incoming Draw callback for the fading-in morphing mesh.
   * @param duration The crossfade duration in frames.
   * @param easing_fn The easing function applied to crossfade progress.
   * @note The cloned meshes and position buffers are arena-allocated with no
   *   per-instance reclamation; the caller must compact the arena between
   *   successive morphs or it grows unbounded (see HankinSolids/MeshFeedback).
   */
  MeshMorph(const MeshState &source, const MeshState &dest, Arena &arena,
            MorphDrawFn draw_outgoing, MorphDrawFn draw_incoming, int duration,
            EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(duration, false), easing_fn(easing_fn),
        draw_outgoing(draw_outgoing), draw_incoming(draw_incoming) {
    HS_CHECK(!source.vertices.is_empty());
    HS_CHECK(!dest.vertices.is_empty());
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();

    MeshOps::clone(source, buf_->mesh_A, arena);
    MeshOps::clone(dest, buf_->mesh_B, arena);

    // step() indexes mesh_B.vertices by dest's vertex count, so trap here if
    // clone() ever welds/dedups rather than writing OOB per-frame.
    HS_CHECK(buf_->mesh_B.vertices.size() == dest.vertices.size());

    buf_->start_pos.bind(arena, dest.vertices.size());
    buf_->end_pos.bind(arena, dest.vertices.size());

    // Symmetry-breaking twist to avoid degenerate nearest-vertex mapping
    Vector twist_axis = Vector(0.0f, 0.0f, 1.0f);
    bool has_poles = false;
    for (const auto &v : source.vertices) {
      if (std::abs(v.z) > 0.99f && std::abs(v.x) < 0.01f && std::abs(v.y) < 0.01f)
        has_poles = true;
    }
    if (has_poles) {
      twist_axis = Vector(1.0f, 1.0f, 1.0f).normalized();
    }
    Quaternion twist = make_rotation(twist_axis, 0.05f);

    // Build nearest-vertex correspondence: an O(V_dest * V_source) brute force,
    // run once at construction. Matched by greatest dot product against the
    // twist-biased dest vertex (the twist breaks ties on symmetric meshes).
    for (size_t i = 0; i < dest.vertices.size(); ++i) {
      Vector v_biased = rotate(dest.vertices[i], twist);
      int best_idx = 0;
      float max_dot = -9999.0f;
      for (size_t j = 0; j < source.vertices.size(); ++j) {
        float d = dot(v_biased, source.vertices[j]);
        if (d > max_dot) {
          max_dot = d;
          best_idx = static_cast<int>(j);
        }
      }
      buf_->start_pos.push_back(source.vertices[best_idx]);
      buf_->end_pos.push_back(dest.vertices[i]);
    }
  }

  // Borrow contract: the draw callbacks are non-owning StoredFunctionRefs read
  // every frame, so they must outlive the timeline; StoredFunctionRef rejects a
  // temporary at the MorphDrawFn parameter, so no `= delete` overload is needed.

  /**
   * @brief Steps the crossfade: interpolates vertices and renders both halves.
   * @param canvas The canvas buffer passed to the draw callbacks.
   */
  void step(Canvas &canvas) override {
    // Increment-first so the final frame lands exactly on the destination mesh
    // (alpha == 1); the skipped progress==0 frame is immaterial.
    AnimationBase::step(canvas);

    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    float alpha = easing_fn(progress);

    for (size_t i = 0; i < buf_->end_pos.size(); ++i) {
      buf_->mesh_B.vertices[i] =
          slerp(buf_->start_pos[i], buf_->end_pos[i], alpha);
    }

    float op_A = 1.0f - alpha;
    if (op_A > 0.01f)
      draw_outgoing(canvas, buf_->mesh_A, op_A);
    if (alpha > 0.01f)
      draw_incoming(canvas, buf_->mesh_B, alpha);
  }

private:
  /**
   * @brief Arena-allocated transient data — keeps MeshMorph inline size small.
   */
  struct Transients {
    MeshState mesh_A;            /**< Outgoing mesh clone. */
    MeshState mesh_B;           /**< Incoming morphing mesh clone. */
    ArenaVector<Vector> start_pos; /**< Per-vertex nearest-source start points. */
    ArenaVector<Vector> end_pos;   /**< Per-vertex dest end points. */
  };

  Transients *buf_;          /**< Pointer to arena-allocated transient state. */
  EasingFn easing_fn;        /**< Easing curve applied to crossfade progress. */
  MorphDrawFn draw_outgoing; /**< Draw callback for the outgoing half. */
  MorphDrawFn draw_incoming; /**< Draw callback for the incoming half. */
};

/**
 * @brief Continuously modulates Mobius parameters to create an evolving warp.
 * @details Uses multiple frequencies for non-repeating chaos.
 *
 * PERPETUAL — never completes: the default AnimationBase ctor gives it
 * duration = -1 with no repeat, so it never reaches done() and a `.then()`
 * callback attached to it NEVER fires (same hazard documented on RandomWalk).
 * Do not chain sequencing logic off its `.then()`, and note the Transformer's
 * slot-recycling `.then()` likewise will not fire for a spawned
 * MobiusWarpEvolving — drive follow-on behavior from a finite animation or
 * cancel() it explicitly.
 */
class MobiusWarpEvolving : public AnimationBase<MobiusWarpEvolving> {
public:
  /**
   * @brief Constructs a MobiusWarpEvolving animation.
   * @param params The params to animate.
   * @param scale Magnitude of modulation.
   * @param speed Speed of the animation.
   * @note `base` snapshots `params` at construction; MobiusParams has no
   * refresh_from, so Transformer::prepare_frame() never re-reads it from
   * template_params. The baseline is latched at spawn — live edits to it
   * require a respawn (live `scale`/`speed` are mutated via the setters).
   */
  MobiusWarpEvolving(MobiusParams &params, float scale = 0.5f,
                     float speed = 0.01f)
      : params_(params), speed_(speed), scale_(scale), base(params),
        seed(hs::random()()) {}

  /** @brief Sets the modulation speed (radians of phase per frame unit). */
  void set_speed(float speed) { speed_ = speed; }

  /** @brief Sets the per-channel modulation magnitude. */
  void set_scale(float scale) { scale_ = scale; }

  /**
   * @brief Derives a per-channel phase offset from the seed and channel index.
   * @param i Channel index.
   * @return Phase offset in [0, 100) for that channel.
   */
  float phase(int i) const {
    uint32_t h = seed ^ (static_cast<uint32_t>(i) * 2654435761u);
    return (h & 0xFFFF) * (100.0f / 65536.0f);
  }

  /**
   * @brief Steps the animation, modulating all eight Mobius params.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    // Accepted limit: past t == 2^24 (~77 h at 60 fps, sooner at higher speed)
    // float can't represent consecutive frames and the phase freezes; this
    // animation is perpetual (duration == -1). A modulo/accumulator fix only
    // trades the slow drift for a sharp artifact.
    float time = t * speed_;
    float s = scale_;

    // Use prime-ish number ratios for frequencies to minimize repetition cycle
    params_.get().a.re = base.a.re + sinf(time * 1.0f + phase(0)) * s;
    params_.get().a.im = base.a.im + cosf(time * 1.13f + phase(1)) * s;

    params_.get().b.re = base.b.re + sinf(time * 1.27f + phase(2)) * s;
    params_.get().b.im = base.b.im + cosf(time * 1.39f + phase(3)) * s;

    params_.get().c.re = base.c.re + sinf(time * 0.71f + phase(4)) * s;
    params_.get().c.im = base.c.im + cosf(time * 0.83f + phase(5)) * s;

    params_.get().d.re = base.d.re + sinf(time * 0.97f + phase(6)) * s;
    params_.get().d.im = base.d.im + cosf(time * 1.09f + phase(7)) * s;
  }

private:
  std::reference_wrapper<MobiusParams> params_; /**< Mobius params to animate. */
  float speed_; /**< Animation speed (radians of phase per frame unit). */
  float scale_; /**< Magnitude of the per-channel modulation. */
  MobiusParams base; /**< Baseline params captured at construction. */
  uint32_t seed;     /**< Seed for the per-channel phase offsets. */
};

/**
 * @brief Parameters for a ripple wave effect.
 */
struct RippleParams {
  Vector center;          /**< Center point of the ripple source. */
  float amplitude = 0.0f; /**< Current height of the wave. */
  float phase = 0.0f;     /**< Current phase offset (time). */
  float decay{5.0};      /**< Spatial decay rate. */
  float thickness{1.0f}; /**< Thickness of the ripple. */

  /** @brief Cached cos(angle) lower fast-reject bound. */
  float cos_threshold_min = 1.0f;
  /** @brief Cached cos(angle) upper fast-reject bound. */
  float cos_threshold_max = -1.0f;

  /**
   * @brief Recomputes the cos(angle) fast-reject bounds for the wavelet's
   * current phase and thickness, so the renderer can skip points outside the
   * ring.
   */
  void prepare_thresholds() {
    float hw = thickness * 0.5f;
    if (hw < 0.001f)
      hw = 0.001f;
    // Clamp into [0,π]: the active ring lies within the sphere's angular range,
    // so cos(clamped) keeps the fast-reject band engaged past phase=π instead of
    // collapsing both bounds to accept-all. cos(0)=1 and cos(π)=-1 reproduce the
    // out-of-range sentinels at the endpoints.
    float d_min = hs::clamp(phase - hw * 2.0f, 0.0f, PI_F);
    float d_max = hs::clamp(phase + hw * 2.0f, 0.0f, PI_F);
    cos_threshold_min = cosf(d_min);
    cos_threshold_max = cosf(d_max);
  }
};

/**
 * @brief Animates a single ripple event: expanding outward and fading away.
 */
class Ripple : public AnimationBase<Ripple> {
public:
  /**
   * @brief Constructs a Ripple animation.
   * @param params Reference to the params struct to animate.
   * @param center The center point of the ripple.
   * @param speed How fast the waves travel.
   * @param duration How long the ripple lasts in frames.
   */
  Ripple(RippleParams &params, const Vector &center, float speed = 0.2f,
         int duration = 100)
      : AnimationBase(duration, false), params(params), speed(speed),
        peak_amplitude(params.amplitude) {
    HS_CHECK(duration >= 0, "Ripple duration must be >= 0");
    this->params.get().center = center;
    this->params.get().phase = 0.0f;
    // Start at 0 to prevent 1-frame singularity before first step()
    this->params.get().amplitude = 0.0f;
  }

  /**
   * @brief Steps the ripple: advances the wave and updates its envelope.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);

    float progress = static_cast<float>(t) / duration;
    float envelope = 0.0f;

    // Past the duration the envelope is pinned to 0, so skip the wave work; the
    // amplitude is still published below so the last frame renders nothing.
    if (t < static_cast<uint32_t>(duration)) {
      params.get().phase += speed;

      // Attack: fast ramp over ~10% of duration, squared for a parabolic ease-in.
      float attack_dur = 0.1f;
      float attack = std::min(progress / attack_dur, 1.0f);
      attack = attack * attack;

      float decay = 1.0f - progress;

      envelope = attack * decay;

      // Re-prepare the reject thresholds against the phase just advanced, so the
      // render never tests the new phase against thresholds cached at the old one.
      params.get().prepare_thresholds();
    }

    params.get().amplitude = peak_amplitude * envelope;
  }

private:
  std::reference_wrapper<RippleParams> params; /**< Ripple params to animate. */
  float speed;          /**< Wave travel speed (phase increment per frame). */
  float peak_amplitude; /**< Peak wave height captured at construction. */
};

/**
 * @brief Parameters for noise transformation.
 */
struct NoiseParams {
  float amplitude = 0.5f;   /**< Noise output amplitude. */
  float speed = 1.0f;       /**< Temporal evolution speed. */
  float frequency = 0.125f; /**< Spatial frequency of the noise. */
  float time = 0.0f;        /**< Current animation time. */
  float scale = 4.0f;       /**< Spatial scale factor. */
  mutable FastNoiseLite noise; /**< Backing generator; mutable for lazy
                                  init/updates. */

  /**
   * @brief Constructs noise params with an OpenSimplex2 generator.
   */
  NoiseParams() { noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2); }

  /**
   * @brief Syncs the generator's frequency with the `frequency` field.
   */
  void sync() const { noise.SetFrequency(frequency); }

  /**
   * @brief Refreshes live-tunable config from a template snapshot.
   * @param t Template params carrying the current slider values.
   * @details Copies the slider-driven fields (amplitude, speed, frequency,
   * scale) but not the animation-advanced `time` axis or the backing generator,
   * so a live edit reaches an already-spawned entity without resetting its
   * phase. prepare_frame() invokes this before sync().
   */
  void refresh_from(const NoiseParams &t) {
    amplitude = t.amplitude;
    speed = t.speed;
    frequency = t.frequency;
    scale = t.scale;
  }
};

/**
 * @brief Animates noise parameters by updating time.
 */
class Noise : public AnimationBase<Noise> {
public:
  /**
   * @brief Constructs a Noise animation.
   * @param params Reference to the NoiseParams to animate.
   * @param duration Duration in frames (-1 for indefinite). Finite durations
   *   are not recommended: this animation repeats, so a finite duration rewinds
   *   t to 0 each cycle and snaps params.time backward — a sawtooth that breaks
   *   the smooth forward flow the time axis exists for. Use -1 (the default).
   */
  Noise(NoiseParams &params, int duration = -1)
      : AnimationBase(duration, true), params(params) {}

  /**
   * @brief Steps the animation, advancing the noise time field.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    // Accepted limit: past t == 2^24 (~77 h at 60 fps, sooner at higher speed)
    // float can't represent consecutive frames and the noise time axis freezes.
    // A modulo/accumulator fix only trades the slow drift for a sharp artifact.
    params.get().time = static_cast<float>(t);
  }

private:
  std::reference_wrapper<NoiseParams> params; /**< Noise params to animate. */
};

} // namespace Animation
