/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Scatters polygon "stars" over a Fibonacci spiral on the sphere and
 *        warps the field with an evolving Möbius transform.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every spiral point is pushed through an evolving Möbius warp so the
 *          field drifts and inflates, while a Languid RandomWalk slowly
 *          reorients the whole field.
 */
template <int W, int H> class GnomonicStars : public Effect {
public:
  /**
   * @brief Live-tunable controls for the star field.
   * @details Star count, per-star radius, polygon side count, and a toggle that
   *          draws each star's bounding box for debugging.
   */
  struct Params {
    float points = 600.0f;    /**< Number of stars scattered on the spiral. */
    float star_radius = 0.02f; /**< Per-star radius in normalized sphere units. */
    float star_sides = 4.0f;  /**< Polygon side count per star. */
    float warp_speed = 0.035f; /**< Möbius warp evolution speed, mirrored into the pinned warp each frame. */
    bool debug_bb = false;    /**< When true, draws each star's bounding box. */
  } params;

  /**
   * @brief Constructs the effect at face resolution W x H.
   * @details Default-initializes the orientation and timeline and binds the
   *          Möbius transformer to the timeline.
   */
  FLASHMEM GnomonicStars()
      : Effect(W, H), orientation(), timeline(), transformer(timeline) {}

  /**
   * @brief Registers the user params and arms the timeline.
   * @details Exposes the user params and arms the timeline with an infinite
   *          Möbius warp whose speed is bound to a slider, plus a Languid
   *          RandomWalk that reorients the field.
   */
  void init() override {
    // Sized to MAX_POINTS so a live "Points" change never reallocates.
    spiral_cache_ = static_cast<Vector *>(persistent_arena.allocate(
        MAX_POINTS * sizeof(Vector), alignof(Vector)));

    registerParam("Points", &params.points, 100.0f, 2000.0f);
    registerParam("Radius", &params.star_radius, 0.01f, 0.1f);
    registerParam("Sides", &params.star_sides, 3.0f, 8.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Args are (scale, speed): fixed 0.5 magnitude; speed is a don't-care here
    // since draw_frame mirrors params.warp_speed into the warp every frame.
    warp_ = transformer.spawn_pinned(0, 0.5f, 0.0f);
    HS_CHECK(warp_, "GnomonicStars: pinned warp spawn must succeed");
    registerParam("Warp Speed", &params.warp_speed, 0.0f, 1.0f);

    baked_palette.bake(persistent_arena, Palettes::mangoPeel);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /// POV column-strobe flag; strobes (see Effect::strobe_columns).
  bool strobe_columns() const override { return true; }

  /**
   * @brief Advances the timeline and renders the warped star field.
   * @details Draws each spiral point as a star whose color is a Y-based gradient
   *          and whose basis carries the current orientation and warp.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Mirror the slider into the warp before the timeline advances it.
    warp_->set_speed(params.warp_speed);

    timeline.step(canvas);

    auto fragment_shader = [this](const Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // Y in [-1, 1] -> gradient t in [0, 1]
      Color4 c = baked_palette.get(t);
      frag.color = c;
    };

    // Clamp to [1, MAX_POINTS]: a sub-1 value would no-op or desync the cache,
    // and MAX_POINTS is the spiral-cache capacity.
    const int points = hs::clamp((int)params.points, 1, MAX_POINTS);
    const float radius = params.star_radius;
    const int sides = (int)params.star_sides;

    // The base spiral depends only on (points, i) — the warp and orientation
    // animate downstream — so rebuild the trig-heavy fib_spiral only when
    // "Points" changes.
    if (points != cached_points_) {
      for (int i = 0; i < points; i++) {
        spiral_cache_[i] = fib_spiral(points, /*eps=*/0.0f, i);
      }
      cached_points_ = points;
    }

    for (int i = 0; i < points; i++) {
      Vector v = transformer.transform(spiral_cache_[i]);

      // make_basis() rotates its normal by the orientation; pass the raw warp
      // output so the orientation is applied exactly once, not twice.
      Basis basis = make_basis(orientation.get(), v);

      Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                             fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  /** @brief Spiral-cache capacity; equals the "Points" slider's upper bound. */
  static constexpr int MAX_POINTS = 2000;

  Orientation<> orientation;        /**< Current field orientation quaternion. */
  FastNoiseLite noise;              /**< Noise source driving the RandomWalk. */
  Timeline timeline;                /**< Animation timeline for warp and walk. */
  Pipeline<W, H> filters;           /**< Render filter pipeline for star scan. */

  Vector *spiral_cache_ = nullptr;  /**< Persistent base lattice, MAX_POINTS slots. */
  int cached_points_ = 0;           /**< Point count the cache holds (0 = unbuilt). */
  BakedPalette baked_palette;       /**< LUT-baked mangoPeel sampled by the shader. */

  MobiusWarpGnomonicTransformer<1> transformer; /**< Evolving Möbius warp applied per point. */
  Animation::MobiusWarpEvolving *warp_ = nullptr; /**< Pinned warp handle; mirrors params.warp_speed each frame. */
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(GnomonicStars)
