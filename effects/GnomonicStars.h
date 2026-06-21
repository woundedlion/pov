/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

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
    registerParam("Points", &params.points, 100.0f, 2000.0f);
    registerParam("Radius", &params.star_radius, 0.01f, 0.1f);
    registerParam("Sides", &params.star_sides, 3.0f, 8.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Spawn the evolving warp and retain its animation reference for the live
    // "Warp Speed" param. The warp is infinite and added before any other
    // timeline event, so pinning it is the safe retained-handle path: if a
    // future edit ever inserts a finite event ahead of it, compaction traps
    // loudly rather than dangling this pointer the GUI reads every frame.
    auto *anim = transformer.spawn_pinned(0, 0.5f, 0.035f);
    // A pinned spawn into a fresh transformer always has a free slot, so a null
    // here is a structural bug, not a runtime condition. Trap rather than
    // silently drop the "Warp Speed" slider.
    HS_CHECK(anim, "GnomonicStars: pinned warp spawn must succeed");
    registerParam("Warp Speed", &anim->speed, 0.0f, 1.0f);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, UP, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /**
   * @brief Reports whether the effect wants the background drawn.
   * @return Always true; this effect renders over a background.
   */
  bool show_bg() const override { return true; }

  /**
   * @brief Advances the timeline and renders the warped star field.
   * @details Draws each spiral point as a star whose color is a Y-based gradient
   *          and whose basis carries the current orientation and warp.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    timeline.step(canvas);

    // Tints each fragment along a Mango Peel gradient keyed to its world-space Y.
    auto fragment_shader = [](const Vector &p, Fragment &frag) {
      float t = (p.y + 1.0f) * 0.5f; // map Y in [-1, 1] to gradient t in [0, 1]
      Color4 c = Palettes::mangoPeel.get(t);
      frag.color = c;
    };

    const int points = (int)params.points;
    const float radius = params.star_radius;
    const int sides = (int)params.star_sides;

    // The base spiral is recomputed from scratch every frame (fib_spiral is
    // trig-heavy and runs up to 2000 times here). That redundancy is deliberate,
    // not an oversight: only the warp + orientation animate, so the raw lattice
    // could be cached and refreshed only when "Points" changes — but the cache is
    // a fixed 2000-Vector (~24 KB) buffer in the RAM-constrained arena, and every
    // base point depends on `points` (so it invalidates on each slider move). This
    // effect is sim-targeted, where the per-frame trig is cheap, so the CPU work
    // is not worth that standing RAM cost. Add the cache here if it ever ships on
    // a device that runs it hot.
    for (int i = 0; i < points; i++) {
      Vector v = fib_spiral(points, 0.0f, i);

      v = transformer.transform(v);

      // Build the basis at the star position. make_basis() rotates its normal
      // by the orientation, so passing the raw warp output applies the latest
      // orientation exactly once — orienting v separately first would rotate it
      // twice (q*q*v) and drift the field at double the RandomWalk rate.
      Basis basis = make_basis(orientation.get(), v);

      Scan::Star::draw<W, H>(filters, canvas, basis, radius, sides,
                             fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  Orientation<> orientation;        /**< Current field orientation quaternion. */
  FastNoiseLite noise;              /**< Noise source driving the RandomWalk. */
  Timeline timeline;                /**< Animation timeline for warp and walk. */
  Pipeline<W, H> filters;           /**< Render filter pipeline for star scan. */

  MobiusWarpGnomonicTransformer<1> transformer; /**< Evolving Möbius warp applied per point. */
};

#include "core/effect_registry.h"
REGISTER_EFFECT(GnomonicStars)
