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
    float warp_speed = 0.5f;  /**< Möbius warp evolution speed, mirrored into the pinned warp each frame. */
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
    // Carve the base-lattice cache once from the persistent arena (only the
    // active effect uses it). Sized to MAX_POINTS so a live "Points" change never
    // reallocates a bump arena. draw_frame() fills slot [0, points) before
    // reading it, so the uninitialized Vector storage here is always overwritten
    // before use. allocate() traps on OOM, so a non-null return needs no check.
    spiral_cache_ = static_cast<Vector *>(persistent_arena.allocate(
        MAX_POINTS * sizeof(Vector), alignof(Vector)));

    registerParam("Points", &params.points, 100.0f, 2000.0f);
    registerParam("Radius", &params.star_radius, 0.01f, 0.1f);
    registerParam("Sides", &params.star_sides, 3.0f, 8.0f);
    registerParam("Debug BB", &params.debug_bb);

    // Spawn the evolving warp and retain its handle so the "Warp Speed" slider —
    // bound to params.warp_speed like every other control — can be mirrored into
    // the animation each frame (see draw_frame). The warp is infinite and added
    // before any other timeline event, so pinning it is the safe retained-handle
    // path: if a future edit ever inserts a finite event ahead of it, compaction
    // traps loudly rather than dangling this pointer.
    warp_ = transformer.spawn_pinned(0, params.warp_speed, 0.035f);
    // A pinned spawn into a fresh transformer always has a free slot, so a null
    // here is a structural bug, not a runtime condition. Trap rather than
    // silently drop the "Warp Speed" slider.
    HS_CHECK(warp_, "GnomonicStars: pinned warp spawn must succeed");
    registerParam("Warp Speed", &params.warp_speed, 0.0f, 1.0f);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
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

    // Mirror the live Warp Speed slider into the pinned warp animation before the
    // timeline advances it, so a mid-run change takes effect this frame.
    warp_->speed = params.warp_speed;

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

    // The base spiral depends only on (points, i): `eps` is fixed at 0 and only
    // the warp + orientation animate, so the raw lattice is identical every frame
    // until "Points" moves. Rebuild it (trig-heavy fib_spiral, up to MAX_POINTS
    // calls) only on a points change and cache it in the persistent arena; the
    // steady state then just reads each base point back before warping it. The
    // cache is a fixed MAX_POINTS-Vector buffer carved once in init() from the
    // persistent arena, which only the active effect uses.
    HS_CHECK(points <= MAX_POINTS,
             "GnomonicStars: Points exceeds spiral-cache capacity");
    if (points != cached_points_) {
      for (int i = 0; i < points; i++) {
        spiral_cache_[i] = fib_spiral(points, 0.0f, i);
      }
      cached_points_ = points;
    }

    for (int i = 0; i < points; i++) {
      Vector v = transformer.transform(spiral_cache_[i]);

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
  /** @brief Spiral-cache capacity; equals the "Points" slider's upper bound. */
  static constexpr int MAX_POINTS = 2000;

  Orientation<> orientation;        /**< Current field orientation quaternion. */
  FastNoiseLite noise;              /**< Noise source driving the RandomWalk. */
  Timeline timeline;                /**< Animation timeline for warp and walk. */
  Pipeline<W, H> filters;           /**< Render filter pipeline for star scan. */

  Vector *spiral_cache_ = nullptr;  /**< Persistent base lattice, MAX_POINTS slots. */
  int cached_points_ = 0;           /**< Point count the cache holds (0 = unbuilt). */

  MobiusWarpGnomonicTransformer<1> transformer; /**< Evolving Möbius warp applied per point. */
  Animation::MobiusWarpEvolving *warp_ = nullptr; /**< Pinned warp handle; mirrors params.warp_speed each frame. */
};

#include "core/effect_registry.h"
REGISTER_EFFECT(GnomonicStars)
