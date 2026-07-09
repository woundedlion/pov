/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

/**
 * @brief Renders concentric rings of a morphing shape in Plot and Scan modes.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Draws nested rings of a polygon/flower/star shape in both the Plot
 * and Scan render modes, cycling the shape type and tumbling the camera.
 */
template <int W, int H> class ShapeShifter : public Effect {
public:
  /**
   * @brief Shape primitives this effect can morph between.
   */
  enum class ShapeType { PlanarPolygon, SphericalPolygon, Flower, Star };
  /**
   * @brief Rasterization paths a ring can be drawn through.
   */
  enum class RenderMode { Plot, Scan };

  /**
   * @brief Maximum rings the raster can resolve: one per row. Rings spaced finer
   * than a row apart collapse, so this bounds both the Count slider and draw_all.
   */
  static constexpr int MAX_RINGS = H > 1 ? H : 1;

  /**
   * @brief Constructs the effect on a WxH canvas.
   * @details Starts on the planar polygon shape.
   */
  FLASHMEM ShapeShifter()
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(plot_filters)::any_crosses_segments ||
                              decltype(scan_filters)::any_crosses_segments}),
        current_shape(ShapeType::PlanarPolygon) {}

  /**
   * @brief Registers tunable params and builds the animation timeline.
   */
  void init() override {
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("Count", &params.num_shapes, 1.0f,
                  static_cast<float>(MAX_RINGS));
    register_param("Radius", &params.radius, 0.1f, 5.0f);
    register_param("Sides", &params.sides, 3.0f, 12.0f);
    register_animated_param("Twist", &params.twist, -5.0f, 5.0f);
    register_param("Debug BB", &params.debug_bb);

    build();
  }

  /**
   * @brief Advances the shape type every 48 frames (unless paused), then steps
   *        the timeline.
   * @details The shape cycle honors anims_paused_ so "Pause Animation" freezes
   *        it alongside the twist/rotation timers (which already carry the
   *        flag). The frame counter is held too, so the cycle resumes mid-phase
   *        rather than jumping on unpause.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    if (!anims_paused_) {
      // Wrap at the 48-frame cycle period to keep the counter bounded.
      frame_count_ = (frame_count_ + 1) % 48;
      if (frame_count_ == 0) {
        int next = (static_cast<int>(current_shape) + 1) % 4;
        current_shape = static_cast<ShapeType>(next);
      }
    }

    timeline.step(canvas);
  }

  /**
   * @brief Builds a fixed five-slot timeline independent of the Count param.
   * @details The slots are the camera tumble, the twist Mutation, two shared
   * orientation tumbles (one per render mode), and a single draw event that
   * walks the live Count each frame. All rings of a mode share one orientation,
   * so Count scales to the raster budget rather than the 64-slot timeline
   * budget.
   */
  FLASHMEM void build() {
    timeline.clear();

    timeline.add(0, Animation::RandomWalk<W>(camera, X_AXIS, noise, {},
                                             hs::rand_int(0, 65536)));

    timeline.add(0, Animation::Mutation(
                        params.twist,
                        [](float t) { return (PI_F / 4.0f) * sinf(t * PI_F); },
                        480, ease_linear, true, &anims_paused_));

    // Shared tumbles: every Plot ring rides plot_orient_ (+X), every Scan ring
    // rides scan_orient_ (-X).
    timeline.add(0, Animation::Rotation<W>(plot_orient_, X_AXIS, 2 * PI_F, 160,
                                           ease_linear, true,
                                           Animation::Space::Local));
    timeline.add(0, Animation::Rotation<W>(scan_orient_, -X_AXIS, 2 * PI_F, 160,
                                           ease_linear, true,
                                           Animation::Space::Local));

    // Single draw event reading Count live. Added last so it steps after the
    // rotations have collapsed.
    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float) { draw_all(canvas); }, -1,
                        0, ease_linear, 0, ease_linear));
  }

  /**
   * @brief Draws all Count rings, outermost first, in both Plot and Scan modes.
   * @param canvas Target canvas the rings are rasterized onto.
   * @details Computes the shared per-mode basis once per frame (all rings of a
   * mode share orientation and a fixed normal) instead of once per ring.
   */
  void draw_all(Canvas &canvas) {
    int count = static_cast<int>(params.num_shapes);
    if (count < 1)
      count = 1;
    if (count > MAX_RINGS)
      count = MAX_RINGS;

    Quaternion cam_q = camera.get();
    Basis plot_basis = make_basis(cam_q * plot_orient_.get(), X_AXIS);
    Basis scan_basis = make_basis(cam_q * scan_orient_.get(), -X_AXIS);

    for (int i = count - 1; i >= 0; --i) {
      // (i+1)/count maps layers to (0, 1]; i/(count-1) would give the i=0 layer
      // radius 0 and render nothing at Count=1.
      float t = static_cast<float>(i + 1) / count;
      Color4 color = Palettes::RICH_SUNSET.get(t);
      draw_ring(canvas, plot_basis, RenderMode::Plot, t, color, i);
      draw_ring(canvas, scan_basis, RenderMode::Scan, t, color, i);
    }
  }

  /**
   * @brief Draws one ring through the current shape's Plot or Scan renderer.
   * @param canvas Target canvas the ring is rasterized onto.
   * @param basis Shared orientation basis for this render mode.
   * @param mode Whether to dispatch to the Plot or Scan renderer.
   * @param scale Radius scale factor in [0, 1] for this layer.
   * @param color Base palette color for the ring, before alpha modulation.
   * @param layer_index Layer ordinal; scales the per-layer twist phase.
   * @details Scales the radius by `scale`, twists the ring by its layer index,
   * and applies a per-fragment alpha shader.
   */
  void draw_ring(Canvas &canvas, const Basis &basis, RenderMode mode, float scale,
                const Color4 &color, int layer_index) {
    auto fragment_shader = [&](const Vector &, Fragment &f) {
      Color4 c = color;
      c.alpha = c.alpha * this->params.alpha;
      f.color = c;
    };

    float phase = layer_index * this->params.twist;
    float r = this->params.radius * scale;
    int sides_int = (int)params.sides;
    if (mode == RenderMode::Plot) {
      dispatch_plot(canvas, basis, r, sides_int, fragment_shader, phase);
    } else {
      dispatch_scan(canvas, basis, r, sides_int, fragment_shader, phase);
    }
  }

private:
  /**
   * @brief Plot-rasterizes the current shape type.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas the shape is rasterized onto.
   * @param basis Orientation basis for the shape.
   * @param r Ring radius in world units.
   * @param sides_int Polygon/flower/star side count.
   * @param fragment_shader Per-fragment shader callable.
   * @param phase Per-layer twist phase in radians.
   * @details Marked noinline to curb code bloat from the per-shape template
   * instantiation.
   */
  template <typename F>
  __attribute__((noinline)) void
  dispatch_plot(Canvas &canvas, const Basis &basis, float r, int sides_int,
               const F &fragment_shader, float phase) {
    switch (current_shape) {
    case ShapeType::Flower:
      Plot::Flower::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                               fragment_shader, {}, phase);
      break;
    case ShapeType::Star:
      Plot::Star::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                             fragment_shader, phase);
      break;
    case ShapeType::PlanarPolygon:
      Plot::PlanarPolygon::draw<W, H>(plot_filters, canvas, basis, r, sides_int,
                                      fragment_shader, phase);
      break;
    default:
      Plot::SphericalPolygon::draw<W, H>(plot_filters, canvas, basis, r,
                                         sides_int, fragment_shader, phase);
      break;
    }
  }

  /**
   * @brief Scan-rasterizes the current shape type.
   * @tparam F Fragment-shader callable type.
   * @param canvas Target canvas the shape is rasterized onto.
   * @param basis Orientation basis for the shape.
   * @param r Ring radius in world units.
   * @param sides_int Polygon/flower/star side count.
   * @param fragment_shader Per-fragment shader callable.
   * @param phase Per-layer twist phase in radians.
   * @details Marked noinline to curb code bloat from the per-shape template
   * instantiation.
   */
  template <typename F>
  __attribute__((noinline)) void
  dispatch_scan(Canvas &canvas, const Basis &basis, float r, int sides_int,
               const F &fragment_shader, float phase) {
    switch (current_shape) {
    case ShapeType::Flower:
      Scan::Flower::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                               fragment_shader, phase, params.debug_bb);
      break;
    case ShapeType::Star:
      Scan::Star::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                             fragment_shader, phase, params.debug_bb);
      break;
    case ShapeType::PlanarPolygon:
      Scan::PlanarPolygon::draw<W, H>(scan_filters, canvas, basis, r, sides_int,
                                      fragment_shader, phase, params.debug_bb);
      break;
    default: // SphericalPolygon
      Scan::SphericalPolygon::draw<W, H>(scan_filters, canvas, basis, r,
                                         sides_int, fragment_shader, phase,
                                         params.debug_bb);
      break;
    }
  }

  FastNoiseLite noise; /**< Noise source driving the camera RandomWalk. */
  Orientation<> camera; /**< Global camera orientation, tumbled each frame. */
  /**
   * @brief Shared Plot-mode tumble orientation; every Plot ring rides it.
   * @details Declared before `timeline` so it outlives the Rotations that point
   * here, which ~Timeline clears on teardown.
   */
  Orientation<> plot_orient_;
  /**
   * @brief Shared Scan-mode tumble orientation; every Scan ring rides it.
   * @details Declared before `timeline` so it outlives the Rotations that point
   * here, which ~Timeline clears on teardown.
   */
  Orientation<> scan_orient_;
  Timeline timeline;   /**< Fixed five-slot animation timeline. */
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> plot_filters; /**< Anti-aliased filter pipeline for Plot mode. */
  Pipeline<W, H> scan_filters; /**< Filter pipeline for Scan mode. */

  /**
   * @brief Tunable rendering parameters exposed to the GUI.
   */
  struct Params {
    float alpha = 0.5f;     /**< Global alpha multiplier in [0, 1]. */
    float num_shapes = 7.0f; /**< Ring count; read live each frame. */
    float radius = 1.0f;     /**< Outermost ring radius in world units. */
    float sides = 5.0f;      /**< Polygon/flower/star side count. */
    float twist = 0.0f;      /**< Per-layer twist phase in radians. */
    bool debug_bb = false;   /**< Draw Scan bounding boxes when true. */
  } params;

  ShapeType current_shape; /**< Shape currently being rendered. */
  int frame_count_ = 0;    /**< Frame phase in [0, 48) — gates shape cycling; never overflows. */
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShapeShifter)
