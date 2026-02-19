/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include "animation.h"
#include "concepts.h"
#include <array>
#include "FastNoiseLite.h"

// Forward declare
class Canvas;

/**
 * @brief A generic manager for state-based geometry transformations.
 * * @tparam W Canvas width (for Timeline).
 * @tparam ParamsT The configuration struct (e.g., RippleParams, MobiusParams).
 * @tparam AnimT The animation class (e.g., Animation::Ripple).
 * @tparam TransformFunc The static function to apply the transformation.
 * @tparam CAPACITY Max number of active transformations.
 */
template <int W, typename ParamsT, typename AnimT, auto TransformFunc,
          int CAPACITY = 32, int TIMELINE_CAPACITY = 32>
  requires TransformerFunction<TransformFunc, ParamsT>
class Transformer {
public:
  struct Entity {
    ParamsT params;
    bool active = false;
  };

  ParamsT params;
  std::array<Entity, CAPACITY> entities;
  Timeline<W, TIMELINE_CAPACITY> &timeline;

  Transformer(Timeline<W, TIMELINE_CAPACITY> &tl) : timeline(tl) {}

  /**
   * @brief Spawns a new transformation animation.
   * @param in_frames Delay before starting.
   * @param args Arguments forwarded to the Animation constructor (after the
   * Params& argument).
   */
  template <typename... Args> void spawn(int in_frames, Args &&...args) {
    // Linear scan for a free slot
    for (auto &e : entities) {
      if (!e.active) {
        e.params = params;
        e.active = true;
        // Create animation with reference to the stable entity params
        auto anim = AnimT(e.params, std::forward<Args>(args)...);
        timeline.add(in_frames, anim.then([&e]() { e.active = false; }));
        return;
      }
    }
    // If we get here, no free slots. Drop the spawn (safe failure).
  }

  /**
   * @brief Applies all active transformations to a vector.
   */
  Vector transform(Vector v) const {
    for (const auto &e : entities) {
      if (e.active) {
        v = TransformFunc(v, e.params);
      }
    }
    return v;
  }
};

/**
 * @brief Generates ripples that warp the sphere
 */
template <int W, int CAPACITY>
using RippleTransformer =
    Transformer<W, RippleParams, Animation::Ripple, ripple_transform, CAPACITY>;

/**
 * @brief Performs Mobius warps that return to the identity.
 */
template <int W, int CAPACITY>
using MobiusWarpTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarp, mobius_transform,
                CAPACITY>;

/**
 * @brief Performs circular Mobius warps stay warped throughout sutiable for
 * repeating animations
 */
template <int W, int CAPACITY>
using MobiusWarpCircularTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarpCircular,
                mobius_transform, CAPACITY>;

/**
 * @brief Performs a changing Mobius warp using gnomonic projection.
 */
template <int W, int CAPACITY>
using MobiusWarpGnomonicTransformer =
    Transformer<W, MobiusParams, Animation::MobiusWarpEvolving,
                gnomonic_mobius_transform, CAPACITY>;

/**
 * @brief Applies 3D noise distortion to vectors.
 */
template <int W, int CAPACITY, int TIMELINE_CAP = 32>
using NoiseTransformer = Transformer<W, NoiseParams, Animation::Noise,
                                     noise_transform, CAPACITY, TIMELINE_CAP>;
