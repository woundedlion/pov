/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "geometry.h"
#include "animation.h"
#include "static_circular_buffer.h"

// Forward declare to avoid circular dependency or missing include issues
class Canvas;

template <int W, int CAPACITY = 32> class RippleGenerator {
public:
  struct RippleEntity {
    RippleParams params;
    bool active = false;
  };

  RippleParams params; // Defaults linked to GUI
  StaticCircularBuffer<RippleEntity, CAPACITY> ripples;
  Timeline<W> &timeline;

  RippleGenerator(Timeline<W> &tl) : timeline(tl) {
    // Default params
    params.amplitude = 1.0f;
    params.thickness = 0.4f;
    params.decay = 0.5f;
    params.phase = 0.0f;
  }

  void set_amplitude(float v) { params.amplitude = v; }
  void set_thickness(float v) { params.thickness = v; }
  void set_decay(float v) { params.decay = v; }

  // Spawn a ripple at origin with a specific duration (frames to reach
  // antipode)
  void ripple(const Vector &origin, int duration = 180) {
    RippleEntity &r = ripples.emplace_back();
    r.params = params; // Copy defaults
    r.params.center = origin;
    r.active = true;

    // Speed: Travel PI radians (half circumference) in 'duration' frames
    float speed = PI_F / static_cast<float>(duration);

    timeline.add(
        0, Animation::Ripple(r.params, origin, speed, duration).then([&r]() {
          r.active = false;
        }));
  }

  // Apply all active ripples to a vector
  Vector apply(Vector v) const {
    for (const auto &r : ripples) {
      if (r.active && r.params.amplitude > 0.001f) {
        v = ripple_transform(v, r.params);
      }
    }
    return v;
  }

  // Get active ripple params for batch processing (optimization)
  // Returns a buffer of valid params to avoid checking 'active' repeatedly per
  // vertex
  StaticCircularBuffer<RippleParams, CAPACITY> get_active() const {
    StaticCircularBuffer<RippleParams, CAPACITY> active_list;
    for (const auto &r : ripples) {
      if (r.active && r.params.amplitude > 0.001f) {
        active_list.push_back(r.params);
      }
    }
    return active_list;
  }
};
