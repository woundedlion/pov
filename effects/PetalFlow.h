/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "../effects_engine.h"

// Helper linear easing
static float ease_linear(float t) { return t; }

template <int W, int H>
class PetalFlow : public Effect {
public:
  PetalFlow() :
    Effect(W, H),
    palette(
      { 0.029f, 0.029f, 0.029f },
      { 0.500f, 0.500f, 0.500f },
      { 0.461f, 0.461f, 0.461f },
      { 0.539f, 0.701f, 0.809f }
    ),
    orientation(),
    filters(
      Filter::World::Orient<W>(orientation),
      Filter::Screen::AntiAlias<W, H>()
    ),
    // Spawner: Manage ring creation based on gap accumulation
    spawner(1, [this](Canvas&) { this->check_spawn(); }, true)
  {
    persist_pixels = false;

    // Initialize Rings
    for (int i = 0; i < MAX_RINGS; ++i) {
      rings[i].active = false;
    }
    
    // Initialize Timeline
    init_timeline();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    
    // Step timeline ONCE per frame
    // Speed control is handled by the magnitude of position updates per frame
    timeline.step(canvas);
  }

private:
  static constexpr int MAX_RINGS = 64;
  static constexpr float START_RHO = -3.75f;
  static constexpr float END_RHO = 3.75f;
  static constexpr float SPACING = 0.3f; // Spacing in rho units

  struct Ring {
    float rho;
    bool active;
  };
  
  Ring rings[MAX_RINGS];
  float gap_accumulator = 0.0f;

  // JS Defaults
  float twist_factor = 2.15f;
  float speed = 8.0f; 
  float alpha = 0.2f;

  ProceduralPalette palette;
  Orientation<W> orientation;
  
  Pipeline<W, H,
    Filter::World::Orient<W>,
    Filter::Screen::AntiAlias<W, H>
  > filters;
  
  Timeline<W> timeline;
  Animation::PeriodicTimer spawner;

  void init_timeline() {
      // Animation for Orientation
      timeline.add(0, Animation::Rotation<W>(orientation, UP, PI_F / 4.0f, 160, ease_mid, true));
      
      // Animation for Twist Mutation
      timeline.add(0, Animation::Mutation(twist_factor, sin_wave(2.0f, 2.5f, 1.0f, 0.0f), 160, ease_mid, true));
      
      // Add Spawner
      gap_accumulator = 0.0f;
      timeline.add(0, spawner);
      
      // Pre-warm: Spawn rings along the path
      // Iterate backwards so the 'newest' (closest to START_RHO) corresponds to recent spawn
      for (float r = END_RHO - 0.01f; r > START_RHO; r -= SPACING) {
          spawn_ring_at_pos(r);
      }
  }

  void check_spawn() {
      // Calculate distance moved this frame based on 16fps target
      // Target rate: speed * 0.015 units/sec
      // Frame rate: 16 fps
      // Per frame: (speed * 0.015) / 16.0 = speed * 0.0009375
      
      float move_dist = speed * 0.0009375f;
      
      gap_accumulator += move_dist;
      
      while (gap_accumulator >= SPACING) {
          gap_accumulator -= SPACING;
          // Spawn at START_RHO + overshoot
          spawn_ring_at_pos(START_RHO + gap_accumulator);
      }
  }
  
  void spawn_ring_at_pos(float initial_rho) {
      for (int i = 0; i < MAX_RINGS; ++i) {
          if (!rings[i].active) {
              Ring& r = rings[i];
              r.active = true;
              r.rho = initial_rho;

              // Use Sprite with effectively infinite duration (1M frames)
              auto sprite_anim = Animation::Sprite(
                  [this, &r](Canvas& c, float) {
                      this->draw_ring(c, r); 
                  },
                  1000000 
              );
              
              timeline.add(0, std::move(sprite_anim));
              return;
          }
      }
  }

  void draw_ring(Canvas& canvas, Ring& ring) {
      if (!ring.active) return;
      
      // Update Position: 
      // Scale down to match JS per-second rate in a 16fps frame context
      // speed 8.0 * 0.0009375 = 0.0075 units/frame
      float move_dist = speed * 0.0009375f;
      ring.rho += move_dist;

      // Lifecycle Check
      if (ring.rho > END_RHO) {
          ring.active = false;
          return;
      }

      // Calculate Opacity (Distance Based)
      // JS: dist = abs(effectiveLogR); if (dist > 2.5) opacity = ...
      float dist = std::abs(ring.rho);
      float opacity = 1.0f;
      if (dist > 2.5f) {
          opacity = std::max(0.0f, 1.0f - (dist - 2.5f) / 1.0f);
      }
      
      // Combined Opacity
      float effective_opacity = opacity * alpha;
      if (effective_opacity <= 0.01f) return;

      const int num_samples = W;
      const float step = 2.0f * PI_F / num_samples;
      
      // Color Logic
      // Index implies cyclic color based on spatial position.
      // We use negative index to maintain color flow direction matching travel
      float color_idx = -(ring.rho / SPACING); 
      float hue = wrap(color_idx * 0.13f, 1.0f);
      Color4 base_col = palette.get(hue).color;
      base_col.alpha = effective_opacity; 

      float twist_angle = (ring.rho / SPACING) * twist_factor;
      
      auto get_shift = [](float t) -> float {
          return 0.6f * std::abs(sinf(3.0f * PI_F * t));
      };

      Fragments fragments;
      fragments.reserve(num_samples + 1);
      
      for (int i = 0; i < num_samples; ++i) {
          float t_norm = static_cast<float>(i) / num_samples;
          float theta = i * step;

          float rho_val = ring.rho + get_shift(t_norm);
          float final_theta = theta + twist_angle;

          float R = expf(rho_val);
          
          float r2 = R * R;
          float denom = 1.0f + r2;
          float x = 2.0f * R * cosf(final_theta) / denom;
          float y = 2.0f * R * sinf(final_theta) / denom;
          float z = (r2 - 1.0f) / denom;

          Fragment f;
          f.pos = Vector(x, y, z);
          f.v0 = t_norm;
          f.age = 0;
          
          fragments.push_back(f);
      }
      
      if (!fragments.empty()) {
          Fragment f = fragments[0];
          f.v0 = 1.0f;
          fragments.push_back(f);
      }

      auto fragment_shader = [&](const Vector&, Fragment& f) {
          f.color = base_col;
      };

      Plot::rasterize<W, H>(filters, canvas, fragments, fragment_shader, true);
  }
};