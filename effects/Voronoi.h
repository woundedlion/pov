/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

#include <cmath>

template <int W, int H> class Voronoi : public Effect {
public:
  struct Site {
    Vector pos;
    Vector axis;
    Color4 color;
    int id;
  };

  FLASHMEM Voronoi() : Effect(W, H) {}

  void init() override {
    registerParam("Num Sites", &params.num_sites, 1.0f,
                  static_cast<float>(MAX_SITES));
    registerParam("Speed", &params.speed, 0.0f, 100.0f);
    registerParam("Smoothness", &params.smoothness, 1.0f, 500.0f);
    registerParam("Border Thick", &params.borderThickness, 0.0f, 0.1f);

    // Allocate the buffer once at its maximum; seed_sites() fills the active
    // count and re-seeds when the slider changes (no further allocation).
    sites_buffer.bind(persistent_arena, MAX_SITES);
    seed_sites();
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);

    // Re-seed when the GUI changes the site count (integer change only, so
    // dragging within a bucket doesn't thrash).
    if (active_site_count() != current_num_sites)
      seed_sites();

    // 1. Animate Sites
    float s = logf(params.speed + 1.0f) * 0.005f;

    for (size_t i = 0; i < sites_buffer.size(); ++i) {
      auto &site = sites_buffer[i];
      // Rotate pos around axis
      Quaternion q = make_rotation(site.axis, s);
      site.pos = rotate(site.pos, q);
    }

    // 2. Render Pixels (Linear Search)
    // No sites (degenerate num_sites == 0): the search below would leave
    // bestIdx == -1 and index sites_buffer[-1] out of bounds for every pixel.
    // Nothing to render, so bail before the pixel loop.
    if (sites_buffer.size() == 0)
      return;

    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        Vector p = pixel_to_vector<W, H>(x, y);

        int bestIdx = -1;
        float bestDot = -2.0f;
        int secondBestIdx = -1;
        float secondBestDot = -2.0f;

        for (size_t i = 0; i < sites_buffer.size(); ++i) {
          float d = dot(p, sites_buffer[i].pos);
          if (d > bestDot) {
            secondBestDot = bestDot;
            secondBestIdx = bestIdx;
            bestDot = d;
            bestIdx = i;
          } else if (d > secondBestDot) {
            secondBestDot = d;
            secondBestIdx = i;
          }
        }

        float maxDot1 = bestDot;
        float maxDot2 = secondBestDot;
        const auto &bestSite = sites_buffer[bestIdx];

        Color4 c = bestSite.color;

        // Smoothing
        if (secondBestIdx != -1 && params.smoothness > 0.0f) {
          const auto &secSite = sites_buffer[secondBestIdx];
          float diff = maxDot1 - maxDot2;
          float factor = std::min(1.0f, diff * params.smoothness);
          factor = quintic_kernel(factor);
          float t = 0.5f + 0.5f * factor;

          // Mix (16-bit linear channels — do NOT truncate to uint8_t)
          c.color.r = static_cast<uint16_t>(
              secSite.color.color.r + (c.color.r - secSite.color.color.r) * t);
          c.color.g = static_cast<uint16_t>(
              secSite.color.color.g + (c.color.g - secSite.color.color.g) * t);
          c.color.b = static_cast<uint16_t>(
              secSite.color.color.b + (c.color.b - secSite.color.color.b) * t);
        }

        // Borders
        if (showBorders && secondBestIdx != -1) {
          float dist1 = acosf(std::min(1.0f, maxDot1));
          float dist2 = acosf(std::min(1.0f, maxDot2));
          if (dist2 - dist1 < params.borderThickness) {
            c = Color4(0, 0, 0, 0); // Black/Transparent
          }
        }

        canvas(x, y) = c.color;
      }
    }
  }

  struct Params {
    float num_sites = 200.0f; // live-tunable (GUI slider), see init/draw_frame
    float speed = 20.0f;
    float borderThickness = 0.0f;
    float smoothness = 100.0f;
  } params;

  bool showBorders = true;

  // The sites buffer is allocated once at MAX_SITES; the active count varies
  // with the "Num Sites" slider and is re-seeded (clear + refill, no realloc)
  // when it changes. current_num_sites tracks the count currently seeded.
  static constexpr int MAX_SITES = 400;
  int current_num_sites = 0;
  ArenaVector<Site> sites_buffer;

  int active_site_count() const {
    return std::clamp(static_cast<int>(params.num_sites), 1, MAX_SITES);
  }

  // (Re)seed the active sites for the current "Num Sites" slider value. The
  // buffer is allocated once at MAX_SITES (in init); here we only clear + refill
  // up to the active count, so changing the slider never re-allocates the arena.
  void seed_sites() {
    const int n = active_site_count();
    sites_buffer.clear();

    for (int i = 0; i < n; i++) {
      float goldenAngle = PI_F * (3.0f - sqrtf(5.0f));
      // Guard n == 1: the denominator would be 0 → NaN y (which then poisons
      // pos/axis). A single site sits at the pole (y = 1). Mirrors the same
      // guard applied to `t` below.
      int span = n > 1 ? n - 1 : 1;
      float y = 1.0f - (i / (float)span) * 2.0f;
      float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
      float theta = goldenAngle * i;

      float x = cosf(theta) * radius;
      float z = sinf(theta) * radius;

      Vector pos = Vector(x, y, z);

      float rx = hs::rand_f() * 2.0f - 1.0f;
      float ry = hs::rand_f() * 2.0f - 1.0f;
      float rz = hs::rand_f() * 2.0f - 1.0f;
      Vector axis = Vector(rx, ry, rz).normalized();

      float t = i / (float)(n > 1 ? n - 1 : 1);
      Color4 color = Palettes::richSunset.get(t);

      sites_buffer.push_back({pos, axis, color, i});
    }
    current_num_sites = n;
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Voronoi)
