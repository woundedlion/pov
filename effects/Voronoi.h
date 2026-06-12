/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"
#include "core/spatial.h"

#include <cmath>
#include <span>

// Spherical Voronoi effect: scatters sites on the unit sphere, animates them,
// and shades each pixel by its nearest site with edge-sharpening and optional
// cell borders. W/H are the framebuffer dimensions.
template <int W, int H> class Voronoi : public Effect {
public:
  // One Voronoi seed: position on the unit sphere, the rotation axis it spins
  // about, and the cell fill color.
  struct Site {
    Vector pos;
    Vector axis;
    Color4 color;
  };

  FLASHMEM Voronoi() : Effect(W, H) {}

  // Configure arenas, register GUI params, allocate the sites buffer, and seed
  // the initial sites.
  void init() override {
    // Persistent holds the sites buffer; scratch_arena_a holds the per-frame
    // KD-tree (positions + nodes + build indices, ~15KB at MAX_SITES). Give it
    // comfortable headroom rather than rely on the 16KB default.
    configure_arenas(GLOBAL_ARENA_SIZE - 64 * 1024, 64 * 1024, 0);

    registerParam("Num Sites", &params.num_sites, 1.0f,
                  static_cast<float>(MAX_SITES));
    registerParam("Speed", &params.speed, 0.0f, 100.0f);
    // "Sharpness", not "Smoothness": a larger value saturates the blend factor
    // sooner, narrowing the border band — i.e. it sharpens cell edges.
    registerParam("Sharpness", &params.sharpness, 1.0f, 500.0f);
    registerParam("Border Thick", &params.borderThickness, 0.0f, 0.1f);

    // Allocate the buffer once at its maximum; seed_sites() fills the active
    // count and re-seeds when the slider changes (no further allocation).
    sites_buffer.bind(persistent_arena, MAX_SITES);
    seed_sites();
  }

  // Opaque effect: fills every pixel, so no background pass is needed.
  bool show_bg() const override { return false; }

  // Animate the sites, build a per-frame KD-tree, and shade each pixel by its
  // nearest site (with edge sharpening and optional borders).
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

    // 2. Render Pixels (KD-tree nearest + second-nearest)
    // No sites (degenerate num_sites == 0): nothing to render, so bail before
    // building the tree / pixel loop (an empty tree would yield no neighbors).
    if (sites_buffer.size() == 0)
      return;

    // Build a KD-tree over the (moving) site positions once per frame. On the
    // unit sphere nearest-by-Euclidean-distance is exactly nearest-by-max-dot
    // (|p−s|² = 2 − 2·p·s for unit p,s), so the k=2 query returns the exact
    // best/second-best site: O(N log N) build + O(log N) per pixel rather than
    // an O(W·H·N) per-pixel scan.
    ScratchScope _scope(scratch_arena_a);
    Vector *positions = static_cast<Vector *>(scratch_arena_a.allocate(
        sites_buffer.size() * sizeof(Vector), alignof(Vector)));
    for (size_t i = 0; i < sites_buffer.size(); ++i)
      positions[i] = sites_buffer[i].pos;
    KDTree tree(scratch_arena_a,
                std::span<const Vector>(positions, sites_buffer.size()));

    // Render via the shared full-screen shader path (clip-aware, instrumented).
    // SAMPLES=1: one sample at pixel center — the per-pixel KD-tree query is too
    // heavy to supersample.
    auto shader = [&](const Vector &p) -> Color4 {
      auto knn = tree.nearest(p, 2);
      // Tree is non-empty (sites_buffer.size() > 0), so there is always a
      // nearest; a second neighbor exists iff there are ≥ 2 sites.
      const auto &bestSite = sites_buffer[knn[0].original_index];
      bool hasSecond = knn.size() > 1;
      float maxDot1 = dot(p, knn[0].point);
      float maxDot2 = hasSecond ? dot(p, knn[1].point) : -2.0f;

      Color4 c = bestSite.color;

      // Border sharpening: a larger sharpness saturates `factor` for smaller
      // nearest/second-nearest gaps, shrinking the cross-cell blend band.
      if (hasSecond && params.sharpness > 0.0f) {
        const auto &secSite = sites_buffer[knn[1].original_index];
        float diff = maxDot1 - maxDot2;
        float factor = std::min(1.0f, diff * params.sharpness);
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

      // Borders — driven entirely by the "Border Thick" slider: a thickness of
      // 0 disables them (and skips the two acosf calls below), any positive
      // value paints the seam between the nearest two sites. knn[0] is the
      // nearest, so maxDot1 >= maxDot2 → dist1 <= dist2 and the cell gap is
      // non-negative.
      if (params.borderThickness > 0.0f && hasSecond) {
        float dist1 = acosf(std::min(1.0f, maxDot1));
        float dist2 = acosf(std::min(1.0f, maxDot2));
        if (dist2 - dist1 < params.borderThickness) {
          c = Color4(0, 0, 0, 0); // Black/Transparent
        }
      }

      return c;
    };
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  struct Params {
    float num_sites = 200.0f; // live-tunable (GUI slider), see init/draw_frame
    float speed = 20.0f;
    float borderThickness = 0.0f;
    float sharpness = 100.0f;
  } params;

  // The sites buffer is allocated once at MAX_SITES; the active count varies
  // with the "Num Sites" slider and is re-seeded (clear + refill, no realloc)
  // when it changes. current_num_sites tracks the count currently seeded.
  static constexpr int MAX_SITES = 400;
  int current_num_sites = 0;
  ArenaVector<Site> sites_buffer;

  // Slider value rounded to an integer and clamped to [1, MAX_SITES].
  int active_site_count() const {
    return std::clamp(static_cast<int>(params.num_sites), 1, MAX_SITES);
  }

  // (Re)seed the active sites for the current "Num Sites" slider value: clear +
  // refill up to the active count (no re-allocation). Sites are placed via the
  // Fibonacci-sphere distribution for an even spread, each with a random spin
  // axis and a palette color.
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

      sites_buffer.push_back({pos, axis, color});
    }
    current_num_sites = n;
  }
};

#include "core/effect_registry.h"
REGISTER_EFFECT(Voronoi)
