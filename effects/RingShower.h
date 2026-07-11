/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/engine/engine.h"

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// pins Ring::radius_at's age+1 endpoint convention (reaches RADIUS_MAX on the
// final visible frame); the smoke harness only proves the rings render.
namespace hs_test {
namespace effects_tests {
struct RingShowerWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Effect that continuously spawns randomly-oriented rings which grow and
 *        fade in over a short lifetime, then recycle their fixed slot.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A RandomTimer drives spawning; each live ring is advanced and drawn
 *          purely from its `age` (see draw_frame).
 */
template <int W, int H> class RingShower : public Effect {
public:
  /**
   * @brief Constructs the effect at the templated canvas dimensions.
   */
  HS_COLD_MEMBER RingShower() : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}) {}

  /**
   * @brief Registers parameters, bakes per-slot palette LUTs, and arms the
   *        spawn timer.
   * @details Allocates each slot's palette LUT once up front; spawn_ring refills
   *          it in place (rebake, no allocation). 16 * 256 entries amortizes
   *          trivially.
   */
  void init() override {
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    for (size_t i = 0; i < MAX_RINGS; ++i)
      rings[i].palette.bake(persistent_arena, dot_keyed(make_palette()));

    timeline.add(0, Animation::RandomTimer(
                        4, 48, [this](Canvas &) { this->spawn_ring(); }, true));
  }

  /**
   * @brief Advances and draws every live ring for the current frame.
   * @details Steps the spawn timer, then advances and draws each live ring from
   *          its slot. Radius, fade-in, and lifetime are pure functions of `age`,
   *          not per-ring animations capturing the slot: a slot is recyclable, so
   *          a stale animation would draw whatever ring later lands in it.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas); // drives the spawn timer only

    for (size_t i = 0; i < MAX_RINGS; ++i) {
      Ring &ring = rings[i];
      if (ring.expired())
        continue; // free slot
      draw_ring(canvas, ring.opacity_at(), i);
      ++ring.age;
    }
  }

private:
  // Test seam: reaches the private Ring type and its radius_at endpoint.
  friend struct ::hs_test::effects_tests::RingShowerWhiteBox;

  /**
   * @brief One ring slot: orientation, baked palette, and age/life timing.
   * @details A slot is free when age >= life and is reused by spawn_ring without
   *          reallocating.
   */
  struct Ring {
    static constexpr int FADE_IN_FRAMES = 4; /**< Fade-in span; opacity uses age+1 so the first draw starts one linear step in (0.25 at 4) rather than 0, then holds at full (no fade-out). */
    static constexpr float RADIUS_MAX = 2.0f; /**< Maximum radius; the ring grows from 0 to this over its whole life. */
    static constexpr int LIFE_MIN = 8;   /**< Minimum lifetime in frames (>0: guards the radius_at divisor). */
    static constexpr int LIFE_SPAN = 72; /**< Width of the random lifetime range, in frames. */

    Vector normal; /**< Plane normal fixing the ring's orientation. */
    /**
     * @brief 256-entry palette LUT, allocated once in init() and rebaked in
     *        place each spawn.
     * @details Per-fragment lookup is then a LUT read, not a GenerativePalette
     *          OKLCH evaluation (the palette is immutable after spawn).
     */
    BakedPalette palette;
    int age = 0;  /**< Frames elapsed since (re)spawn. */
    int life = 0; /**< Total visible frames; the slot is free once age >= life. */

    /**
     * @brief Constructs a free slot with a random orientation.
     */
    Ring() : normal(random_vector()) {}

    /**
     * @brief Linear radius for the frame currently being drawn.
     * @return Radius in world units, in [0, RADIUS_MAX].
     * @details Uses age + 1 (not age) so the first draw renders one linear step
     *          in rather than radius 0, and the ring reaches RADIUS_MAX on its
     *          final visible frame (age + 1 == life).
     */
    float radius_at() const {
      HS_CHECK(life > 0, "RingShower::Ring::radius_at: life must be > 0");
      float t = hs::clamp(static_cast<float>(age + 1) / life, 0.0f, 1.0f);
      return RADIUS_MAX * t;
    }

    /**
     * @brief Whether the slot has outlived its life and is recyclable/free.
     */
    bool expired() const { return age >= life; }

    /**
     * @brief Opacity multiplier for the frame currently being drawn.
     * @return Opacity in [0, 1]: eased fade-in over the first FADE_IN_FRAMES
     *         frames, then a held 1.0; no fade-out.
     */
    float opacity_at() const {
      if (age + 1 < FADE_IN_FRAMES)
        return ease_linear(static_cast<float>(age + 1) / FADE_IN_FRAMES);
      return 1.0f;
    }
  };

  /**
   * @brief Builds a fresh random palette for a spawning ring.
   * @return A circular, analogous, flat-brightness GenerativePalette.
   * @details Each construction reseeds, so every ring draws from its own
   *          distinct palette.
   */
  static GenerativePalette make_palette() {
    return GenerativePalette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS,
                             BrightnessProfile::FLAT);
  }

  /**
   * @brief Bake-time adapter that maps the LUT coordinate from the cos domain
   *        into the source's angle parameter.
   * @details Baking through this folds the d -> acos(d)/PI radial mapping into
   *          the bake (256 acos per spawn, not one per fragment): the fragment
   *          lookup keys the LUT by the raw dot product via dot_key(d). dot_key
   *          inverts DotKeyed: u -> d = 1 - 2u, get(u) returns the source at
   *          acos(d)/PI.
   */
  template <typename Source> struct DotKeyed {
    const Source &source;
    Color4 get(float u) const {
      float d = hs::clamp(1.0f - 2.0f * u, -1.0f, 1.0f);
      return source.get(fast_acos(d) / PI_F);
    }
  };
  template <typename Source>
  static DotKeyed<Source> dot_keyed(const Source &source) {
    return DotKeyed<Source>{source};
  }
  /// LUT coordinate for cos-value d = dot(X, v); inverse of DotKeyed's mapping.
  static float dot_key(float d) {
    return (1.0f - hs::clamp(d, -1.0f, 1.0f)) * 0.5f;
  }

  /**
   * @brief Reinitializes the first free ring slot with a new orientation, life,
   *        and palette.
   * @details Scans for a slot whose age has reached its life; if none is free,
   *          the spawn is silently dropped. The pool can briefly fill, so a
   *          dropped spawn is an EXPECTED transient (a missed ring is invisible
   *          against the shower), NOT an invariant violation — do not convert
   *          this to an HS_CHECK.
   */
  HS_COLD_MEMBER void spawn_ring() {
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].expired()) { // free slot
        Ring &ring = rings[i];
        ring.normal = random_vector();
        ring.life = static_cast<int>(hs::rand_f() * Ring::LIFE_SPAN +
                                     Ring::LIFE_MIN);
        ring.age = 0;
        ring.palette.rebake(dot_keyed(make_palette()));
        return;
      }
    }
  }

  /**
   * @brief Draws one ring slot at its current radius and opacity.
   * @param canvas Target canvas for this frame.
   * @param opacity Per-frame opacity multiplier in [0, 1] from opacity_at().
   * @param index Slot index into the rings array.
   * @details Rings are never re-oriented after spawn, so the basis comes
   *          straight from the ring's own normal (identity rotation) and the
   *          palette's radial axis is the fixed X axis — both constant across
   *          the ring's fragments.
   */
  void draw_ring(Canvas &canvas, float opacity, size_t index) {
    Ring &ring = rings[index];
    Basis basis = make_basis(Quaternion(), ring.normal);
    // v is unit (the rasterizer renormalizes every shaded position), so
    // dot(X, v) is just v.x; the palette is baked in this cos domain (dot_keyed),
    // folding the acos radial mapping into the bake.
    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      f.color = ring.palette.get(dot_key(v.x));
      f.color.alpha *= opacity * params.alpha;
    };
    Plot::Ring::draw<W, H>(filters, canvas, basis, ring.radius_at(),
                           fragment_shader);
  }

  static constexpr size_t MAX_RINGS = 16; /**< Fixed pool of ring slots. */
  Ring rings[MAX_RINGS];                  /**< The recyclable ring slots. */
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters; /**< Screen-space anti-alias pipeline. */
  Timeline timeline;    /**< Drives the spawn timer. */

  /**
   * @brief User-tunable parameters for the effect.
   */
  struct Params {
    float alpha = 0.2f; /**< Global opacity scale in [0, 1]. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(RingShower)
