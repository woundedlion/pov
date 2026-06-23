/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/engine.h"

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
  FLASHMEM RingShower() : Effect(W, H) {}

  /**
   * @brief Registers parameters, bakes per-slot palette LUTs, and arms the
   *        spawn timer.
   * @details Allocates each slot's palette LUT once up front; spawn_ring refills
   *          it in place (rebake, no allocation). 16 * 256 entries amortizes
   *          trivially.
   */
  void init() override {
    registerParam("Alpha", &params.alpha, 0.0f, 1.0f);

    // Allocate each slot's palette LUT once up front; spawn_ring refills it in
    // place.
    for (size_t i = 0; i < MAX_RINGS; ++i)
      rings[i].palette.bake(persistent_arena, make_palette());

    timeline.add(0, Animation::RandomTimer(
                        4, 48, [this](Canvas &) { this->spawn_ring(); }, true));
  }

  /**
   * @brief Reports whether the engine should clear to the background each frame.
   * @return false; rings accumulate over the prior frame's contents.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Advances and draws every live ring for the current frame.
   * @details Steps the spawn timer, then advances and draws each live ring
   *          directly from its slot. Radius growth, fade-in and lifetime are
   *          pure functions of `age` computed here rather than by per-ring
   *          Sprite/Transition animations capturing the slot: a slot is
   *          recyclable, so an animation outliving its reuse would draw whatever
   *          ring later lands in it. Driving everything from `age` makes the
   *          slot's lifetime explicit and self-contained (mirrors Thrusters).
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas); // drives the spawn timer only

    for (size_t i = 0; i < MAX_RINGS; ++i) {
      Ring &ring = rings[i];
      if (ring.age >= ring.life)
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
    static constexpr int FADE_IN_FRAMES = 4; /**< Frames spent fading in from 0 to full opacity; opacity then holds (no fade-out). */
    static constexpr float RADIUS_MAX = 2.0f; /**< Maximum radius; the ring grows from 0 to this over its whole life. */
    // Lifetime is drawn uniformly from [LIFE_MIN, LIFE_MIN + LIFE_SPAN) frames.
    // LIFE_MIN is strictly positive *by design*: it is the sole guarantee that
    // radius_at()'s `/ life` divisor is never zero (spawn_ring() is the only
    // writer of life). Keep LIFE_MIN >= 1 if these are ever retuned.
    static constexpr int LIFE_MIN = 8;   /**< Minimum lifetime in frames (>0: guards the radius_at divisor). */
    static constexpr int LIFE_SPAN = 72; /**< Width of the random lifetime range, in frames. */

    Vector normal; /**< Plane normal fixing the ring's orientation. */
    /**
     * @brief 256-entry palette LUT, allocated once in init() and rebaked in
     *        place each spawn.
     * @details The per-fragment lookup is then a cheap LUT read rather than a
     *          full GenerativePalette OKLCH evaluation (the palette is immutable
     *          after spawn).
     */
    BakedPalette palette;
    int life = 0; /**< Total visible frames. */
    int age = 0;  /**< Frames elapsed; age >= life marks the slot free. */

    /**
     * @brief Constructs a free slot with a random orientation.
     */
    Ring() : normal(random_vector()) {}

    /**
     * @brief Eased radius for the frame currently being drawn.
     * @return Radius in world units, in [0, RADIUS_MAX].
     * @details Uses age + 1 (not age) so the first draw renders one eased step
     *          in rather than radius 0, and the ring reaches RADIUS_MAX on its
     *          final visible frame (age + 1 == life).
     */
    float radius_at() const {
      // `/ life` is divisor-safe because spawn_ring() always sets life >= LIFE_MIN
      // (> 0); a slot is never drawn before being spawned (age >= life marks it free).
      float t = hs::clamp(static_cast<float>(age + 1) / life, 0.0f, 1.0f);
      return RADIUS_MAX * ease_linear(t);
    }

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
   * @brief Reinitializes the first free ring slot with a new orientation, life,
   *        and palette.
   * @details Scans for a slot whose age has reached its life; if none is free,
   *          the spawn is silently dropped. life is drawn uniformly in
   *          [LIFE_MIN, LIFE_MIN + LIFE_SPAN) = [8, 80) frames. With MAX_RINGS
   *          slots, lives up to 80 frames, and a
   *          4-48 frame spawn cadence the pool can briefly fill, so a dropped
   *          spawn is an EXPECTED transient (bounded soft handling — a missed
   *          ring is invisible against the shower), NOT an invariant violation:
   *          it intentionally drops rather than trapping. Do not convert this to
   *          an HS_CHECK.
   */
  void spawn_ring() {
    for (size_t i = 0; i < MAX_RINGS; ++i) {
      if (rings[i].age >= rings[i].life) { // free slot
        Ring &ring = rings[i];
        ring.normal = random_vector();
        ring.life = static_cast<int>(hs::rand_f() * Ring::LIFE_SPAN +
                                     Ring::LIFE_MIN);
        ring.age = 0;
        ring.palette.rebake(make_palette());
        return;
      }
    }
  }

  /**
   * @brief Draws one ring slot at its current eased radius and opacity.
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
    const Vector z = X_AXIS;
    auto fragment_shader = [&](const Vector &v, Fragment &f) {
      f.color = ring.palette.get(angle_between(z, v) / PI_F);
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

#include "core/effect_registry.h"
REGISTER_EFFECT(RingShower)
