/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Growing and warped ring white-box cases. */

/** @brief White-box accessor for the Thrusters warp curve and fire path. */
struct ThrustersWhiteBox {
  using FX = Thrusters<DEFAULT_W, DEFAULT_H>;
  using Slot = typename FX::ThrusterContext;

  /** @brief Checks the column-time main ring and spawned rings. */
  static void check_temporal_rings();

  /**
   * @brief Verifies warp_decay peaks at 0.7 at t=0 and relaxes to exactly 0 at t=1.
   */
  static void check_warp_endpoints() {
    HS_EXPECT_NEAR(FX::warp_decay(0.0f), 0.7f, 1e-6f); // peak at fire
    HS_EXPECT_NEAR(FX::warp_decay(1.0f), 0.0f, 1e-6f); // full relaxation by end
  }

  /**
   * @brief Verifies one fire spawns two thrusters a half turn apart on the ring.
   * @details The pair's polar angles differ (the warp shifts each radially), so
   *          the opposition is in azimuth about the ring axis, not antipodality.
   */
  static void check_fire_spawns_opposed_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.on_fire_thruster();
    HS_EXPECT_SIZE_OR_RETURN(effect.thrusters, 2);

    const math::Vector a =
        azimuth_dir(effect.thrusters[0].point, effect.ring_vec);
    const math::Vector b =
        azimuth_dir(effect.thrusters[1].point, effect.ring_vec);
    HS_EXPECT_NEAR(math::dot(a, b), -1.0f, 1e-4f);
  }

  /**
   * @brief Verifies a collapsed ring still fires.
   * @details theta_eq is radius * pi/2, so a sub-microradian ring puts both
   *          thrust points on ring_vec and drops the spin-axis cross product
   *          below EPS_NORMALIZE_SQ. Reaching the assertions is the check.
   */
  static void check_collapsed_ring_falls_back_to_a_spin_axis() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.params.radius = 1e-7f;
    effect.on_fire_thruster();
    HS_EXPECT_SIZE_OR_RETURN(effect.thrusters, 2);
    HS_EXPECT_NEAR(math::dot(effect.thrusters[0].point, effect.ring_vec), 1.0f,
                   1e-6f);
    HS_EXPECT_NEAR(math::dot(effect.thrusters[1].point, effect.ring_vec), 1.0f,
                   1e-6f);
  }

  /**
   * @brief Verifies the slot pool saturates and evicts whole fires in order.
   */
  static void check_fifo_evicts_the_oldest_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    constexpr size_t CAPACITY = effect.thrusters.capacity();
    constexpr int FIRES = CAPACITY / 2 + 1;
    math::Vector lead[FIRES];
    for (int fire = 0; fire < FIRES; ++fire) {
      effect.on_fire_thruster();
      const size_t spawned = 2u * static_cast<size_t>(fire + 1);
      const size_t expected = spawned < CAPACITY ? spawned : CAPACITY;
      HS_EXPECT_EQ(effect.thrusters.size(), expected);
      // The fire's first spawn: back() is its opposed twin.
      lead[fire] = effect.thrusters[effect.thrusters.size() - 2].point;
    }
    // The last fire pushed two slots into a full pool, retiring fire 0's pair
    // and leaving fire 1's first slot at the front.
    HS_EXPECT_NEAR(math::dot(effect.thrusters.front().point, lead[1]), 1.0f,
                   1e-6f);
  }

  /**
   * @brief Verifies draw_frame retires the expired prefix and ages the rest.
   */
  static void check_expired_slots_retire_by_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.on_fire_thruster();
    effect.on_fire_thruster();
    HS_EXPECT_SIZE_OR_RETURN(effect.thrusters, 4);

    // The first fire's pair is spent; the second's has one frame left.
    effect.thrusters[0].age = Slot::LIFE;
    effect.thrusters[1].age = Slot::LIFE;
    effect.thrusters[2].age = Slot::LIFE - 1;
    effect.thrusters[3].age = Slot::LIFE - 1;

    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_SIZE_OR_RETURN(effect.thrusters, 2);
    HS_EXPECT_EQ(effect.thrusters.front().age, Slot::LIFE);

    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_TRUE(effect.thrusters.is_empty());
  }

private:
  /** @brief Unit azimuth of @p p about @p axis, with the radial warp removed. */
  static math::Vector azimuth_dir(const math::Vector &p,
                                  const math::Vector &axis) {
    return (p - axis * math::dot(p, axis)).normalized();
  }
};

/** @brief White-box accessor for RingShower's radius easing endpoints. */
struct RingShowerWhiteBox {
  /** @brief Checks growth, fade, first draws and recycled slots at column time. */
  static void check_temporal_rings();

  /**
   * @brief Verifies radius_at hits RADIUS_MAX on the final frame and is non-zero
   *        on the first.
   */
  static void check_radius_endpoints() {
    using RS = RingShower<DEFAULT_W, DEFAULT_H>;
    typename RS::Ring ring;
    ring.life = 50;
    ring.age = ring.life - 1; // final visible frame: age+1 == life -> t == 1
    HS_EXPECT_NEAR(ring.radius_at(), RS::Ring::RADIUS_MAX, 1e-5f);
    ring.age = 0; // first frame: one eased step in, not radius 0
    HS_EXPECT_GT(ring.radius_at(), 0.0f);
    HS_EXPECT_LT(ring.radius_at(), RS::Ring::RADIUS_MAX);
  }
};
