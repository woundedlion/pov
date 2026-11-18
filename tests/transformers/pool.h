/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Pool test cases. */

#include "tests/transformers/fixture.h"

namespace hs_test::transformers_tests {

// ============================================================================
// Transformer<> active slots — composition, live edits and pinned repeats
// ============================================================================

/**
 * @brief Verifies the compact active-slot list: a spawned entity is applied by
 *        transform(), and two active entities compose while staying on the unit
 *        sphere.
 * @details A spawned Noise entity displaces immediately, with no timeline
 *          stepping.
 */
inline void test_transformer_spawn_applies_and_composes() {
  Timeline tl;
  NoiseTransformer<4> nt(tl);
  nt.init_storage(persistent_arena);
  nt.template_params.amplitude = 0.6f;
  nt.template_params.scale = 4.0f;
  nt.template_params.time = 1.0f;
  nt.template_params.sync();

  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 0, 1),
                                  math::Vector(0.4f, 0.6f, 0.7f).normalized()};

  for (const math::Vector &v : samples) {
    math::Vector r = nt.transform(v);
    HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
    HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
    HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
  }

  HS_EXPECT_TRUE(nt.spawn_pinned(0) != nullptr);
  constexpr size_t N = sizeof(samples) / sizeof(samples[0]);
  math::Vector single[N];
  float total_moved = 0.0f;
  for (size_t i = 0; i < N; ++i) {
    const math::Vector &v = samples[i];
    single[i] = nt.transform(v);
    HS_EXPECT_TRUE(finite_vec(single[i]));
    HS_EXPECT_NEAR(single[i].length(), 1.0f, 1e-3f);
    total_moved += displacement(single[i], v);
  }
  HS_EXPECT_GT(total_moved, 1e-2f);

  HS_EXPECT_TRUE(nt.spawn_pinned(0) != nullptr);
  float total_divergence = 0.0f;
  for (size_t i = 0; i < N; ++i) {
    const math::Vector &v = samples[i];
    math::Vector r = nt.transform(v);
    HS_EXPECT_TRUE(finite_vec(r));
    HS_EXPECT_NEAR(r.length(), 1.0f, 1e-3f);
    total_divergence += displacement(r, single[i]);
  }
  HS_EXPECT_GT(total_divergence, 1e-3f);
}

/** @brief Spawning noise synchronizes the template frequency before any frame preparation. */
inline void test_noise_transformer_spawn_synchronizes_frequency() {
  Timeline timeline;
  NoiseTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  transformer.template_params.frequency = .4f;
  transformer.template_params.set_seed(42);
  HS_EXPECT_TRUE(transformer.spawn_pinned(0) != nullptr);
  const auto &active = transformer.active_params(0);
  Animation::NoiseParams expected;
  expected.frequency = .4f;
  expected.set_seed(42);
  expected.sync();
  for (float coordinate : {.13f, .3f, .7f, 1.2f})
    HS_EXPECT_EQ(active.noise.GetNoise(coordinate, .7f, 1.2f),
                 expected.noise.GetNoise(coordinate, .7f, 1.2f));
}

/** @brief Pins transformer live edit preserves instance phase and seed. */
inline void test_transformer_live_edit_preserves_instance_phase_and_seed() {
  Timeline timeline;
  NoiseTransformer<4> transformer(timeline);
  transformer.init_storage(persistent_arena);
  transformer.template_params.time = 3.0f;
  transformer.template_params.set_seed(42);
  HS_EXPECT_TRUE(transformer.spawn_pinned(0) != nullptr);
  transformer.template_params.amplitude = 0.7f;
  transformer.template_params.speed = 0.2f;
  transformer.template_params.frequency = 0.4f;
  transformer.template_params.scale = 6.0f;
  transformer.template_params.time = 99.0f;
  transformer.template_params.set_seed(99);
  transformer.prepare_frame();
  const auto &active = transformer.active_params(0);
  HS_EXPECT_EQ(active.amplitude, 0.7f);
  HS_EXPECT_EQ(active.speed, 0.2f);
  HS_EXPECT_EQ(active.frequency, 0.4f);
  HS_EXPECT_EQ(active.scale, 6.0f);
  HS_EXPECT_EQ(active.time, 3.0f);
  HS_EXPECT_EQ(active.seed, 42);
  Animation::NoiseParams expected;
  expected.frequency = 0.4f;
  expected.set_seed(42);
  expected.sync();
  HS_EXPECT_EQ(active.noise.GetNoise(0.3f, 0.7f, 1.2f),
               expected.noise.GetNoise(0.3f, 0.7f, 1.2f));
  RippleTransformer<1> ripple(timeline);
  ripple.init_storage(persistent_arena);
  HS_EXPECT_TRUE(ripple.spawn(0, math::Y_AXIS, 0.4f, 24) != nullptr);
  StubEffect effect(8, 8);
  Canvas canvas(effect);
  for (int frame = 0; frame < 5; ++frame)
    timeline.step(canvas);
  const auto before = ripple.active_params(0);
  ripple.template_params.amplitude = 0.9f;
  ripple.template_params.decay = 2.0f;
  ripple.template_params.thickness = 0.3f;
  ripple.prepare_frame();
  const auto &refreshed = ripple.active_params(0);
  HS_EXPECT_EQ(refreshed.decay, 2.0f);
  HS_EXPECT_EQ(refreshed.thickness, 0.3f);
  HS_EXPECT_EQ(refreshed.amplitude, before.amplitude);
  HS_EXPECT_EQ(refreshed.phase, before.phase);
  HS_EXPECT_TRUE(vec_identical(refreshed.center, before.center));
  auto expected_ripple = before;
  expected_ripple.decay = 2.0f;
  expected_ripple.thickness = 0.3f;
  expected_ripple.sync();
  HS_EXPECT_EQ(refreshed.cos_near, expected_ripple.cos_near);
  HS_EXPECT_EQ(refreshed.cos_far, expected_ripple.cos_far);
}

/** @brief A pinned repeating warp retains its slot across cycles until cancellation. */
inline void test_transformer_pinned_repeat_retains_slot() {
  Timeline timeline;
  StubEffect effect(8, 8);
  Canvas canvas(effect);
  MobiusWarpTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  auto *animation = transformer.spawn_pinned(
      0, 0.4f, 4,
      Animation::MobiusWarpOptions{.repeat = true,
                                   .easing = math::ease_linear});
  HS_EXPECT_TRUE(animation != nullptr);
  if (!animation)
    return;
  for (int frame = 0; frame < 12; ++frame) {
    HS_CONTEXT("repeat frame", frame);
    timeline.step(canvas);
    HS_EXPECT_EQ(transformer.active_count(), 1);
    if (transformer.active_count() != 1)
      return;
    if (frame % 4 == 0)
      HS_EXPECT_GT(std::abs(transformer.active_params(0).b.im), 0.1f);
  }
  animation->cancel();
  timeline.step(canvas);
  HS_EXPECT_EQ(transformer.active_count(), 0);
  HS_EXPECT_TRUE(
      transformer.spawn_pinned(
          0, 0.4f, 4, Animation::MobiusWarpOptions{.repeat = true}) != nullptr);
}

/** @brief The circular Mobius alias applies the authored Mobius kernel. */
inline void test_mobius_circular_alias_kernel() {
  Timeline timeline;
  StubEffect effect(8, 8);
  Canvas canvas(effect);
  MobiusWarpCircularTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  HS_EXPECT_TRUE(transformer.spawn_pinned(0, 1.0f, 160) != nullptr);
  timeline.step(canvas);
  HS_EXPECT_EQ(transformer.active_count(), 1);
  const math::Vector input = math::Vector(0.3f, 0.5f, 0.7f).normalized();
  HS_EXPECT_TRUE(vec_identical(
      transformer.transform(input),
      math::mobius_transform(input, transformer.active_params(0))));
}

/** @brief The gnomonic Mobius alias applies the authored gnomonic kernel. */
inline void test_mobius_gnomonic_alias_kernel() {
  static_assert(std::is_same_v<
                MobiusWarpGnomonicTransformer<1>,
                Transformer<math::MobiusParams, Animation::MobiusWarpEvolving,
                            math::gnomonic_mobius_transform, 1>>);
  Timeline timeline;
  StubEffect effect(8, 8);
  Canvas canvas(effect);
  MobiusWarpGnomonicTransformer<1> transformer(timeline);
  transformer.init_storage(persistent_arena);
  HS_EXPECT_TRUE(transformer.spawn_pinned(0, 0.5f, 0.0f) != nullptr);
  timeline.step(canvas);
  HS_EXPECT_EQ(transformer.active_count(), 1);
  const math::Vector input = math::Vector(0.3f, 0.5f, 0.7f).normalized();
  const auto actual = transformer.transform(input);
  const auto expected =
      math::gnomonic_mobius_transform(input, transformer.active_params(0));
  const auto wrong =
      math::mobius_transform(input, transformer.active_params(0));
  HS_EXPECT_EQ(actual.x, expected.x);
  HS_EXPECT_EQ(actual.y, expected.y);
  HS_EXPECT_EQ(actual.z, expected.z);
  HS_EXPECT_FALSE(values_equal(wrong.x, expected.x) &&
                  values_equal(wrong.y, expected.y) &&
                  values_equal(wrong.z, expected.z));
}

// ============================================================================
// Transformer<> spawn_pausable() gates the whole event, start delay included
// ============================================================================

/**
 * @brief Verifies spawn_pausable() freezes the spawned event's pending start
 *        delay, not just the animation once it has begun.
 * @details An unpaused pool is the control. Its completion relocates the
 *          gated event, so the pause gate must survive compaction.
 */
inline void test_transformer_spawn_pausable_freezes_start_delay() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  constexpr int DELAY = 20;
  constexpr int DURATION = 4;
  const math::Vector CENTER(0, 1, 0);

  RippleTransformer<1> control(tl);
  control.init_storage(persistent_arena);
  RippleTransformer<1> gated(tl);
  gated.init_storage(persistent_arena);

  bool paused = false;
  HS_EXPECT_TRUE(control.spawn(DELAY, CENTER, 0.2f, DURATION) != nullptr);
  HS_EXPECT_TRUE(gated.spawn_pausable(&paused, DELAY, CENTER, 0.2f, DURATION) !=
                 nullptr);
  HS_EXPECT_EQ(gated.active_count(), 1);

  constexpr int ELAPSED = 5;
  for (int i = 0; i < ELAPSED; ++i)
    tl.step(cv);
  paused = true;
  for (int i = 0; i < DELAY + DURATION + 4; ++i)
    tl.step(cv);
  HS_EXPECT_EQ(control.active_count(), 0);
  HS_EXPECT_EQ(gated.active_count(), 1);

  paused = false;
  for (int i = 0; i < DELAY - ELAPSED - 1; ++i)
    tl.step(cv);
  HS_EXPECT_EQ(gated.active_count(), 1);
  HS_EXPECT_NEAR(gated.active_params(0).amplitude, 0.0f, 1e-6f);

  for (int i = 0; i < DURATION + 1; ++i)
    tl.step(cv);
  HS_EXPECT_EQ(gated.active_count(), 0);
}

// ============================================================================
// Transformer<> slot release on rejection and compaction
// ============================================================================

/** @brief A full timeline releases a rejected transformer's pool slot. */
inline void test_transformer_full_timeline_releases_slot() {
  Timeline tl;
  hs_test::StubEffect effect(8, 8);
  Canvas canvas(effect);
  RippleTransformer<1> transformer(tl);
  transformer.init_storage(persistent_arena);
  float value = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(value, 1.0f, 2, math::ease_linear));
  HS_EXPECT_TRUE(transformer.spawn(0, math::Y_AXIS, 0.2f, 4) == nullptr);
  HS_EXPECT_EQ(transformer.active_count(), 0);
  tl.step(canvas);
  tl.step(canvas);
  HS_EXPECT_TRUE(transformer.spawn(0, math::Y_AXIS, 0.2f, 4) != nullptr);
}

/**
 * @brief Verifies a non-pinned spawned transform's pool slot is reclaimed after
 *        the timeline relocates its event during compaction.
 * @details With CAPACITY 1, the re-spawn succeeds only if the relocated
 *          then() callback reclaimed the slot.
 */
inline void test_transformer_nonpinned_slot_reclaimed_after_compaction() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  RippleTransformer<1> rt(tl);
  rt.init_storage(persistent_arena);

  float dummy = 0.0f;
  tl.add(0, Animation::Transition(dummy, 1.0f, 2, math::ease_linear));

  Animation::Ripple *p = rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4);
  HS_EXPECT_TRUE(p != nullptr);
  // Slot occupied: with CAPACITY 1 a second spawn finds nothing free.
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4) == nullptr);

  // Step past the earlier event (relocates the ripple) and the ripple's own
  // completion, so the then() callback fires through the relocation.
  for (int i = 0; i < 12; ++i)
    tl.step(cv);

  Animation::Ripple *reclaimed = rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 4);
  HS_EXPECT_TRUE(reclaimed != nullptr);
}

/**
 * @brief A transformer applied with substitute params composes them in spawn
 *        order, as the live params compose.
 */
inline void test_transformer_substitute_params_compose_in_spawn_order() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);
  RippleTransformer<4> rt(tl);
  rt.init_storage(persistent_arena);
  rt.template_params.amplitude = 0.15f;
  rt.template_params.thickness = 0.7f;
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 1, 0), 0.05f, 30) != nullptr);
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(1, 0, 0.3f), 0.07f, 30) != nullptr);
  for (int i = 0; i < 6; ++i)
    tl.step(cv);
  rt.prepare_frame();
  Animation::RippleParams live[2] = {rt.active_params(0), rt.active_params(1)};
  for (int i = 0; i < 50; ++i) {
    const math::Vector v = math::fib_spiral(50, 0.5f, i);
    const math::Vector a = rt.transform(v);
    const math::Vector b = rt.transform(v, live);
    HS_EXPECT_EQ(a.x, b.x);
    HS_EXPECT_EQ(a.y, b.y);
    HS_EXPECT_EQ(a.z, b.z);
  }
}

// ============================================================================
// Transformer<> completion callbacks outliving their pool
// ============================================================================

/**
 * @brief Verifies removal of a canceled animation after its pool is destroyed
 *        leaves the replacement pool's storage alone.
 * @details A second pool is built at the freed address, so the stale
 *          callback must be rejected on identity, not address.
 */
inline void test_transformer_callback_after_pool_destroyed() {
  using Pool = RippleTransformer<1>;
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  alignas(Pool) static uint8_t pool_storage[sizeof(Pool)];
  Pool *first = new (pool_storage) Pool(tl);
  first->init_storage(persistent_arena);
  HS_EXPECT_TRUE(first->spawn(0, math::Vector(0, 1, 0), 0.2f, 4) != nullptr);
  first->~Pool();
  HS_EXPECT_TRUE(global_timeline_events[0].animation()->is_canceled());

  Pool *second = new (pool_storage) Pool(tl);
  second->init_storage(persistent_arena);
  HS_EXPECT_TRUE(second->spawn(0, math::Vector(1, 0, 0), 0.2f, 120) != nullptr);
  HS_EXPECT_EQ(second->active_count(), 1);

  // Removes the first pool's canceled ripple, firing its now-stale callback.
  for (int i = 0; i < 8; ++i)
    tl.step(cv);

  HS_EXPECT_EQ(second->active_count(), 1);
  HS_EXPECT_TRUE(second->spawn(0, math::Vector(0, 0, 1), 0.2f, 4) == nullptr);
  second->~Pool();
}

// ============================================================================
// Transformer<> slots are reclaimed by Timeline::clear()
// ============================================================================

/**
 * @brief Verifies Timeline::clear() frees the pool slots of the animations it
 *        destroys.
 * @details clear() destroys events without running their then() callbacks, so
 *          the pool's clear hook must reclaim the slots.
 */
inline void test_transformer_slots_released_by_timeline_clear() {
  Timeline tl;

  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);

  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(1, 0, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_EQ(rt.active_count(), 2);
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 0, 1), 0.2f, 60) == nullptr);

  tl.clear();
  HS_EXPECT_EQ(rt.active_count(), 0);
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 0, 1), 0.2f, 60) != nullptr);
  HS_EXPECT_TRUE(rt.spawn(0, math::Vector(0, 1, 0), 0.2f, 60) != nullptr);
  HS_EXPECT_EQ(rt.active_count(), 2);
}

/**
 * @brief Verifies a recycled freed slot composes its new warp in spawn order
 *        (last), not slot-index order (first).
 * @details Spawn order composes tag 1 then tag 2, so x reads 12; slot-index
 *          order would read 21.
 */
inline void test_transformer_recycled_slot_composes_in_spawn_order() {
  Timeline tl;
  hs_test::StubEffect fx(8, 8);
  Canvas cv(fx);

  Transformer<OrderParams, TagAnim, order_transform, 2> tr(tl);
  tr.init_storage(persistent_arena);

  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/9, /*duration=*/2) != nullptr); // slot 0
  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/1, /*duration=*/100) !=
                 nullptr); // slot 1

  // Step until the short warp completes and frees slot 0; the long one survives.
  for (int i = 0; i < 5; ++i)
    tl.step(cv);

  HS_EXPECT_TRUE(tr.spawn(0, /*order=*/2, /*duration=*/100) !=
                 nullptr); // recycles slot 0

  const math::Vector r = tr.transform(math::Vector(0, 0, 0));
  HS_EXPECT_NEAR(r.x, 12.0f, 1e-4f);
  const OrderParams subs[2] = {{3}, {4}};
  HS_EXPECT_NEAR(tr.transform(math::Vector(0, 0, 0), subs).x, 34.0f, 1e-4f);
}

// ============================================================================
// Transformer<> empty-pool identity
// ============================================================================

/**
 * @brief Verifies a manager with nothing spawned (all entity slots inactive)
 *        is the identity.
 */
inline void test_transformer_no_entities_is_identity() {
  Timeline tl;
  RippleTransformer<8> rt(tl);
  math::Vector v = math::Vector(0.6f, 0.2f, 0.77f).normalized();
  math::Vector r = rt.transform(v);
  HS_EXPECT_NEAR(r.x, v.x, 1e-6f);
  HS_EXPECT_NEAR(r.y, v.y, 1e-6f);
  HS_EXPECT_NEAR(r.z, v.z, 1e-6f);
}

} // namespace hs_test::transformers_tests
