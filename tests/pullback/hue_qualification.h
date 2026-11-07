/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file @brief Cached hue modulation qualification over cycling palettes. */

/** @brief Cached hue modulation meets channel bounds across reachable palette intervals. */
inline void test_pullback_cached_hue_qualification() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Layer = Color::Modulate<DisplacementLayerProvider>;
  for (uint32_t offset : {0u, 17u, 83u})
    for (float chroma : {.314f, .718f, .788f, 1.f})
      for (auto mode :
           {Cycle::TRIADIC, Cycle::COMPLEMENTARY, Cycle::ANALOGOUS}) {
        HS_CONTEXT("cycle-offset", offset);
        HS_CONTEXT("chroma-milli", int(chroma * 1000));
        HS_CONTEXT("harmony", int(mode));
        std::array<uint8_t, GeneratedPaletteBank::required_arena_bytes() +
                                2 * BakedPalette::required_arena_bytes()>
            bytes{};
        Arena arena(bytes.data(), bytes.size());
        GeneratedPaletteBank bank;
        bank.init(arena, chroma, math::ease_in_out_sin);
        auto snapshot = bank.snapshot();
        for (size_t cycle = 0; cycle < snapshot.cycles.size(); ++cycle) {
          snapshot.cycles[cycle].next_sequence = 2 + offset;
          snapshot.hues[cycle] = (1 + offset) * GeneratedPaletteBank::HUE_STEP;
        }
        HS_EXPECT_TRUE(GeneratedPaletteBank::valid_snapshot(snapshot));
        bank.restore_snapshot(snapshot);

        BakedPaletteStorage from, to;
        from.clone_from(bank.palette(mode), arena);
        to.clone_from(bank.palette(mode), arena);
        std::array<PaletteLabLut, 2> caches;
        std::array<PaletteRgbHueLut, 2> rgb;
        PaletteFrame frame;
        frame.palette = &to.view();
        frame.temporal.palette = &from.view();
        frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.f};
        const auto prepared_layer = Layer::prepare(frame);
        int elapsed = 0, last = -1;
        uint32_t maximum = 0;
        uint64_t sum = 0, channels = 0;
        size_t samples = 0;
        for (int target :
             {0, 1, 150, 300, 600, 899, 900, 901, 902, 1200, 1800}) {
          while (elapsed < target) {
            from.rebake_copy(bank.palette(mode));
            bank.step(mode);
            ++elapsed;
          }
          to.rebake_copy(bank.palette(mode));
          const int slot = (last + 1) % 2;
          auto &cache = caches[slot];
          cache.bake(from.view(), to.view());
          rgb[slot].bake(cache, last < 0 ? nullptr : &caches[last],
                         last < 0 ? nullptr : &rgb[last]);
          frame.lab = &cache;
          HS_EXPECT_TRUE(cache.rgb_nodes && cache.rgb_masks);
          for (float fraction : {.17f, .5f, .83f}) {
            frame.temporal.fraction = fraction;
            const auto prepared = Policy::prepare(frame);
            const int LOW_ROW =
                std::min(int(fraction * (PaletteLabLut::TIMES - 1)),
                         int(PaletteLabLut::TIMES) - 2);
            const auto *const LOW_MASK =
                &rgb[slot].masks[cache.rgb_midpoint ? cache.rgb_midpoint_row
                                                    : LOW_ROW];
            const auto *const HIGH_MASK =
                &rgb[slot].masks[cache.rgb_midpoint ? cache.rgb_midpoint_row
                                                    : LOW_ROW + 1];
            HS_EXPECT_TRUE(prepared.lab.rgb_mask_low == LOW_MASK);
            HS_EXPECT_TRUE(prepared.lab.rgb_mask_high ==
                           (LOW_MASK == HIGH_MASK ? nullptr : HIGH_MASK));
            for (int value = 0; value < 256; ++value) {
              const float coordinate = (value + .5f) / 256;
              for (int hue = 0; hue < 128; ++hue) {
                const float scalar = 2.f * (hue + .5f) / 128 - 1.f;
                const FieldSample field{coordinate, .75f, math::Y_AXIS, scalar};
                ColorSample cached;
                HS_EXPECT_TRUE(
                    Policy::apply_shared_base(field, frame, prepared, cached));
                const auto pair = Policy::apply_pair(field, frame, prepared);
                auto direct =
                    blend_color_samples(pair.first, pair.second, fraction);
                const auto actual = finish_color_sample(
                    Layer::apply(cached, frame, prepared_layer));
                const auto expected = finish_color_sample(
                    Layer::apply(direct, frame, prepared_layer));
                HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-7f);
                for (auto error :
                     {channel_gap(actual.color.r, expected.color.r),
                      channel_gap(actual.color.g, expected.color.g),
                      channel_gap(actual.color.b, expected.color.b)}) {
                  maximum = std::max(maximum, uint32_t(error));
                  sum += error;
                  ++channels;
                }
                ++samples;
              }
            }
          }
          last = slot;
        }
        const double mean = double(sum) / channels;
        std::printf(
            "  cached hue offset=%u chroma=%.3f mode=%d max=%u mean=%.9f channels=%llu samples=%zu\n",
            offset, double(chroma), int(mode), maximum, mean,
            (unsigned long long)channels, samples);
        HS_EXPECT_LE(float(maximum), Layer::METRICS[0].limit);
        HS_EXPECT_LE(mean, Layer::METRICS[1].limit);
      }
}
