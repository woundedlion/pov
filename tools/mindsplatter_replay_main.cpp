/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#include "core/engine/memory.h"
#include "effects/mindsplatter_replay_corpus.h"
#include "tools/mindsplatter_whitebox.h"
#include "tools/mindsplatter_replay_metrics.h"

#include <cstdio>
#include <span>

namespace {

constexpr int WIDTH = 288;
constexpr int HEIGHT = 144;
constexpr int REPLAY_FRAMES = 3;
using ReplayEffect = MindSplatter<WIDTH, HEIGHT>;
using WhiteBox = hs_test::effects_tests::MindSplatterWhiteBox;

} // namespace

int main() {
  const auto &corpus = mindsplatter_replay::DENSE_PRESET0_F136;
  hs::random().seed(1337u);
  GenerativePalette::reset_hue_seed(0);
  configure_arenas_default();
  ReplayEffect effect;
  effect.init();
  const uint16_t particle_count = WhiteBox::restore_render(
      effect, std::span(corpus.state, corpus.state_size));

  std::printf("replay corpus id=%s source=%s preset=%d traits=%u particles=%u "
              "bytes=%zu hash=%llu\n",
              corpus.id, corpus.source, corpus.preset, corpus.traits,
              particle_count, corpus.state_size,
              static_cast<unsigned long long>(corpus.corpus_hash));
  bool exact = true;
  for (int frame = 0; frame < REPLAY_FRAMES; ++frame) {
    mindsplatter_replay::FrameStats stats;
    WhiteBox::draw_particles_inspect(
        effect, [&](Canvas &canvas, const ClipRegion &clip) {
          stats = mindsplatter_replay::compare_frame<WIDTH>(
              canvas, clip, corpus.framebuffer, corpus.framebuffer_entries);
        });
    effect.advance_display();
    std::printf(
        "replay frame=%d hash=%llu expected=%llu changed=%u "
        "channels=%u max=%u abs=%llu\n",
        frame + 1, static_cast<unsigned long long>(stats.framebuffer_hash),
        static_cast<unsigned long long>(stats.expected_hash),
        stats.changed_pixels, stats.changed_channels, stats.max_channel_error,
        static_cast<unsigned long long>(stats.total_absolute_error));
    exact &= stats.changed_pixels == 0 &&
             stats.expected_hash == corpus.framebuffer_hash;
  }
  return exact ? 0 : 1;
}
