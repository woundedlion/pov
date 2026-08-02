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
#include <vector>

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
  std::vector<Pixel> reference(static_cast<size_t>(WIDTH) * HEIGHT);
  for (int frame = 0; frame < REPLAY_FRAMES; ++frame) {
    mindsplatter_replay::FrameStats stats;
    mindsplatter_replay::FrameStats host_stats;
    {
      Canvas canvas(effect);
      const ClipRegion clip = canvas.clip();
      WhiteBox::draw_particles_replay_reference(effect, canvas);
      const size_t reference_count = mindsplatter_replay::capture_frame<WIDTH>(
          canvas.data(), clip, reference.data(), reference.size());
      mindsplatter_replay::clear_frame(canvas, clip);
      WhiteBox::draw_particles_candidate(effect, canvas);
      stats = mindsplatter_replay::compare_frame_reference<WIDTH>(
          canvas, clip, reference.data(), reference_count);
      host_stats = mindsplatter_replay::compare_frame<WIDTH>(
          canvas, clip, corpus.framebuffer, corpus.framebuffer_entries);
    }
    effect.advance_display();
    std::printf(
        "replay frame=%d candidate=%llu reference=%llu changed=%u "
        "channels=%u max=%u abs=%llu\n",
        frame + 1, static_cast<unsigned long long>(stats.framebuffer_hash),
        static_cast<unsigned long long>(stats.expected_hash),
        stats.changed_pixels, stats.changed_channels, stats.max_channel_error,
        static_cast<unsigned long long>(stats.total_absolute_error));
    std::printf("replay host frame=%d candidate=%llu host=%llu changed=%u\n",
                frame + 1,
                static_cast<unsigned long long>(host_stats.framebuffer_hash),
                static_cast<unsigned long long>(host_stats.expected_hash),
                host_stats.changed_pixels);
    exact &= stats.changed_pixels == 0 && host_stats.changed_pixels == 0 &&
             host_stats.expected_hash == corpus.framebuffer_hash;
  }
  return exact ? 0 : 1;
}
