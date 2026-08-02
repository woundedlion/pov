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

/**
 * @brief Production-corpus visual bounds for single-pass sample phasing.
 * @details At most 18.1% of pixels and 17.3% of channels may differ. Mean
 * absolute error is capped at 73 channel counts and 61 luminance counts per
 * framebuffer pixel; added and dropped AA fringe pixels are capped at 0.2%
 * and 0.5% of the frame, with their luminance error capped at five counts per
 * framebuffer pixel and 9,000 for any one fringe pixel.
 */
struct VisualGate {
  static constexpr uint32_t PIXELS = WIDTH * HEIGHT;
  static constexpr uint32_t MAX_CHANGED_PIXELS = PIXELS * 181 / 1000;
  static constexpr uint32_t MAX_CHANGED_CHANNELS = PIXELS * 3 * 173 / 1000;
  static constexpr uint16_t MAX_CHANNEL_ERROR = 14000;
  static constexpr uint64_t MAX_TOTAL_ERROR = PIXELS * 3ull * 73;
  static constexpr uint64_t MAX_LUMINANCE_ERROR = PIXELS * 61ull;
  static constexpr int64_t MAX_LUMINANCE_BIAS = PIXELS * 2ll;
  static constexpr uint32_t MAX_ADDED_PIXELS = PIXELS * 2 / 1000;
  static constexpr uint32_t MAX_DROPPED_PIXELS = PIXELS * 5 / 1000;
  static constexpr uint64_t MAX_COVERAGE_LUMINANCE_ERROR = PIXELS * 5ull;
  static constexpr uint16_t MAX_COVERAGE_LUMINANCE = 9000;
};

} // namespace

int main() {
  const auto &corpus = mindsplatter_replay::HEAVY_SEARCH_V1;
  hs::random().seed(1337u);
  GenerativePalette::reset_hue_seed(0);
  configure_arenas_default();
  ReplayEffect effect;
  effect.init();
  const uint16_t particle_count = WhiteBox::restore_render(
      effect, std::span(corpus.state, corpus.state_size));

  std::printf("replay corpus id=%s source=%s preset=%d traits=%u particles=%u "
              "bytes=%zu hash=%llu revision=%s search_frame=%u "
              "search_score=%llu adaptive=%u long=%u peak_clip=%u\n",
              corpus.id, corpus.source, corpus.preset, corpus.traits,
              particle_count, corpus.state_size,
              static_cast<unsigned long long>(corpus.corpus_hash),
              corpus.source_revision, corpus.search_frame,
              static_cast<unsigned long long>(corpus.selection_score),
              corpus.search_adaptive_samples, corpus.search_long_edges,
              corpus.peak_clip);
  bool accepted = true;
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
        "replay frame=%d hash=%llu expected=%llu changed=%u "
        "channels=%u max=%u abs=%llu luma_abs=%llu luma_bias=%lld "
        "added=%u dropped=%u coverage_luma=%llu coverage_max=%u\n",
        frame + 1, static_cast<unsigned long long>(stats.framebuffer_hash),
        static_cast<unsigned long long>(stats.expected_hash),
        stats.changed_pixels, stats.changed_channels, stats.max_channel_error,
        static_cast<unsigned long long>(stats.total_absolute_error),
        static_cast<unsigned long long>(stats.total_luminance_error),
        static_cast<long long>(stats.luminance_bias), stats.added_pixels,
        stats.dropped_pixels,
        static_cast<unsigned long long>(stats.coverage_luminance_error),
        stats.max_coverage_luminance);
    const bool frame_accepted =
        stats.changed_pixels <= VisualGate::MAX_CHANGED_PIXELS &&
        stats.changed_channels <= VisualGate::MAX_CHANGED_CHANNELS &&
        stats.max_channel_error <= VisualGate::MAX_CHANNEL_ERROR &&
        stats.total_absolute_error <= VisualGate::MAX_TOTAL_ERROR &&
        stats.total_luminance_error <= VisualGate::MAX_LUMINANCE_ERROR &&
        stats.luminance_bias <= VisualGate::MAX_LUMINANCE_BIAS &&
        stats.luminance_bias >= -VisualGate::MAX_LUMINANCE_BIAS &&
        stats.added_pixels <= VisualGate::MAX_ADDED_PIXELS &&
        stats.dropped_pixels <= VisualGate::MAX_DROPPED_PIXELS &&
        stats.coverage_luminance_error <=
            VisualGate::MAX_COVERAGE_LUMINANCE_ERROR &&
        stats.max_coverage_luminance <=
            VisualGate::MAX_COVERAGE_LUMINANCE &&
        host_stats.expected_hash == corpus.framebuffer_hash;
    accepted = accepted && frame_accepted;
  }
  return accepted ? 0 : 1;
}
