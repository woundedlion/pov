/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "core/engine/memory.h"
#include "tests/mindsplatter_replay_corpus.h"
#include "tests/mindsplatter_whitebox.h"
#include "tests/mindsplatter_replay_metrics.h"

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
 * @brief Visual bounds for the two replay oracles.
 * @details The reference terms bound single-pass sample phasing: the candidate
 * against a reference render produced by the same binary. Measured over the
 * corpus frame, 2,954 changed pixels (7.1%), 5,718 changed channels (4.6%),
 * peak channel error 273, total absolute error 23,009 (0.18 counts per
 * channel). Counts retain at least 2x headroom, aggregate error 5x, and the
 * isolated overlap-order peak 7x.
 * @details The corpus terms bound the candidate against a golden recorded by
 * the generic-reference renderer, which differs by construction across a third
 * of the frame. Measured: 13,602 changed pixels (32.8%), 32,577 changed
 * channels (26.2%), peak channel error 13,799, total absolute error 5,994,146
 * (48.2 counts per channel), luminance error 1,850,744 (44.6 counts per pixel),
 * bias -59,060 (1.42 per pixel), 25 added and 66 dropped fringe pixels,
 * coverage error 2,919 with a 2,320 peak. Area and magnitude terms carry 1.5x,
 * the single-pixel peak 2x, covering the drift between the compiler recorded
 * in the corpus header and the one running the replay.
 * @details Every area bound is a whole-frame budget, so a clipped pass spends
 * it over fewer pixels. The per-pixel densities do not scale down with the
 * region: the peak-workload quadrant carries the frame's densest splats and
 * measures 2.06 luminance-bias counts per pixel against the frame's 1.42.
 */
struct VisualGate {
  static constexpr uint32_t PIXELS = WIDTH * HEIGHT;
  static constexpr uint32_t MAX_CHANGED_PIXELS = PIXELS * 16 / 100;
  static constexpr uint32_t MAX_CHANGED_CHANNELS = PIXELS * 3 * 11 / 100;
  static constexpr uint16_t MAX_CHANNEL_ERROR = 2000;
  static constexpr uint64_t MAX_TOTAL_ERROR = PIXELS * 3ull;
  static constexpr uint32_t MAX_CORPUS_CHANGED_PIXELS = PIXELS * 50 / 100;
  static constexpr uint32_t MAX_CORPUS_CHANGED_CHANNELS = PIXELS * 3 * 40 / 100;
  static constexpr uint16_t MAX_CORPUS_CHANNEL_ERROR = 28000;
  static constexpr uint64_t MAX_CORPUS_TOTAL_ERROR = PIXELS * 3ull * 72;
  static constexpr uint64_t MAX_LUMINANCE_ERROR = PIXELS * 61ull;
  static constexpr int64_t MAX_LUMINANCE_BIAS = PIXELS * 2ll;
  static constexpr uint32_t MAX_ADDED_PIXELS = PIXELS * 2 / 1000;
  static constexpr uint32_t MAX_DROPPED_PIXELS = PIXELS * 5 / 1000;
  static constexpr uint64_t MAX_COVERAGE_LUMINANCE_ERROR = PIXELS * 5ull;
  static constexpr uint16_t MAX_COVERAGE_LUMINANCE = 9000;
};

/** @brief Collects budget checks for one pass, naming each that fails. */
struct GateReport {
  const char *label;
  bool accepted = true;

  void check(const char *name, uint64_t achieved, uint64_t budget) {
    if (achieved <= budget)
      return;
    std::printf("replay %s over budget %s=%llu limit=%llu\n", label, name,
                static_cast<unsigned long long>(achieved),
                static_cast<unsigned long long>(budget));
    accepted = false;
  }
};

/**
 * @brief Renders one reference/candidate pass over the effect's current clip
 *        and gates the difference metrics.
 * @param effect Restored replay effect; its clip selects the compared region.
 * @param corpus Golden corpus the candidate is also measured against.
 * @param reference Scratch buffer holding the reference frame, WIDTH*HEIGHT.
 * @param label Pass label for the emitted line.
 * @param whole_frame True when the clip covers the whole frame, where the
 *        corpus's expanded pixels can be rehashed against the recorded hash to
 *        confirm the corpus itself is intact; the golden is recorded unclipped,
 *        so a partial pass rehashes only its own region.
 * @return True when every bound holds.
 */
bool run_pass(ReplayEffect &effect, const mindsplatter_replay::Corpus &corpus,
              std::vector<Pixel> &reference, const char *label,
              bool whole_frame) {
  mindsplatter_replay::ReferenceStats stats;
  mindsplatter_replay::FrameStats host_stats;
  ClipRegion clip;
  {
    Canvas canvas(effect);
    clip = canvas.clip();
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
      "replay %s clip=[%d,%d)x[%d,%d) hash=%llu expected=%llu changed=%u "
      "channels=%u max=%u abs=%llu corpus_pixels_hash=%llu corpus_changed=%u "
      "corpus_channels=%u corpus_max=%u corpus_abs=%llu luma_abs=%llu "
      "luma_bias=%lld added=%u dropped=%u coverage_luma=%llu "
      "coverage_max=%u\n",
      label, clip.x_start, clip.x_end, clip.y_start, clip.y_end,
      static_cast<unsigned long long>(stats.framebuffer_hash),
      static_cast<unsigned long long>(stats.expected_hash),
      stats.changed_pixels, stats.changed_channels, stats.max_channel_error,
      static_cast<unsigned long long>(stats.total_absolute_error),
      static_cast<unsigned long long>(host_stats.expected_hash),
      host_stats.changed_pixels, host_stats.changed_channels,
      host_stats.max_channel_error,
      static_cast<unsigned long long>(host_stats.total_absolute_error),
      static_cast<unsigned long long>(host_stats.total_luminance_error),
      static_cast<long long>(host_stats.luminance_bias),
      host_stats.added_pixels, host_stats.dropped_pixels,
      static_cast<unsigned long long>(host_stats.coverage_luminance_error),
      host_stats.max_coverage_luminance);

  GateReport gate{label};
  gate.check("changed", stats.changed_pixels, VisualGate::MAX_CHANGED_PIXELS);
  gate.check("channels", stats.changed_channels,
             VisualGate::MAX_CHANGED_CHANNELS);
  gate.check("max", stats.max_channel_error, VisualGate::MAX_CHANNEL_ERROR);
  gate.check("abs", stats.total_absolute_error, VisualGate::MAX_TOTAL_ERROR);
  gate.check("corpus_changed", host_stats.changed_pixels,
             VisualGate::MAX_CORPUS_CHANGED_PIXELS);
  gate.check("corpus_channels", host_stats.changed_channels,
             VisualGate::MAX_CORPUS_CHANGED_CHANNELS);
  gate.check("corpus_max", host_stats.max_channel_error,
             VisualGate::MAX_CORPUS_CHANNEL_ERROR);
  gate.check("corpus_abs", host_stats.total_absolute_error,
             VisualGate::MAX_CORPUS_TOTAL_ERROR);
  gate.check("luma_abs", host_stats.total_luminance_error,
             VisualGate::MAX_LUMINANCE_ERROR);
  gate.check("luma_bias_abs",
             static_cast<uint64_t>(host_stats.luminance_bias < 0
                                       ? -host_stats.luminance_bias
                                       : host_stats.luminance_bias),
             VisualGate::MAX_LUMINANCE_BIAS);
  gate.check("added", host_stats.added_pixels, VisualGate::MAX_ADDED_PIXELS);
  gate.check("dropped", host_stats.dropped_pixels,
             VisualGate::MAX_DROPPED_PIXELS);
  gate.check("coverage_luma", host_stats.coverage_luminance_error,
             VisualGate::MAX_COVERAGE_LUMINANCE_ERROR);
  gate.check("coverage_max", host_stats.max_coverage_luminance,
             VisualGate::MAX_COVERAGE_LUMINANCE);
  if (whole_frame && host_stats.expected_hash != corpus.framebuffer_hash) {
    std::printf("replay %s corpus integrity failed: expanded=%llu "
                "recorded=%llu\n",
                label,
                static_cast<unsigned long long>(host_stats.expected_hash),
                static_cast<unsigned long long>(corpus.framebuffer_hash));
    gate.accepted = false;
  }
  return gate.accepted;
}

} // namespace

int main() {
  const auto &corpus = mindsplatter_replay::HEAVY_SEARCH_V1;
  hs::random().seed(1337u);
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
    char label[16];
    std::snprintf(label, sizeof(label), "frame=%d", frame + 1);
    accepted = run_pass(effect, corpus, reference, label, true) && accepted;
  }

  // The corpus winner was scored under quadrant clips; replay the peak one so
  // the clip-boundary splat path the search stresses is actually rendered.
  const ClipRegion peak =
      mindsplatter_replay::search_clip<WIDTH, HEIGHT>(corpus.peak_clip);
  effect.set_clip(peak.y_start, peak.y_end, peak.x_start, peak.x_end);
  accepted =
      run_pass(effect, corpus, reference, "peak_clip", false) && accepted;
  effect.set_clip(0, HEIGHT, 0, WIDTH);

  return accepted ? 0 : 1;
}
