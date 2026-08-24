/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "core/engine/memory.h"
#include "tests/mindsplatter_replay_corpus.h"
#include "tests/mindsplatter_whitebox.h"
#include "tests/mindsplatter_replay_metrics.h"

#include <cstdio>
#include <cstring>
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
 * corpus frame, 6,767 changed pixels (16.3%), 11,907 changed channels (9.6%),
 * peak channel error 2,298, total absolute error 120,849 (0.97 counts per
 * channel). The sparse trail's longer edges amplify the optimized transform's
 * numerical drift; the bounds retain roughly 1.5-2x headroom.
 * @details The corpus terms bound the candidate against the recorded golden,
 * so they measure toolchain drift alone. Across the supported host toolchains,
 * the largest total error is 15,544, luminance error is 4,006, and luminance
 * bias is 1,056. Their bounds retain at least 2x headroom. The fringe and
 * coverage terms measure zero, so they keep a small absolute budget.
 * @details Every area bound is a whole-frame budget, so a clipped pass spends
 * it over fewer pixels. The per-pixel densities do not scale down with the
 * region: the peak-workload quadrant carries the frame's densest splats.
 */
struct VisualGate {
  static constexpr uint32_t PIXELS = WIDTH * HEIGHT;
  static constexpr uint32_t MAX_CHANGED_PIXELS = PIXELS * 25 / 100;
  static constexpr uint32_t MAX_CHANGED_CHANNELS = PIXELS * 3 * 11 / 100;
  static constexpr uint16_t MAX_CHANNEL_ERROR = 5000;
  static constexpr uint64_t MAX_TOTAL_ERROR = PIXELS * 6ull;
  static constexpr uint32_t MAX_CORPUS_CHANGED_PIXELS = PIXELS * 21 / 100;
  static constexpr uint32_t MAX_CORPUS_CHANGED_CHANNELS = PIXELS * 3 * 9 / 100;
  static constexpr uint16_t MAX_CORPUS_CHANNEL_ERROR = 216;
  static constexpr uint64_t MAX_CORPUS_TOTAL_ERROR = PIXELS * 4ull / 5;
  static constexpr uint64_t MAX_LUMINANCE_ERROR = PIXELS / 5;
  static constexpr int64_t MAX_LUMINANCE_BIAS = PIXELS / 16;
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
 * @brief Rehashes the corpus's own bytes the way mindsplatter_replay_gen
 *        recorded them: the restore blob, then every expanded golden channel.
 * @param corpus Corpus under replay.
 * @return FNV-1a over the state bytes followed by the whole framebuffer.
 * @details framebuffer_hash covers only the golden pixels; this also covers the
 * restore blob, so the two together pin every byte the corpus carries.
 */
uint64_t rehash_corpus(const mindsplatter_replay::Corpus &corpus) {
  uint64_t hash = 1469598103934665603ull;
  for (size_t i = 0; i < corpus.state_size; ++i)
    hash = (hash ^ corpus.state[i]) * 1099511628211ull;
  size_t entry = 0;
  for (size_t index = 0; index < static_cast<size_t>(WIDTH) * HEIGHT; ++index) {
    while (entry < corpus.framebuffer_entries &&
           corpus.framebuffer[entry].index < index)
      ++entry;
    const bool lit = entry < corpus.framebuffer_entries &&
                     corpus.framebuffer[entry].index == index;
    const uint16_t channels[] = {
        lit ? corpus.framebuffer[entry].r : uint16_t{0},
        lit ? corpus.framebuffer[entry].g : uint16_t{0},
        lit ? corpus.framebuffer[entry].b : uint16_t{0}};
    for (uint16_t channel : channels)
      hash = mindsplatter_replay::hash_channel(hash, channel);
  }
  return hash;
}

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
  if (std::strcmp(corpus.source_revision,
                  mindsplatter_replay::SOURCE_REVISION) != 0) {
    std::printf("replay corpus revision mismatch: corpus=%s expected=%s\n",
                corpus.source_revision, mindsplatter_replay::SOURCE_REVISION);
    accepted = false;
  }
  const uint64_t rehashed = rehash_corpus(corpus);
  if (rehashed != corpus.corpus_hash) {
    std::printf("replay corpus hash mismatch: rehashed=%llu recorded=%llu\n",
                static_cast<unsigned long long>(rehashed),
                static_cast<unsigned long long>(corpus.corpus_hash));
    accepted = false;
  }
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
