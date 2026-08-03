/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#include "core/engine/memory.h"
#include "tools/mindsplatter_replay_metrics.h"
#include "tools/mindsplatter_whitebox.h"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace {

constexpr int WIDTH = 288;
constexpr int HEIGHT = 144;
constexpr int FIRST_SEARCH_FRAME = 136;
constexpr int LAST_SEARCH_FRAME = 384;
constexpr int SEARCH_FRAME_STRIDE = 8;
constexpr unsigned long FRAME_MS = 16;
constexpr unsigned long FRAME_US = 16000;
constexpr uint32_t SEARCH_SEED = 1337;
constexpr const char *SOURCE_REVISION = "msp-heavy-search-v1";
constexpr uint32_t TRAIT_SATURATED = 1u << 0;
constexpr uint32_t TRAIT_LONG_EDGE = 1u << 3;
constexpr uint32_t TRAIT_MEASURED_WORST = 1u << 5;

using ReplayEffect = MindSplatter<WIDTH, HEIGHT>;
using WhiteBox = hs_test::effects_tests::MindSplatterWhiteBox;
using Snapshot = WhiteBox::ReplaySnapshot<WIDTH, HEIGHT>;

struct Workload {
  uint64_t score = 0;
  uint64_t adaptive_samples = 0;
  uint64_t long_edges = 0;
  uint64_t fragment_shader_calls = 0;
  uint64_t tap_writes = 0;

  Workload &operator+=(const Workload &other) {
    score += other.score;
    adaptive_samples += other.adaptive_samples;
    long_edges += other.long_edges;
    fragment_shader_calls += other.fragment_shader_calls;
    tap_writes += other.tap_writes;
    return *this;
  }
};

struct SearchResult {
  Snapshot snapshot;
  Workload aggregate;
  Workload peak_clip_workload;
  uint16_t frame = 0;
  uint8_t preset = 0;
  uint8_t peak_clip = 0;
};

struct GoldenPixel {
  uint16_t index;
  uint16_t r;
  uint16_t g;
  uint16_t b;
};

uint64_t fnv1a(uint64_t hash, uint8_t byte) {
  return (hash ^ byte) * 1099511628211ull;
}

Workload read_workload() {
  const hs::MindSplatterCounts &counts = hs::g_mindsplatter_counts;
  Workload workload;
  workload.adaptive_samples = counts.adaptive_samples;
  workload.long_edges = counts.long_edges;
  workload.fragment_shader_calls = counts.fragment_shader_calls;
  for (size_t taps = 1; taps < std::size(counts.aa_tap_masks); ++taps)
    workload.tap_writes += taps * counts.aa_tap_masks[taps];
  workload.score = workload.adaptive_samples * 64 + workload.long_edges * 512 +
                   workload.fragment_shader_calls * 8 + workload.tap_writes;
  return workload;
}

std::optional<SearchResult> search_corpus() {
  std::optional<SearchResult> best;
  for (uint8_t preset = 0; preset < 4; ++preset) {
    Workload preset_peak;
    uint16_t preset_peak_frame = 0;
    uint8_t preset_peak_clip = 0;
    GenerativePalette::reset_hue_seed(0);
    hs::random().seed(SEARCH_SEED);
    hs::set_mock_time(0, 0);
    configure_arenas_default();
    ReplayEffect effect;
    effect.init();
    WhiteBox::select_preset(effect, preset);
    effect.setAnimationsPaused(true);

    for (int frame = 1; frame <= LAST_SEARCH_FRAME; ++frame) {
      hs::set_mock_time(static_cast<unsigned long>(frame - 1) * FRAME_MS,
                        static_cast<unsigned long>(frame - 1) * FRAME_US);
      WhiteBox::step_state_without_render(effect);
      effect.advance_display();
      if (frame < FIRST_SEARCH_FRAME ||
          (frame - FIRST_SEARCH_FRAME) % SEARCH_FRAME_STRIDE != 0)
        continue;

      Workload aggregate;
      Workload peak_clip_workload;
      uint8_t peak_clip = 0;
      for (uint8_t clip_index = 0;
           clip_index < mindsplatter_replay::SEARCH_CLIP_COUNT; ++clip_index) {
        const ClipRegion clip =
            mindsplatter_replay::search_clip<WIDTH, HEIGHT>(clip_index);
        effect.set_clip(clip.y_start, clip.y_end, clip.x_start, clip.x_end);
        hs::g_mindsplatter_counts.reset();
        {
          Canvas canvas(effect);
          WhiteBox::draw_particles_replay_reference(effect, canvas);
        }
        effect.advance_display();
        const Workload workload = read_workload();
        aggregate += workload;
        if (workload.score > peak_clip_workload.score) {
          peak_clip_workload = workload;
          peak_clip = clip_index;
        }
      }
      effect.set_clip(0, HEIGHT, 0, WIDTH);

      if (!best || aggregate.score > best->aggregate.score) {
        SearchResult result;
        result.snapshot = WhiteBox::capture(effect);
        result.aggregate = aggregate;
        result.peak_clip_workload = peak_clip_workload;
        result.frame = static_cast<uint16_t>(frame);
        result.preset = preset;
        result.peak_clip = peak_clip;
        best = std::move(result);
      }
      if (aggregate.score > preset_peak.score) {
        preset_peak = aggregate;
        preset_peak_frame = static_cast<uint16_t>(frame);
        preset_peak_clip = peak_clip;
      }
    }
    std::printf("search preset=%u peak_frame=%u score=%llu adaptive=%llu "
                "long=%llu peak_clip=%u\n",
                static_cast<unsigned>(preset),
                static_cast<unsigned>(preset_peak_frame),
                static_cast<unsigned long long>(preset_peak.score),
                static_cast<unsigned long long>(preset_peak.adaptive_samples),
                static_cast<unsigned long long>(preset_peak.long_edges),
                static_cast<unsigned>(preset_peak_clip));
  }
  return best;
}

template <typename T>
void emit_array(std::ostream &out, const char *type, const char *name,
                const std::vector<T> &values, int columns) {
  out << "inline const " << type << ' ' << name << "[] PROGMEM = {\n";
  for (size_t i = 0; i < values.size(); ++i) {
    if (i % columns == 0)
      out << "    ";
    out << static_cast<uint64_t>(values[i]) << ',';
    if (i % columns == static_cast<size_t>(columns - 1) ||
        i + 1 == values.size())
      out << '\n';
    else
      out << ' ';
  }
  out << "};\n\n";
}

void emit_golden(std::ostream &out, const std::vector<GoldenPixel> &pixels) {
  out << "inline const GoldenPixel HEAVY_SEARCH_V1_FRAMEBUFFER[] PROGMEM = {\n";
  for (const GoldenPixel &pixel : pixels)
    out << "    {" << pixel.index << ", " << pixel.r << ", " << pixel.g << ", "
        << pixel.b << "},\n";
  out << "};\n\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: mindsplatter_replay_gen <output-header>\n");
    return 2;
  }

  std::optional<SearchResult> selected = search_corpus();
  if (!selected) {
    std::fprintf(stderr, "MindSplatter replay search produced no candidates\n");
    return 1;
  }

  const std::vector<unsigned char> state =
      WhiteBox::serialize_render(selected->snapshot);
  GenerativePalette::reset_hue_seed(0);
  hs::random().seed(SEARCH_SEED);
  configure_arenas_default();
  ReplayEffect effect;
  effect.init();
  WhiteBox::restore(effect, selected->snapshot);
  effect.set_clip(0, HEIGHT, 0, WIDTH);
  WhiteBox::draw_particles(effect);
  effect.advance_display();

  const Pixel *pixels = effect.display_buffer();
  std::vector<uint16_t> framebuffer;
  framebuffer.reserve(static_cast<size_t>(WIDTH) * HEIGHT * 3);
  std::vector<GoldenPixel> golden;
  uint64_t framebuffer_hash = 1469598103934665603ull;
  for (int i = 0; i < WIDTH * HEIGHT; ++i) {
    if (pixels[i].r | pixels[i].g | pixels[i].b)
      golden.push_back(
          {static_cast<uint16_t>(i), pixels[i].r, pixels[i].g, pixels[i].b});
    for (uint16_t channel : {pixels[i].r, pixels[i].g, pixels[i].b}) {
      framebuffer.push_back(channel);
      framebuffer_hash = fnv1a(framebuffer_hash, channel & 0xffu);
      framebuffer_hash = fnv1a(framebuffer_hash, channel >> 8);
    }
  }
  uint64_t corpus_hash = 1469598103934665603ull;
  for (unsigned char byte : state)
    corpus_hash = fnv1a(corpus_hash, byte);
  for (uint16_t channel : framebuffer) {
    corpus_hash = fnv1a(corpus_hash, channel & 0xffu);
    corpus_hash = fnv1a(corpus_hash, channel >> 8);
  }

  const ClipRegion peak_clip =
      mindsplatter_replay::search_clip<WIDTH, HEIGHT>(selected->peak_clip);
  const std::string corpus_id = "heavy_search_v1_p" +
                                std::to_string(selected->preset) + "_f" +
                                std::to_string(selected->frame);
  const std::string source =
      "seed=1337 presets=0..3 frames=136..384/8 clips=quadrants "
      "renderer=generic-reference "
      "score=64*adaptive+512*long+8*shader+taps compiler=" __VERSION__;
  uint32_t traits = TRAIT_LONG_EDGE | TRAIT_MEASURED_WORST;
  if (selected->snapshot.particles.size() ==
      WhiteBox::particle_capacity(effect))
    traits |= TRAIT_SATURATED;

  std::ofstream out(argv[1], std::ios::binary | std::ios::trunc);
  if (!out) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 1;
  }
  out << "/*\n"
         " * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.\n"
         " * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit\n"
         " * permission.\n"
         " * Generated by tools/mindsplatter_replay_gen.cpp.\n"
         " */\n"
         "#pragma once\n\n"
         "#include \"core/engine/platform.h\"\n\n"
         "#include <cstddef>\n"
         "#include <cstdint>\n\n"
         "namespace mindsplatter_replay {\n\n"
         "struct GoldenPixel {\n"
         "  uint16_t index;\n"
         "  uint16_t r;\n"
         "  uint16_t g;\n"
         "  uint16_t b;\n"
         "};\n\n"
         "enum CorpusTraits : uint32_t {\n"
         "  CORPUS_SATURATED = 1u << 0,\n"
         "  CORPUS_POLE = 1u << 1,\n"
         "  CORPUS_SEAM = 1u << 2,\n"
         "  CORPUS_LONG_EDGE = 1u << 3,\n"
         "  CORPUS_ONE_DOT_HEAVY = 1u << 4,\n"
         "  CORPUS_MEASURED_WORST = 1u << 5,\n"
         "};\n\n"
         "struct Corpus {\n"
         "  const char *id;\n"
         "  const char *source;\n"
         "  const char *source_revision;\n"
         "  int8_t preset;\n"
         "  uint32_t traits;\n"
         "  const unsigned char *state;\n"
         "  size_t state_size;\n"
         "  const GoldenPixel *framebuffer;\n"
         "  size_t framebuffer_entries;\n"
         "  uint16_t particle_count;\n"
         "  uint16_t search_frame;\n"
         "  uint8_t peak_clip;\n"
         "  uint64_t selection_score;\n"
         "  uint32_t search_adaptive_samples;\n"
         "  uint32_t search_long_edges;\n"
         "  uint64_t corpus_hash;\n"
         "  uint64_t framebuffer_hash;\n"
         "};\n\n";
  emit_array(out, "unsigned char", "HEAVY_SEARCH_V1_STATE", state, 16);
  emit_golden(out, golden);
  out << "inline const Corpus HEAVY_SEARCH_V1 = {\n"
      << "    \"" << corpus_id << "\",\n"
      << "    \"" << source << "\",\n"
      << "    \"" << SOURCE_REVISION << "\",\n"
      << "    " << static_cast<unsigned>(selected->preset) << ",\n"
      << "    " << traits << ",\n"
      << "    HEAVY_SEARCH_V1_STATE,\n"
      << "    sizeof(HEAVY_SEARCH_V1_STATE),\n"
      << "    HEAVY_SEARCH_V1_FRAMEBUFFER,\n"
      << "    sizeof(HEAVY_SEARCH_V1_FRAMEBUFFER) / sizeof(GoldenPixel),\n"
      << "    " << selected->snapshot.particles.size() << ",\n"
      << "    " << selected->frame << ",\n"
      << "    " << static_cast<unsigned>(selected->peak_clip) << ",\n"
      << "    " << selected->aggregate.score << "ull,\n"
      << "    " << selected->aggregate.adaptive_samples << ",\n"
      << "    " << selected->aggregate.long_edges << ",\n"
      << "    " << corpus_hash << "ull,\n"
      << "    " << framebuffer_hash << "ull};\n\n"
      << "inline const Corpus *const CORPUS_MANIFEST[] = {&HEAVY_SEARCH_V1};\n\n"
      << "} // namespace mindsplatter_replay\n";
  hs::clear_mock_time();
  std::printf(
      "%s revision=%s preset=%u frame=%u particles=%zu score=%llu "
      "adaptive=%llu long=%llu peak_clip=%u[%d,%d,%d,%d] "
      "peak_score=%llu state=%zu framebuffer=%zu hash=%llu\n",
      corpus_id.c_str(), SOURCE_REVISION,
      static_cast<unsigned>(selected->preset),
      static_cast<unsigned>(selected->frame),
      selected->snapshot.particles.size(),
      static_cast<unsigned long long>(selected->aggregate.score),
      static_cast<unsigned long long>(selected->aggregate.adaptive_samples),
      static_cast<unsigned long long>(selected->aggregate.long_edges),
      static_cast<unsigned>(selected->peak_clip), peak_clip.x_start,
      peak_clip.x_end, peak_clip.y_start, peak_clip.y_end,
      static_cast<unsigned long long>(selected->peak_clip_workload.score),
      state.size(), framebuffer.size() * sizeof(uint16_t),
      static_cast<unsigned long long>(corpus_hash));
  return out ? 0 : 1;
}
