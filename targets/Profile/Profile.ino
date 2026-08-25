/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Profile — single-effect on-device profiling harness (288×144, segmented)
 *
 * Target: one Teensy 4.0 running as segment 0 of the shipping 4-segment
 * Phantasm configuration. Runs exactly one effect (selected at build time via
 * -D HS_PROFILE_TARGET=<EffectClass>, default DisplacementField) under the
 * real POVSegmented driver — flywheel ISR, DMA LED output, and segment
 * clipping all live — and periodically dumps the HS_PROFILE cycle-counter
 * tree plus exact per-frame wall-clock stats over USB serial.
 *
 * Build with the `profile` PlatformIO env (which defines HS_PROFILE_ENABLE);
 * drive it via `just profile <EffectClass>`.
 */

// Effect class to profile; overridden per-run by `just profile <EffectClass>`.
#ifndef HS_PROFILE_TARGET
#define HS_PROFILE_TARGET DisplacementField
#endif

// Frames per readout window; override (via PLATFORMIO_BUILD_FLAGS) for finer
// phase resolution at the cost of more serial traffic.
#ifndef HS_PROFILE_WINDOW
#define HS_PROFILE_WINDOW 128
#endif

#ifndef HS_PROFILE_ENABLE
#error "Profile builds require HS_PROFILE_ENABLE (use the `profile` env)"
#endif

#ifndef HS_PROFILE_CONFIG_TAG
#error "Profile builds require HS_PROFILE_CONFIG_TAG"
#endif

// Further per-run knobs consumed elsewhere (all via PLATFORMIO_BUILD_FLAGS):
//   HS_PROFILE_EPOCH_REVS    epoch length override (pov_segmented.h) so one
//                            instance covers a full preset cycle
//   HS_PROFILE_ORDERED_CYCLE random-next cyclers advance in order instead
//                            (HankinSolids, SphericalHarmonics)
//   HS_PROFILE_TRANS_SPEED   "Trans Speed" applied after init (below)
//   HS_PROFILE_PRESET        zero-based fixed preset selected after init
//   HS_PROFILE_AUTO_RESUME   resume choreography after selecting that preset
//   HS_SCAN_METRICS          compiles in the per-pixel hs::g_scan_metrics probe
//                            counters and adds the "scan totals" window line.
//                            Every probe pays a non-atomic global increment, so
//                            read the COUNTS from such a build, not the times.
//   HS_PROBE_BREAKDOWN       compiles in the per-probe stage cycle buckets and
//                            adds the "probe cycles"/"probe counts" lines. Every
//                            stage boundary is a cycle-counter read, so read
//                            RATIOS from such a build, discounted by `tick`.
//   HS_PLOT_COUNTS           compiles in count-only Plot workload attribution
//                            and adds one "plot counts" line per window.
//   HS_MINDSPLATTER_REPLAY    restores the selected frozen MindSplatter corpus
//                            and renders it repeatedly without physics.
//   HS_MINDSPLATTER_REPLAY_AB renders the generic baseline and candidate from
//                            one state and reports their pixel error.
//   HS_MINDSPLATTER_REPLAY_CORPUS
//                            corpus symbol selected from the generated header.
//   HS_PROFILE_MINDSPLATTER_COUNTS
//                            adds count-only MindSplatter path attribution.
//   HS_PROFILE_MINDSPLATTER_STALLS
//                            adds short-batch DWT cycle/stall attribution.
//   HS_PROFILE_SHADER_WORKBENCH_STAGES
//                            adds raw ShaderWorkbench pipeline cycle totals.

#if defined(HS_PROFILE_MINDSPLATTER_COUNTS) &&                                 \
    defined(HS_PROFILE_MINDSPLATTER_STALLS)
#error "MindSplatter count and stall captures require separate images"
#endif

#define HS_PROFILE_STR2(x) #x
#define HS_PROFILE_STR(x) HS_PROFILE_STR2(x)

#ifdef HS_PROFILE_PULLBACK_TELEMETRY
#ifndef HS_PULLBACK_ARM
#define HS_PULLBACK_ARM LANDED
#endif
#ifndef HS_PULLBACK_SHORT_SHA
#error "Pullback telemetry requires HS_PULLBACK_SHORT_SHA"
#endif
#endif

#include "../Phantasm/phantasm_target.h"

#ifdef HS_MINDSPLATTER_REPLAY
#include "tests/mindsplatter_replay_corpus.h"
#include "tests/mindsplatter_whitebox.h"
#include "tests/mindsplatter_replay_metrics.h"
#ifndef HS_MINDSPLATTER_REPLAY_CORPUS
#define HS_MINDSPLATTER_REPLAY_CORPUS HEAVY_SEARCH_V1
#endif
#endif

#if defined(HS_MINDSPLATTER_REPLAY_AB) && !defined(HS_MINDSPLATTER_REPLAY)
#error "HS_MINDSPLATTER_REPLAY_AB requires HS_MINDSPLATTER_REPLAY"
#endif

// The .ino -> .cpp converter injects prototypes for every function it detects
// immediately before the first one, which here sits inside the anonymous
// namespace below — giving the injected setup/loop internal linkage they never
// get a definition for. Declaring them at global scope suppresses the injection.
FLASHMEM void setup();
void loop();

namespace {

#ifdef HS_MINDSPLATTER_REPLAY_AB
static Pixel replay_reference_pixels[MAX_W * MAX_H / NUM_SEGMENTS];
#endif

/**
 * @brief Wraps the profiled effect: roots the HS_PROFILE tree at one `frame`
 * counter per draw_frame() and accumulates exact wall-clock stats, dumping
 * both every WINDOW_FRAMES frames.
 */
template <int W, int H> class ProfiledEffect : public HS_PROFILE_TARGET<W, H> {
#ifdef HS_PROFILE_EFFECT_HEAP_BYTES
  alignas(16) static inline uint8_t
      effect_storage[HS_PROFILE_EFFECT_HEAP_BYTES]{};
#endif

public:
#ifdef HS_PROFILE_EFFECT_HEAP_BYTES
  static void *operator new(size_t size, const std::nothrow_t &) noexcept {
    return size <= sizeof(effect_storage) ? effect_storage : nullptr;
  }

  static void operator delete(void *) noexcept {}
  static void operator delete(void *, const std::nothrow_t &) noexcept {}
#endif

#ifdef HS_MINDSPLATTER_REPLAY
  void init() override {
    using Target = HS_PROFILE_TARGET<W, H>;
    static_assert(std::is_same_v<Target, MindSplatter<W, H>>,
                  "HS_MINDSPLATTER_REPLAY requires MindSplatter");
    Target::init();
    const auto &corpus = replay_corpus();
    HS_CHECK(corpus.framebuffer_entries <= static_cast<size_t>(W) * H,
             "MindSplatter replay framebuffer entries differ");
    const uint16_t particles = ReplayWhiteBox::restore_render(
        *this, {corpus.state, corpus.state_size});
    HS_CHECK(particles == corpus.particle_count,
             "MindSplatter replay particle count differs");
    char hash[21], score[21];
    hs::log("replay corpus: id=%s preset=%d traits=%lu particles=%u bytes=%lu "
            "hash=%s",
            corpus.id, (int)corpus.preset, (unsigned long)corpus.traits,
            (unsigned)particles, (unsigned long)corpus.state_size,
            hs::u64_dec(corpus.corpus_hash, hash));
    hs::log("replay source: %s", corpus.source);
    hs::log("replay search: revision=%s frame=%u peak_clip=%u score=%s "
            "adaptive=%lu long=%lu",
            corpus.source_revision, (unsigned)corpus.search_frame,
            (unsigned)corpus.peak_clip,
            hs::u64_dec(corpus.selection_score, score),
            (unsigned long)corpus.search_adaptive_samples,
            (unsigned long)corpus.search_long_edges);
#ifdef HS_MINDSPLATTER_REPLAY_AB
    hs::log("replay compare: candidate vs generic device reference; "
            "timing=combined");
#else
    hs::log("replay compare: candidate vs host-generated fingerprint");
#endif
  }
#endif

  void draw_frame() override {
    const uint64_t bw0 = buffer_wait ? buffer_wait->cycles : 0;
    const unsigned long t0 = micros();
#ifdef HS_MINDSPLATTER_REPLAY
    const Pixel *pixels = nullptr;
#endif
    {
      HS_PROFILE(frame);
#ifdef HS_MINDSPLATTER_REPLAY
#ifdef HS_MINDSPLATTER_REPLAY_AB
      Canvas canvas(*this);
      replay_clip = canvas.clip();
      ReplayWhiteBox::draw_particles_replay_reference(*this, canvas);
      const size_t reference_count = mindsplatter_replay::capture_frame<W>(
          canvas.data(), replay_clip, replay_reference_pixels,
          std::size(replay_reference_pixels));
      mindsplatter_replay::clear_frame(canvas, replay_clip);
      ReplayWhiteBox::draw_particles_candidate(*this, canvas);
      pixels = canvas.data();
      replay_stats = mindsplatter_replay::compare_frame_reference<W>(
          pixels, replay_clip, replay_reference_pixels, reference_count);
#else
      ReplayWhiteBox::draw_particles_inspect(
          *this, [&](Canvas &canvas, const ClipRegion &clip) {
            pixels = &canvas(0, 0);
            replay_clip = clip;
          });
#endif
#else
      HS_PROFILE_TARGET<W, H>::draw_frame();
#endif
    }
    const unsigned long dt = micros() - t0;
#ifdef HS_MINDSPLATTER_REPLAY
#ifndef HS_MINDSPLATTER_REPLAY_AB
    replay_stats = mindsplatter_replay::compare_frame<W>(
        pixels, replay_clip, replay_corpus().framebuffer,
        replay_corpus().framebuffer_entries);
#endif
#endif
    // Per-frame render = wall minus this frame's display-sync wait, read as
    // the effect's *_buffer_wait counter delta. The counter self-registers on
    // the first draw_frame, so the lookup retries until it appears.
    if (!buffer_wait)
      buffer_wait = hs::CycleCounter::find_suffix("_buffer_wait");
    const uint64_t bw1 = buffer_wait ? buffer_wait->cycles : 0;
    const unsigned long wait_us =
        (unsigned long)((bw1 - bw0) / hs::CycleCounter::CYCLES_PER_US);
    const unsigned long render = dt > wait_us ? dt - wait_us : 0;
    // One compact line per frame (the full counter tree still dumps per
    // window — per-frame log_all would perturb the frames it measures).
    hs::log("f %lu w=%lu r=%lu", total_frames + 1, dt, render);
#ifdef HS_SCAN_METRICS
    drain_scan_metrics();
#endif
#ifdef HS_PROBE_BREAKDOWN
    drain_probe_breakdown();
#endif
#ifdef HS_PROFILE_SHADER_WORKBENCH_STAGES
    drain_shader_workbench_stages();
#endif
    render_sum += render;
    if (render > render_max)
      render_max = render;
    wall_sum += dt;
    if (dt < wall_min)
      wall_min = dt;
    if (dt > wall_max)
      wall_max = dt;
    ++total_frames;
    if (++window_frames == WINDOW_FRAMES)
      dump();
  }

private:
  static constexpr int WINDOW_FRAMES =
      HS_PROFILE_WINDOW; /**< Frames per readout window. */

#ifdef HS_MINDSPLATTER_REPLAY
  using ReplayWhiteBox = hs_test::effects_tests::MindSplatterWhiteBox;

  static const mindsplatter_replay::Corpus &replay_corpus() {
    return mindsplatter_replay::HS_MINDSPLATTER_REPLAY_CORPUS;
  }

  void dump_replay_stats() const {
    char actual[21], expected[21], total_error[21];
#ifdef HS_MINDSPLATTER_REPLAY_AB
    static constexpr const char *EXPECTED_LABEL = "reference";
#else
    static constexpr const char *EXPECTED_LABEL = "host";
#endif
    hs::log("replay frame: clip=%d,%d,%d,%d candidate=%s %s=%s "
            "changed=%lu channels=%lu max=%u abs=%s",
            replay_stats.clip.x_start, replay_stats.clip.x_end,
            replay_stats.clip.y_start, replay_stats.clip.y_end,
            hs::u64_dec(replay_stats.framebuffer_hash, actual), EXPECTED_LABEL,
            hs::u64_dec(replay_stats.expected_hash, expected),
            (unsigned long)replay_stats.changed_pixels,
            (unsigned long)replay_stats.changed_channels,
            (unsigned)replay_stats.max_channel_error,
            hs::u64_dec(replay_stats.total_absolute_error, total_error));
  }
#endif

  void dump() {
    const unsigned long now = micros();
    hs::log("=== profile %s [%dx%d] frames %lu-%lu window=%lu us ===",
            HS_PROFILE_STR(HS_PROFILE_TARGET), W, H,
            total_frames - window_frames + 1, total_frames, now - window_start);
    hs::log("frame wall us: min=%lu avg=%lu max=%lu sum=%lu (%d frames)",
            wall_min, wall_sum / WINDOW_FRAMES, wall_max, wall_sum,
            WINDOW_FRAMES);
    hs::log("frame render us: avg=%lu max=%lu", render_sum / WINDOW_FRAMES,
            render_max);
    hs::CycleCounter::log_all();
#ifdef HS_SCAN_METRICS
    dump_scan_totals();
#endif
#ifdef HS_PROBE_BREAKDOWN
    dump_probe_breakdown();
#endif
#ifdef HS_PLOT_COUNTS
    dump_plot_counts();
#endif
#ifdef HS_MINDSPLATTER_REPLAY
    dump_replay_stats();
#endif
#ifdef HS_PROFILE_MINDSPLATTER_COUNTS
    dump_mindsplatter_counts();
#endif
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    dump_mindsplatter_stalls();
#endif
#ifdef HS_PROFILE_SHADER_WORKBENCH_STAGES
    dump_shader_workbench_stages();
#endif
    dump_isr_stats(now - window_start);
    hs::CycleCounter::reset_all();
    window_frames = 0;
    wall_sum = 0;
    wall_min = ~0ul;
    wall_max = 0;
    render_sum = 0;
    render_max = 0;
    window_start = micros();
  }

  /**
   * @brief Prints and resets the column-ISR accumulators for this window.
   * @param window_us Window wall-clock span, for the CPU-share figure.
   * @details Copy + reset under a brief IRQ-off window (the ISR is the sole
   * writer). Per-call figures are exact cycles; us/ns derive from
   * CycleCounter::CYCLES_PER_US.
   */
  static void dump_isr_stats(unsigned long window_us) {
    hs::IsrCycleStats wake, pack, submit;
    const uint32_t primask = hs::save_disable_interrupts();
    wake = hs::g_flywheel_wake_cycles;
    pack = hs::g_column_pack_cycles;
    submit = hs::g_dma_submit_cycles;
    hs::g_flywheel_wake_cycles.reset();
    hs::g_column_pack_cycles.reset();
    hs::g_dma_submit_cycles.reset();
    hs::restore_interrupts(primask);
    log_isr("isr_wake", wake, window_us);
    log_isr("isr_pack", pack, window_us);
    log_isr("isr_dma_submit", submit, window_us);
  }

  static_assert((uint64_t)hs::CycleCounter::CYCLES_PER_US * 1000000u ==
                    (uint64_t)F_CPU,
                "CycleCounter::CYCLES_PER_US disagrees with F_CPU");

  /**
   * @brief Converts a cycle count to nanoseconds at the core clock.
   * @details The 64-bit intermediate keeps cycles * 1000 from wrapping 32 bits.
   */
  static constexpr uint32_t cycles_to_ns(uint32_t cycles) {
    return (uint32_t)((uint64_t)cycles * 1000u /
                      hs::CycleCounter::CYCLES_PER_US);
  }

  static void log_isr(const char *name, const hs::IsrCycleStats &s,
                      unsigned long window_us) {
    if (!s.count) {
      hs::log("%-14s n=0", name);
      return;
    }
    const uint32_t avg = (uint32_t)(s.cycles / s.count);
    const uint32_t total_us =
        (uint32_t)(s.cycles / hs::CycleCounter::CYCLES_PER_US);
    // CPU share in hundredths of a percent; window_us is far below the
    // 32-bit ceiling of total_us * 10000.
    const uint32_t share_c =
        window_us ? (uint32_t)((uint64_t)total_us * 10000u / window_us) : 0;
    hs::log("%-14s n=%lu cyc min/avg/max=%lu/%lu/%lu "
            "ns min/avg/max=%lu/%lu/%lu total=%lu us cpu=%lu.%02lu%%",
            name, (unsigned long)s.count, (unsigned long)s.min,
            (unsigned long)avg, (unsigned long)s.max,
            (unsigned long)cycles_to_ns(s.min),
            (unsigned long)cycles_to_ns(avg),
            (unsigned long)cycles_to_ns(s.max), (unsigned long)total_us,
            (unsigned long)(share_c / 100u), (unsigned long)(share_c % 100u));
  }

#ifdef HS_PROFILE_SHADER_WORKBENCH_STAGES
  struct ShaderWorkbenchStageTotals {
    uint64_t lens = 0;
    uint64_t surface_noise = 0;
    uint64_t projection = 0;
    uint64_t planar_warp = 0;
    uint64_t source = 0;
    uint64_t material = 0;
    uint64_t color = 0;
    uint64_t mirror_tile = 0;
    uint64_t polyhedral_pixels = 0;
    uint64_t polyhedral_reflections = 0;
    uint32_t polyhedral_max_reflections = 0;

    void reset() { *this = {}; }
  };

  void drain_shader_workbench_stages() {
    const hs::ShaderWorkbenchStageCycles &frame = hs::g_shader_workbench_stage_cycles;
    shader_workbench_stage_totals.lens += frame.lens;
    shader_workbench_stage_totals.surface_noise += frame.surface_noise;
    shader_workbench_stage_totals.projection += frame.projection;
    shader_workbench_stage_totals.planar_warp += frame.planar_warp;
    shader_workbench_stage_totals.source += frame.source;
    shader_workbench_stage_totals.material += frame.material;
    shader_workbench_stage_totals.color += frame.color;
    shader_workbench_stage_totals.mirror_tile += frame.mirror_tile;
    shader_workbench_stage_totals.polyhedral_pixels += frame.polyhedral_pixels;
    shader_workbench_stage_totals.polyhedral_reflections +=
        frame.polyhedral_reflections;
    shader_workbench_stage_totals.polyhedral_max_reflections =
        std::max(shader_workbench_stage_totals.polyhedral_max_reflections,
                 frame.polyhedral_max_reflections);
    hs::g_shader_workbench_stage_cycles.reset();
  }

  void dump_shader_workbench_stages() {
    char c0[21], c1[21], c2[21], c3[21], c4[21], c5[21], c6[21];
    hs::log("sb stages: lens=%s surface_noise=%s projection=%s warp=%s "
            "source=%s material=%s color=%s",
            hs::u64_dec(shader_workbench_stage_totals.lens, c0),
            hs::u64_dec(shader_workbench_stage_totals.surface_noise, c1),
            hs::u64_dec(shader_workbench_stage_totals.projection, c2),
            hs::u64_dec(shader_workbench_stage_totals.planar_warp, c3),
            hs::u64_dec(shader_workbench_stage_totals.source, c4),
            hs::u64_dec(shader_workbench_stage_totals.material, c5),
            hs::u64_dec(shader_workbench_stage_totals.color, c6));
    char c7[21], c8[21], c9[21];
    hs::log("sb detail: mirror=%s poly_pixels=%s poly_reflections=%s "
            "poly_max=%lu",
            hs::u64_dec(shader_workbench_stage_totals.mirror_tile, c7),
            hs::u64_dec(shader_workbench_stage_totals.polyhedral_pixels, c8),
            hs::u64_dec(shader_workbench_stage_totals.polyhedral_reflections, c9),
            (unsigned long)shader_workbench_stage_totals.polyhedral_max_reflections);
    shader_workbench_stage_totals.reset();
  }

  ShaderWorkbenchStageTotals shader_workbench_stage_totals;
#endif

#ifdef HS_SCAN_METRICS
  /** @brief 64-bit window accumulators for the per-pixel scan counters. */
  struct ScanTotals {
    uint64_t tested = 0; /**< Face::distance probes. */
    uint64_t culled =
        0; /**< Probes rejected by the back-face / radius guards. */
    uint64_t exact =
        0; /**< Probes taking a full evaluation (convex + sector + walk). */
    uint64_t convex = 0; /**< Full evaluations on the convex half-plane path. */
    uint64_t sector = 0; /**< Full evaluations on the concave sector walk. */
    uint64_t lut = 0;    /**< Probes served by the class-LUT bilinear fetch. */
    uint64_t cand = 0;   /**< Pixels passing the scan's d < pixel_width test. */
    uint64_t backstop = 0; /**< plot() steps_cache capacity-backstop trips. */
    void reset() {
      tested = culled = exact = convex = sector = lut = cand = backstop = 0;
    }
  };

  /**
   * @brief Folds this frame's scan counters into the window totals, rezeroed.
   * @details The source counters are uint32; a long window of a many-faced
   *          solid probes enough pixels to wrap one, so they drain per frame.
   */
  void drain_scan_metrics() {
    const hs::ScanMetrics &m = hs::g_scan_metrics;
    scan_totals.tested += m.pixels_tested;
    scan_totals.culled += m.pixels_culled;
    scan_totals.exact += m.exact_hits;
    scan_totals.convex += m.convex_hits;
    scan_totals.sector += m.sector_hits;
    scan_totals.lut += m.lut_hits;
    scan_totals.cand += m.shade_candidates;
    scan_totals.backstop += m.plot_backstop_hits;
    hs::g_scan_metrics.reset();
  }

  /**
   * @brief Prints and resets this window's scan-probe totals.
   * @details Window totals, like the counter tree; divide by the header's frame
   *          count for per-frame figures. walk = exact - convex - sector is the
   *          residual exact-edge walk. Alpha survivors are the raster_shade
   *          scope's call count, already in the tree above.
   */
  void dump_scan_totals() {
    const ScanTotals &t = scan_totals;
    const uint64_t walk = t.exact - t.convex - t.sector;
    char b0[21], b1[21], b2[21], b3[21], b4[21], b5[21], b6[21], b7[21];
    hs::log("scan totals: tested=%s culled=%s lut=%s convex=%s sector=%s "
            "walk=%s cand=%s backstop=%s",
            hs::u64_dec(t.tested, b0), hs::u64_dec(t.culled, b1),
            hs::u64_dec(t.lut, b2), hs::u64_dec(t.convex, b3),
            hs::u64_dec(t.sector, b4), hs::u64_dec(walk, b5),
            hs::u64_dec(t.cand, b6), hs::u64_dec(t.backstop, b7));
    scan_totals.reset();
  }

  ScanTotals scan_totals; /**< This window's drained scan counters. */
#endif

#ifdef HS_PLOT_COUNTS
  static void dump_plot_counts() {
    const hs::PlotCounts &t = hs::g_plot_counts;
    hs::log(
        "plot counts: r=%lu,e=%lu,p=%lu,g=%lu,d=%lu,o=%lu,t=%lu,c=%lu,"
        "s=%lu,y=%lu,u=%lu,a=%lu,n=%lu,h=%lu,x=%lu,k=%lu,b=%lu",
        (unsigned long)t.rings, (unsigned long)t.edges, (unsigned long)t.planar,
        (unsigned long)t.geodesic, (unsigned long)t.degenerate,
        (unsigned long)t.one_dot, (unsigned long)t.cull_tests,
        (unsigned long)t.culled, (unsigned long)t.sim_samples,
        (unsigned long)t.replay_samples, (unsigned long)t.planar_unprojects,
        (unsigned long)t.planar_arc_samples, (unsigned long)t.normalizations,
        (unsigned long)t.shader_calls, (unsigned long)t.plotted_samples,
        (unsigned long)t.steps_peak, (unsigned long)t.backstops);
    hs::g_plot_counts.reset();
  }
#endif

#ifdef HS_PROFILE_MINDSPLATTER_COUNTS
  static void dump_mindsplatter_counts() {
    const hs::MindSplatterCounts &t = hs::g_mindsplatter_counts;
    hs::log("msp counts particles: resident=%lu live=%lu full=%lu partial=%lu "
            "draining=%lu",
            (unsigned long)t.resident_particles,
            (unsigned long)t.live_particles, (unsigned long)t.full_histories,
            (unsigned long)t.partial_histories,
            (unsigned long)t.draining_histories);
    hs::log("msp counts gate: cart_lat=%lu cart_mer=%lu cart_fallback=%lu "
            "row=%lu col=%lu edge=%lu visible=%lu exact=%lu",
            (unsigned long)t.cartesian_latitude_rejects,
            (unsigned long)t.cartesian_meridian_rejects,
            (unsigned long)t.cartesian_fallbacks,
            (unsigned long)t.prologue_row_rejects,
            (unsigned long)t.prologue_column_rejects,
            (unsigned long)t.edge_rejects, (unsigned long)t.visible_trails,
            (unsigned long)t.exact_gate_fallbacks);
    hs::log("msp counts render: dot=%lu long=%lu adaptive=%lu shader=%lu "
            "pal_end=%lu pal_lerp=%lu hole_early=%lu",
            (unsigned long)t.one_dot_edges, (unsigned long)t.long_edges,
            (unsigned long)t.adaptive_samples,
            (unsigned long)t.fragment_shader_calls,
            (unsigned long)t.palette_endpoints,
            (unsigned long)t.palette_interpolated,
            (unsigned long)t.hole_early_outs);
    hs::log("msp counts aa: tap0=%lu tap1=%lu tap2=%lu tap3=%lu tap4=%lu "
            "interior=%lu boundary=%lu",
            (unsigned long)t.aa_tap_masks[0], (unsigned long)t.aa_tap_masks[1],
            (unsigned long)t.aa_tap_masks[2], (unsigned long)t.aa_tap_masks[3],
            (unsigned long)t.aa_tap_masks[4], (unsigned long)t.interior_splats,
            (unsigned long)t.clip_boundary_splats);
    hs::g_mindsplatter_counts.reset();
  }
#endif

#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  static void log_mindsplatter_stall(const char *stage,
                                     const hs::DwtStallBucket &b) {
    char c0[21], c1[21], c2[21], c3[21];
    hs::log("msp stall: stage=%s batches=%lu cyc=%s cpi=%s lsu=%s exc=%s",
            stage, (unsigned long)b.batches, hs::u64_dec(b.cycles, c0),
            hs::u64_dec(b.cpi, c1), hs::u64_dec(b.lsu, c2),
            hs::u64_dec(b.exc, c3));
    if (b.wrapped)
      hs::log("msp stall aliased: stage=%s intervals=%lu/%lu ran 256+ cycles; "
              "cpi/lsu/exc understated",
              stage, (unsigned long)b.wrapped, (unsigned long)b.batches);
  }

  static void dump_mindsplatter_stalls() {
    const hs::MindSplatterStalls &t = hs::g_mindsplatter_stalls;
    log_mindsplatter_stall("history_vertex", t.history_vertex);
    log_mindsplatter_stall("trail_gate", t.trail_gate);
    log_mindsplatter_stall("edge_setup", t.edge_setup);
    log_mindsplatter_stall("adaptive_sim", t.adaptive_sim);
    log_mindsplatter_stall("normalized_replay", t.normalized_replay);
    log_mindsplatter_stall("shade_palette", t.shade_palette);
    log_mindsplatter_stall("projection", t.projection);
    log_mindsplatter_stall("aa_weights", t.aa_weights);
    log_mindsplatter_stall("framebuffer_blend", t.framebuffer_blend);
    log_mindsplatter_stall("signed_axis_physics", t.signed_axis_physics);
    hs::g_mindsplatter_stalls.reset();
  }
#endif

#ifdef HS_PROBE_BREAKDOWN
  /** @brief 64-bit window accumulators for the per-probe stage buckets. */
  struct ProbeTotals {
    uint64_t point = 0, project = 0, lut = 0, convex = 0, sector = 0, exact = 0;
    uint64_t pack = 0, alpha = 0, tick = 0;
    uint64_t n_probe = 0, n_cull_cos = 0, n_cull_r = 0, n_lut = 0, n_convex = 0;
    uint64_t n_sector = 0, n_exact = 0, n_alpha = 0;
    void reset() {
      point = project = lut = convex = sector = exact = pack = alpha = tick = 0;
      n_probe = n_cull_cos = n_cull_r = n_lut = n_convex = n_sector = n_exact =
          n_alpha = 0;
    }
  };

  /**
   * @brief Folds this frame's probe buckets into the window totals, rezeroed.
   * @details Same uint32 wrap argument as drain_scan_metrics: a frame of a
   *          many-faced solid accumulates enough cycles to wrap a bucket.
   */
  void drain_probe_breakdown() {
    const hs::ProbeBreakdown &b = hs::g_probe_breakdown;
    probe_totals.point += b.point;
    probe_totals.project += b.project;
    probe_totals.lut += b.edge_lut;
    probe_totals.convex += b.edge_convex;
    probe_totals.sector += b.edge_sector;
    probe_totals.exact += b.edge_exact;
    probe_totals.pack += b.pack;
    probe_totals.alpha += b.alpha;
    probe_totals.tick += b.tick;
    probe_totals.n_probe += b.n_probe;
    probe_totals.n_cull_cos += b.n_cull_cos;
    probe_totals.n_cull_r += b.n_cull_r;
    probe_totals.n_lut += b.n_lut;
    probe_totals.n_convex += b.n_convex;
    probe_totals.n_sector += b.n_sector;
    probe_totals.n_exact += b.n_exact;
    probe_totals.n_alpha += b.n_alpha;
    hs::g_probe_breakdown.reset();
  }

  /**
   * @brief Prints and resets this window's probe-stage cycle buckets.
   * @details Window totals in cycles beside their event counts, as the scan
   *          totals do. Each bucket carries one counter read; tick is the summed
   *          cost of a back-to-back read pair per probe, so tick/n_probe/2 is
   *          the per-read inflation to subtract from each bucket's mean.
   */
  void dump_probe_breakdown() {
    const ProbeTotals &t = probe_totals;
    char b0[21], b1[21], b2[21], b3[21], b4[21], b5[21], b6[21], b7[21], b8[21];
    hs::log("probe cycles: point=%s project=%s lut=%s convex=%s sector=%s "
            "exact=%s pack=%s alpha=%s tick=%s",
            hs::u64_dec(t.point, b0), hs::u64_dec(t.project, b1),
            hs::u64_dec(t.lut, b2), hs::u64_dec(t.convex, b3),
            hs::u64_dec(t.sector, b4), hs::u64_dec(t.exact, b5),
            hs::u64_dec(t.pack, b6), hs::u64_dec(t.alpha, b7),
            hs::u64_dec(t.tick, b8));
    hs::log("probe counts: probe=%s cull_cos=%s cull_r=%s lut=%s convex=%s "
            "sector=%s exact=%s alpha=%s",
            hs::u64_dec(t.n_probe, b0), hs::u64_dec(t.n_cull_cos, b1),
            hs::u64_dec(t.n_cull_r, b2), hs::u64_dec(t.n_lut, b3),
            hs::u64_dec(t.n_convex, b4), hs::u64_dec(t.n_sector, b5),
            hs::u64_dec(t.n_exact, b6), hs::u64_dec(t.n_alpha, b7));
    probe_totals.reset();
  }

  ProbeTotals probe_totals; /**< This window's drained probe buckets. */
#endif

  unsigned long total_frames =
      0; /**< Frames since this effect instance began. */
  unsigned long window_frames = 0; /**< Frames in the current readout window. */
  unsigned long wall_sum =
      0; /**< Summed draw_frame wall time this window (µs). */
  unsigned long wall_min = ~0ul; /**< Fastest draw_frame this window (µs). */
  unsigned long wall_max = 0;    /**< Slowest draw_frame this window (µs). */
  unsigned long render_sum =
      0; /**< Summed render (wall − sync wait) this window (µs). */
  unsigned long render_max = 0; /**< Slowest render this window (µs). */
  hs::CycleCounter *buffer_wait =
      nullptr; /**< The effect's *_buffer_wait counter. */
  unsigned long window_start = micros(); /**< Window wall-clock start (µs). */
#ifdef HS_MINDSPLATTER_REPLAY
#ifdef HS_MINDSPLATTER_REPLAY_AB
  mindsplatter_replay::ReferenceStats replay_stats;
#else
  mindsplatter_replay::FrameStats replay_stats;
#endif
  ClipRegion replay_clip{};
#endif
};

// Bytes the harness wrapper adds over the effect it profiles, derived so every
// conditionally-compiled member (replay stats, scan totals, probe buckets) is
// counted and the budget below stays a measurement of the effect alone.
static constexpr size_t PROFILE_WRAPPER_BYTES =
    sizeof(ProfiledEffect<CANVAS_W, CANVAS_H>) -
    sizeof(HS_PROFILE_TARGET<CANVAS_W, CANVAS_H>);
#ifdef HS_PROFILE_EFFECT_HEAP_BYTES
static constexpr size_t MAX_EFFECT_HEAP_BYTES = HS_PROFILE_EFFECT_HEAP_BYTES;
#else
static constexpr size_t MAX_EFFECT_HEAP_BYTES =
    HS_PHANTASM_EFFECT_HEAP_BYTES + PROFILE_WRAPPER_BYTES;
#endif

#ifdef HS_PROFILE_PRESET
template <typename E> void select_profile_preset(E &effect) {
  static_assert(
      requires(E &e) { e.profile_select_preset(size_t{}); },
      "HS_PROFILE_PRESET target does not support preset selection");
  effect.profile_select_preset(static_cast<size_t>(HS_PROFILE_PRESET));
}
#endif

Effect *construct_profiled() {
  using Target = ProfiledEffect<CANVAS_W, CANVAS_H>;
  auto *e =
      static_cast<Target *>(construct_effect<Target, MAX_EFFECT_HEAP_BYTES>());
#ifdef HS_PROFILE_PRESET
  select_profile_preset(*e);
#ifdef HS_PROFILE_AUTO_RESUME
  e->setAnimationsPaused(false);
#endif
#endif
#ifdef HS_PROFILE_TRANS_SPEED
  // Per-run knob (e.g. IslamicStars carousel speed-up so a single epoch walks the
  // whole shape roster). No-op for effects that don't register "Trans Speed".
  e->updateParameter("Trans Speed", (float)(HS_PROFILE_TRANS_SPEED));
#endif
  return e;
}

const POV::EffectFactory EFFECT_FACTORIES[] = {&construct_profiled};

// One hour at RPM. The epoch only rebuilds the same effect here, under the
// K-revolution commit budget the flywheel ISR traps on, and it resets the
// profile counters mid-capture; park the boundary past any capture.
// HS_PROFILE_EPOCH_REVS still overrides this (pov_segmented.h).
constexpr uint32_t PROFILE_REVOLUTIONS[] = {RPM * 60};

// Config::effect_revolutions must span the whole roster: valid() and
// revolutions_for_effect() index it over [0, effect_count).
static_assert(std::size(PROFILE_REVOLUTIONS) == std::size(EFFECT_FACTORIES));

// The same per-effect RNG identity the shipping playlist builds, derived from
// the profiled class rather than the ProfiledEffect wrapper. Without it the
// driver falls back to hs::epoch_seed(0) and the harness measures a different
// random realization than the show renders.
constexpr uint64_t PROFILE_SEEDS[] = {hs::stable_effect_seed(
    hs::stable_effect_id<HS_PROFILE_TARGET<CANVAS_W, CANVAS_H>>(
        HS_PROFILE_STR(HS_PROFILE_TARGET)))};

static_assert(std::size(PROFILE_SEEDS) == std::size(EFFECT_FACTORIES));

constexpr pov::sync::Config profile_config() {
  auto cfg = pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W, 1);
  cfg.effect_revolutions = PROFILE_REVOLUTIONS;
  return cfg;
}

static_assert(profile_config().valid() == nullptr,
              "Profile pov::sync::Config invariants violated");
} // namespace

FLASHMEM void setup() {
  POV::park_sync_out();
  boot_serial();
  hs::log("profile harness: effect=%s config=%s segments=%d rpm=%u f_cpu=%lu",
          HS_PROFILE_STR(HS_PROFILE_TARGET),
          HS_PROFILE_STR(HS_PROFILE_CONFIG_TAG), NUM_SEGMENTS, RPM,
          (unsigned long)F_CPU);
#ifdef HS_PROFILE_PULLBACK_TELEMETRY
  hs::log("Pullback arm: %s sha=%s", HS_PROFILE_STR(HS_PULLBACK_ARM),
          HS_PULLBACK_SHORT_SHA);
#endif
  log_reset_cause();
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  hs::enable_mindsplatter_stall_counters();
#endif
  create_pov();
  // A strapped-nonzero board comes up downstream: dark in ACQUIRE, waiting on a
  // sync wire nothing drives, so the whole capture runs black and counter-empty.
  HS_CHECK(POV::segment_index() == 0,
           "profile harness needs segment 0 (master); ID straps read %d",
           POV::segment_index());
}

void loop() {
  // Never returns: runs the single-entry playlist forever.
  g_pov->run_show(EFFECT_FACTORIES, &PROFILE_REVOLUTIONS, &PROFILE_SEEDS);
}
