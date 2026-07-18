/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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

// Select the DMA HD107S output path (same guard dance as Phantasm.ino: the
// profile env also passes -D USE_DMA_LEDS).
#ifndef USE_DMA_LEDS
#define USE_DMA_LEDS
#endif

#ifndef PHANTASM_NUM_SEGMENTS
#define PHANTASM_NUM_SEGMENTS 4
#endif

// Effect class to profile; overridden per-run by `just profile <EffectClass>`.
#ifndef HS_PROFILE_TARGET
#define HS_PROFILE_TARGET DisplacementField
#endif

// Frames per readout window; override (via PLATFORMIO_BUILD_FLAGS) for finer
// phase resolution at the cost of more serial traffic.
#ifndef HS_PROFILE_WINDOW
#define HS_PROFILE_WINDOW 128
#endif

// Further per-run knobs consumed elsewhere (all via PLATFORMIO_BUILD_FLAGS):
//   HS_PROFILE_EPOCH_REVS    epoch length override (pov_segmented.h) so one
//                            instance covers a full preset cycle
//   HS_PROFILE_ORDERED_CYCLE random-next cyclers advance in order instead
//                            (HankinSolids, SphericalHarmonics)
//   HS_PROFILE_TRANS_SPEED   "Trans Speed" applied after init (below)
//   HS_SCAN_METRICS          compiles in the per-pixel hs::g_scan_metrics probe
//                            counters and adds the "scan totals" window line.
//                            Every probe pays a non-atomic global increment, so
//                            read the COUNTS from such a build, not the times.

#define HS_PROFILE_STR2(x) #x
#define HS_PROFILE_STR(x) HS_PROFILE_STR2(x)

#include <FastLED.h>
#include <SPI.h>
#include <new> // std::nothrow — fail-fast OOM check on the POV allocation below

#include "pov_segmented.h"
#include "engine/effects.h"

static constexpr int TOTAL_PIXELS = 288;
static constexpr int NUM_SEGMENTS = PHANTASM_NUM_SEGMENTS;
static constexpr unsigned int RPM = 480;

using POV = POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>;

#if defined(USE_DMA_LEDS)
// Out-of-line definition for this target's controller, emitted as the required
// DMAMEM explicit specialization (see pov_segmented.h).
HS_DEFINE_POV_SEGMENTED_LED_CONTROLLER(TOTAL_PIXELS, NUM_SEGMENTS, RPM);
#endif

namespace {

/**
 * @brief Wraps the profiled effect: roots the HS_PROFILE tree at one `frame`
 * counter per draw_frame() and accumulates exact wall-clock stats, dumping
 * both every WINDOW_FRAMES frames.
 */
template <int W, int H> class ProfiledEffect : public HS_PROFILE_TARGET<W, H> {
public:
  void draw_frame() override {
    const uint64_t bw0 = buffer_wait_ ? buffer_wait_->cycles : 0;
    const unsigned long t0 = micros();
    {
      HS_PROFILE(frame);
      HS_PROFILE_TARGET<W, H>::draw_frame();
    }
    const unsigned long dt = micros() - t0;
    // Per-frame render = wall minus this frame's display-sync wait, read as
    // the effect's *_buffer_wait counter delta. The counter self-registers on
    // the first draw_frame, so the lookup retries until it appears.
    if (!buffer_wait_)
      buffer_wait_ = hs::CycleCounter::find_suffix("_buffer_wait");
    const uint64_t bw1 = buffer_wait_ ? buffer_wait_->cycles : 0;
    const unsigned long wait_us = (unsigned long)(
        (bw1 - bw0) / hs::CycleCounter::CYCLES_PER_US);
    const unsigned long render = dt > wait_us ? dt - wait_us : 0;
    // One compact line per frame (the full counter tree still dumps per
    // window — per-frame log_all would perturb the frames it measures).
    hs::log("f %lu w=%lu r=%lu", total_frames_ + 1, dt, render);
#ifdef HS_SCAN_METRICS
    drain_scan_metrics();
#endif
    render_sum_ += render;
    if (render > render_max_) render_max_ = render;
    wall_sum_ += dt;
    if (dt < wall_min_) wall_min_ = dt;
    if (dt > wall_max_) wall_max_ = dt;
    ++total_frames_;
    if (++window_frames_ == WINDOW_FRAMES) dump();
  }

private:
  static constexpr int WINDOW_FRAMES = HS_PROFILE_WINDOW; /**< Frames per readout window. */

  void dump() {
    const unsigned long now = micros();
    hs::log("=== profile %s [%dx%d] frames %lu-%lu window=%lu us ===",
            HS_PROFILE_STR(HS_PROFILE_TARGET), W, H,
            total_frames_ - window_frames_ + 1, total_frames_,
            now - window_start_);
    hs::log("frame wall us: min=%lu avg=%lu max=%lu sum=%lu (%d frames)",
            wall_min_, wall_sum_ / WINDOW_FRAMES, wall_max_, wall_sum_,
            WINDOW_FRAMES);
    hs::log("frame render us: avg=%lu max=%lu",
            render_sum_ / WINDOW_FRAMES, render_max_);
    dump_isr_stats(now - window_start_);
    hs::CycleCounter::log_all();
#ifdef HS_SCAN_METRICS
    dump_scan_totals();
#endif
    hs::CycleCounter::reset_all();
    window_frames_ = 0;
    wall_sum_ = 0;
    wall_min_ = ~0ul;
    wall_max_ = 0;
    render_sum_ = 0;
    render_max_ = 0;
    window_start_ = micros();
  }

  /**
   * @brief Prints and resets the column-ISR accumulators for this window.
   * @param window_us Window wall-clock span, for the CPU-share figure.
   * @details Copy + reset under a brief IRQ-off window (the ISR is the sole
   * writer). Per-call figures are exact cycles; ns = cycles * 5 / 3 at 600 MHz.
   */
  static void dump_isr_stats(unsigned long window_us) {
    hs::IsrCycleStats wake, pack, submit;
    __disable_irq();
    wake = hs::g_flywheel_wake_cycles;
    pack = hs::g_column_pack_cycles;
    submit = hs::g_dma_submit_cycles;
    hs::g_flywheel_wake_cycles.reset();
    hs::g_column_pack_cycles.reset();
    hs::g_dma_submit_cycles.reset();
    __enable_irq();
    log_isr("isr_wake", wake, window_us);
    log_isr("isr_pack", pack, window_us);
    log_isr("isr_dma_submit", submit, window_us);
  }

  static void log_isr(const char *name, const hs::IsrCycleStats &s,
                      unsigned long window_us) {
    if (!s.count) {
      hs::log("%-14s n=0", name);
      return;
    }
    const uint32_t avg = (uint32_t)(s.cycles / s.count);
    const uint32_t total_us = (uint32_t)(s.cycles / 600u);
    // CPU share in hundredths of a percent; window_us is far below the
    // 32-bit ceiling of total_us * 10000.
    const uint32_t share_c =
        window_us ? (uint32_t)((uint64_t)total_us * 10000u / window_us) : 0;
    hs::log("%-14s n=%lu cyc min/avg/max=%lu/%lu/%lu "
            "ns min/avg/max=%lu/%lu/%lu total=%lu us cpu=%lu.%02lu%%",
            name, (unsigned long)s.count, (unsigned long)s.min,
            (unsigned long)avg, (unsigned long)s.max,
            (unsigned long)(s.min * 5u / 3u), (unsigned long)(avg * 5u / 3u),
            (unsigned long)(s.max * 5u / 3u), (unsigned long)total_us,
            (unsigned long)(share_c / 100u), (unsigned long)(share_c % 100u));
  }

#ifdef HS_SCAN_METRICS
  /** @brief 64-bit window accumulators for the per-pixel scan counters. */
  struct ScanTotals {
    uint64_t tested = 0;    /**< Face::distance probes. */
    uint64_t culled = 0;    /**< Probes rejected by the back-face / radius guards. */
    uint64_t exact = 0;     /**< Probes taking a full evaluation (convex + sector + walk). */
    uint64_t convex = 0;    /**< Full evaluations on the convex half-plane path. */
    uint64_t sector = 0;    /**< Full evaluations on the concave sector walk. */
    uint64_t lut = 0;       /**< Probes served by the class-LUT bilinear fetch. */
    uint64_t cand = 0;      /**< Pixels passing the scan's d < pixel_width test. */
    uint64_t backstop = 0;  /**< plot() steps_cache capacity-backstop trips. */
    void reset() { tested = culled = exact = convex = sector = lut = cand = backstop = 0; }
  };

  /**
   * @brief Folds this frame's scan counters into the window totals, rezeroed.
   * @details The source counters are uint32; a long window of a many-faced
   *          solid probes enough pixels to wrap one, so they drain per frame.
   */
  void drain_scan_metrics() {
    const hs::ScanMetrics &m = hs::g_scan_metrics;
    scan_totals_.tested += m.pixels_tested;
    scan_totals_.culled += m.pixels_culled;
    scan_totals_.exact += m.exact_hits;
    scan_totals_.convex += m.convex_hits;
    scan_totals_.sector += m.sector_hits;
    scan_totals_.lut += m.lut_hits;
    scan_totals_.cand += m.shade_candidates;
    scan_totals_.backstop += m.plot_backstop_hits;
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
    const ScanTotals &t = scan_totals_;
    const uint64_t walk = t.exact - t.convex - t.sector;
    char b0[21], b1[21], b2[21], b3[21], b4[21], b5[21], b6[21], b7[21];
    hs::log("scan totals: tested=%s culled=%s lut=%s convex=%s sector=%s "
            "walk=%s cand=%s backstop=%s",
            hs::u64_dec(t.tested, b0), hs::u64_dec(t.culled, b1),
            hs::u64_dec(t.lut, b2), hs::u64_dec(t.convex, b3),
            hs::u64_dec(t.sector, b4), hs::u64_dec(walk, b5),
            hs::u64_dec(t.cand, b6), hs::u64_dec(t.backstop, b7));
    scan_totals_.reset();
  }

  ScanTotals scan_totals_; /**< This window's drained scan counters. */
#endif

  unsigned long total_frames_ = 0;  /**< Frames since this effect instance began. */
  unsigned long window_frames_ = 0; /**< Frames in the current readout window. */
  unsigned long wall_sum_ = 0;      /**< Summed draw_frame wall time this window (µs). */
  unsigned long wall_min_ = ~0ul;   /**< Fastest draw_frame this window (µs). */
  unsigned long wall_max_ = 0;      /**< Slowest draw_frame this window (µs). */
  unsigned long render_sum_ = 0;    /**< Summed render (wall − sync wait) this window (µs). */
  unsigned long render_max_ = 0;    /**< Slowest render this window (µs). */
  hs::CycleCounter* buffer_wait_ = nullptr; /**< The effect's *_buffer_wait counter. */
  unsigned long window_start_ = micros(); /**< Window wall-clock start (µs). */
};

POV *g_pov;  // g_-prefixed: a bare `pov` collides with the hardware `namespace pov`

// Slightly above Phantasm's shipping per-effect budget: the wrapper adds its
// profiling bookkeeping on top of the wrapped effect.
static constexpr size_t MAX_EFFECT_HEAP_BYTES = 3584 + 64;

Effect *construct_profiled() {
  using E = ProfiledEffect<288, 144>;
  static_assert(sizeof(E) <= MAX_EFFECT_HEAP_BYTES,
                "profiled effect exceeds the heap-object budget");
  // Eager-fill the scanline LUTs before the first frame so the flywheel ISR
  // never observes a half-filled table.
  GeometryResolution<E>::init();
  E *e = new (std::nothrow) E();
  HS_CHECK(e != nullptr, "effect allocation failed (OOM)");
  configure_arenas_default(); // Reset before init so effects can override
  e->init();
#ifdef HS_PROFILE_TRANS_SPEED
  // Per-run knob (e.g. IslamicStars carousel speed-up so a single epoch walks the
  // whole shape roster). No-op for effects that don't register "Trans Speed".
  e->updateParameter("Trans Speed", (float)(HS_PROFILE_TRANS_SPEED));
#endif
  return e;
}

const POV::EffectFactory EFFECT_FACTORIES[] = {&construct_profiled};

static_assert(pov::sync::phantasm_config(F_CPU, RPM, CANVAS_W, 1).valid(),
              "Profile pov::sync::Config invariants violated");
} // namespace

void setup() {
  Serial.begin(9600); // baud inert on Teensy USB-CDC; initializes Serial only
  delay(1000);        // USB-CDC enumeration settle so early output isn't lost
  hs::log("profile harness: effect=%s segments=%d rpm=%u f_cpu=%lu",
          HS_PROFILE_STR(HS_PROFILE_TARGET), NUM_SEGMENTS, RPM,
          (unsigned long)F_CPU);
  g_pov = new (std::nothrow) POV();
  HS_CHECK(g_pov != nullptr, "POV allocation failed (OOM)");
}

void loop() {
  // Never returns: runs the single-entry playlist forever.
  g_pov->run_show(EFFECT_FACTORIES, 1);
}
