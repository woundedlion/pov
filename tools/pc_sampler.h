/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pc_sampler.h
 * @brief Statistical program-counter sampler for the Profile harness: the
 *        exception-frame decoding and fixed-capacity PC histogram, plus on
 *        Teensy the GPT2 interrupt that feeds them.
 */
#pragma once

#include "core/platform/profiling.h"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>

namespace hs {

/** @brief EXC_RETURN bit set when the interrupted context stacked on PSP. */
inline constexpr uint32_t EXC_RETURN_PSP = 1u << 2;
/** @brief EXC_RETURN bit set when the exception returns to Thread mode. */
inline constexpr uint32_t EXC_RETURN_THREAD = 1u << 3;
/** @brief EXC_RETURN bit clear when the stacked frame carries FP state. */
inline constexpr uint32_t EXC_RETURN_BASIC_FRAME = 1u << 4;

/**
 * @brief Word index of the stacked return address in an exception frame.
 * @details The FP-extended frame appends its registers after xPSR, so the
 *          index is the same for the basic and the extended frame, whether or
 *          not lazy stacking has filled the FP slots yet.
 */
inline constexpr size_t STACKED_PC_WORD = 6;

/**
 * @brief Return address the hardware stacked for the interrupted context.
 * @param msp Main stack pointer at handler entry, before any push.
 * @param psp Process stack pointer at handler entry.
 * @param exc_return The handler's entry LR.
 * @return The address the interrupted context resumes at.
 */
inline uint32_t interrupted_pc(const uint32_t *msp, const uint32_t *psp,
                               uint32_t exc_return) {
  const uint32_t *frame = (exc_return & EXC_RETURN_PSP) ? psp : msp;
  return frame[STACKED_PC_WORD];
}

/**
 * @brief Timer periods that elapsed without a serviced tick.
 * @param elapsed_cycles Core cycles since the previous serviced tick.
 * @param period_cycles Nominal tick period in core cycles (nonzero).
 * @return Elapsed periods, rounded to nearest, beyond the one that fired.
 */
constexpr uint32_t pc_lost_ticks(uint32_t elapsed_cycles,
                                 uint32_t period_cycles) {
  const uint32_t half = period_cycles / 2;
  if (elapsed_cycles < period_cycles + half)
    return 0;
  const uint32_t whole = elapsed_cycles / period_cycles;
  const bool round_up = elapsed_cycles % period_cycles >= period_cycles - half;
  return round_up ? whole : whole - 1;
}

/**
 * @brief Open-addressing histogram of PCs, keyed by aligned address bin.
 * @tparam CAPACITY Slot count; a power of two of at least 2.
 * @tparam BIN_SHIFT log2 of the bin size in bytes.
 * @details Zero-initialized storage is the empty table, so a static instance
 *          needs no constructor run. Insertion probes linearly from a
 *          Fibonacci hash and gives up after MAX_PROBES slots.
 */
template <size_t CAPACITY, unsigned BIN_SHIFT> class PcHistogram {
  static_assert(CAPACITY >= 2 && (CAPACITY & (CAPACITY - 1)) == 0,
                "PcHistogram capacity must be a power of two >= 2");
  static_assert(BIN_SHIFT >= 1 && BIN_SHIFT <= 16,
                "PcHistogram bins span 2 to 64 KiB");

public:
  static constexpr uint32_t BIN_BYTES = 1u << BIN_SHIFT; /**< Bin size. */
  static constexpr size_t MAX_PROBES =
      CAPACITY < 64 ? CAPACITY : 64; /**< Slots tried per insert. */

  /**
   * @brief Counts one sample at @p pc.
   * @param pc Sampled address.
   * @return False when no slot was free within MAX_PROBES.
   */
  bool add(uint32_t pc) {
    const uint32_t key = OCCUPIED | (pc >> BIN_SHIFT);
    size_t i = (size_t)((key * 0x9E3779B1u) >> (32 - LOG2_CAPACITY));
    for (size_t n = 0; n < MAX_PROBES; ++n, i = (i + 1) & (CAPACITY - 1)) {
      Slot &s = slots[i];
      if (s.key == key) {
        ++s.count;
        return true;
      }
      if (s.key == 0) {
        s.key = key;
        s.count = 1;
        ++used;
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Visits every occupied bin in slot order.
   * @param f Called as f(bin_start_address, count).
   */
  template <typename F> void for_each(F &&f) const {
    for (const Slot &s : slots)
      if (s.key)
        f((s.key & ~OCCUPIED) << BIN_SHIFT, s.count);
  }

  /** @brief Occupied slots. @return The distinct bins counted. */
  size_t size() const { return used; }

private:
  static constexpr uint32_t OCCUPIED = 0x80000000u;
  static constexpr unsigned LOG2_CAPACITY = std::countr_zero(CAPACITY);

  /** @brief One bin: its tagged key and its sample count. */
  struct Slot {
    uint32_t key;   /**< OCCUPIED | (pc >> BIN_SHIFT); 0 marks a free slot. */
    uint32_t count; /**< Samples in the bin. */
  };

  Slot slots[CAPACITY]{};
  size_t used = 0;
};

/** @brief Number of log2 latency buckets PcSampleCounters keeps. */
inline constexpr size_t PC_LATENCY_BUCKETS = 12;

/**
 * @brief Log2 bucket for a tick's service latency.
 * @param latency Timer counts between the compare match and the read.
 * @return 0 for 0, else bit_width(latency), capped at the last bucket.
 */
constexpr size_t pc_latency_bucket(uint32_t latency) {
  const size_t b = (size_t)std::bit_width(latency);
  return b < PC_LATENCY_BUCKETS ? b : PC_LATENCY_BUCKETS - 1;
}

/** @brief Tick classification totals. */
struct PcSampleCounters {
  uint32_t ticks;    /**< Ticks serviced. */
  uint32_t idle;     /**< Ticks, and their lost periods, with the gate shut. */
  uint32_t deferred; /**< Open-gate ticks serviced late or after a loss. */
  uint32_t missed;   /**< Lost periods charged to an open gate. */
  uint32_t handler;  /**< Binned ticks that interrupted a handler. */
  uint32_t recorded; /**< Open-gate ticks binned by PC. */
  uint32_t dropped;  /**< Open-gate ticks the histogram had no slot for. */
  uint32_t latency[PC_LATENCY_BUCKETS]; /**< Ticks by pc_latency_bucket. */
  uint64_t service_cycles;              /**< Summed handler body cycles. */
  uint32_t service_max;                 /**< Longest handler body in cycles. */
};

/** @brief Fixed sampler parameters echoed in the dump header. */
struct PcSampleConfig {
  uint32_t hz;            /**< Nominal tick rate. */
  uint32_t defer_counts;  /**< Latency above which a tick is deferred. */
  uint32_t counts_per_us; /**< Timer counts per microsecond. */
};

/**
 * @brief Classifies timer ticks against the render gate and bins the PCs of
 *        the ones taken promptly inside it.
 * @tparam CAPACITY Histogram slot count.
 * @tparam BIN_SHIFT log2 of the histogram bin size in bytes.
 * @details A late tick's interrupted PC is where the thread resumed after a
 *          higher-priority handler or a masked section, not where the time
 *          went, so it is counted as deferred and not binned.
 */
template <size_t CAPACITY, unsigned BIN_SHIFT> class PcSampler {
public:
  /** @brief Histogram type the sampler bins into. */
  using Histogram = PcHistogram<CAPACITY, BIN_SHIFT>;

  /**
   * @brief Sets the parameters the classification and dump use.
   * @param c Rate, deferral threshold and timer scale.
   */
  void configure(const PcSampleConfig &c) { config = c; }

  /**
   * @brief Classifies one serviced tick.
   * @param open Whether the render gate was open at service.
   * @param latency Timer counts from the compare match to service.
   * @param lost Periods that elapsed with no serviced tick.
   * @param exc_return The handler's entry LR.
   * @param pc The interrupted context's return address.
   */
  void tick(bool open, uint32_t latency, uint32_t lost, uint32_t exc_return,
            uint32_t pc) {
    ++counters.ticks;
    ++counters.latency[pc_latency_bucket(latency)];
    if (!open) {
      counters.idle += 1 + lost;
      return;
    }
    if (lost || latency > config.defer_counts) {
      counters.missed += lost;
      ++counters.deferred;
      return;
    }
    if (!histogram.add(pc)) {
      ++counters.dropped;
      return;
    }
    ++counters.recorded;
    if (!(exc_return & EXC_RETURN_THREAD))
      ++counters.handler;
  }

  /**
   * @brief Accumulates one handler body's cycle cost.
   * @param cycles Cycles from handler body entry to exit.
   */
  void note_service(uint32_t cycles) {
    counters.service_cycles += cycles;
    if (cycles > counters.service_max)
      counters.service_max = cycles;
  }

  /** @brief Tick totals. @return The counters. */
  const PcSampleCounters &totals() const { return counters; }

  /** @brief Binned samples. @return The histogram. */
  const Histogram &bins() const { return histogram; }

  /**
   * @brief Writes the dump `tools/pc_profile.py` reads, one line per call.
   * @param emit Called with each NUL-terminated line, no newline.
   * @param elapsed_ms Wall time the sampler ran.
   */
  template <typename Emit> void dump(Emit &&emit, uint32_t elapsed_ms) const {
    const PcSampleCounters &c = counters;
    char line[192];
    std::snprintf(line, sizeof(line),
                  "pc profile: version=1 hz=%lu bin_bytes=%lu capacity=%lu "
                  "defer_counts=%lu counts_per_us=%lu elapsed_ms=%lu",
                  (unsigned long)config.hz, (unsigned long)Histogram::BIN_BYTES,
                  (unsigned long)CAPACITY, (unsigned long)config.defer_counts,
                  (unsigned long)config.counts_per_us,
                  (unsigned long)elapsed_ms);
    emit(line);
    const uint32_t service_avg =
        c.ticks ? (uint32_t)(c.service_cycles / c.ticks) : 0;
    std::snprintf(line, sizeof(line),
                  "pc profile: ticks=%lu idle=%lu deferred=%lu missed=%lu "
                  "handler=%lu recorded=%lu dropped=%lu bins=%lu "
                  "service_avg=%lu service_max=%lu",
                  (unsigned long)c.ticks, (unsigned long)c.idle,
                  (unsigned long)c.deferred, (unsigned long)c.missed,
                  (unsigned long)c.handler, (unsigned long)c.recorded,
                  (unsigned long)c.dropped, (unsigned long)histogram.size(),
                  (unsigned long)service_avg, (unsigned long)c.service_max);
    emit(line);
    int n = std::snprintf(line, sizeof(line), "pc profile latency:");
    for (size_t b = 0; b < PC_LATENCY_BUCKETS; ++b)
      n += std::snprintf(line + n, sizeof(line) - (size_t)n, " %lu",
                         (unsigned long)c.latency[b]);
    emit(line);
    histogram.for_each([&](uint32_t address, uint32_t count) {
      std::snprintf(line, sizeof(line), "pc %08lx %lu", (unsigned long)address,
                    (unsigned long)count);
      emit(line);
    });
    emit("pc profile end");
  }

private:
  Histogram histogram;
  PcSampleCounters counters{};
  PcSampleConfig config{};
};

} // namespace hs

#ifdef CORE_TEENSY
#include <Arduino.h>

#ifndef HS_PROFILE_PC_HZ
/** @brief PC sample rate in Hz. */
#define HS_PROFILE_PC_HZ 10000
#endif
#ifndef HS_PROFILE_PC_BINS
/** @brief PC histogram slot count; a power of two. */
#define HS_PROFILE_PC_BINS 4096
#endif
#ifndef HS_PROFILE_PC_BIN_SHIFT
/** @brief log2 of the PC histogram bin size in bytes. */
#define HS_PROFILE_PC_BIN_SHIFT 4
#endif
#ifndef HS_PROFILE_PC_MS
/** @brief Milliseconds from sampler start to the one-shot dump. */
#define HS_PROFILE_PC_MS 30000
#endif
#ifndef HS_PROFILE_PC_DEFER_US
/** @brief Service latency in microseconds past which a tick is deferred. */
#define HS_PROFILE_PC_DEFER_US 2
#endif

namespace hs {

/** @brief GPT2 input clock: PERCLK, which the Teensy core runs at 24 MHz. */
inline constexpr uint32_t PC_TIMER_HZ = 24000000;
/** @brief NVIC priority of the sampler: below the core's default 128. */
inline constexpr uint8_t PC_SAMPLER_IRQ_PRIORITY = 192;

static_assert(PC_TIMER_HZ % HS_PROFILE_PC_HZ == 0 && HS_PROFILE_PC_HZ >= 1000 &&
                  HS_PROFILE_PC_HZ <= 50000,
              "HS_PROFILE_PC_HZ must divide 24 MHz within 1-50 kHz");

/** @brief The sampler's histogram and counters (zero-initialized .bss). */
inline PcSampler<HS_PROFILE_PC_BINS, HS_PROFILE_PC_BIN_SHIFT> g_pc_sampler;
/** @brief Cycle counter at the previous serviced tick. */
inline uint32_t g_pc_last_tick_cycles;
/** @brief Tick period in core cycles. */
inline uint32_t g_pc_period_cycles;
/** @brief millis() when the sampler started. */
inline uint32_t g_pc_start_ms;

} // namespace hs

/**
 * @brief Body of the GPT2 handler; reached by tail branch from
 *        hs_pc_sampler_isr with LR still holding EXC_RETURN.
 * @param msp Main stack pointer at handler entry.
 * @param psp Process stack pointer at handler entry.
 * @param exc_return The handler's entry LR.
 */
extern "C" inline __attribute__((used)) void
hs_pc_sampler_tick(const uint32_t *msp, const uint32_t *psp,
                   uint32_t exc_return) {
  const uint32_t start = ARM_DWT_CYCCNT;
  const uint32_t latency = GPT2_CNT;
  GPT2_SR = GPT_SR_OF1;
  const uint32_t lost = hs::pc_lost_ticks(start - hs::g_pc_last_tick_cycles,
                                          hs::g_pc_period_cycles);
  hs::g_pc_last_tick_cycles = start;
  hs::g_pc_sampler.tick(hs::pc_sample_gate() != 0, latency, lost, exc_return,
                        hs::interrupted_pc(msp, psp, exc_return));
  hs::g_pc_sampler.note_service(ARM_DWT_CYCCNT - start);
  // The status clear is a posted write; retire it before returning.
  __asm__ volatile("dsb" ::: "memory");
}

/**
 * @brief GPT2 vector: captures both stack pointers before any push and
 *        tail-branches to hs_pc_sampler_tick.
 */
extern "C" inline __attribute__((naked, used)) void hs_pc_sampler_isr() {
  __asm__ volatile("mrs r0, msp\n\t"
                   "mrs r1, psp\n\t"
                   "mov r2, lr\n\t"
                   "b hs_pc_sampler_tick\n\t");
}

namespace hs {

/** @brief Starts GPT2 ticking at HS_PROFILE_PC_HZ in restart mode. */
HS_FLASH_INLINE inline void pc_sampler_begin() {
  g_pc_sampler.configure({HS_PROFILE_PC_HZ,
                          HS_PROFILE_PC_DEFER_US * (PC_TIMER_HZ / 1000000),
                          PC_TIMER_HZ / 1000000});
  CCM_CCGR0 |=
      CCM_CCGR0_GPT2_BUS(CCM_CCGR_ON) | CCM_CCGR0_GPT2_SERIAL(CCM_CCGR_ON);
  GPT2_CR = 0;
  GPT2_PR = 0;
  GPT2_SR = 0x3F;
  GPT2_IR = GPT_IR_OF1IE;
  GPT2_OCR1 = PC_TIMER_HZ / HS_PROFILE_PC_HZ - 1;
  GPT2_CR = GPT_CR_CLKSRC(1) | GPT_CR_ENMOD;
  attachInterruptVector(IRQ_GPT2, hs_pc_sampler_isr);
  NVIC_SET_PRIORITY(IRQ_GPT2, PC_SAMPLER_IRQ_PRIORITY);
  g_pc_period_cycles = F_CPU_ACTUAL / HS_PROFILE_PC_HZ;
  g_pc_start_ms = millis();
  g_pc_last_tick_cycles = ARM_DWT_CYCCNT;
  NVIC_ENABLE_IRQ(IRQ_GPT2);
  GPT2_CR |= GPT_CR_EN;
}

/**
 * @brief Once @p after_ms have passed since pc_sampler_begin(), stops the
 *        timer and logs the dump; later calls do nothing.
 * @param after_ms Sampling duration in milliseconds.
 */
HS_FLASH_INLINE inline void pc_sampler_dump_when_due(uint32_t after_ms) {
  static bool done = false;
  const uint32_t elapsed = millis() - g_pc_start_ms;
  if (done || elapsed < after_ms)
    return;
  done = true;
  NVIC_DISABLE_IRQ(IRQ_GPT2);
  GPT2_CR &= ~GPT_CR_EN;
  GPT2_SR = 0x3F;
  g_pc_sampler.dump([](const char *line) { hs::log("%s", line); }, elapsed);
}

} // namespace hs
#endif
