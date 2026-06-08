---
name: Teensy Hardware Cycle Profiling
description: Instrument Holosphere effects with ARM DWT cycle counter profiling for the Teensy 4.0 (600MHz Cortex-M7)
---

# Teensy Hardware Cycle Profiling

## Overview

Use the ARM DWT (Data Watchpoint and Trace) hardware cycle counter to measure exact cycle costs of code sections on the Teensy 4.0. The counter is free-running at 600MHz (1 cycle = 1.667ns, 600 cycles = 1µs).

## Platform Macro

The project defines a cross-platform macro in `core/platform.h`:

```cpp
#ifdef CORE_TEENSY
#define HS_OS_CYCLES() ARM_DWT_CYCCNT
#else
#define HS_OS_CYCLES() 0
#endif
```

**Always use `HS_OS_CYCLES()`** instead of `ARM_DWT_CYCCNT` directly. This ensures the code compiles on both Teensy and WASM targets (WASM returns 0).

## Logging

`hs::log()` in `core/platform.h` supports printf-style format strings on Teensy via `vsnprintf`. Use `%lu` for `uint32_t` values.

```cpp
hs::log("PROF phase_a:%lu phase_b:%lu", cycles_a / 600, cycles_b / 600);
```

Division by 600 converts cycles to microseconds for readability.

## Instrumentation Patterns

### Pattern 1: Phase Timing in draw_frame()

Wrap each major phase with cycle counter reads:

```cpp
void draw_frame() override {
    Canvas canvas(*this);
    
    uint32_t c0 = HS_OS_CYCLES();
    timeline.step(canvas);
    uint32_t c1 = HS_OS_CYCLES();
    
    // ... more phases ...
    uint32_t c2 = HS_OS_CYCLES();
    
    hs::log("PROF timeline:%lu draw:%lu", (c1-c0)/600, (c2-c1)/600);
}
```

### Pattern 2: Accumulating Lambda Costs

For lambdas called many times per frame (shaders, callbacks), accumulate cycle counts into member variables:

```cpp
// Member variables
uint32_t prof_vs_cycles_ = 0, prof_vs_calls_ = 0;

// In draw_frame, reset before use:
prof_vs_cycles_ = 0; prof_vs_calls_ = 0;

// In lambda:
auto vertex_shader = [&](Fragment &f) {
    uint32_t t0 = HS_OS_CYCLES();
    // ... work ...
    prof_vs_cycles_ += HS_OS_CYCLES() - t0;
    prof_vs_calls_++;
};

// After: compute per-call average
// prof_vs_cycles_ / prof_vs_calls_ = avg cycles per call
```

### Pattern 3: Global Counters for Core Functions

For profiling inside shared core code (plot.h, scan.h), use global counters:

```cpp
// In a header, before namespace:
struct MyProfile {
    uint32_t phase_a = 0;
    uint32_t phase_b = 0;
    uint32_t count = 0;
    void reset() { phase_a = phase_b = count = 0; }
};
inline MyProfile g_my_prof;

// In hot path:
uint32_t t0 = HS_OS_CYCLES();
// ... work ...
g_my_prof.phase_a += HS_OS_CYCLES() - t0;
g_my_prof.count++;

// In effect's draw_frame (reset before, log after):
g_my_prof.reset();
// ... call hot path ...
hs::log("PROF a:%lu(%lu)", g_my_prof.phase_a/600, g_my_prof.count);
```

## Overhead

Each `HS_OS_CYCLES()` read is a single register load (~1 cycle). For tight inner loops with 100K+ iterations, the accumulation overhead is ~0.3ms. Acceptable for profiling.

## Converting Results

| Cycles | Microseconds | Milliseconds |
|--------|-------------|-------------|
| 600 | 1 µs | 0.001 ms |
| 60,000 | 100 µs | 0.1 ms |
| 600,000 | 1,000 µs | 1 ms |
| 60,000,000 | 100,000 µs | 100 ms |

## Existing Effect Profiling

Some effects have their own `OS_CYCLES()` macros (e.g., `IslamicStars.h`, `HankinSolids.h`). These are equivalent to `HS_OS_CYCLES()`. Prefer the global `HS_OS_CYCLES()` for consistency.

## Cleanup Checklist

After profiling, always remove:
- [ ] Cycle counter reads (`HS_OS_CYCLES()` calls)
- [ ] Profiling member variables  
- [ ] Global profiling structs
- [ ] `hs::log()` profiling output
- [ ] Any accumulated counter variables

## Key Learnings from Past Profiling

### MindSplatter (290ms → 64ms)
- **Pipeline tween amplification**: Orient filter's `tween(orientation, ...)` multiplied ALL downstream work by N orientation frames. Removing it saved 47ms.
- **Sub-pixel overdraw**: 89% of rasterize segments were sub-pixel, wasting full pipeline.plot cost per segment.
- **pipeline.plot cost**: ~1,500-1,900 cycles per call (vector_to_pixel + AntiAlias 4× writes + DMAMEM blend). This is the dominant per-pixel cost.
- **Baked palette**: Replacing virtual `GenerativePalette::get()` with LUT lookup saved 130ms.
- **Fragment shader** is cheap (~115 cycles) when using baked palettes.
- **Physics** is cheap (~2ms for 800 particles) — not worth optimizing.
