/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

// platform.h defines NDEBUG on device; include before <cassert> so assert
// stripping does not depend on include order.
#include "platform/platform.h"
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <utility>
#include "platform/constants.h"
#include "control/preset_host.h"
#include "render/clip.h"
#include "color/color.h"
#if HS_ENABLE_TEST_HOOKS
#include <cstdlib>
#endif

/**
 * @file canvas.h
 * @brief EffectConfig, the Effect base class, and the Canvas pixel buffer.
 * @details Effect's parameter registry and preset controller live in
 * control/param_host.h and control/preset_host.h.
 */

class Canvas;

/**
 * @brief Construction-time flags for an Effect (see the accessors of the same
 *        name). Defaults suit a plain, non-strobing, band-clippable effect.
 * @details pipeline_config (render/filter/pipeline.h) folds a filter
 * pipeline's segment traits into this config.
 */
struct EffectConfig {
  bool strobe = false;  /**< POV column strobe (Effect::strobe_columns). */
  bool persist = false; /**< Copy previous frame forward (persists_pixels). */
  bool full_frame = false; /**< Force full-canvas render (needs_full_frame). */
  bool reads_outside_band =
      false; /**< Frame generation samples pixels outside the display band. */
  /**
   * @brief Render-bound expansion the effect's filters need, in pixels
   *        (ClipRegion::margin). Raised to the ClipRegion default when lower.
   */
  int margin = ClipRegion{}.margin;
};

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for
 * drawing a frame. The parameter registry (ParamHost) and preset controller
 * (PresetHost) come in through the base chain.
 */
class Effect : public PresetHost {
  friend class Canvas;

public:
  using BufferReadyHook = void (*)(Effect &);

  bool debug_visuals = false; /**< Flag to enable visual debugging overlays. */
#if HS_ENABLE_TEST_HOOKS
  /** @brief Forces the whole-buffer clear, for clip-clear parity tests. */
  bool force_full_buffer_clear = false;
#endif

  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect, in [1, MAX_W].
   * @param H The height (resolution) of the effect, in [1, MAX_H].
   * @param cfg Construction-time behavior flags; see EffectConfig. Defaults to
   *        a plain band-clippable effect.
   */
  HS_COLD_MEMBER Effect(int W, int H, EffectConfig cfg = {})
      : persist_pixels(cfg.persist), full_frame(cfg.full_frame),
        reads_outside_band(cfg.reads_outside_band), strobe(cfg.strobe),
        frame_width(W), frame_height(H) {
    HS_CHECK(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H,
             "Effect dimensions %d x %d are outside 1..%d x 1..%d", W, H, MAX_W,
             MAX_H);
    // Single-live-Effect precondition: every Effect aliases the same two static
    // buffers, so a second live instance corrupts both frames.
    HS_CHECK(
        !s_alive,
        "Effect: a second Effect was constructed while one is still alive; "
        "buffer_a/buffer_b are shared static storage (one live Effect only)");
    s_alive = true;
    // Point bufs at the shared static storage and clear both buffers.
    clear_buffers();
    clip_region.w = W;
    clip_region.h = H;
    clip_region.y_end = H;
    clip_region.x_end = W;
    // Only widen: a filter fold of 0 must not shrink the coverage every effect
    // already gets from the ClipRegion default. Routed through set_margin so an
    // unrepresentable folded requirement traps instead of being weakened.
    const int widened =
        cfg.margin > clip_region.margin ? cfg.margin : clip_region.margin;
    // A one-column canvas is already fully covered and has no representable
    // positive cylindrical margin.
    set_margin(W == 1 ? 0 : widened);
  }

  /**
   * @brief Destroys the Effect instance.
   * @details Clears the single-live-Effect guard so the next construction is
   * admitted; every effect-swap path destroys the outgoing instance first.
   */
  virtual ~Effect() { s_alive = false; }

  /**
   * @brief Post-construction initialization. Override to move heavy
   * setup logic here (avoids GCC C1/C2 constructor duplication).
   */
  virtual void __attribute__((noinline)) init() {}

  /**
   * @brief Abstract method to be implemented by derived classes to generate a
   * frame of graphics data.
   */
  virtual void draw_frame() = 0;

  /**
   * @brief POV-display strobe control: whether each LED column is blanked to
   *        black immediately after it is shown. Governs inter-column strip
   *        behavior, not framebuffer contents (that is `persist_pixels`).
   * @return true to strobe (each column a sharp slice with dark gaps); false to
   *         persist the lit column until the next overwrites it.
   */
  [[nodiscard]] bool strobe_columns() const { return strobe; }

  /**
   * @brief Whether this effect must render the FULL canvas per simulator worker
   *        rather than be clipped to a segment band.
   * @return True for an effect whose output can cross segment boundaries (for
   *         example MeshFeedback's unbounded warp or World::Trails);
   *         false (the default) when each segment can render independently.
   * @details Read by the segment drivers to leave the clip at full canvas for
   *          stateful effects. Set once at construction from the filter
   *          pipeline's `any_crosses_segments` fold; defaults to false.
   */
  [[nodiscard]] bool needs_full_frame() const { return full_frame; }

  /**
   * @brief Whether this effect copies its previous frame forward (trails/decay).
   * @return True when persist_pixels is set.
   * @details The segmented device driver leaves a persisting effect at full
   *          canvas, else its per-frame arm-half alternation breaks trail
   *          continuity.
   */
  [[nodiscard]] bool persists_pixels() const { return persist_pixels; }

  /**
   * @brief The effect's current clip region (display band + render margin).
   * @return Const reference to the clip region; mutate only through set_clip(),
   *         set_clip_x(), and set_margin(), which enforce its invariants.
   */
  [[nodiscard]] const ClipRegion &clip() const { return clip_region; }

  /**
   * @brief Driver sets display bounds (which segment this Teensy owns).
   * @param y0 Inclusive start row of the owned segment.
   * @param y1 Exclusive end row of the owned segment.
   * @param x0 Inclusive start column of the owned segment.
   * @param x1 Exclusive end column of the owned segment.
   */
  void set_clip(int y0, int y1, int x0, int x1) {
    check_clip_mutable();
    HS_CHECK(y0 >= 0 && y0 <= y1 && y1 <= clip_region.h && x0 >= 0 &&
                 x0 <= x1 && x1 <= clip_region.w,
             "set_clip band must be non-inverted and within canvas bounds");
    clip_region.y_start = y0;
    clip_region.y_end = y1;
    clip_region.x_start = x0;
    clip_region.x_end = x1;
    update_render_clip();
  }
  /**
   * @brief Update only the horizontal clip band, leaving the y bounds intact.
   *
   * For callers that retune just the x-band (canvas unit tests); no device
   * driver calls this.
   * @param x0 Inclusive start column of the horizontal clip band.
   * @param x1 Exclusive end column of the horizontal clip band.
   */
  void set_clip_x(int x0, int x1) {
    check_clip_mutable();
    HS_CHECK(x0 >= 0 && x0 <= x1 && x1 <= clip_region.w,
             "set_clip_x band must be non-inverted and within canvas width");
    clip_region.x_start = x0;
    clip_region.x_end = x1;
    update_render_clip();
  }
  /**
   * @brief Effect sets render margin for stateful filters.
   * @param m Render margin width in pixels.
   * @details ClipRegion's cylindrical wrap only corrects a single period of
   *          underflow, so its [0, w) contract holds only while margin < w. Trap
   *          a wrapping margin here rather than let a negative column leak into
   *          the per-fragment clip predicates.
   */
  void set_margin(int m) {
    check_clip_mutable();
    HS_CHECK(m >= 0 && m < clip_region.w,
             "render margin must be in [0, canvas width)");
    clip_region.margin = m;
    update_render_clip();
  }

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The Pixel color reference.
   * @note An override that applies a per-pixel transform (rather than a plain
   *       buffer read) must also override overrides_get_pixel() to return true;
   *       otherwise ISR/readback fast paths that index display_buffer() directly
   *       bypass the transform.
   */
  virtual const Pixel &get_pixel(int x, int y) const {
    // Debug-only bounds guard, matching the write-path accessors (stripped on
    // device, catches an out-of-range display read in test/sim).
    assert(x >= 0 && x < frame_width && y >= 0 && y < frame_height);
    return bufs[prev.load(std::memory_order_relaxed)][y * frame_width + x];
  }

  /**
   * @brief Raw base pointer to the currently displayed buffer (row-major,
   *        `width()` stride): `display_buffer()[y * width() + x]` equals
   *        `get_pixel(x, y)` for any effect that does not override `get_pixel`.
   *
   * ISR fast path: index pixels directly, skipping the virtual `get_pixel`
   * dispatch. Valid only until the next `advance_display()` flip. Effects that
   * override `get_pixel` must NOT use this (it bypasses their transform); check
   * `overrides_get_pixel()` first.
   */
  [[nodiscard]] const Pixel *display_buffer() const {
    return bufs[prev.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Whether this effect overrides get_pixel with a per-pixel transform
   *        that display_buffer() does NOT reflect.
   *
   * Base effects return false: `display_buffer()[y * width() + x]` equals
   * `get_pixel(x, y)`, so ISR fast paths may index the buffer directly. An
   * override that reads through a per-pixel transform (e.g. the RingTwist
   * scroller) returns true, forcing those paths back to virtual dispatch.
   */
  [[nodiscard]] virtual bool overrides_get_pixel() const { return false; }

  /** @brief Sets the complete-output transition envelope in [0,1]. */
  void set_output_envelope(float value) {
    HS_AUDIT_CHECK(std::isfinite(value) && value >= 0.0f && value <= 1.0f,
                   "output envelope must be finite and in [0,1]");
    output_envelope = static_cast<uint16_t>(value * 65535.0f + 0.5f);
  }

  [[nodiscard]] uint16_t output_envelope_u16() const { return output_envelope; }

  [[nodiscard]] Pixel apply_output_envelope(const Pixel &pixel) const {
    auto scale = [this](uint16_t channel) {
      return static_cast<uint16_t>(
          (static_cast<uint32_t>(channel) * output_envelope + 32767u) / 65535u);
    };
    return {scale(pixel.r), scale(pixel.g), scale(pixel.b)};
  }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  [[nodiscard]] inline int width() const { return frame_width; }
  /**
   * @brief Gets the height of the effect.
   * @return The height.
   */
  [[nodiscard]] inline int height() const { return frame_height; }
  /**
   * @brief Checks whether the queued frame has been picked up for display.
   * @return True when `prev == next` (no frame still waiting to be shown), so
   *         the writer is free to claim the other buffer.
   * @details The acquire load pairs with `advance_display()`'s release store:
   * everything the display ISR stored before the flip — the segmented driver's
   * display-window half — is visible to the writer this gate releases.
   */
  [[nodiscard]] inline bool buffer_free() const {
    return prev.load(std::memory_order_acquire) ==
           next.load(std::memory_order_relaxed);
  }
  /**
   * @brief Installs a callback run once per frame from the Canvas constructor.
   * @param hook Callback to install, or null to disable it.
   * @details Runs after the buffer_free() wait returns and before
   * `advance_buffer()` and the stale-pixel clear. That slot is the only one that
   * sees the segment window the ISR settled on during the wait *and* still
   * precedes the clear, so it is where the segmented driver sets the clip that
   * governs which band gets cleared.
   */
  void set_buffer_ready_hook(BufferReadyHook hook) { buffer_ready_hook = hook; }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   * @details The acquire load pairs with `queue_frame()`'s release store, so
   * the queued frame's pixel writes are visible before any read through `prev`.
   * The release store publishes the display flip and the caller's preceding
   * display-window update to `buffer_free()`.
   */
  inline void advance_display() {
    int n = next.load(std::memory_order_acquire);
    prev.store(n, std::memory_order_release);
  }

protected:
  /** @brief Whole-output transition envelope; 65535 is unattenuated. */
  uint16_t output_envelope = 65535u;
  /**
   * @brief Flag indicating if the previous frame's pixels should be copied to
   * the new buffer (for trails/decay).
   */
  bool persist_pixels;
  /**
   * @brief Full-canvas render gate (see needs_full_frame()); set once at
   * construction from the filter pipeline `any_crosses_segments` trait.
   */
  bool full_frame;
  /** @brief Whether frame generation samples pixels outside the display band. */
  bool reads_outside_band;
  /**
   * @brief POV column-strobe flag (see strobe_columns()); set at construction.
   */
  bool strobe;

private:
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content
   * to the new buffer.
   */
  inline void advance_buffer() {
    int c = cur.load(std::memory_order_relaxed) ? 0 : 1;
    // The new write buffer must not be the one the ISR is scanning out (prev);
    // with two buffers this holds only if buffer_free() gated the advance. Trap
    // it here (once per frame, cold) instead of tearing.
    HS_CHECK(c != prev.load(std::memory_order_relaxed),
             "advance_buffer: new write buffer is the one the display ISR is "
             "scanning out");
    cur.store(c, std::memory_order_relaxed);
    if (persist_pixels) {
      // The trail base is the last COMPLETED frame (next). The buffer_free()
      // gate forces prev == next, so copy from next and assert the equality
      // rather than depend on the gate silently across methods.
      int last = next.load(std::memory_order_relaxed);
      HS_CHECK(last == prev.load(std::memory_order_relaxed),
               "advance_buffer: trail base is not the last completed frame");
      memcpy(bufs[c], bufs[last], sizeof(Pixel) * frame_width * frame_height);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   * @details Publishes `cur` as the new `next`. The release store orders the
   * frame's pixel writes before the publish, pairing with the acquire load in
   * `advance_display()`; the IRQ-off bracket keeps the publish atomic against
   * the on-device display ISR.
   */
  inline void queue_frame() {
    hs::disable_interrupts();
    next.store(cur.load(std::memory_order_relaxed), std::memory_order_release);
    hs::enable_interrupts();
  }

  /**
   * @brief Points bufs at the shared static storage and zeroes both buffers.
   * @details noinline so the two full-frame fills are emitted once, not inlined
   * into both GCC constructor variants (C1/C2). Invoked from the ctor (not init())
   * because derived init() overrides do not chain to Effect::init().
   */
  HS_FLASH_MEMBER void __attribute__((noinline)) clear_buffers() {
    bufs[0] = buffer_a;
    std::fill_n(bufs[0], MAX_W * MAX_H, Pixel(0, 0, 0));
    bufs[1] = buffer_b;
    std::fill_n(bufs[1], MAX_W * MAX_H, Pixel(0, 0, 0));
  }

  inline void notify_buffer_ready() {
    if (buffer_ready_hook)
      buffer_ready_hook(*this);
  }

  std::atomic<int> prev{0}; /**< Buffer the ISR is currently reading. */
  std::atomic<int> cur{0};  /**< Buffer the main loop is currently writing. */
  std::atomic<int> next{0}; /**< Last completed frame, queued for display. */
  int frame_width;          /**< The width of the effect. */
  int frame_height;         /**< The height of the effect. */
  ClipRegion clip_region; /**< Segment clip region (display + render margin). */
  ClipRegion::XClip render_x_clip{};
  int render_y_start = 0;
  int render_y_end = 0;
  bool canvas_active = false; /**< A Canvas is currently drawing a frame. */

  void update_render_clip() {
    render_x_clip = clip_region.x_clip();
    render_y_start = clip_region.render_y_start();
    render_y_end = clip_region.render_y_end();
  }

  void check_clip_mutable() const {
    HS_CHECK(!canvas_active, "clip cannot change while a frame is active");
  }
  // Shared static storage for the double buffer. PRECONDITION: at most one Effect
  // live at a time (s_alive guard); a second would alias these arrays and the
  // prev/cur/next indices.
  static DMAMEM Pixel
      buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A (shared). */
  static DMAMEM Pixel
      buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B (shared). */
  Pixel *bufs[2]; /**< Pointers to the two buffer storage locations. */
  // True while an Effect is constructed-but-not-destroyed. Guards the
  // single-live-Effect precondition on the shared buffer_a/buffer_b: the ctor
  // traps if already set, the dtor clears it.
  static bool s_alive;
  BufferReadyHook buffer_ready_hook = nullptr;
};

/**
 * @brief Context class providing a safe, scoped interface to the current
 * drawing buffer.
 * @details Ensures the buffer advances correctly when the object is constructed
 * and queued when destroyed.
 */
class Canvas {
public:
  /**
   * @brief Constructs the Canvas, advancing the effect buffer and clearing
   * whatever of it can still hold stale pixels.
   * @param owner The effect instance owning the buffer.
   */
  explicit Canvas(Effect &owner) : effect(owner) {
    HS_CHECK(!effect.canvas_active, "Canvas already active for this Effect");
    wait_for_free_buffer();
    // Ordering is load-bearing: the hook runs post-wait but pre-clear, so the
    // clip it sets is the one clear_stale_pixels() honours.
    effect.notify_buffer_ready();
    effect.advance_buffer();
    if (!effect.persist_pixels) {
      HS_PROFILE(canvas_clear);
      clear_stale_pixels();
    }
    effect.canvas_active = true;
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() {
    effect.queue_frame();
    effect.canvas_active = false;
  }

  Canvas(const Canvas &) = delete;
  Canvas(Canvas &&) = delete;
  Canvas &operator=(const Canvas &) = delete;
  Canvas &operator=(Canvas &&) = delete;

  /**
   * @brief Accesses a pixel in the current drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int x, int y) {
    // assert, NOT HS_CHECK — the one hot-loop exception: a debug-only bounds
    // guard (stripped on device), since an always-on branch on every pixel access
    // is the single place HS_CHECK's contract forbids it.
    assert(x >= 0 && x < effect.frame_width && y >= 0 &&
           y < effect.frame_height);
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)]
                      [y * effect.frame_width + x];
  }

  /**
   * @brief Accesses a pixel in the previous drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Copy of the Pixel from the previous frame.
   * @details Returns by value, unlike `operator()`: the previous frame is
   *          read-only, so a reference into it would invite an accidental write
   *          to a buffer about to be recycled. A Pixel is small enough the copy
   *          is free.
   */
  inline Pixel prev(int x, int y) const {
    assert(x >= 0 && x < effect.frame_width && y >= 0 &&
           y < effect.frame_height);
    return effect.bufs[effect.prev.load(std::memory_order_relaxed)]
                      [y * effect.frame_width + x];
  }

  /**
   * @brief Raw base pointer to the previous-frame buffer (row-major,
   *        `width()` stride): `prev_data()[y * width() + x]` equals
   *        `prev(x, y)`.
   * @return Const base pointer, valid until this Canvas is destroyed.
   * @details Const because `prev()` returns by value specifically so a
   *          reference cannot be used to write into a buffer about to be
   *          recycled; the const pointer preserves that guarantee.
   */
  [[nodiscard]] inline const Pixel *prev_data() const {
    return effect.bufs[effect.prev.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Raw base pointer to the current drawing buffer (row-major,
   *        `width()` stride): `data()[y * width() + x]` equals
   *        `(*this)(x, y)`.
   * @return Mutable base pointer, valid until this Canvas is destroyed.
   */
  [[nodiscard]] inline Pixel *data() {
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int xy) {
    assert(xy >= 0 && xy < effect.frame_width * effect.frame_height);
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)][xy];
  }

  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    int c = effect.cur.load(std::memory_order_relaxed);
    std::fill_n(effect.bufs[c], effect.frame_width * effect.frame_height,
                Pixel(0, 0, 0));
  }

  /**
   * @brief Gets the width of the underlying effect.
   * @return The width in pixels.
   */
  [[nodiscard]] inline int width() const { return effect.width(); }
  /**
   * @brief Gets the height of the underlying effect.
   * @return The height in pixels.
   */
  [[nodiscard]] inline int height() const { return effect.height(); }
  /**
   * @brief Gets the effect's current clip region.
   * @return Const reference to the clip region (display + render margin).
   */
  [[nodiscard]] inline const ClipRegion &clip() const {
    return effect.clip_region;
  }

  /** @brief Tests a column against the frame's cached render clip. */
  [[nodiscard]] inline bool clip_contains_x(int x) const {
    return !effect.render_x_clip.clipped(x);
  }

  /** @brief Tests a row against the frame's cached render clip. */
  [[nodiscard]] inline bool clip_contains_y(int y) const {
    return y >= effect.render_y_start && y < effect.render_y_end;
  }

  /**
   * @brief Checks if debug visuals are enabled.
   * @return True if debugging is active.
   * @details Read by Scan::rasterize, Scan::rasterize_solid and
   * Scan::RingGroup. The fused walks — Scan::rasterize_face (and so
   * Scan::Mesh) and Scan::DistortedRingStack — ignore it and render
   * identically either way; a per-shape fallback costs them the whole shared
   * rasterizer in ITCM.
   */
  inline bool debug() const { return effect.debug_visuals; }

#if HS_ENABLE_TEST_HOOKS
  /**
   * @brief Test-only count of buffer_free() spin iterations across all Canvas
   *        ctors in this process.
   * @return The running spin-iteration count.
   * @details Lets a test detect that a ctor has actually entered the wait loop
   *          and release it on observed progress, rather than racing a fixed
   *          sleep. Compiled out of the device/sim image.
   */
  static unsigned long buffer_free_spin_count() {
    return s_buffer_free_spins.load(std::memory_order_relaxed);
  }
#endif

private:
  /**
   * @brief Spins until the effect has a free back buffer.
   * @details Not a TOCTOU race: single-core, strict index ownership — the main
   * loop writes cur/next, the ISR only sets prev = next, so the gate can't be
   * falsified between check and flip. buffer_free() is re-satisfied only when the
   * display ISR advances at a frame boundary, so an unbounded wait means that ISR
   * stalled and the watchdog traps rather than hanging. The outer buffer_free()
   * guard skips the micros() reads on the no-wait path. The counter name is
   * load-bearing: Profile.ino derives each frame's render as wall minus the
   * counter whose name ends in _buffer_wait.
   */
  void wait_for_free_buffer() {
    HS_PROFILE(canvas_buffer_wait);
    if (effect.buffer_free())
      return;
#if HS_ENABLE_TEST_HOOKS
    const unsigned long WATCHDOG_US = buffer_free_watchdog_us();
#else
    constexpr unsigned long WATCHDOG_US = BUFFER_FREE_WATCHDOG_US;
#endif
    const unsigned long wait_start = micros();
    while (!effect.buffer_free()) {
      HS_CHECK(micros() - wait_start < WATCHDOG_US,
               "buffer_free watchdog timeout — display ISR stalled");
#if HS_ENABLE_TEST_HOOKS
      s_buffer_free_spins.fetch_add(1, std::memory_order_relaxed);
#endif
    }
  }

  /**
   * @brief Clears whatever of the freshly acquired buffer can still show stale
   *        pixels from the frame that last wrote it.
   * @details For an effect that does not read outside its display clip, the
   *          leftovers there are invisible and only the display band has to
   *          be cleared. It does draw into the margin-expanded render bounds,
   *          but that band is write-only scratch: nothing samples or displays
   *          it. A filter that samples across the band edge needs the whole
   *          buffer cleared, independently of whether its output crosses
   *          segment boundaries. Unsegmented targets clip to the full canvas,
   *          where the two are the same fill.
   */
  void clear_stale_pixels() {
#if HS_ENABLE_TEST_HOOKS
    if (effect.force_full_buffer_clear) {
      clear_buffer();
      return;
    }
#endif
    if (effect.reads_outside_band)
      clear_buffer();
    else
      clear_display_clip_buffer();
  }

  /**
   * @brief Clears the current display clip, excluding its render margin.
   * @details Device builds keep this helper out of line and outside ITCM.
   */
  HS_COLD_MEMBER void clear_display_clip_buffer() {
    const int c = effect.cur.load(std::memory_order_relaxed);
    const ClipRegion &clip = effect.clip_region;
    const int span = clip.x_end - clip.x_start;
    Pixel *const buffer = effect.bufs[c];

    if (span == effect.frame_width) {
      std::fill_n(buffer + clip.y_start * effect.frame_width,
                  span * (clip.y_end - clip.y_start), Pixel(0, 0, 0));
      return;
    }

    for (int y = clip.y_start; y < clip.y_end; ++y) {
      std::fill_n(buffer + y * effect.frame_width + clip.x_start, span,
                  Pixel(0, 0, 0));
    }
  }

  /** Watchdog bound for the ctor buffer_free() spin (µs). One display
   *  revolution is tens-to-hundreds of ms even at low RPM; 2 s is well above
   *  that, so only a genuinely stalled display ISR trips it. */
  static constexpr unsigned long BUFFER_FREE_WATCHDOG_US = 2000000UL;

#if HS_ENABLE_TEST_HOOKS
  /**
   * @brief Test-only watchdog bound, overridable via HS_BUFFER_FREE_WATCHDOG_US.
   * @return The bound in µs; BUFFER_FREE_WATCHDOG_US unless the environment
   *         supplies a positive override.
   * @details The watchdog is a trap, so tripping it kills the whole test
   *          process rather than failing one test; a sanitizer job whose
   *          threads deschedule past the shipping bound raises it here. Read
   *          once per process, outside the spin loop.
   */
  static unsigned long buffer_free_watchdog_us() {
    static const unsigned long RESOLVED = [] {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      const char *e = std::getenv("HS_BUFFER_FREE_WATCHDOG_US");
#pragma clang diagnostic pop
      if (e) {
        const unsigned long v = std::strtoul(e, nullptr, 10);
        if (v > 0)
          return v;
      }
      return BUFFER_FREE_WATCHDOG_US;
    }();
    return RESOLVED;
  }
#endif

  Effect &effect; /**< Reference to the owning Effect instance. */
#if HS_ENABLE_TEST_HOOKS
  inline static std::atomic<unsigned long> s_buffer_free_spins{0};
#endif
};
