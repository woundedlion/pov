/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstring>
#include <cstdio>
#include <algorithm>
#include "platform.h"
#include "constants.h"
#include "color.h"

class Canvas;

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for drawing a frame.
 */
class Effect {
  friend class Canvas;

public:
  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect.
   * @param H The height (resolution) of the effect.
   */
  Effect(int W, int H) :
    persist_pixels(true),
    width_(W),
    height_(H)
  {
    bufs_[0] = buffer_a;
    memset(bufs_[0], 0, sizeof(Pixel) * MAX_W * MAX_H);
    bufs_[1] = buffer_b;
    memset(bufs_[1], 0, sizeof(Pixel) * MAX_W * MAX_H);
  }

  virtual ~Effect() {
  };

  /**
   * @brief Abstract method to be implemented by derived classes to generate a frame of graphics data.
   */
  virtual void draw_frame() = 0;
  /**
   * @brief Indicates whether the display should show black between POV sweeps.
   * @return True if the background should be shown (often black).
   */
  virtual bool show_bg() const = 0;

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The Pixel color reference.
   */
  virtual const Pixel& get_pixel(int x, int y) const {
    return bufs_[prev_][x * height_ + y];
  }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  inline int width() const { return width_; }
  /**
   * @brief Gets the height of the effect.
   * @return The height.
   */
  inline int height() const { return height_; }
  /**
   * @brief Checks if the display buffer and drawing buffer are synced.
   * @return True if the current frame is ready to be shown.
   */
  inline bool buffer_free() const { return prev_ == next_; }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   */
  inline void advance_display() { prev_ = next_; }
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content to the new buffer.
   */
  inline void advance_buffer() {
    hs::disable_interrupts();
    cur_ = cur_ ? 0 : 1;
    hs::enable_interrupts();
    if (persist_pixels) {
      memcpy(bufs_[cur_], bufs_[prev_], sizeof(Pixel) * width_ * height_);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   */
  inline void queue_frame() {
    hs::disable_interrupts();
    next_ = cur_;
    hs::enable_interrupts();
  }

protected:
  /**
   * @brief Flag indicating if the previous frame's pixels should be copied to the new buffer (for trails/decay).
   */
  bool persist_pixels;

private:
  volatile int prev_ = 0, cur_ = 0, next_ = 0; /**< Pointers/indices for double buffering logic. */
  int width_; /**< The width of the effect. */
  int height_; /**< The height of the effect. */
  inline static Pixel buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A. */
  inline static Pixel buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B. */
  Pixel* bufs_[2]; /**< Pointers to the two buffer storage locations. */
};

/**
 * @brief Context class providing a safe, scoped interface to the current drawing buffer.
 * @details Ensures the buffer advances correctly when the object is constructed and queued when destroyed.
 */
class Canvas {
public:
  /**
   * @brief Constructs the Canvas, advancing the effect buffer and optionally clearing it.
   * @param effect The effect instance owning the buffer.
   */
  Canvas(Effect& effect) : effect_(effect) {
    while (!effect_.buffer_free()) {}
    //    start_time = millis();
    effect_.advance_buffer();
    if (!effect_.persist_pixels) {
      clear_buffer();
    }
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() {
    //    Serial.printf("draw_frame_duration: %d\n", (millis() - start_time));
    effect_.queue_frame();
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Reference to the Pixel.
   */
  inline Pixel& operator()(int x, int y) {
    return effect_.bufs_[effect_.cur_][x * effect_.height_ + y];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the Pixel.
   */
  inline Pixel& operator()(int xy) {
    return effect_.bufs_[effect_.cur_][xy];
  }

  /**
   * @brief Gets the width of the canvas.
   * @return The width.
   */

  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    memset(effect_.bufs_[effect_.cur_], 0, sizeof(Pixel) * effect_.width_ * effect_.height_);
  }

  inline int width() const { return effect_.width(); }
  inline int height() const { return effect_.height(); }

private:
  Effect& effect_; /**< Reference to the owning Effect instance. */
  size_t start_time; /**< Tracks frame drawing duration (debug). */
};
