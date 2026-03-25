/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_CANVAS_H_
#define HOLOSPHERE_CORE_CANVAS_H_

#include <cstring>
#include <cstdio>
#include <algorithm>
#include <variant>
#include <atomic>
#include "platform.h"
#include "constants.h"
#include "color.h"
#include <array>

class Canvas;

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for
 * drawing a frame.
 */
class Effect {
  friend class Canvas;

public:
  bool debug_visuals = false; /**< Flag to enable visual debugging overlays. */

  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect.
   * @param H The height (resolution) of the effect.
   */
  Effect(int W, int H) : persist_pixels(false), width_(W), height_(H) {
    bufs_[0] = buffer_a;
    std::fill_n(bufs_[0], MAX_W * MAX_H, Pixel(0, 0, 0));
    bufs_[1] = buffer_b;
    std::fill_n(bufs_[1], MAX_W * MAX_H, Pixel(0, 0, 0));
  }

  virtual ~Effect() {};

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
  virtual const Pixel &get_pixel(int x, int y) const {
    return bufs_[prev_.load(std::memory_order_relaxed)][y * width_ + x];
  }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  [[nodiscard]] inline int width() const { return width_; }
  /**
   * @brief Gets the height of the effect.
   * @return The height.
   */
  [[nodiscard]] inline int height() const { return height_; }
  /**
   * @brief Checks if the display buffer and drawing buffer are synced.
   * @return True if the current frame is ready to be shown.
   */
  [[nodiscard]] inline bool buffer_free() const {
    return prev_.load(std::memory_order_relaxed) ==
           next_.load(std::memory_order_relaxed);
  }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   */
  inline void advance_display() {
    prev_.store(next_.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
  }
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content
   * to the new buffer.
   */
  inline void advance_buffer() {
    int c = cur_.load(std::memory_order_relaxed) ? 0 : 1;
    cur_.store(c, std::memory_order_relaxed);
    if (persist_pixels) {
      memcpy(bufs_[c], bufs_[prev_.load(std::memory_order_relaxed)],
             sizeof(Pixel) * width_ * height_);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   */
  inline void queue_frame() {
    hs::disable_interrupts();
    next_.store(cur_.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
    hs::enable_interrupts();
  }

  /**
   * @brief Defines a runtime-adjustable parameter.
   */
  struct ParamDef {
    const char *name; /**< Parameter name. */
    std::variant<float *, bool *>
        target;             /**< Type-safe pointer to the variable. */
    float min = 0;          /**< Minimum value (for floats). */
    float max = 1;          /**< Maximum value (for floats). */
    float defaultValue = 0; /**< Default value. */

    /** @brief Read the current value as float (bool maps to 0/1). */
    float get() const {
      return std::visit(
          [](auto *p) -> float {
            using T = std::remove_pointer_t<decltype(p)>;
            if constexpr (std::is_same_v<T, bool>)
              return *p ? 1.0f : 0.0f;
            else
              return *p;
          },
          target);
    }

    /** @brief Write a float value (bool threshold at 0.5). */
    void set(float v) {
      std::visit(
          [v](auto *p) {
            using T = std::remove_pointer_t<decltype(p)>;
            if constexpr (std::is_same_v<T, bool>)
              *p = (v > 0.5f);
            else
              *p = v;
          },
          target);
    }

    /** @brief Check if this parameter targets a bool. */
    bool is_bool() const { return std::holds_alternative<bool *>(target); }
  };

  struct ParamList {
    std::array<ParamDef, 32> elements;
    size_t count = 0;

    const ParamDef *begin() const { return elements.data(); }
    const ParamDef *end() const { return elements.data() + count; }
    ParamDef *begin() { return elements.data(); }
    ParamDef *end() { return elements.data() + count; }

    ParamDef *find(const char *name) {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
    const ParamDef *find(const char *name) const {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
    size_t size() const { return count; }
  };

  // Parameter System
  /**
   * @brief Updates a parameter's value by name.
   * @param name The name of the parameter.
   * @param value The new value (mapped to bool if necessary).
   */
  void updateParameter(const char *name, float value) {
    auto *def = parameters.find(name);
    if (def != nullptr) {
      def->set(value);
    }
  }

  /**
   * @brief Retrieves the list of registered parameters.
   * @return Const reference to the parameter list.
   */
  const ParamList &getParameters() const { return parameters; }

protected:
  /**
   * @brief Flag indicating if the previous frame's pixels should be copied to
   * the new buffer (for trails/decay).
   */
  bool persist_pixels;
  ParamList parameters; /**< List of parameters. */

  /**
   * @brief Registers a floating-point parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the float variable.
   * @param min Minimum value.
   * @param max Maximum value.
   */
  void registerParam(const char *name, float *ptr, float min = 0.0f,
                     float max = 1.0f) {
    if (parameters.count < parameters.elements.size()) {
      parameters.elements[parameters.count++] = {name, ptr, min, max, *ptr};
    }
  }

  /**
   * @brief Registers a boolean parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable.
   * @param defaultValue Initial value.
   */
  void registerParam(const char *name, bool *ptr, bool defaultValue = false) {
    *ptr = defaultValue;
    if (parameters.count < parameters.elements.size()) {
      parameters.elements[parameters.count++] = {name, ptr, 0.0f, 1.0f,
                                                 (float)defaultValue};
    }
  }

private:
  std::atomic<int> prev_{0}; /**< Buffer the ISR is currently reading. */
  std::atomic<int> cur_{0};  /**< Buffer the main loop is currently writing. */
  std::atomic<int> next_{0}; /**< Last completed frame, queued for display. */
  int width_;             /**< The width of the effect. */
  int height_;            /**< The height of the effect. */
  static DMAMEM Pixel
      buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A. */
  static DMAMEM Pixel
      buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B. */
  Pixel *bufs_[2]; /**< Pointers to the two buffer storage locations. */
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
   * @brief Constructs the Canvas, advancing the effect buffer and optionally
   * clearing it.
   * @param effect The effect instance owning the buffer.
   */
  Canvas(Effect &effect) : effect_(effect) {
    while (!effect_.buffer_free()) {
    }
    start_time = micros();
    effect_.advance_buffer();
    if (!effect_.persist_pixels) {
      clear_buffer();
    }
    clear_time = micros() - start_time;
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() {
#ifdef ARDUINO
    unsigned long render_time = micros() - start_time - clear_time;
    Serial.print("clear ");
    Serial.print(clear_time);
    Serial.print(" render ");
    Serial.println(render_time);
#endif
    effect_.queue_frame();
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int x, int y) {
    assert(x >= 0 && x < effect_.width_ && y >= 0 && y < effect_.height_);
    return effect_.bufs_[effect_.cur_.load(std::memory_order_relaxed)][y * effect_.width_ + x];
  }

  /**
   * @brief Accesses a pixel in the previous drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Copy of the Pixel from the previous frame.
   */
  inline Pixel prev(int x, int y) const {
    assert(x >= 0 && x < effect_.width_ && y >= 0 && y < effect_.height_);
    return effect_.bufs_[effect_.prev_.load(std::memory_order_relaxed)][y * effect_.width_ + x];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int xy) {
    assert(xy >= 0 && xy < effect_.width_ * effect_.height_);
    return effect_.bufs_[effect_.cur_.load(std::memory_order_relaxed)][xy];
  }


  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    int c = effect_.cur_.load(std::memory_order_relaxed);
    std::fill_n(effect_.bufs_[c], effect_.width_ * effect_.height_,
                Pixel(0, 0, 0));
  }

  [[nodiscard]] inline int width() const { return effect_.width(); }
  [[nodiscard]] inline int height() const { return effect_.height(); }
  /**
   * @brief Checks if debug visuals are enabled.
   * @return True if debugging is active.
   */
  inline bool debug() const { return effect_.debug_visuals; }

private:
  Effect &effect_;   /**< Reference to the owning Effect instance. */
  unsigned long start_time; /**< Tracks frame drawing duration (debug). */
  unsigned long clear_time = 0; /**< Tracks clear phase duration. */
};

#endif // HOLOSPHERE_CORE_CANVAS_H_
