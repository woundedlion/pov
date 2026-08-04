/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 *
 * @file effects_legacy.h
 * @brief Legacy effects written before this engine, plus the timers,
 *        oscillators and color sequences they share.
 */
#include "render/canvas.h"
#include "math/rotate.h"
#include "engine/util.h"
#include "render/led.h"

///////////////////////////////////////////////////////////////////////////////
// Legacy effects built before this engine
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Ages a pixel one step along a fading rainbow trail.
 * @param c Pixel to age, modified in place.
 * @param hue_falloff Hue advance per step, in 8-bit hue units.
 * @param dim_falloff Value (brightness) reduction per step, saturating.
 * @details Forces full saturation so already-trailing pixels keep vivid color
 * as they decay.
 */
void trail_rainbow(Pixel &c, uint8_t hue_falloff = 32,
                   uint8_t dim_falloff = 32) {
  auto p = rgb2hsv_approximate(static_cast<CRGB>(c));
  p.s = 255;
  p.h += hue_falloff;
  p.v = qsub8(p.v, dim_falloff);
  c = p;
}

/**
 * @brief Ages a pixel along a rainbow trail with slower decay for dim pixels.
 * @param c Pixel to age, modified in place.
 * @param hue_falloff Hue advance per step, in 8-bit hue units.
 * @param dim_falloff Value (brightness) reduction per step, saturating.
 * @details Halves the hue/dim step for dim pixels (value < 128), stretching the
 * tail so faint trail ends fade and rotate more slowly.
 */
void trail_rainbow_lin(Pixel &c, uint8_t hue_falloff = 32,
                       uint8_t dim_falloff = 32) {
  auto p = rgb2hsv_approximate(static_cast<CRGB>(c));
  p.s = 255;
  if (p.v < 128) {
    p.h += hue_falloff / 2;
    p.v = qsub8(p.v, dim_falloff / 2);
  } else {
    p.h += hue_falloff;
    p.v = qsub8(p.v, dim_falloff);
  }
  c = p;
}

/**
 * @brief Randomly kills a pixel to produce a sparse dissolve.
 * @param p Color to maybe extinguish, modified in place.
 * @param prob Kill probability numerator over 255 (chance the value is zeroed).
 * @details With probability prob/255, zeroes the pixel's value while keeping its
 * hue and saturation.
 */
void disintegrate(CHSV &p, int prob) {
  if (random8() <= prob) {
    p = CHSV(p.h, p.s, 0);
  }
}

/**
 * @brief Returns a saturating linear-blend functor.
 * @param a Blend weight in [0, 1] for the second color.
 * @return A functor (c1, c2) -> CRGB mixing c1 and c2 per channel with weights
 * (1-a) and a respectively, saturating each channel.
 */
auto blend(float a) {
  return [a](const Pixel &src, const CRGB &c2) {
    const CRGB c1 = static_cast<CRGB>(src);
    return CRGB(qadd8(c1.r * (1 - a), c2.r * a),
                qadd8(c1.g * (1 - a), c2.g * a),
                qadd8(c1.b * (1 - a), c2.b * a));
  };
}

/**
 * @brief Plots a color at a fractional position with anti-aliasing.
 * @param cv Canvas to draw into.
 * @param x Fractional column coordinate; wraps around the ring (cylinder).
 * @param y Fractional row coordinate; clamped so the top row never bleeds past
 * the last row.
 * @param c Color to deposit.
 * @details Distributes c over the four surrounding pixels weighted by sub-pixel
 * coverage.
 */
void plot_aa(Canvas &cv, const float &x, const float &y, const CHSV &c) {
  int x_i = floorf(x);
  int y_i = floorf(y);
  float x_m = x - x_i;
  float y_m = y - y_i;

  auto p = CRGB(c);

  float v = (1 - x_m) * (1 - y_m);
  cv(x_i, y_i) = blend(v)(cv(x_i, y_i), p);

  v = x_m * (1 - y_m);
  cv((x_i + 1) % cv.width(), y_i) =
      blend(v)(cv((x_i + 1) % cv.width(), y_i), p);

  if (y_i < cv.height() - 1) {
    v = (1 - x_m) * y_m;
    cv(x_i, y_i + 1) = blend(v)(cv(x_i, y_i + 1), p);

    v = x_m * y_m;
    cv((x_i + 1) % cv.width(), y_i + 1) =
        blend(v)(cv((x_i + 1) % cv.width(), y_i + 1), p);
  }
}

/**
 * @brief Timer that fires at randomized intervals; poll elapsed() each frame.
 * @details Each fire rearms with a fresh interval drawn from [min_ms, max_ms).
 */
class PollingRandomTimer {
public:
  /**
   * @brief Constructs the timer and arms the first random interval.
   * @param min_ms Minimum interval between fires, in milliseconds.
   * @param max_ms Maximum interval between fires, in milliseconds.
   */
  PollingRandomTimer(uint32_t min_ms, uint32_t max_ms)
      : min_ms(min_ms), max_ms(max_ms) {
    randomSeed(analogRead(PIN_RANDOM));
    next = millis() + random(min_ms, max_ms);
  }

  /**
   * @brief Reports whether the current interval has elapsed.
   * @return True once per random interval; rearms with a fresh interval on each
   * fire.
   */
  bool elapsed() {
    if (millis() >= next) {
      next = millis() + random(min_ms, max_ms);
      return true;
    }
    return false;
  }

private:
  uint32_t min_ms; /**< Minimum interval between fires, in milliseconds. */
  uint32_t max_ms; /**< Maximum interval between fires, in milliseconds. */
  uint32_t next;   /**< Absolute time (ms) of the next scheduled fire. */
};

/**
 * @brief Random-walk color generator over a complementary palette.
 * @details Mostly returns the complement (hue + 64) of the previous color,
 * occasionally nudging the hue +/- a small step instead.
 */
class ComplementaryColorSequence {
public:
  /**
   * @brief Constructs the sequence with a random initial color.
   */
  ComplementaryColorSequence() : last(0, random(128, 255), 255) {}

  /**
   * @brief Produces the next color in the sequence.
   * @return The next CHSV color; updates internal state.
   */
  CHSV get() {
    CHSV r;
    if (random8() < 32) {
      last = CHSV(last.h + random(16, 32), random(128, 255), 255);
      return last;
    } else if (random8() < 64) {
      last = CHSV(last.h - random(16, 32), random(128, 255), 255);
      return last;
    }
    r = CHSV(last.h + 64, last.s, last.v);
    last = r;
    return r;
  }

private:
  CHSV last; /**< Color returned by the previous get() call. */
};

///////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Integer triangle-wave counter bouncing between [min, max].
 * @details Optionally pauses end_delay polls at each endpoint before reversing
 * direction.
 */
class Oscillator {
public:
  /**
   * @brief Constructs the oscillator starting at min, moving toward max.
   * @param min Lower bound of the bounce range (inclusive).
   * @param max Upper bound of the bounce range (inclusive).
   * @param dwell_polls Number of extra polls to dwell at each endpoint before
   * reversing.
   */
  Oscillator(int min, int max, int dwell_polls = 0)
      : t(min), lower_bound(min), upper_bound(max), end_delay(dwell_polls),
        end_count(0), direction(1) {}

  /**
   * @brief Returns the current value, then advances one step.
   * @return The value before this step (honoring endpoint delay).
   */
  int get() {
    int r = t;
    if ((t == upper_bound && direction == 1) ||
        (t == lower_bound && direction == -1)) {
      if (end_count++ >= end_delay) {
        direction = direction * -1;
        end_count = 0;
      } else {
        return r;
      }
    }
    t += direction;
    return r;
  }

  /**
   * @brief Returns the current travel direction.
   * @return +1 when moving toward max, -1 when moving toward min.
   */
  int dir() const { return direction; }

  /**
   * @brief Sets the endpoint dwell count.
   * @param d Number of extra polls to dwell at each endpoint before reversing.
   */
  void set_delay(int d) { end_delay = d; }

private:
  int t;           /**< Current counter value. */
  int lower_bound; /**< Lower bound of the bounce range (inclusive). */
  int upper_bound; /**< Upper bound of the bounce range (inclusive). */
  int end_delay;   /**< Polls to dwell at each endpoint before reversing. */
  int end_count;   /**< Polls already dwelled at the current endpoint. */
  int direction;   /**< Current travel direction, +1 or -1. */
};

/**
 * @brief Chain of beads (one dot per row) dragged across the ring.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The driven row's motion drags the others like a string of beads,
 * each lagging within a varying gap and leaving rainbow trails; the chain is
 * replicated around the ring for symmetry.
 */
template <int W, int H> class ChainWiggle : public Effect {
public:
  /**
   * @brief Constructs the effect and seeds one dot per row at column 0.
   */
  FLASHMEM ChainWiggle() : Effect(W, H, {.persist = true}), osc(1, 4) {
    for (size_t i = 0; i < (sizeof(dots) / sizeof(Dot)); ++i) {
      dots[i] = Dot(0, i, 0);
    }
  }

  /**
   * @brief Renders one frame of the chained-bead animation.
   * @details Ages all trails, then periodically flips/randomizes speed and gap
   * before pulling the chain one step.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 24, 24);
      }
    }

    EVERY_N_MILLIS(1000) { speed *= -1; }

    EVERY_N_MILLIS(4000) {
      gap = osc.get();
      speed = dir(speed) * random(1, 2);
    }

    pull(c, dots, 0, speed);
  };

private:
  /**
   * @brief One chain bead: ring position x on row y, signed velocity, draw hue.
   */
  struct Dot {
    /**
     * @brief Constructs a zeroed dot at the origin with red hue.
     */
    Dot() : x(0), y(0), v(0), hue(HUE_RED) {}

    /**
     * @brief Constructs a dot at a given position and velocity.
     * @param x Ring column position.
     * @param y Row index.
     * @param v Signed velocity in columns per step.
     */
    Dot(int x, int y, int v) : x(x), y(y), v(v), hue(HUE_RED) {}

    /**
     * @brief Copy-constructs a dot.
     * @param d Dot to copy.
     */
    Dot(const Dot &d) : x(d.x), y(d.y), v(d.v), hue(d.hue) {}

    int x;       /**< Ring column position. */
    int y;       /**< Row index. */
    int v;       /**< Signed velocity in columns per step. */
    uint8_t hue; /**< Hue used to draw the bead. */
  };

  /**
   * @brief Draws a color at x and at n-1 evenly spaced copies around the ring.
   * @param c Canvas to draw into.
   * @param x Ring column position of the first copy.
   * @param y Row index.
   * @param color Color to draw.
   * @param n Total number of evenly spaced copies around the ring.
   */
  void replicate(Canvas &c, int x, int y, const CHSV &color, int n) {
    for (int i = 0; i < W; i += W / n) {
      c(wrap(x + i, W), y) = color;
    }
  }

  /**
   * @brief Signed shortest ring distance from a to b on a width-m circle.
   * @param a Start column.
   * @param b End column.
   * @param m Ring circumference in columns.
   * @return Number of columns from a to b; positive if b is reached going
   * forward (+x), negative going backward.
   */
  int distance(int a, int b, int m) {
    int fwd = wrap(b - a, m);
    int rev = wrap(a - b, m);
    if (fwd <= rev) {
      return fwd;
    }
    return -rev;
  }

  /**
   * @brief Returns the sign of a velocity.
   * @param v Velocity to test.
   * @return -1 when v is negative, +1 otherwise (treats 0 as +1).
   */
  int dir(int v) { return v < 0 ? -1 : 1; }

  /**
   * @brief Reverses the velocity of the dot on row y.
   * @param dots Array of per-row dots.
   * @param y Row index whose dot to reverse.
   */
  void reverse(Dot *dots, int y) { dots[y].v *= -1; }

  /**
   * @brief Drives row y at velocity v, then drags every other row toward it.
   * @param c Canvas to draw into.
   * @param dots Array of per-row dots.
   * @param y Driven row index.
   * @param v Signed velocity for the driven row, in columns per step.
   */
  void pull(Canvas &c, Dot *dots, int y, int v) {
    dots[y].v = v;
    move(c, dots[y]);
    for (int i = y - 1; i >= 0; --i) {
      drag(c, dots[i + 1], dots[i]);
    }
    for (int i = y + 1; i < H; ++i) {
      drag(c, dots[i - 1], dots[i]);
    }
  }

  /**
   * @brief Advances a follower dot toward its leader, keeping the chain taut.
   * @param c Canvas to draw into.
   * @param leader Dot being followed.
   * @param follower Dot to advance, modified in place.
   * @details If keeping its own velocity would let it drift more than gap
   * behind, snaps it onto the leader's path and locks its speed to the
   * leader's so the chain stays taut.
   */
  void drag(Canvas &c, Dot &leader, Dot &follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (abs(distance(dest, leader.x, W)) > gap) {
      // Sweep out to the far end of the leader's just-traveled path...
      dest = wrap(leader.x - dir(leader.v) * (abs(leader.v) - 1 + gap), W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      // ...then settle one gap-length behind the leader.
      dest = wrap(leader.x - dir(leader.v) * gap, W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      follower.v = leader.v;
    } else {
      move(c, follower);
    }
  }

  /**
   * @brief Steps a dot from its current x to x+v, drawing the swept path.
   * @param c Canvas to draw into.
   * @param dot Dot to move, modified in place.
   * @details Draws replicated copies along the whole swept path so no gaps
   * appear between frames.
   */
  void move(Canvas &c, Dot &dot) {
    int dest = wrap(dot.x + dot.v, W);
    for (int i = dot.x;; i = wrap(i + dir(dot.v), W)) {
      replicate(c, i, dot.y, CHSV(dot.hue, 255, 255), replicas);
      if (i == dest)
        break;
    }
    dot.x = dest;
  }

  Dot dots[H];    /**< Per-row chain beads. */
  Oscillator osc; /**< Drives the periodic gap value. */
  int gap = 1;    /**< Maximum lag (columns) a follower may trail its leader. */
  int speed = 1;  /**< Current signed drive velocity in columns per step. */
  int replicas = 2; /**< Number of symmetric copies drawn around the ring. */
};

/**
 * @brief Seed spokes that each row twists around the ring at its own offset.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details COUNT seed spokes are drawn; a leader row drives the twist and
 * drags the rest, with recursive trails filling the swept arcs. The effect
 * periodically halts, then re-randomizes direction and leader.
 */
template <int W, int H> class RingTwist : public Effect {
public:
  /**
   * @brief Constructs the effect and paints COUNT seed spokes.
   */
  FLASHMEM RingTwist()
      : Effect(W, H, {.persist = true}), seed(CHSV(104, 0, 255)) {
    randomSeed(analogRead(PIN_RANDOM));
    Canvas c(*this);
    for (int x = 0; x < W; x += W / COUNT) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = seed;
      }
    }
  }

  /**
   * @brief Samples a pixel with each row scrolled by its own offset.
   * @param x Logical column index.
   * @param y Row index.
   * @return The stored pixel from column wrap(x - pos[y], W) on row y.
   */
  const Pixel &get_pixel(int x, int y) const override {
    return Effect::get_pixel(wrap(x - pos[y], W), y);
  }

  /**
   * @brief Declares that this effect overrides per-pixel sampling.
   * @return True; the per-row scroll offset lives only in get_pixel(), so ISR
   * fast paths must fall back to virtual get_pixel dispatch.
   */
  bool overrides_get_pixel() const override { return true; }

  /**
   * @brief Renders one frame of the ring-twist animation.
   * @details Decays all non-seed pixels, drives the leader row, and
   * periodically latches a stop position so the twist eventually parks.
   */
  void draw_frame() {
    Canvas c(*this);
    canvas = &c;
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (c(x, y) != seed) {
          decay(c(x, y));
        }
      }
    }

    pull(leader);
    EVERY_N_MILLISECONDS(STOP_TIMER) {
      if (stop < 0) {
        stop = pos[leader];
      }
    }
  }

private:
  /**
   * @brief Ages one pixel by one decay step.
   * @param p Pixel to age, modified in place.
   * @details A still-seed pixel takes a small first step off the seed color; an
   * already-decaying one rides the rainbow trail.
   */
  void decay(Pixel &p) {
    if (p == seed) {
      auto c = rgb2hsv_approximate(static_cast<CRGB>(p));
      c.h = seed.h - 8;
      c.s = 255;
      c.v = qsub8(seed.v, 24);
      p = c;
    } else {
      trail_rainbow(p, 8, 24);
    }
  }

  /**
   * @brief Advances the leader row and drags the rest of the rows behind it.
   * @param y Leader row index.
   * @details Advances leader row y one step, then pulls each other row toward
   * its neighbor, never letting the lag exceed lead_length. Hands off to
   * block() once the leader reaches a latched stop position.
   */
  void pull(int y) {
    if (stop != -1 && pos[y] == stop) {
      block(y);
      return;
    }
    move(y, pos[y], wrap(pos[y] + dir, W), dir);
    for (int i = y - 1; i >= 0; --i) {
      if (distance(pos[i], pos[i + 1], dir) > lead_length) {
        move(i, pos[i], wrap(pos[i + 1] - lead_length * dir, W), dir);
      }
    }
    for (int i = y + 1; i < H; ++i) {
      if (distance(pos[i], pos[i - 1], dir) > lead_length) {
        move(i, pos[i], wrap(pos[i - 1] - lead_length * dir, W), dir);
      }
    }
  }

  /**
   * @brief Drains trailing rows after the leader has parked, then re-randomizes.
   * @param y Leader row index that has reached the latched stop position.
   * @details Keeps advancing every row still trailing the leader until all
   * align, then releases and re-randomizes direction, leader, and lead length
   * for the next twist.
   */
  void block(int y) {
    bool all_stop = true;
    for (int i = 0; i < H; ++i) {
      if (distance(pos[i], pos[y], dir) != 0) {
        move(i, pos[i], wrap(pos[i] + dir, W), dir);
        all_stop = false;
      }
    }
    if (all_stop) {
      stop = -1;
      dir = random(2) < 1 ? 1 : -1;
      leader = random(0, H);
      lead_length = lead_lengths[random(1, sizeof(lead_lengths) / sizeof(int))];
    }
  }

  /**
   * @brief Unsigned ring distance from a to b measured in a given direction.
   * @param a Start column.
   * @param b End column.
   * @param direction Travel direction; positive measures forward, otherwise
   * backward.
   * @return Number of columns from a to b along the ring in that direction.
   */
  int distance(int a, int b, int direction) {
    return direction > 0 ? wrap(b - a, W) : wrap(a - b, W);
  }

  /**
   * @brief Sets a row's offset to x1, painting fresh trails behind each spoke.
   * @param y Row index to move.
   * @param x0 Starting offset of the row.
   * @param x1 New offset of the row.
   * @param direction Travel direction (+1 forward, -1 backward).
   * @details Paints a fresh trail behind each spoke for every step swept from
   * x0 to x1.
   */
  void move(int y, int x0, int x1, int direction) {
    pos[y] = x1;
    for (int x = 0; x < W; x += W / COUNT) {
      for (int d = 0; d < distance(x0, x1, direction); ++d) {
        add_trail(*canvas, y, wrap(x - direction, W), direction);
      }
    }
  }

  /**
   * @brief Recursively extends a decaying trail backward from a spoke.
   * @param c Canvas to draw into.
   * @param y Row index.
   * @param x Column at which to write the copied pixel.
   * @param direction Travel direction (+1 forward, -1 backward).
   * @details Recursively copies the pixel ahead (in direction) backward into x
   * and decays it, extending the trail until it hits a black pixel or the seed.
   */
  void add_trail(Canvas &c, int y, int x, int direction) {
    if (c(wrap(x + direction, W), y) == CRGB(0, 0, 0) || c(x, y) == seed) {
      return;
    }
    add_trail(c, y, wrap(x - direction, W), direction);
    c(x, y) = c(wrap(x + direction, W), y);
    decay(c(x, y));
  }

  NoTempCorrection _;  /**< Disables temperature correction for this effect. */
  const CHSV seed;     /**< Color of the seed spokes. */
  const int COUNT = 2; /**< Number of seed spokes around the ring. */
  const int STOP_TIMER = 15000; /**< Interval (ms) between stop latches. */
  int pos[H] = {0};             /**< Per-row scroll offset. */
  int leader = 0; /**< Index of the row currently driving the twist. */
  int dir = 1;    /**< Current twist direction (+1 or -1). */
  int stop = -1;  /**< Latched stop offset, or -1 when running. */
  int lead_length =
      1; /**< Maximum lag (columns) a row may trail its neighbor. */
  int lead_lengths[7] = {1, 2, 3, 4, 6, 12, 16}; /**< Candidate lead lengths. */
  Canvas *canvas = NULL; /**< Canvas for the frame in progress. */
};

/**
 * @brief Falling "digital rain" columns.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam HUE Base hue of the rain gradient, in 8-bit hue units.
 * @details Each pixel holds an intensity that scrolls down, mapped to a
 * HUE-based gradient (bright head, dim tail).
 */
template <int W, int H, uint8_t HUE> class TheMatrix : public Effect {
public:
  /**
   * @brief Constructs the effect and clears the intensity field.
   */
  FLASHMEM TheMatrix() : Effect(W, H, {.strobe = true}) {
    random16_add_entropy(random());
    memset(pixels, 0, sizeof(pixels));
  }

  /**
   * @brief Renders one frame of digital rain.
   * @details Every 125 ms scrolls each column down, seeds new drops, and
   * renders.
   */
  void draw_frame() {
    Canvas c(*this);
    EVERY_N_MILLIS(125) {
      for (int x = 0; x < W; ++x) {
        fall(x);
        generate(x);
        for (int y = 0; y < H; ++y) {
          c(x, y) =
              blend(CHSV(HUE + (3 * y), 255, 40),
                    CHSV(HUE + (3 * (H - 1 - y)), 255, 255), pixels[x][y]);
        }
      }
    }
  }

private:
  /**
   * @brief Shifts a column down one row, clearing the top.
   * @param x Column index to scroll.
   */
  void fall(int x) {
    memmove8(&pixels[x][1], &pixels[x][0], H - 1);
    pixels[x][0] = 0;
  }

  /**
   * @brief Randomly spawns a new drop at the top of a column.
   * @param x Column index to maybe seed.
   * @details Spawns a bright head fading over three rows.
   */
  void generate(int x) {
    if (random8() < 15) {
      pixels[x][2] = 255;
      pixels[x][1] = 100;
      pixels[x][0] = 30;
    }
  }

  uint8_t pixels[W][H]; /**< Per-cell rain intensity. */
  NoColorCorrection _;  /**< Disables color correction for this effect. */
};

/**
 * @brief Mirror-paired curved strokes sweeping the ring.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each pixel's brightness is remapped through a warm gradient palette
 * so the fading trails shift hue with their decay.
 */
template <int W, int H> class Curves : public Effect {
public:
  /**
   * @brief Constructs the effect and builds the warm gradient palette.
   */
  FLASHMEM Curves() : Effect(W, H, {.persist = true}) {
    fill_gradient<CHSV>(palette, sizeof(palette) / sizeof(CHSV),
                        rgb2hsv_approximate(CRGB(6, 4, 47)),
                        rgb2hsv_approximate(CRGB(162, 84, 84)),
                        rgb2hsv_approximate(CRGB(252, 114, 0)), SHORTEST_HUES);
  }

  /**
   * @brief Renders one frame of mirrored sweeping strokes.
   * @details Fades trails and remaps their hue by brightness, then draws COUNT
   * mirrored stroke pairs and advances the slope, color, and ring offset.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 0, 32);
        CHSV h = rgb2hsv_approximate(static_cast<CRGB>(c(x, y)));
        h.h = palette[h.v].h;
        c(x, y) = h;
      }
    }

    for (int i = 0; i < COUNT; ++i) {
      CHSV color = palette[map8(triwave8(color_offset), 100, 140)];
      plot(c, num, 255, (i * W / COUNT + offset) % W, color);
      plot(c, -num, 255, (i * W / COUNT + W - 1 - offset) % W, color);
    }

    num = wrap(num - 1, 1024);
    color_offset++;
    offset = (offset + 1) % W;
  }

private:
  /**
   * @brief Draws a straight stroke of slope n/d anchored at a column.
   * @param cv Canvas to draw into.
   * @param n Slope numerator.
   * @param d Slope denominator.
   * @param c Anchor column at row 0.
   * @param color Color to draw the stroke.
   * @details The column advances y*n/d per row, wrapping around the ring.
   */
  void plot(Canvas &cv, int n, int d, int c, const CHSV &color) {
    for (int y = 0; y < H; ++y) {
      cv(wrap(y * n / d + c, W), y) = color;
    }
  }

  int COUNT = 4;            /**< Number of mirrored stroke pairs. */
  int offset = 0;           /**< Current ring offset of the strokes. */
  uint8_t color_offset = 0; /**< Phase into the palette color cycle. */
  uint16_t num = 1024;      /**< Current stroke slope numerator. */
  uint8_t falloff = 48;     /**< Trail fade rate. */
  CHSVPalette256 palette;   /**< Warm gradient palette. */
  NoTempCorrection _; /**< Disables temperature correction for this effect. */
};

/**
 * @brief Sparse stars that wink on and fade out along rainbow trails.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Stars wink on at the current hue and fade out along rainbow trails;
 * the spawn hue drifts one step per frame.
 */
template <int W, int H> class StarsFade : public Effect {
private:
  uint8_t hue; /**< Hue used to spawn the next star; drifts per frame. */

public:
  /**
   * @brief Constructs the effect with persistent pixels.
   */
  FLASHMEM StarsFade()
      : Effect(W, H, {.strobe = true, .persist = true}), hue(0) {}

  /**
   * @brief Renders one frame of fading stars.
   * @details Fades lit pixels; rarely lights a dark pixel with a fresh star.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (c(x, y) != CRGB(0, 0, 0)) {
          trail_rainbow(c(x, y), 16, 32);
        } else if (random8() < 3) {
          c(x, y) = CRGB(CHSV(hue, 255, 255));
        }
      }
    }
    hue++;
  }
};

/**
 * @brief Static rainbow spiral of diagonal arms.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam SPREAD Spacing control; every (SPREAD+1)th diagonal cell is lit.
 * @details Lights every (SPREAD+1)th cell along the x+y diagonal, coloring it
 * by position so the lit cells form diagonal spiral arms.
 */
template <int W, int H, int SPREAD> class Spiral : public Effect {
private:
  CHSVPalette16 palette; /**< Rainbow palette indexed by diagonal position. */
  NoTempCorrection _; /**< Disables temperature correction for this effect. */

public:
  /**
   * @brief Constructs the effect, fills the palette, and paints the pattern.
   */
  FLASHMEM Spiral() : Effect(W, H, {.strobe = true}) {
    fill_rainbow(palette.entries, 16, 0, 256 / 16);
    draw_frame();
  }

  /**
   * @brief Paints the diagonal spiral pattern.
   * @details Output is identical every frame.
   */
  FASTRUN void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if ((x + y) % (SPREAD + 1) == 0) {
          c(x, y) = palette[(x + y) % 16];
        } else {
          c(x, y) = CHSV(0, 0, 0);
        }
      }
    }
  }
};

/**
 * @brief Three sine-wave traces sweeping up and down each column.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The blue, green, and aqua traces sweep at phase offsets set per
 * column, leaving rainbow trails behind them.
 */
template <int W, int H> class WaveTrails : public Effect {
public:
  /**
   * @brief Constructs the effect with persistent pixels.
   */
  FLASHMEM WaveTrails() : Effect(W, H, {.persist = true}) {}

  /**
   * @brief Renders one frame of the three wave trails.
   * @details Fades trails, then plots each column's three wave samples.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow(c(x, y), 5, 24);
      }
      c(x, beatsin16(5, H / 15, H * 14 / 15, 0, map(x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_BLUE, 100, 255);
      c(x, beatsin16(5, 0, H - 1, 0, map(W - 1 - x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_GREEN, 100, 255);
      c(x, beatsin16(5, 0, H - 1, 32767, map(W - 1 - x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_AQUA, 100, 255);
    }
  }
};

/**
 * @brief A horizontal ring at mid-height tumbled through a 3D rotation.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The ring is projected through a tumbling 3D rotation and drawn
 * anti-aliased so it sweeps across the cylinder leaving rainbow trails; the
 * rotation rates beat on three independent sine cycles.
 */
template <int W, int H> class RingTrails : public Effect {
public:
  /**
   * @brief Constructs the effect, builds the palette, and clears the TTL field.
   */
  FLASHMEM RingTrails() : Effect(W, H, {.persist = true}), dot(0) {
    fill_gradient<CHSV>(palette, sizeof(palette) / sizeof(CHSV),
                        rgb2hsv_approximate(CRGB(6, 4, 47)),
                        rgb2hsv_approximate(CRGB(162, 84, 84)),
                        rgb2hsv_approximate(CRGB(252, 114, 0)), SHORTEST_HUES);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        ttl[x][y] = 0;
      }
    }
  }

  /**
   * @brief Renders one frame of the tumbling ring.
   * @details Fades trails, plots the projected ring, then advances the rotation
   * one step.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 10, 24);
      }
    }

    for (float x = 0; x < W; x += 1) {
      if (static_cast<int>(x) % 2 == 0) {
      }
      Point p = projection.project(x, H / 2);
      CHSV color;
      color = CHSV(HUE_YELLOW, 0, 255);
      plot_aa(c, p.x, p.y, color);
    }

    uint8_t dl = beatsin8(1, 1, 3);
    uint8_t dg = beatsin8(2, 1, 3, 16384);
    uint8_t dp = beatsin8(3, 1, 3, 32768);
    projection.rotate(dl, dg, dp);
    dot = (dot + 1) % W;
  }

private:
  typedef typename Projection<W, H>::Point Point; /**< Projected point type. */
  Projection<W, H> projection; /**< 3D rotation/projection of the ring. */
  uint8_t hue = HUE_RED;       /**< Base hue. */
  CHSVPalette256 palette;      /**< Warm gradient palette. */
  int ttl[W][H];               /**< Per-cell time-to-live for trails. */
  int dot;                     /**< Current ring scan position. */
  NoColorCorrection _; /**< Disables color correction for this effect. */
};

/**
 * @brief Point-symmetric kaleidoscope of palette-colored stamps.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each frame dims the canvas, then stamps a palette-colored point and
 * its central-symmetry mirror, scanning the stamp position across the ring (x)
 * and bouncing it up and down (y). The number of symmetric arms cycles every
 * few seconds.
 */
template <int W, int H> class Kaleidoscope : public Effect {
public:
  /**
   * @brief Constructs the effect and fills the rainbow palette.
   */
  FLASHMEM Kaleidoscope() : Effect(W, H, {.persist = true}) {
    fill_rainbow(palette.entries, 16, 0, 256 / 16);
  }

  /**
   * @brief Renders one frame of the kaleidoscope.
   * @details Fades the canvas, draws the symmetric arms, then steps the scan
   * position and (at the vertical turnaround) advances the arm count and
   * palette.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = c(x, y) * (159.0f / 255.0f);
      }
    }

    for (uint8_t i = 0; i < counts[count]; ++i) {
      uint8_t x = addmod8(i * offset, tx, W);
      uint8_t y = ty;
      uint8_t x2 = W - 1 - x;
      uint8_t y2 = H - 1 - y;
      c(x, y) = palette[addmod8(palette_offset, i * (16 / counts[count]), 16)];
      c(x2, y2) =
          palette[addmod8(palette_offset, 16 - (i * (16 / counts[count])), 16)];
    }

    tx = addmod8(tx, 2, W);
    ty = addmod8(ty, inc_y, H);
    if (ty == H - 1 || ty == 0) {
      inc_y = 0 - inc_y;
      tx = addmod8(tx, 2, W);
      EVERY_N_MILLISECONDS(3000) {
        count = addmod8(count, 1, sizeof(counts));
        offset = W / counts[count];
      }
    }
    palette_offset = addmod8(palette_offset, 1, 16);
  }

private:
  uint8_t counts[9] = {1, 2, 3, 4, 6, 8, 12, 16, 20}; /**< Arm-count cycle. */
  uint8_t count = 0;                                  /**< Index into counts. */
  uint8_t offset = W / counts[count]; /**< Angular spacing between arms. */
  uint8_t tx = 0;                     /**< Horizontal scan position. */
  uint8_t ty = 0;                     /**< Vertical scan position. */
  int inc_y = 1;              /**< Vertical scan direction (+1 or -1). */
  CHSVPalette16 palette;      /**< Rainbow palette. */
  uint8_t palette_offset = 0; /**< Phase into the palette color cycle. */
  NoTempCorrection _; /**< Disables temperature correction for this effect. */
};

/**
 * @brief Three rainbow rings stacked up the cylinder, each tumbled in 3D.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each ring is tumbled through one of two 3D projections (alternating
 * rings counter-rotate) and drawn anti-aliased. The rotation rates swap every
 * 10 s; the rainbow also scrolls around each ring.
 */
template <int W, int H> class RingRotate : public Effect {
public:
  /**
   * @brief Constructs the effect and fills the rainbow palette.
   */
  FLASHMEM RingRotate() : Effect(W, H) {
    fill_rainbow(pal.entries, 256, HUE_RED, 1);
  }

  /**
   * @brief Renders one frame of the three tumbling rings.
   * @details Clears, draws the three projected rings, advances the rotations
   * and the rainbow scroll offset, and cycles the rotation-rate set every 10 s.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = CHSV(0, 0, 0);
      }
    }

    const float rings[] = {4.5f, 9.5f, 15.5f};
    for (size_t i = 0; i < sizeof(rings) / sizeof(float); ++i) {
      float y = rings[i];
      for (float x = 0; x < W; ++x) {
        Point p(x, y);
        int dir = 1;
        if (i % 2 == 0) {
          dir = 0;
          p = projection.project(p);
        } else {
          p = projection2.project(p);
        }
        CHSV color = pal[(static_cast<uint8_t>(x) + c_off * dir) * 255 / W];
        plot_aa(c, p.x, p.y, color);
      }
    }

    EVERY_N_SECONDS(10) {
      r_off = (r_off + 1) % (sizeof(rotations) / sizeof(uint16_t) / 3 / 2);
    }

    projection.rotate(rotations[r_off][0][0], rotations[r_off][0][1],
                      rotations[r_off][0][2]);

    projection2.rotate(rotations[r_off][1][0], rotations[r_off][1][1],
                       rotations[r_off][1][2]);

    c_off += 2;
  }

private:
  typedef typename Projection<W, H>::Point Point; /**< Projected point type. */
  Projection<W, H> projection;  /**< First projection (even rings). */
  Projection<W, H> projection2; /**< Second projection (odd rings). */

  CHSVPalette256 pal; /**< Rainbow palette. */
  uint8_t c_off = 0;  /**< Rainbow scroll offset around each ring. */
  bool bg = 0;        /**< Background-draw flag. */

  /** @brief Rotation-rate sets: [set][projection][axis], swapped every 10 s. */
  uint16_t rotations[4][2][3] = {
      {{0, 5, 0}, {0, 5, 0}},
      {{0, 1, 3}, {0, 2, 5}},
      {{0, 5, 0}, {0, 355, 0}},
      {{0, 355, 0}, {0, 3, 0}},
  };

  uint8_t r_off = 0; /**< Index of the active rotation-rate set. */
};

/**
 * @brief Pixels ignite one at a time in random order, then fall as embers.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam HUE Base hue of the glow, in 8-bit hue units.
 * @tparam BURNRATE Ignition-rate control for the burn order.
 * @details Each lit pixel glows bright at hue HUE, then its ember falls down
 * its column row by row until it drops off the bottom, leaving the field dark
 * behind it.
 */
template <int W, int H, uint8_t HUE, uint8_t BURNRATE>
class Burnout : public Effect {
public:
  /**
   * @brief Constructs the effect and shuffles the pixel ignition order.
   */
  FLASHMEM Burnout() : Effect(W, H, {.strobe = true}), timer(30), burn_idx(0) {
    random16_add_entropy(random());
    memset8(pixels, INIT, sizeof(pixels));
    fill_seq(burn, W * H);
    shuffle(burn, W * H);
  }

  /**
   * @brief Renders one frame of the burnout animation.
   * @details Renders each pixel by state (burning bright, unlit dim, or
   * burnt-out dark), lets embers fall, and ignites the next pixel in the
   * shuffled order.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (pixels[x][y] & BURN) {
          c(x, y) = CHSV(HUE + (3 * y), 255, 255);
        } else if (pixels[x][y] == INIT) {
          c(x, y) = CHSV(HUE + (3 * y), 255, 40);
        } else {
          c(x, y) = CHSV(0, 0, 0);
        }
      }
      fall(x);
    }
    if (timer) {
      burnout();
    }
  }

private:
  /**
   * @brief Per-pixel state.
   * @details NONE is burnt-out dark, INIT is unlit dim, BURN (high bit) is
   * currently burning.
   */
  enum State {
    NONE,
    INIT,
    BURN = 0x80,
  };

  /**
   * @brief Ignites the next pixel in the shuffled burn order.
   */
  void burnout() {
    if (burn_idx < W * H) {
      *((uint8_t *)pixels + burn[burn_idx++]) = BURN;
    }
  }

  /**
   * @brief Drops every burning ember in a column down one row.
   * @param x Column index to advance.
   * @details An ember dropping off the bottom clears it.
   */
  void fall(int x) {
    for (int y = 0; y < H; ++y) {
      if (pixels[x][y] & BURN) {
        pixels[x][y] ^= BURN;
        if (y < (H - 1)) {
          pixels[x][y + 1] |= BURN;
          ++y;
        }
      }
    }
  }

  /**
   * @brief Fills an array with the sequence 0..len-1.
   * @param arr Destination array of flat pixel indices to be burned.
   * @param len Number of entries to fill.
   */
  void fill_seq(short *arr, int len) {
    for (int i = 0; i < len; ++i) {
      arr[i] = i;
    }
  }

  /**
   * @brief In-place Fisher-Yates shuffle so pixels ignite in random order.
   * @param arr Array to shuffle, modified in place.
   * @param len Number of entries in arr.
   */
  void shuffle(short *arr, int len) {
    for (int i = 0; i < len - 1; ++i) {
      short j = random(i, len);
      short temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }

  CEveryNMillis timer; /**< Gates how often a new pixel ignites. */

  uint8_t pixels[W][H]; /**< Per-pixel State. */
  short burn[W * H];    /**< Shuffled flat pixel indices, burn order. */
  int burn_idx;         /**< Index of the next pixel to ignite. */
};

/**
 * @brief Classic per-column heat-diffusion fire.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam COOL Cooling-rate control; larger values cool faster.
 * @tparam SPARK Spark probability numerator over 255 for igniting the base.
 * @details Each column cools, heat rises by averaging the cells below, sparks
 * randomly ignite the base, and the heat map is rendered through FastLED's
 * HeatColor palette.
 */
template <int W, int H, uint8_t COOL, uint8_t SPARK>
class Fire : public Effect {
public:
  /**
   * @brief Constructs the effect and seeds the entropy pool.
   */
  FLASHMEM Fire() : Effect(W, H) { random16_add_entropy(random()); }

  /**
   * @brief Renders one frame of fire.
   * @details Every 125 ms cools, rises, and sparks each column, then maps heat
   * to color.
   */
  void draw_frame() {
    EVERY_N_MILLIS(125) {
      Canvas c(*this);
      for (int x = 0; x < W; ++x) {
        cool(x);
        rise(x);
        spark(x);
        for (int y = 0; y < H; ++y) {
          c(x, y) = rgb2hsv_approximate(HeatColor(heat[x][y]));
        }
      }
    }
  }

private:
  /**
   * @brief Subtracts a small random amount of heat from every cell in a column.
   * @param x Column index to cool.
   */
  inline void cool(uint8_t x) {
    for (uint8_t y = 0; y < H; ++y) {
      heat[x][y] = qsub8(heat[x][y], random(0, ((COOL * 10) / H) + 2));
    }
  }

  /**
   * @brief Propagates heat upward by averaging each cell with the two below it.
   * @param x Column index to update.
   */
  inline void rise(uint8_t x) {
    for (uint8_t y = H - 1; y >= 2; --y) {
      heat[x][H - 1 - y] = (heat[x][H - 1 - y + 1] + heat[x][H - 1 - y + 2] +
                            heat[x][H - 1 - y + 2]) /
                           3;
    }
  }

  /**
   * @brief Randomly injects a hot spark into one of the bottom rows.
   * @param x Column index to maybe spark.
   * @details Fires with probability SPARK/255.
   */
  inline void spark(uint8_t x) {
    if (random8() < SPARK) {
      uint8_t y = random8(3);
      heat[x][H - 1 - y] = qadd8(heat[x][H - 1 - y], random8(160, 255));
    }
  }

  uint8_t heat[W][H]; /**< Per-cell heat map. */
};

/**
 * @brief One white dot per row circling the ring with fading trails.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Alternating rows run in opposite directions; each dot leaves a
 * fading rainbow trail and teleports to a new random spot once its
 * time-to-live expires.
 */
template <int W, int H> class DotTrails : public Effect {
public:
  /**
   * @brief Constructs the effect and seeds one dot per row.
   */
  FLASHMEM DotTrails() : Effect(W, H, {.persist = true}) {
    random16_add_entropy(random());
    for (int i = 0; i < H - 1; ++i) {
      bool rev = i % 2 == 0;
      uint8_t x = random8() % W;
      dots[i] = Dot(x, rev, (random8() % 30) + 20);
    }
  }

  /**
   * @brief Renders one frame of the circling dots.
   * @details Draws each row's dot and advances it; ages every other cell along
   * its trail.
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (dots[y].x == x) {
          if (dots[y].drawn) {
            dots[y].drawn = false;
          } else {
            c(x, y) = CHSV(hue - 12, 0, 255);
            move(dots[y]);
          }
        } else {
          trail_rainbow(c(x, y), 12, 12);
        }
      }
    }
  }

private:
  /**
   * @brief A moving dot circling the ring on one row.
   * @details Holds ring position x, direction (rev), frames-to-live ttl, and a
   * drawn flag that defers one frame so a freshly painted cell isn't
   * re-trailed.
   */
  struct Dot {
    /**
     * @brief Constructs a zeroed dot.
     */
    Dot() : x(0), rev(0), ttl(0), drawn(0) {}

    /**
     * @brief Constructs a dot with a position, direction, and lifetime.
     * @param x Ring column position.
     * @param rev True to travel in reverse (backward) around the ring.
     * @param ttl Frames-to-live before the dot respawns.
     */
    Dot(uint8_t x, bool rev, int ttl) : x(x), rev(rev), ttl(ttl), drawn(0) {}

    uint8_t x; /**< Ring column position. */
    bool rev;  /**< True when travelling backward around the ring. */
    int ttl;   /**< Frames remaining before respawn. */
    bool
        drawn; /**< Defers one frame so a freshly painted cell isn't re-trailed. */
  };

  /**
   * @brief Steps a dot one cell, respawning it when its lifetime expires.
   * @param dot Dot to advance, modified in place.
   * @details When ttl hits 0, respawns it at a random position with a fresh
   * lifetime; otherwise steps one cell respecting direction.
   */
  void move(Dot &dot) {
    if (dot.ttl == 0) {
      dot.x = random8() % W;
      dot.ttl = (random8() % 30) + 20;
    } else if (dot.rev) {
      dot.x = (dot.x - 1 + W) % W;
      dot.ttl--;
    } else {
      dot.x = addmod8(dot.x, 1, W);
      dot.drawn = true;
      dot.ttl--;
    }
  }

  Dot dots[H];           /**< Per-row dots. */
  uint8_t hue = HUE_RED; /**< Hue used to draw the dot heads. */
  NoTempCorrection _; /**< Disables temperature correction for this effect. */
};

/**
 * @brief Horizontal bands of two two-color palettes scrolling oppositely.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The red/yellow and aqua/blue palettes scroll oppositely around the
 * ring; the band thickness cycles through p[] every few seconds and the two
 * palettes periodically swap rows.
 */
template <int W, int H> class Spinner : public Effect {
public:
  /**
   * @brief Constructs the effect and seeds the two two-color palettes.
   */
  FLASHMEM Spinner() : Effect(W, H), spin_timer(62), pos(0) {
    memset(palette1, 0, sizeof(palette1));
    memset(palette2, 0, sizeof(palette1));
    palette1[0] = palette1[8] = CHSV(HUE_RED, 255, 255);
    palette1[1] = palette1[9] = CHSV(HUE_YELLOW, 255, 255);
    palette2[0] = palette2[8] = CHSV(HUE_AQUA, 255, 255);
    palette2[1] = palette2[9] = CHSV(HUE_BLUE, 255, 255);
  }

  /**
   * @brief Renders one frame of the scrolling bands.
   * @details Paints the two scrolling band sets, advances the scroll, and every
   * 5 s steps the band-thickness index (reversing and swapping palettes at the
   * ends).
   */
  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (!swap) {
          if (y % (p[i] * 2) < p[i]) {
            c(x, y) = palette1[addmod8(x, pos, 16)];
          } else {
            c(x, y) = palette2[(uint8_t)(pos - x) % 16];
          }
        } else {
          if (y % (p[i] * 2) < p[i]) {
            c(x, y) = palette2[(uint8_t)(pos - x) % 16];
          } else {
            c(x, y) = palette1[addmod8(x, pos, 16)];
          }
        }
      }
      if (spin_timer) {
        pos = addmod8(pos, 1, W);
      }
    }
    EVERY_N_SECONDS(5) {
      if ((!swap && (i == sizeof(p) - 1)) || (swap && (i == 0))) {
        swap = !swap;
      } else {
        i += swap ? -1 : 1;
      }
    }
  }

private:
  CHSVPalette16 palette1;   /**< Red/yellow band palette. */
  CHSVPalette16 palette2;   /**< Aqua/blue band palette. */
  CEveryNMillis spin_timer; /**< Gates the scroll advance. */
  uint8_t pos;              /**< Current scroll offset. */
  bool swap = false;        /**< When true, the two palettes swap rows. */
  uint8_t p[5] = {20, 10, 4, 2, 1}; /**< Band-thickness cycle. */
  uint8_t i = 0;                    /**< Index into p[]. */
  NoColorCorrection _; /**< Disables color correction for this effect. */
};
