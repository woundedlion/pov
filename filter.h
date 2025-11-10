#pragma once
template <int W>
class Filter {
public:
  Filter() : next(nullptr) {}

  virtual Filter& chain(Filter& filter) {
    next = &filter;
    return *this;
  }

  virtual void plot(Canvas& canvas, float_t x, float_t y,
    const Pixel& c, float_t age, float_t alpha) = 0;

protected:

  void pass(Canvas& canvas, float_t x, float_t y,
    const Pixel& c, float_t age, float_t alpha)
  {
    if (next == nullptr) {
      auto xy = XY(static_cast<int>(x + 0.5f), static_cast<int>(y + 0.5f));
      auto p = blend_alpha(alpha)(canvas(xy), c);
      canvas(xy) = p;
    } else {
      next->plot(canvas, x, y, c, age, alpha);
    }
  }

  Filter* next;
};

template <int W>
class FilterRaw : public Filter<W> {
  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& c, float_t age, float_t alpha) {
    this->pass(canvas, x, y, c, age, alpha);
  }
};

template <int W>
class FilterAntiAlias : public Filter<W> {
public:
  FilterAntiAlias() {}

  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& c, float_t age, float_t alpha) {
    float_t x_i = 0;
    float_t x_m = modf(x, &x_i);
    float_t y_i = 0;
    float_t y_m = modf(y, &y_i);

    float_t v = (1 - x_m) * (1 - y_m);
    this->pass(canvas, x_i, y_i, c, age, alpha * v);

    v = x_m * (1 - y_m);
    this->pass(canvas, (static_cast<int>(x_i + 1)) % W, y_i, c, age, alpha * v);

    if (y_i < H - 1) {
      v = (1 - x_m) * y_m;
      this->pass(canvas, x_i, y_i + 1, c, age, alpha * v);

      v = x_m * y_m;
      this->pass(canvas, static_cast<int>(x_i + 1) % W, y_i + 1, c, age, alpha * v);
    }
  }
};

template <int W, int MAX_PIXELS>
class FilterDecay : public Filter<W> {
public:

  FilterDecay(int lifetime) :
    lifetime(lifetime),
    num_pixels(0)
  {}

  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls[i].ttl < 0.00001) {
        num_pixels--;
        if (i < num_pixels) {
          ttls[i] = std::move(ttls[num_pixels]);
          i--;
        }
      }
    }
  }

  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& color, float_t age, float_t alpha) {
    if (age >= 0) {
      if (num_pixels < MAX_PIXELS) {
        ttls[num_pixels++] = { static_cast<float>(x), static_cast<float>(y), static_cast<float>(lifetime - age) };
      }
      else {
        Serial.println("FilterDecay full!");
      }
    }
    if (age <= 0) {
      this->pass(canvas, x, y, color, age, alpha);
    }
  }

  void trail(Canvas& canvas, TrailFn trailFn, float_t alpha) {
    for (int i = 0; i < num_pixels; ++i) {
      auto color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      this->pass(canvas, ttls[i].x, ttls[i].y, color, lifetime - ttls[i].ttl, alpha);
    }
  }

private:

  struct DecayPixel {
    float_t x;
    float_t y;
    float_t ttl;
  };

  int lifetime;
  std::array<DecayPixel, MAX_PIXELS> ttls;
  int num_pixels;
};

template<int W>
class FilterChromaticShift : public Filter<W> {
public:

  FilterChromaticShift()
  {
  }

  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& color, float_t age, float_t alpha) {
    CRGB r(color.r, 0, 0);
    CRGB g(0, color.g, 0);
    CRGB b(0, 0, color.b);
    this->pass(canvas, x, y, color, age, alpha);
    this->pass(canvas, wrap(x + 1, W), y, r, age, alpha);
    this->pass(canvas, wrap(x + 2, W), y, g, age, alpha);
    this->pass(canvas, wrap(x + 3, W), y, b, age, alpha);
  }
};

template <int W>
class FilterOrient : public Filter<W> {
public:
  FilterOrient(Orientation& orientation) :
    orientation(orientation)
  {}

  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& color, float_t age, float_t alpha) {
    auto v = pixel_to_vector<W>(x, y);
    auto r = vector_to_pixel<W>(orientation.orient(v));
    orientation.collapse();
    pass(canvas, r.x, r.y, color, age, alpha);
  }

private:

  Orientation& orientation;
};

template <int W>
class FilterReplicate : public Filter<W> {
public:

  FilterReplicate(int count) :
    count(std::clamp(count, 1, W))
  {}

  void plot(Canvas& canvas, float_t x, float_t y, const Pixel& color, float_t age, float_t alpha) {
    for (int i = 0; i < W; i += W / count) {
      this->pass(canvas, wrap(x + i, W), y, color, age, alpha);
    }
  }

private:

  int count;
};
