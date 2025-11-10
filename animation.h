#pragma once
static const int FPS = 16;

float_t ease_in_out_bicubic(float_t t) {
  return t < 0.5 ? 4 * pow(t, 3) : 1 - pow(-2 * t + 2, 3) / 2;
}

float_t ease_in_out_sin(float_t t) {
  return -(cos(PI_F * t) - 1) / 2;
}

float_t ease_in_sin(float_t t) {
  return 1 - cos((t * PI_F) / 2);
}

float_t ease_out_sin(float_t t) {
  return sin((t * PI_F) / 2);
}

float_t ease_in_cubic(float_t t) {
  return pow(t, 3);
}

float_t ease_in_circ(float_t t) {
  return 1 - sqrt(1 - pow(t, 2));
}

float_t ease_mid(float_t t) {
  return t;
}

float_t ease_out_expo(float_t t) {
  return t == 1 ? 1 : 1 - pow(2, -10 * t);
}

float_t ease_out_circ(float_t t) {
  return sqrt(1 - pow(t - 1, 2));
}

template <typename Derived>
class Animation {
public:

  void cancel() { canceled = true; }
  virtual bool done() { return canceled || (duration >= 0 && t >= duration); }

  virtual void step(Canvas& canvas) {
    t++;
    if (done()) {
      if (repeat) {
        t = 0;
      }
    }
  }

  Derived& then(std::function<void()> callback) {
    post = callback;
    return static_cast<Derived&>(*this);
  }

  void post_callback() {
    post();
  }

protected:

  Animation(int duration, bool repeat) :
    duration(duration),
    repeat(repeat),
    canceled(false),
    post([]() {})
  {
  }

  int duration;
  bool repeat;
  int t = 0;

private:

  bool canceled;
  std::function<void()> post;
};

class NullAnimation : public Animation<NullAnimation> {
public:
  NullAnimation() :
    Animation(0, false)
  {}
  void step(Canvas&) {}
  bool done() { return true; }
};

class RandomTimer : public Animation<RandomTimer> {
public:

  RandomTimer(int min, int max, TimerFn f, bool repeat = false) :
    Animation(-1, repeat),
    min(min),
    max(max),
    f(f),
    next(0)
  {
    reset();
  }

  void reset() {
    next = t + static_cast<int>(std::round(hs::rand_int(min, max)));
  }

  void step(Canvas& canvas) {
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step(canvas);
  }

private:

  int min;
  int max;
  TimerFn f;
  int next;
};

class PeriodicTimer : public Animation<PeriodicTimer> {
public:
  PeriodicTimer(int period, TimerFn f, bool repeat = false) :
    Animation(-1, repeat),
    period(period),
    f(f)
  {
    reset();
  }

  void reset() {
    next = t + period;
  }

  void step(Canvas& canvas) {
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step(canvas);
  }

private:

  int period;
  TimerFn f;
  int next;
};

class Transition : public Animation<Transition> {
public:
  Transition(float_t& mutant, float_t to, int duration, EasingFn easing_fn, bool quantized = false, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    to(to),
    easing_fn(easing_fn),
    quantized(quantized)
  {
  }

  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    Animation::step(canvas);
    auto t = std::min(1.0f, static_cast<float_t>(this->t) / duration);
    auto n = easing_fn(t) * (to - from) + from;
    if (quantized) {
      n = std::floor(n);
    }
    mutant.get() = n;
  }

  void rebind_mutant(float_t& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<float_t> mutant;
  float_t from;
  float_t to;
  EasingFn easing_fn;
  bool quantized;
};

class Mutation : public Animation<Mutation> {
public:

  Mutation(float_t& mutant, MutateFn f, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    f(f),
    easing_fn(easing_fn)
  {}

  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    Animation::step(canvas);
    auto t = std::min(1.0f, static_cast<float_t>(this->t) / duration);
    mutant.get() = f(easing_fn(t));
  }

  void rebind_mutant(float_t& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<float_t> mutant;
  float_t from;
  MutateFn f;
  EasingFn easing_fn;
};

class Sprite : public Animation<Sprite> {
public:

  Sprite(SpriteFn draw_fn, int duration,
    int fade_in_duration = 0, EasingFn fade_in_easing_fn = ease_mid,
    int fade_out_duration = 0, EasingFn fade_out_easing_fn = ease_mid) :
    Animation(duration, false),
    draw_fn(draw_fn),
    fader(fade_in_duration > 0 ? 0 : 1),
    fade_in_duration(fade_in_duration),
    fade_out_duration(fade_out_duration),
    fade_in(fader, 1, fade_in_duration, fade_in_easing_fn),
    fade_out(fader, 0, fade_out_duration, fade_out_easing_fn)
  {}

  Sprite(Sprite&& other) noexcept
    : Animation(std::move(other)),
    draw_fn(std::move(other.draw_fn)),
    fader(other.fader),
    fade_in_duration(other.fade_in_duration),
    fade_out_duration(other.fade_out_duration),
    fade_in(std::move(other.fade_in)),
    fade_out(std::move(other.fade_out))
  {
    fade_in.rebind_mutant(this->fader);
    fade_out.rebind_mutant(this->fader);
  }

  Sprite& operator=(Sprite&& other) noexcept {
    if (this == &other) return *this;

    Animation::operator=(std::move(other));
    draw_fn = std::move(other.draw_fn);
    fader = other.fader;
    fade_in_duration = other.fade_in_duration;
    fade_out_duration = other.fade_out_duration;
    fade_in = std::move(other.fade_in);
    fade_in.rebind_mutant(this->fader);
    fade_out = std::move(other.fade_out);
    fade_out.rebind_mutant(this->fader);

    return *this;
  }

  void step(Canvas& canvas) {
    if (!fade_in.done()) {
      fade_in.step(canvas);
    }
    else if (duration >= 0 && t >= (duration - fade_out_duration)) {
      fade_out.step(canvas);
    }
    draw_fn(canvas, fader);
    Animation::step(canvas);
  }

private:

  SpriteFn draw_fn;
  float_t fader;
  int fade_in_duration;
  int fade_out_duration;
  Transition fade_in;
  Transition fade_out;
};

template <int W>
class Motion : public Animation<Motion<W>> {
public:

  Motion(Orientation& orientation, const Path<W>& path, int duration, bool repeat = false) :
    Animation<Motion<W>>(duration, repeat),
    orientation(orientation),
    path(path),
    to(this->path.get_point(0))
  {}

  void step(Canvas& canvas) {
    Animation<Motion<W>>::step(canvas);
    orientation.get().collapse();
    from = to;
    to = path.get().get_point(static_cast<float_t>(this->t) / this->duration);
    if (from != to) {
      Vector axis = cross(from, to).normalize();
      auto angle = angle_between(from, to);
      auto step_angle = angle / std::ceil(angle / MAX_ANGLE);
      auto& origin = orientation.get().get();
      for (auto a = step_angle; angle - a > 0.0001; a += step_angle) {
        orientation.get().push(make_rotation(axis, a) * origin);
      }
      orientation.get().push(make_rotation(axis, angle) * origin);
    }
  }

private:

  static constexpr float_t MAX_ANGLE = 2 * PI_F / W;
  std::reference_wrapper<Orientation> orientation;
  std::reference_wrapper<const Path<W>> path;
  Vector from;
  Vector to;
};

template <int W>
class Rotation : public Animation<Rotation<W>> {
public:

  Rotation(Orientation& orientation, const Vector& axis, float_t angle, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation<Rotation<W>>(duration, repeat),
    orientation(orientation),
    axis(axis),
    total_angle(angle),
    easing_fn(easing_fn),
    from(0),
    to(0)
  {
  }

  void step(Canvas& canvas) {
    Animation<Rotation<W>>::step(canvas);
    orientation.get().collapse();
    from = to;
    to = easing_fn(static_cast<float_t>(this->t) / this->duration) * total_angle;
    auto angle = fwd_distance(from, to, total_angle);
    if (angle > 0.00001) {
      auto step_angle = angle / std::ceil(angle / MAX_ANGLE);
      auto origin = orientation.get().get();
      for (auto a = step_angle; angle - a > 0.00001; a += step_angle) {
        orientation.get().push(make_rotation(axis, a) * origin);
      }
      orientation.get().push(make_rotation(axis, angle) * origin);
    }
  }

  static void animate(Canvas& canvas, Orientation& orientation, const Vector& axis, float_t angle, EasingFn easing_fn) {
    Rotation<W> r(orientation, axis, angle, 1, easing_fn, false);
    r.step(canvas);
  }

private:

  static constexpr float_t MAX_ANGLE = 2 * PI_F / W;
  std::reference_wrapper<Orientation> orientation;
  Vector axis;
  float_t total_angle;
  EasingFn easing_fn;
  float_t from;
  float_t to;
};

template<int W>
class RandomWalk : public Animation<RandomWalk<W>> {
public:
  RandomWalk(Orientation& orientation, const Vector& v_start) :
    Animation<RandomWalk<W>>(-1, false),
    orientation(orientation),
    v(Vector(v_start).normalize())
  {
    Vector u = X_AXIS;
    if (std::abs(dot(v, u)) > 0.99) {
      u = Y_AXIS;
    }
    direction = cross(v, u).normalize();
    noiseGenerator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noiseGenerator.SetFrequency(NOISE_SCALE);
    noiseGenerator.SetSeed(hs::rand_int(0, 65535));
  }

  void step(Canvas& canvas) override {
    Animation<RandomWalk<W>>::step(canvas);
    float_t pivotAngle = noiseGenerator.GetNoise(this->t * NOISE_SCALE, 0.0f) * PIVOT_STRENGTH;
    direction = rotate(direction, make_rotation(v, pivotAngle)).normalize();
    Vector walk_axis = cross(v, direction).normalize();
    v = rotate(v, make_rotation(walk_axis, WALK_SPEED)).normalize();
    direction = rotate(direction, make_rotation(walk_axis, WALK_SPEED)).normalize();
    Rotation<W>::animate(canvas, orientation, walk_axis, WALK_SPEED, ease_mid);
  }

private:

  static constexpr float_t WALK_SPEED = 0.12f;
  static constexpr float_t PIVOT_STRENGTH = 1.5f;
  static constexpr float_t NOISE_SCALE = 0.05f;

  FastNoiseLite noiseGenerator;
  std::reference_wrapper<Orientation> orientation;
  Vector v;
  Vector direction;
};

///////////////////////////////////////////////////////////////////////////////

using AnimationVariant = std::variant<
  NullAnimation,
  Sprite,
  Transition,
  Mutation,
  RandomTimer,
  PeriodicTimer,
  Rotation<MAX_W>,
  Motion<MAX_W>,
  RandomWalk<MAX_W>
>;

struct TimelineEvent {
  int start;
  AnimationVariant animation;
};

class Timeline {
public:

  Timeline() :
    num_events(0)
  {}

  template <typename A>
  Timeline& add(float_t in_frames, A animation) {
    if (num_events >= MAX_EVENTS) {
      Serial.println("Timeline full, failed to add animation!");
      return *this;
    }
    TimelineEvent& e = events[num_events++];
    e.start = t + in_frames;
    e.animation = std::move(animation);
    return *this;
  }

  void step(Canvas& canvas) {
    t++;
    for (int i = 0; i < num_events; ++i) {
      if (t < events[i].start) {
        continue;
      }
      std::visit([&](auto& a) { a.step(canvas); }, events[i].animation);
      bool done = std::visit([](auto& a) { return a.done(); }, events[i].animation);
      if (done) {
        std::visit([](auto& a) { a.post_callback(); }, events[i].animation);
        num_events--;
        if (i < num_events) {
          events[i] = std::move(events[num_events]);
          i--;
        }
      }
    }
  }

  int t = 0;

private:

  static constexpr int MAX_EVENTS = 32;
  std::array<TimelineEvent, MAX_EVENTS> events;
  int num_events;
};
