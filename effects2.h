#pragma once

#include "effects_infra.h"
#include <functional>
#include <memory>
#include <vector>
#include <map>



/*
template <int>
class Thrusters;

template <int W>
class Thruster : public Animation {
public:

  Thruster(Thrusters<W>& parent, const Orientation& orientation, const Vector& thrust_point) :
    Animation(0, false),
    orientation(orientation),
    thrust_point(thrust_point),
    exhaust_radius(0),
    exhaust_motion(exhaust_radius, 0.3, 8, ease_mid),
    exhaust_sprite([this, &parent](auto opacity) {
        parent.draw_thruster(this->orientation, this->thrust_point, this->exhaust_radius, opacity); 
      },
      32, 0, ease_mid, 32, ease_out_expo)
  {
    Serial.println("Thruster() Create");
  }

  ~Thruster() {
    Serial.println("Thruster() Destroy");
  }

  bool done() {
    return exhaust_motion.done()
      && exhaust_sprite.done();
    ;
  }

  void step() {
    exhaust_sprite.step();
    exhaust_motion.step();
  }

private:

  Orientation orientation;
  Vector thrust_point;
  double exhaust_radius;
  Transition exhaust_motion;
  Sprite exhaust_sprite;
};

template <int W>
class Thrusters : public Effect {
public:
  Thrusters() :
    Effect(W),
    palette(
      { 0.5, 0.5, 0.5 },
      { 0.5, 0.5, 0.5 },
      { 0.3, 0.3, 0.3 },
      { 0.0, 0.2, 0.6 }),
    ring(Vector(0.5, 0.5, 0.5).normalize())
  {
    timeline
      .add(0,
        std::make_shared<Sprite>(
          [this](auto o) { return draw_ring(o); },
          -1,
          16, ease_in_sin,
          16, ease_out_sin)
      )
      
      .add(0,
        std::make_shared<RandomTimer>(
          8, 48,
          [this]() {  this->on_fire_thruster(); },
          true)
      );
    Serial.println("Thrusters()");
  }

  void draw_thruster(const Orientation& orientation, const Vector& thrust_point, double radius, double opacity) {
    auto dots = ::draw_ring<W>(orientation, thrust_point, radius,
      [opacity](auto&, auto) { return CHSV(0, 0, 255 * opacity); });
    plot_dots<W>(dots, filters, pixels);
  }

  void on_fire_thruster() {
    Serial.println("on_fire_thruster()");
    auto warp_phase = hs::rand_dbl() * 2 * PI;
    auto thrust_point = fn_point(
      [this](auto t) { return ring_fn(t); },
      ring, 1, warp_phase);
    auto thrust_opp = fn_point(
      [this](auto t) { return ring_fn(t); },
      ring, 1, 
      warp_phase + PI);

    // warp ring
    if (warp && !warp->done()) {
      warp->cancel();
    }
    
    warp.reset(new Mutation(
      amplitude, 
      [](auto t) { return 0.7 * exp(-2.0 * t); },
      32, 
      ease_mid));
    timeline.add(1.0 / 16,
      warp
    );
   
    // Spin ring
    auto thrust_axis = cross(
      orientation.orient(thrust_point),
      orientation.orient(ring))
      .normalize();
    timeline.add(0,
      std::make_shared<Rotation<W>>(
        orientation, thrust_axis, 2 * PI, 8 * 16, ease_out_expo)
    );
    
    // show thruster
    timeline.add(0,
      std::make_shared<Thruster<W>>(
        *this,
        orientation,
        thrust_point));
    timeline.add(0,
      std::make_shared<Thruster<W>>(
        *this,
        orientation,
        thrust_opp));
    Serial.println("end on_fire_thruster()");
  }

  double ring_fn(double t) {
    return sin_wave(-1, 1, 2, warp_phase)(t) // ring
      * sin_wave(-1, 1, 3, 0)((this->timeline.t % 32) / 32.0) // oscillation
      * amplitude;
  }

  void draw_ring(double opacity) {
    Serial.println("draw_ring");
    rotate_between<W>(orientation, to);
    orientation.collapse();
    to.collapse();
    auto dots = draw_fn<W>(orientation, ring, radius,
      [this](auto t) -> auto { 
        return this->ring_fn(t); 
      },
      [this, opacity](auto& v, auto t) -> auto {
        auto z = orientation.orient(X_AXIS);
        auto color = rgb2hsv_approximate((palette.get(angle_between(z, v) / PI)));
        color.v = dim8_lin(color.v * opacity);
        return color;
      }
    );
    plot_dots<W>(dots, filters, pixels);
    Serial.println("end draw_ring");
  }

  void draw_frame() {
    pixels.clear();
    timeline.step();
    Canvas canvas(*this);
    canvas.clear_buffer();
    for (auto& [xy, p] : pixels) {
      canvas(xy) = p;
    }
  }

  bool show_bg() const { return false; }

private:

  ProceduralPalette palette;
  FilterAntiAlias<W> filters;
  Vector ring;
  Orientation orientation;
  Orientation to;
  std::vector<Thruster<W>> thrusters;
  double amplitude = 0;
  double warp_phase = 0;
  double radius = 1;
  std::shared_ptr<Mutation> warp;
  Timeline timeline;
  Pixels pixels;
};

///////////////////////////////////////////////////////////////////////////////

template <int W>
class Wormhole : public Effect {
public:

  Wormhole() :
    Effect(W),
    palette(
      { 0.5, 0.5, 0.5 },
      { 1.0, 0.2, 0.5 },
      { 0.5, 0.5, 0.5 },
      { 0.3, 0.5, 0.0 }),
    normal(Z_AXIS)
  {
    timeline.add(0,
      std::make_shared<Sprite>([this](double opacity) { this->draw_rings(opacity); },
        -1, 8, ease_mid, 0, ease_mid)
    );

    // T1: Spin everything
    on_thrust_rings(1);
    on_spin_rings(1);
    on_mutate_duty_cyle(1);
    on_mutate_twist(1);
  }

  void on_mutate_duty_cyle(double in_secs = 0) {
    return;
    timeline.add(in_secs,
      std::make_shared<Mutation>(duty_cycle, sin_wave((2 * PI) / W, (8 * 2 * PI) / W, 1, PI / 2),
        160, ease_mid, true)
    );
  }

  void on_mutate_twist(double in_secs = 0) {
    return;
    timeline.add(in_secs,
      std::make_shared<Mutation>(twist, sin_wave(3.0 / W, 10.0 / W, 1, PI / 2),
        64, ease_mid, true)
    );
  }

  void on_thrust_rings(double in_secs = 0) {
    timeline.add(in_secs,
      std::make_shared<Rotation<W>>(
        orientation,
        ring_point(orientation.orient(normal), 1, hs::rand_dbl() * 2 * PI),
        4 * PI,
        16 * 16, ease_in_out_sin, false)
    );

    timeline.add(in_secs,
      std::make_shared<RandomTimer>(4 * 16, 8 * 16, [this]() { this->on_thrust_rings(); }, false)
    );
  }

  void on_spin_rings(double in_secs = 0) {
    timeline.add(in_secs,
      std::make_shared<Transition>(phase, 2 * PI, 32, ease_mid, false, true)
    );
  }

  void calc_ring_spread() {
    radii.resize(num_rings);
    for (int i = 0; i < num_rings; ++i) {
      double x = (static_cast<double>(i + 1) / (num_rings + 1)) * 2 - 1;
      double r = sqrt(pow(1 - x, 2));
      radii[i] = lerp(home_radius, r, spread_factor);
    }
  }

 void draw_rings(double opacity) {
    calc_ring_spread();
    Serial.println(phase);
    orientation.collapse();
    for (unsigned int i = 0; i < radii.size(); ++i) {
      auto dots = draw_ring<W>(orientation, normal, radii[i],
        [=](auto& v, auto t) {
          double idx = num_rings == 1 ? 0 : (1 - (static_cast<double>(i) / (num_rings - 1)));
          double darken = pow(1 - abs(radii[i] - 1), 3);
          auto color = dim(palette.get(idx), darken);
          auto r = dotted_brush(dim(color, opacity), freq,
            duty_cycle, twist, t);
          return r;
        },
        twist * i + phase);
      plot_dots(dots, filters, pixels);
    }
  }

  void draw_frame() {
    pixels.clear();
    timeline.step();
    Canvas canvas(*this);
    for (auto& [xy, p] : pixels) {
      canvas(xy) = p;
    }
  }

  bool show_bg() const { return false; }

  private:

    ProceduralPalette palette;
    FilterAntiAlias<W> filters;
    Vector normal;
    Orientation orientation;

    // TODO: int template for Mutation
    double num_rings = W;
    double spread_factor = 1;
    double home_radius = 1;
    double duty_cycle = (2 * PI) / W;
    double freq = 2;
    double twist = 10.0 / W;
    double phase = 0;
    std::vector<double> radii;

    Timeline timeline;
    Pixels pixels;
};
*/
///////////////////////////////////////////////////////////////////////////////

template <int W>
class RingSpin : public Effect {
public:

  struct Ring {
    Ring(const Vector& normal, Filter& filters, const Palette& palette, uint8_t trail_length) :
      normal(normal),
      palette(palette),
      trails(trail_length)
    {    
      trails.chain(filters);
    }

    Vector normal;
    const Palette& palette;
    Orientation orientation;
    FilterDecay<W> trails;
  };

  RingSpin() :
    Effect(W),
    alpha(0.2),
    trail_length(15)
  {
    persist_pixels = false;
    spawn_ring(X_AXIS, *palettes[0]);
    spawn_ring(Y_AXIS, *palettes[1]);
    spawn_ring(Z_AXIS, *palettes[2]);
    spawn_ring(X_AXIS, *palettes[3]);
    spawn_ring(Y_AXIS, *palettes[4]);
    spawn_ring(Z_AXIS, *palettes[5]);
  }

  bool show_bg() const { return false; }

  void spawn_ring(const Vector& normal, Palette& palette) {
    auto ring_index = rings.size();
    rings.emplace_back(normal, filters, palette, trail_length);
    timeline.add(0,
      Sprite(
        [=](Canvas& canvas, double opacity) { draw_ring(canvas, opacity, ring_index); },
        -1,
        4, easeMid,
        0, easeMid
      ));
    timeline.add(0,
      RandomWalk<W>(rings[ring_index].orientation, rings[ring_index].normal));
  }

  void draw_ring(Canvas& canvas, double opacity, size_t ring_index) {
    Dots dots;
    auto& ring = rings[ring_index];
    double s = ring.orientation.length();
    for (int i = 0; i < s; ++i) {
      dots = draw_ring(ring.orientation.orient(ring.normal, i), 1,
        [=](auto& v, auto t) { return ring.palette.get(0); });
      plot_dots(dots, ring.trails, canvas, 
       (s - 1 - i) / s,
       alpha * opacity);
    }
    ring.trails.trail(canvas,
      [=](double x, double y, double t) { return vignette(ring.palette)(1 - t); },
      alpha * opacity);
    ring.trails.decay();
    ring.orientation.collapse();
  }

  void draw_frame() {
    Canvas canvas(*this);
    timeline.step(canvas);
  }

private:

  std::vector<Ring> rings;
  double alpha;
  double trail_length;
  FilterAntiAlias<W> filters;
  std::array<const Palette*> palettes = { &richSunset, &undersea, &mangoPeel, &lemonLime, &algae, &lateSunset };
  Timeline timeline;
};

template <int W>
class Dynamo : public Effect {
public:

  struct Node {
    Node() :
      x(0), y(0), v(0)
    {}

    int x;
    int y;
    int v;
  };

  Dynamo() :
    Effect(W),
    palettes({ GenerativePalette(GradientShape::VIGNETTE, HarmonyType::ANALOGOUS) }),
    palette_normal(X_AXIS),
    speed(2),
    gap(4),
    trail_length(10),
    filters(4),
    trails(trail_length),
    orient(orientation)
  {
    persist_pixels = false;

    for (int i = 0; i < NUM_NODES; ++i) {
      nodes[i].y = i;
    }

    filters
      .chain(trails)
      .chain(orient)
      .chain(aa);

    timeline
      .add(0,
        RandomTimer(4, 64, [=](auto&) { reverse(); }, true))
      .add(0,
        RandomTimer(20, 64, [=](auto&) { color_wipe(); }, true))
      .add(0,
        RandomTimer(80, 150, [=](auto&) { rotate(); }, true));
  }

  bool show_bg() const { return false; }

  void reverse() {
    speed *= -1;
  }

  void rotate() {
    timeline.add(0,
      Rotation<W>(
        orientation, random_vector(), PI, 40, ease_in_out_sin, false));
  }

  void color_wipe() {
    palettes.push(GenerativePalette(GradientShape::VIGNETTE));
    palette_boundaries.push(0);
    timeline.add(0,
      Transition(palette_boundaries[0], PI, 20, ease_mid)
      .then([this]() {
          palette_boundaries.pop();
          palettes.pop();
        });
  }

  Pixel color(const Vector& v, double t) {
    constexpr double blend_width = PI / 8.0;
    double a = angle_between(v, palette_normal);

    for (int i = 0; i < palette_boundaries.size(); ++i) {
      double boundary = palette_boundaries[i];
      auto lower_edge = boundary - blend_width;
      auto upper_edge = boundary + blend_width;

      if (a < lower_edge) {
        return std::visit([=](auto& p) { return p.get(t); },  palettes[i]);
      }

      if (a >= lower_edge && a <= upper_edge) {
        auto blend_factor = (a - lower_edge) / (2 * blend_width);
        auto clamped_blend_factor = std::max(0, std::min(blend_factor, 1));
        auto c1 = std::visit([=](auto& p) { return p.get(t); }, palettes[i]);
        auto c2 = std::visit([=](auto& p) { return p.get(t); }, palettes[i + 1]);
        return c1.lerp16(c2, to_short(clamped_blend_factor));
      }

      auto next_boundary_lower_edge = (i + 1 < palette_boundaries.size()
        ? palette_boundaries[i + 1] - blend_width
        : PI;
      if (a > upper_edge && a < next_boundary_lower_edge) {
        return std::visit([=](auto& p) { return p.get(t); }, palettes[i + 1]);
      }
    }

    return std::visit([=](auto& p) { return p.get(t); }, palettes[0]);
  }

  void draw_frame() {
    Canvas canvas(*this);

    for (int i = std::abs(speed) - 1; i >= 0; --i) {
      pull(0);
      draw_nodes(canvas, static_cast<double>(i) / std::abs(speed));
    }
    trails.trail(
      canvas,
      [this](double x, double y, double t) { return color(pixel_to_vector(x, y), t); },
      0.5
    );

    trails.decay();
    timeline.step(canvas);
  }

private:

  double node_y(const Node& node) const {
    return (node.y / (nodes.size() - 1)) * (H - 1);
  }

  void draw_nodes(Canvas& canvas, double age) {
    Dots dots;
    for (int i = 0; i < nodes.size(); ++i) {
      if (i == 0) {
        auto from = pixel_to_vector(nodes[i].x, node_y(nodes[i]));
        auto drawn = draw_vector(from, [this](auto& v) { return color(v, 0); });
        dots.insert(dots.end(), drawn.begin(), drawn.end());
      } else {
        auto from = pixel_to_vector(nodes[i - 1].x, node_y(nodes[i - 1]));
        auto to = pixel_to_vector(nodes[i].x, node_y(nodes[i]));
        auto drawn = draw_line<W>(from, to, [this](auto& v) { return color(v, 0); });
        dots.insert(dots.end(), drawn.begin(), drawn.end());
      }
    }
    plot_dots(dots, filters, canvas, age, 0.5);
  }

  void pull(int leader) {
    nodes[leader].v = dir(speed);
    move(nodes[leader]);
    for (int i = leader - 1; i >= 0; --i) {
      drag(nodes[i + 1], nodes[i]);
    }
    for (int i = leader + 1; i < nodes.size(); ++i) {
      drag(nodes[i - 1], nodes[i]);
    }
  }

  void drag(Node& leader, Node& follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (distance(dest, leader.x, W) > gap) {
      follower.v = leader.v;
      while (distance(follower.x, leader.x, W) > gap) {
        move(follower);
      }
    } else {
      move(follower);
    }
  }
  
  void move(Node& node) {
    node.x = wrap(node.x + node.v, W);
  }

  int dir(int speed) {
    return speed < 0 ? -1 : 1;
  }

  Timeline timeline;
 
  static constexpr MAX_PALETTES = 16;
  static constexpr NUM_NODES = H;
  StaticCircularBuffer <PaletteVariant, MAX_PALETTES> palettes;
  StaticCircularBuffer <double, MAX_PALETTES - 1> palette_boundaries; 

  Vector palette_normal;
  std::array<Node, NUM_NODES> nodes;
  int speed;
  int gap;
  uint32_t trail_length;
  Orientation orientation;
  
  FilterReplicate<W> filters;
  FilterDecay<W> trails;
  FilterOrient<W> orient;
  FilterAntiAlias<W> aa;
};