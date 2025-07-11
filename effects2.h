#pragma once

#include "effects_infra.h"
#include <functional>
#include <memory>
#include <vector>
#include <map>

DEFINE_GRADIENT_PALETTE(lucid_dream_p) {
  0, 6, 4, 47,  
  128, 162, 84, 84,
  255, 252, 114, 0 
};

DEFINE_GRADIENT_PALETTE(red_orange_p) {
  0, 255, 96, 0,
  255, 255, 0, 0 
};

DEFINE_GRADIENT_PALETTE(blue_purple_p) {
  0, 0, 0, 255,
  255, 135, 0, 135
};

/*

struct EventHandlers {
  std::function<void()> enter;
  std::function<void()> draw;
  std::function<void()> animate;
  std::function<void()> exit;
};

template <int W>
class PolyRot : public Effect {
public:
  PolyRot() :
    Effect(W),
    poly_mask_mask(4),
  {
    poly_mask.chain(poly_mask_mask);
    transition();
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    states[state].draw();
    states[state].animate();
  }

private:

  enum State {
    INIT,
    GEN_POLY,
    SPIN_RING,
    SPIN_POLY,
    SPLIT_POLY
  };

  int t = 0;
  int state_index = -1;
  State state = INIT;

  Dodecahedron poly;
  Orientation poly_top_orientation;
  Orientation poly_bottom_orientation;

  Vector ring = { 0, 0, 1 };
  Orientation ring_orientation;

  Vector spin_axis = { 0, 0, 1 };
  Orientation spin_axis_orientation;

  std::unique_ptr<Path> path;
  std::unique_ptr<Motion> motion;
  std::unique_ptr<Rotation> rotation;
  std::unique_ptr<Rotation> rotation_rev;

  FilterAntiAlias<W, H> poly_mask;
  FilterDecayMask<W, H> poly_mask_mask;
  FilterAntiAlias<W, H> output;

  void transition() {
    state_index = (state_index + 1) % sequence.size();
    transitionTo(sequence[state_index]);
  }

  void transitionTo(State next) {
    Serial.printf("State transition to %d\n", next);
    if (state != INIT) {
      states[state].exit();
    }
    t = 0;
    state = next;
    states[state].enter();
  }

  std::map<State, EventHandlers> states = {
    {GEN_POLY, {
      [this]() { // enter
        Serial.printf("GEN_POLY Enter\n");
        double duration = 160;
        poly = Dodecahedron();
        path.reset(new Path());
        path->append_segment(
          [this](double t) {
            return Vector(ring_orientation.orient(lissajous(10, 0.5, 0, t)));
          },
          2 * PI,
          duration,
          ease_in_out_sin
        );
        motion.reset(new Motion(*path, duration));
      },
      [this]() { // draw
        Serial.printf("GEN_POLY Draw\n");
        poly_mask_mask.decay();
        Pixels pixels;
        Dots dots;
        Canvas canvas(*this);

        for (int i = 0; i < W * H; ++i) {
          canvas(i) = CHSV(0, 0, 0);
        }

        // Draw ring trail into polygon mask
        auto n = ring_orientation.length();
        for (int i = 0; i < n; ++i) {
          auto normal = ring_orientation.orient(ring, i);
          plot_dots(
            poly_mask,
            draw_ring(normal, 1, [](auto& v) { return CHSV(0, 0, 0); }),
            n == 1 ? 0 : (n - 1 - i) / (n - 1));
        }
        ring_orientation.collapse();
        Vector normal = ring_orientation.orient(ring);

        // Draw polyhedron
        auto vertices = poly_top_orientation.orient(poly.vertices);
        dots = draw_polyhedron(
          vertices,
          poly.euler_path,
          [this, &normal](auto& v) { return distance_gradient(v, normal, red_orange_p, blue_purple_p); }
        );
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), poly_mask_mask.mask(xy, p.color));
        }

        // Draw current ring position
        dots = draw_ring(normal, 1, [](auto& v) { return CHSV(HUE_BLUE, 0, 64); });
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }
    },
      [this]() { // animate
        Serial.printf("GEN_POLY Animate\n");
        if (motion->done()) {
          transition();
        } else {
          motion->move(ring_orientation);
        }
      },
      [this]() { // exit

      },
    }},
    {SPIN_POLY, {
      [this]() { // enter
        Serial.printf("SPIN_POLY Enter\n");
        double duration = 192;
        poly = Dodecahedron();
        spin_axis = spin_axis_orientation.orient(spin_axis);
        rotation.reset(new Rotation(spin_axis, 4 * PI, duration));
        path.reset(new Path());
        path->append_segment(
          [this](double t) {
            return Vector(lissajous(12.8, 2 * PI, 0, t));
          },
          1,
          duration,
          ease_mid
        );
        motion.reset(new Motion(*path, duration));
      },
      [this]() { // draw
        Serial.printf("SPIN_POLY Draw\n");
        Canvas canvas(*this);
        Dots dots;
        Pixels pixels;

        for (int i = 0; i < W * H; ++i) {
          canvas(i) = CHSV(0, 0, 0);
        }

        // Draw polyhedron
        auto normal = ring_orientation.orient(ring);
        auto vertices = poly_top_orientation.orient(poly.vertices);
        dots = draw_polyhedron(
          vertices,
          poly.euler_path,
          [this, &normal](auto& v) { return distance_gradient(v, normal, red_orange_p, blue_purple_p); }
        );
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }

        // Draw Ring
        dots = draw_ring(normal, 1, [](auto& v) { return CHSV(HUE_BLUE, 0, 128); });
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }
      },
      [this]() { // animate
        Serial.printf("SPIN_POLY Animate\n");
        if (rotation->done()) {
          transition();
        } else {
          motion->move(spin_axis_orientation);
          rotation->set_axis(spin_axis_orientation.orient(spin_axis))
            .rotate(poly_top_orientation, ease_in_out_sin);
        }
      },
      [this]() { // exit

      },
    }},
    {SPLIT_POLY, {
      [this]() { // enter
        Serial.printf("SPLIT_POLY Enter\n");
        double duration = 96;
        poly = Dodecahedron();
        auto normal = ring_orientation.orient(ring);
        bisect(poly, poly_top_orientation, normal);
        poly_bottom_orientation.set(poly_top_orientation.get());
        rotation.reset(new Rotation(normal, 4 * PI, duration));
        rotation_rev.reset(new Rotation(normal.inverse(), 4 * PI, duration));
      },
      [this]() { // draw
        Serial.printf("SPLIT_POLY Draw\n");
        Dots dots;
        Pixels pixels;
        Canvas canvas(*this);
        for (int i = 0; i < W * H; ++i) {
          canvas(i) = CHSV(0, 0, 0);
        }

        auto normal = ring_orientation.orient(ring);
        VertexList vertices;
        std::transform(poly.vertices.begin(), poly.vertices.end(), std::back_inserter(vertices),
          [this, &normal](auto& v) {
            auto u = poly_top_orientation.orient(v);
            if (is_over(u, normal)) {
              return u;
            } else {
              return poly_bottom_orientation.orient(v);
            }
          });
        
        // Draw polyhedron
        dots = draw_polyhedron(
          vertices,
          poly.euler_path,
          [this, &normal](auto& v) { return distance_gradient(v, normal, red_orange_p, blue_purple_p); }
        );
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }

        // Draw Ring
        dots = draw_ring(normal, 1, [](auto& v) { return CHSV(HUE_BLUE, 0, 128); });
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }

      },
      [this]() { // animate
        Serial.printf("SPLIT_POLY Animate\n");
        if (rotation->done()) {
          transition();
        } else {
          rotation->rotate(poly_top_orientation, ease_in_out_sin);
          rotation_rev->rotate(poly_bottom_orientation, ease_in_out_sin);
        }
      },
      [this]() { // exit

      },
    }},
    {SPIN_RING, {
      [this]() { // enter
        srand(time(NULL));
        Serial.printf("SPIN_RING Enter\n");
        double duration = 16;
        poly = Dodecahedron();
        Vector from_normal = ring_orientation.orient(ring);
        Vector to_normal = ring_orientation.orient(poly.vertices[rand() % poly.vertices.size()]);
        path.reset(new Path());
        path->append_line(from_normal, to_normal, true);
        motion.reset(new Motion(*path, duration));

      },
      [this]() { // draw
        Serial.printf("SPIN_RING Draw\n");
        Canvas canvas(*this);
        for (int i = 0; i < W * H; ++i) {
          canvas(i) = CHSV(0, 0, 0);
        }

        Dots dots;
        Pixels pixels;

        // Draw polyhedron
        auto normal = ring_orientation.orient(ring);
        auto vertices = poly_top_orientation.orient(poly.vertices);
        dots = draw_polyhedron(
          vertices,
          poly.euler_path,
          [this, &normal](auto& v) { return distance_gradient(v, normal, red_orange_p, blue_purple_p); }
        );
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }
        
        // Draw Ring
        dots = draw_ring(normal, 1, [](auto& v) { return CHSV(HUE_BLUE, 0, 128); });
        pixels = plot_dots(output, dots);
        for (auto& [xy, p] : pixels) {
          canvas(xy) = blend_over(canvas(xy), p.color);
        }
      },
      [this]() { // animate
        Serial.printf("SPIN_RING Animate\n");
        if (motion->done()) {
          transition();
        } else {
          motion->move(ring_orientation);
        }
      },
      [this]() { // exit

      },
    }},
  };

  std::vector<State> sequence = {
   GEN_POLY,
    SPIN_RING,
    SPIN_POLY,
    SPIN_RING,
    SPLIT_POLY,
    SPIN_RING,
  };

};
*/

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
      16, 0, ease_mid, 16, ease_out_expo)
  {
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
  }

  bool show_bg() const { return false; }

  void draw_thruster(const Orientation& orientation, const Vector& thrust_point, double radius, double opacity) {
    auto dots = ::draw_ring<W>(orientation, thrust_point, radius,
      [opacity](auto&, auto) { return CHSV(255, 255, 255 * opacity); });
    plot_dots(dots, filters, pixels);
  }

  void on_fire_thruster() {
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
      [](double t, double m) { return 0.7 * exp(-2 * t); },
      32, 
      ease_mid));
    timeline.add(1 / 16,
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
  }

  double ring_fn(double t) {
    return sin_wave(-1, 1, 2, warp_phase)(t) // ring
      * sin_wave(-1, 1, 3, 0)((this->timeline.t % 32) / 32) // oscillation
      * amplitude;
  }

  void draw_ring(double opacity) {
    rotate_between<W>(orientation, to);
    orientation.collapse();
    to.collapse();
    auto dots = draw_fn<W>(orientation, ring, radius,
      [this](auto t) -> auto { 
        return this->ring_fn(t); 
      },
      [this, opacity](auto& v, auto t) -> auto {
        auto z = orientation.orient(X_AXIS);
        auto color = palette.get(angle_between(z, v) / PI);
        color.v = dim8_lin(color.v * opacity);
        return color;
      }
    );
    plot_dots(dots, filters, pixels);
  }

  void draw_frame() {
      pixels.clear();
      timeline.step();

      Canvas canvas(*this);
      for (auto& [xy, p] : pixels) {
        canvas(xy) = blend_over(canvas(xy), p.color);
      }
  }

private:

  ProceduralPalette palette;
  FilterAntiAlias<W, H> filters;
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