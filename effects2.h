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

typedef ()

class Thruster {
public:

  Thruster(DrawFn drawFn, const Orientation& orientation, const Vector& thrustPoint) :
    drawFn(drawFn),
    orientation(orientation),
    thrustPoint(thrustPoint),
    exhaustRadius(0),
    exhaustMotion(exhaustRadius, 0.3, 8, easeMid),
    exhaustSprite(drawFn, 16, 0, easeMid, 16, easeOutExpo)
  {
  }

  done() {
    return exhaustMotion.done()
      && exhaustSprite.done();
    ;
  }

  step() {
    exhaustSprite.step();
    exhaustMotion.step();
  }

private:

  DrawFn drawFn;
  Orientation orientation;
  Vector thrustPoint;
  double exhaustRadius;
  Transition exhaustMotion;
  Sprite exhaustSprite;
}

template <int W>
class Thruster : public Effect {
public:
  Thrusters() :
    Effect(W),
    palette(
      [0.5, 0.5, 0.5],
      [0.5, 0.5, 0.5],
      [0.3, 0.3, 0.3],
      [0.0, 0.2, 0.6]),

  {
    ringOutput.chain(new FilterAntiAlias());
    timeline.add(0,
      new Sprite(this.drawRing.bind(this), -1,
        16, easeInSin,
        16, easeOutSin)
    );
  }

  bool show_bg() const { return false; }

  void draw_frame() {
  }

private:

  ProceduralPalette palette;
  FilterRaw ringOutput;
  int t = 0;
  Vector ring;
  Orientation orientation;
  Orientation to;
  std::vector<Thruster> thrusters;
  double amplitude = 0;
  double warp_phase = 0;
  double radius = 1;

  Timeline timeline;

}