/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include <map>
#include <stack>
#include <array>

struct SphericalTurtle {
  Vector pos;
  Vector heading;

  SphericalTurtle(Vector p, Vector h)
      : pos(p.normalize()), heading(h.normalize()) {}

  // Move forward by dist radians
  std::pair<Vector, Vector> forward(float dist) {
    Vector axis = cross(pos, heading).normalize();
    Quaternion rot = make_rotation(axis, dist);

    Vector start = pos;
    pos = rotate(pos, rot);
    heading = rotate(heading, rot);

    return {start, pos};
  }

  // Turn by angle radians around local normal (pos)
  void turn(float angle) {
    Quaternion rot = make_rotation(pos, angle);
    heading = rotate(heading, rot);
  }
};

template <int W, int H> class LSystem : public Effect {
public:
  struct Rule {
    char init;
    const char *result;
  };

  struct Ruleset {
    const char *name;
    const char *axiom;
    std::array<Rule, 4> rules;
    int num_rules;
    float angle; // Degrees
    float step;
    int iterations;
  };

  FLASHMEM LSystem()
      : Effect(W, H), filters(Filter::World::Orient<W>(orientation),
                              Filter::Screen::AntiAlias<W, H>()) {
    registerParam("Rule", &params.rule_idx, 0.0f, 2.0f);
    registerParam("Angle", &params.angle_mod, -15.0f, 15.0f);
    registerParam("Step", &params.step_mod, 0.5f, 2.0f);

    setup_rules();
    set_ruleset(0); // Tree
    timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 2400,
                                           ease_mid, true));
  }

  bool show_bg() const override { return false; }

  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    // Check for ruleset change
    static int last_rule = -1;
    if ((int)params.rule_idx != last_rule) {
      set_ruleset((int)params.rule_idx);
      last_rule = (int)params.rule_idx;
    }

    // Check for live parameter updates (requires regeneration if changed)
    // Optimization: only regenerate if changed? For now, we regenerate or check
    // diff But regenerate() is heavy? L-System string gen is fast for small N.
    // Geometry gen is okay. Let's check diffs.
    static float last_angle = 0;
    static float last_step = 0;
    if (std::abs(params.angle_mod - last_angle) > 0.01f ||
        std::abs(params.step_mod - last_step) > 0.01f) {
      regenerate();
      last_angle = params.angle_mod;
      last_step = params.step_mod;
    }

    for (const auto &seg : segments) {
      // Color by y position
      auto fragment_shader = [&](const Vector &v, Fragment &f) {
        f.color = palette.get((v.j + 1.0f) * 0.5f);
      };
      Plot::Line::draw<W, H>(filters, canvas, Fragment(seg.first),
                             Fragment(seg.second), fragment_shader);
    }
  }

  // Params
  void set_ruleset(int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= rulesets.size())
      return;
    current_ruleset = rulesets[idx]; // Copy
    regenerate();
  }

  void regenerate() {
    char s[2048];
    strncpy(s, current_ruleset.axiom, sizeof(s));
    s[sizeof(s) - 1] = '\0';

    char next_s[2048];

    for (int i = 0; i < current_ruleset.iterations; ++i) {
      next_s[0] = '\0';
      int next_len = 0;
      for (int j = 0; s[j] != '\0'; ++j) {
        char c = s[j];
        bool found = false;
        for (int r = 0; r < current_ruleset.num_rules; ++r) {
          if (current_ruleset.rules[r].init == c) {
            int r_len = strlen(current_ruleset.rules[r].result);
            if (next_len + r_len < sizeof(next_s) - 1) {
              strcpy(&next_s[next_len], current_ruleset.rules[r].result);
              next_len += r_len;
            }
            found = true;
            break;
          }
        }
        if (!found) {
          if (next_len < sizeof(next_s) - 1) {
            next_s[next_len++] = c;
            next_s[next_len] = '\0';
          }
        }
      }
      strcpy(s, next_s);
    }

    float step = current_ruleset.step * params.step_mod;
    float angle = (current_ruleset.angle + params.angle_mod) * PI_F / 180.0f;

    Vector pos(0, -1, 0);
    Vector heading(0, 0, 1);

    SphericalTurtle turtle(pos, heading);
    std::stack<SphericalTurtle> stack;
    segments.clear();

    for (int i = 0; s[i] != '\0'; ++i) {
      char c = s[i];
      if (c == 'F') {
        segments.push_back(turtle.forward(step));
      } else if (c == '+') {
        turtle.turn(angle);
      } else if (c == '-') {
        turtle.turn(-angle);
      } else if (c == '[') {
        stack.push(turtle);
      } else if (c == ']') {
        if (!stack.empty()) {
          turtle = stack.top();
          stack.pop();
        }
      }
    }
  }

private:
  Orientation<W> orientation;
  Timeline<W> timeline;
  Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>>
      filters;
  ProceduralPalette palette = Palettes::richSunset;

  std::vector<Ruleset> rulesets;
  Ruleset current_ruleset;
  std::vector<std::pair<Vector, Vector>> segments;

  void setup_rules() {
    rulesets.push_back({"Tree",
                        "X",
                        {{{'X', "F[+X][-X]FX"}, {'F', "FF"}}},
                        2,
                        35.0f,
                        0.25f,
                        4});

    rulesets.push_back(
        {"Bush", "F", {{{'F', "FF-[-F+F+F]+[+F-F-F]"}}}, 1, 25.0f, 0.1f, 4});

    rulesets.push_back(
        {"Mosaic", "F++F++F", {{{'F', "F-F++F-F"}}}, 1, 79.0f, 0.33f, 3});
  }

  struct Params {
    float rule_idx = 0.0f;
    float angle_mod = 0.0f;
    float step_mod = 1.0f;
  } params;
};
