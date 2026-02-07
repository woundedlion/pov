/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../plot.h"
#include "../animation.h"
#include "../filter.h"
#include "../palettes.h"
#include <vector>
#include <string>
#include <map>
#include <stack>

struct SphericalTurtle {
    Vector pos;
    Vector heading;

    SphericalTurtle(Vector p, Vector h) : pos(p.normalize()), heading(h.normalize()) {}

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

template <int W>
class LSystem : public Effect {
public:
    struct Ruleset {
        std::string name;
        std::string axiom;
        std::map<char, std::string> rules;
        float angle; // Degrees
        float step;
        int iterations;
    };

    LSystem() : Effect(W), filters(FilterOrient<W>(orientation), FilterAntiAlias<W>()) {
        setup_rules();
        set_ruleset(0); // Tree
        timeline.add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 2400, ease_mid, true));
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        for (const auto& seg : segments) {
            // Color by y position
            auto fragment_shader = [&](const Vector& v, const Fragment& f) -> Fragment {
                Fragment out = f;
                out.color = palette.get((v.j + 1.0f) * 0.5f);
                return out;
            };
            Plot::Line::draw<W>(filters, canvas, seg.first, seg.second, fragment_shader);
        }
    }
    
    // Params
    void set_ruleset(int idx) {
        if (idx < 0 || static_cast<size_t>(idx) >= rulesets.size()) return;
        current_ruleset = rulesets[idx]; // Copy
        regenerate();
    }
    
    void regenerate() {
        std::string s = current_ruleset.axiom;
        for (int i = 0; i < current_ruleset.iterations; ++i) {
             std::string next_s;
             for (char c : s) {
                 if (current_ruleset.rules.count(c)) {
                     next_s += current_ruleset.rules[c];
                 } else {
                     next_s += c;
                 }
             }
             s = next_s;
        }
        
        float step = current_ruleset.step;
        float angle = current_ruleset.angle * PI_F / 180.0f;
        
        Vector pos(0, -1, 0);
        Vector heading(0, 0, 1);
        
        SphericalTurtle turtle(pos, heading);
        std::stack<SphericalTurtle> stack;
        segments.clear();
        
        for (char c : s) {
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
    Orientation orientation;
    Timeline timeline;
    Pipeline<W, FilterOrient<W>, FilterAntiAlias<W>> filters;
    ProceduralPalette palette = Palettes::richSunset;
    
    std::vector<Ruleset> rulesets;
    Ruleset current_ruleset;
    std::vector<std::pair<Vector, Vector>> segments;
    
    void setup_rules() {
        rulesets.push_back({
            "Tree", "X", 
            {{'X', "F[+X][-X]FX"}, {'F', "FF"}},
            35.0f, 0.25f, 4
        });
        
        rulesets.push_back({
            "Bush", "F",
            {{'F', "FF-[-F+F+F]+[+F-F-F]"}},
            25.0f, 0.1f, 4
        });
        
        rulesets.push_back({
            "Mosaic", "F++F++F",
            {{'F', "F-F++F-F"}},
            79.0f, 0.33f, 3
        });
    }
};
