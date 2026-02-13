/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include <vector>
#include "../effects_engine.h"

template <int W, int H>
class SpinShapes : public Effect {
public:
    SpinShapes() : Effect(W, H), filters(Filter::World::Orient<W>(camera), Filter::Screen::AntiAlias<W, H>()) {
        rebuild();
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        
        // Step animations (manual management to avoid Timeline limits)
        for (auto& rot : rotations) {
            rot.step(canvas);
        }
        
        // Draw shapes
        for (auto& shape : shapes) {
            draw_shape(canvas, shape);
        }
    }

private:
    struct Shape {
        Vector normal;
        Orientation<W> orientation;
        int layer;
    };
    
    // Parameters
    int sides = 3;
    float radius = 0.2f;
    int count = 40;
    
    std::vector<Shape> shapes;
    std::vector<Animation::Rotation<W>> rotations;
    
    Orientation<W> camera;
    Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>> filters;

    void rebuild() {
        shapes.clear();
        rotations.clear();
        
        // Ensure stability for references
        shapes.reserve(count * 2); 
        rotations.reserve(count * 2);

        for (int i = 0; i < count; i++) {
            Vector normal = fib_spiral(count, 0, i);

            // Layer 1
            {
                shapes.push_back({normal, Orientation<W>(), 0});
                Shape& s = shapes.back();
                // Rotation(orientation, axis, angle, duration, easing, repeat, space)
                rotations.emplace_back(s.orientation, s.normal, 2 * PI_F, 300 + i * 2, ease_mid, true, Animation::Space::Local);
            }

            // Layer 2
            {
                shapes.push_back({normal, Orientation<W>(), 1});
                Shape& s = shapes.back();
                rotations.emplace_back(s.orientation, s.normal, -2 * PI_F, 300 + i * 2, ease_mid, true, Animation::Space::Local);
            }
        }
    }

    void draw_shape(Canvas& canvas, Shape& shape) {
        float t = (shape.normal.j + 1.0f) / 2.0f;
        auto c = Palettes::richSunset.get(t);

        auto fragment_shader = [&](const Vector& p, const Fragment& f) -> Fragment {
            Fragment out = f;
            out.color = Color4(c, 0.6f); // alpha fixed
            return out;
        };

        Basis basis = make_basis(shape.orientation.get(), shape.normal);
        float phase = (shape.layer == 0) ? 0.0f : PI_F / sides;
        
        // Scan::SphericalPolygon::draw(pipeline, canvas, basis, radius, sides, fragment_shader, phase)
        Scan::SphericalPolygon::draw<W, H>(filters, canvas, basis, radius, sides, fragment_shader, phase);
    }
};
