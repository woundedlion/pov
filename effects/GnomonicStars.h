/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"
#include <vector>

template <int W, int H>
class GnomonicStars : public Effect {
public:
    GnomonicStars() : Effect(W, H),
        filters(Filter::World::Orient<W>(orientation), Filter::Screen::AntiAlias<W, H>()),
        orientation(),
        timeline(),
        params(1,0,0,0, 0,0,1,0)
    {
        this->persist_pixels = false;
        timeline.add(0, Animation::MobiusGenerate(params, 0.5f, 0.05f));
        timeline.add(0, Animation::RandomWalk<W>(orientation, UP, Animation::RandomWalk<W>::Options::Energetic()));
    }

    bool show_bg() const override { return true; } 

    void draw_frame() override {
        Canvas canvas(*this);
        
        timeline.step(canvas);
        
        // Fragment Shader for Stars (Mango Peel Gradient based on Y)
        auto fragment_shader = [](const Vector& p, const Fragment& frag) -> Fragment {
            Fragment f_out = frag;
            float t = (p.j + 1.0f) * 0.5f; // Map [-1, 1] to [0, 1]
            Color4 c = Palettes::mangoPeel.get(t);
            f_out.color = c;
            f_out.blend = BLEND_OVER; 
            return f_out;
        };

        const int points = spiral_params.points;
        const float radius = spiral_params.star_radius;
        const int sides = spiral_params.star_sides;
        
        for (int i = 0; i < points; i++) {
            Vector v = fib_spiral(points, 0.0f, i);

            if (enable_transform) {
                v = gnomonic_mobius_transform(v, params);
            }

            // Orient the star position (using latest orientation)
            v = orientation.orient(v);
            
            // Create Basis at the star position (aligned with normal v)
            Basis basis = make_basis(Quaternion(), v); 
            
            Scan::Star::draw<W, H>(
                filters, 
                canvas, 
                basis, 
                radius, 
                sides, 
                fragment_shader
            );
        }
    }

private:
    Orientation<W> orientation;
    Timeline<W> timeline;
    Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>> filters;
    
    MobiusParams params;
    bool enable_transform = true;
    
    struct SpiralParams {
        int points = 600;
        float star_radius = 0.02f;
        int star_sides = 4;
    } spiral_params;
};
