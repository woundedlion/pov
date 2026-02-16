/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

template <int W, int H>
class FlamingMesh : public Effect {
public:
    FlamingMesh() : Effect(W, H),
        orientation(),
        noise(),
        filters(
            Filter::World::Orient<W>(orientation),
            Filter::Screen::Temporal<W, 100000>([this](float x, float y) { return calculate_delay(x, y); }, 2.0f),
            Filter::Screen::AntiAlias<W, H>()
        ),
        palette(Palettes::richSunset)
    {
        persist_pixels = false;
        
        // Initialize with JS defaults
        params.temporalEnabled = true;
        params.windowSize = 8;
        params.delayBase = 10;
        params.delayAmp = 10;
        params.speed = 1.0f;
        params.noiseFreq = 0.125f;

        // Initialize noise
        noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
        noise.SetFrequency(params.noiseFreq);

        // Initialize mesh (Dodecahedron matches Solids.get(3) or 'dodecahedron')
        mesh = Solids::Platonic::dodecahedron();

        timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS));
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        noise.SetFrequency(params.noiseFreq);
        filters.next.set_window_size(params.windowSize);

        Plot::Mesh::draw<W, H>(filters, canvas, mesh, [&](const Vector& v, Fragment& f) {
            float t_val = (v.j + 1.0f) * 0.5f;
            f.color = palette.get(t_val);
        });

        filters.flush(canvas, [](float x, float y, float t) { return Color4(0, 0, 0, 0); }, 1.0f);
        
        t += 1.0f;
    }

    // Parameters matching JS
    struct Params {
        bool temporalEnabled;
        int windowSize;
        float delayBase;
        float delayAmp;
        float speed;
        float noiseFreq;
    } params;

private:
    float t = 0.0f;

    Orientation<W> orientation;
    Timeline<W> timeline;
    FastNoiseLite noise;
    PolyMesh mesh;
    ProceduralPalette palette;

    // Temporal filter delay calculation 
    float calculate_delay(float x, float y) {
        if (!params.temporalEnabled) return 0.0f;
        float noiseVal = noise.GetNoise(x, y, t * params.speed);
        return std::max(0.0f, params.delayBase + (noiseVal * 0.5f + 0.5f) * params.delayAmp);
    }

    Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::Temporal<W, 100000>, Filter::Screen::AntiAlias<W, H>> filters;
};
