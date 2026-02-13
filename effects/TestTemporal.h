/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include <vector>
#include <functional>
#include "../effects_engine.h"

template <int W, int H>
class TestTemporal : public Effect {
public:
    enum DelayMode {
        VerticalWave,
        DiagonalSpiral,
        LiquidTime,
        QuantumTunnel,
        Datamosh
    };

    TestTemporal() : Effect(W, H),
        orientation(),
        noise(),
        filters(
             Filter::World::Orient<W>(orientation),
             Filter::Screen::Temporal<W, 200>([this](float x, float y) { return calculate_delay(x, y); }),
             Filter::Screen::AntiAlias<W, H>()
        ),
        circular_source(Palettes::richSunset),
        palette(circular_source),
        modifier(0.02f)
    {
          this->persist_pixels = false;

        noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    
        timeline.add(0, Animation::RandomWalk<W>(orientation, Vector(0, 1, 0)));
        palette.add(&modifier);
        timeline.add(0, Animation::PaletteAnimation(modifier));

        rebuild_mesh();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        t += 0.01f; // internal time
        
        
        noise.SetFrequency(liquid_params.noiseFreq);
        
        Plot::Mesh::draw<W, H>(filters, canvas, mesh, [&](const Vector& v, const Fragment& f) -> Fragment {
             float t = (v.j + 1.0f) * 0.5f;
             Fragment out = f;
             out.color = palette.get(t);
             return out;
        });
        
        filters.flush(canvas, [](float x, float y, float t) { return Color4(0,0,0,0); }, 1.0f);
    }
    
    // Params
    struct GlobalParams {
        bool temporalEnabled = true;
        int windowSize = 2; // Unused in C++ buffer logic? 
        float speed = 0.01f;
    } global;
    
    struct VerticalWaveParams {
        float delayBase = 8.0f;
        float delayAmp = 4.0f;
        float frequency = 0.3f;
    } vertical_params;
    
    struct DiagonalSpiralParams {
        float delayBase = 8.0f;
        float delayAmp = 4.0f;
        float xSpirals = 2.0f;
        float yFreq = 0.3f;
    } diagonal_params;

    struct LiquidTimeParams {
        float delayBase = 8.0f;
        float delayAmp = 4.0f;
        float noiseFreq = 0.03f;
        float timeScale = 2.0f;
    } liquid_params;

    struct QuantumTunnelParams {
        float delayBase = 8.0f;
        float delayAmp = 4.0f;
        float spiralTightness = 10.0f;
        float spiralAngle = 5.0f;
    } quantum_params;

    struct DatamoshParams {
        float delayBase = 8.0f;
        float delayAmp = 4.0f;
        float flowSpeed = 0.1f;
        float glitchScale = 15.0f;
    } datamosh_params;


private:
    DelayMode current_mode = VerticalWave;
    float t = 0;

    CircularPalette circular_source;
    AnimatedPalette palette;
    CycleModifier modifier;
  
    Orientation<W> orientation;
    Timeline<W> timeline;
    FastNoiseLite noise;
    Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::Temporal<W, 200>, Filter::Screen::AntiAlias<W, H>> filters;
    PolyMesh mesh;
    

    void rebuild_mesh() {
        mesh = Solids::Archimedean::truncatedIcosidodecahedron();
    }
    
    float calculate_delay(float x, float y) {
        if (!global.temporalEnabled) return 0.0f;
        
        switch (current_mode) {
            case VerticalWave: {
                float phase = y * vertical_params.frequency + t * global.speed * 10.0f; 
                return std::max(0.0f, vertical_params.delayBase + sinf(phase) * vertical_params.delayAmp);
            }
            case DiagonalSpiral: {
                float xPhase = (x / W) * PI_F * 2.0f * diagonal_params.xSpirals;
                float yPhase = y * diagonal_params.yFreq;
                float phase = xPhase + yPhase + t * global.speed * 10.0f;
                return std::max(0.0f, diagonal_params.delayBase + sinf(phase) * diagonal_params.delayAmp);
            }
            case LiquidTime: {
                // 3D Noise
                float noiseVal = noise.GetNoise(x, y, t * liquid_params.timeScale * 100.0f);
                return liquid_params.delayBase + (noiseVal + 1.0f) * 0.5f * liquid_params.delayAmp;
            }
            case QuantumTunnel: {
                // Polar
                float u = (x / W) * 2.0f - 1.0f;
                float v = (y / H) * 2.0f - 1.0f;
                float radius = sqrtf(u*u + v*v);
                float angle = atan2f(v, u);
                float spiral = sinf(radius * quantum_params.spiralTightness - angle * quantum_params.spiralAngle + t * global.speed * 10.0f);
                return std::max(0.0f, quantum_params.delayBase + (spiral + 1.0f) * quantum_params.delayAmp);
            }
            case Datamosh: {
                float flow = sinf(y * datamosh_params.flowSpeed + t * 0.05f * 10.0f);
                int blockSize = 8;
                int column = (int)(x / blockSize);
                float glitchOffset = sinf(column * 12.9898f) * datamosh_params.glitchScale;
                float total = flow * 10.0f + glitchOffset;
                return std::max(0.0f, fmodf(datamosh_params.delayBase + std::abs(total), datamosh_params.delayAmp));
            }
        }
        return 0.0f;
    }
};
