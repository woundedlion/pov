/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../scan.h"
#include "../plot.h"
#include "../animation.h"
#include "../filter.h"
#include "../solids.h"
#include "../palettes.h"
#include "../FastNoiseLite.h"
#include <vector>
#include <functional>

template <int W>
class TestTemporal : public Effect {
public:
    enum DelayMode {
        VerticalWave,
        DiagonalSpiral,
        LiquidTime,
        QuantumTunnel,
        Datamosh
    };

    TestTemporal() : Effect(W),
        orientation(),
        noise(),
        filters(
             FilterOrient<W>(orientation),
             FilterTemporal<W, 200>([this](float x, float y) { return calculate_delay(x, y); }),
             FilterAntiAlias<W>()
        )
    {
        noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
        
        // Random Walk
        timeline.add(0, Motion<W>(orientation, random_path, 2000, true)); // Placeholder
        
        rebuild_mesh();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        t += 0.01f; // internal time
        
        noise.SetFrequency(liquid_params.noiseFreq);
        
        // Refresh proxy
        // The lambda capture handles it.
        
        Plot::Mesh::draw<W>(filters, canvas, mesh, [&](const Vector& v, const Fragment& f) -> Fragment {
             // renderMesh
             // JS: return colors.get((v.y + 1) * 0.5);
             float t = (v.j + 1.0f) * 0.5f;
             Fragment out = f;
             out.color = Palettes::richSunset.get(t);
             return out;
        });
        
        // Flush not needed for Plot (it draws immediately), but FilterTemporal processes trails in 'trail'.
        // Plot pipeline handles `trail`?
        // `Effect::draw` usually calls `pipeline.trail` if needed? 
        // No, `plot.h` draws primitives.
        // `FilterTemporal` buffer is filled during `plot`.
        // To drain/animate the buffer, we MUST call `filters.trail(...)`.
        // JS: `this.filters.flush(null, 1.0)`.
        // C++: `filters.trail(canvas, trailFn, 1.0f)`.
        // `trailFn` logic? 
        // JS `filters.flush` default uses standard decay? 
        // `FilterTemporal` C++ `trail` decrements delay and plots when <=0.
        // It needs a `trailFn` to generate color?? 
        // `FilterTemporal::plot` stores the color!
        // `FilterTemporal::trail` in C++ line 508: `trailFn` is passed but...
        // Line 514: `pass(... item.c ...)`. It uses stored color. 
        // It ignores `trailFn` output? 
        // Line 508: `trail(TrailFn auto trailFn...)`.
        // It doesn't seem to use `trailFn`.
        // So I can pass a dummy lambda.
        
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

    DelayMode current_mode = VerticalWave;
    float t = 0;

private:
    Orientation orientation;
    Timeline timeline;
    FastNoiseLite noise;
    
    // Using simple ProceduralPath for motion
    ProceduralPath<std::function<Vector(float)>> random_path{
        [](float t) { return Vector(0,1,0); } // Placeholder
    };
    
    struct SimpleMesh {
        std::vector<Vector> vertices;
        std::vector<std::vector<int>> faces;
    };
    SimpleMesh mesh; 
    
    Pipeline<W, FilterOrient<W>, FilterTemporal<W, 200>, FilterAntiAlias<W>> filters;

    void rebuild_mesh() {
        // Load Icosahedron
        using S = Icosahedron;
        // Vertices
        for(const auto& v : S::vertices) mesh.vertices.push_back(v);
        // Faces
        // Faces
        int offset = 0;
        for(uint8_t count : S::face_counts) {
            std::vector<int> face;
            for(int i = 0; i < count; ++i) {
                face.push_back(S::faces[offset + i]);
            }
            mesh.faces.push_back(face);
            offset += count;
        }
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
