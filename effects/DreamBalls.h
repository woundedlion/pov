/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include <map>
#include <memory> 

template <int W, int H>
class DreamBalls : public Effect {
public:
    enum class SolidType {
        Rhombicuboctahedron,
        Rhombicosidodecahedron,
        TruncatedCuboctahedron,
        Icosidodecahedron
    };

    struct Params {
        SolidType solid_type;
        int num_copies;
        float offset_radius;
        float offset_speed;
        float warp_scale;
        const Palette* palette;
        float alpha;
        bool enable_slice = false;
    };

    DreamBalls() : Effect(W, H), 
        filters(Filter::World::OrientSlice<W>(orientations, Y_AXIS), Filter::Screen::AntiAlias<W, H>()),
        slice_filter(filters) // Filters inherits Head (FilterOrientSlice)
    {
        persist_pixels = false;

        // Initialize Orientations
        orientations.resize(2); // 2 slices

        // Initialize Presets
        setup_presets();
        
        // Start Sequence
        timeline.add(0, Animation::PeriodicTimer(160, [this](Canvas& c) { this->spin_slices(); }, true));
        timeline.add(9, Animation::RandomWalk<W>(global_orientation, Y_AXIS, Animation::RandomWalk<W>::Options::Languid()));
        
        spawn_sprite(0); // Start first preset
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        
        slice_filter.get().enabled = false; 
        
        t += 0.01f;
        timeline.step(canvas);
    }
    
    // Exposed for consistency
    Params params; 

private:
    float t = 0;

    struct Tangent { Vector u; Vector v; };
    
    struct PresetData {
        // Base Mesh Data
        MeshState mesh_state; 
        std::vector<Tangent> tangents;
    };
    
    std::vector<PresetData> loaded_presets; 

    // Timeline
    Timeline<W> timeline;
    
    // Slices
    std::vector<Orientation<W>> orientations;
    Orientation<W> global_orientation; 
    
    // Pipeline
    Pipeline<W, H, Filter::World::OrientSlice<W>, Filter::Screen::AntiAlias<W, H>> filters;
    std::reference_wrapper<Filter::World::OrientSlice<W>> slice_filter;

    // Presets Definition
    std::vector<Params> presets;
    
    AlphaFalloffPalette bloodStreamFalloff {
        [](float t){ return 1.0f - t; }, 
        Palettes::bloodStream
    };

    void setup_presets() {
        // Define Params
        presets.push_back({SolidType::Rhombicuboctahedron, 18, 0.3f, 0.4f, 0.3f, &bloodStreamFalloff, 0.7f});
        presets.push_back({SolidType::Rhombicosidodecahedron, 6, 0.05f, 1.0f, 1.8f, &bloodStreamFalloff, 0.7f});
        presets.push_back({SolidType::TruncatedCuboctahedron, 6, 0.16f, 1.0f, 2.0f, &Palettes::richSunset, 0.3f});
        presets.push_back({SolidType::Icosidodecahedron, 10, 0.16f, 1.0f, 0.5f, &Palettes::lavenderLake, 0.3f});
        
        // Pre-load Geometry
        loaded_presets.reserve(presets.size());
        
        for(const auto& p : presets) {
            loaded_presets.emplace_back();
            auto& data = loaded_presets.back();
            PolyMesh m;
            switch(p.solid_type) {
                case SolidType::Rhombicuboctahedron: m = Solids::Archimedean::rhombicuboctahedron(); break;
                case SolidType::Rhombicosidodecahedron: m = Solids::Archimedean::rhombicosidodecahedron(); break;
                case SolidType::TruncatedCuboctahedron: m = Solids::Archimedean::truncatedCuboctahedron(); break;
                case SolidType::Icosidodecahedron: m = Solids::Archimedean::icosidodecahedron(); break;
            }
            
            // Store Verts
            data.mesh_state.vertices = m.vertices;
            
            // Store Faces
            data.mesh_state.faces.clear();
            data.mesh_state.face_counts.clear();
            for(const auto& f : m.faces) {
                data.mesh_state.face_counts.push_back((uint8_t)f.size());
                for(int idx : f) data.mesh_state.faces.push_back(idx);
            }
            
            // Compute Tangents
            for(const auto& v : data.mesh_state.vertices) {
                Vector axis = (std::abs(v.j) > 0.99f) ? X_AXIS : Y_AXIS;
                Vector u = cross(v, axis).normalize();
                Vector frame_v = cross(v, u).normalize();
                data.tangents.push_back({u, frame_v});
            }
        }
    }

    void spawn_sprite(int idx) {
        int safe_idx = idx % presets.size();
        const auto& preset_params = presets[safe_idx];
        const auto& preset_data = loaded_presets[safe_idx]; 
        
        auto mobius = std::make_shared<MobiusParams>(); 
        auto warp = std::make_shared<Animation::MobiusWarp>(*mobius, 1.0f, 200, true);
        warp->scale = preset_params.warp_scale;
        
        auto draw_fn = [this, preset_params, preset_data, mobius, warp](Canvas& canvas, float opacity) mutable {
            warp->step(canvas); 
            
            // Stack allocated MeshState
            MeshState target_mesh = preset_data.mesh_state; // Deep copy
            
            this->draw_scene(canvas, preset_params, opacity, preset_data.mesh_state, target_mesh, preset_data.tangents, *mobius);
        };
        
        timeline.add(0, Animation::Sprite(draw_fn, 320, 32, ease_mid, 32, ease_mid));
        
        // Queue next
        timeline.add(320 - 32, Animation::PeriodicTimer(0, [this, idx](Canvas& c) { this->spawn_sprite(idx + 1); }, false));
    }

    void update_displaced_mesh(const MeshState& base, MeshState& target, const std::vector<Tangent>& tangents, const Params& p, float angle_offset) {
        size_t count = base.vertices.size();
        float r = p.offset_radius;
        float speed = p.offset_speed;

        for(size_t i=0; i<count; ++i) {
            const Vector& v = base.vertices[i];
            const auto& tan = tangents[i];

            float phase = i * 0.1f;
            float angle = t * speed * 2 * PI_F + phase + angle_offset;

            float cosA = cosf(angle);
            float sinA = sinf(angle);

            Vector disp = v + (tan.u * cosA + tan.v * sinA) * r;
            target.vertices[i] = disp.normalize();
        }
    }

    void draw_scene(Canvas& canvas, const Params& p, float opacity, const MeshState& base, MeshState& target, const std::vector<Tangent>& tangents, const MobiusParams& m_params) {
        
        auto vertex_shader = [&](Fragment& f) {
            f.pos = mobius_transform(f.pos, m_params);
            f.pos = global_orientation.orient(f.pos);
        };

        auto fragment_shader = [&](const Vector& v, Fragment& f) {
            Color4 c = p.palette->get(f.v0); 
            c.alpha *= p.alpha * opacity;
            f.color = c;
        };

        for(int i=0; i<p.num_copies; ++i) {
            float offset = (static_cast<float>(i) / p.num_copies) * 2 * PI_F;
            update_displaced_mesh(base, target, tangents, p, offset);
            Plot::Mesh::draw<W, H>(filters, canvas, target, fragment_shader, vertex_shader);
        }
    }

    void spin_slices() {
        Vector axis = random_vector();
        slice_filter.get().axis = axis; 
        
        for (size_t i = 0; i < orientations.size(); ++i) {
            float direction = (i % 2 == 0) ? 1.0f : -1.0f;
            timeline.add(0, Animation::Rotation<W>(orientations[i], axis, direction * 2 * PI_F, 80, ease_in_out_sin, false));
        }
    }
};
