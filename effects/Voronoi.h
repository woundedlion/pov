/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../animation.h"
#include "../palettes.h"
#include "../filter.h"
#include <vector>
#include <cmath>
#include <cstdlib>

template <int W>
class Voronoi : public Effect {
public:
    struct Site {
        Vector pos;
        Vector axis;
        Color4 color;
        int id;
    };

    Voronoi() : Effect(W),
        sites_buffer(200) 
    {
        init_sites();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        
        // 1. Animate Sites
        float s = logf(speed + 1.0f) * 0.005f;
        
        for(size_t i=0; i<static_cast<size_t>(num_sites); ++i) {
            auto& site = sites_buffer[i];
            // Rotate pos around axis
            Quaternion q = make_rotation(site.axis, s);
            site.pos = rotate(site.pos, q);
        }
        
        // 2. Render Pixels (Linear Search)
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                Vector p = pixel_to_vector<W>(x, y);
                
                int bestIdx = -1;
                float bestDot = -2.0f;
                int secondBestIdx = -1;
                float secondBestDot = -2.0f;
                
                for(size_t i=0; i<static_cast<size_t>(num_sites); ++i) {
                     float d = dot(p, sites_buffer[i].pos);
                     if (d > bestDot) {
                         secondBestDot = bestDot;
                         secondBestIdx = bestIdx;
                         bestDot = d;
                         bestIdx = i;
                     } else if (d > secondBestDot) {
                         secondBestDot = d;
                         secondBestIdx = i;
                     }
                }
                
                float maxDot1 = bestDot;
                float maxDot2 = secondBestDot;
                const auto& bestSite = sites_buffer[bestIdx];
                
                Color4 c = bestSite.color;
                
                // Smoothing
                if (secondBestIdx != -1 && smoothness > 0.0f) {
                    const auto& secSite = sites_buffer[secondBestIdx];
                    float diff = maxDot1 - maxDot2;
                    float factor = std::min(1.0f, diff * smoothness);
                    factor = quintic_kernel(factor);
                    float t = 0.5f + 0.5f * factor;
                    
                    // Mix
                    c.color.r = static_cast<uint8_t>(secSite.color.color.r + (c.color.r - secSite.color.color.r) * t);
                    c.color.g = static_cast<uint8_t>(secSite.color.color.g + (c.color.g - secSite.color.color.g) * t);
                    c.color.b = static_cast<uint8_t>(secSite.color.color.b + (c.color.b - secSite.color.color.b) * t);
                }
                
                // Borders
                if (showBorders && secondBestIdx != -1) {
                    float dist1 = acosf(std::min(1.0f, maxDot1));
                    float dist2 = acosf(std::min(1.0f, maxDot2));
                    if (dist2 - dist1 < borderThickness) {
                        c = Color4(0,0,0,0); // Black/Transparent
                    }
                }
                
                canvas(x, y) = c.color; 
            }
        }
    }
    
    // Params
    int num_sites = 20; 
    float speed = 20.0f;
    bool showBorders = true;
    float borderThickness = 0.0f;
    float smoothness = 100.0f;
    bool showSites = false;

private:
    std::vector<Site> sites_buffer;
    
    void init_sites() {
        sites_buffer.clear();
        sites_buffer.resize(num_sites);
        
        for (int i = 0; i < num_sites; i++) {
            float goldenAngle = PI_F * (3.0f - sqrtf(5.0f));
            float y = 1.0f - (i / (float)(num_sites - 1)) * 2.0f;
            float radius = sqrtf(std::max(0.0f, 1.0f - y * y));
            float theta = goldenAngle * i;

            float x = cosf(theta) * radius;
            float z = sinf(theta) * radius;

            sites_buffer[i].pos = Vector(x, y, z);
            
            float rx = (rand() % 1000) / 500.0f - 1.0f;
            float ry = (rand() % 1000) / 500.0f - 1.0f;
            float rz = (rand() % 1000) / 500.0f - 1.0f;
            sites_buffer[i].axis = Vector(rx, ry, rz).normalize();
            
            float t = i / (float)(num_sites > 1 ? num_sites - 1 : 1);
            sites_buffer[i].color = Palettes::richSunset.get(t);
            sites_buffer[i].id = i;
        }
    }
};
