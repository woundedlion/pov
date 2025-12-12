#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

template <int W>
class MobiusGrid : public Effect {
public:
    MobiusGrid() :
        Effect(W),
        palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
        next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
        holeN(Z_AXIS),
        holeS(-Z_AXIS),
        filters(
            FilterHole<W>(holeN, 1.2f),
            FilterHole<W>(holeS, 1.2f),
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;

        timeline
            .add(0, MobiusWarp(params, num_rings, 160, true))
            .add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400, ease_mid, true))
            .add(0, PeriodicTimer(120, [this](auto&) { wipe_palette(); }, true))
            .add(0, Mutation(num_rings, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true))
            .add(160, Mutation(num_lines, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true));
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        dots.clear();

        float phase = fmodf(static_cast<float>(timeline.t), 120.0f) / 120.0f;

        draw_axis_rings(Z_AXIS, num_rings, phase);
        draw_longitudes(num_lines, phase);

        // Calculate Stabilizing Counter-Rotation
        Vector n_in = Z_AXIS;
        Vector n_trans = inv_stereo(mobius(stereo(n_in), params));
        Vector s_in = -Z_AXIS;
        Vector s_trans = inv_stereo(mobius(stereo(s_in), params));
        Vector mid = (n_trans + s_trans).normalize();
        Quaternion q = make_rotation(mid, Z_AXIS);

        // Apply counter-rotation
        for (auto& dot : dots) {
            dot.position = rotate(dot.position, q).normalize();
        }

        // Update hole origins to match the rotated geometry
        holeN = rotate(n_trans, q).normalize();
        holeS = rotate(s_trans, q).normalize();

        plot_dots<W>(dots, filters, canvas, 0, alpha);
    }

private:
    void wipe_palette() {
        next_palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT);
        timeline.add(0, ColorWipe(palette, next_palette, 60, ease_mid));
    }

    void draw_axis_rings(const Vector& normal, float num, float phase) {
        const float log_min = -2.5f;
        const float log_max = 2.5f;
        const float range = log_max - log_min;
        int count = static_cast<int>(std::ceil(num));

        for (int i = 0; i < count; ++i) {
            float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
            float log_r = log_min + t * range;
            float r_val = expf(log_r);
            float radius = (4.0f / PI_F) * atanf(1.0f / r_val);

            Points points;
            sample_ring<W>(points, Quaternion(), normal, radius);

            Points transformed_points;
            for (const auto& p : points) {
                transformed_points.push_back(inv_stereo(mobius(stereo(p), params)));
            }

            float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
            rasterize<W>(dots, transformed_points, [&](const Vector&, float) {
                return palette.get(static_cast<float>(i) / num) * opacity;
                }, true);
        }
    }

    void draw_longitudes(float num, float phase) {
        int count = static_cast<int>(std::ceil(num));
        for (int i = 0; i < count; ++i) {
            float theta = (static_cast<float>(i) / num) * PI_F;
            Vector normal(cosf(theta), sinf(theta), 0.0f);

            Points points;
            sample_ring<W>(points, Quaternion(), normal, 1.0f);
            Points transformed_points;
            for (const auto& p : points) {
                transformed_points.push_back(inv_stereo(mobius(stereo(p), params)));
            }

            float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
            rasterize<W>(dots, transformed_points, [&](const Vector&, float t_line) {
                // Approximate original Z to calculate log-gradient
                float original_idx = t_line * points.size();
                int idx1 = static_cast<int>(original_idx) % points.size();
                int idx2 = (idx1 + 1) % points.size();
                float f = original_idx - std::floor(original_idx);
                float z = points[idx1].k * (1.0f - f) + points[idx2].k * f;

                float R = sqrtf((1.0f + z) / (1.0f - z));
                float log_r = logf(R);
                const float log_min = -2.5f;
                const float log_max = 2.5f;
                float t = (log_r - log_min) / (log_max - log_min);

                return palette.get(wrap(t - phase, 1.0f)) * opacity;
                }, true);
        }
    }

    float alpha = 0.2f;
    float num_rings = 0;
    float num_lines = 0;
    GenerativePalette palette;
    GenerativePalette next_palette;
    MobiusParams params;
    Orientation orientation;
    Timeline timeline;
    Dots dots;

    Vector holeN;
    Vector holeS;

    Pipeline<W,
        FilterHoleRef<W>,
        FilterHoleRef<W>,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
};
