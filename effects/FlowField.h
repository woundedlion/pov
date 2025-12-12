#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"
#include "../FastNoiseLite.h"

template <int W>
class FlowField : public Effect {
public:
    FlowField() :
        Effect(W),
        palette(iceMelt),
        trails(k_trail_length),
        filters(
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;
        noise_generator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise_generator.SetSeed(hs::rand_int(0, 65535));
        reset_particles();
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        dots.clear();

        for (Particle& p : particles) {
            // 1. Get acceleration from the 3D noise field (simulating 4D)
            Vector accel = get_noise_force(p.pos, timeline.t);

            // 2. Apply a gravity force pulling the particle to the center
            Vector gravity_force = p.pos * -k_gravity;
            accel = accel + gravity_force;

            // 3. Update velocity
            p.vel = p.vel + accel;

            // 4. Clamp velocity to max speed (manual implementation of clampLength)
            float speed = p.vel.length();
            if (speed > k_max_speed) {
                p.vel = (p.vel / speed) * k_max_speed;
            }

            // 5. Update position and re-normalize to keep it on the sphere
            p.pos = p.pos + p.vel;
            p.pos.normalize();

            // 6. Get color based on speed (ratio of current speed to max)
            float speed_ratio = speed / k_max_speed;
            Pixel color = palette.get(speed_ratio);

            // 7. Add the particle's "head" to the dot buffer
            dots.emplace_back(Dot(p.pos, color));
        }

        trails.record(dots, 0, k_alpha);

        trails.render(canvas, filters,
            [this](const Vector& v, float t) {
                return palette.get(1.0f - t) * (1.0f - t);
            }
        );
    }

private:

    struct Particle {
        Vector pos;
        Vector vel;

        Particle() : pos(random_vector()), vel(0, 0, 0) {}
    };

    static constexpr int k_num_particles = 100;
    static constexpr int k_trail_length = 8;
    static constexpr float k_noise_scale = 100;
    static constexpr float k_force_scale = 0.001;
    static constexpr float k_gravity = 0.001;
    static constexpr float k_max_speed = 0.1;
    static constexpr float k_alpha = 0.7;

    static constexpr float k_time_scale = 0.01;
    static constexpr int k_max_trail_dots = 1024;
    Timeline timeline;
    FastNoiseLite noise_generator;
    const ProceduralPalette& palette;
    std::array<Particle, k_num_particles> particles;

    Orientation orientation;
    DecayBuffer<W, k_max_trail_dots> trails;

    Pipeline<W,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
    Dots dots;

    void reset_particles() {
        for (int i = 0; i < k_num_particles; ++i) {
            particles[i] = Particle();
        }
    }

    Vector get_noise_force(const Vector& pos, float t) {
        float t_scaled = t * k_time_scale;
        Vector n_pos = pos * k_noise_scale;

        float x_force = noise_generator.GetNoise(n_pos.i, n_pos.j, t_scaled);
        float y_force = noise_generator.GetNoise(n_pos.j, n_pos.k, t_scaled + 100.0f);
        float z_force = noise_generator.GetNoise(n_pos.k, n_pos.i, t_scaled + 200.0f);

        return Vector(x_force, y_force, z_force) * k_force_scale;
    }
};
