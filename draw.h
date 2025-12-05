#pragma once
#include <vector>
#include <algorithm>

/**
 * @brief Draws a geodesic line (arc) between two vectors on the sphere with adaptive sampling.
 * @tparam W The width of the display (used for step calculation).
 * @param dots The buffer to which the resulting Dots will be added.
 * @param v1 The start vector.
 * @param v2 The end vector.
 * @param color_fn Function to determine the color (takes vector and normalized progress t).
 * @param start Starting angle multiplier for drawing the line arc.
 * @param end Ending multiplier for the total arc angle.
 * @param long_way If true, draws the longer arc.
 * @param omit_last If true, omits the last point (useful for connecting segments).
 */
template <int W>
void draw_line(Dots& dots, const Vector& v1, const Vector& v2, ColorFn auto color_fn,
    float start, float end, bool long_way, bool omit_last)
{
    Vector u(v1);
    Vector v(v2);
    float a = angle_between(u, v);
    Vector w;

    if (std::abs(a) < TOLERANCE) {
        if (!omit_last) {
            dots.emplace_back(Dot(u, color_fn(u, 0.0f)));
        }
        return;
    }
    else if (std::abs(PI_F - a) < TOLERANCE) {
        if (std::abs(dot(v, X_AXIS)) > 0.9999f) {
            w = cross(u, Y_AXIS).normalize();
        }
        else {
            w = cross(u, X_AXIS).normalize();
        }
    }
    else {
        w = cross(u, v).normalize();
    }

    if (long_way) {
        a = 2 * PI_F - a;
        w = -w;
    }

    if (std::abs(start) > TOLERANCE) {
        Quaternion q = make_rotation(w, start * a);
        u = rotate(u, q).normalize();
    }
    a *= std::abs(end - start);

    // Simulation P\phase
    Vector sim_u = u;
    float sim_angle = 0;
    float steps[W * 4]; // Safety margin to accomodate worst case of meridian ring
    size_t step_count = 0;
    const float base_step = 2 * PI_F / W;
    while (sim_angle < a) {
        float scale_factor = std::max(0.05f, sqrtf(std::max(0.0f, 1.0f - sim_u.j * sim_u.j)));
        float step = base_step * scale_factor;
        if (step_count < W * 4) {
            steps[step_count++] = step;
        }
        else {
			Serial.println("draw_line: Exceeded maximum step count during simulation phase.");
            break;
        }
        sim_angle += step;
        Quaternion q = make_rotation(w, step);
        sim_u = rotate(sim_u, q).normalize();
    }

    // Drawing phase
    float scale = a / sim_angle;
    if (step_count == 0) {
        return;
    }
    float current_angle = 0;
    dots.emplace_back(Dot(u, color_fn(u, 0)));
    size_t loop_limit = omit_last ? step_count - 1 : step_count;
    for (size_t i = 0; i < loop_limit; i++) {
        float step = steps[i] * scale;
        Quaternion q = make_rotation(w, step);
        u = rotate(u, q).normalize();
        current_angle += step;
        float t = (a > 0) ? (current_angle / a) : 1;
        dots.emplace_back(Dot(u, color_fn(u, t)));
    }
}

/**
 * @brief Rasterizes a list of points into Dot objects by connecting them with geodesic lines.
 * @param dots The buffer to which the resulting Dots will be added.
 * @param points The buffer of points.
 * @param color_fn Function to determine color (takes vector and normalized progress t).
 * @param close_loop If true, connects the last point to the first.
 */
template <int W>
void rasterize(Dots& dots, const Points& points, ColorFn auto color_fn, bool close_loop = false) {
    size_t len = points.size();
    if (len == 0) return;

    size_t count = close_loop ? len : len - 1;
    for (size_t i = 0; i < count; i++) {
        const Vector& p1 = points[i];
        const Vector& p2 = points[(i + 1) % len];

        auto segment_color_fn = [&](const Vector& p, float sub_t) {
            float global_t = (i + sub_t) / count;
            return color_fn(p, global_t);
         };

        draw_line<W>(dots, p1, p2, segment_color_fn, 0.0f, 1.0f, false, true);
    }
}

/**
 * @brief Draws a single vector (point) on the unit sphere.
 * @tparam W The width of the display (used for projection context).
 * @param dots The buffer to which the Dot will be added.
 * @param v The 3D vector position (will be normalized).
 * @param color_fn The color function to determine the pixel's color.
 */
template <int W>
void draw_vector(Dots& dots, const Vector& v, ColorFn auto color_fn) {
    Vector u(v);
    u.normalize();
    dots.emplace_back(Dot(u, color_fn(u, 0)));
}

/**
 * @brief Draws a set of vertices (points) without connecting them.
 * @param dots The buffer to which the resulting Dots will be added.
 * @param points The list of vertex vectors.
 * @param color_fn The color function.
 */
void draw_vertices(Dots& dots, const Points& points, ColorFn auto color_fn) {
    for (Vector v : points) {
        dots.emplace_back(Dot(v.normalize(), color_fn(v, 0)));
    }
}

/**
 * @brief Draws the wireframe of a polyhedron by connecting vertices based on adjacency.
 * @tparam W The width of the display.
 * @param dots The buffer to which the resulting Dots will be added.
 * @param vertices The list of vertex vectors.
 * @param edges The adjacency list defining which vertices are connected.
 * @param color_fn The color function.
 */
template <int W>
void draw_polyhedron(Dots& dots, const VertexList& vertices, const AdjacencyList& edges, ColorFn auto color_fn) {
    for (size_t i = 0; i < edges.size(); ++i) {
        Vector a(vertices[i]);
        for (auto j : edges[i]) {
            if (i < j) {
                Vector b(vertices[j]);
                draw_line<W>(dots, a, b, color_fn, 0, 1, false, true);
            }
        }
    }
}

/**
 * @brief Calculates a point on a circle that lies on the surface of the unit sphere.
 * Used internally by drawing functions.
 */
Vector calc_ring_point(float a, float radius, const Vector& u, const Vector& v, const Vector& w) {
    auto d = sqrtf((1 - radius) * (1 - radius));
    return Vector(
        d * v.i + radius * u.i * cosf(a) + radius * w.i * sinf(a),
        d * v.j + radius * u.j * cosf(a) + radius * w.j * sinf(a),
        d * v.k + radius * u.k * cosf(a) + radius * w.k * sinf(a)
    ).normalize();
}

/**
 * @brief Calculates a single point on a sphere distorted by a function, often for an oscillating ring.
 */
Vector fn_point(ScalarFn auto f, const Vector& normal, float radius, float angle) {
    Vector v(normal);
    if (radius > 1) {
        v = -v;
        radius = 2 - radius;
    }
    Vector u;
    if (std::abs(dot(v, X_AXIS)) > 0.99995f) {
        u = cross(v, Y_AXIS).normalize();
    }
    else {
        u = cross(v, X_AXIS).normalize();
    }
    Vector w(cross(v, u));
    auto d = sqrtf((1 - radius) * (1 - radius));

    auto vi = calc_ring_point(angle, radius, u, v, w);
    auto vp = calc_ring_point(angle, 1, u, v, w);
    Vector axis = cross(v, vp).normalize();
    auto shift = make_rotation(axis, f(angle * PI_F / 2));
    return rotate(vi, shift);
};

/**
 * @brief Samples points for a circular ring on the sphere surface with adaptive sampling.
 * @param points The buffer to which points are added.
 * @param orientationQuaternion The orientation of the ring.
 * @param normal The normal vector defining the ring plane.
 * @param radius The radius of the ring.
 * @param phase Starting phase.
 */
template<int W>
void sample_ring(Points& points, const Quaternion& orientationQuaternion, const Vector& normal, float radius, float phase = 0) {
    // Basis
    Vector ref_axis = X_AXIS;
    if (std::abs(dot(normal, ref_axis)) > 0.9999f) {
        ref_axis = Y_AXIS;
    }
    Vector v = rotate(normal, orientationQuaternion).normalize();
    Vector ref = rotate(ref_axis, orientationQuaternion).normalize();
    Vector u = cross(v, ref).normalize();
    Vector w = cross(v, u).normalize();

    // Backside rings
    Vector v_dir = v;
    if (radius > 1.0f) {
        v_dir = -v_dir;
        radius = 2.0f - radius;
    }

    const float theta_eq = radius * (PI_F / 2.0f);
    const float r = sinf(theta_eq);
    const float d = cosf(theta_eq);

    // Calculate Samples
    const int num_samples = W / 4;
    const float step = 2.0f * PI_F / num_samples;
    Vector u_temp;

    for (int i = 0; i < num_samples; i++) {
        float theta = i * step;
        float t = theta + phase;
        float cos_ring = cosf(t);
        float sin_ring = sinf(t);
        u_temp = (u * cos_ring) + (w * sin_ring);
        Vector p = ((v_dir * d) + (u_temp * r)).normalize();
        points.push_back(p);
    }
}

/**
 * @brief Draws a circular ring on the sphere surface with adaptive sampling
 * to prevent artifacts near the poles.
 */
template<int W>
void draw_ring(Dots& dots, const Quaternion& orientation, const Vector& normal,
    float radius, ColorFn auto color_fn, float phase = 0)
{
    Points points;
    sample_ring<W>(points, orientation, normal, radius, phase);
    rasterize<W>(dots, points, color_fn, true);
}

/**
 * @brief Samples points for a function-distorted ring with adaptive sampling.
 * @param points The buffer to which points are added.
 * @param orientationQuaternion Orientation of the base ring.
 * @param normal Normal of the base ring.
 * @param radius Base radius (0-1).
 * @param shift_fn Function(t) returning angle offset.
 * @param phase Starting phase offset.
 */
template <int W>
void sample_fn(Points& points, const Quaternion& orientationQuaternion, const Vector& normal, float radius, ScalarFn auto shift_fn, float phase = 0) {
    // Basis
    Vector ref_axis = X_AXIS;
    if (std::abs(dot(normal, ref_axis)) > 0.9999f) {
        ref_axis = Y_AXIS;
    }
    Vector v = rotate(normal, orientationQuaternion).normalize();
    Vector ref = rotate(ref_axis, orientationQuaternion).normalize();
    Vector u = cross(v, ref).normalize();
    Vector w = cross(v, u).normalize();

    // Backside rings
    float v_sign = 1.0f;
    if (radius > 1.0f) {
        v_sign = -1.0f;
        radius = 2.0f - radius;
    }

    // Equidistant projection
    const float theta_eq = radius * (PI_F / 2.0f);
    const float r = sinf(theta_eq);
    const float d = cosf(theta_eq);

    // Calculate Samples
    const int num_samples = W;
    const float step = 2.0f * PI_F / num_samples;
    Vector u_temp;

    for (int i = 0; i < num_samples; i++) {
        float theta = i * step;
        float t = theta + phase;
        float cos_ring = cosf(t);
        float sin_ring = sinf(t);
        u_temp = (u * cos_ring) + (w * sin_ring);

        // Apply Shift
        float shift = shift_fn(theta / (2.0f * PI_F));
        float cos_shift = cosf(shift);
        float sin_shift = sinf(shift);
        float v_scale = (v_sign * d) * cos_shift - r * sin_shift;
        float u_scale = r * cos_shift + (v_sign * d) * sin_shift;
        Vector p = ((v * v_scale) + (u_temp * u_scale)).normalize();

        points.push_back(p);
    }
}

/**
*@brief Draws a function-distorted ring with adaptive sampling.
*/
template <int W>
void draw_fn(Dots& dots, const Quaternion& orientationQuaternion, const Vector& normal,
    float radius, ScalarFn auto shift_fn, ColorFn auto color_fn, float phase = 0)
{
    Points points;
    sample_fn<W>(points, orientationQuaternion, normal, radius, shift_fn, phase);
    rasterize<W>(dots, points, color_fn, true);
}

/**
 * @brief Calculates a single point on a ring based on angle and phase.
 */
Vector ring_point(const Vector& normal, float radius, float angle, float phase = 0) {
    Vector v(normal);
    if (radius > 1) {
        v = -v;
    }
    Vector u;
    if (std::abs(dot(v, X_AXIS)) > 0.99995f) {
        u = cross(v, Y_AXIS).normalize();
    }
    else {
        u = cross(v, X_AXIS).normalize();
    }
    Vector w(cross(v, u));
    if (radius > 1) {
        w = -w;
        radius = 2 - radius;
    }
    auto d = sqrtf((1 - radius) * (1 - radius));
    return Vector(
        d * v.i + radius * u.i * cosf(angle + phase) + radius * w.i * sinf(angle + phase),
        d * v.j + radius * u.j * cosf(angle + phase) + radius * w.j * sinf(angle + phase),
        d * v.k + radius * u.k * cosf(angle + phase) + radius * w.k * sinf(angle + phase)
    ).normalize();
};

/**
 * @brief Represents a customizable path or trajectory in 3D space.
 * @tparam W The width of the display (used for samples/resolution).
 */
template <int W>
class Path {
public:
    Path() {}

    /**
     * @brief Appends a great-circle arc segment to the path.
     */
    Path& append_line(const Vector& v1, const Vector& v2, bool long_way = false) {
        if (points.size() > 0) {
            points.pop_back(); // Overlap previous segment
        }
        Dots seg;
        draw_line<W>(seg, v1, v2, [](auto&, auto) { return Pixel(); }, 0.0f, 1.0f, long_way, false);
        std::transform(seg.begin(), seg.end(), std::back_inserter(points),
            [](auto& d) { return d.position; });
        return *this;
    }

    /**
     * @brief Appends a segment generated by a plotting function.
     */
    Path& append_segment(PlotFn auto plot, float domain, float samples, ScalarFn auto easing) {
        if (points.size() > 0) {
            points.pop_back(); // Overlap previous segment
        }
        for (float t = 0; t <= samples; t++) {
            points.push_back(plot(easing(t / samples) * domain));
        }
        return *this;
    }

    /**
     * @brief Retrieves a point along the path based on a normalized factor.
     */
    Vector get_point(float t) const {
        return points[static_cast<int>(t * (points.size() - 1))];
    }

    /**
     * @brief Gets the total number of discrete points in the path.
     */
    size_t num_points() const { return points.size(); }

    /**
     * @brief Reduces the path to a single point: the last recorded position.
     */
    void collapse() {
        if (points.size() > 1) {
            points = { points.back() };
        }
    }

private:

    StaticCircularBuffer<Vector, 1024> points; /**< The discrete points making up the path. */
};

/**
 * @brief Draws a Path object by sampling its points and adding them as dots.
 */
template <int W>
void draw_path(Dots& dots, const Path<W>& path, ColorFn auto color) {
    size_t samples = path.num_points();
    for (size_t i = 0; i < samples; ++i) {
        auto v = path.get_point(static_cast<float>(i) / samples);
        dots.push_back(Dot(v, color(v, i / (samples - 1))));
    }
}

/**
 * @brief Applies a discrete square wave mask to a color, simulating a dotted or striped brush.
 */
Pixel dotted_brush(const Pixel& color, float freq, float duty_cycle, float phase, float t) {
    return dim(color, square_wave(0, 1, freq, duty_cycle, phase)(t));
}

/**
 * @brief Plots a buffer of Dots to the Canvas, routing the data through the filter chain.
 */
template<int W, typename Pipeline>
void plot_dots(const Dots& dots, Pipeline& filters, Canvas& canvas, float age, float alpha) {
    for (auto& dot : dots) {
        auto p = vector_to_pixel<W>(dot.position);
        filters.plot(canvas, p.x, p.y, gamma_correct(dot.color), age, alpha);
    }
}

/**
 * @brief Iterates over the Orientation's history and calls a draw function for each historical step.
 */
void tween(const Orientation& orientation, TweenFn auto draw_fn) {
    size_t s = orientation.length();
    size_t start = (s > 1) ? 1 : 0;
    for (size_t i = start; i < s; ++i) {
        draw_fn(orientation.get(i),
            static_cast<float>((s - 1 - i)) / s);
    }
}