/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief White-box access to MindSplatter state.
 */
#pragma once

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include "effects/MindSplatter.h"
#include "tests/temporal_oracle.h"
#include "tests/test_fixture.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace effects_tests {

/** @brief White-box accessor for MindSplatter internals. */
struct MindSplatterWhiteBox {
  using MS = MindSplatter<288, 144>;

  /** @brief Renders production column trails or their frozen geometry.
   * @param timed Uses the production temporal draw.
   * @param fraction Frozen frame fraction.
   * @param displacement Vertical frozen-reference offset in pixels.
   * @param turn Each head's turn over the frame, radians.
   * @param shift Warp translation change over the frame.
   * @param spacing Columns between a trail's anchors.
   * @param setups Receives the production draw's column setups, when given.
   * @return Row-major linear luminance.
   * @details Four particles on the equator carry full trails under a moving
   * warp and orientation. The frozen reference draws each trail through the
   * column walk without motion: anchors under the warp at @p fraction, the
   * head turned back by the rest of its step, under the pose at @p fraction.
   * Anchors a column apart are steady single samples, which the point route
   * plots to its own tolerance, so their displaced reference also moves a
   * column's @p displacement along the rows.
   */
  static std::vector<float> render_temporal_trails(bool timed, float fraction,
                                                   float displacement,
                                                   float turn, float shift,
                                                   float spacing,
                                                   uint32_t *setups = nullptr) {
    reset_effect_globals();
    MS fx;
    fx.init();
    fx.history_ready = true;
    fx.particle_system.attractors.clear();
    fx.active_base_mesh = MS::BaseMesh::TETRAHEDRON;
    const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
    const auto start_pose = math::make_rotation(math::X_AXIS, 0.27f);
    fx.orientation.set(start_pose);
    fx.orientation.push(math::make_rotation(axis, 0.012f) * start_pose);
    fx.previous_warp = math::Complex(0.15f, 0.06f);
    fx.mobius.b = math::Complex(0.15f + shift, 0.06f + 0.4f * shift);
    constexpr uint16_t LIFE = 148;
    for (float x : {1.0f, 60.0f, 143.0f, 222.0f}) {
      const math::Vector head = math::pixel_to_vector<288, 144>(x, 70.0f);
      const math::Vector east = math::cross(math::Y_AXIS, head).normalized();
      const math::Vector north = math::cross(head, east);
      fx.particle_system.spawn(head, (east + north * 0.4f).normalized() * turn,
                               0);
      auto &p = fx.particle_system.pool[fx.particle_system.active() - 1];
      p.life = LIFE;
      for (int j = 0; j < 8; ++j)
        p.history.record(math::pixel_to_vector<288, 144>(
            x - spacing * float(8 - j), 70.0f + 0.3f * spacing * j));
    }
    {
      Canvas canvas(fx);
      if (timed) {
        Temporal::geometry_counters = {};
        fx.draw_particles(canvas);
        if (setups)
          *setups = Temporal::geometry_counters.column_draw_setups;
      } else {
        const math::Complex warp_step(fx.mobius.b.re - fx.previous_warp.re,
                                      fx.mobius.b.im - fx.previous_warp.im);
        math::MobiusParams warp = fx.mobius;
        warp.b = math::Complex(fx.previous_warp.re + warp_step.re * fraction,
                               fx.previous_warp.im + warp_step.im * fraction);
        const auto pose = fx.orientation.get(fraction);
        fx.orientation.set(pose);
        ScratchScope guard(scratch_arena_a);
        constexpr size_t POINTS = 9;
        const size_t count = POINTS * fx.particle_system.active();
        auto *vertices = scratch_arena_a.make_n<math::Vector>(count);
        auto *progress = scratch_arena_a.make_n<float>(count);
        ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_a, count);
        for (int i = 0; i < fx.particle_system.active(); ++i) {
          const auto &p = fx.particle_system.pool[i];
          for (size_t j = 0; j < POINTS; ++j) {
            math::Vector source =
                j + 1 < POINTS ? p.history.get(j) : p.position;
            if (j + 1 == POINTS && turn != 0.0f) {
              const float speed = p.velocity.magnitude();
              const float back = fminf(speed, fx.particle_system.motion_cap) *
                                 (1.0f - fraction);
              source = source * cosf(back) - p.velocity * (sinf(back) / speed);
            }
            math::Vector warped = math::mobius_transform(source, warp);
            if (displacement != 0.0f) {
              math::Vector world = hs_test::displace_reference_rows<144>(
                  math::rotate(warped, pose), displacement);
              if (spacing <= 1.0f)
                world = math::rotate(
                    world, math::make_rotation(
                               math::Y_AXIS,
                               displacement * math::RADIANS_PER_COLUMN<288>));
              warped = math::rotate(world, pose.conjugate());
            }
            vertices[i * POINTS + j] = warped;
            progress[i * POINTS + j] = static_cast<float>(j) / (POINTS - 1);
          }
          for (size_t j = POINTS - 1; j > 0; --j)
            edges.push_back({static_cast<uint16_t>(i * POINTS + j),
                             static_cast<uint16_t>(i * POINTS + j - 1)});
        }
        const float fade =
            static_cast<float>(LIFE) / fx.particle_system.max_life;
        auto shader = [&](const math::Vector &, Fragment &f) {
          const auto &edge = edges[static_cast<size_t>(f.v2)];
          const float t =
              progress[edge.u] + (progress[edge.v] - progress[edge.u]) * f.v0;
          f.color = Color4(
              MS::sample_trail_palette(MINDSPLATTER_PALETTES[0], 1.0f - t),
              fmaxf(0.0f, fminf(fmaxf(0.0f, t), fade)));
        };
        Plot::GeometryCallback geometry;
        geometry.source_positions = false;
        geometry.chained = true;
        Plot::Mesh::draw_column_edges<288, 144>(
            fx.timed_filters, canvas, typename MS::TrailMesh{{vertices, count}},
            edges, shader, geometry, 0.0f);
      }
    }
    fx.advance_display();
    std::vector<float> luminance(288 * 144);
    for (int y = 0; y < 144; ++y)
      for (int x = 0; x < 288; ++x) {
        const Pixel pixel = fx.get_pixel(x, y);
        luminance[y * 288 + x] =
            (pixel.r + pixel.g + pixel.b) / (3.0f * 65535.0f);
      }
    return luminance;
  }

  static constexpr uint32_t REPLAY_MAGIC = 0x3152534du;
  static constexpr uint16_t REPLAY_VERSION = 3;
  static constexpr size_t REPLAY_HEADER_BYTES = 4 + 8 * 2 + 1;

  template <typename T>
  using ObjectBytes = std::array<unsigned char, sizeof(T)>;

  template <int W, int H> struct ReplaySnapshot {
    using EffectType = MindSplatter<W, H>;
    using Particle = Animation::Particle<EffectType::TRAIL_LEN>;
    using Attractor = typename EffectType::ParticleSystem::Attractor;
    using BaseMesh = typename EffectType::BaseMesh;
    using Params = typename EffectType::Params;
    using Transition = typename EffectType::Transition;

    std::vector<ObjectBytes<Particle>> particles;
    std::vector<ObjectBytes<Attractor>> attractors;
    ObjectBytes<Params> params{};
    ObjectBytes<Transition> transition{};
    ObjectBytes<math::Orientation<>> orientation{};
    ObjectBytes<math::MobiusParams> mobius{};
    std::array<float, EffectType::MAX_EMITTERS> emit_phases{};
    uint8_t palette_sequence = 0;
    ClipRegion clip{};
    float friction = 0.0f;
    float gravity = 0.0f;
    uint16_t max_life = 0;
    BaseMesh active_base_mesh = BaseMesh::CUBE;
  };

  template <typename T>
  static void append_object(std::vector<unsigned char> &bytes,
                            const ObjectBytes<T> &object) {
    bytes.insert(bytes.end(), object.begin(), object.end());
  }

  static void append_u16(std::vector<unsigned char> &bytes, uint16_t value) {
    bytes.push_back(static_cast<unsigned char>(value));
    bytes.push_back(static_cast<unsigned char>(value >> 8));
  }

  static void append_u32(std::vector<unsigned char> &bytes, uint32_t value) {
    append_u16(bytes, static_cast<uint16_t>(value));
    append_u16(bytes, static_cast<uint16_t>(value >> 16));
  }

  template <int W, int H>
  static std::vector<unsigned char>
  serialize_render(const ReplaySnapshot<W, H> &snapshot) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    HS_CHECK(snapshot.particles.size() <= UINT16_MAX);
    static_assert(sizeof(Particle) <= UINT16_MAX);
    static_assert(sizeof(math::Orientation<>) <= UINT16_MAX);
    static_assert(sizeof(math::MobiusParams) <= UINT16_MAX);

    std::vector<unsigned char> bytes;
    bytes.reserve(REPLAY_HEADER_BYTES + sizeof(math::Orientation<>) +
                  sizeof(math::MobiusParams) +
                  snapshot.particles.size() * sizeof(Particle));
    append_u32(bytes, REPLAY_MAGIC);
    append_u16(bytes, REPLAY_VERSION);
    append_u16(bytes, W);
    append_u16(bytes, H);
    append_u16(bytes, static_cast<uint16_t>(sizeof(Particle)));
    append_u16(bytes, static_cast<uint16_t>(sizeof(math::Orientation<>)));
    append_u16(bytes, static_cast<uint16_t>(sizeof(math::MobiusParams)));
    append_u16(bytes, static_cast<uint16_t>(snapshot.particles.size()));
    append_u16(bytes, snapshot.max_life);
    bytes.push_back(static_cast<unsigned char>(snapshot.active_base_mesh));
    append_object<math::Orientation<>>(bytes, snapshot.orientation);
    append_object<math::MobiusParams>(bytes, snapshot.mobius);
    for (const ObjectBytes<Particle> &particle : snapshot.particles)
      append_object<Particle>(bytes, particle);
    return bytes;
  }

  class ReplayReader {
  public:
    explicit ReplayReader(std::span<const unsigned char> bytes)
        : bytes(bytes) {}

    uint16_t read_u16() {
      require(2);
      const uint16_t value = static_cast<uint16_t>(bytes[offset]) |
                             static_cast<uint16_t>(bytes[offset + 1]) << 8;
      offset += 2;
      return value;
    }

    uint32_t read_u32() {
      const uint32_t low = read_u16();
      return low | static_cast<uint32_t>(read_u16()) << 16;
    }

    uint8_t read_u8() {
      require(1);
      return bytes[offset++];
    }

    template <typename T> ObjectBytes<T> read_object() {
      require(sizeof(T));
      ObjectBytes<T> object;
      std::memcpy(object.data(), bytes.data() + offset, sizeof(T));
      offset += sizeof(T);
      return object;
    }

    bool empty() const { return offset == bytes.size(); }

  private:
    void require(size_t count) const {
      HS_CHECK(count <= bytes.size() - offset,
               "MindSplatter replay corpus is truncated");
    }

    std::span<const unsigned char> bytes;
    size_t offset = 0;
  };

  template <typename T> static ObjectBytes<T> object_bytes(const T &source) {
    static_assert(std::is_trivially_copyable_v<T>);
    ObjectBytes<T> bytes;
    std::memcpy(bytes.data(), &source, sizeof(T));
    return bytes;
  }

  template <typename T>
  static void restore_object(T &target, const ObjectBytes<T> &bytes) {
    static_assert(std::is_trivially_copyable_v<T>);
    std::memcpy(&target, bytes.data(), sizeof(T));
  }

  template <int W, int H>
  static ReplaySnapshot<W, H> capture(const MindSplatter<W, H> &ms) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    Snapshot snapshot;
    snapshot.particles.reserve(ms.particle_system.active());
    for (size_t i = 0; i < ms.particle_system.active(); ++i)
      snapshot.particles.push_back(object_bytes(ms.particle_system.pool[i]));
    snapshot.attractors.reserve(ms.particle_system.attractors.size());
    for (size_t i = 0; i < ms.particle_system.attractors.size(); ++i)
      snapshot.attractors.push_back(
          object_bytes(ms.particle_system.attractors[i]));
    snapshot.params = object_bytes(ms.params);
    snapshot.transition = object_bytes(ms.transition);
    snapshot.orientation = object_bytes(ms.orientation);
    snapshot.mobius = object_bytes(ms.mobius);
    snapshot.emit_phases = ms.emit_phases;
    snapshot.palette_sequence = ms.palette_sequence;
    snapshot.clip = ms.clip();
    snapshot.friction = ms.particle_system.friction;
    snapshot.gravity = ms.particle_system.gravity;
    snapshot.max_life = ms.particle_system.max_life;
    snapshot.active_base_mesh = ms.active_base_mesh;
    static_assert(std::is_trivially_copyable_v<Particle>);
    return snapshot;
  }

  template <int W, int H>
  static void restore(MindSplatter<W, H> &ms,
                      const ReplaySnapshot<W, H> &snapshot) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    HS_CHECK(ms.particle_system.active() == 0);
    HS_CHECK(snapshot.particles.size() <= ms.particle_system.pool.capacity());
    ms.configure_particle_geometry(snapshot.active_base_mesh);
    HS_CHECK(snapshot.attractors.size() ==
             ms.particle_system.attractors.size());

    restore_object(ms.params, snapshot.params);
    restore_object(ms.transition, snapshot.transition);
    restore_object(ms.orientation, snapshot.orientation);
    restore_object(ms.mobius, snapshot.mobius);
    ms.emit_phases = snapshot.emit_phases;
    ms.palette_sequence = snapshot.palette_sequence;
    ms.set_clip(snapshot.clip.y_start, snapshot.clip.y_end,
                snapshot.clip.x_start, snapshot.clip.x_end);
    ms.set_margin(snapshot.clip.margin);
    ms.particle_system.friction = snapshot.friction;
    ms.particle_system.gravity = snapshot.gravity;
    ms.particle_system.max_life = snapshot.max_life;

    for (size_t i = 0; i < snapshot.attractors.size(); ++i)
      restore_object(ms.particle_system.attractors[i], snapshot.attractors[i]);
    for (const ObjectBytes<Particle> &bytes : snapshot.particles) {
      const size_t i = ms.particle_system.active();
      ms.particle_system.spawn(math::Vector(), math::Vector(), 0);
      restore_object(ms.particle_system.pool[i], bytes);
    }
  }

  template <int W, int H>
  static uint16_t restore_render(MindSplatter<W, H> &ms,
                                 std::span<const unsigned char> bytes) {
    using Snapshot = ReplaySnapshot<W, H>;
    using Particle = typename Snapshot::Particle;
    ReplayReader reader(bytes);
    HS_CHECK(reader.read_u32() == REPLAY_MAGIC,
             "MindSplatter replay corpus magic differs");
    HS_CHECK(reader.read_u16() == REPLAY_VERSION,
             "MindSplatter replay corpus version differs");
    HS_CHECK(reader.read_u16() == W && reader.read_u16() == H,
             "MindSplatter replay corpus resolution differs");
    HS_CHECK(reader.read_u16() == sizeof(Particle),
             "MindSplatter replay particle layout differs");
    HS_CHECK(reader.read_u16() == sizeof(math::Orientation<>),
             "MindSplatter replay orientation layout differs");
    HS_CHECK(reader.read_u16() == sizeof(math::MobiusParams),
             "MindSplatter replay Mobius layout differs");
    const uint16_t particle_count = reader.read_u16();
    const uint16_t max_life = reader.read_u16();
    const auto active_base_mesh =
        static_cast<typename Snapshot::BaseMesh>(reader.read_u8());
    HS_CHECK(static_cast<size_t>(active_base_mesh) <
                 Solids::PLATONIC_BASE_MESH_COUNT,
             "MindSplatter replay base mesh %u is out of range",
             static_cast<unsigned>(active_base_mesh));
    HS_CHECK(ms.particle_system.active() == 0);
    HS_CHECK(particle_count <= ms.particle_system.pool.capacity());

    ms.configure_particle_geometry(active_base_mesh);
    restore_object(ms.orientation, reader.read_object<math::Orientation<>>());
    restore_object(ms.mobius, reader.read_object<math::MobiusParams>());
    ms.particle_system.max_life = max_life;
    for (uint16_t i = 0; i < particle_count; ++i) {
      ms.particle_system.spawn(math::Vector(), math::Vector(), 0);
      restore_object(ms.particle_system.pool[i],
                     reader.read_object<Particle>());
    }
    HS_CHECK(reader.empty(), "MindSplatter replay corpus has trailing bytes");

    ms.params.active_count = static_cast<float>(particle_count);
    return particle_count;
  }

  template <int W, int H> static void step_physics(MindSplatter<W, H> &ms) {
    Canvas canvas(ms);
    ms.particle_system.step(canvas);
    ms.params.active_count = static_cast<float>(ms.particle_system.active());
  }

  template <int W, int H>
  static void select_preset(MindSplatter<W, H> &ms, size_t index) {
    HS_CHECK(ms.selectPreset(index),
             "MindSplatter replay preset is out of range");
  }

  template <int W, int H>
  static void step_state_without_render(MindSplatter<W, H> &ms) {
    Canvas canvas(ms);
    ms.step_state(canvas);
  }

  template <int W, int H>
  static bool same_snapshot(const ReplaySnapshot<W, H> &a,
                            const ReplaySnapshot<W, H> &b) {
    if (a.particles != b.particles || a.attractors != b.attractors ||
        a.params != b.params || a.transition != b.transition ||
        a.orientation != b.orientation || a.mobius != b.mobius)
      return false;
    if (std::memcmp(a.emit_phases.data(), b.emit_phases.data(),
                    sizeof(a.emit_phases)) != 0 ||
        std::memcmp(&a.clip, &b.clip, sizeof(a.clip)) != 0 ||
        std::memcmp(&a.friction, &b.friction, sizeof(a.friction)) != 0 ||
        std::memcmp(&a.gravity, &b.gravity, sizeof(a.gravity)) != 0 ||
        a.max_life != b.max_life || a.palette_sequence != b.palette_sequence ||
        a.active_base_mesh != b.active_base_mesh)
      return false;
    return true;
  }

  template <int W, int H>
  static size_t num_emitters(const MindSplatter<W, H> &ms) {
    return MindSplatter<W, H>::emitter_vertices(ms.active_base_mesh).size();
  }
  template <int W, int H>
  static size_t active_emitters(const MindSplatter<W, H> &ms) {
    return ms.particle_system.emitters.size();
  }
  template <int W, int H>
  static size_t active_attractors(const MindSplatter<W, H> &ms) {
    return ms.particle_system.attractors.size();
  }
  template <int W, int H>
  static auto active_base_mesh(const MindSplatter<W, H> &ms) {
    return ms.active_base_mesh;
  }
  template <int W, int H>
  static auto preset_base_mesh(const MindSplatter<W, H> &ms, size_t index) {
    (void)ms;
    return MindSplatter<W, H>::PRESETS[index].params.base_mesh;
  }
  static float emit_phase(const MS &ms, size_t i) { return ms.emit_phases[i]; }
  static float event_horizon() { return MS::EVENT_HORIZON; }
  static float hole_alpha(const math::Vector &p) {
    return MS::octahedral_hole_alpha(p, math::fast_cosf(MS::EVENT_HORIZON));
  }
  static float reference_hole_alpha(const math::Vector &p) {
    const float cos_event_horizon = math::fast_cosf(MS::EVENT_HORIZON);
    float alpha = 1.0f;
    for (const math::Vector &axis : Solids::Octahedron::vertices) {
      const float cos_d = math::dot(p, axis);
      if (cos_d < cos_event_horizon)
        continue;
      const float d = math::fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
      alpha *= math::quintic_kernel(d / MS::EVENT_HORIZON);
    }
    return alpha;
  }
  template <int W, int H>
  static float attractor_hole_alpha(const MindSplatter<W, H> &ms,
                                    const math::Vector &p) {
    return ms.attractor_hole_alpha(
        p, math::fast_cosf(MindSplatter<W, H>::EVENT_HORIZON));
  }
  template <int W, int H>
  static math::Vector attractor_position(const MindSplatter<W, H> &ms,
                                         size_t i) {
    return ms.particle_system.attractors[i].position;
  }
  static math::Vector matrix_vertex(const math::Vector &v,
                                    const math::MobiusParams &mobius,
                                    const math::Quaternion &orientation) {
    math::RotationMatrix rotation(orientation);
    return rotation.apply(math::mobius_transform(v, mobius));
  }
  static math::Vector reference_vertex(const math::Vector &v,
                                       const math::MobiusParams &mobius,
                                       const math::Quaternion &orientation) {
    return math::rotate(math::mobius_transform(v, mobius), orientation);
  }
  static constexpr int trail_length() { return MS::TRAIL_LEN; }
  template <int W, int H>
  static void use_reference_orientation(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_orientation = enabled;
  }
  template <int W, int H>
  static void use_reference_vertex_pass(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_vertex_pass = enabled;
  }
  template <int W, int H>
  static void use_reference_hole_kernel(MindSplatter<W, H> &ms, bool enabled) {
    ms.reference_hole_kernel = enabled;
  }
  /** @brief Drives one choreography-origin preset advance. */
  template <int W, int H> static void advance_preset(MindSplatter<W, H> &ms) {
    HS_CHECK(try_advance_preset(ms),
             "MindSplatter automatic preset advance failed");
  }
  template <int W, int H>
  static bool try_advance_preset(MindSplatter<W, H> &ms) {
    return ms.advance_preset();
  }
  template <int W, int H>
  static bool transition_active(const MindSplatter<W, H> &ms) {
    return ms.transition.active;
  }
  template <int W, int H>
  static void tick_choreography(MindSplatter<W, H> &ms) {
    ms.step_choreography();
  }
  template <int W, int H>
  static void saturate_timeline(MindSplatter<W, H> &ms, float &sink) {
    while (Timeline::remaining() > 0)
      ms.timeline.add(0,
                      Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  }
  template <int W, int H> static void clear_timeline(MindSplatter<W, H> &ms) {
    ms.timeline.clear();
  }
  template <int W, int H>
  static const typename MindSplatter<W, H>::Params &
  live_params(const MindSplatter<W, H> &ms) {
    return ms.params;
  }
  template <int W, int H>
  static const typename MindSplatter<W, H>::Params &
  preset_params(const MindSplatter<W, H> &ms, size_t index) {
    (void)ms;
    return MindSplatter<W, H>::PRESETS[index].params;
  }
  template <int W, int H>
  static size_t preset_index(const MindSplatter<W, H> &ms) {
    return ms.getPresetIndex();
  }
  template <int W, int H> static constexpr size_t preset_count() {
    return MindSplatter<W, H>::PRESETS.size();
  }
  template <int W, int H>
  static size_t preset_count(const MindSplatter<W, H> &) {
    return preset_count<W, H>();
  }
  template <int W, int H>
  static void use_reference_signed_axis_physics(MindSplatter<W, H> &ms,
                                                bool enabled) {
    ms.particle_system.reference_signed_axis_physics = enabled;
  }
  template <int W, int H>
  static uint16_t active_particles(const MindSplatter<W, H> &ms) {
    return ms.particle_system.active();
  }
  template <int W, int H>
  static void fill_particle_capacity(MindSplatter<W, H> &ms) {
    const size_t source_count = ms.particle_system.active();
    HS_CHECK(source_count > 0);
    while (ms.particle_system.active() < ms.particle_system.pool.capacity()) {
      const size_t destination = ms.particle_system.active();
      const auto source = ms.particle_system.pool[destination % source_count];
      ms.particle_system.spawn(math::Vector(), math::Vector(), 0);
      ms.particle_system.pool[destination] = source;
    }
    ms.params.active_count = static_cast<float>(ms.particle_system.active());
  }
  template <int W, int H>
  static uint16_t particle_color_seed(const MindSplatter<W, H> &ms,
                                      size_t index) {
    HS_CHECK(index < ms.particle_system.active());
    return ms.particle_system.pool[index].color_seed;
  }
  template <int W, int H> static void draw_particles(MindSplatter<W, H> &ms) {
    Canvas canvas(ms);
    ms.draw_particles(canvas);
  }
  template <int W, int H>
  static void draw_particles_candidate(MindSplatter<W, H> &ms, Canvas &canvas) {
    ms.draw_particles(canvas);
  }
  /** @brief Enables reference screen stepping, restoring the prior value. */
  struct ReferenceScreenStepScope {
    const bool previous = Plot::g_reference_screen_step;

    ReferenceScreenStepScope() { Plot::g_reference_screen_step = true; }
    ~ReferenceScreenStepScope() { Plot::g_reference_screen_step = previous; }

    ReferenceScreenStepScope(const ReferenceScreenStepScope &) = delete;
    ReferenceScreenStepScope &
    operator=(const ReferenceScreenStepScope &) = delete;
  };

  template <int W, int H>
  static void draw_particles_replay_reference(MindSplatter<W, H> &ms,
                                              Canvas &canvas) {
    const ReferenceScreenStepScope scope;
    ms.draw_particles_with(ms.filters, canvas);
  }
  template <int W, int H>
  static void draw_particles_reference(MindSplatter<W, H> &ms) {
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> sink;
    Canvas canvas(ms);
    ms.draw_particles_with(sink, canvas);
  }
  template <int W, int H>
  static size_t particle_capacity(const MindSplatter<W, H> &ms) {
    return ms.particle_system.pool.capacity();
  }
};

} // namespace effects_tests
} // namespace hs_test
