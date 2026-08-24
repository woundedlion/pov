/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "effects/MindSplatter.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace effects_tests {

/** @brief White-box accessor shared by MindSplatter tests and replay tools. */
struct MindSplatterWhiteBox {
  using MS = MindSplatter<288, 144>;

  static constexpr uint32_t REPLAY_MAGIC = 0x3152534du;
  static constexpr uint16_t REPLAY_VERSION = 3;

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
    ObjectBytes<Orientation<>> orientation{};
    ObjectBytes<MobiusParams> mobius{};
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
    HS_CHECK(sizeof(Particle) <= UINT16_MAX);
    HS_CHECK(sizeof(Orientation<>) <= UINT16_MAX);
    HS_CHECK(sizeof(MobiusParams) <= UINT16_MAX);

    std::vector<unsigned char> bytes;
    bytes.reserve(24 + sizeof(Orientation<>) + sizeof(MobiusParams) +
                  snapshot.particles.size() * sizeof(Particle));
    append_u32(bytes, REPLAY_MAGIC);
    append_u16(bytes, REPLAY_VERSION);
    append_u16(bytes, W);
    append_u16(bytes, H);
    append_u16(bytes, static_cast<uint16_t>(sizeof(Particle)));
    append_u16(bytes, static_cast<uint16_t>(sizeof(Orientation<>)));
    append_u16(bytes, static_cast<uint16_t>(sizeof(MobiusParams)));
    append_u16(bytes, static_cast<uint16_t>(snapshot.particles.size()));
    append_u16(bytes, snapshot.max_life);
    bytes.push_back(static_cast<unsigned char>(snapshot.active_base_mesh));
    append_object<Orientation<>>(bytes, snapshot.orientation);
    append_object<MobiusParams>(bytes, snapshot.mobius);
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
      ms.particle_system.spawn(Vector(), Vector(), 0);
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
    HS_CHECK(reader.read_u16() == sizeof(Orientation<>),
             "MindSplatter replay orientation layout differs");
    HS_CHECK(reader.read_u16() == sizeof(MobiusParams),
             "MindSplatter replay Mobius layout differs");
    const uint16_t particle_count = reader.read_u16();
    const uint16_t max_life = reader.read_u16();
    const auto active_base_mesh =
        static_cast<typename Snapshot::BaseMesh>(reader.read_u8());
    HS_CHECK(static_cast<size_t>(active_base_mesh) <
                 Solids::PLATONIC_BASE_MESH_COUNT,
             "MindSplatter replay base mesh differs");
    HS_CHECK(ms.particle_system.active() == 0);
    HS_CHECK(particle_count <= ms.particle_system.pool.capacity());

    ms.configure_particle_geometry(active_base_mesh);
    restore_object(ms.orientation, reader.read_object<Orientation<>>());
    restore_object(ms.mobius, reader.read_object<MobiusParams>());
    ms.particle_system.max_life = max_life;
    for (uint16_t i = 0; i < particle_count; ++i) {
      ms.particle_system.spawn(Vector(), Vector(), 0);
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
    ms.timeline.step(canvas);
    if (!ms.particle_geometry_ready ||
        ms.params.base_mesh != ms.active_base_mesh)
      ms.configure_particle_geometry(ms.params.base_mesh);
    ms.particle_system.friction = ms.params.friction;
    for (size_t i = 0; i < ms.particle_system.attractors.size(); ++i)
      ms.particle_system.attractors[i].strength = ms.params.well_strength;
    ms.particle_system.step(canvas);
    ms.params.active_count = static_cast<float>(ms.particle_system.active());
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
  static float hole_alpha(const Vector &p) {
    return MS::octahedral_hole_alpha(p, fast_cosf(MS::EVENT_HORIZON));
  }
  static float reference_hole_alpha(const Vector &p) {
    const float cos_event_horizon = fast_cosf(MS::EVENT_HORIZON);
    float alpha = 1.0f;
    for (const Vector &axis : Solids::Octahedron::vertices) {
      const float cos_d = dot(p, axis);
      if (cos_d < cos_event_horizon)
        continue;
      const float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
      alpha *= quintic_kernel(d / MS::EVENT_HORIZON);
    }
    return alpha;
  }
  static Vector matrix_vertex(const Vector &v, const MobiusParams &mobius,
                              const Quaternion &orientation) {
    RotationMatrix rotation(orientation);
    return rotation.apply(mobius_transform(v, mobius));
  }
  static Vector reference_vertex(const Vector &v, const MobiusParams &mobius,
                                 const Quaternion &orientation) {
    return rotate(mobius_transform(v, mobius), orientation);
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
  template <int W, int H>
  static auto trail_palette(const MindSplatter<W, H> &ms, uint16_t seed) {
    (void)ms;
    std::array<Pixel, MINDSPLATTER_PALETTE_LUT_SIZE> colors;
    for (int i = 0; i < MINDSPLATTER_PALETTE_LUT_SIZE; ++i)
      colors[i] = MINDSPLATTER_PALETTES[seed >> 8][i];
    return colors;
  }
  /** @brief Drives one choreography-origin preset advance. */
  template <int W, int H> static void advance_preset(MindSplatter<W, H> &ms) {
    HS_CHECK(ms.advancePreset(),
             "MindSplatter automatic preset advance failed");
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
  template <int W, int H>
  static size_t preset_count(const MindSplatter<W, H> &ms) {
    return ms.PRESETS.size();
  }
  template <int W, int H>
  static void use_reference_signed_axis_physics(MindSplatter<W, H> &ms,
                                                bool enabled) {
    ms.particle_system.reference_signed_axis_physics = enabled;
  }
  template <int W, int H>
  static void use_clip_clear(MindSplatter<W, H> &ms, bool enabled) {
    ms.force_full_buffer_clear = !enabled;
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
      ms.particle_system.spawn(Vector(), Vector(), 0);
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
  template <int W, int H>
  static void draw_particles(MindSplatter<W, H> &ms, float opacity = 1.0f) {
    Canvas canvas(ms);
    ms.draw_particles(canvas, opacity);
  }
  template <int W, int H, typename Inspect>
  static void draw_particles_inspect(MindSplatter<W, H> &ms, Inspect inspect) {
    Canvas canvas(ms);
    ms.draw_particles(canvas);
    inspect(canvas, ms.clip());
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
