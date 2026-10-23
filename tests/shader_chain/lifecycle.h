/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- lifecycle test models ------------------------------------------------

struct CountLifecycle {
  static inline int inits = 0;
  static inline int migrates = 0;
  static inline int destroys = 0;
  static inline bool fail_migrate = false;
  static void reset() {
    inits = 0;
    migrates = 0;
    destroys = 0;
    fail_migrate = false;
  }
};

struct CountParams {
  float gain = 0.5f;
  static constexpr auto FIELDS = std::array{
      PB::Field<CountParams>{"gain", &CountParams::gain, "Gain", 0.0f, 1.0f,
                             PB::FieldCurve::LERP},
  };
};

struct CountingState {
  float accumulator = 0.0f;
  bool live = false;
  CountingState() = default;
  CountingState(const CountingState &) = delete;
  CountingState &operator=(const CountingState &) = delete;
  ~CountingState() {
    if (live)
      ++CountLifecycle::destroys;
  }
};

template <int Tag> struct CountModel {
  static constexpr const char *ID =
      Tag == 0 ? "test.count-a.v2" : "test.count-b.v2";
  static constexpr const char *NAME = "Count";
  using Input = PB::SphereSample;
  using Output = PB::SphereSample;
  using Params = CountParams;
  using State = CountingState;
  struct Prepared {};

  static void init(State &state, In::InstanceId) {
    state.live = true;
    ++CountLifecycle::inits;
  }
  static In::Status migrate(State &dst, const State &src, In::InstanceId) {
    ++CountLifecycle::migrates;
    if (CountLifecycle::fail_migrate)
      return In::Status::FAILED;
    dst.accumulator = src.accumulator;
    dst.live = true;
    return In::Status::OK;
  }
  static void advance(State &state, const Params &) {
    state.accumulator += 1.0f;
  }
  static Prepared prepare(const In::FrameContext &, const Params &,
                          const State &) {
    return {};
  }
  static PB::SphereSample run(const PB::SphereSample &input,
                              const In::FrameContext &, const Params &,
                              const Prepared &) {
    return input;
  }
};

/** @brief Addresses of every constructed LifetimeState, and the runs that
    read a state no longer in the set. */
struct LifetimeRegistry {
  static inline std::vector<const void *> live;
  static inline int dead_reads = 0;
  static void reset() {
    live.clear();
    dead_reads = 0;
  }
  static bool alive(const void *state) {
    for (const void *entry : live)
      if (entry == state)
        return true;
    return false;
  }
  static void remove(const void *state) {
    for (size_t index = 0; index < live.size(); ++index)
      if (live[index] == state) {
        live.erase(live.begin() + static_cast<std::ptrdiff_t>(index));
        return;
      }
  }
};

/** @brief Temporal state that registers its own lifetime. */
struct LifetimeState {
  float phase = 0.0f;
  LifetimeState() { LifetimeRegistry::live.push_back(this); }
  LifetimeState(const LifetimeState &other) : phase(other.phase) {
    LifetimeRegistry::live.push_back(this);
  }
  LifetimeState &operator=(const LifetimeState &) = default;
  ~LifetimeState() { LifetimeRegistry::remove(this); }
};

/** @brief LifetimeState's clock between two endpoints. */
struct LifetimeInterval {
  float from;
  float to;
};

/** @brief Prepares the clock between consecutive endpoints. */
inline LifetimeInterval temporal_interval(const LifetimeState &from,
                                          const LifetimeState &to) {
  return {from.phase, to.phase};
}

/** @brief Samples the clock at @p fraction. */
template <typename Params>
void temporal_sample(LifetimeState &out, const LifetimeInterval &interval,
                     const Params &, float fraction) {
  out.phase = interval.from + (interval.to - interval.from) * fraction;
}

/** @brief Pass-through operator whose prepared block borrows its state and
    whose run counts reads of a destroyed state. */
struct LifetimeModel {
  static constexpr const char *ID = "test.lifetime.v2";
  static constexpr const char *NAME = "Lifetime";
  using Input = PB::SphereSample;
  using Output = PB::SphereSample;
  using Params = CountParams;
  using State = LifetimeState;
  struct Prepared {
    const LifetimeState *state;
  };

  static void init(State &, In::InstanceId) {}
  static In::Status migrate(State &dst, const State &src, In::InstanceId) {
    dst.phase = src.phase;
    return In::Status::OK;
  }
  static void advance(State &state, const Params &) { state.phase += 0.25f; }
  static Prepared prepare(const In::FrameContext &, const Params &,
                          const State &state) {
    return {&state};
  }
  static PB::SphereSample run(const PB::SphereSample &input,
                              const In::FrameContext &, const Params &,
                              const Prepared &prepared) {
    if (!LifetimeRegistry::alive(prepared.state))
      ++LifetimeRegistry::dead_reads;
    return input;
  }
};

inline constexpr auto OUT_OF_RANGE_DEFAULT_SCHEMA = [] {
  auto schema = In::SCHEMA<CountModel<0>>;
  schema[0].def = 1.5f;
  return schema;
}();
static_assert(In::defaults_in_range(In::SCHEMA<CountModel<0>>));
static_assert(!In::defaults_in_range(OUT_OF_RANGE_DEFAULT_SCHEMA));

consteval In::OperatorDescriptor table_entry(std::string_view operator_id) {
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE)
    if (operator_id == op.operator_id)
      return op;
  throw "unknown operator id";
}

template <size_t N>
constexpr std::array<In::ParamFieldInfo, N> filler_schema() {
  std::array<In::ParamFieldInfo, N> out{};
  for (auto &field : out)
    field = In::ParamFieldInfo{
        "fat", nullptr, 0.0f, 1.0f,    0.0f,    PB::FieldCurve::LERP,
        false, 0,       0,    nullptr, nullptr, 0};
  return out;
}

inline constexpr auto FAT_SCHEMA = filler_schema<In::MAX_CHAIN_PARAMS + 1>();
/** With the default chain's project/sample/colorize tail, exactly
    MAX_CHAIN_PARAMS fields. */
inline constexpr auto EXACT_FIT_SCHEMA =
    filler_schema<In::MAX_CHAIN_PARAMS -
                  table_entry("project.stereographic.v2").schema_count -
                  table_entry("sample.grid.v3").schema_count -
                  table_entry("colorize.generated-palette.v3").schema_count>();

/** Never run; exercises only the schema-field budget. */
inline constexpr In::OperatorDescriptor
make_filler_descriptor(const char *operator_id,
                       std::span<const In::ParamFieldInfo> schema) {
  In::OperatorDescriptor descriptor =
      In::make_operator_descriptor<CountModel<0>>();
  descriptor.operator_id = operator_id;
  descriptor.schema = schema.data();
  descriptor.schema_count = static_cast<uint16_t>(schema.size());
  return descriptor;
}

/** @brief A colour chain's operators plus LifetimeModel. */
inline std::span<const In::OperatorDescriptor> lifetime_table() {
  static constexpr std::array<In::OperatorDescriptor, 4> TABLE{
      table_entry("project.stereographic.v2"),
      table_entry("sample.grid.v3"),
      table_entry("colorize.generated-palette.v3"),
      In::make_operator_descriptor<LifetimeModel>(),
  };
  return TABLE;
}

inline std::span<const In::OperatorDescriptor> extended_table() {
  static constexpr std::array<In::OperatorDescriptor, 8> TABLE{
      table_entry("sphere.rotate.v2"),
      table_entry("project.stereographic.v2"),
      table_entry("sample.grid.v3"),
      table_entry("colorize.generated-palette.v3"),
      In::make_operator_descriptor<CountModel<0>>(),
      In::make_operator_descriptor<CountModel<1>>(),
      make_filler_descriptor("test.fat.v2", FAT_SCHEMA),
      make_filler_descriptor("test.exact-fit.v2", EXACT_FIT_SCHEMA),
  };
  return TABLE;
}
