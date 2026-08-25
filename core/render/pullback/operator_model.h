/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <array>
#include <cstdint>
#include <new>
#include <span>
#include <string_view>
#include <type_traits>

#include "render/pullback/color.h"
#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file operator_model.h
 * @brief The chain-interpreter operator model: the erased runtime ABI, the
 *        catalog-facing parameter schema, and the factory that derives one
 *        OperatorDescriptor from a typed operator model.
 */

namespace Pullback {

namespace Interp {

/** @brief Carrier identity of a chain endpoint; the CarrierList rank. */
enum class CarrierId : uint8_t {
  SPHERE = static_cast<uint8_t>(family_of<SphereSample>),
  PLANE = static_cast<uint8_t>(family_of<PlaneSample>),
  FIELD = static_cast<uint8_t>(family_of<FieldSample>),
  COLOR = static_cast<uint8_t>(family_of<Color4>),
  COUNT = static_cast<uint8_t>(CarrierList::SIZE)
};

/** @brief Result of a fallible lifecycle callback. */
enum class Status : uint8_t { OK, FAILED };

/** @brief CarrierId of a canonical carrier type. */
template <CanonicalCarrier T> consteval CarrierId carrier_id_of() {
  return static_cast<CarrierId>(family_of<T>);
}

/** @brief Catalog spelling of each carrier, indexed by CarrierId. */
inline constexpr auto CARRIER_NAMES =
    std::to_array<const char *>({"sphere", "plane", "field", "color"});
static_assert(CARRIER_NAMES.size() == CarrierList::SIZE);
static_assert(CARRIER_NAMES.size() == static_cast<size_t>(CarrierId::COUNT));

/**
 * @brief The (instance_id, operator_id) pair identity of one chain entry.
 * @details `stable_hash` seeds instance-owned resources, so a fresh init of
 * the same pair reconstructs the same resources deterministically.
 */
struct InstanceId {
  const char *instance;
  const char *operator_id;
  uint32_t stable_hash;
};

/** @brief FNV-1a over the instance label, a separator, and the operator id. */
constexpr uint32_t instance_hash(std::string_view instance,
                                 std::string_view operator_id) {
  uint32_t hash = 2166136261u;
  const auto mix = [&hash](char c) {
    hash ^= static_cast<uint8_t>(c);
    hash *= 16777619u;
  };
  for (char c : instance)
    mix(c);
  mix('\0');
  for (char c : operator_id)
    mix(c);
  return hash;
}

/**
 * @brief The interpreter's per-frame snapshot: shared engine-owned resources
 *        and transforms the operators read.
 * @details Pointers alias engine-owned storage and are valid only within the
 * draw_frame() that built the context — the same lifetime contract as
 * ComposedEffect's FrameState.
 */
struct FrameContext {
  uint32_t frame = 0;
  float time = 0.0f;
  /** Base orientation composed under the projection spin/wander frame. */
  Quaternion projection_base;
  /** Generated palette bakes, indexed by the colorize palette-mode value. */
  std::array<const BakedPalette *, 3> palettes{};
  const Pixel *hue_rotation_lut = nullptr;
  const int8_t *hue_noise_lut = nullptr;
};

/**
 * @brief One topology enum8 of a parameter family: a structural switch that
 *        selects among eagerly constructed variants at run time.
 * @tparam Owner The family struct the enum8 lives in.
 */
template <typename Owner> struct TopologyField {
  const char *id; /**< Stable machine id, kebab-case. */
  uint8_t Owner::*member;
  const char *const *value_ids; /**< Kebab-case value spellings. */
  uint8_t value_count;
  uint8_t def; /**< Default value index. */
};

/**
 * @brief One entry of an operator's registered parameter schema — the full
 *        union across topology variants, topology enum8s included.
 */
struct ParamFieldInfo {
  const char *id;
  const char *display_name; /**< Nullable; cosmetic. */
  float min;
  float max;
  float def;
  FieldCurve curve;
  bool topology;
  uint8_t enum_count; /**< 0 marks a float field. */
  uint8_t enum_def;
  const char *const *enum_ids;
  const char *const *enum_names; /**< Nullable; cosmetic. */
};

/** @brief Size and alignment of one arena-allocated block. */
struct BlockLayout {
  uint32_t size;
  uint32_t align;
};

/**
 * @brief The erased runtime ABI of one operator-table entry.
 * @details Blocks are arena-allocated by chain compilation from the declared
 * layouts. `init` is infallible and constructs owned resources eagerly for
 * every topology variant; `migrate` clones `src` into unconstructed `dst`
 * storage (on failure `dst` is left unconstructed and `src` untouched);
 * `advance` steps per-frame clocks; `prepare` derives the frame's prepared
 * block from context, params and state; `run` placement-constructs the output
 * carrier, licensed by validated adjacency.
 */
struct OperatorRuntime {
  BlockLayout param;
  BlockLayout prepared;
  BlockLayout state;
  void (*construct_params)(void *params);
  void (*init)(void *state, InstanceId id);
  Status (*migrate)(void *dst, const void *src, InstanceId id);
  void (*destroy)(void *state);
  void (*advance)(void *state, const uint8_t *params);
  void (*prepare)(const FrameContext &ctx, const uint8_t *params,
                  const void *state, uint8_t *prepared);
  void (*run)(const void *in, void *out, const FrameContext &ctx,
              const uint8_t *params, const uint8_t *prepared);
  /** Address of one schema field inside a param block, by schema index. */
  void *(*param_address)(void *params, uint16_t schema_index);
};

/**
 * @brief One operator's record: the single authority the interpreter table,
 *        the tool catalog and (later) promotion all read.
 * @details Always produced by make_operator_descriptor(), never aggregate
 * literals at call sites, so promotion fields can be added additively.
 */
struct OperatorDescriptor {
  const char *operator_id;
  const char *display_name;
  CarrierId input;
  CarrierId output;
  const ParamFieldInfo *schema;
  uint16_t schema_count;
  OperatorRuntime runtime;
  bool approximate;
  ApproximationOracleId oracle;
  const ApproximationMetric *metrics;
  uint8_t metric_count;

  std::span<const ParamFieldInfo> schema_span() const {
    return {schema, schema_count};
  }
};

/**
 * @brief Rebuilds a base family's field table over a derived family.
 * @details A `float Base::*` converts implicitly to `float Combined::*`, so
 * the concatenation stays a homogeneous Field<Combined> array.
 */
template <typename Combined, typename BaseOwner, size_t N, size_t M>
consteval std::array<Field<Combined>, N + M>
concat_fields(const std::array<Field<BaseOwner>, N> &base,
              const std::array<Field<Combined>, M> &extra) {
  std::array<Field<Combined>, N + M> out{};
  for (size_t index = 0; index < N; ++index)
    out[index] = Field<Combined>{
        base[index].id,  base[index].member, base[index].name, base[index].min,
        base[index].max, base[index].curve,  base[index].gate};
  for (size_t index = 0; index < M; ++index)
    out[N + index] = extra[index];
  return out;
}

namespace Detail {

/** @brief The family's TOPOLOGY table, or an empty one when it declares none. */
template <typename Params> consteval auto topology_of() {
  if constexpr (requires { Params::TOPOLOGY; })
    return Params::TOPOLOGY;
  else
    return std::array<TopologyField<Params>, 0>{};
}

} // namespace Detail

/** @brief Whether topology catalog defaults match the parameter defaults. */
template <typename Params> consteval bool topology_defaults_match() {
  constexpr auto TOPOLOGY = Detail::topology_of<Params>();
  constexpr Params DEFAULTS{};
  for (const auto &field : TOPOLOGY)
    if (DEFAULTS.*(field.member) != field.def)
      return false;
  return true;
}

namespace Detail {

template <typename Model> consteval auto make_schema() {
  using Params = typename Model::Params;
  static_assert(topology_defaults_match<Params>());
  constexpr auto TOPOLOGY = topology_of<Params>();
  constexpr size_t FIELD_COUNT = Params::FIELDS.size();
  std::array<ParamFieldInfo, FIELD_COUNT + TOPOLOGY.size()> out{};
  constexpr Params DEFAULTS{};
  for (size_t index = 0; index < FIELD_COUNT; ++index) {
    const auto &field = Params::FIELDS[index];
    out[index] = ParamFieldInfo{
        field.id,    field.name, field.min, field.max, DEFAULTS.*(field.member),
        field.curve, false,      0,         0,         nullptr,
        nullptr};
  }
  for (size_t index = 0; index < TOPOLOGY.size(); ++index) {
    const auto &topo = TOPOLOGY[index];
    out[FIELD_COUNT + index] =
        ParamFieldInfo{topo.id,  nullptr,          0.0f,   0.0f,
                       0.0f,     FieldCurve::SNAP, true,   topo.value_count,
                       topo.def, topo.value_ids,   nullptr};
  }
  return out;
}

} // namespace Detail

/** @brief Whether every schema id is non-null and unique across the schema. */
template <size_t N>
consteval bool schema_ids_unique(const std::array<ParamFieldInfo, N> &schema) {
  for (size_t i = 0; i < N; ++i) {
    if (schema[i].id == nullptr)
      return false;
    for (size_t j = 0; j < i; ++j)
      if (std::string_view(schema[i].id) == schema[j].id)
        return false;
  }
  return true;
}

/** @brief Every topology enum8 is a genuine switch with well-formed values. */
template <size_t N>
consteval bool
topology_wellformed(const std::array<ParamFieldInfo, N> &schema) {
  for (const ParamFieldInfo &field : schema) {
    if (!field.topology)
      continue;
    if (field.enum_count < 2 || field.enum_ids == nullptr ||
        field.enum_def >= field.enum_count)
      return false;
    for (uint8_t value = 0; value < field.enum_count; ++value)
      if (field.enum_ids[value] == nullptr)
        return false;
  }
  return true;
}

/** @brief Whether every float default lies within its declared range. */
template <size_t N>
consteval bool defaults_in_range(const std::array<ParamFieldInfo, N> &schema) {
  for (const ParamFieldInfo &field : schema)
    if (!field.topology && !(field.def >= field.min && field.def <= field.max))
      return false;
  return true;
}

/** @brief The registered schema of @p Model: family fields, then topology. */
template <typename Model>
inline constexpr auto SCHEMA = Detail::make_schema<Model>();

namespace Detail {

/** Typed-to-erased trampolines over one operator model. */
template <typename Model> struct ErasedAdapter {
  using Params = typename Model::Params;
  using State = typename Model::State;
  using Prepared = typename Model::Prepared;
  using Input = typename Model::Input;
  using Output = typename Model::Output;

  static void construct_params(void *params) { ::new (params) Params{}; }

  static void init(void *state, InstanceId id) {
    Model::init(*::new (state) State{}, id);
  }

  static Status migrate(void *dst, const void *src, InstanceId id) {
    State &clone = *::new (dst) State{};
    const Status status =
        Model::migrate(clone, *static_cast<const State *>(src), id);
    if (status != Status::OK)
      clone.~State();
    return status;
  }

  static void destroy(void *state) { static_cast<State *>(state)->~State(); }

  static void advance(void *state, const uint8_t *params) {
    Model::advance(*std::launder(static_cast<State *>(state)),
                   *std::launder(reinterpret_cast<const Params *>(params)));
  }

  static void prepare(const FrameContext &ctx, const uint8_t *params,
                      const void *state, uint8_t *prepared) {
    ::new (static_cast<void *>(prepared)) Prepared{Model::prepare(
        ctx, *std::launder(reinterpret_cast<const Params *>(params)),
        *std::launder(static_cast<const State *>(state)))};
  }

  static void run(const void *in, void *out, const FrameContext &ctx,
                  const uint8_t *params, const uint8_t *prepared) {
    ::new (out) Output{Model::run(
        *std::launder(static_cast<const Input *>(in)), ctx,
        *std::launder(reinterpret_cast<const Params *>(params)),
        *std::launder(reinterpret_cast<const Prepared *>(prepared)))};
  }

  static void *param_address(void *params, uint16_t schema_index) {
    Params &block = *static_cast<Params *>(params);
    constexpr size_t FIELD_COUNT = Params::FIELDS.size();
    if (schema_index < FIELD_COUNT)
      return &(block.*(Params::FIELDS[schema_index].member));
    constexpr auto TOPOLOGY = topology_of<Params>();
    const size_t topo_index = schema_index - FIELD_COUNT;
    HS_CHECK(topo_index < TOPOLOGY.size(),
             "operator param_address: schema index out of range");
    return &(block.*(TOPOLOGY[topo_index].member));
  }
};

template <typename Model> consteval bool model_approximate() {
  if constexpr (requires { Model::APPROXIMATE; })
    return Model::APPROXIMATE;
  else
    return false;
}

template <typename Model> consteval ApproximationOracleId model_oracle() {
  if constexpr (requires { Model::ORACLE; })
    return Model::ORACLE;
  else
    return ApproximationOracleId::NONE;
}

template <typename Model> constexpr const ApproximationMetric *model_metrics() {
  if constexpr (requires { Model::METRICS; })
    return Model::METRICS.data();
  else
    return nullptr;
}

template <typename Model> consteval uint8_t model_metric_count() {
  if constexpr (requires { Model::METRICS; })
    return static_cast<uint8_t>(Model::METRICS.size());
  else
    return 0;
}

template <typename Model> consteval bool model_non_floating_fields_exact() {
  if constexpr (requires { Model::NON_FLOATING_FIELDS_EXACT; })
    return Model::NON_FLOATING_FIELDS_EXACT;
  else
    return true;
}

template <typename Model> consteval bool model_final_framebuffer_metric() {
  if constexpr (requires { Model::METRICS; })
    return has_final_framebuffer_metric<Model>();
  else
    return false;
}

} // namespace Detail

/** @brief Instance state of a model that owns nothing. */
struct EmptyState {};

/**
 * @brief State and defaulted lifecycle for a model whose instance state is a
 *        plain value.
 * @details `init` leaves the default-constructed state; `migrate` clones by
 * assignment, so `dst` is fully constructed and `src` untouched. A model that
 * seeds instance-owned resources declares its own `init`.
 * @tparam S The model's State type.
 */
template <typename S> struct ValueStateModel {
  using State = S;
  static void init(State &, InstanceId) {}
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
};

/** @brief State and defaulted lifecycle for a model with no instance state. */
struct StatelessModel : ValueStateModel<EmptyState> {
  template <typename Params> static void advance(State &, const Params &) {}
};

/**
 * @brief Derives one OperatorDescriptor from a typed operator model.
 * @details Every derivable field comes from the member that owns it: carriers
 * from Input/Output, block layouts from sizeof/alignof, the schema from the
 * family FIELDS and TOPOLOGY tables, the callbacks as typed-to-erased
 * trampolines, and approximation metadata from the model's declarations.
 * @tparam Model The operator model; see the models in operators.h.
 */
template <typename Model>
constexpr OperatorDescriptor make_operator_descriptor() {
  using Params = typename Model::Params;
  using State = typename Model::State;
  using Prepared = typename Model::Prepared;
  static_assert(CanonicalCarrier<typename Model::Input> &&
                    CanonicalCarrier<typename Model::Output>,
                "operator model: carriers must be canonical");
  static_assert(family_of<typename Model::Input> <=
                    family_of<typename Model::Output>,
                "operator model: family rank may not decrease");
  static_assert(std::is_trivially_destructible_v<Params> &&
                    std::is_trivially_destructible_v<Prepared>,
                "operator model: param and prepared blocks must be trivially "
                "destructible");
  static_assert(std::is_default_constructible_v<State>,
                "operator model: state must be default-constructible");
  static_assert(alignof(Params) <= alignof(std::max_align_t) &&
                    alignof(Prepared) <= alignof(std::max_align_t) &&
                    alignof(State) <= alignof(std::max_align_t),
                "operator model: block alignment exceeds the arena guarantee");
  static_assert(field_ids_unique<Params>());
  static_assert(schema_ids_unique(SCHEMA<Model>));
  static_assert(topology_wellformed(SCHEMA<Model>));
  static_assert(defaults_in_range(SCHEMA<Model>),
                "operator model: parameter default outside declared range");
  static_assert(
      Detail::model_approximate<Model>()
          ? (Detail::model_oracle<Model>() != ApproximationOracleId::NONE &&
             Detail::model_metric_count<Model>() > 0 &&
             Detail::model_non_floating_fields_exact<Model>() &&
             Detail::model_final_framebuffer_metric<Model>())
          : (Detail::model_oracle<Model>() == ApproximationOracleId::NONE &&
             Detail::model_metric_count<Model>() == 0),
      "operator model: malformed approximation metadata — an approximate "
      "operator declares an oracle, exact non-floating fields, and metrics "
      "including a final framebuffer bound; an exact one declares none");
  using Adapter = Detail::ErasedAdapter<Model>;
  return OperatorDescriptor{
      Model::ID,
      Model::NAME,
      carrier_id_of<typename Model::Input>(),
      carrier_id_of<typename Model::Output>(),
      SCHEMA<Model>.data(),
      static_cast<uint16_t>(SCHEMA<Model>.size()),
      OperatorRuntime{
          {sizeof(Params), alignof(Params)},
          {sizeof(Prepared), alignof(Prepared)},
          {sizeof(State), alignof(State)},
          &Adapter::construct_params,
          &Adapter::init,
          &Adapter::migrate,
          &Adapter::destroy,
          &Adapter::advance,
          &Adapter::prepare,
          &Adapter::run,
          &Adapter::param_address,
      },
      Detail::model_approximate<Model>(),
      Detail::model_oracle<Model>(),
      Detail::model_metrics<Model>(),
      Detail::model_metric_count<Model>(),
  };
}

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
