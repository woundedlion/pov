/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Derivation reach of the composed layer over the operator catalog
// ============================================================================

namespace In = Pullback::Interp;

using ReachBinding = Pullback::ComposedDetail::DiscoveryBinding;

template <typename Family>
concept DerivableSource = requires {
  typename Pullback::SourcePolicyFor<Family, ReachBinding>::Type;
};

// Derivation-side half of the reach table.
static_assert(DerivableSource<Pullback::GridSourceParams>);
// The shared noise field group is not itself a family: it carries the fields
// for both noise sources, so it names no single policy.
static_assert(!DerivableSource<Pullback::Source::NoiseSourceParams>);
static_assert(DerivableSource<Pullback::ProjectedNoiseSourceParams>);
static_assert(DerivableSource<Pullback::SphericalNoiseSourceParams>);
// SphericalRings stays chain-only: its policy needs a prepared axis and
// phase the composed frame does not carry.
static_assert(!DerivableSource<Pullback::Source::SphericalRingsSourceParams>);
static_assert(
    std::is_same_v<
        typename Pullback::SourcePolicyFor<Pullback::SphericalNoiseSourceParams,
                                           ReachBinding>::Type,
        Pullback::Source::SphericalNoise<
            Pullback::SourceProvider<ReachBinding,
                                     Pullback::SphericalNoiseSourceParams>,
            math::NoiseBasis::SIMPLEX>>);
static_assert(
    std::is_same_v<
        typename Pullback::SourcePolicyFor<Pullback::ProjectedNoiseSourceParams,
                                           ReachBinding>::Type,
        Pullback::Source::ProjectedNoise<
            Pullback::SourceProvider<ReachBinding,
                                     Pullback::ProjectedNoiseSourceParams>,
            math::NoiseBasis::SIMPLEX>>);
static_assert(
    std::is_same_v<
        typename Pullback::WarpPolicyFor<Pullback::PolarParams, ReachBinding,
                                         "outer_warp", false>::Type,
        Pullback::Warp::PolarChart<
            Pullback::WarpProvider<ReachBinding, "outer_warp",
                                   Pullback::PolarParams, false>,
            Pullback::Warp::LinearPolar, 1>>);
static_assert(std::is_same_v<
              typename Pullback::WarpPolicyFor<Pullback::VectorNoiseParams,
                                               ReachBinding, "inner_warp",
                                               false>::Type,
              Pullback::Warp::VectorNoise<
                  Pullback::WarpProvider<ReachBinding, "inner_warp",
                                         Pullback::VectorNoiseParams, false>,
                  math::NoiseBasis::SIMPLEX, Pullback::Warp::FlatEnvelope>>);

/**
 * @brief How much of one catalog operator's vocabulary a ComposedEffect
 *        specialization can emit.
 * @details A null `topology_id` marks an operator no Spec or parameter family
 * selects, so none of its topology values are reachable either. Otherwise
 * `reachable` lists exactly the values the derivation layer emits for that
 * enum8 and the rest are workbench-only, so a value the catalog gains counts
 * as unreachable until the row is revisited.
 */
struct DerivationReach {
  const char *operator_id;
  const char *topology_id;
  std::array<const char *, 9>
      reachable; /**< Reachable IDs with nullable unused slots. */
};

constexpr DerivationReach DERIVATION_REACH[] = {
    // ProjectionKind has no spelling for these four.
    {"project.peirce.v3", nullptr, {}},
    {"project.peirce-square-fast.v3", nullptr, {}},
    {"project.bonne.v3", nullptr, {}},
    {"project.airocean.v3", nullptr, {}},
    // SourcePolicyFor has no policy for these samplers; spherical-rings needs
    // a prepared axis and phase.
    {"sample.rings.v2", nullptr, {}},
    {"sample.spherical-rings.v3", nullptr, {}},
    {"sample.fractal.v2", nullptr, {}},
    {"sample.tessellation.v2", nullptr, {}},
    // No WarpPolicyFor specialization carries these warp families.
    {"warp.vortex.v2", nullptr, {}},
    {"warp.curl-flow.v2", nullptr, {}},
    // TransferKind is NONE or ISO_CONTOUR.
    {"field.transfer.ridge.v2", nullptr, {}},
    {"field.transfer.smooth-bands.v2", nullptr, {}},
    // ProjectionKind names the folded gnomonic alone.
    {"project.gnomonic.v2", "hemisphere", {{"folded"}}},
    {"project.stereographic.v2", "frame", {{"identity", "spin-wander"}}},
    {"project.folded-sinusoidal.v2", "frame", {{"identity", "spin-wander"}}},
    {"project.equirectangular.v2", "frame", {{"identity", "spin-wander"}}},
    {"project.gnomonic.v2", "frame", {{"identity", "spin-wander"}}},
    // The displacement policies pin the basis and the integrator.
    {"sphere.displace.curl.v2", "basis", {{"simplex"}}},
    {"sphere.displace.curl.v2", "integrator", {{"euler"}}},
    {"sphere.displace.direct.v2", "basis", {{"simplex"}}},
    // Spec::LensPolicy names one lens policy per symmetry.
    {"sphere.lens.kaleidoscope.v2",
     "symmetry",
     {{"azimuthal", "tetrahedral", "octahedral", "dodecahedral",
       "triangular-prism", "square-prism", "pentagonal-prism",
       "hexagonal-prism", "octagonal-prism"}}},
    // WarpPolicyFor pins the flat envelope, the simplex basis, the linear
    // polar chart and its first harmonic.
    {"warp.wave-shear.v2", "envelope", {{"flat"}}},
    {"warp.vector-noise.v2", "basis", {{"simplex"}}},
    {"warp.vector-noise.v2", "envelope", {{"flat"}}},
    {"warp.polar-chart.v2", "mode", {{"linear"}}},
    {"warp.polar-chart.v2", "harmonic", {{"h1"}}},
    // SampleStage pins Weight::Projection; no shipped Spec selects the
    // none coverage.
    {"sample.grid.v3", "weight-mode", {{"projection"}}},
    {"sample.grid.v3",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.twin-wave.v3", "weight-mode", {{"projection"}}},
    {"sample.twin-wave.v3",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.spiral.v2", "weight-mode", {{"projection"}}},
    {"sample.spiral.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.projected-noise.v2", "weight-mode", {{"projection"}}},
    {"sample.projected-noise.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.projected-noise.v2", "basis", {{"simplex"}}},
    {"sample.spherical-noise.v3", "basis", {{"simplex"}}},
    {"sample.lattice.v2", "weight-mode", {{"projection"}}},
    {"sample.lattice.v2",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    // Shipped Specs select only the none and cup brightness envelopes.
    {"colorize.generated-palette.v3",
     "palette-mode",
     {{"triadic", "complementary", "analogous"}}},
    {"colorize.generated-palette.v3",
     "palette-mapping",
     {{"cup", "bell", "linear", "reverse"}}},
    {"colorize.generated-palette.v3",
     "hue-shift-mode",
     {{"none", "noise", "path-length"}}},
    {"colorize.generated-palette.v3", "brightness-envelope", {{"none", "cup"}}},
};

/** @brief The topology enum8 @p topology_id of @p op, or null. */
inline const In::ParamFieldInfo *find_topology(const In::OperatorDescriptor &op,
                                               const char *topology_id) {
  for (const In::ParamFieldInfo &field : op.schema_span())
    if (field.topology && std::string_view(field.id) == topology_id)
      return &field;
  return nullptr;
}

/** @brief Rows of DERIVATION_REACH naming @p topology_id of @p operator_id. */
inline size_t reach_rows(const char *operator_id, const char *topology_id) {
  size_t rows = 0;
  for (const DerivationReach &row : DERIVATION_REACH)
    if (std::string_view(row.operator_id) == operator_id &&
        row.topology_id != nullptr &&
        std::string_view(row.topology_id) == topology_id)
      ++rows;
  return rows;
}

/** @brief Reports whether derivation emits an operator topology value. */
inline bool derivation_value_reachable(std::string_view operator_id,
                                       std::string_view field_id,
                                       std::string_view value) {
  for (const auto &row : DERIVATION_REACH) {
    if (row.operator_id != operator_id || row.topology_id == nullptr ||
        row.topology_id != field_id)
      continue;
    for (const char *reachable : row.reachable)
      if (reachable != nullptr && reachable == value)
        return true;
  }
  return false;
}

/**
 * @brief Pins the catalog vocabulary no ComposedEffect specialization emits.
 * @details DERIVATION_REACH records how far the derivation layer falls short of
 * the chain interpreter's OPERATOR_TABLE and is resolved against the live
 * table, so a renamed or removed operator, topology enum8 or value reds, as
 * does a catalog addition the table has not classified.
 */
inline void test_composed_derivation_reach() {
  static_assert(AshCloudSpec::FIELD_COVERAGE ==
                Pullback::FieldCoverageKind::VALUE_CUTOUT);
  size_t unreachable_operators = 0;
  size_t unreachable_values = 0;
  for (const DerivationReach &row : DERIVATION_REACH) {
    HS_CONTEXT(row.operator_id);
    const In::OperatorDescriptor *op = In::find_operator(row.operator_id);
    HS_EXPECT_TRUE(op != nullptr);
    if (op == nullptr)
      continue;
    if (row.topology_id == nullptr) {
      ++unreachable_operators;
      for (const In::ParamFieldInfo &field : op->schema_span())
        if (field.topology)
          unreachable_values += field.enum_count;
      continue;
    }
    HS_CONTEXT(row.topology_id);
    const In::ParamFieldInfo *field = find_topology(*op, row.topology_id);
    HS_EXPECT_TRUE(field != nullptr);
    if (field == nullptr)
      continue;
    size_t reachable = 0;
    for (const char *value : row.reachable) {
      if (value == nullptr)
        break;
      HS_CONTEXT(value);
      bool found = false;
      for (uint8_t index = 0; index < field->enum_count; ++index)
        found = found || std::string_view(field->enum_ids[index]) == value;
      HS_EXPECT_TRUE(found);
      reachable += found ? 1 : 0;
    }
    HS_EXPECT_LE(reachable, static_cast<size_t>(field->enum_count));
    unreachable_values += field->enum_count - reachable;
  }

  size_t catalog_values = 0;
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    bool whole_operator = false;
    for (const DerivationReach &row : DERIVATION_REACH)
      whole_operator = whole_operator ||
                       (std::string_view(row.operator_id) == op.operator_id &&
                        row.topology_id == nullptr);
    for (const In::ParamFieldInfo &field : op.schema_span()) {
      if (!field.topology)
        continue;
      HS_CONTEXT(field.id);
      catalog_values += field.enum_count;
      HS_EXPECT_EQ(reach_rows(op.operator_id, field.id),
                   whole_operator ? 0u : 1u);
    }
  }

  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 38u);
  HS_EXPECT_EQ(unreachable_operators, 12u);
  HS_EXPECT_EQ(catalog_values, 150u);
  HS_EXPECT_EQ(unreachable_values, 90u);
}
