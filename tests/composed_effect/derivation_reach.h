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
// The shared noise field group is not itself a family, so it names no single
// policy.
static_assert(!DerivableSource<Pullback::Source::NoiseSourceParams>);
static_assert(DerivableSource<Pullback::ProjectedNoiseSourceParams>);
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
 * @brief Catalog vocabulary emitted by shipped ComposedEffect specializations.
 * @details A null `topology_id` marks an operator no Spec or parameter family
 * selects, so none of its topology values are reachable either. Otherwise
 * `reachable` lists values emitted by a shipped specialization through a
 * policy pin, Spec choice or authored preset. Selecting another value extends
 * its row and updates the unreachable-count pin.
 */
struct DerivationReach {
  const char *operator_id;
  const char *topology_id;
  std::array<const char *, 9>
      reachable; /**< Reachable IDs with nullable unused slots. */
};

constexpr DerivationReach DERIVATION_REACH[] = {
    {"project.peirce.v3", nullptr, {}},
    {"project.peirce-square-fast.v3", nullptr, {}},
    {"project.bonne.v3", nullptr, {}},
    {"project.airocean.v3", nullptr, {}},
    {"project.tetrahedral.v1", nullptr, {}},
    {"project.cubic.v1", "frame", {{"spin-wander"}}},
    {"project.octahedral.v1", nullptr, {}},
    {"project.dodecahedral.v1", nullptr, {}},
    {"project.icosahedral.v1", nullptr, {}},
    {"warp.curl-flow.v2", nullptr, {}},
    {"sphere.displace.potential.v1", "integrator", {{"euler"}}},
    {"sphere.displace.potential.v1", "term-1-type", {{"mode", "dodecahedral"}}},
    {"sphere.displace.potential.v1", "term-1-kind", {{"toroidal", "poloidal"}}},
    {"sphere.displace.potential.v1", "term-1-l", {{"one", "eight"}}},
    {"sphere.displace.potential.v1", "term-1-m", {{"one", "four"}}},
    {"sphere.displace.potential.v1", "term-2-type", {{"mode"}}},
    {"sphere.displace.potential.v1", "term-2-kind", {{"toroidal"}}},
    {"sphere.displace.potential.v1", "term-2-l", {{"one"}}},
    {"sphere.displace.potential.v1", "term-2-m", {{"zero"}}},
    {"sphere.displace.potential.v1", "term-3-type", {{"mode"}}},
    {"sphere.displace.potential.v1", "term-3-kind", {{"toroidal"}}},
    {"sphere.displace.potential.v1", "term-3-l", {{"one"}}},
    {"sphere.displace.potential.v1", "term-3-m", {{"zero"}}},
    {"sphere.displace.potential.v1", "term-4-type", {{"mode"}}},
    {"sphere.displace.potential.v1", "term-4-kind", {{"toroidal"}}},
    {"sphere.displace.potential.v1", "term-4-l", {{"one"}}},
    {"sphere.displace.potential.v1", "term-4-m", {{"zero"}}},
    {"sphere.displace.rossby-haurwitz.v1", nullptr, {}},
    {"sphere.displace.symmetric.v2", nullptr, {}},
    {"field.transfer.ridge.v2", nullptr, {}},
    {"field.transfer.smooth-bands.v2", nullptr, {}},
    // ProjectionKind names the folded gnomonic alone.
    {"project.gnomonic.v2", "hemisphere", {{"folded"}}},
    {"project.stereographic.v2", "frame", {{"identity", "spin-wander"}}},
    {"project.folded-sinusoidal.v2", "frame", {{"spin-wander"}}},
    {"project.equirectangular.v2", "frame", {{"identity", "spin-wander"}}},
    {"project.gnomonic.v2", "frame", {{"identity", "spin-wander"}}},
    // The displacement policies pin the basis and the integrator.
    {"sphere.displace.curl.v2", "basis", {{"simplex"}}},
    {"sphere.displace.curl.v2", "integrator", {{"euler"}}},
    {"sphere.displace.curl.v2", "limit", {{"smooth"}}},
    {"sphere.displace.direct.v2", "basis", {{"simplex"}}},
    {"sphere.displace.ripple.v2", "burst-size", {{"four"}}},
    {"sphere.displace.ripple.v2", "random-center", {{"on"}}},
    {"sphere.lens.kaleidoscope.v2",
     "symmetry",
     {{"azimuthal", "tetrahedral", "octahedral", "dodecahedral",
       "hexagonal-prism"}}},
    // The affine TrackPath flag pins its path length.
    {"warp.affine.v4", "path-length", {{"off", "on"}}},
    // WarpPolicyFor pins the flat envelope and simplex basis.
    {"warp.wave-shear.v2", "envelope", {{"flat"}}},
    {"warp.vector-noise.v2", "basis", {{"simplex"}}},
    {"warp.vector-noise.v2", "envelope", {{"flat"}}},
    {"warp.polar-chart.v2", nullptr, {}},
    // SampleStage pins Weight::Projection; no shipped Spec selects the
    // none coverage.
    {"sample.grid.v4", "weight-mode", {{"projection"}}},
    {"sample.grid.v4",
     "coverage-mode",
     {{"weight", "weight-squared", "edge-fade"}}},
    {"sample.twin-wave.v3", "weight-mode", {{"projection"}}},
    {"sample.twin-wave.v3", "coverage-mode", {{"weight-squared"}}},
    {"sample.bands.v1", nullptr, {}},
    {"sample.spiral.v3", "weight-mode", {{"projection"}}},
    {"sample.spiral.v3", "coverage-mode", {{"weight-squared"}}},
    {"sample.projected-noise.v2", "weight-mode", {}},
    {"sample.projected-noise.v2", "coverage-mode", {}},
    {"sample.projected-noise.v2", "basis", {}},
    {"sample.tessellation.v3", "weight-mode", {{"projection"}}},
    {"sample.tessellation.v3", "coverage-mode", {{"weight", "edge-fade"}}},
    // Shipped Specs pin square cells.
    {"sample.tessellation.v3", "kind", {{"square"}}},
    {"sample.tessellation.v3", "fill", {{"cells"}}},
    {"colorize.generated-palette.v4",
     "palette-mode",
     {{"triadic", "complementary", "analogous"}}},
    {"colorize.generated-palette.v4",
     "palette-type",
     {{"cycling-generated", "generated", "procedural"}}},
    {"colorize.generated-palette.v4", "palette-harmony", {{"triadic"}}},
    {"colorize.generated-palette.v4",
     "palette-lightness-curve",
     {{"cup", "constant"}}},
    {"colorize.generated-palette.v4", "palette-name", {{"embers"}}},
    {"colorize.generated-palette.v4", "palette-name-b", {{"fire-and-ice"}}},
    {"colorize.generated-palette.v4",
     "palette-mapping",
     {{"cup", "bell", "linear"}}},
    {"color.modulate.v1", "modulation-enabled", {{"enabled"}}},
    {"color.modulate.v1",
     "modulation-source",
     {{"noise", "displacement", "sphere-smooth-grid", "directional-sine"}}},
    {"color.modulate.v1", "modulation-evaluation", {{"direct", "baked-cube"}}},
    {"color.modulate.v1",
     "modulation-mechanic",
     {{"hue-rotation", "positive-lightness-lift"}}},
    {"color.modulate.v1", "pattern-noise-basis", {{"simplex"}}},
    {"color.modulate.v1", "pattern-noise-variant", {{"configured-single"}}},
    // Shipped Specs select the none, cup and bell brightness envelopes.
    {"colorize.generated-palette.v4",
     "brightness-envelope",
     {{"none", "cup", "bell"}}},
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
 * @details DERIVATION_REACH classifies shipped policy pins and Spec choices
 * against the live catalog. A Spec selecting another value extends its row
 * and updates the unreachable-count pin.
 */
inline void test_composed_derivation_reach() {
  static_assert(Pullback::MetadataFor<AshCloudSpec>::FIELD_COVERAGE ==
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

  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 43u);
  HS_EXPECT_EQ(unreachable_operators, 15u);
  HS_EXPECT_EQ(catalog_values, 372u);
  HS_EXPECT_EQ(unreachable_values, 285u);
}
