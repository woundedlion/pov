/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Authored fixed-pipeline policy death cases.

/** @brief Frame binding for one keyed authored stage.
 * @tparam Family Parameter family.
 * @tparam Kind Resource kind.
 * @tparam Key Resource key. */
template <typename Family, Pullback::ResourceKind Kind,
          Pullback::ResourceKey Key = "probe">
struct AuthoredStageFrame {
  using Params = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<
          Pullback::ParameterResource<Key, Family, Kind>>>;
  using Frame = Pullback::FrameIntervalState<Params>;
  using Binding = Pullback::Binding<Frame>;
};

/** @brief Prepares an authored warp after corrupting one topology selector.
 * @tparam Family Authored warp family.
 * @tparam Corrupt Selector writer.
 * @param corrupt Writes the invalid selector. */
template <typename Family, typename Corrupt>
inline void authored_warp_prepare(Corrupt corrupt) {
  using Fixture = AuthoredStageFrame<Family, Pullback::ResourceKind::WARP>;
  using Policy =
      typename Pullback::WarpPolicyFor<Family, typename Fixture::Binding,
                                       "probe", true>::Type;
  typename Fixture::Frame frame{};
  corrupt(frame.params.template get<"probe">());
  (void)Policy::prepare(frame);
}

/** @brief Prepares an authored source after corrupting one topology selector.
 * @tparam Family Authored source family.
 * @tparam Corrupt Selector writer.
 * @param corrupt Writes the invalid selector. */
template <typename Family, typename Corrupt>
inline void authored_source_prepare(Corrupt corrupt) {
  using Fixture =
      AuthoredStageFrame<Family, Pullback::ResourceKind::SOURCE, "source">;
  using Policy =
      typename Pullback::SourcePolicyFor<Family,
                                         typename Fixture::Binding>::Type;
  typename Fixture::Frame frame{};
  corrupt(frame.params.template get<"source">());
  (void)Policy::prepare(frame);
}

/** @brief Death case: authored wave shear rejects an unknown envelope. */
inline void case_authored_wave_shear_invalid_envelope() {
  authored_warp_prepare<Pullback::Warp::AuthoredWaveShearParams>(
      [](auto &params) { params.envelope = opaque<uint8_t>(0xff); });
}

/** @brief Death case: authored vector noise rejects an unknown basis. */
inline void case_authored_vector_noise_invalid_basis() {
  authored_warp_prepare<Pullback::Warp::AuthoredVectorNoiseParams>(
      [](auto &params) { params.basis = opaque<uint8_t>(0xff); });
}

/** @brief Death case: an authored polar chart rejects an unknown mode. */
inline void case_authored_polar_chart_invalid_mode() {
  authored_warp_prepare<Pullback::Warp::AuthoredPolarParams>(
      [](auto &params) { params.mode = opaque<uint8_t>(0xff); });
}

/** @brief Death case: authored curl flow rejects an unknown integrator. */
inline void case_authored_curl_flow_invalid_integrator() {
  authored_warp_prepare<Pullback::Warp::AuthoredCurlFlowParams>(
      [](auto &params) { params.integrator = opaque<uint8_t>(0xff); });
}

/** @brief Death case: authored projected noise rejects an unknown basis. */
inline void case_authored_projected_noise_invalid_basis() {
  authored_source_prepare<Pullback::Source::AuthoredProjectedNoiseSourceParams>(
      [](auto &params) { params.basis = opaque<uint8_t>(0xff); });
}

/** @brief Applies an authored sample crossing after corrupting one mode.
 * @tparam Policy Crossing policy over a source provider.
 * @tparam Corrupt Mode writer.
 * @param corrupt Writes the invalid mode. */
template <template <typename> class Policy, typename Corrupt>
inline void authored_sample_apply(Corrupt corrupt) {
  using Family = Pullback::Source::AuthoredRingsSampleParams;
  using Fixture =
      AuthoredStageFrame<Family, Pullback::ResourceKind::SOURCE, "source">;
  using Provider = Pullback::SourceProvider<typename Fixture::Binding, Family>;
  typename Fixture::Frame frame{};
  corrupt(frame.params.template get<"source">());
  Pullback::ProjectionProvenance provenance{};
  float result;
  if constexpr (requires { Policy<Provider>::apply(provenance, frame); })
    result = Policy<Provider>::apply(provenance, frame);
  else
    result = Policy<Provider>::apply(1.0f, provenance, frame);
  if (result != 0.0f)
    std::printf("x");
}

/** @brief Death case: authored signal weighting rejects an unknown mode. */
inline void case_authored_sample_invalid_weight_mode() {
  authored_sample_apply<Pullback::Weight::Authored>(
      [](auto &params) { params.weight_mode = opaque<uint8_t>(0xff); });
}

/** @brief Death case: authored projection coverage rejects an unknown mode. */
inline void case_authored_sample_invalid_coverage_mode() {
  authored_sample_apply<Pullback::ProjectionCoverage::Authored>(
      [](auto &params) { params.coverage_mode = opaque<uint8_t>(0xff); });
}
