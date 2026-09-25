/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "animation.h"
#include "mesh/recipe.h"

namespace Animation {

/**
 * @brief Schedules swept recipe legs and carries mesh palettes across compaction.
 * @tparam Host Owner providing timeline, carousel, palette storage, drawing and
 * persistent-storage reclamation.
 * @tparam MAX_BUILD_STEPS Maximum lowered primitive count.
 * @tparam MAX_BUILD_FACES Maximum mesh face count.
 */
template <class Host, size_t MAX_BUILD_STEPS, size_t MAX_BUILD_FACES>
class RecipeBuild {
  Host &host() { return static_cast<Host &>(*this); }
  const Host &host() const { return static_cast<const Host &>(*this); }

protected:
  // Recipe-build leg lengths; divided by the Trans Speed divisor like every
  // other stage.
  static constexpr int HANKIN_LEG_FRAMES = 32;
  static constexpr int SWEEP_LEG_FRAMES =
      24; /**< ambo / truncate / snub / chamfer. */
  static constexpr int RELAX_LEG_FRAMES = 16;
  static constexpr int RECONCILE_LEG_FRAMES =
      24; /**< identity-mesh -> authored kis/needle slerp. */
  /** Identity-mesh truncate depth of the smooth kis/needle path: the "uniform"
   * Conway depth at which dual(truncate(X)) matches kis(dual(X)) exactly on
   * regular seeds (docs/specs/opchain_morph_spec.md, smooth kis/needle). */
  static constexpr float MACRO_TRUNCATE_T = 1.0f / 3.0f;
  // Build-chain state (entries with a non-null recipe): the shape is built op
  // by op between the fade-in and the still hold. Null-recipe entries never
  // touch any of it.
  bool build_active = false; /**< Legs draw; the sprite draw_fn is muted. */
  /** Entry the resident shape was spawned from; must outlive its build (the
   * registries are static). The build chain reads its recipe rather than
   * re-indexing a registry, so a spawn is not tied to a registry slot. */
  const Solids::Entry *build_entry = nullptr;
  Solids::OpStep
      build_step_chain[MAX_BUILD_STEPS];      /**< Lowered primitive chain. */
  size_t build_step_count = 0;                /**< Lowered step count. */
  size_t build_step = 0;                      /**< Current leg index. */
  int build_leg_frames[MAX_BUILD_STEPS] = {}; /**< Per-leg frame budget. */
  int build_total_frames = 0;                 /**< Sum of leg frames. */
  PolyMesh build_seed;                        /**< Leg-k seed (persistent). */
  PolyMesh build_next_seed;  /**< Clean endpoint seed_{k+1}: built eagerly at
                                 leg start, or from the leg's own topology at
                                 its end on a hankin step. */
  PolyMesh dual_bridge_ambo; /**< ambo(P) held across a DUAL bridge: leg 1's
                                 arrival grouping and leg 2's departed mesh. */
  size_t dual_bridge_ambo_faces =
      0; /**< ambo(P) face count, kept for leg 3's handoff length after the mesh
            itself is dropped at the medial leg (persistent-budget relief). */
  /** Device persistent budget of the current shape's split, set by spawn_shape
   * before any read; the host arena is over-provisioned, so gates check the
   * resident persistent high-water against this. */
  size_t device_persistent_budget = 0;
  const Animation::OpLeg::Landing *build_landing =
      nullptr; /**< Latest leg's arrival data (leg-arena backed). */
  const uint8_t *build_from_pal =
      nullptr; /**< Per-face palette the previous leg landed on; survives the
                  leg-boundary compaction that drops its landing. */
  size_t build_from_faces = 0; /**< Length of build_from_pal. */
  int dual_bridges_built = 0;  /**< DUAL bridges scheduled (test coverage). */
  int build_macro_sweep_frames = SWEEP_LEG_FRAMES; /**< Truncate leg of a smooth
                                                       kis/needle macro. */
  int build_reconcile_frames =
      RECONCILE_LEG_FRAMES; /**< Reconcile leg length. */
  enum class BuildContinuation : uint8_t {
    FINISH,
    DUAL_MEDIAL,
    DUAL_UNTRUNCATE,
    DUAL_DONE,
    DT_AFTER_TRUNCATE,
    DT_AFTER_BRIDGE,
    DTD_AFTER_BRIDGE1,
    DTD_AFTER_TRUNCATE,
    DTD_AFTER_BRIDGE2,
  };
  BuildContinuation dual_bridge_done = BuildContinuation::FINISH;

  /**
   * @brief Draw callback for build-leg frames.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, MeshState &, const Animation::OpLeg::Shading &), 8>
      draw_build_fn{
          [this](Canvas &c, MeshState &m, const Animation::OpLeg::Shading &sh) {
            host().draw_build_mesh(c, m, sh);
          }};

  /**
   * @brief Leg frame budget of one lowered primitive step, before the Trans
   *        Speed divisor.
   * @param op Lowered primitive op.
   * @return Frames the leg runs for.
   */
  static int leg_frames(Solids::Op op) {
    switch (op) {
    case Solids::Op::HANKIN:
      return HANKIN_LEG_FRAMES;
    case Solids::Op::RELAX:
      return RELAX_LEG_FRAMES;
    case Solids::Op::DUAL:
      // The smooth dual is a three-leg bridge (truncate to ambo, medial slerp,
      // truncate down to the dual), each a normal single-mesh sweep.
      return 3 * SWEEP_LEG_FRAMES;
    default:
      return SWEEP_LEG_FRAMES;
    }
  }

  /**
   * @brief Whether a lowered DUAL step pairs with a trailing KIS (needle = kd =
   *        dt): the pair builds as the smooth dt macro over both steps.
   * @param k Lowered step index.
   */
  HS_COLD_MEMBER bool dt_pair_at(size_t k) const {
    return k + 1 < build_step_count &&
           build_step_chain[k].op == Solids::Op::DUAL &&
           build_step_chain[k + 1].op == Solids::Op::KIS;
  }

  /**
   * @brief Whether a lowered KIS step stands alone (not the tail of a dt pair);
   *        such a kis builds as the dtd macro (kis = dtd).
   * @param k Lowered step index.
   */
  HS_COLD_MEMBER bool standalone_kis_at(size_t k) const {
    return build_step_chain[k].op == Solids::Op::KIS &&
           !(k > 0 && build_step_chain[k - 1].op == Solids::Op::DUAL);
  }

  /**
   * @brief Whether the lowered chain runs a smooth kis/needle bridge (a dt pair
   *        or a standalone kis), which spawns on the scratch_a-heavy split.
   * @details Scanned at spawn to pick the per-shape arena split before the build
   * grows the arenas. Every such shape fits the bridge split (needle is the
   * largest).
   */
  bool build_uses_smooth_bridge() const {
    for (size_t k = 0; k < build_step_count; ++k)
      if (dt_pair_at(k) || standalone_kis_at(k))
        return true;
    return false;
  }

  /**
   * @brief Log label of a lowered primitive step.
   * @param op Lowered primitive op.
   * @return Static name string.
   */
  static const char *leg_name(Solids::Op op) {
    switch (op) {
    case Solids::Op::HANKIN:
      return "hankin";
    case Solids::Op::AMBO:
      return "ambo";
    case Solids::Op::TRUNCATE:
      return "truncate";
    case Solids::Op::SNUB:
      return "snub";
    case Solids::Op::CHAMFER:
      return "chamfer";
    case Solids::Op::RELAX:
      return "relax";
    case Solids::Op::KIS:
      return "kis";
    case Solids::Op::DUAL:
      return "dual";
    default:
      return "?";
    }
  }

  /**
   * @brief Lays out the build chain's per-leg frame budgets and resets the
   *        build cursor to the chain's first leg.
   * @param sp Trans Speed divisor, >= 1.
   * @return Frames the whole chain occupies.
   * @details A smooth kis/needle macro spans more legs than its lowered step
   * count: a trailing dual,kis (dt) is truncate + dual bridge + reconcile; a
   * standalone kis (dtd) is dual bridge + truncate + dual bridge + reconcile.
   * build_leg_frames[k] carries the dual-bridge budget the bridge splits by
   * three; the truncate and reconcile legs draw fixed member budgets.
   */
  HS_COLD_MEMBER int plan_build_legs(float sp) {
    build_step = 0;
    build_from_pal = nullptr;
    build_from_faces = 0;
    build_total_frames = 0;
    build_macro_sweep_frames =
        std::max(1, static_cast<int>(SWEEP_LEG_FRAMES / sp));
    build_reconcile_frames =
        std::max(1, static_cast<int>(RECONCILE_LEG_FRAMES / sp));
    const int bridge_frames =
        std::max(1, static_cast<int>(leg_frames(Solids::Op::DUAL) / sp));
    // A step this plan skips (a dt pair's trailing kis) must read 0, not the
    // previous build's budget, so a cursor that lands there is caught.
    for (int &frames : build_leg_frames)
      frames = 0;
    for (size_t k = 0; k < build_step_count; ++k) {
      const Solids::Op op = build_step_chain[k].op;
      if (dt_pair_at(k)) {
        build_leg_frames[k] = bridge_frames; // dual bridge (DUAL step)
        build_total_frames +=
            build_macro_sweep_frames + bridge_frames + build_reconcile_frames;
        ++k; // the trailing kis is consumed by the dt macro
        continue;
      }
      if (standalone_kis_at(k)) {
        build_leg_frames[k] = bridge_frames; // both dtd bridges split this
        build_total_frames += bridge_frames + build_macro_sweep_frames +
                              bridge_frames + build_reconcile_frames;
        continue;
      }
      const int frames = std::max(1, static_cast<int>(leg_frames(op) / sp));
      build_leg_frames[k] = frames;
      build_total_frames += frames;
    }
    return build_total_frames;
  }

  /**
   * @brief Constructs and schedules build leg build_step: eagerly builds the
   *        clean endpoint seed and its bookend classification, derives the
   *        palette handoff, and chains completion into finish_build_leg.
   */
  HS_COLD_MEMBER void start_build_leg() {
    const size_t k = build_step;
    const Solids::OpStep &step = build_step_chain[k];

    // Smooth kis/needle macros (docs/specs/opchain_morph_spec.md, smooth
    // kis/needle): a trailing dual,kis is the dt macro (spanning both steps),
    // a standalone kis is the dtd macro; each ends on a reconcile leg onto the
    // exact authored mesh. The recipe/expand_to_primitives still lower needle
    // to {DUAL,KIS}.
    if (dt_pair_at(k)) {
      schedule_dt_macro();
      return;
    }
    if (standalone_kis_at(k)) {
      schedule_dtd_macro();
      return;
    }

    // A DUAL reaching here is a lone one: dt_pair_at() carries no eligibility
    // predicate, so every DUAL,KIS pair took the macro above. The lone DUAL is
    // the smooth three-leg bridge; it builds its own endpoints and chains its
    // legs, then rejoins at finish_build_leg.
    if (step.op == Solids::Op::DUAL) {
      schedule_dual_bridge(BuildContinuation::FINISH);
      return;
    }

    // Eager clean endpoint seed_{k+1}: the mesh the leg lands on and the next
    // leg sweeps from. Runs first — generate() resets the scratch arenas the
    // handoff arrays below live in. A hankin leg builds none; its arrival is
    // the mesh its baked topology already carries. Colours re-key per leg: the
    // arrival's classification maps to a freshly shuffled palette set, and
    // every face crossfades from the previous leg's landing over the leg
    // (Animation::OpLeg palette handoff, core/animation/opleg.h).
    Animation::OpLeg::BookendClasses bookend;
    if (step.op != Solids::Op::HANKIN) {
      hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        build_next_seed =
            Solids::finalize_solid(clean_endpoint(step, a, b), target);
      });
      bookend = next_seed_bookend();
    }

    // Handoff arrays are ctor-scoped: scratch-backed under this scope, alive
    // through the OpLeg constructor's own LIFO-stacked scratch scopes.
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);

    const int frames = build_leg_frames[k];
    HS_CHECK(frames > 0, "RecipeBuild: build leg on a step plan_build_legs "
                         "never budgeted");
    hs::log("Build leg: %s (%d frames)", leg_name(step.op), frames);

    // The swept operator's endpoints: every inflate leg opens at the clamped
    // zero-area birth limit and lands on the step's own parameter; ambo is a
    // truncate swept to the short-circuit point.
    switch (step.op) {
    case Solids::Op::HANKIN:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                            .theta_end = step.param,
                                            .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::RELAX:
      schedule_build_leg(
          Animation::OpLeg(build_seed,
                           Animation::OpLeg::RelaxSpec{
                               .iterations = static_cast<int>(step.param),
                               .bake = step.bake,
                               .sweep_frames = frames},
                           persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::AMBO:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                           .t_start = 0.0f,
                                           .t_end = 0.5f,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::TRUNCATE:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::SNUB:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .twist_end = step.twist,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::CHAMFER:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::CHAMFER,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    default:
      // Neither DUAL nor KIS reaches here: both route through the smooth
      // bridges above.
      HS_CHECK(false, "RecipeBuild: unsweepable primitive op reached a leg");
      break;
    }
  }

  HS_COLD_MEMBER void
  schedule_build_leg(Animation::OpLeg &&leg,
                     BuildContinuation next = BuildContinuation::FINISH) {
    check_build_budget();
    build_landing = &leg.landing();
    Animation::OpLeg::require_event_slot();
    host().timeline.add(
        0, std::move(leg).then([this, next] { continue_build(next); }));
  }

  void check_build_budget() const {
    HS_CHECK(persistent_arena.get_offset() <= device_persistent_budget,
             "RecipeBuild: build leg exceeds the device persistent budget");
  }

  HS_COLD_MEMBER void continue_build(BuildContinuation next) {
    switch (next) {
    case BuildContinuation::FINISH:
      finish_build_leg();
      break;
    case BuildContinuation::DUAL_MEDIAL:
      schedule_dual_medial();
      break;
    case BuildContinuation::DUAL_UNTRUNCATE:
      schedule_dual_untruncate();
      break;
    case BuildContinuation::DUAL_DONE:
      continue_build(dual_bridge_done);
      break;
    case BuildContinuation::DT_AFTER_TRUNCATE:
      dt_after_truncate();
      break;
    case BuildContinuation::DT_AFTER_BRIDGE:
      dt_after_bridge();
      break;
    case BuildContinuation::DTD_AFTER_BRIDGE1:
      dtd_after_bridge1();
      break;
    case BuildContinuation::DTD_AFTER_TRUNCATE:
      dtd_after_truncate();
      break;
    case BuildContinuation::DTD_AFTER_BRIDGE2:
      dtd_after_bridge2();
      break;
    }
  }

  /**
   * @brief Classifies the eagerly built endpoint into topology groups and
   * returns the bookend grouping keyed on it.
   * @return Bookend classes over build_next_seed's fresh classification.
   * @details ScratchScope-guarded so the caller's prior allocations in the
   * shared scratch arenas survive; a bare reset() would drop them.
   */
  HS_COLD_MEMBER Animation::OpLeg::BookendClasses next_seed_bookend() {
    MeshOps::classify_faces_by_topology(build_next_seed, scratch_arena_a,
                                        scratch_arena_b, persistent_arena);
    const size_t faces = build_next_seed.face_counts.size();
    HS_CHECK(faces <= MAX_BUILD_FACES,
             "RecipeBuild: leg endpoint exceeds MAX_BUILD_FACES");
    return {.topology = build_next_seed.topology.data(), .faces = faces};
  }

  /**
   * @brief Captures palette provenance for a leg departing the current seed.
   * @param scratch Arena the centroid and palette arrays live in; must outlive
   * the OpLeg constructor that reads them.
   * @param correspondence Departed-to-swept face-order relationship.
   */
  HS_COLD_MEMBER Animation::OpLeg::PaletteHandoff
  seed_handoff(Arena &scratch,
               Animation::OpLeg::FaceCorrespondence correspondence =
                   Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
    const size_t prev_faces = build_seed.face_counts.size();
    HS_CHECK(prev_faces <= MAX_BUILD_FACES,
             "RecipeBuild: leg seed exceeds MAX_BUILD_FACES");
    math::Vector *prev_centroid = nullptr;
    if (correspondence == Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
      prev_centroid = scratch.allocate_n<math::Vector>(prev_faces);
      Animation::OpLeg::face_centroids_into(build_seed, prev_centroid);
    }
    const uint8_t *prev_pal;
    if (!build_from_pal) {
      // The build's first leg departs the carousel seed slot: its per-face
      // spawn colours are the chain's FROM state. Keyed on build_from_pal
      // (reset to null per build) rather than build_step == 0, which a smooth
      // kis/needle macro's later sub-legs can still sit on.
      HS_CHECK(prev_faces <= host().carousel.current().topology.size(),
               "RecipeBuild: spawn palette does not cover the leg seed");
      prev_pal = host().slot_face_palette[host().carousel.front_index()];
    } else {
      // Depart from the palette the previous leg landed on.
      HS_CHECK(build_from_pal && build_from_faces == prev_faces,
               "RecipeBuild: carried palette does not cover the leg seed");
      prev_pal = build_from_pal;
    }
    return {.bank = &host().palette_bank.bank,
            .prev_face_palette = prev_pal,
            .prev_faces = prev_faces,
            .prev_face_centroid = prev_centroid,
            .correspondence = correspondence};
  }

  /**
   * @brief Captures palette provenance for a leg departing a prior landing.
   * @param departed Mesh the next leg departs from, in the previous leg's
   * landing face order (so its face f carries build_landing face f's palette).
   * @param scratch Arena the arrays live in.
   * @param correspondence Departed-to-swept face-order relationship.
   */
  HS_COLD_MEMBER Animation::OpLeg::PaletteHandoff
  landing_handoff(const PolyMesh &departed, Arena &scratch,
                  Animation::OpLeg::FaceCorrespondence correspondence =
                      Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
    const size_t nf = departed.face_counts.size();
    HS_CHECK(nf <= MAX_BUILD_FACES && build_landing &&
                 build_landing->faces >= nf,
             "RecipeBuild: landing does not cover the departed mesh");
    math::Vector *cen = nullptr;
    uint8_t *pal = scratch.allocate_n<uint8_t>(nf);
    if (correspondence == Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
      cen = scratch.allocate_n<math::Vector>(nf);
      Animation::OpLeg::face_centroids_into(departed, cen);
    }
    for (size_t f = 0; f < nf; ++f)
      pal[f] = build_landing->landed_palette(f);
    return {.bank = &host().palette_bank.bank,
            .prev_face_palette = pal,
            .prev_faces = nf,
            .prev_face_centroid = cen,
            .correspondence = correspondence};
  }

  /**
   * @brief Frame count of DUAL bridge sub-leg `sub` (0/1/2), summing to the
   * step's budget.
   */
  int dual_sub_frames(int sub) const {
    const int total = build_leg_frames[build_step];
    HS_CHECK(total > 0, "RecipeBuild: dual bridge on a step plan_build_legs "
                        "never budgeted");
    const int third = std::max(1, total / 3);
    return sub < 2 ? third : std::max(1, total - 2 * third);
  }

  /**
   * @brief Schedules the smooth dual as three legs: truncate P -> ambo(P), a
   * medial slerp to ambo(dual(P)), and truncate dual(P) back down to dual(P).
   * @details Only ambo(P) (leg 1's arrival) is built up front; dual(P) is
   * deferred to leg 3 so it never co-resides with the medial leg's peak. Each
   * leg's scheduler compacts the arena before it runs, reclaiming the finished
   * legs and the endpoints they no longer need -- the heaviest gyro and
   * ambo_dual seeds run the whole bridge co-resident ~21 KB over budget.
   * @param done Build stage entered after the closing leg.
   */
  HS_COLD_MEMBER void schedule_dual_bridge(BuildContinuation done) {
    HS_CHECK(done == BuildContinuation::FINISH ||
                 done == BuildContinuation::DT_AFTER_BRIDGE ||
                 done == BuildContinuation::DTD_AFTER_BRIDGE1 ||
                 done == BuildContinuation::DTD_AFTER_BRIDGE2,
             "RecipeBuild: invalid dual bridge continuation");
    dual_bridge_done = done;
    ++dual_bridges_built;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      dual_bridge_ambo =
          Solids::finalize_solid(MeshOps::ambo(build_seed, a, b), target);
    });
    HS_CHECK(dual_bridge_ambo.face_counts.size() <= MAX_BUILD_FACES,
             "RecipeBuild: dual bridge ambo exceeds MAX_BUILD_FACES");

    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);
    const int frames = dual_sub_frames(0);
    hs::log("Build leg: dual bridge 1/3 truncate->ambo (%d frames)", frames);
    Animation::OpLeg leg(
        Animation::OpLeg::SweepSeed::borrow(build_seed),
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.0f,
                                         .t_end = 0.5f,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true},
        persistent_arena, draw_build_fn, handoff);
    schedule_build_leg(std::move(leg), BuildContinuation::DUAL_MEDIAL);
  }

  /**
   * @brief Schedules the dual bridge's medial-slerp leg (ambo(P) ->
   * ambo(dual(P))), departing from leg 1's landing. No bookend: the medial
   * arrival's own classification is the target leg 3 departs from.
   */
  HS_COLD_MEMBER void schedule_dual_medial() {
    ScratchScope handoff_guard(scratch_arena_a);
    // Copy leg 1's landed palette into scratch_a before compacting away that
    // finished leg; IDENTITY takes no centroids, since the medial leg sweeps
    // ambo(P)'s faces in place. The persistent reset below leaves scratch_a
    // intact.
    Animation::OpLeg::PaletteHandoff handoff =
        landing_handoff(dual_bridge_ambo, scratch_arena_a,
                        Animation::OpLeg::FaceCorrespondence::IDENTITY);
    const size_t medial_faces = dual_bridge_ambo.face_counts.size();
    HS_CHECK(build_landing && build_landing->faces == medial_faces,
             "RecipeBuild: medial bookend does not match the landing");
    uint16_t *medial_topology =
        scratch_arena_a.allocate_n<uint16_t>(medial_faces);
    std::copy_n(build_landing->topology, medial_faces, medial_topology);
    build_landing = nullptr;

    // Only ambo(P)'s face count survives to leg 3 (its handoff length); the mesh
    // is dead once the landing palette and topology above are copied.
    dual_bridge_ambo_faces = medial_faces;
    dual_bridge_ambo = PolyMesh();

    // Compact: keep only the seed P (leg 2 rebuilds its medial from it, leg 3 its
    // dual), drop leg 1 and the ambo endpoint.
    {
      Persist<PolyMesh> ps(build_seed, scratch_arena_b, persistent_arena);
      host().carousel.compact_drop_all(
          [this](Arena &arena) { host().reclaim_persistent(arena); });
    }

    const int frames = dual_sub_frames(1);
    hs::log("Build leg: dual bridge 2/3 medial (%d frames)", frames);
    Animation::OpLeg leg(
        build_seed, Animation::OpLeg::MedialSpec{.sweep_frames = frames},
        persistent_arena, draw_build_fn, handoff,
        Animation::OpLeg::BookendClasses{.topology = medial_topology,
                                         .faces = medial_faces});
    schedule_build_leg(std::move(leg), BuildContinuation::DUAL_UNTRUNCATE);
  }

  /**
   * @brief Schedules the dual bridge's closing leg: truncate dual(P) from the
   * ambo point down to dual(P), landing on the macro's clean endpoint. The
   * departed palettes come from leg 2's landing in its face order.
   */
  HS_COLD_MEMBER void schedule_dual_untruncate() {
    ScratchScope a_guard(scratch_arena_a);
    const size_t nf = dual_bridge_ambo_faces;
    HS_CHECK(nf <= MAX_BUILD_FACES && build_landing &&
                 build_landing->faces >= nf,
             "RecipeBuild: dual bridge landing does not cover ambo(P)");
    uint8_t *pal = scratch_arena_a.allocate_n<uint8_t>(nf);
    for (size_t f = 0; f < nf; ++f)
      pal[f] = build_landing->landed_palette(f);
    build_landing = nullptr;

    // Compact: keep the seed P (leg 3 builds dual(P) from it), drop leg 2 (the
    // ambo(P) endpoint was already dropped at the medial leg).
    {
      Persist<PolyMesh> ps(build_seed, scratch_arena_b, persistent_arena);
      host().carousel.compact_drop_all(
          [this](Arena &arena) { host().reclaim_persistent(arena); });
    }

    // Build the deferred leg-3 seed dual(P) into the compacted arena. Not
    // generate(): its depth-0 scratch_a reset would drop the palette snapshot
    // above, so scope the pipeline off the live frame instead.
    {
      ScratchScope da(scratch_arena_a);
      ScratchScope db(scratch_arena_b);
      build_next_seed = Solids::finalize_solid(
          MeshOps::dual(build_seed, scratch_arena_a, scratch_arena_b),
          persistent_arena);
    }
    Animation::OpLeg::BookendClasses bookend = next_seed_bookend();

    Animation::OpLeg::PaletteHandoff handoff{
        .bank = &host().palette_bank.bank,
        .prev_face_palette = pal,
        .prev_faces = nf,
        .correspondence = Animation::OpLeg::FaceCorrespondence::DUAL_CLOSING};
    const int frames = dual_sub_frames(2);
    hs::log("Build leg: dual bridge 3/3 truncate->dual (%d frames)", frames);
    Animation::OpLeg leg(
        Animation::OpLeg::SweepSeed::borrow(build_next_seed),
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.5f,
                                         .t_end = 0.0f,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true},
        persistent_arena, draw_build_fn, handoff, bookend);
    schedule_build_leg(std::move(leg), BuildContinuation::DUAL_DONE);
  }

  /**
   * @brief Adopts the eagerly built endpoint (build_next_seed) as the next
   * leg's seed, snapshots the palette the finished leg landed on so its
   * successor departs continuously, and compacts the finished leg's storage.
   * @details The shared middle of finish_build_leg's non-last path, reused by
   * the smooth kis/needle macro stages, which chain their own next leg rather
   * than advancing build_step.
   */
  HS_COLD_MEMBER void carry_landing_to_seed() {
    const size_t landed_faces = build_next_seed.face_counts.size();
    HS_CHECK(build_landing && landed_faces <= build_landing->faces,
             "RecipeBuild: next seed larger than the leg landing");
    HS_CHECK(landed_faces <= MAX_BUILD_FACES,
             "RecipeBuild: next seed exceeds the slot palette capacity");
    // The carry lives in the idle slot's face-palette array: the outgoing
    // shape was dropped at the recipe spawn, and the array's same-address
    // re-claim keeps the bytes across every boundary compaction.
    uint8_t *carry =
        host().slot_face_palette[1 - host().carousel.front_index()];
    for (size_t f = 0; f < landed_faces; ++f)
      carry[f] = build_landing->landed_palette(f);
    build_landing = nullptr;

    {
      Persist<PolyMesh> pn(build_next_seed, scratch_arena_b, persistent_arena);
      build_seed = PolyMesh();
      host().carousel.compact_drop_all(
          [this](Arena &arena) { host().reclaim_persistent(arena); });
    }
    build_seed = std::move(build_next_seed);

    build_from_pal = carry;
    build_from_faces = landed_faces;
  }

  /**
   * @brief Schedules a plain truncate sweep of the current build seed to
   * MACRO_TRUNCATE_T, landing on build_next_seed = truncate(seed, 1/3).
   * @param log Log label ("dt truncate" / "dtd truncate").
   * @param next Build stage entered after the leg.
   */
  HS_COLD_MEMBER void schedule_macro_truncate(const char *log,
                                              BuildContinuation next) {
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      build_next_seed = Solids::finalize_solid(
          MeshOps::truncate(build_seed, a, b, MACRO_TRUNCATE_T), target);
    });
    Animation::OpLeg::BookendClasses bookend = next_seed_bookend();
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);
    const int frames = build_macro_sweep_frames;
    hs::log("Build leg: %s (%d frames)", log, frames);
    Animation::OpLeg leg(
        Animation::OpLeg::SweepSeed::borrow(build_seed),
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.0f,
                                         .t_end = MACRO_TRUNCATE_T,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true},
        persistent_arena, draw_build_fn, handoff, bookend);
    schedule_build_leg(std::move(leg), next);
  }

  /**
   * @brief Schedules the smooth dt macro for a trailing dual,kis (needle = kd =
   * dt): truncate(X, 1/3) sweep, dual bridge on it, then reconcile onto the
   * exact kis(dual(X)) mesh. Covers both the DUAL step and its trailing KIS.
   */
  HS_COLD_MEMBER void schedule_dt_macro() {
    schedule_macro_truncate("dt truncate",
                            BuildContinuation::DT_AFTER_TRUNCATE);
  }
  HS_COLD_MEMBER void dt_after_truncate() {
    carry_landing_to_seed(); // build_seed = truncate(X, 1/3)
    schedule_dual_bridge(BuildContinuation::DT_AFTER_BRIDGE);
  }
  HS_COLD_MEMBER void dt_after_bridge() {
    carry_landing_to_seed(); // build_seed = dual(truncate(X, 1/3)) (identity)
    HS_CHECK(dt_pair_at(build_step),
             "RecipeBuild: dt macro left its dual,kis pair");
    ++build_step; // advance onto the KIS index the reconcile finishes at
    schedule_reconcile(build_step - 1, /*kis_of_dual=*/true);
  }

  /**
   * @brief Schedules the smooth dtd macro for a standalone kis (kis = dtd):
   * dual bridge on X, truncate(dual(X), 1/3) sweep, dual bridge on it, then
   * reconcile onto the exact kis(X) mesh. Runs entirely on the KIS step.
   */
  HS_COLD_MEMBER void schedule_dtd_macro() {
    schedule_dual_bridge(BuildContinuation::DTD_AFTER_BRIDGE1);
  }
  HS_COLD_MEMBER void dtd_after_bridge1() {
    carry_landing_to_seed(); // build_seed = dual(X)
    schedule_macro_truncate("dtd truncate",
                            BuildContinuation::DTD_AFTER_TRUNCATE);
  }
  HS_COLD_MEMBER void dtd_after_truncate() {
    carry_landing_to_seed(); // build_seed = truncate(dual(X), 1/3)
    schedule_dual_bridge(BuildContinuation::DTD_AFTER_BRIDGE2);
  }
  HS_COLD_MEMBER void dtd_after_bridge2() {
    carry_landing_to_seed(); // build_seed = dual(truncate(dual(X), 1/3))
    HS_CHECK(standalone_kis_at(build_step),
             "RecipeBuild: dtd macro left its kis step");
    schedule_reconcile(build_step, /*kis_of_dual=*/false);
  }

  /**
   * @brief Schedules the reconcile leg closing a smooth kis/needle macro: a
   * per-vertex great-circle slerp from the identity mesh (build_seed) onto the
   * exact authored kis/needle mesh, landing on the generator's mesh.
   * @param x_prefix Lowered-step count replayed to rebuild X, the mesh the macro
   * departed from (the generator's exact intermediate).
   * @param kis_of_dual True for a dt macro (authored = kis(dual(X)) = needle(X)),
   * false for a dtd macro (authored = kis(X)).
   */
  HS_COLD_MEMBER void schedule_reconcile(size_t x_prefix, bool kis_of_dual) {
    const uint8_t seed = build_entry->recipe->seed;
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      PolyMesh X = Solids::build_steps(seed, build_step_chain, x_prefix,
                                       scratch_arena_a, scratch_arena_b);
      PolyMesh Xc;
      MeshOps::clone(X, Xc, scratch_arena_a);
      PolyMesh authored =
          kis_of_dual ? MeshOps::needle(Xc, scratch_arena_a, scratch_arena_b)
                      : MeshOps::kis(Xc, scratch_arena_a, scratch_arena_b);
      MeshOps::reconcile_vertices(build_seed, authored, build_next_seed,
                                  persistent_arena, scratch_arena_a);
    }
    Animation::OpLeg::BookendClasses bookend = next_seed_bookend();
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(
        scratch_arena_a, Animation::OpLeg::FaceCorrespondence::IDENTITY);
    const int frames = build_reconcile_frames;
    hs::log("Build leg: reconcile (%d frames)", frames);
    Animation::OpLeg leg(build_seed,
                         Animation::OpLeg::ReconcileSpec{
                             .to_positions = build_next_seed.vertices.data(),
                             .to_count = build_next_seed.vertices.size(),
                             .sweep_frames = frames},
                         persistent_arena, draw_build_fn, handoff, bookend);
    schedule_build_leg(std::move(leg));
  }

  /**
   * @brief Builds the clean endpoint mesh a leg lands on.
   * @param step Lowered primitive step; a hankin step has no entry, since its
   * leg rebuilds its own arrival from the topology it swept.
   * @param a Output arena for even pipeline stages.
   * @param b Scratch arena for odd pipeline stages.
   * @return The op applied to the current build seed at its exact parameter.
   */
  HS_COLD_MEMBER PolyMesh clean_endpoint(const Solids::OpStep &step, Arena &a,
                                         Arena &b) {
    switch (step.op) {
    case Solids::Op::AMBO:
      return MeshOps::ambo(build_seed, a, b);
    case Solids::Op::TRUNCATE:
      return MeshOps::truncate(build_seed, a, b, step.param);
    case Solids::Op::SNUB:
      return MeshOps::snub(build_seed, a, b, step.param, step.twist);
    case Solids::Op::CHAMFER:
      return MeshOps::chamfer(build_seed, a, b, step.param);
    case Solids::Op::RELAX:
      if (step.bake)
        return MeshOps::relax_baked(build_seed, a, *step.bake);
      return MeshOps::relax(build_seed, a, b, static_cast<int>(step.param));
    default:
      HS_CHECK(false, "RecipeBuild: step builds no eager endpoint");
      return PolyMesh{};
    }
  }

  /**
   * @brief Leg completion: adopts the eagerly built clean endpoint as the next
   *        leg's seed, then starts the next leg or finishes the build.
   */
  HS_COLD_MEMBER void finish_build_leg() {
    // Reclaim the finished leg. Only the endpoint the next leg sweeps from
    // crosses the reset. A hankin leg's endpoint is rebuilt here from the
    // topology it swept, into the scratch the evacuation below reads; the other
    // kinds carry the endpoint start_build_leg built eagerly. The palette the
    // leg landed on is snapshotted too, since the landing does not survive.
    {
      ScratchScope a_guard(scratch_arena_a);
      HS_CHECK(build_step < build_step_count,
               "RecipeBuild: build cursor ran past the lowered chain");
      if (build_step_chain[build_step].op == Solids::Op::HANKIN) {
        HS_CHECK(build_landing, "RecipeBuild: finished leg has no landing");
        Animation::OpLeg::arrival_mesh(*build_landing, build_next_seed,
                                       scratch_arena_a);
      }

      if (build_step + 1 >= build_step_count) {
        // finish_build consumes the last leg's landing, so that one is reclaimed
        // by the closing compaction instead.
        build_seed = std::move(build_next_seed);
        ++build_step;
        finish_build();
        return;
      }

      // Carry the emission-order prefix the next leg departs from (the whole
      // landing for a face-count-preserving leg, the survivor prefix where a leg's
      // arrival has more faces than its clean endpoint -- the DUAL bridge's
      // closing truncate lands V+F faces but hands off the V dual faces), then
      // advance to the next lowered step.
      carry_landing_to_seed();
      ++build_step;
    }
    start_build_leg();
  }

  /**
   * @brief Build completion: recompiles the finished solid into the front slot
   *        and hands the last leg's per-face colours to the sprite, so its
   *        next frame is pixel-equal to the leg's last one.
   */
  HS_COLD_MEMBER void finish_build() {
    // The finished solid is build_seed (the last leg's clean endpoint): its
    // face count is the emission-order prefix of the landing (the whole landing
    // for a normal leg, the surviving dual faces for the DUAL bridge's closing
    // truncate, whose zero-area corner births drop at the compile below).
    const size_t landed_faces = build_seed.face_counts.size();
    HS_CHECK(build_landing && landed_faces <= build_landing->faces,
             "RecipeBuild: finished solid larger than the leg landing");
    HS_CHECK(landed_faces <= MAX_BUILD_FACES,
             "RecipeBuild: finished solid exceeds the slot palette capacity");
    HS_CHECK(build_landing->topology,
             "RecipeBuild: finished leg has no topology");
    ScratchScope topology_guard(scratch_arena_b);
    uint16_t *landed_topology =
        scratch_arena_b.allocate_n<uint16_t>(landed_faces);
    std::copy_n(build_landing->topology, landed_faces, landed_topology);
    const int front = host().carousel.front_index();
    // Per-face sprite handoff: copied before the compaction below, whose
    // same-address re-claim keeps the array's bytes.
    for (size_t f = 0; f < landed_faces; ++f)
      host().slot_face_palette[front][f] = build_landing->landed_palette(f);
    build_landing = nullptr;
    build_from_pal = nullptr;
    build_from_faces = 0;

    MeshState &slot = host().carousel.slot(front);
    {
      ScratchScope a_guard(scratch_arena_a);
      slot.clear();
      // The seed the closing compile reads is evacuated to scratch and never
      // restored; the compiled slot is the only form the effect renders. It is
      // dropped again before the classification, which needs scratch_arena_b in
      // full.
      {
        ScratchScope seed_guard(scratch_arena_b);
        PolyMesh built;
        MeshOps::clone(build_seed, built, scratch_arena_b);
        build_next_seed = PolyMesh();
        build_seed = PolyMesh();
        host().carousel.compact_drop_all(
            [this](Arena &arena) { host().reclaim_persistent(arena); });
        MeshOps::compile(built, slot, persistent_arena, scratch_arena_a);
      }
      slot.topology.bind(persistent_arena, landed_faces);
      for (size_t f = 0; f < landed_faces; ++f)
        slot.topology.push_back(landed_topology[f]);
    }

    // The per-face colours index the compiled slot by emission order, and
    // compile() strips faces with fewer than 3 sides, so the counts must agree
    // exactly.
    HS_CHECK(landed_faces == slot.num_faces(),
             "RecipeBuild: compiled face count differs from the leg landing");
    build_active = false;

    hs::log("Built Shape: %s (V=%d, E=%d, F=%d, I=%d)", build_entry->name,
            (int)slot.vertices.size(), (int)(slot.faces.size() / 2),
            (int)slot.face_counts.size(), (int)slot.faces.size());
  }
};

} // namespace Animation
