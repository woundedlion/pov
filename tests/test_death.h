/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Death tests for the fail-fast (HS_CHECK / __builtin_trap) seams, plus
 * cross-process cold/warm effect determinism.
 *
 * Each case runs in a re-exec'd child (HS_DEATH_CHILD=harness,
 * HS_DEATH_CASE=<name>). The parent requires the illegal-instruction trap
 * status (SIGILL / STATUS_ILLEGAL_INSTRUCTION) and the HS_CHECK breadcrumb of
 * the case's exact guard: under -fsanitize-trap=undefined, UB lowers to the
 * same instruction.
 */
#pragma once

#include "core/animation/orientation.h"
#include <array>
#include <cerrno>
#include <chrono>
#include <thread>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <string>
#include <utility>

#include "death_guard_sites.h" // generated HS_CHECK census
#include "tests/test_fixture.h"
#include "tests/test_pullback.h"
#include "tests/test_effects.h"
#include "tests/test_shapeshifter_oracle.h"
#include "tests/test_harness.h"

#include "core/math/3dmath.h"
#include "core/math/4dmath.h"
#include "core/animation/animation.h"
#include "core/animation/carousel.h"
#include "core/render/canvas.h"
#include "core/color/color.h"
#include "core/color/noise_hue_palette.h"
#include "core/color/noise_shimmer_palette.h"
#include "core/math/geometry.h"
#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/control/choreography.h"
#include "core/control/presets.h"
#include "core/control/registry.h"
#include "core/memory.h"
#include "core/platform/led.h"
#include "core/math/lenses.h"
#include "core/mesh/hankin.h"
#include "core/mesh/mesh.h"
#include "core/mesh/recipe.h"
#include "core/render/plot.h"
#include "core/render/pullback/interpreter.h"
#include "core/render/scan.h"
#include "core/render/sdf.h"
#include "core/render/sdf/volume.h"
#include "core/mesh/solids.h"
#include "core/math/spherical_field.h"
#include "core/math/spherical_harmonics.h"
#include "core/spatial/kd_tree.h"
#include "core/spatial/reaction_graph.h"
#include "core/containers/triangular_bitset.h"
#include "core/containers/static_circular_buffer.h"
#include "core/animation/transformer.h"
#include "hardware/pov_sync.h"
#include "targets/wasm/param_marshal.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <fcntl.h>    // open / O_WRONLY for the capture file
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#include <unistd.h>   // fork / execv / dup2 / close / _exit — shell-free spawn
#else
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <process.h> // _getpid
#endif

namespace hs_test {
namespace death_tests {

/**
 * @brief Accessor for argv[0] of the running test binary, used to re-exec self.
 * @return Reference to the static char-pointer slot; captured in main().
 */
inline const char *&self_exe() {
  static const char *s = nullptr;
  return s;
}

/**
 * @brief Launders a value through a volatile to defeat constant-folding.
 * @tparam T Value type to pass through opaquely.
 * @param v Value to make opaque to the optimizer.
 * @return A copy of @p v the optimizer cannot prove constant.
 * @details Keeps the compiler from proving a trap is taken at compile time and
 *          reshaping the case; each case must trap at run time.
 */
template <typename T> inline T opaque(T v) {
  volatile T x = v;
  return x;
}

#include "tests/death/fixtures.h"
#include "tests/death/animation_cases.h"
#include "tests/death/callables_cases.h"
#include "tests/death/color_cases.h"
#include "tests/death/effects_cases.h"
#include "tests/death/math_spatial_cases.h"
#include "tests/death/memory_cases.h"
#include "tests/death/mesh_cases.h"
#include "tests/death/params_canvas_cases.h"
#include "tests/death/particles_cases.h"
#include "tests/death/plot_filter_cases.h"
#include "tests/death/pullback_cases.h"
#include "tests/death/recipes_cases.h"
#include "tests/death/registry_cases.h"
#include "tests/death/sdf_cases.h"
#include "tests/death/transformers_cases.h"
/**
 * @brief A named death case selected by HS_DEATH_CASE in the child process.
 */
struct Case {
  const char *name;        /**< Case selector matched against HS_DEATH_CASE. */
  void (*fn)();            /**< The trap-triggering case body to run. */
  const char *guard_file;  /**< Repo-relative path of the guard source. */
  const char *guard_text;  /**< Expected "(condition) message" tail of that
                               guard's breadcrumb line. */
  bool debug_only = false; /**< Guard absent in NDEBUG builds. */
};

inline bool case_enabled(const Case &entry) {
#ifdef NDEBUG
  return !entry.debug_only;
#else
  (void)entry;
  return true;
#endif
}

/**
 * @brief Returns the full death-case table.
 * @param n Out-param set to the number of cases in the table.
 * @return Pointer to the static case array.
 */
inline const Case *all_cases(int &n) {
  static const Case cases[] = {
      {"sdf_smooth_union_zero_radius", case_sdf_smooth_union_zero_radius,
       "core/render/sdf/csg.h",
       "(k > 0.0f) SDF CSG: smoothness must be positive"},
      {"sdf_angular_repeat_zero_copies", case_sdf_angular_repeat_zero_copies,
       "core/render/sdf/csg.h",
       "(reps > 0) SDF CSG: repetition count must be positive"},
      {"square_wave_invalid_duty", case_square_wave_invalid_duty,
       "core/math/waves.h",
       "(duty_cycle >= 0.0f && duty_cycle <= 1.0f) square_wave: duty_cycle must be in [0,1]"},
      {"random_timer_max_int", case_random_timer_max_int,
       "core/animation/timers.h",
       "(max < std::numeric_limits<int>::max()) RandomTimer max must be < INT_MAX (reset adds 1)"},
      {"effect_margin_equal_width", case_effect_margin_equal_width,
       "core/render/canvas.h",
       "(m >= 0 && m < clip_region.w) render margin must be in [0, canvas width)"},
      {"fib_spiral_zero_points", case_fib_spiral_zero_points,
       "core/math/spherical.h", "(n > 0) fib_spiral: n must be positive"},
      {"noise_hue_bake_invalid_scale", case_noise_hue_bake_invalid_scale,
       "core/color/noise_hue_palette.h",
       "(std::isfinite(bake_scale) && bake_scale > 0.0f) HueNoiseBakeCache: scale must be finite and positive"},
      {"direct_sink_wrong_dimensions", case_direct_sink_wrong_dimensions,
       "core/render/filter/screen_direct_aa_sink.h",
       "(cv.width() == W && cv.height() == H) DirectAntiAliasSink: framebuffer dimensions do not match"},
      {"recipe_bake_wrong_op", case_recipe_bake_wrong_op, "core/mesh/recipe.h",
       "(!step.bake || step.op == Op::RELAX) apply_step: only RELAX accepts a bake"},
      {"vertex_replicate_short_input", case_vertex_replicate_short_input,
       "core/render/filter/world_vertex_replicate.h",
       "(std::size(vertices) >= static_cast<size_t>(N)) VertexReplicate: "
       "vertex array is smaller than replica count"},
      {"recipe_twist_wrong_op", case_recipe_twist_wrong_op,
       "core/mesh/recipe.h",
       "(step.twist == 0.0f || step.op == Op::SNUB) apply_step: only SNUB accepts a twist"},
      {"star_mismatched_radius_cache", case_star_mismatched_radius_cache,
       "core/render/plot/shapes.h",
       "(cache_matches) Star: cached trigonometry does not match inputs"},
      {"star_mismatched_step_cache", case_star_mismatched_step_cache,
       "core/render/plot/shapes.h",
       "(cache_matches) Star: cached trigonometry does not match inputs"},
      {"reaction_lattice_uninitialized", case_reaction_lattice_uninitialized,
       "effects/ReactionDiffusionBase.h",
       "(nodes != nullptr) ReactionDiffusion: lattice is not initialized"},
      {"recipe_bake_live_iterations", case_recipe_bake_live_iterations,
       "core/mesh/recipe.h",
       "(!step.bake || step.param == 0.0f) apply_step: a baked RELAX step must not specify live iterations"},
      {"pullback_mobius_degenerate", case_pullback_mobius_degenerate,
       "core/render/pullback/operators/sphere.h",
       "(Lens::MobiusLensParams::nondegenerate(mobius)) sphere.lens.mobius: degenerate coefficients"},
      {"shapeshifter_count_over_capacity",
       case_shapeshifter_count_over_capacity, "effects/ShapeShifter.h",
       "(count >= 1 && count <= MAX_SHAPES) ShapeShifter: contour count"},
      {"planar_chords_empty_storage", case_planar_chords_empty_storage,
       "core/render/plot/chords.h",
       "(max_vertices >= 1) PlanarChords: max_vertices"},
      {"planar_chords_over_capacity", case_planar_chords_over_capacity,
       "core/render/plot/chords.h",
       "(vertices >= 1 && vertices <= capacity) PlanarChords:"},
      {"planar_chords_unprepared", case_planar_chords_unprepared,
       "core/render/plot/chords.h",
       "(prepared && clip_stamp == canvas.clip()) PlanarChords:"},
      {"planar_chords_stale_clip", case_planar_chords_stale_clip,
       "core/render/plot/chords.h",
       "(prepared && clip_stamp == canvas.clip()) PlanarChords:"},
      {"planar_band_split_empty_storage", case_planar_band_split_empty_storage,
       "core/render/plot/chords.h",
       "(max_points >= 2) PlanarBandSplit: max_points"},
      {"planar_band_split_nonempty_output",
       case_planar_band_split_nonempty_output, "core/render/plot/chords.h",
       "(out.empty()) PlanarBandSplit::split: out must be empty"},
      {"planar_band_split_over_capacity", case_planar_band_split_over_capacity,
       "core/render/plot/chords.h",
       "(max_points(edges, pieces) <= capacity) PlanarBandSplit:"},
      {"lattice_shells_oob", case_lattice_shells_oob,
       "core/render/sdf/lattice.h",
       "(static_cast<uint8_t>(settings.shells) < MAX_SHELLS) lattice shell count exceeds crossing capacity"},
      {"lattice_zero_softness", case_lattice_zero_softness,
       "core/render/sdf/lattice.h",
       "(settings.softness > 0 && settings.cell_size > 0 && "
       "settings.aa_strength >= 0) lattice requires positive softness"},
      {"lattice_zero_cell_size", case_lattice_zero_cell_size,
       "core/render/sdf/lattice.h",
       "(settings.softness > 0 && settings.cell_size > 0 && "
       "settings.aa_strength >= 0) lattice requires positive softness"},
      {"lattice_negative_aa", case_lattice_negative_aa,
       "core/render/sdf/lattice.h",
       "(settings.softness > 0 && settings.cell_size > 0 && "
       "settings.aa_strength >= 0) lattice requires positive softness"},
      {"hyperlattice_pattern_defaults_invalid",
       case_hyperlattice_pattern_defaults_invalid, "effects/HyperLattice.h",
       "(false) HyperLattice: unsupported pattern"},
      {"hyperlattice_frame_without_crossings",
       case_hyperlattice_frame_without_crossings, "effects/HyperLattice.h",
       "(frame.crossings) HyperLattice: frame has no crossing list"},
      {"hyperlattice_cubic_trace_geometry",
       case_hyperlattice_cubic_trace_geometry, "effects/HyperLattice.h",
       "(false) HyperLattice: pattern has no trace geometry"},
      {"mindsplatter_profile_preset_oob", case_mindsplatter_profile_preset_oob,
       "core/control/choreography.h",
       "(index < authored_preset_count()) profile preset index out of range"},

      {"direct_sink_unprepared_plot", case_direct_sink_unprepared_plot,
       "core/render/filter/screen_direct_aa_sink.h",
       "(prepared_for(cv)) DirectAntiAliasSink: prepare current canvas before plotting"},
      {"timeline_add_into_live_slot", case_timeline_add_into_live_slot,
       "core/animation/timeline.h",
       "(!e.manager) add_get would overwrite a live animation"},
      {"feedback_storage_twice", case_feedback_storage_twice,
       "core/render/filter/pixel_feedback.h",
       "(!cached_warp_x || !stamp.block_alive(cached_warp_x, CACHE_CELLS * sizeof(int16_t))) feedback filter: storage already initialized",
       true},
      {"world_storage_twice", case_world_storage_twice,
       "core/render/filter/world_trails.h",
       "(!items || !stamp.block_alive(items, STORAGE_BYTES)) world filter: storage already initialized",
       true},
      {"screen_storage_twice", case_screen_storage_twice,
       "core/render/filter/screen_trails.h",
       "(!points || !stamp.block_alive(points, STORAGE_BYTES)) screen filter: storage already initialized",
       true},
      {"sample_sphere_nan", case_sample_sphere_nan,
       "core/render/pullback/contract.h",
       "(value == value) unit clamp: NaN input"},
      {"sample_plane_nan", case_sample_plane_nan,
       "core/render/pullback/contract.h",
       "(value == value) unit clamp: NaN input"},
      {"projection_invalid_frame_prepare",
       case_projection_invalid_frame_prepare,
       "core/render/pullback/operators/project.h",
       "(params.frame == static_cast<uint8_t>(ProjectionFrame::IDENTITY) || params.frame == static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER)) projection operator: invalid frame policy"},
      {"projection_invalid_frame_advance",
       case_projection_invalid_frame_advance,
       "core/render/pullback/operators/project.h",
       "(params.frame == static_cast<uint8_t>(ProjectionFrame::IDENTITY) || params.frame == static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER)) projection operator: invalid frame policy"},
      {"peirce_invalid_layout", case_peirce_invalid_layout,
       "core/math/projections.h",
       "(folded_layout || strip_layout) Peirce projection: invalid layout"},
      {"ball_drop_nonfinite_azimuth", case_ball_drop_nonfinite_azimuth,
       "core/animation/params.h",
       "(std::isfinite(azimuth)) BallDrop azimuth must be finite"},
      {"bump_offset_outside_cap_distance",
       case_bump_offset_outside_cap_distance, "core/animation/transformer.h",
       "(std::abs(y) <= d + 1e-3f) bump offset exceeds the angular distance to its center"},
      {"field_transfer_outside_range", case_field_transfer_outside_range,
       "core/render/pullback/stage.h",
       "(value >= 0.0f && value <= 1.0f) field transfer must remain in [0, 1]"},
      {"field_coverage_increases", case_field_coverage_increases,
       "core/render/pullback/stage.h",
       "(factor >= 0.0f && factor <= 1.0f) field coverage factor must remain in [0, 1]"},
      {"field_sample_coverage_outside_range",
       case_field_sample_coverage_outside_range, "core/render/pullback/stage.h",
       "(coverage >= 0.0f && coverage <= 1.0f && input.provenance.domain_coverage >= 0.0f && input.provenance.domain_coverage <= 1.0f) field coverage factors must remain in [0, 1]"},
      {"reconcile_aliased_arenas", case_reconcile_aliased_arenas,
       "core/mesh/conway.h",
       "(&target != &scratch) reconcile_vertices: target and scratch must differ"},
      {"opleg_no_event_slot", case_opleg_no_event_slot,
       "core/animation/opleg.h",
       "(Timeline::remaining() > 0) OpLeg requires a free timeline event"},
      {"arena_oom", case_arena_oom, "core/memory.cpp",
       "(false) Arena::allocate: out of memory"},
      {"arena_make_oom", case_arena_make_oom, "core/memory.cpp",
       "(false) Arena::allocate: out of memory"},
      {"arena_zero_size_alloc", case_arena_zero_size_alloc,
       "core/memory/arena.h", "(size > 0) Arena::allocate: zero-size request"},
      {"arena_allocate_n_overflow", case_arena_allocate_n_overflow,
       "core/memory/arena.h",
       "(n <= SIZE_MAX / sizeof(T)) Arena::allocate_n element count overflows "
       "size_t"},
      {"arena_bad_alignment", case_arena_bad_alignment, "core/memory/arena.h",
       "(align != 0 && (align & (align - 1)) == 0) Arena::allocate: "
       "alignment "},
      {"resplit_scratch_not_empty", case_resplit_scratch_not_empty,
       "core/memory.cpp",
       "(scratch_arena_a.get_offset() == 0 && scratch_arena_b.get_offset() "
       "== 0) resplit_arenas: both scratch arenas must be empty"},
      {"resplit_persistent_strands", case_resplit_persistent_strands,
       "core/memory/arena.h",
       "(offset <= new_capacity) Arena::rebind_capacity below the live offset "
       "would strand content"},
      {"arena_set_offset_forward", case_arena_set_offset_forward,
       "core/memory/arena.h", "(new_offset <= offset) Arena::set_offset: "},
      {"scratch_scope_non_lifo", case_scratch_scope_non_lifo,
       "core/memory/scratch.h",
       "(arena.get_offset() >= saved_offset) ScratchScope: non-LIFO teardown"},
      {"arena_rewind_history_overflow", case_arena_rewind_history_overflow,
       "core/memory/arena.h",
       "(rewind_history_size < REWIND_HISTORY_CAPACITY) Arena: debug rewind history capacity exceeded",
       true},
      {"scratch_scope_reset", case_scratch_scope_reset, "core/memory/scratch.h",
       "(arena.get_generation() == saved_generation) ScratchScope: arena reset during scope lifetime",
       true},
      {"arena_vector_overflow", case_arena_vector_overflow,
       "core/memory/vector.h",
       "(element_count < element_capacity) ArenaVector push_back exact "
       "capacity exceeded!"},
      {"arena_vector_emplace_overflow", case_arena_vector_emplace_overflow,
       "core/memory/vector.h",
       "(element_count < element_capacity) ArenaVector emplace_back exact "
       "capacity exceeded!"},
      {"generate_target_is_scratch", case_generate_target_is_scratch,
       "core/memory/generate.h",
       "(&target != &scratch_arena_a && &target != &scratch_arena_b) "
       "generate: target must not alias an engine scratch arena"},
      {"generate_recursion_too_deep", case_generate_recursion_too_deep,
       "core/memory/generate.h",
       "(depth <= MAX_GENERATE_DEPTH) generate: recursion too deep"},
      {"normalize_zero", case_normalize_zero, "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Vector: zero length"},
      {"vector_normalize_in_place_zero", case_vector_normalize_in_place_zero,
       "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Vector: zero length"},
      {"rotate_plane_degenerate", case_rotate_plane_degenerate,
       "core/math/4dmath.h",
       "(a >= 0 && a < VEC4_DIMENSIONS && b >= 0 && b < VEC4_DIMENSIONS && a "
       "!= b) rotate_plane: "},
      {"angle_between_zero", case_angle_between_zero, "core/math/3dmath.h",
       "(m1 >= math::EPS_LEN_SQ && m2 >= math::EPS_LEN_SQ) "
       "angle_between: degenerate vector"},
      {"normalize_nan", case_normalize_nan, "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Vector: zero length"},
      {"solids_index_oob", case_solids_index_oob, "core/mesh/solids.h",
       "(index < static_cast<size_t>(NUM_ENTRIES)) Solids::get_entry: index "
       "out of range"},
      {"solids_unknown_name", case_solids_unknown_name, "core/mesh/solids.h",
       "(entry) Solids::get_by_name: unknown solid name"},
      {"circular_buffer_oob", case_circular_buffer_oob,
       "core/containers/static_circular_buffer.h",
       "(index < count) StaticCircularBuffer::operator[]: index 5 is "
       "outside [0, 2)"},
      {"circular_buffer_front_empty", case_circular_buffer_front_empty,
       "core/containers/static_circular_buffer.h",
       "(!is_empty()) front() on empty StaticCircularBuffer"},
      {"circular_buffer_back_empty", case_circular_buffer_back_empty,
       "core/containers/static_circular_buffer.h",
       "(!is_empty()) back() on empty StaticCircularBuffer"},
      {"circular_buffer_const_back_empty",
       case_circular_buffer_const_back_empty,
       "core/containers/static_circular_buffer.h",
       "(!is_empty()) back() on empty StaticCircularBuffer"},
      {"arena_vector_append_bulk_overflow",
       case_arena_vector_append_bulk_overflow, "core/memory/vector.h",
       "(count <= element_capacity - element_count) ArenaVector bulk append "
       "exceeds capacity!"},
      {"spatial_knn_over_max", case_spatial_knn_over_max,
       "core/spatial/kd_tree.h",
       "(k <= static_cast<size_t>(MAX_K)) KDTree::nearest k exceeds MAX_K"},
      {"reaction_graph_node_index_out_of_range",
       case_reaction_graph_node_index_out_of_range,
       "core/spatial/reaction_graph.h",
       "(i >= 0 && i < RD_N) node() index outside the lattice"},
      {"reaction_graph_slot_out_of_range",
       case_reaction_graph_slot_out_of_range, "core/spatial/reaction_graph.h",
       "(table[i][k] >= 0 && table[i][k] < RD_N) neighbors[] slot is not a "
       "lattice node index"},
      {"gs_neighbor_exceeds_history", case_gs_neighbor_exceeds_history,
       "effects/GSReactionDiffusion.h",
       "(delta >= -PHYSICS_NEIGHBOR_REACH) GS neighbor exceeds delayed-write history"},
      {"gs_color_noise_zero_scale", case_gs_color_noise_zero_scale,
       "core/color/noise_hue_palette.h",
       "(std::isfinite(bake_scale) && bake_scale > 0.0f) HueNoiseBakeCache: scale must be finite and positive"},
      {"gs_color_noise_nan_scale", case_gs_color_noise_nan_scale,
       "core/color/noise_hue_palette.h",
       "(std::isfinite(bake_scale) && bake_scale > 0.0f) HueNoiseBakeCache: scale must be finite and positive"},
      {"arena_oversubscribed", case_arena_oversubscribed, "core/memory.cpp",
       "(total <= GLOBAL_ARENA_SIZE) split_bases: "},
      {"arena_split_scratch_too_large", case_arena_split_scratch_too_large,
       "core/memory/arena.h",
       "(scratch_a <= total && scratch_b <= total - scratch_a) ArenaSplit: "},
      {"arena_partition_too_large", case_arena_partition_too_large,
       "core/memory.cpp",
       "(persistent <= GLOBAL_ARENA_SIZE && scratch_a <= GLOBAL_ARENA_SIZE && "
       "scratch_b <= GLOBAL_ARENA_SIZE) split_bases: "},
      {"persist_forgot_reset", case_persist_forgot_reset,
       "core/memory/persist.h",
       "(persistent.get_offset() <= persistent_offset_at_ctor) Persist: "
       "restore grew the persistent arena past its construction watermark — "
       "the caller did not rewind/reset it during the scope, so the restore "
       "appended a duplicate instead of reconstructing"},
      {"persist_same_arena", case_persist_same_arena, "core/memory/persist.h",
       "(&scratch_arena != &restore_arena) Persist: scratch and persistent "
       "must be distinct arenas — the dtor's watermark restore assumes the "
       "backup lives in a different arena than the one it restores into"},
      {"triangular_bitset_unordered_pair",
       case_triangular_bitset_unordered_pair,
       "core/containers/triangular_bitset.h",
       "(small >= 0 && small < large && large < MAX_V) "
       "TriangularBitset::index: pair "},
      {"timeline_pinned_relocation", case_timeline_pinned_relocation,
       "core/animation/timeline.h",
       "(!pinned) move_into would dangle a pinned animation's retained "
       "pointer"},
      {"timeline_move_into_live_destination",
       case_timeline_move_into_live_destination, "core/animation/timeline.h",
       "(!dst.manager) move_into would leak the destination's live animation"},
      {"timeline_negative_delay", case_timeline_negative_delay,
       "core/animation/timeline.h",
       "(in_frames >= 0) Timeline delay must be non-negative"},
      {"finite_param_perpetual_duration", case_finite_param_perpetual_duration,
       "core/animation/params.h",
       "(duration >= 0) finite parameter animation duration must be >= 0"},
      {"transition_nonfinite_target", case_transition_nonfinite_target,
       "core/animation/params.h",
       "(std::isfinite(to)) Transition target must be finite"},
      {"timeline_pinned_completion", case_timeline_pinned_completion,
       "core/animation/timeline.h",
       "(!e.pinned || anim->is_canceled()) pinned animation completed; only "
       "cancel() may destroy a pinned event"},
      {"timeline_pinned_finite_animation",
       case_timeline_pinned_finite_animation, "core/animation/timeline.h",
       "(!animation.is_finite() || animation.repeats()) pinned animation "
       "must be infinite or repeating"},
      {"timeline_pinned_add_on_full_timeline",
       case_timeline_pinned_add_on_full_timeline, "core/animation/timeline.h",
       "(pin == Pin::UNPINNED) Timeline full, dropped a pinned animation"},
      {"timeline_pinned_after_cancelled_perpetual",
       case_timeline_pinned_after_cancelled_perpetual,
       "core/animation/timeline.h",
       "(!prev || (!prev->is_canceled() && (!prev->is_finite() || prev->repeats()))) "
       "pinned animation added after a retiring predecessor"},
      {"timeline_pinned_one_shot_timer", case_timeline_pinned_one_shot_timer,
       "core/animation/timeline.h",
       "(!e.pinned || anim->is_canceled()) pinned animation completed; only "
       "cancel() may destroy a pinned event"},
      {"timeline_clear_pinned", case_timeline_clear_pinned,
       "core/animation/timeline.h",
       "(!global_timeline_events[i].pinned) clear() would destroy a pinned "
       "animation"},
      {"timeline_clear_during_step", case_timeline_clear_during_step,
       "core/animation/timeline.h",
       "(!stepping) clear() from inside step() would destroy the animation "
       "whose callback is running"},
      {"timeline_clear_hook_adds_event", case_timeline_clear_hook_adds_event,
       "core/animation/timeline.h",
       "(global_timeline_num_events == event_count) clear hook added or "
       "removed timeline events"},
      {"segue_sprite_no_slot", case_segue_sprite_no_slot,
       "core/animation/segue.h",
       "(Timeline::remaining() >= 1) segue: the transition sprite needs a free "
       "timeline slot"},
      {"mesh_carousel_unflipped_slot", case_mesh_carousel_unflipped_slot,
       "core/animation/carousel.h",
       "(slot == front) MeshCarousel segue scheduled before incoming slot "
       "flip"},
      {"timeline_double_construct", case_timeline_double_construct,
       "core/animation/timeline.h",
       "(!global_timeline_live) a second live Timeline would stomp the shared "
       "global events"},
      {"transformer_pool_init_storage_twice",
       case_transformer_pool_init_storage_twice, "core/animation/transformer.h",
       "(!entities) TransformerPool: init_storage() called twice"},
      {"transformer_unpinned_full", case_transformer_unpinned_full,
       "core/animation/transformer.h",
       "(pin == Timeline::Pin::PINNED || (anim.is_finite() && !anim.repeats())) "
       "Transformer::spawn needs a finite, non-repeating animation"},
      {"transformer_pool_spawn_before_init",
       case_transformer_pool_spawn_before_init, "core/animation/transformer.h",
       "(entities) TransformerPool: call init_storage() before spawn"},
      {"transformer_pool_prepare_frame_before_init",
       case_transformer_pool_prepare_frame_before_init,
       "core/animation/transformer.h",
       "(entities) TransformerPool: call init_storage() before prepare_frame"},
      {"transformer_pool_spawn_pausable_null_flag",
       case_transformer_pool_spawn_pausable_null_flag,
       "core/animation/transformer.h",
       "(paused != nullptr) pausable spawn needs a pause flag"},
      {"transformer_pool_active_index_oob",
       case_transformer_pool_active_index_oob, "core/animation/transformer.h",
       "(k >= 0 && k < active_slot_count) TransformerPool: active index out of "
       "range"},
      {"transformer_pool_reclaim_storage_moved",
       case_transformer_pool_reclaim_storage_moved,
       "core/animation/transformer.h",
       "(e == entities && s == active_slots) TransformerPool: reclaimed "
       "storage moved"},
      {"transformer_pool_arena_reclaimed",
       case_transformer_pool_arena_reclaimed, "core/animation/transformer.h",
       "(storage_arena->get_offset() >= storage_end) TransformerPool: arena "
       "reclaimed under a live pool; init_storage() runs after "
       "configure_arenas()"},
      {"transformer_pinned_owner_order", case_transformer_pinned_owner_order,
       "core/animation/timeline.h",
       "(!retiring_predecessor || event.owner == nullptr || !event.pinned || event.iface->is_canceled()) retire later pinned owners before their predecessors"},
      {"transformer_pool_outlives_timeline",
       case_transformer_pool_outlives_timeline, "core/animation/transformer.h",
       "(global_timeline_live) TransformerPool outlived its Timeline: declare "
       "the Timeline before the pools that schedule on it"},
      {"effect_margin_below_pipeline", case_effect_margin_below_pipeline,
       "core/render/canvas.h",
       "(m >= required_margin) render margin is smaller than the pipeline requirement"},
      {"effect_double_construct", case_effect_double_construct,
       "core/render/canvas.h",
       "(!s_alive) Effect: a second Effect was constructed while one is "
       "still alive; buffer_a/buffer_b are shared static storage (one live "
       "Effect only)"},
      {"effect_width_zero", case_effect_width_zero, "core/render/canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 0 x "
       "16 are outside 1..288 x 1..144"},
      {"effect_height_zero", case_effect_height_zero, "core/render/canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 32 x "
       "0 are outside 1..288 x 1..144"},
      {"effect_width_over_max", case_effect_width_over_max,
       "core/render/canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 289 x "
       "16 are outside 1..288 x 1..144"},
      {"effect_height_over_max", case_effect_height_over_max,
       "core/render/canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 32 x "
       "145 are outside 1..288 x 1..144"},
      {"particle_lifetime_zero", case_particle_lifetime_zero,
       "core/animation/sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"particle_friction_nan", case_particle_friction_nan,
       "core/animation/sprites.h",
       "(std::isfinite(friction) && std::isfinite(gravity)) ParticleSystem "
       "friction and gravity must be finite"},
      {"particle_gravity_nan", case_particle_gravity_nan,
       "core/animation/sprites.h",
       "(std::isfinite(friction) && std::isfinite(gravity)) ParticleSystem "
       "friction and gravity must be finite"},
      {"particle_attractor_softening_negative",
       case_particle_attractor_softening_negative, "core/animation/sprites.h",
       "(std::isfinite(softening) && softening >= 0.0f) ParticleSystem "
       "attractor softening must be nonnegative"},
      {"particle_lifetime_nan", case_particle_lifetime_nan,
       "core/animation/sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"random_walk_nonfinite_options", case_random_walk_nonfinite_options,
       "core/animation/motion.h",
       "(std::isfinite(options.speed) && std::isfinite(options.pivot_strength) "
       "&& std::isfinite(options.noise_scale) && "
       "std::isfinite(options.smoothing) && std::isfinite(options.drift)) "
       "RandomWalk: options must be finite"},
      {"particle_lifetime_over_max", case_particle_lifetime_over_max,
       "core/animation/sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"particle_render_zero_lifetime", case_particle_render_zero_lifetime,
       "core/render/plot/particles.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem render max_life must be finite and in [1, "
       "65535]"},
      {"correction_guard_double_construct",
       case_correction_guard_double_construct, "core/platform/led.h",
       "(!correction_guard_live()) NoColorCorrection and NoTempCorrection "
       "guards cannot overlap"},
      {"correction_guard_cross_type", case_correction_guard_cross_type,
       "core/platform/led.h",
       "(!correction_guard_live()) NoColorCorrection and NoTempCorrection "
       "guards cannot overlap"},
      {"mesh_narrow_index", case_mesh_narrow_index, "core/mesh/mesh.h",
       "(i <= MeshLimits::MAX_VERTEX_INDEX) mesh index exceeds topology "
       "vertex range (oversized mesh?)"},
      {"medial_aliases_input", case_medial_aliases_input, "core/mesh/conway.h",
       "(&mesh != &out_a) medial input mesh must not alias output mesh"},
      {"needle_aliased_arenas", case_needle_aliased_arenas,
       "core/mesh/conway.h",
       "(&target != &temp) needle: target and temp must differ"},
      {"zip_aliased_arenas", case_zip_aliased_arenas, "core/mesh/conway.h",
       "(&target != &temp) zip: target and temp must differ"},
      {"gyro_aliased_arenas", case_gyro_aliased_arenas, "core/mesh/conway.h",
       "(&target != &temp) gyro: target and temp must differ"},
      {"bevel_aliased_arenas", case_bevel_aliased_arenas, "core/mesh/conway.h",
       "(&target != &temp) bevel: target and temp must differ"},
      {"ambo_aliased_arenas", case_ambo_aliased_arenas, "core/mesh/conway.h",
       "(&target != &temp) ambo: target and temp must differ"},
      {"mesh_compile_aliased_arenas", case_mesh_compile_aliased_arenas,
       "core/mesh/mesh.h",
       "(&geom_arena != &scratch) MeshOps::compile geom_arena must not alias "
       "scratch"},
      {"mesh_ops_clone_aliases_dst", case_mesh_ops_clone_aliases_dst,
       "core/mesh/mesh.h",
       "(&src != &dst) MeshOps::clone src must not alias dst"},
      {"mesh_state_clone_aliases_dst", case_mesh_state_clone_aliases_dst,
       "core/mesh/mesh_state.h",
       "(&src != &dst) MeshState::clone src must not alias "
       "dst"},
      {"mesh_transform_aliases_source", case_mesh_transform_aliases_source,
       "core/mesh/conway.h",
       "(&mesh != &transformed) MeshOps::transform source mesh must not alias "
       "the destination"},
      {"chamfer_collapsed_endpoint", case_chamfer_collapsed_endpoint,
       "core/mesh/conway.h", "(t >= 0.0f && t < 1.0f) chamfer: t out of [0,1)"},
      {"snub_collapsed_endpoint", case_snub_collapsed_endpoint,
       "core/mesh/conway.h", "(t >= 0.0f && t < 1.0f) snub: t out of [0,1)"},
      {"conway_empty_mesh", case_conway_empty_mesh, "core/mesh/mesh.h",
       "(total_indices > 0) half-edge mesh requires at least one face index"},
      {"conway_degenerate_mesh", case_conway_degenerate_mesh,
       "core/mesh/mesh.h",
       "(he_mesh.half_edges[i].pair != HE_NONE) MeshOps::truncate requires a "
       "closed manifold (unpaired half-edge)"},
      {"conway_target_exhausted", case_conway_target_exhausted,
       "core/memory.cpp", "(false) Arena::allocate: out of memory"},
      {"relax_baked_dimension_mismatch", case_relax_baked_dimension_mismatch,
       "core/mesh/conway.h",
       "(V == bake.vertex_count && F == bake.face_count && I == "
       "bake.index_count) relax_baked: source dimensions differ"},
      {"relax_baked_topology_mismatch", case_relax_baked_topology_mismatch,
       "core/mesh/conway.h",
       "(relax_topology_hash(mesh) == bake.topology_hash) relax_baked: source "
       "topology differs"},
      {"relax_baked_source_mismatch", case_relax_baked_source_mismatch,
       "core/mesh/relax_bake.h",
       "(relax_source_hash(mesh) == bake.source_hash) relax_baked: source "
       "vertices differ"},
      {"relax_baked_output_hash_mismatch",
       case_relax_baked_output_hash_mismatch, "core/mesh/conway.h",
       "(output_hash == bake.output_hash) relax_baked: output hash differs"},
      {"half_edge_face_counts_short", case_half_edge_face_counts_short,
       "core/mesh/mesh.h",
       "(counted_indices == total_indices) mesh face counts do not span flat "
       "index length"},
      {"half_edge_face_counts_long", case_half_edge_face_counts_long,
       "core/mesh/mesh.h",
       "(count <= total_indices - counted_indices) mesh face counts exceed "
       "flat index length"},
      {"mesh_compile_face_counts_short", case_mesh_compile_face_counts_short,
       "core/mesh/mesh.h",
       "(counted_indices == total_indices) mesh face counts do not span flat "
       "index length"},
      {"mesh_compile_face_counts_long", case_mesh_compile_face_counts_long,
       "core/mesh/mesh.h",
       "(count <= total_indices - counted_indices) mesh face counts exceed "
       "flat index length"},
      {"mesh_compile_face_span_over_16bit",
       case_mesh_compile_face_span_over_16bit, "core/mesh/mesh.h",
       "(static_cast<size_t>(current_offset) + count <= MeshLimits::MAX_HALF_EDGES) mesh face_offsets exceeds "
       "16-bit index range"},
      {"update_hankin_stale_topology", case_update_hankin_stale_topology,
       "core/mesh/hankin.h",
       "(prior_topology_size == 0 || out_mesh.topology_key == "
       "compiled.topology_key) update_hankin: reused out_mesh carries "
       "a topology from a different compiled pattern (clear it first)"},
      {"update_hankin_dual_seed_topology",
       case_update_hankin_dual_seed_topology, "core/mesh/hankin.h",
       "(prior_topology_size == 0 || out_mesh.topology_key == "
       "compiled.topology_key) update_hankin: reused out_mesh carries "
       "a topology from a different compiled pattern (clear it first)"},
      {"update_hankin_borrowed_stale_topology",
       case_update_hankin_borrowed_stale_topology, "core/mesh/hankin.h",
       "(prior_topology_size == 0 || out_mesh.topology_key == "
       "compiled.topology_key) update_hankin: reused out_mesh carries "
       "a topology from a different compiled pattern (clear it first)"},
      {"update_hankin_nonfinite_angle", case_update_hankin_nonfinite_angle,
       "core/mesh/hankin.h",
       "(std::isfinite(angle)) update_hankin: contact angle must be finite"},
      {"hankin_clone_aliases_dst", case_hankin_clone_aliases_dst,
       "core/mesh/hankin.h",
       "(&src != &dst) CompiledHankin::clone src must not alias dst"},
      {"mesh_state_set_borrowed_offsets_count_mismatch",
       case_mesh_state_set_borrowed_offsets_count_mismatch,
       "core/mesh/mesh_state.h",
       "(face_offsets_span.size() == face_counts_span.size()) "
       "MeshState::set_borrowed: one face offset per face count required"},
      {"mesh_state_set_borrowed_offsets_short_span",
       case_mesh_state_set_borrowed_offsets_short_span,
       "core/mesh/mesh_state.h",
       "(static_cast<size_t>(face_offsets_span[last]) + "
       "face_counts_span[last] == faces_span.size()) MeshState::set_borrowed: "
       "face offsets do not span faces"},
      {"mesh_state_set_borrowed_keyed_empty_topology",
       case_mesh_state_set_borrowed_keyed_empty_topology,
       "core/mesh/mesh_state.h",
       "(!topology_span.is_empty() || key == 0) MeshState::set_borrowed: an "
       "empty topology span requires a zero topology key"},
      {"mesh_state_set_borrowed_offsets_not_prefix_sum",
       case_mesh_state_set_borrowed_offsets_not_prefix_sum,
       "core/mesh/mesh_state.h",
       "(offsets_are_prefix_sum(face_counts_span, face_offsets_span)) "
       "MeshState::set_borrowed: face offsets are not the prefix sum of the "
       "face counts"},
      {"half_edge_zero_side_face", case_half_edge_zero_side_face,
       "core/mesh/mesh.h", "(count > 0) half-edge mesh face has zero sides"},
      {"half_edge_non_manifold_edge", case_half_edge_non_manifold_edge,
       "core/mesh/mesh.h",
       "(j - i <= 2) non-manifold edge: >2 half-edges share an edge"},
      {"half_edge_inconsistent_winding", case_half_edge_inconsistent_winding,
       "core/mesh/mesh.h",
       "(out.half_edges[a].vertex != out.half_edges[b].vertex) half-edge "
       "mesh faces are inconsistently wound"},
      {"mesh_narrow_face_count", case_mesh_narrow_face_count,
       "core/mesh/mesh.h",
       "(count >= 0 && count <= MeshLimits::MAX_FACE_DEGREE) mesh face side count exceeds "
       "uint8_t range"},
      {"mesh_require_closed_manifold", case_mesh_require_closed_manifold,
       "core/mesh/mesh.h",
       "(he_mesh.half_edges[i].pair != HE_NONE) MeshOps::death requires a "
       "closed manifold (unpaired half-edge)"},
      {"mesh_require_vertex_manifold", case_mesh_require_vertex_manifold,
       "core/mesh/mesh.h",
       "(walked == fan_size[origin]) MeshOps::death requires a vertex "
       "manifold (split vertex fan)"},
      {"mesh_require_matching_face_sides",
       case_mesh_require_matching_face_sides, "core/mesh/mesh.h",
       "(static_cast<size_t>(he_mesh.faces[fi].half_edge) == face_offset) "
       "MeshOps::death: half-edge mesh face sides differ from the source mesh"},
      {"mesh_require_matching_half_edge_census",
       case_mesh_require_matching_half_edge_census, "core/mesh/mesh.h",
       "(face_offset == he_mesh.half_edges.size()) "
       "MeshOps::death: half-edge mesh census differs from the source mesh"},
      {"mesh_require_matching_half_edge_loops",
       case_mesh_require_matching_half_edge_loops, "core/mesh/mesh.h",
       "(head == faces[face_offset + (k + 1) % count]) "
       "MeshOps::death: half-edge mesh loops different faces than the "
       "source mesh"},
      {"apply_step_hankin_no_angle", case_apply_step_hankin_no_angle,
       "core/mesh/recipe.h",
       "(step.param > 0.0f) apply_step: HANKIN step has no contact angle"},
      {"apply_step_bevel_no_depth", case_apply_step_bevel_no_depth,
       "core/mesh/recipe.h",
       "(step.param > 0.0f) apply_step: BEVEL step has no depth"},
      {"slerp_nan", case_slerp_nan, "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Vector: zero length"},
      {"make_rotation_vectors_nan", case_make_rotation_vectors_nan,
       "core/math/3dmath.h",
       "(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ && "
       "std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ) "
       "make_rotation(from, to): inputs must be unit vectors"},
      {"make_rotation_angle_nan", case_make_rotation_angle_nan,
       "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Quaternion: zero magnitude"},
      {"make_rotation_nonunit", case_make_rotation_nonunit,
       "core/math/3dmath.h",
       "(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ && "
       "std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ) "
       "make_rotation(from, to): inputs must be unit vectors"},
      {"make_basis_nan", case_make_basis_nan, "core/math/3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) Vector: zero length"},
      {"noise_transform_nan", case_noise_transform_nan,
       "core/animation/transformer.h",
       "(std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z)) "
       "noise_transform: non-finite direction"},
      {"param_def_unknown_get_target_type",
       case_param_def_unknown_get_target_type, "core/control/params.h",
       "(false) ParamDef::get_from: unknown target type "},
      {"param_def_unknown_set_target_type",
       case_param_def_unknown_set_target_type, "core/control/params.h",
       "(false) ParamDef::write_unchecked: unknown target type "},
      {"driver_null_speed_src", case_driver_null_speed_src,
       "core/animation/params.h",
       "(speed_src != nullptr) Driver: live speed_src is null"},
      {"path_append_zero_samples", case_path_append_zero_samples,
       "core/animation/motion.h",
       "(samples >= 1) Path: samples must be positive"},
      {"motion_empty_path_origin_sample", case_motion_empty_path_origin_sample,
       "core/animation/motion.h",
       "(math::dot(current_v, current_v) >= math::EPS_LEN_SQ && math::dot(target_v, "
       "target_v) >= math::EPS_LEN_SQ) Motion: path sampled at the origin "
       "(empty or origin-crossing path)"},
      {"hue_wobble_depth_out_of_range", case_hue_wobble_depth_out_of_range,
       "core/color/composition.h",
       "(fabsf(depth) * (2.0f * math::PI_F) < PALETTE_PHASE_ARG_LIMIT) "
       "HueWobbleShade: depth must stay inside the fast-trig argument range"},
      {"iridescent_weight_negative", case_iridescent_weight_negative,
       "core/color/composition.h",
       "(weight >= 0.0f) IridescentShade: weight must be non-negative"},
      {"alpha_falloff_null", case_alpha_falloff_null,
       "core/color/composition.h",
       "(fn != nullptr) AlphaFalloffShade: falloff function must not be null"},
      {"noise_hue_palette_direct_null_source",
       case_noise_hue_palette_direct_null_source,
       "core/color/noise_hue_palette.h",
       "(source != nullptr) NoiseHuePalette direct mode bound to null source"},
      {"noise_hue_palette_direct_null_noise_lut",
       case_noise_hue_palette_direct_null_noise_lut,
       "core/color/noise_hue_palette.h",
       "(hue_noise_lut != nullptr) NoiseHuePalette direct mode bound to null hue-noise LUT"},
      {"noise_shimmer_palette_null_source",
       case_noise_shimmer_palette_null_source,
       "core/color/noise_shimmer_palette.h",
       "(source != nullptr) NoiseShimmerPalette bound to null source"},
      {"noise_shimmer_palette_null_noise_lut",
       case_noise_shimmer_palette_null_noise_lut,
       "core/color/noise_shimmer_palette.h",
       "(noise_lut != nullptr) NoiseShimmerPalette bound to null noise LUT"},
      {"noise_hue_palette_null_source", case_noise_hue_palette_null_source,
       "core/color/noise_hue_palette.h",
       "(source != nullptr) NoiseHuePalette bound to null source"},
      {"noise_hue_palette_null_rotation_lut",
       case_noise_hue_palette_null_rotation_lut,
       "core/color/noise_hue_palette.h",
       "(hue_rotation_lut != nullptr) NoiseHuePalette bound to null "
       "hue-rotation LUT"},
      {"noise_hue_palette_null_noise_lut",
       case_noise_hue_palette_null_noise_lut, "core/color/noise_hue_palette.h",
       "(hue_noise_lut != nullptr) NoiseHuePalette bound to null hue-noise "
       "LUT"},
      {"palette_cycler_mutated_policy", case_palette_cycler_mutated_policy,
       "core/color/palette_cycler.h",
       "(entries[current].generative->morph_compatible( "
       "*entries[next_of(current)].generative)) "
       "PaletteCycler borrowed palette policy changed after init"},
      {"generated_palette_bank_unknown_mode",
       case_generated_palette_bank_unknown_mode, "core/color/palette_cycler.h",
       "(false) GeneratedPaletteBank::palette: unknown palette mode"},
      {"palette_cycler_restore_without_generated_init",
       case_palette_cycler_restore_without_generated_init,
       "core/color/palette_cycler.h",
       "(from_slot != nullptr && to_slot != nullptr && morph != nullptr) "
       "PaletteCycler restore needs a generated cycle"},
      {"baked_palette_clone_from_self", case_baked_palette_clone_from_self,
       "core/color/baked_palette.h",
       "(&src != &table) BakedPaletteStorage::clone_from from itself"},
      {"baked_palette_bake_blend_self", case_baked_palette_bake_blend_self,
       "core/color/baked_palette.h",
       "(&from != &table && &to != &table) BakedPaletteStorage::bake_blend endpoint is "
       "the output"},
      {"float_options_missing_labels", case_float_options_missing_labels,
       "core/control/param_host.h",
       "((options == nullptr) == (option_count == 0)) "
       "register_param: inconsistent options and count"},
      {"float_options_missing_count", case_float_options_missing_count,
       "core/control/param_host.h",
       "((options == nullptr) == (option_count == 0)) "
       "register_param: inconsistent options and count"},
      {"float_options_wrong_range", case_float_options_wrong_range,
       "core/control/param_host.h",
       "(option_count > 0 && min == 0 && max == option_count - 1) "
       "register_param: option range does not match labels"},
      {"register_param_overflow", case_register_param_overflow,
       "core/control/param_host.h",
       "(parameters.count < parameters.capacity()) register_param: "
       "exceeded ParamList capacity"},
      {"wasm_param_capacity_exceeded", case_wasm_param_capacity_exceeded,
       "targets/wasm/param_marshal.h",
       "(count <= ParamStreams::CAPACITY) effect exposes 257 params, past the "
       "256 reserved"},
      {"register_param_duplicate", case_register_param_duplicate,
       "core/control/param_host.h",
       "(parameters.find(name) == nullptr) register_param: duplicate parameter name name=duplicate"},
      {"register_param_default_outside_range",
       case_register_param_default_outside_range, "core/control/param_host.h",
       "(*ptr >= min && *ptr <= max) register_param: default *ptr outside [min,max] name=outside value_bits=40000000 min_bits=00000000 max_bits=3f800000"},
      {"restore_parameters_unknown_name", case_restore_parameters_unknown_name,
       "core/control/param_host.h",
       "(def != nullptr && !def->readonly) replay_parameter_writes: unknown or readonly parameter"},
      {"restore_parameters_readonly_name",
       case_restore_parameters_readonly_name, "core/control/param_host.h",
       "(def != nullptr && !def->readonly) replay_parameter_writes: unknown or readonly parameter"},
      {"restore_parameters_singular_mobius",
       case_restore_parameters_singular_mobius, "core/control/param_host.h",
       "(parameter_write_admitted(def, def.get_requested())) replay_parameter_writes: inadmissible final state"},
      {"register_enum_param_range", case_register_enum_param_range,
       "core/control/param_host.h",
       "(static_cast<int64_t>(option_count - 1) <= "
       "static_cast<int64_t>(std::numeric_limits<Integer>::max())) "
       "register_param: options must fit the target enum type"},
      {"register_enum_param_bound_inexact",
       case_register_enum_param_bound_inexact, "core/control/param_host.h",
       "(static_cast<int64_t>(static_cast<float>(option_count - 1)) == "
       "option_count - 1) register_param: enum bound must be exactly representable as float"},
      {"register_int_param_range", case_register_int_param_range,
       "core/control/param_host.h",
       "(range_fits) register_param: [min,max] must fit the target "
       "integer type"},
      {"register_int_param_max_inexact", case_register_int_param_max_inexact,
       "core/control/param_host.h",
       "(bounds_exact) register_param: bounds must be exactly representable as float"},
      {"register_int_param_min_inexact", case_register_int_param_min_inexact,
       "core/control/param_host.h",
       "(bounds_exact) register_param: bounds must be exactly representable as float"},
      {"param_spec_invalid_option_values",
       case_param_spec_invalid_option_values, "core/control/param_host.h",
       "(spec.valid_option_values(*ptr)) register_param: invalid explicit option values"},
      {"param_spec_uint32_bound_outside_storage",
       case_param_spec_uint32_bound_outside_storage,
       "core/control/param_host.h",
       "(range_fits) register_param: [min,max] must fit the target integer type"},
      {"param_spec_uint32_bound_inexact", case_param_spec_uint32_bound_inexact,
       "core/control/param_host.h",
       "(bounds_exact) register_param: bounds must be exactly representable as float"},
      {"param_spec_integer_preserve_policy",
       case_param_spec_integer_preserve_policy, "core/control/param_host.h",
       "(spec.initial_value == ParamInitialValue::REQUIRE_IN_RANGE) register_param: "
       "preserving requested values requires a non-enum float"},
      {"param_spec_float_bound_nonfinite",
       case_param_spec_float_bound_nonfinite, "core/control/param_host.h",
       "(std::isfinite(min) && std::isfinite(max)) register_param: bounds must be finite"},
      {"param_spec_name_null", case_param_spec_name_null,
       "core/control/param_host.h",
       "(name != nullptr) register_param: null parameter name"},
      {"param_spec_requested_nonfinite", case_param_spec_requested_nonfinite,
       "core/control/param_host.h",
       "(std::isfinite(*ptr)) register_param: requested value must be finite"},
      {"param_spec_option_label_null", case_param_spec_option_label_null,
       "core/control/param_host.h",
       "(options[i] != nullptr) register_param: null option label"},
      {"param_spec_export_label_null", case_param_spec_export_label_null,
       "core/control/param_host.h",
       "(spec.export_options == nullptr || spec.export_options[i] != nullptr) "
       "register_param: null export label"},
      {"set_clip_out_of_bounds", case_set_clip_out_of_bounds,
       "core/render/canvas.h",
       "(y0 >= 0 && y0 <= y1 && y1 <= clip_region.h && x0 >= 0 && x0 <= x1 "
       "&& x1 <= clip_region.w) set_clip band must be non-inverted and "
       "within canvas bounds"},
      {"set_clip_mid_frame", case_set_clip_mid_frame, "core/render/canvas.h",
       "(!canvas_active) clip cannot change while a frame is active"},
      {"output_envelope_out_of_range", case_output_envelope_out_of_range,
       "core/render/canvas.h",
       "(std::isfinite(value) && value >= 0.0f && value <= 1.0f) output "
       "envelope must be finite and in [0,1]"},
      {"set_clip_x_out_of_bounds", case_set_clip_x_out_of_bounds,
       "core/render/canvas.h",
       "(x0 >= 0 && x0 <= x1 && x1 <= clip_region.w) set_clip_x band must be "
       "non-inverted and within canvas width"},
      {"arcs_overlap_start_out_of_range", case_arcs_overlap_start_out_of_range,
       "core/render/clip.h", "(s1 >= 0 && s1 < w && s2 >= 0 && s2 < w) "},
      {"scan_block_coherent_zero_block", case_scan_block_coherent_zero_block,
       "core/render/scan/shader.h",
       "(block > 0) block-coherent shader requires a positive block size"},
      {"scan_clip_out_of_bounds", case_scan_clip_out_of_bounds,
       "core/render/scan/shader.h",
       "(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 && "
       "cr.render_y_end() <= H) scan clip region outside the canvas and trig "
       "LUT domain"},
      {"scan_clip_rows_out_of_bounds", case_scan_clip_rows_out_of_bounds,
       "core/render/scan/shader.h",
       "(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 && "
       "cr.render_y_end() <= H) scan clip region outside the canvas and trig "
       "LUT domain"},
      {"face_scratch_retargeted", case_face_scratch_retargeted,
       "core/render/sdf/face_geometry.h",
       "(!scratch_owner || scratch_owner->claim_seq == scratch_claim) "
       "SDF::Face scanned after a later Face claimed its scratch buffer"},
      {"face_virtual_height_mismatched_geometry",
       case_face_virtual_height_mismatched_geometry, "core/render/sdf/face.h",
       "(build_geometry.row_to_phi(0) == DISPLAY_GEOMETRY.row_to_phi(0) && "
       "build_geometry.row_to_phi(height - 1) == "
       "DISPLAY_GEOMETRY.row_to_phi(height - 1)) "
       "Face: virtual height must match the display geometry"},
      {"face_scratch_retargeted_by_culled_face",
       case_face_scratch_retargeted_by_culled_face,
       "core/render/sdf/face_geometry.h",
       "(!scratch_owner || scratch_owner->claim_seq == scratch_claim) "
       "SDF::Face scanned after a later Face claimed its scratch buffer"},
      {"scan_canvas_dim_mismatch", case_scan_canvas_dim_mismatch,
       "core/render/scan/raster.h",
       "(canvas.width() == W && canvas.height() == H) canvas size differs from "
       "the scan's W/H"},
      {"plot_canvas_dim_mismatch", case_plot_canvas_dim_mismatch,
       "core/render/plot/raster.h",
       "(canvas.width() == W && canvas.height() == H) canvas size differs from "
       "the plot's W/H"},
      {"plot_window_multi_segment", case_plot_window_multi_segment,
       "core/render/plot/raster.h",
       "(!plot_window || count == 1) a plot window requires a single-segment "
       "polyline"},
      {"planar_chords_sink_unprepared", case_planar_chords_sink_unprepared,
       "core/render/plot/chords.h",
       "(pipeline.prepared_for(canvas)) direct raster pipeline not prepared for this canvas"},
      {"scan_pipeline_not_prepared", case_scan_pipeline_not_prepared,
       "core/render/scan/raster.h",
       "(pipeline.prepared_for(canvas)) direct raster pipeline not prepared "
       "for this canvas"},
      {"pipeline_ref_erase_not_prepared", case_pipeline_ref_erase_not_prepared,
       "core/engine/concepts.h",
       "(t.prepared_for(cv)) direct raster pipeline not prepared for this "
       "canvas"},
      {"scan_mesh_face_index_out_of_range",
       case_scan_mesh_face_index_out_of_range, "core/render/scan/mesh.h",
       "(static_cast<size_t>(faces[k]) < num_verts) mesh face index exceeds "
       "the vertex pool"},
      {"scan_mesh_missing_offsets", case_scan_mesh_missing_offsets,
       "core/render/scan/mesh.h",
       "(mesh.get_face_offsets_size() == num_f) solid mesh scan requires one "
       "face offset per face"},
      {"scan_mesh_class_id_out_of_range", case_scan_mesh_class_id_out_of_range,
       "core/render/scan/mesh.h",
       "(rec.class_id < bake->classes.size()) mesh class bake face record "
       "names an unknown class"},
      {"plot_mesh_vertex_over_capacity", case_plot_mesh_vertex_over_capacity,
       "core/render/plot/mesh.h",
       "(large < DEDUP_CAPACITY) Mesh edge dedup: vertex index "},
      {"plot_four_regular_open_mesh", case_plot_four_regular_open_mesh,
       "core/render/plot/mesh.h",
       "(he.pair != HE_NONE) extract_four_regular_edges: mesh is not closed"},
      {"plot_four_regular_non_bipartite", case_plot_four_regular_non_bipartite,
       "core/render/plot/mesh.h",
       "(colors[adjacent] == adjacent_color) four-regular mesh dual is not bipartite"},
      {"plot_find_missing_edge", case_plot_find_missing_edge,
       "core/render/plot/mesh.h",
       "(false) find_edge_index: face edge missing from the edge list"},
      {"plot_extract_edges_vertex_over_capacity",
       case_plot_extract_edges_vertex_over_capacity, "core/render/plot/mesh.h",
       "(large < DEDUP_CAPACITY) Mesh edge dedup: vertex index "},
      {"feedback_downsample_indivisible", case_feedback_downsample_indivisible,
       "core/render/filter/pixel_feedback.h",
       "(downsample > 0 && W % downsample == 0) feedback downsample 5 must "
       "be > 0 and divide width 32"},
      {"feedback_uncached_scratch_budget",
       case_feedback_uncached_scratch_budget,
       "core/render/filter/pixel_feedback.h",
       "(uncached_scratch_bytes(grid.downsample) <= scratch.get_capacity() - scratch.get_offset()) uncached feedback needs more scratch:"},
      {"screen_trails_set_lifetime_nonpositive",
       case_screen_trails_set_lifetime_nonpositive,
       "core/render/filter/screen_trails.h",
       "(new_lifetime > 0) Screen::Trails: lifetime 0 must be positive"},
      {"world_trails_lifetime_over_max", case_world_trails_lifetime_over_max,
       "core/render/filter/world_trails.h",
       "(lifetime > 0 && lifetime <= 255) World::Trails: lifetime 256 outside "
       "[1, 255]"},
      {"screen_trails_plot_without_storage",
       case_screen_trails_plot_without_storage,
       "core/render/filter/screen_trails.h",
       "(points) Screen::Trails needs init_storage() from effect init()"},
      {"world_trails_plot_without_storage",
       case_world_trails_plot_without_storage,
       "core/render/filter/world_trails.h",
       "(items) World::Trails needs init_storage() from effect init()"},
      {"plot_open_loop_seam", case_plot_open_loop_seam,
       "core/render/plot/shapes.h",
       "(params.loop_seam == nullptr || params.close_loop) a raster seam fragment requires a closed loop"},
      {"raster_point_projection_pair_mismatch",
       case_raster_point_projection_pair_mismatch, "core/render/plot/raster.h",
       "(rows.size() == cols.size()) hoisted point projection rows and columns differ in length"},
      {"raster_nonfinite_point", case_raster_nonfinite_point,
       "core/render/plot/raster.h",
       "(std::isfinite(point.pos.x) && std::isfinite(point.pos.y) && std::isfinite(point.pos.z)) rasterize control points must be finite"},
      {"raster_single_nonfinite_point", case_raster_single_nonfinite_point,
       "core/render/plot/raster.h",
       "(std::isfinite(point.pos.x) && std::isfinite(point.pos.y) && std::isfinite(point.pos.z)) rasterize control points must be finite"},
      {"raster_empty_path_null_shader", case_raster_empty_path_null_shader,
       "core/render/plot/raster.h",
       "(fragment_shader) rasterize requires a non-null fragment_shader"},
      {"raster_edge_flags_short", case_raster_edge_flags_short,
       "core/render/plot/raster.h",
       "(edge_flags == nullptr || opts.projection.flags().size() == count) edge_flags length must match the rasterized edge count"},
      {"raster_point_projections_short", case_raster_point_projections_short,
       "core/render/plot/raster.h",
       "(point_rows == nullptr || opts.point_projections.size() == len) hoisted "
       "point projections need one entry per polyline point"},
      {"spherical_field_ring_index_oob", case_spherical_field_ring_index_oob,
       "core/math/spherical_field.h",
       "(y < H - 1) SphericalFieldLayout: ring index 5 out of range"},
      {"spherical_field_populate_ring_end_oob",
       case_spherical_field_populate_ring_end_oob,
       "core/math/spherical_field.h",
       "(ring_end < layout.ring_count()) SphericalField::populate: ring_end 5 "
       "past the last ring 4"},
      {"spherical_harmonic_order_over_degree",
       case_spherical_harmonic_order_over_degree,
       "core/math/spherical_harmonics.h",
       "(l >= 0 && abs_m <= l) spherical harmonic normalization: order 3 is "
       "outside [-2, 2]"},
      {"spherical_harmonic_decode_negative_index",
       case_spherical_harmonic_decode_negative_index,
       "core/math/spherical_harmonics.h",
       "(idx >= 0) decode_lm: flat index -1 is negative"},
      {"spherical_field_infill_over_domain",
       case_spherical_field_infill_over_domain, "core/math/spherical_field.h",
       "(north_infill >= 0 && south_infill >= 0 && north_infill + south_infill "
       "<= H) SphericalFieldLayout: infills 0 + 17 must be non-negative and "
       "fit within H = 16"},
      {"spherical_field_negative_equator_samples",
       case_spherical_field_negative_equator_samples,
       "core/math/spherical_field.h",
       "(equator_samples >= 0) SphericalFieldLayout: equator_samples -1 must "
       "be >= 0"},
      {"feedback_negative_fade", case_feedback_negative_fade,
       "core/render/filter/feedback_style.h",
       "(fade >= 0.0f && fade <= 1.0f) Feedback::Style::fade must be in "
       "[0, 1]"},
      {"feedback_infinite_fade", case_feedback_infinite_fade,
       "core/render/filter/feedback_style.h",
       "(fade >= 0.0f && fade <= 1.0f) Feedback::Style::fade must be in "
       "[0, 1]"},
      {"gradient_no_stops", case_gradient_no_stops,
       "core/color/palette_sources.h",
       "(points.size() > 0) Gradient requires at least one stop"},
      {"gradient_stop_out_of_range", case_gradient_stop_out_of_range,
       "core/color/palette_sources.h",
       "(stop.first >= 0.0f && stop.first <= 1.0f) Gradient stop position "
       "out of [0,1]"},
      {"gradient_stops_unsorted", case_gradient_stops_unsorted,
       "core/color/palette_sources.h",
       "(stop.first >= prev_check) Gradient stops must be sorted ascending"},
      {"random_timer_inverted_range", case_random_timer_inverted_range,
       "core/animation/timers.h",
       "(min >= 0 && min <= max) RandomTimer: invalid frame range"},
      {"empty_fn_call", case_empty_fn_call, "core/engine/static_storage.cpp",
       "(vtable != empty) empty hs::inplace_function called"},
      {"empty_function_ref_call", case_empty_function_ref_call,
       "core/engine/static_storage.cpp",
       "(thunk != empty_thunk) empty FunctionRef called"},
      {"effect_registry_duplicate_name", case_effect_registry_duplicate_name,
       "core/control/registry.h",
       "(registration_names_unique(entries)) duplicate effect registration identity"},
      {"effect_registry_duplicate_stable_id",
       case_effect_registry_duplicate_stable_id, "core/control/registry.h",
       "(registration_names_unique(entries)) duplicate effect registration identity"},
      {"effect_registry_stable_id_matches_name",
       case_effect_registry_stable_id_matches_name, "core/control/registry.h",
       "(registration_names_unique(entries)) duplicate effect registration identity"},
      {"effect_registry_name_matches_stable_id",
       case_effect_registry_name_matches_stable_id, "core/control/registry.h",
       "(registration_names_unique(entries)) duplicate effect registration identity"},
      {"flywheel_period_zero", case_flywheel_period_zero,
       "hardware/pov_sync_flywheel.h",
       "(p > 0 && p <= static_cast<uint32_t>(INT32_MAX) / MIN_SAFE_HALF_REVS) "
       "Flywheel: cycles_per_half_rev outside the range position()'s int32 "
       "elapsed window holds for MIN_SAFE_HALF_REVS of coast"},
      {"latitude_geometry_degenerate_height",
       case_latitude_geometry_degenerate_height, "core/math/display_geometry.h",
       "(height > 1 && north >= 0.0f && south <= PI_F && north < south) "
       "Invalid latitude geometry"},
      {"latitude_geometry_reversed_span", case_latitude_geometry_reversed_span,
       "core/math/display_geometry.h",
       "(height > 1 && north >= 0.0f && south <= PI_F && north < south) "
       "Invalid latitude geometry"},
      {"y_to_phi_degenerate_height", case_y_to_phi_degenerate_height,
       "core/math/pixel_mapping.h",
       "(h_virt > 1) y_to_phi_virtual: h_virt must be > 1"},
      {"orientation_frame_index_oob", case_orientation_frame_index_oob,
       "core/animation/orientation.h",
       "(i >= 0 && i < num_frames) Orientation: frame index out of range"},
      {"make_basis_nonunit_quaternion", case_make_basis_nonunit_quaternion,
       "core/math/spherical.h",
       "(std::abs(orientation_norm_sq - 1.0f) < "
       "math::EPS_UNIT_QUAT_SQ) make_basis: orientation |q|^2 is "},
      {"parallel_transport_antipodal", case_parallel_transport_antipodal,
       "core/math/spherical.h",
       "(denominator > 1.0f || dot(cross(from, to), cross(from, to)) > "
       "MIN_TRANSPORT_CROSS_SQ) parallel_transport: antipodal endpoints"},
      {"polyhedral_kaleidoscope_no_converge",
       case_polyhedral_kaleidoscope_no_converge, "core/math/lenses.h",
       "(false) polyhedral kaleidoscope fold did not converge"},
      {"framework_invalid_geometry", case_framework_invalid_geometry,
       "core/render/sdf/framework.h",
       "(valid) framework event geometry must be valid"},
      {"framework_invalid_ray", case_framework_invalid_ray,
       "core/render/sdf/framework.h",
       "(ray.valid()) framework event ray must be valid"},
      {"octet_invalid_geometry", case_octet_invalid_geometry,
       "core/render/sdf/framework.h",
       "(geometry.valid()) octet event geometry must be valid"},
      {"octet_invalid_ray", case_octet_invalid_ray,
       "core/render/sdf/framework.h",
       "(ray.valid()) octet event ray must be valid"},
      {"octet4_invalid_geometry", case_octet4_invalid_geometry,
       "core/render/sdf/framework.h",
       "(geometry.valid()) octet4 event geometry must be valid"},
      {"octet4_invalid_interval", case_octet4_invalid_interval,
       "core/render/sdf/framework.h",
       "(interval.valid()) octet4 event interval must be valid"},
      {"octet4_nonfinite_ray", case_octet4_nonfinite_ray,
       "core/render/sdf/framework.h",
       "(Raycast::finite(origin[i]) && Raycast::finite(direction[i])) octet4 event ray components must be finite"},
      {"octet4_nonunit_direction", case_octet4_nonunit_direction,
       "core/render/sdf/framework.h",
       "(fabsf(length2 - 1.0f) < 1e-4f) octet4 event direction must be unit length"},
      {"sdf_polygon_side_count", case_sdf_polygon_side_count,
       "core/render/sdf/shapes.h",
       "(sides >= 3) SDF PlanarPolygon: sides must be at least 3"},
      {"sdf_spherical_polygon_radius_over_hemisphere",
       case_sdf_spherical_polygon_radius_over_hemisphere,
       "core/render/sdf/shapes.h",
       "(radius <= 1.0f) SDF SphericalPolygon: radius exceeds a hemisphere (build inverted)"},
      {"sdf_flower_radius_over_hemisphere",
       case_sdf_flower_radius_over_hemisphere, "core/render/sdf/shapes.h",
       "(radius <= 1.0f) SDF Flower: radius exceeds a hemisphere (build inverted)"},
      {"sdf_flower_zero_radius", case_sdf_flower_zero_radius,
       "core/render/sdf/shapes.h",
       "(radius > 0.0f) SDF Flower: radius must be positive"},
      {"sdf_class_lut_too_few_vertices", case_sdf_class_lut_too_few_vertices,
       "core/render/sdf/face_lut.h",
       "(count >= 3) build_canonical_distance_lut requires at least 3 polygon "
       "vertices"},
      {"sdf_class_lut_grid_too_small", case_sdf_class_lut_grid_too_small,
       "core/render/sdf/face_lut.h",
       "(n >= 2) build_canonical_distance_lut requires a grid resolution of at "
       "least 2"},
      {"sdf_bind_class_lut_offset_out_of_range",
       case_sdf_bind_class_lut_offset_out_of_range,
       "core/render/sdf/face_geometry.h",
       "(vert_offset >= 0 && vert_offset < count) bind_class_lut: vertex "
       "offset outside the face"},
      {"star_mismatched_chart", case_star_mismatched_chart,
       "core/render/plot/shapes.h",
       "(CHART_MATCHES) Star: edge chart does not match inputs"},
      {"sdf_distorted_ring_negative_distortion",
       case_sdf_distorted_ring_negative_distortion, "core/render/sdf/rings.h",
       "(md >= 0.0f) DistortedRing: negative maximum distortion"},
      {"chain_overlapping_storage", case_chain_overlapping_storage,
       "core/render/pullback/interpreter.h",
       "(A < B ? B - A >= block_capacity : A - B >= block_capacity) "
       "ChainProgram::bind_storage: overlapping blocks"},
      {"chain_identical_storage", case_chain_identical_storage,
       "core/render/pullback/interpreter.h",
       "(A < B ? B - A >= block_capacity : A - B >= block_capacity) "
       "ChainProgram::bind_storage: overlapping blocks"},
      {"sdf_distorted_ring_null_shift", case_sdf_distorted_ring_null_shift,
       "core/render/sdf/rings.h",
       "(sf) DistortedRing: shift_fn must be non-null"},
      {"sdf_ring_radius_past_antipode", case_sdf_ring_radius_past_antipode,
       "core/render/sdf/rings.h",
       "(radius >= 0.0f && radius <= 2.0f) Ring: radius outside [0, 2]"},
      {"sdf_ring_negative_thickness", case_sdf_ring_negative_thickness,
       "core/render/sdf/rings.h",
       "(thickness >= 0.0f) Ring: negative stroke half-width"},
      {"sdf_distorted_ring_radius_past_antipode",
       case_sdf_distorted_ring_radius_past_antipode, "core/render/sdf/rings.h",
       "(radius >= 0.0f && radius <= 2.0f) DistortedRing: radius outside "
       "[0, 2]"},
      {"sdf_distorted_ring_negative_thickness",
       case_sdf_distorted_ring_negative_thickness, "core/render/sdf/rings.h",
       "(thickness >= 0.0f) DistortedRing: negative stroke half-width"},
      {"sdf_line_negative_thickness", case_sdf_line_negative_thickness,
       "core/render/sdf/shapes.h",
       "(thickness >= 0.0f) Line: negative stroke half-width"},
      {"chain_zero_alignment",
       case_chain_zero_alignment, "core/render/pullback/interpreter.h", "(layout.align != 0 && (layout.align & (layout.align - 1)) == 0 && layout.align <= alignof(std::max_align_t)) ChainProgram::bind_storage: invalid block alignment"},
      {"chain_non_power_alignment", case_chain_non_power_alignment,
       "core/render/pullback/interpreter.h",
       "(layout.align != 0 && (layout.align & (layout.align - 1)) == 0 && layout.align <= alignof(std::max_align_t)) ChainProgram::bind_storage: invalid block alignment"},
      {"chain_overaligned_block", case_chain_overaligned_block,
       "core/render/pullback/interpreter.h",
       "(layout.align != 0 && (layout.align & (layout.align - 1)) == 0 && layout.align <= alignof(std::max_align_t)) ChainProgram::bind_storage: invalid block alignment"},
      {"chain_zero_size", case_chain_zero_size,
       "core/render/pullback/interpreter.h",
       "(layout.size > 0 && layout.size % layout.align == 0) ChainProgram::bind_storage: invalid block size"},
      {"chain_misaligned_size", case_chain_misaligned_size,
       "core/render/pullback/interpreter.h",
       "(layout.size > 0 && layout.size % layout.align == 0) ChainProgram::bind_storage: invalid block size"},
      {"chain_capacity_overflow", case_chain_capacity_overflow,
       "core/render/pullback/interpreter.h",
       "(block_capacity < std::numeric_limits<uint32_t>::max()) ChainProgram::bind_storage: capacity exceeds offset range"},
      {"chain_table_rank_decreases", case_chain_table_rank_decreases,
       "core/render/pullback/interpreter.h",
       "(entry.input <= entry.output) ChainProgram::bind_storage: operator "
       "family rank decreases"},

      {"dreamballs_woven_owner_vertex_oob",
       case_dreamballs_woven_owner_vertex_oob, "effects/DreamBalls.h",
       "(vertex < vertex_count) DreamBalls: woven edge start vertex "
       "outside the owner table"},
      {"dreamballs_woven_owner_edge_oob", case_dreamballs_woven_owner_edge_oob,
       "effects/DreamBalls.h",
       "(edge_index < edges.size()) DreamBalls: woven owner query edge index "
       "out of range"},
      {"raymarch_placement_solid_oob", case_raymarch_placement_solid_oob,
       "effects/Raymarch.h",
       "(placement_index < PLACEMENT_SOLID_COUNT) Raymarch placement solid is "
       "out of range"},
      {"spherical_harmonics_invalid_morph_mode",
       case_spherical_harmonics_invalid_morph_mode,
       "effects/SphericalHarmonics.h",
       "(synchronized) SphericalHarmonics preset synchronization failed"},
      {"hankinsolids_missing_topology", case_hankinsolids_missing_topology,
       "effects/HankinSolids.h",
       "(topology_faces == rotated_mesh.num_faces()) Hankin topology must cover every face"},
      {"islamicstars_build_budget", case_islamicstars_build_budget,
       "core/animation/recipe_build.h",
       "(persistent_arena.get_offset() <= device_persistent_budget) "
       "RecipeBuild: build leg exceeds the device persistent budget"},
      {"islamicstars_bridge_continuation",
       case_islamicstars_bridge_continuation, "core/animation/recipe_build.h",
       "(done == BuildContinuation::FINISH || "
       "done == BuildContinuation::DT_AFTER_BRIDGE || "
       "done == BuildContinuation::DTD_AFTER_BRIDGE1 || "
       "done == BuildContinuation::DTD_AFTER_BRIDGE2) "
       "RecipeBuild: invalid dual bridge continuation"},
      {"islamicstars_hankin_eager_endpoint",
       case_islamicstars_hankin_eager_endpoint, "core/animation/recipe_build.h",
       "(false) RecipeBuild: step builds no eager endpoint"},
      {"reconcile_vertices_size_mismatch",
       case_reconcile_vertices_size_mismatch, "core/mesh/conway.h",
       "(authored.vertices.size() == V) reconcile_vertices: endpoints differ "
       "in vertex count"},
      {"reconcile_vertices_empty", case_reconcile_vertices_empty,
       "core/mesh/conway.h",
       "(V > 0) reconcile_vertices: endpoint pair has no vertices"},
      {"reconcile_vertices_aliased_output",
       case_reconcile_vertices_aliased_output, "core/mesh/conway.h",
       "(&out != &identity && &out != &authored) reconcile_vertices: output must not alias either input"},
      {"sdf_angular_repeat_nonunit_axis", case_sdf_angular_repeat_nonunit_axis,
       "core/render/sdf/csg.h",
       "(fabsf(ax.length() - 1.0f) < 1e-3f) SDF CSG: repetition axis must be unit length"},
      {"sdf_distorted_ring_zero_knots", case_sdf_distorted_ring_zero_knots,
       "core/render/sdf/rings.h",
       "(kn != nullptr && n >= 3) DistortedRing: knot storage requires at least three knots"},
      {"sdf_distorted_ring_one_knots", case_sdf_distorted_ring_one_knots,
       "core/render/sdf/rings.h",
       "(kn != nullptr && n >= 3) DistortedRing: knot storage requires at least three knots"},
      {"sdf_distorted_ring_two_knots", case_sdf_distorted_ring_two_knots,
       "core/render/sdf/rings.h",
       "(kn != nullptr && n >= 3) DistortedRing: knot storage requires at least three knots"},
      {"gamut_lut_scratch_a", case_gamut_lut_scratch_a,
       "core/color/color_space.h",
       "(&arena != &scratch_arena_a && &arena != &scratch_arena_b) init_gamut_lut: global LUT requires non-scratch storage"},
      {"gamut_lut_scratch_b", case_gamut_lut_scratch_b,
       "core/color/color_space.h",
       "(&arena != &scratch_arena_a && &arena != &scratch_arena_b) init_gamut_lut: global LUT requires non-scratch storage"},
      {"scan_ring_stack_too_many_slots", case_scan_ring_stack_too_many_slots,
       "core/render/scan/shapes.h",
       "(n_slots <= INT8_MAX) ring stack exceeds the signed slot index range"},
      {"scan_ring_stack_too_many_rings", case_scan_ring_stack_too_many_rings,
       "core/render/scan/shapes.h",
       "(n_rings <= Table::MAX_RINGS) ring stack exceeds the candidate table's ring index range"},
      {"scan_ring_stack_callback_ring", case_scan_ring_stack_callback_ring,
       "core/render/scan/shapes.h",
       "(shapes[s].knots != nullptr) ring stack rings must be knot rings"},
      {"sdf_twist_zero_major_radius", case_sdf_twist_zero_major_radius,
       "core/render/sdf/volume.h",
       "(R > 0.0f) SDF Volume: radius must be positive"},
      {"transformed_torus_invalid_minor_radius",
       case_transformed_torus_invalid_minor_radius, "core/render/sdf/volume.h",
       "(base.r <= base.R * 0.5f) WarpedVolume Torus minor radius must not exceed half its major radius"},
      {"pick_next_edge_unknown_node", case_pick_next_edge_unknown_node,
       "core/mesh/conway_graph.h",
       "(n > 0) pick_next_edge: node outside the graph"},
      {"opleg_edge_sweep_no_edge", case_opleg_edge_sweep_no_edge,
       "core/animation/opleg.h",
       "(spec.edge) OpLeg: edge sweep carries no graph edge"},
      {"opleg_rewind_refill", case_opleg_rewind_refill,
       "core/animation/opleg.h",
       "(stamp.block_alive(buf, live_bytes)) OpLeg: leg arena storage reclaimed under a live leg",
       true},
      {"opleg_zero_sweep_frames", case_opleg_zero_sweep_frames,
       "core/animation/opleg.h",
       "(spec.sweep_frames >= 1) OpLeg: parameter sweep needs a positive sweep length"},
      {"opleg_incomplete_palette_handoff",
       case_opleg_incomplete_palette_handoff, "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_edge_settle_mismatch", case_opleg_edge_settle_mismatch,
       "core/animation/opleg.h",
       "(spec.settle_frames >= 0 && edge.settle == (spec.settle_frames > 0)) "
       "OpLeg: settle frames disagree with the edge"},
      {"opleg_edge_sweep_incomplete_handoff",
       case_opleg_edge_sweep_incomplete_handoff, "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_hankin_incomplete_handoff", case_opleg_hankin_incomplete_handoff,
       "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_hankin_backward_theta", case_opleg_hankin_backward_theta,
       "core/animation/opleg.h",
       "(spec.theta_start <= spec.theta_end) OpLeg: hankin leg sweeps back to "
       "a smaller contact angle"},
      {"opleg_relax_no_iterations", case_opleg_relax_no_iterations,
       "core/animation/opleg.h",
       "(spec.bake || spec.iterations >= 1) OpLeg: relax leg needs a positive "
       "iteration count"},
      {"opleg_relax_incomplete_handoff", case_opleg_relax_incomplete_handoff,
       "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_medial_incomplete_handoff", case_opleg_medial_incomplete_handoff,
       "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_reconcile_endpoint_count", case_opleg_reconcile_endpoint_count,
       "core/animation/opleg.h",
       "(spec.to_count == from_mesh.vertices.size()) "
       "OpLeg: reconcile endpoint count must match seed vertices"},
      {"opleg_reconcile_no_endpoints", case_opleg_reconcile_no_endpoints,
       "core/animation/opleg.h",
       "(spec.to_positions) OpLeg: reconcile leg carries no "
       "endpoints"},
      {"opleg_reconcile_incomplete_handoff",
       case_opleg_reconcile_incomplete_handoff, "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_gated_swap_zero_gate_frames",
       case_opleg_gated_swap_zero_gate_frames, "core/animation/opleg.h",
       "(spec.gate_frames >= 1) OpLeg needs a positive gate length"},
      {"opleg_gated_swap_incomplete_handoff",
       case_opleg_gated_swap_incomplete_handoff, "core/animation/opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: incomplete palette handoff"},
      {"opleg_shading_face_out_of_range", case_opleg_shading_face_out_of_range,
       "core/animation/opleg.h",
       "(face < faces) OpLeg::Shading: ramp face out of range"},
      {"pullback_project_nonunit_direction",
       case_pullback_project_nonunit_direction, "core/render/pullback/stage.h",
       "(fabsf(input.dir.x * input.dir.x + input.dir.y * input.dir.y + "
       "input.dir.z * input.dir.z - 1.0f) <= 0.004f) "
       "Project requires a unit direction within lens approximation error"},
      {"pullback_operator_invalid_coverage_mode",
       case_pullback_operator_invalid_coverage_mode,
       "core/render/pullback/operators/sample.h",
       "(coverage_mode <= static_cast<uint8_t>("
       "ProjectionCoverageMode::EDGE_FADE)) "
       "sample operator: invalid projection coverage mode"},
      {"pullback_operator_invalid_weight_mode",
       case_pullback_operator_invalid_weight_mode,
       "core/render/pullback/operators/sample.h",
       "(weight_mode <= static_cast<uint8_t>(WeightMode::PROJECTION)) "
       "sample operator: invalid weight mode"},
      {"pullback_operator_invalid_warp_envelope",
       case_pullback_operator_invalid_warp_envelope,
       "core/render/pullback/operators/warp.h",
       "(envelope <= static_cast<uint8_t>(WarpEnvelope::EDGE_FADE)) "
       "warp operator: invalid envelope"},
      {"pullback_operator_invalid_polar_mode",
       case_pullback_operator_invalid_polar_mode,
       "core/render/pullback/operators/warp.h",
       "(params.mode <= static_cast<uint8_t>(PolarMode::LOGARITHMIC)) "
       "warp.polar-chart: invalid polar mode"},
      {"pullback_operator_invalid_polar_harmonic",
       case_pullback_operator_invalid_polar_harmonic,
       "core/render/pullback/operators/warp.h",
       "(params.harmonic < Warp::MAX_POLAR_HARMONIC) "
       "warp.polar-chart: invalid harmonic"},
      {"pullback_operator_invalid_surface_integrator",
       case_pullback_operator_invalid_surface_integrator,
       "core/render/pullback/operators/sphere.h",
       "(params.integrator <= static_cast<uint8_t>(Surface::Integrator::MIDPOINT_2X)) sphere.displace.curl: invalid integrator"},
      {"pullback_operator_invalid_bonne_hemisphere",
       case_pullback_operator_invalid_bonne_hemisphere,
       "core/render/pullback/operators/project.h",
       "(params.hemisphere < std::size(BONNE_HEMISPHERE_IDS)) project.bonne: invalid hemisphere"},
      {"pullback_operator_invalid_airocean_layout",
       case_pullback_operator_invalid_airocean_layout,
       "core/render/pullback/operators/project.h",
       "(params.layout < std::size(AIROCEAN_LAYOUT_IDS)) project.airocean: invalid layout"},
      {"pullback_operator_invalid_peirce_layout",
       case_pullback_operator_invalid_peirce_layout,
       "core/render/pullback/operators/project.h",
       "(params.layout < std::size(PEIRCE_LAYOUT_IDS)) project.peirce: invalid layout"},
      {"pullback_operator_invalid_gnomonic_hemisphere",
       case_pullback_operator_invalid_gnomonic_hemisphere,
       "core/render/pullback/operators/project.h",
       "(params.hemisphere <= static_cast<uint8_t>(Projection::GnomonicHemisphere::BACK)) project.gnomonic: invalid hemisphere"},
      {"pullback_operator_invalid_curl_integrator",
       case_pullback_operator_invalid_curl_integrator,
       "core/render/pullback/operators/warp.h",
       "(params.integrator <= static_cast<uint8_t>(CurlIntegrator::MIDPOINT4)) warp.curl-flow: invalid integrator"},
      {"pullback_operator_invalid_noise_basis",
       case_pullback_operator_invalid_noise_basis,
       "core/render/pullback/operators/common.h",
       "(basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3)) "
       "pullback operator: invalid noise basis"},
      {"pullback_operator_invalid_tessellation_kind",
       case_pullback_operator_invalid_tessellation_kind,
       "core/render/pullback/operators/sample.h",
       "(params.kind <= "
       "static_cast<uint8_t>(Source::TessellationKind::HEXAGONAL)) "
       "sample.tessellation: invalid kind"},
      {"pullback_operator_invalid_kaleidoscope_symmetry",
       case_pullback_operator_invalid_kaleidoscope_symmetry,
       "core/render/pullback/operators/sphere.h",
       "(params.symmetry <= "
       "static_cast<uint8_t>(KaleidoscopeSymmetry::OCTAGONAL_PRISM)) "
       "sphere.lens.kaleidoscope: invalid symmetry"},
      {"pullback_operator_invalid_hue_mode",
       case_pullback_operator_invalid_hue_mode,
       "core/render/pullback/operators.h",
       "(params.hue_mode <= static_cast<uint8_t>(HueShiftMode::PATH_LENGTH)) "
       "colorize.generated-palette: invalid hue shift mode"},
      {"pullback_operator_invalid_palette_mode",
       case_pullback_operator_invalid_palette_mode,
       "core/render/pullback/operators.h",
       "(params.palette_mode < ctx.palettes.size()) "
       "colorize.generated-palette: invalid palette mode"},
      {"pullback_operator_invalid_palette_mapping",
       case_pullback_operator_invalid_palette_mapping,
       "core/render/pullback/operators.h",
       "(params.mapping_mode <= static_cast<uint8_t>("
       "Color::PaletteMapping::REVERSE)) "
       "colorize.generated-palette: invalid palette mapping"},
      {"pullback_operator_invalid_brightness_envelope",
       case_pullback_operator_invalid_brightness_envelope,
       "core/render/pullback/operators.h",
       "(params.envelope_mode <= static_cast<uint8_t>("
       "EnvelopeMode::DESCENDING)) "
       "colorize.generated-palette: invalid brightness envelope"},
  };
  n = static_cast<int>(sizeof(cases) / sizeof(cases[0]));
  return cases;
}

/**
 * @brief Dedicated always-trapping case proves the trap is observable.
 */
inline constexpr const char *SHAPE_PROBE_CASE = "__shape_probe__";
inline constexpr const char *DETERMINISM_PROBE_CASE =
    "__cross_process_determinism__";

/**
 * @brief Child entry point: runs exactly one named death case, then returns.
 * @param name Case selector; an unknown name (e.g. the "__spawn_check__"
 *             control) simply returns, so the child exits 0.
 * @details A case is expected to trap; returning means the child exits 0 and
 *          the parent flags it. The determinism selector instead emits cold/warm
 *          per-effect frame folds.
 */
inline void run_child_case(const char *name) {
#if defined(_WIN32)
  if (std::strcmp(name, "__literal_trap_exit__") == 0)
    std::exit(static_cast<int>(EXCEPTION_ILLEGAL_INSTRUCTION));
#endif
  if (std::strcmp(name, "__timeout_check__") == 0) {
    for (;;)
      std::this_thread::sleep_for(std::chrono::seconds(1));
  }
  if (std::strcmp(name, SHAPE_PROBE_CASE) == 0) {
    HS_CHECK(false, "death-harness trap-shape probe"); // always traps
    return;
  }
  if (std::strcmp(name, DETERMINISM_PROBE_CASE) == 0) {
    std::vector<Pixel> frame;
    uint64_t fold = 0;
    const auto CAPTURE = [&]<template <int, int> class E>() {
      effects_tests::render_capture<E, effects_tests::SMALL_W,
                                    effects_tests::SMALL_H>(frame, 8, &fold);
      const uint64_t COLD = fold;
      effects_tests::render_capture<E, effects_tests::SMALL_W,
                                    effects_tests::SMALL_H>(frame, 8, &fold);
      std::printf("capture %016llx %016llx\n",
                  static_cast<unsigned long long>(COLD),
                  static_cast<unsigned long long>(fold));
    };
#define HS_CAPTURE_COLD_WARM(effect) CAPTURE.operator()<effect>();
    HS_EFFECT_LIST(HS_CAPTURE_COLD_WARM)
#undef HS_CAPTURE_COLD_WARM
    return;
  }
  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i)
    if (case_enabled(cs[i]) && std::strcmp(cs[i].name, name) == 0) {
      cs[i].fn();
      return;
    }
}

/** @brief Bytes of child output kept for the guard-identity check. */
inline constexpr size_t CHILD_OUTPUT_CAP = 4096;

/**
 * @brief Accessor for the text captured from the most recent child spawn.
 * @return Pointer to the NUL-terminated capture buffer.
 */
inline char *child_output() {
  static char buf[CHILD_OUTPUT_CAP];
  return buf;
}

/**
 * @brief Path the spawned child's stdout/stderr is redirected to.
 * @return Stable NUL-terminated path, built once per process.
 * @details Named from the parent's pid so two test binaries running
 *          concurrently cannot overwrite each other's capture.
 */
inline const char *child_capture_path() {
  static char path[512];
  if (path[0] == '\0') {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#if defined(_WIN32)
    const char *dir = std::getenv("TEMP");
    if (!dir || dir[0] == '\0')
      dir = ".";
    std::snprintf(path, sizeof(path), "%s\\hs_death_%d.out", dir, _getpid());
#else
    const char *dir = std::getenv("TMPDIR");
    if (!dir || dir[0] == '\0')
      dir = "/tmp";
    std::snprintf(path, sizeof(path), "%s/hs_death_%d.out", dir,
                  static_cast<int>(getpid()));
#endif
#pragma clang diagnostic pop
  }
  return path;
}

/**
 * @brief Loads the tail of the capture file into child_output().
 * @details Keeps the last CHILD_OUTPUT_CAP-1 bytes, where a trapping child's
 *          breadcrumb lands.
 */
inline void load_child_output() {
  char *buf = child_output();
  buf[0] = '\0';
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *f = std::fopen(child_capture_path(), "rb");
#pragma clang diagnostic pop
  if (!f)
    return;
  if (std::fseek(f, 0, SEEK_END) == 0) {
    const long size = std::ftell(f);
    if (size > 0) {
      const long cap = static_cast<long>(CHILD_OUTPUT_CAP) - 1;
      const long keep = size < cap ? size : cap;
      if (std::fseek(f, size - keep, SEEK_SET) == 0)
        buf[std::fread(buf, 1, static_cast<size_t>(keep), f)] = '\0';
    }
  }
  std::fclose(f);
}

/**
 * @brief Tests whether captured child output carries one guard's breadcrumb.
 * @param out Text captured from the child.
 * @param file Repo-relative path of the guard source.
 * @param text Expected "(condition) message" tail of the breadcrumb line.
 * @return The matched source line, or zero when no breadcrumb matches.
 */
inline int breadcrumb_names_guard(const char *out, const char *file,
                                  const char *text) {
  std::string normalized(out);
  std::replace(normalized.begin(), normalized.end(), '\\', '/');
  const char *prefix = "HS_CHECK failed: ";
  const size_t file_len = std::strlen(file);
  const size_t text_len = std::strlen(text);
  for (const char *p = std::strstr(normalized.c_str(), prefix); p;
       p = std::strstr(p + 1, prefix)) {
    const char *path = p + std::strlen(prefix);
    const char *q = std::strchr(path, '\n');
    const char *end = q ? q : path + std::strlen(path);
    for (q = path; q < end; ++q) {
      if (*q != ':' || q + 1 == end || q[1] < '0' || q[1] > '9')
        continue;
      const size_t path_len = static_cast<size_t>(q - path);
      if (path_len < file_len ||
          std::strncmp(q - file_len, file, file_len) != 0 ||
          (path_len > file_len &&
           q[-static_cast<ptrdiff_t>(file_len) - 1] != '/'))
        break;
      ++q;
      int line = 0;
      while (q < end && *q >= '0' && *q <= '9') {
        line = line * 10 + (*q - '0');
        ++q;
      }
      if (end - q >= 2 && q[0] == ':' && q[1] == ' ' &&
          static_cast<size_t>(end - q - 2) >= text_len &&
          std::strncmp(q + 2, text, text_len) == 0)
        return line;
      break;
    }
  }
  return 0;
}

/**
 * @brief Sets HS_DEATH_CASE in this process's environment.
 * @param name Case selector to publish; the spawned child inherits it through
 *             the environment.
 */
inline void set_case_env(const char *name) {
#if defined(_WIN32)
  _putenv_s("HS_DEATH_CASE", name);
  _putenv_s("HS_DEATH_CHILD", name[0] ? "harness" : "");
#else
  setenv("HS_DEATH_CASE", name, 1);
  setenv("HS_DEATH_CHILD", name[0] ? "harness" : "", 1);
#endif
}

#if defined(_WIN32)
/** @brief Whether the latest child raised an unhandled illegal instruction. */
inline bool &child_unhandled_illegal_instruction() {
  static bool observed = false;
  return observed;
}
#endif

/**
 * @brief Spawns the test binary as a child running the given death case.
 * @param name Case selector passed to the child via HS_DEATH_CASE.
 * @param timeout_ms Maximum child runtime in milliseconds.
 * @return The child's raw status (debugged process on Windows, fork+execv wait status
 *         on POSIX). -1 on a spawn failure.
 * @details Child stdout/stderr are redirected to child_capture_path() and the
 *          tail is loaded into child_output(), so the caller can require the
 *          HS_CHECK breadcrumb of the guard the case is supposed to fire.
 */
inline int spawn_child(const char *name, unsigned timeout_ms = 10000) {
  set_case_env(name);
  const char *capture_path = child_capture_path();
  // A spawn that fails before the redirect must leave an EMPTY capture, not
  // the last child's breadcrumb.
  child_output()[0] = '\0';
  std::remove(capture_path);
#if defined(_WIN32)
  child_unhandled_illegal_instruction() = false;
  SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
  HANDLE capture = CreateFileA(capture_path, GENERIC_WRITE,
                               FILE_SHARE_READ | FILE_SHARE_WRITE, &security,
                               CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (capture == INVALID_HANDLE_VALUE)
    return -1;
  STARTUPINFOA startup{};
  startup.cb = sizeof(startup);
  startup.dwFlags = STARTF_USESTDHANDLES;
  startup.hStdOutput = capture;
  startup.hStdError = capture;
  startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
  PROCESS_INFORMATION process{};
  std::string command = std::string("\"") + self_exe() + "\"";
  const BOOL started =
      CreateProcessA(self_exe(), command.data(), nullptr, nullptr, TRUE,
                     DEBUG_ONLY_THIS_PROCESS | CREATE_NO_WINDOW, nullptr,
                     nullptr, &startup, &process);
  CloseHandle(capture);
  if (!started)
    return -1;
  int rc = -1;
  bool timed_out = false;
  bool stopping = false;
  bool initial_breakpoint = true;
  auto deadline =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout_ms);
  for (;;) {
    const auto now = std::chrono::steady_clock::now();
    if (now >= deadline) {
      if (stopping) {
        DebugActiveProcessStop(process.dwProcessId);
        WaitForSingleObject(process.hProcess, 5000);
        break;
      }
      timed_out = true;
      stopping = true;
      TerminateProcess(process.hProcess, 1);
      deadline = now + std::chrono::seconds(5);
    }
    const auto remaining =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            deadline - std::chrono::steady_clock::now())
            .count();
    DEBUG_EVENT event{};
    if (!WaitForDebugEvent(
            &event, static_cast<DWORD>(std::max<int64_t>(1, remaining)))) {
      if (GetLastError() == ERROR_SEM_TIMEOUT)
        continue;
      TerminateProcess(process.hProcess, 1);
      DebugActiveProcessStop(process.dwProcessId);
      WaitForSingleObject(process.hProcess, 5000);
      break;
    }
    DWORD continuation = DBG_CONTINUE;
    bool exited = false;
    switch (event.dwDebugEventCode) {
    case CREATE_PROCESS_DEBUG_EVENT:
      if (event.u.CreateProcessInfo.hFile)
        CloseHandle(event.u.CreateProcessInfo.hFile);
      break;
    case LOAD_DLL_DEBUG_EVENT:
      if (event.u.LoadDll.hFile)
        CloseHandle(event.u.LoadDll.hFile);
      break;
    case EXCEPTION_DEBUG_EVENT: {
      const auto &exception = event.u.Exception;
      const DWORD code = exception.ExceptionRecord.ExceptionCode;
      if (initial_breakpoint && exception.dwFirstChance &&
          code == EXCEPTION_BREAKPOINT) {
        initial_breakpoint = false;
      } else {
        continuation = DBG_EXCEPTION_NOT_HANDLED;
        if (!exception.dwFirstChance && code == EXCEPTION_ILLEGAL_INSTRUCTION)
          child_unhandled_illegal_instruction() = true;
      }
      break;
    }
    case EXIT_PROCESS_DEBUG_EVENT:
      if (!stopping)
        rc = static_cast<int>(event.u.ExitProcess.dwExitCode);
      exited = true;
      break;
    }
    if (!ContinueDebugEvent(event.dwProcessId, event.dwThreadId,
                            continuation)) {
      TerminateProcess(process.hProcess, 1);
      DebugActiveProcessStop(process.dwProcessId);
      WaitForSingleObject(process.hProcess, 5000);
      rc = -1;
      break;
    }
    if (exited) {
      if (WaitForSingleObject(process.hProcess, 5000) != WAIT_OBJECT_0)
        rc = -1;
      break;
    }
  }
  CloseHandle(process.hThread);
  CloseHandle(process.hProcess);
  load_child_output();
  if (timed_out)
    std::fprintf(stderr, "death child timed out after %u ms: %s\n", timeout_ms,
                 name);
  return rc;
#else
  // Shell-free spawn: fork + execv with stdout/stderr redirected to the
  // capture file; returns the raw wait status.
  std::fflush(stdout);
  std::fflush(stderr);
  const char *exe = self_exe();
  pid_t pid = fork();
  if (pid < 0)
    return -1;
  if (pid == 0) {
    int capture = open(capture_path, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (capture >= 0) {
      dup2(capture, 1);
      dup2(capture, 2);
      if (capture > 2)
        close(capture);
    }
    const char *argv[] = {exe, nullptr};
    execv(exe, const_cast<char *const *>(argv));
    _exit(127); // exec failed — never returns to the harness
  }
  int status = 0;
  const auto DEADLINE =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout_ms);
  for (;;) {
    const pid_t WAITED = waitpid(pid, &status, WNOHANG);
    if (WAITED == pid)
      break;
    if (WAITED < 0 && errno != EINTR)
      return -1;
    if (std::chrono::steady_clock::now() >= DEADLINE) {
      kill(pid, SIGKILL);
      while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {
      }
      load_child_output();
      std::fprintf(stderr, "death child timed out after %u ms: %s\n",
                   timeout_ms, name);
      return -1;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  load_child_output();
  return status;
#endif
}

#if defined(_WIN32)
/**
 * @brief The Windows EXCEPTION_ILLEGAL_INSTRUCTION process exit code.
 * @details Accepted only with a second-chance illegal-instruction debug event.
 */
inline constexpr int TRAP_STATUS = static_cast<int>(0xC000001D);
#endif

/**
 * @brief Tests whether a child died from an observable illegal instruction.
 * @param rc The raw spawn_child() return value to interpret.
 * @return True iff the child died by SIGILL or an unhandled Windows trap.
 */
inline bool child_trapped(int rc) {
#if defined(_WIN32)
  return rc == TRAP_STATUS && child_unhandled_illegal_instruction();
#else
  return rc != -1 && WIFSIGNALED(rc) && WTERMSIG(rc) == SIGILL;
#endif
}

/**
 * @brief Tests whether the child exited cleanly (exit code 0).
 * @param rc The raw spawn_child() return value to interpret.
 * @return True iff the child exited normally with status 0.
 */
inline bool child_exited_clean(int rc) {
#if defined(_WIN32)
  return rc == 0;
#else
  return rc != -1 && WIFEXITED(rc) && WEXITSTATUS(rc) == 0;
#endif
}

/**
 * @brief Reports that the death suite could not run.
 * @param why Human-readable reason the suite is unrunnable.
 * @param rc The associated child return code, for diagnostics.
 * @details An unrunnable death tier is a test failure on every host.
 */
inline void report_unrunnable(const char *why, int rc) {
  std::printf("  [FAIL] death tests: %s (rc=%d)\n", why, rc);
  HS_EXPECT_TRUE(false && "death suite must run");
}

/**
 * @brief Counts distinct fired guard lines in @p file.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param file Repo-relative source path.
 * @param lines Matched source line for each successfully trapped case, or zero.
 * @return Number of distinct covered source lines.
 */
inline int pinned_guards_in(const Case *cs, int n, const char *file,
                            const int *lines) {
  int pinned = 0;
  for (int i = 0; i < n; ++i) {
    if (!lines[i] || std::strcmp(cs[i].guard_file, file) != 0)
      continue;
    bool duplicate = false;
    for (int j = 0; j < i && !duplicate; ++j)
      duplicate =
          std::strcmp(cs[j].guard_file, file) == 0 && lines[j] == lines[i];
    if (!duplicate)
      ++pinned;
  }
  return pinned;
}

/** @brief One file's approved count of guard sites no case pins. */
struct GuardGapAllowance {
  const char *file; /**< Repo-relative source path from the census. */
  int gap;          /**< Unpinned HS_CHECK sites tolerated in that file. */
};

/**
 * @brief Per-file allowances for the sites the case table leaves unpinned.
 * @details Every census file's gap is gated against its row here, and a file
 *          with no row must be fully pinned. Each row is exact in both
 *          directions, so pinning or deleting a guard forces the row down.
 */
inline constexpr GuardGapAllowance GUARD_GAP_ALLOW[] = {
    {"core/animation/animation.h", 2},
    {"core/animation/carousel.h", 4},
    {"core/animation/motion.h", 6},
    {"core/animation/opleg.h", 41},
    {"core/animation/orientation.h", 12},
    {"core/animation/params.h", 9},
    {"core/animation/recipe_build.h", 21},
    {"core/animation/segue.h", 1},
    {"core/animation/sprites.h", 10},
    {"core/animation/timeline.h", 7},
    {"core/animation/transformer.h", 3},
    {"core/color/baked_palette.h", 9},
    {"core/color/color_space.h", 1},
    {"core/color/composition.h", 22},
    {"core/color/generative_palette.h", 4},
    {"core/color/palette_cycler.h", 8},
    {"core/containers/static_circular_buffer.h", 2},
    {"core/control/param_host.h", 15},
    {"core/control/preset_host.h", 2},
    {"core/memory/arena.h", 1},
    {"core/math/3dmath.h", 2},
    {"core/math/lenses.h", 1},
    {"core/math/pixel_mapping.h", 2},
    {"core/math/spherical_field.h", 2},
    {"core/mesh/conway.h", 31},
    {"core/mesh/conway_graph.h", 1},
    {"core/mesh/hankin.h", 8},
    {"core/mesh/mesh.h", 9},
    {"core/mesh/mesh_state.h", 2},
    {"core/mesh/recipe.h", 13},
    {"core/mesh/solid_builder.h", 5},
    {"core/mesh/solids.h", 1},
    {"core/render/canvas.h", 4},
    {"core/render/filter/pixel_feedback.h", 7},
    {"core/render/filter/screen_trails.h", 2},
    {"core/render/filter/world_trails.h", 2},
    {"core/render/plot/cull/spans.h", 1},
    {"core/render/plot/mesh.h", 4},
    {"core/render/plot/raster.h", 4},
    {"core/render/plot/shapes.h", 13},
    {"core/render/pullback/interpreter.h", 13},
    {"core/render/pullback/operators/model.h", 1},
    {"core/render/pullback/operators.h", 1},
    {"core/render/scan/mesh.h", 4},
    {"core/render/scan/raster.h", 1},
    {"core/render/scan/shader.h", 2},
    {"core/render/scan/shapes.h", 7},
    {"core/render/scan/volume.h", 4},
    {"core/render/sdf/common.h", 4},
    {"core/render/sdf/face.h", 1},
    {"core/render/sdf/face_geometry.h", 2},
    {"core/render/sdf/face_class_bake.h", 6},
    {"core/render/sdf/shapes.h", 8},
    {"core/render/sdf/volume.h", 2},
    {"core/render/shading.h", 1},
    {"core/spatial/kd_tree.h", 2},
    {"core/spatial/reaction_graph.h", 1},
    {"effects/Comets.h", 1},
    {"effects/DisplacementField.h", 1},
    {"effects/DreamBalls.h", 9},
    {"effects/Dynamo.h", 1},
    {"effects/Fishbowl.h", 1},
    {"effects/GnomonicStars.h", 1},
    {"effects/HankinSolids.h", 13},
    {"effects/HyperLattice.h", 2},
    {"effects/IslamicStars.h", 5},
    {"effects/MeshFeedback.h", 1},
    {"effects/MindSplatter.h", 5},
    {"effects/MobiusRings.h", 1},
    {"effects/ReactionDiffusionBase.h", 1},
    {"effects/RingShower.h", 1},
    {"effects/ShapeShifter.h", 2},
    {"hardware/dma_led.h", 4},
    {"hardware/pov_segmented.h", 9},
    {"hardware/pov_single.h", 10},
    {"targets/Phantasm/phantasm_target.h", 1},
    {"targets/Profile/Profile.ino", 4},
    // WASM-only bootstrap and reconstruction invariants; exercised by engine contracts.
    {"targets/wasm/engine_bindings.h", 8},
    {"targets/wasm/mesh_ops_bindings.h", 2},
    {"targets/wasm/workbench_bindings.h", 2},
    {"workbench/shader/chain_host.h", 1},
};

/**
 * @brief Looks up a file's approved unpinned-site count.
 * @param file Repo-relative source path from the census.
 * @return The approved gap, or 0 for a file with no allowance row.
 */
inline int allowed_guard_gap(const char *file) {
  int debug_gap = 0;
#ifdef NDEBUG
  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i) {
    if (!cs[i].debug_only || std::strcmp(cs[i].guard_file, file) != 0)
      continue;
    bool duplicate = false;
    for (int j = 0; j < i; ++j)
      duplicate |= cs[j].debug_only &&
                   std::strcmp(cs[j].guard_file, file) == 0 &&
                   std::strcmp(cs[j].guard_text, cs[i].guard_text) == 0;
    if (!duplicate)
      ++debug_gap;
  }
#endif
  for (const GuardGapAllowance &a : GUARD_GAP_ALLOW)
    if (std::strcmp(a.file, file) == 0)
      return a.gap + debug_gap;
  return debug_gap;
}

/**
 * @brief Prints what fraction of the engine's fail-fast surface is pinned.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param lines Matched source line for each successfully trapped case, or zero.
 * @details The denominator is the generated HS_CHECK census
 *          (death_guard_sites.h); the numerator is the distinct source lines
 *          observed in trapped cases. A case naming a file outside the census
 *          fails the module. Per-file gaps are gated against GUARD_GAP_ALLOW;
 *          the ratio is not gated.
 */
inline void report_guard_coverage(const Case *cs, int n, const int *lines) {
  int covered = 0;
  int off_census = 0;
  int unapproved_gaps = 0;
  int stale_allowances = 0;
  constexpr int GAPS = 5;
  const GuardSiteCount *worst[GAPS] = {};
  int worst_gap[GAPS] = {};
  for (const GuardSiteCount &f : GUARD_SITE_COUNTS) {
    int pinned = pinned_guards_in(cs, n, f.file, lines);
    HS_EXPECT_LE(pinned, f.sites);
    covered += pinned;
    int gap = f.sites - pinned;
    const int allowed = allowed_guard_gap(f.file);
    if (gap > allowed) {
      std::printf("  [FAIL] %s leaves %d HS_CHECK site(s) unpinned, %d "
                  "approved — add a death case, or write the gap down as "
                  "{\"%s\", %d} in GUARD_GAP_ALLOW\n",
                  f.file, gap, allowed, f.file, gap);
      ++unapproved_gaps;
    } else if (gap < allowed) {
      std::printf("  [FAIL] {\"%s\", %d} over-approves — the file's gap is %d; "
                  "lower the GUARD_GAP_ALLOW row to {\"%s\", %d}\n",
                  f.file, allowed, gap, f.file, gap);
      ++stale_allowances;
    }
    for (int slot = 0; slot < GAPS; ++slot) {
      if (gap <= worst_gap[slot])
        continue;
      for (int k = GAPS - 1; k > slot; --k) {
        worst[k] = worst[k - 1];
        worst_gap[k] = worst_gap[k - 1];
      }
      worst[slot] = &f;
      worst_gap[slot] = gap;
      break;
    }
  }
  for (int i = 0; i < n; ++i) {
    bool in_census = false;
    for (const GuardSiteCount &f : GUARD_SITE_COUNTS)
      if (std::strcmp(cs[i].guard_file, f.file) == 0) {
        in_census = true;
        break;
      }
    if (!in_census)
      ++off_census;
  }
  std::printf(
      "  guard coverage: %d/%d fail-fast sites pinned by a case (%d%%), "
      "%d case(s) outside the census\n",
      covered, GUARD_SITE_TOTAL,
      GUARD_SITE_TOTAL ? covered * 100 / GUARD_SITE_TOTAL : 0, off_census);
  std::printf("  widest gaps:");
  for (int slot = 0; slot < GAPS && worst[slot]; ++slot)
    std::printf(" %s %d/%d", worst[slot]->file,
                worst[slot]->sites - worst_gap[slot], worst[slot]->sites);
  std::printf("\n");
  for (const GuardGapAllowance &a : GUARD_GAP_ALLOW) {
    bool in_census = false;
    for (const GuardSiteCount &f : GUARD_SITE_COUNTS)
      if (std::strcmp(a.file, f.file) == 0) {
        in_census = true;
        break;
      }
    if (!in_census) {
      std::printf("  [FAIL] allowance \"%s\" names no guard-bearing file\n",
                  a.file);
      ++stale_allowances;
    }
  }
  HS_EXPECT_EQ(off_census, 0);
  HS_EXPECT_EQ(unapproved_gaps, 0);
  HS_EXPECT_EQ(stale_allowances, 0);
}

/**
 * @brief Parent entry point for the death module.
 * @return The module's failure count.
 * @details Spawn-checks the harness, then runs every case in a child and asserts
 *          each died by the exact trap status. Two fresh processes also compare
 *          per-effect frame folds across cold/warm runs and process boundaries;
 *          their complete records must fit in CHILD_OUTPUT_CAP.
 */
inline int run_death_tests() {
  hs_test::ModuleFixture fixture("death");

  if (!self_exe() || self_exe()[0] == '\0') {
    report_unrunnable("no argv[0] to re-exec", 0);
    return fixture.result();
  }

  HS_EXPECT_EQ(spawn_child("__timeout_check__", 100), -1);

  // Control: a child given an unknown case must exit cleanly.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    report_unrunnable("cannot re-exec self", control);
    set_case_env("");
    return fixture.result();
  }

  int n;
  const Case *cs = all_cases(n);

  // The sentinel proves the trap is observable independently of real cases.
  const int probe = spawn_child(SHAPE_PROBE_CASE);
  if (!child_trapped(probe)) {
    report_unrunnable("trap sentinel did not trap; trap status is unobservable",
                      probe);
    set_case_env("");
    return fixture.result();
  }

  HS_EXPECT_EQ(
      breadcrumb_names_guard(
          "HS_CHECK failed: C:\\tree\\core\\render\\sdf\\shapes.h:12: (false) probe\n",
          "core/render/sdf/shapes.h", "(false) probe"),
      12);
  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: core/other/shapes.h:12: (false) probe\n",
      "core/render/sdf/shapes.h", "(false) probe"));
  HS_EXPECT_FALSE(breadcrumb_names_guard(
      "HS_CHECK failed: notcore/render/sdf/shapes.h:12: (false) probe\n",
      "core/render/sdf/shapes.h", "(false) probe"));

  // The same sentinel proves the capture channel.
  if (!breadcrumb_names_guard(child_output(), "tests/test_death.h",
                              "(false) death-harness trap-shape probe")) {
    report_unrunnable("cannot capture child output; which guard fired is "
                      "unverifiable",
                      0);
    set_case_env("");
    return fixture.result();
  }

#if defined(_WIN32)
  const int literal_exit = spawn_child("__literal_trap_exit__");
  HS_EXPECT_EQ(literal_exit, TRAP_STATUS);
  HS_EXPECT_FALSE(child_trapped(literal_exit));
#endif

  const int determinism_a = spawn_child(DETERMINISM_PROBE_CASE);
  const std::string fold_a = child_output();
  const int determinism_b = spawn_child(DETERMINISM_PROBE_CASE);
  const std::string fold_b = child_output();
  HS_EXPECT_TRUE(child_exited_clean(determinism_a));
  HS_EXPECT_TRUE(child_exited_clean(determinism_b));
  constexpr size_t EFFECT_COUNT = static_cast<size_t>(HS_EFFECT_COUNT);
  constexpr size_t PREFIX_BYTES = sizeof("capture ") - 1;
  constexpr size_t PAYLOAD_BYTES = 33;
  size_t record_count = 0;
  size_t offset = 0;
  while ((offset = fold_a.find("capture ", offset)) != std::string::npos) {
    HS_CONTEXT("effect", record_count);
    offset += PREFIX_BYTES;
    if (offset + PAYLOAD_BYTES >= fold_a.size()) {
      HS_EXPECT_TRUE(false);
      break;
    }
    const std::string COLD = fold_a.substr(offset, 16);
    const std::string WARM = fold_a.substr(offset + 17, 16);
    HS_EXPECT_EQ(COLD.find_first_not_of("0123456789abcdef"), std::string::npos);
    HS_EXPECT_EQ(WARM.find_first_not_of("0123456789abcdef"), std::string::npos);
    HS_EXPECT_EQ(fold_a[offset + 16], ' ');
    HS_EXPECT_TRUE(fold_a[offset + PAYLOAD_BYTES] == '\r' ||
                   fold_a[offset + PAYLOAD_BYTES] == '\n');
    HS_EXPECT_EQ(COLD, WARM);
    ++record_count;
    offset += PAYLOAD_BYTES;
  }
  HS_EXPECT_EQ(record_count, EFFECT_COUNT);
  HS_EXPECT_EQ(fold_a, fold_b);

#if defined(_WIN32)
  DWORD handles_before = 0;
  HS_EXPECT_TRUE(GetProcessHandleCount(GetCurrentProcess(), &handles_before));
#endif

  std::vector<int> lines(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (!case_enabled(cs[i]))
      continue;
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc);
    // The child must die at this case's guard, not at another trap.
    int line = breadcrumb_names_guard(child_output(), cs[i].guard_file,
                                      cs[i].guard_text);
    bool at_guard = line != 0;
    lines[i] = trapped ? line : 0;
    HS_EXPECT_TRUE(trapped);
    HS_EXPECT_TRUE(at_guard);
    std::printf("  [%s] trap fires: %-26s (child rc=%d)\n",
                trapped && at_guard ? "ok" : "FAIL", cs[i].name, rc);
    if (trapped && !at_guard)
      std::printf("      expected %s: %s\n      child logged: %s\n",
                  cs[i].guard_file, cs[i].guard_text, child_output());
  }

#if defined(_WIN32)
  DWORD handles_after = 0;
  HS_EXPECT_TRUE(GetProcessHandleCount(GetCurrentProcess(), &handles_after));
  HS_EXPECT_EQ(handles_after, handles_before);
#endif

  std::remove(child_capture_path());
  set_case_env(""); // leave the env clean for anything that runs after us
  report_guard_coverage(cs, n, lines.data());
  return fixture.result();
}

} // namespace death_tests
} // namespace hs_test
