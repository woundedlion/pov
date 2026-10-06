/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file mesh_ops_bindings.h
 * @brief JS-facing mesh editor bridge: the MeshOpsWrapper class and its embind
 *        registration.
 *
 * Owns the lazily-allocated tooling arenas, the wipe generation counter and the
 * re-entrancy guard the Conway/Goldberg operators run under, plus the boundary
 * guards that keep a JS-driven operator chain out of the engine's fail-fast
 * traps.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/platform/platform.h"
#include "core/mesh/solids.h"
#include "targets/wasm/arena_metrics.h"
#include "targets/wasm/mesh_op_bounds.h"
#include "targets/wasm/wasm_predicates.h"
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>

// Arenas for the JS mesh-editor tools, malloc'd lazily on first MeshOps use
// (capacity 0 until then); the block lives for the module's lifetime.
inline constexpr size_t TOOLING_ARENA_BYTES = 8 * 1024 * 1024;
inline constexpr size_t TOOLING_SCRATCH_BYTES = 4 * 1024 * 1024;
static Arena tooling_arena(nullptr, 0);
// Scratch-using MeshOps entry points reset both arenas and hold a ToolingOpGuard.
// Scratch is valid only within one synchronous call.
static Arena tooling_scratch_a(nullptr, 0);
static Arena tooling_scratch_b(nullptr, 0);

// Largest element count any operator stage may reach: face/index counts narrow
// to uint16_t behind an always-on HS_CHECK, so a larger mesh is rejected at the
// JS boundary.
inline constexpr size_t MAX_MESH_CONNECTIVITY_ELEMENTS =
    MeshLimits::MAX_HALF_EDGES;
static_assert(MAX_MESH_CONNECTIVITY_ELEMENTS <=
                  TOOLING_SCRATCH_BYTES /
                      hs_wasm::TOOLING_BYTES_PER_MESH_ELEMENT,
              "a stage at the 16-bit ceiling must still fit a scratch arena");

static_assert(sizeof(math::Vector) + sizeof(uint8_t) + 2 * sizeof(uint16_t) <
                  hs_wasm::TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
              "finalized mesh element must fit its predicted arena bytes");

// Widest face a mesh can hold: per-face side counts are uint8_t and
// narrow_face_count traps past this.
inline constexpr size_t MAX_MESH_FACE_DEGREE = MeshLimits::MAX_FACE_DEGREE;

// Bumped on every clearToolingMemory(). Each wrapper records the generation it
// was built under and rejects via wrapper_live() if a wipe reclaimed its storage.
static uint32_t tooling_generation = 0;

// Module-global scratch permits one active MeshOps call; traps leave this
// latched because wasm unreachable does not unwind and the module is terminal.
static bool tooling_op_active = false;
struct ToolingOpGuard {
  ToolingOpGuard() {
    HS_CHECK(!tooling_op_active,
             "re-entrant MeshOps call aliases module-global tooling scratch");
    tooling_op_active = true;
  }
  ~ToolingOpGuard() { tooling_op_active = false; }
};

/**
 * @brief Allocates and binds the tooling arenas on first MeshOps use.
 * @return True once the three arenas are bound, false when the block
 *         allocation failed and none of them is usable.
 * @details A no-op once bound. Reading an unbound arena's metrics is safe and
 *          reports 0/0/0.
 */
static bool ensure_tooling_arenas() {
  if (tooling_arena.get_capacity() != 0)
    return true;
  const size_t total = TOOLING_ARENA_BYTES + 2 * TOOLING_SCRATCH_BYTES;
  uint8_t *block = static_cast<uint8_t *>(std::malloc(total));
  if (block == nullptr) {
    hs::log("WASM: tooling arena block allocation of %zu bytes failed — "
            "ignored",
            total);
    return false;
  }
  tooling_arena.rebind(block, TOOLING_ARENA_BYTES);
  tooling_scratch_a.rebind(block + TOOLING_ARENA_BYTES, TOOLING_SCRATCH_BYTES);
  tooling_scratch_b.rebind(block + TOOLING_ARENA_BYTES + TOOLING_SCRATCH_BYTES,
                           TOOLING_SCRATCH_BYTES);
  return true;
}

/**
 * @brief Builds a {usage, high_water_mark, lifetime_high_water_mark, capacity}
 *        report for the engine and tooling arenas.
 * @return JS object mapping each arena name to its {usage, high_water_mark,
 *         lifetime_high_water_mark, capacity} metrics, in bytes.
 */
static emscripten::val collect_arena_metrics() {
  emscripten::val metrics = collect_engine_arena_metrics();
  add_arena_metrics(metrics, "tooling_arena", tooling_arena);
  add_arena_metrics(metrics, "tooling_scratch_a", tooling_scratch_a);
  add_arena_metrics(metrics, "tooling_scratch_b", tooling_scratch_b);
  return metrics;
}

/**
 * @brief Why the most recent MeshOps call answered null.
 * @details Read back via MeshOps.getLastResult(). Exposed to JS as the
 *          Module.MeshOpResult embind enum; compare against its values, never by
 *          truthiness.
 */
enum class MeshOpResult {
  OK,                    /**< No rejection was recorded. */
  UNKNOWN_NAME,          /**< No registry entry carries that name. */
  CONNECTIVITY_OVERFLOW, /**< A stage would pass the 16-bit element ceiling. */
  FACE_DEGREE_OVERFLOW,  /**< A stage would emit a face past the 8-bit side
                              count. */
  ARENA_EXHAUSTED,       /**< The result would not fit tooling_arena's
                              remaining bytes. */
  NON_FINITE_ARG,        /**< An operator argument was NaN or infinite. */
  ANGLE_OUT_OF_DOMAIN,   /**< An angle argument sat outside its operator's
                              domain. */
  STALE_WRAPPER,         /**< The wrapper's storage was reclaimed by a
                              clearToolingMemory(). */
  ARENA_UNAVAILABLE,     /**< The tooling arena block could not be allocated;
                              arena-backed mesh operations are unavailable. */
};

// Outcome of the most recent MeshOps call that could answer null
// (getLastResult()).
static MeshOpResult last_mesh_op_result = MeshOpResult::OK;

// Whether that call saturated an argument into its operator's domain
// (getLastAdjusted()).
static bool last_mesh_op_adjusted = false;

/**
 * @brief JS-facing wrapper around a PolyMesh and the Conway/Goldberg operators.
 * @details Each wrapper's mesh is built into the tooling arena and records the generation it was
 *          built under, so a wipe via clearToolingMemory() is detected by
 *          wrapper_live().
 */
struct MeshOpsWrapper {
private:
  PolyMesh mesh; /**< The wrapped mesh, stored in the tooling arena. */
  /** Generation of the tooling arena this mesh was built into. */
  uint32_t generation = tooling_generation;

public:
  /**
   * @brief Constructs a wrapper taking ownership of an existing mesh.
   * @param m Mesh to move into this wrapper.
   */
  MeshOpsWrapper(PolyMesh &&m) : mesh(std::move(m)) {}

private:
  /**
   * @brief Opens a MeshOps entry point by clearing the previous call's outcome.
   * @details Result getters, getRegistry() and getArenaMetrics() preserve the
   *          previous outcome.
   */
  static void begin_mesh_op() {
    last_mesh_op_result = MeshOpResult::OK;
    last_mesh_op_adjusted = false;
  }

  /**
   * @brief Reports whether this wrapper outlived a clearToolingMemory() wipe.
   * @return true while its mesh still owns live arena storage; false — after
   *         logging and recording STALE_WRAPPER — once a wipe reclaimed it.
   * @details A stale wrapper's mesh aliases reclaimed arena storage. Rejects
   *          rather than traps.
   */
  bool wrapper_live() const {
    if (generation == tooling_generation)
      return true;
    hs::log("WASM: MeshOps wrapper used after clearToolingMemory() — ignored");
    last_mesh_op_result = MeshOpResult::STALE_WRAPPER;
    return false;
  }

  /**
   * @brief Runs the tooling boundary sequence one mesh must clear.
   * @param verts Mesh vertex count.
   * @param faces Mesh face count.
   * @param indices Mesh flat face-index count.
   * @param expansion Largest multiple of the mesh's biggest element count any
   *        stage reaches; 1 when the call only stores the mesh it was given.
   * @param context Call-site label prefixed to a rejection log.
   * @param bytes_per_element Tooling-arena bytes the call retains per output
   *        element; defaults to a whole finalized mesh.
   * @return true when the caller must answer null; last_mesh_op_result carries
   *         ARENA_UNAVAILABLE, CONNECTIVITY_OVERFLOW or ARENA_EXHAUSTED.
   * @details Binds the arenas first.
   */
  static bool
  tooling_bounds_reject(size_t verts, size_t faces, size_t indices,
                        size_t expansion, const char *context,
                        size_t bytes_per_element =
                            hs_wasm::TOOLING_ARENA_BYTES_PER_MESH_ELEMENT) {
    if (!ensure_tooling_arenas()) {
      last_mesh_op_result = MeshOpResult::ARENA_UNAVAILABLE;
      return true;
    }
    if ((expansion != 0 && verts > MeshLimits::MAX_VERTICES / expansion) ||
        hs_wasm::mesh_op_expansion_over_ceiling(
            verts, faces, indices, expansion, MAX_MESH_CONNECTIVITY_ELEMENTS)) {
      hs::log("WASM: %s: mesh of %zu verts / %zu faces / %zu indices expands "
              "%zux, past the %zu-element 16-bit connectivity range — ignored",
              context, verts, faces, indices, expansion,
              MAX_MESH_CONNECTIVITY_ELEMENTS);
      last_mesh_op_result = MeshOpResult::CONNECTIVITY_OVERFLOW;
      return true;
    }
    if (hs_wasm::mesh_op_output_over_arena(
            verts, faces, indices, expansion, bytes_per_element,
            tooling_arena.get_offset(), tooling_arena.get_capacity())) {
      hs::log("WASM: %s: result does not fit the tooling arena (%zu of %zu "
              "bytes used) — ignored; call clearToolingMemory() to reclaim it, "
              "which invalidates every live mesh",
              context, tooling_arena.get_offset(),
              tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return true;
    }
    return false;
  }

public:
  /**
   * @brief Resets all tooling arenas to empty and invalidates live wrappers.
   * @details Bumps the generation so any wrapper built before this wipe is
   *          rejected on next use. Does NOT shrink linear memory: only the arena
   *          bump-pointers reset. Not a post-trap recovery: the re-entrancy latch
   *          (tooling_op_active) stays set.
   */
  static void clearToolingMemory() {
    begin_mesh_op();
    tooling_arena.reset();
    tooling_arena.reset_high_water_mark();
    tooling_scratch_a.reset();
    tooling_scratch_a.reset_high_water_mark();
    tooling_scratch_b.reset();
    tooling_scratch_b.reset_high_water_mark();
    ++tooling_generation;
  }

  /**
   * @brief Builds a wrapper for the named base solid.
   * @param name Solid name to look up in the Solids registry.
   * @return Owning pointer to the new wrapper, or null for an unknown name, an
   *         unallocatable tooling arena, or a solid that would not fit what is
   *         left of tooling_arena; getLastResult() names which.
   * @details Generates into the scratch arenas and prices the finalized copy
   *          against tooling_arena's remaining bytes before committing it.
   */
  static std::unique_ptr<MeshOpsWrapper>
  fromSolidName(const std::string &name) {
    begin_mesh_op();
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: fromSolidName unknown solid '%s' — ignored", name.c_str());
      last_mesh_op_result = MeshOpResult::UNKNOWN_NAME;
      return nullptr;
    }
    ToolingOpGuard guard;
    if (!ensure_tooling_arenas()) {
      last_mesh_op_result = MeshOpResult::ARENA_UNAVAILABLE;
      return nullptr;
    }
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    const PolyMesh generated =
        entry->generate(tooling_scratch_a, tooling_scratch_b);
    const std::string context = "fromSolidName '" + name + "'";
    if (tooling_bounds_reject(generated.vertices.size(),
                              generated.get_face_counts_size(),
                              generated.get_faces_size(), 1, context.c_str()))
      return nullptr;
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(generated, tooling_arena));
  }

  /**
   * @brief Returns the mesh vertices as a JS Float32Array.
   * @return Float32Array of flattened [x,y,z] triples, copied out of the mesh,
   *         or null if a clearToolingMemory() reclaimed this wrapper's storage;
   *         getLastResult() then reports STALE_WRAPPER.
   * @details Copies directly from packed arena storage without allocating in
   *          WASM, so this readback cannot detach outstanding memory views.
   */
  emscripten::val getVertices() const {
    begin_mesh_op();
    if (!wrapper_live())
      return emscripten::val::null();
    static_assert(std::is_standard_layout_v<math::Vector> &&
                      offsetof(math::Vector, x) == 0 &&
                      offsetof(math::Vector, y) == sizeof(float) &&
                      offsetof(math::Vector, z) == 2 * sizeof(float) &&
                      sizeof(math::Vector) == 3 * sizeof(float) &&
                      alignof(math::Vector) == alignof(float),
                  "flat [x,y,z] view requires tightly packed vertices");
    return emscripten::val::global("Float32Array")
        .new_(emscripten::val(emscripten::typed_memory_view(
            mesh.vertices.size() * 3,
            reinterpret_cast<const float *>(mesh.vertices.data()))));
  }

  /**
   * @brief Returns the mesh faces as flat index + per-face side-count buffers.
   * @return JS object `{ indices: Uint16Array, counts: Uint8Array }`; JS
   *         unflattens the parallel arrays into per-face index lists. Both are
   *         copied out of WASM memory (same tooling-arena lifetime contract as
   *         getVertices()), so they are safe to hold across later mesh ops.
   *         Null if a clearToolingMemory() reclaimed this wrapper's storage;
   *         getLastResult() then reports STALE_WRAPPER.
   */
  emscripten::val getFaces() const {
    begin_mesh_op();
    if (!wrapper_live())
      return emscripten::val::null();
    emscripten::val out = emscripten::val::object();
    out.set("indices", emscripten::val::global("Uint16Array")
                           .new_(emscripten::val(emscripten::typed_memory_view(
                               mesh.get_faces_size(), mesh.get_faces_data()))));
    out.set("counts", emscripten::val::global("Uint8Array")
                          .new_(emscripten::val(emscripten::typed_memory_view(
                              mesh.get_face_counts_size(),
                              mesh.get_face_counts_data()))));
    return out;
  }

  /**
   * @brief Classifies faces by topology and returns the per-face codes.
   * @return JS Int32Array of one topology code per face, copied out of the
   *         mesh's now-populated topology buffer, or null when this wrapper's
   *         storage was reclaimed, the mesh is past
   *         MAX_MESH_CONNECTIVITY_ELEMENTS, the tooling arena could not be
   *         allocated, or its topology block would not fit what is left of
   *         tooling_arena; getLastResult() names which.
   * @details The mesh's topology buffer lives in tooling_arena until
   *          clearToolingMemory(). The returned Int32Array is a JS-owned copy,
   *          valid across later mesh operations and arena resets.
   */
  emscripten::val classifyFaces() {
    begin_mesh_op();
    if (!wrapper_live())
      return emscripten::val::null();
    // Only the topology block lands in tooling_arena: one uint16_t per face.
    if (tooling_bounds_reject(mesh.vertices.size(), mesh.get_face_counts_size(),
                              mesh.get_faces_size(), 1, "classifyFaces",
                              sizeof(uint16_t)))
      return emscripten::val::null();
    ToolingOpGuard guard;
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    MeshOps::classify_faces_by_topology(mesh, tooling_scratch_a,
                                        tooling_scratch_b, tooling_arena);
    return emscripten::val::global("Int32Array")
        .new_(emscripten::val(emscripten::typed_memory_view(
            mesh.topology.size(), mesh.topology.data())));
  }

  // --- Conway/Goldberg operators -------------------------------------------

private:
  /**
   * @brief Highest vertex valence in this wrapper's mesh.
   * @return Faces meeting at the most-incident vertex, or 0 for an empty mesh.
   * @details Counts incidences in tooling_scratch_a, released before return.
   */
  size_t max_vertex_valence() {
    const size_t verts = mesh.vertices.size();
    if (verts == 0)
      return 0;
    tooling_scratch_a.reset();
    uint32_t *incidence = tooling_scratch_a.allocate_n<uint32_t>(verts);
    const size_t valence = hs_wasm::mesh_max_vertex_valence(
        mesh.get_faces_data(), mesh.get_faces_size(), incidence, verts);
    tooling_scratch_a.reset();
    return valence;
  }

  /**
   * @brief Runs a mesh operator and wraps the result.
   * @tparam Op Callable of signature (const PolyMesh&, Arena&, Arena&) ->
   *         PolyMesh.
   * @param bounds Operator's growth factors (see MESHOP_LIST).
   * @param op Operator to run against this wrapper's mesh.
   * @return Owning pointer to a new wrapper holding the finalized result mesh, or
   *         null if this wrapper's storage was reclaimed, some stage of this
   *         operator would pass MAX_MESH_CONNECTIVITY_ELEMENTS or
   *         MAX_MESH_FACE_DEGREE, the tooling arena could not be allocated, or
   *         its output would not fit what is left of tooling_arena;
   *         getLastResult() names which.
   * @details Runs the op into the scratch arenas and finalizes the result into
   *          tooling_arena. Does not call begin_mesh_op(): the caller's argument
   *          clamp record must survive.
   */
  template <typename Op>
  std::unique_ptr<MeshOpsWrapper> apply(hs_wasm::MeshOpBounds bounds, Op &&op) {
    if (!wrapper_live())
      return nullptr;
    if (tooling_bounds_reject(mesh.vertices.size(), mesh.get_face_counts_size(),
                              mesh.get_faces_size(), bounds.elements,
                              "MeshOps operator"))
      return nullptr;
    ToolingOpGuard guard;
    const size_t face_degree = hs_wasm::mesh_max_face_degree(
        mesh.get_face_counts_data(), mesh.get_face_counts_size());
    const size_t valence = bounds.valence == 0 ? 0 : max_vertex_valence();
    if (hs_wasm::mesh_op_face_degree_overflows(
            face_degree, valence, bounds.face_degree, bounds.valence,
            MAX_MESH_FACE_DEGREE)) {
      hs::log("WASM: MeshOps input mesh (widest face %zu sides x%zu, highest "
              "valence %zu x%zu) would emit a face past the %zu-side limit — "
              "ignored",
              face_degree, bounds.face_degree, valence, bounds.valence,
              MAX_MESH_FACE_DEGREE);
      last_mesh_op_result = MeshOpResult::FACE_DEGREE_OVERFLOW;
      return nullptr;
    }
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        op(mesh, tooling_scratch_a, tooling_scratch_b), tooling_arena));
  }

  /**
   * @brief Validates that an operator argument is finite.
   * @param arg Operator argument crossing the untrusted JS boundary. Taken as a
   *        double so a count-valued argument is tested at its own width.
   * @param op Operator name, for the rejection log message.
   * @return true if arg is finite; false (after logging) otherwise.
   */
  bool finite_arg(double arg, const char *op) const {
    if (std::isfinite(arg))
      return true;
    hs::log("WASM: MeshOps::%s got a non-finite argument — ignored", op);
    last_mesh_op_result = MeshOpResult::NON_FINITE_ARG;
    return false;
  }

  /**
   * @brief Saturates a fraction argument into its operator's domain.
   * @param arg Fraction crossing the untrusted JS boundary.
   * @param out_of_domain Whether @p arg sits outside that domain.
   * @param clamped @p arg saturated into it.
   * @param op Operator name, for the adjustment log.
   * @param domain Domain notation for the adjustment log.
   * @return @p clamped.
   * @details Logs the adjustment and records it for getLastAdjusted().
   */
  static float note_clamped_arg(double arg, bool out_of_domain, float clamped,
                                const char *op, const char *domain) {
    last_mesh_op_adjusted = out_of_domain;
    if (out_of_domain)
      hs::log("WASM: MeshOps::%s clamped %g to %s", op, arg, domain);
    return clamped;
  }

public:
/**
 * @brief Defines a zero-argument Conway/Goldberg operator method.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details The generated method runs MeshOps::name(mesh) via apply() and returns
 *          a new wrapper holding the result.
 */
#define MESHOP_0(name, elements, degree, valence)                              \
  std::unique_ptr<MeshOpsWrapper> name() {                                     \
    begin_mesh_op();                                                           \
    return apply({elements, degree, valence},                                  \
                 [](const PolyMesh &m, Arena &a, Arena &b) {                   \
                   return MeshOps::name(m, a, b);                              \
                 });                                                           \
  }
/**
 * @brief Defines a one-number-argument operator whose argument is a [0,1]
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Rejects a non-finite arg, then clamps to [0,1]; the clamp is
 *          recorded for getLastAdjusted() and logged.
 */
#define MESHOP_1U(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(double arg) {                           \
    begin_mesh_op();                                                           \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    const float t =                                                            \
        note_clamped_arg(arg, hs_wasm::unit_fraction_out_of_range(arg),        \
                         hs_wasm::clamp_unit_fraction(arg), #name, "[0,1]");   \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

/**
 * @brief Defines a one-number-argument operator whose argument is a [0,1)
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Like MESHOP_1U, for operators that assert `t < 1.0f`.
 */
#define MESHOP_1H(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(double arg) {                           \
    begin_mesh_op();                                                           \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    const float t = note_clamped_arg(                                          \
        arg, hs_wasm::half_open_fraction_out_of_range(arg),                    \
        hs_wasm::clamp_half_open_fraction(arg), #name, "[0,1)");               \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

  MESHOP_LIST(MESHOP_0, MESHOP_1U, MESHOP_1H)

#undef MESHOP_0
#undef MESHOP_1U
#undef MESHOP_1H

  /** Ceiling on relax smoothing passes. */
  static constexpr int MAX_RELAX_ITERATIONS = 1000;

  /**
   * @brief Applies relax smoothing passes to the mesh.
   * @param iterations Number of smoothing passes, truncated toward zero;
   *        floored at 0 and clamped to MAX_RELAX_ITERATIONS.
   * @return Owning pointer to a new wrapper holding the relaxed mesh, or null if
   *         the input, arena, or allocation is invalid; getLastResult() names why.
   * @details Non-finite counts are rejected. The clamp is reported by
   *          getLastAdjusted() and logged.
   */
  std::unique_ptr<MeshOpsWrapper> relax(double iterations) {
    begin_mesh_op();
    if (!finite_arg(iterations, "relax"))
      return nullptr;
    const int clamped =
        hs_wasm::clamp_relax_iterations(iterations, MAX_RELAX_ITERATIONS);
    last_mesh_op_adjusted = static_cast<double>(clamped) != iterations;
    if (last_mesh_op_adjusted)
      hs::log("WASM: MeshOps::relax clamped %g iterations to %d", iterations,
              clamped);
    return apply(hs_wasm::RELAX_BOUNDS,
                 [clamped](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::relax(m, a, b, clamped);
                 });
  }

  /**
   * Inclusive upper bound of the Hankin contact-angle domain; past pi/2 an angle
   * mirrors one already in domain.
   */
  static constexpr float MAX_HANKIN_ANGLE = math::PI_F / 2.0f;

  /**
   * @brief Applies the Hankin interlace operator to the mesh.
   * @param radians Interlace angle in radians (the unit MeshOps::hankin
   *        expects), in the operator's [0, MAX_HANKIN_ANGLE] domain.
   * @return Owning pointer to a new wrapper holding the result, or null if the
   *         angle is non-finite or out of domain; getLastResult() names which.
   * @details An out-of-domain angle is rejected, not clamped.
   */
  std::unique_ptr<MeshOpsWrapper> hankin(double radians) {
    begin_mesh_op();
    if (!finite_arg(radians, "hankin"))
      return nullptr;
    if (hs_wasm::hankin_angle_out_of_range(radians, MAX_HANKIN_ANGLE)) {
      hs::log("WASM: MeshOps::hankin angle %g outside [0, %g] — ignored",
              radians, MAX_HANKIN_ANGLE);
      last_mesh_op_result = MeshOpResult::ANGLE_OUT_OF_DOMAIN;
      return nullptr;
    }
    return apply(hs_wasm::HANKIN_BOUNDS,
                 [radians](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::hankin(m, a, b, static_cast<float>(radians));
                 });
  }

  /**
   * @brief Applies the chiral snub operator with explicit inset and twist.
   * @param t Inset factor of each face toward its centroid, clamped to [0, 1)
   *          (its documented domain) at the JS boundary.
   * @param twist Per-face rotation about the face normal, in radians (0 = none);
   *          finite values saturate to the engine float range before narrowing.
   * @return Owning pointer to a new wrapper, or null on invalid input, stale wrapper, or allocation failure.
   * @details A clamped inset is recorded for getLastAdjusted() and logged.
   */
  std::unique_ptr<MeshOpsWrapper> snub(double t, double twist) {
    begin_mesh_op();
    if (!finite_arg(t, "snub") || !finite_arg(twist, "snub"))
      return nullptr;
    const float ct =
        note_clamped_arg(t, hs_wasm::half_open_fraction_out_of_range(t),
                         hs_wasm::clamp_half_open_fraction(t), "snub", "[0,1)");
    const bool inset_adjusted = last_mesh_op_adjusted;
    const float clamped_twist = hs_wasm::clamp_finite_float(twist);
    const float rotation =
        note_clamped_arg(twist,
                         twist != static_cast<double>(clamped_twist) &&
                             (twist > std::numeric_limits<float>::max() ||
                              twist < -std::numeric_limits<float>::max()),
                         clamped_twist, "snub twist", "finite float range");
    last_mesh_op_adjusted = last_mesh_op_adjusted || inset_adjusted;
    return apply(hs_wasm::SNUB_BOUNDS,
                 [ct, rotation](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::snub(m, a, b, ct, rotation);
                 });
  }
  /**
   * @brief Lists all available solids for the editor's solid picker.
   * @return JS array of {name, category} objects, one per registered solid.
   */
  static emscripten::val getRegistry() {
    emscripten::val registry = emscripten::val::array();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      const auto &entry = Solids::get_entry(i);
      emscripten::val item = emscripten::val::object();
      item.set("name", emscripten::val(entry.name));
      item.set("category", entry.category == Solids::Category::Simple
                               ? "Simple"
                               : "Complex");
      registry.set(i, item);
    }
    return registry;
  }

  /**
   * @brief Maps an authored Solids::Op to the editor's lowercase op string.
   * @param op Authored operator from a recipe step.
   * @return The op name matching the solids editor's vocabulary.
   */
  static const char *op_name(Solids::Op op) {
    switch (op) {
    case Solids::Op::TRUNCATE:
      return "truncate";
    case Solids::Op::EXPAND:
      return "expand";
    case Solids::Op::SNUB:
      return "snub";
    case Solids::Op::CHAMFER:
      return "chamfer";
    case Solids::Op::HANKIN:
      return "hankin";
    case Solids::Op::RELAX:
      return "relax";
    case Solids::Op::KIS:
      return "kis";
    case Solids::Op::DUAL:
      return "dual";
    case Solids::Op::AMBO:
      return "ambo";
    case Solids::Op::BEVEL:
      return "bevel";
    case Solids::Op::GYRO:
      return "gyro";
    case Solids::Op::META:
      return "meta";
    case Solids::Op::NEEDLE:
      return "needle";
    case Solids::Op::ZIP:
      return "zip";
    }
    HS_CHECK(false, "unhandled Solids::Op %d in op_name", static_cast<int>(op));
  }

  /**
   * @brief Returns a solid's authored recipe chain for the editor.
   * @param name Registry name to look up.
   * @return JS object {seed: string, ops: [{op: string, param, twist}]}, or
   *         null for an unknown name (getLastResult() then reports
   *         UNKNOWN_NAME) or for a known entry without a recipe (OK).
   * @details Pure table read: no arenas, no wrapper. Params cross in
   *          engine-native units (radians for hankin, raw t, relax iteration
   *          count). A recipe-less known entry returns null without logging.
   */
  static emscripten::val getRecipe(const std::string &name) {
    begin_mesh_op();
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: getRecipe unknown solid '%s' — ignored", name.c_str());
      last_mesh_op_result = MeshOpResult::UNKNOWN_NAME;
      return emscripten::val::null();
    }
    if (!entry->recipe)
      return emscripten::val::null();
    const Solids::Recipe &recipe = *entry->recipe;
    // recipe.seed indexes simple_registry, not get_entry's combined index.
    emscripten::val out = emscripten::val::object();
    out.set("seed", emscripten::val(Solids::simple_registry[recipe.seed].name));
    emscripten::val ops = emscripten::val::array();
    for (size_t i = 0; i < recipe.count; ++i) {
      const Solids::OpStep &step = recipe.steps[i];
      emscripten::val item = emscripten::val::object();
      item.set("op", emscripten::val(op_name(step.op)));
      item.set("param", step.param);
      item.set("twist", step.twist);
      ops.set(i, item);
    }
    out.set("ops", ops);
    return out;
  }

#ifdef HS_WASM_DEV_BINDINGS
  /**
   * @brief Measures the maximum vertex/face/index counts across all solids.
   * @return JS object with {max_v, v_name, max_f, f_name, max_i, i_name} giving
   *         the largest counts and the solids that produce them, or null if the
   *         tooling arenas could not be allocated.
   * @details Dev-only, compiled with the HS_WASM_DEV_BINDINGS CMake option.
   *          Measures in the scratch arenas only, never tooling_arena.
   */
  static emscripten::val getMaxBounds() {
    begin_mesh_op();
    int max_v = 0;
    int max_f = 0;
    int max_i = 0;
    const char *mv_name = "";
    const char *mf_name = "";
    const char *mi_name = "";

    ToolingOpGuard guard;
    if (!ensure_tooling_arenas()) {
      last_mesh_op_result = MeshOpResult::ARENA_UNAVAILABLE;
      return emscripten::val::null();
    }
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      tooling_scratch_a.reset();
      tooling_scratch_b.reset();
      PolyMesh temp =
          Solids::get_entry(i).generate(tooling_scratch_a, tooling_scratch_b);

      int v = static_cast<int>(temp.vertices.size());
      int f = static_cast<int>(temp.get_face_counts_size());
      int idxs = static_cast<int>(temp.get_faces_size());

      if (v > max_v) {
        max_v = v;
        mv_name = Solids::get_entry(i).name;
      }
      if (f > max_f) {
        max_f = f;
        mf_name = Solids::get_entry(i).name;
      }
      if (idxs > max_i) {
        max_i = idxs;
        mi_name = Solids::get_entry(i).name;
      }
    }

    tooling_scratch_a.reset();
    tooling_scratch_b.reset();

    emscripten::val stats = emscripten::val::object();
    stats.set("max_v", max_v);
    stats.set("v_name", emscripten::val(mv_name));
    stats.set("max_f", max_f);
    stats.set("f_name", emscripten::val(mf_name));
    stats.set("max_i", max_i);
    stats.set("i_name", emscripten::val(mi_name));
    return stats;
  }
#endif

  /**
   * @brief Reports the engine and tooling arena metrics.
   * @return JS object of {usage, high_water_mark, lifetime_high_water_mark,
   *         capacity} metrics per arena, in bytes.
   */
  static emscripten::val getArenaMetrics() { return collect_arena_metrics(); }

  /**
   * @brief Reports why the most recent checked mesh operation answered null.
   * @return OK when no rejection was recorded, otherwise the rejection reason.
   * @details Read it immediately after the null; the next checked call
   *          overwrites it.
   */
  static MeshOpResult getLastResult() { return last_mesh_op_result; }

  /**
   * @brief Reports whether the most recent checked mesh operation saturated an argument
   *        into its operator's domain.
   * @return true when the mesh that call produced was rendered from a value
   *         other than the one passed in. Meaningful only when the call
   *         returned a mesh.
   * @details A clamped argument leaves getLastResult() OK; a caller exporting
   *          the argument it passed must check this, since the raw value would
   *          trip the engine's HS_CHECK. The next checked operation overwrites
   *          it.
   */
  static bool getLastAdjusted() { return last_mesh_op_adjusted; }
};

/** @brief Registers the mesh editor bridge's enum and class with Embind. */
static void bind_mesh_ops() {
  emscripten::enum_<MeshOpResult>("MeshOpResult")
      .value("OK", MeshOpResult::OK)
      .value("UNKNOWN_NAME", MeshOpResult::UNKNOWN_NAME)
      .value("CONNECTIVITY_OVERFLOW", MeshOpResult::CONNECTIVITY_OVERFLOW)
      .value("FACE_DEGREE_OVERFLOW", MeshOpResult::FACE_DEGREE_OVERFLOW)
      .value("ARENA_EXHAUSTED", MeshOpResult::ARENA_EXHAUSTED)
      .value("NON_FINITE_ARG", MeshOpResult::NON_FINITE_ARG)
      .value("ANGLE_OUT_OF_DOMAIN", MeshOpResult::ANGLE_OUT_OF_DOMAIN)
      .value("STALE_WRAPPER", MeshOpResult::STALE_WRAPPER)
      .value("ARENA_UNAVAILABLE", MeshOpResult::ARENA_UNAVAILABLE);

  // No public .constructor<>(): construction goes through fromSolidName.
  emscripten::class_<MeshOpsWrapper>("MeshOps")
      .class_function("clearToolingMemory", &MeshOpsWrapper::clearToolingMemory)
      .class_function("getLastResult", &MeshOpsWrapper::getLastResult)
      .class_function("getLastAdjusted", &MeshOpsWrapper::getLastAdjusted)
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName)
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
      .class_function("getRecipe", &MeshOpsWrapper::getRecipe)
#ifdef HS_WASM_DEV_BINDINGS
      .class_function("getMaxBounds", &MeshOpsWrapper::getMaxBounds)
#endif
      .class_function("getArenaMetrics", &MeshOpsWrapper::getArenaMetrics)
      .function("getVertices", &MeshOpsWrapper::getVertices)
      .function("getFaces", &MeshOpsWrapper::getFaces)
      .function("classifyFaces", &MeshOpsWrapper::classifyFaces)
// Binds both MESHOP_LIST metadata entries and MESHOP_IRREGULAR_LIST bare names.
#define MESHOP_BIND(name, ...) .function(#name, &MeshOpsWrapper::name)
          MESHOP_LIST(MESHOP_BIND, MESHOP_BIND, MESHOP_BIND)
              MESHOP_IRREGULAR_LIST(MESHOP_BIND);
#undef MESHOP_BIND
}
