/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file mesh_ops_bindings.h
 * @brief JS-facing mesh editor bridge: the MeshOpsWrapper class and its embind
 *        registration.
 *
 * Owns the lazily-allocated tooling arenas, the wipe generation counter and the
 * re-entrancy guard the Conway/Goldberg operators run under, plus the boundary
 * guards that keep a JS-driven operator chain out of the engine's fail-fast
 * traps. Included only by targets/wasm/wasm.cpp.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/engine/platform.h"
#include "core/mesh/solids.h"
#include "targets/wasm/arena_metrics.h"
#include "targets/wasm/mesh_op_bounds.h" // pure, host-tested operator roster
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <cstdlib> // std::malloc for the lazily-allocated tooling arenas
#include <cmath>   // std::isfinite — validate MeshOps args at the JS boundary
#include <memory>
#include <string>
#include <vector>

using namespace emscripten;

// Arenas for the JS mesh-editor tools (8 MB build + two 4 MB scratch), used
// only by MeshOpsWrapper. malloc'd lazily on first MeshOps use (start at
// capacity 0) so engine/worker instances that never touch MeshOps don't reserve
// 16 MB; the block lives for the module's lifetime (reset via clearToolingMemory).
static constexpr size_t TOOLING_ARENA_BYTES = 8 * 1024 * 1024;
static constexpr size_t TOOLING_SCRATCH_BYTES = 4 * 1024 * 1024;
static Arena tooling_arena(nullptr, 0);
// Transient single-op scratch, shared module-globally. Every MeshOps entry
// point reset()s both at its head; valid only within one synchronous call. A
// ToolingOpGuard at each entry traps any re-entrant or interleaved use (e.g. a
// future worker or async refactor) before it can alias this scratch.
static Arena tooling_scratch_a(nullptr, 0);
static Arena tooling_scratch_b(nullptr, 0);

// Worst-case scratch bytes one operator allocates per element of its largest
// stage, counting vertices, faces and flat face indices alike. The heaviest
// operator (meta = kis of dual of ambo) keeps an intermediate mesh, its output and
// its index scratch live in one arena, at 12 bytes per vertex and 6 per index;
// this bound covers that with margin.
static constexpr size_t TOOLING_BYTES_PER_MESH_ELEMENT = 64;

// Largest element count any operator stage may reach. Half-edge construction and
// the topology classifier both narrow face/index counts to uint16_t behind an
// always-on HS_CHECK, so a mesh past this must be rejected at the JS boundary
// rather than allowed to reach one. The scratch arenas hold exactly this many
// elements, so the same ceiling also keeps Arena::allocate's trap out of reach.
static constexpr size_t MAX_MESH_CONNECTIVITY_ELEMENTS = UINT16_MAX;
static_assert(MAX_MESH_CONNECTIVITY_ELEMENTS <=
                  TOOLING_SCRATCH_BYTES / TOOLING_BYTES_PER_MESH_ELEMENT,
              "a stage at the 16-bit ceiling must still fit a scratch arena");

// Tooling-arena bytes a finalized mesh retains per element: one Vector of
// vertex, one uint8_t of side count, one uint16_t of index and the uint16_t
// topology code classifyFaces() later binds into the same arena, with slack for
// the four blocks' alignment padding.
static constexpr size_t TOOLING_ARENA_BYTES_PER_MESH_ELEMENT = 20;
static_assert(sizeof(Vector) + sizeof(uint8_t) + 2 * sizeof(uint16_t) <
                  TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
              "finalized mesh element must fit its predicted arena bytes");

// Widest face a mesh can hold: PolyMesh stores per-face side counts as uint8_t
// and narrow_face_count traps past this, so an operator that would emit a wider
// face must be rejected at the JS boundary.
static constexpr size_t MAX_MESH_FACE_DEGREE = UINT8_MAX;

// Bumped on every clearToolingMemory(). Each wrapper records the generation it
// was built under and rejects via wrapper_live() if a wipe reclaimed its storage.
static uint32_t tooling_generation = 0;

// Set for the duration of one MeshOps entry point. The tooling scratch arenas
// are module-global and reset() at each entry's head, so a second op entered
// before the first returns would alias the first's scratch and corrupt its
// geometry. Synchronous single-threaded calls never overlap today; the guard
// makes that contract enforced rather than implicit, so a future worker/async
// refactor traps here instead of silently aliasing.
//
// A trap inside an op compiles to wasm `unreachable`, which unwinds nothing, so
// ~ToolingOpGuard() never runs and this stays set. clearToolingMemory() clears
// it, so a JS caller that catches the RuntimeError can recover through there.
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
 * @details A no-op once bound, so it is cheap to call at the head of every
 *          MeshOps entry point. Reading an unbound arena's metrics
 *          (collect_arena_metrics) is safe and reports 0/0/0, so engine
 *          instances that never call MeshOps never trigger this allocation.
 */
static void ensure_tooling_arenas() {
  if (tooling_arena.get_capacity() != 0)
    return;
  const size_t total = TOOLING_ARENA_BYTES + 2 * TOOLING_SCRATCH_BYTES;
  uint8_t *block = static_cast<uint8_t *>(std::malloc(total));
  HS_CHECK(block != nullptr,
           "tooling arena block allocation of %zu bytes failed", total);
  tooling_arena.rebind(block, TOOLING_ARENA_BYTES);
  tooling_scratch_a.rebind(block + TOOLING_ARENA_BYTES, TOOLING_SCRATCH_BYTES);
  tooling_scratch_b.rebind(block + TOOLING_ARENA_BYTES + TOOLING_SCRATCH_BYTES,
                           TOOLING_SCRATCH_BYTES);
}

/**
 * @brief Builds a {usage, high_water_mark, capacity} report for the three engine
 *        arenas and the three tooling arenas.
 * @return JS object mapping each arena name to its {usage, high_water_mark,
 *         capacity} metrics, in bytes.
 * @details Read on demand through MeshOps.getArenaMetrics(). The tooling scratch
 *          arenas are the regions TOOLING_BYTES_PER_MESH_ELEMENT is sized
 *          against, and an operator that overruns one takes the module down, so
 *          they are reported alongside the rest rather than left unobservable.
 */
static val collect_arena_metrics() {
  val metrics = collect_engine_arena_metrics();
  add_arena_metrics(metrics, "tooling_arena", tooling_arena);
  add_arena_metrics(metrics, "tooling_scratch_a", tooling_scratch_a);
  add_arena_metrics(metrics, "tooling_scratch_b", tooling_scratch_b);
  return metrics;
}

/**
 * @brief Why the most recent MeshOps call answered null.
 * @details A bare null collapses reasons that demand opposite caller actions —
 *          shrinking the op chain versus calling clearToolingMemory() — so the
 *          reason is recorded here and read back via MeshOps.getLastResult().
 *          Exposed to JS as the Module.MeshOpResult embind enum; compare
 *          against its values, never by truthiness.
 */
enum class MeshOpResult {
  OK,                    /**< The call produced a result. */
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
};

// Outcome of the most recent MeshOps call that could answer null
// (getLastResult()).
static MeshOpResult last_mesh_op_result = MeshOpResult::OK;

/**
 * @brief JS-facing wrapper around a PolyMesh and the Conway/Goldberg operators.
 * @details Named to avoid collision with the MeshOps namespace. Each wrapper's
 *          mesh is built into the tooling arena and records the generation it was
 *          built under, so a wipe via clearToolingMemory() is detected by
 *          wrapper_live().
 */
struct MeshOpsWrapper {
private:
  PolyMesh mesh; /**< The wrapped mesh, stored in the tooling arena. */
  /**
   * Generation of the tooling arena this mesh was built into; compared against
   * the live counter on every use (see wrapper_live()).
   */
  uint32_t generation = tooling_generation;

public:
  /**
   * @brief Constructs a wrapper taking ownership of an existing mesh.
   * @param m Mesh to move into this wrapper.
   */
  MeshOpsWrapper(PolyMesh &&m) : mesh(std::move(m)) {}

private:
  /**
   * @brief Reports whether this wrapper outlived a clearToolingMemory() wipe.
   * @return true while its mesh still owns live arena storage; false — after
   *         logging and recording STALE_WRAPPER — once a wipe reclaimed it.
   * @details A stale wrapper's mesh aliases reclaimed arena storage, which
   *          release builds would otherwise read back as silently wrong
   *          geometry. Called at every entry point that touches `mesh`, which
   *          rejects rather than traps: a JS caller holding a wrapper across an
   *          interleaved clearToolingMemory() is an ordering slip the page can
   *          recover from, not an engine invariant violation.
   */
  bool wrapper_live() const {
    if (generation == tooling_generation)
      return true;
    hs::log("WASM: MeshOps wrapper used after clearToolingMemory() — ignored");
    last_mesh_op_result = MeshOpResult::STALE_WRAPPER;
    return false;
  }

public:
  /**
   * @brief Resets all tooling arenas to empty and invalidates live wrappers.
   * @details Reclaims the storage behind every live wrapper and bumps the
   *          generation so any wrapper built before this wipe is rejected on
   *          next use (wrapper_live). JS-callable.
   *
   *          Despite the name, this does NOT shrink the module's linear memory:
   *          the 16 MB tooling block is retained for the module's lifetime and
   *          only its arena bump-pointers are reset. A JS caller will not see
   *          memory usage drop; this frees the arenas for reuse, not the OS.
   *
   *          Also clears the re-entrancy flag, which a trap inside an op leaves
   *          set (see tooling_op_active). No op is in flight once JS regains
   *          control, so this is the recovery entry point after a trap.
   */
  static void clearToolingMemory() {
    tooling_op_active = false;
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
   * @return Owning pointer to the new wrapper, or null for an unknown name or a
   *         solid that would not fit what is left of tooling_arena;
   *         getLastResult() names which.
   * @details Rejects an unknown name at the untrusted JS boundary rather than
   *          tripping get_by_name()'s fail-fast HS_CHECK and aborting the module.
   *          Generates into the scratch arenas and prices the finalized copy
   *          against tooling_arena's remaining bytes before committing it, since
   *          that arena accumulates one finalized mesh per live wrapper until
   *          clearToolingMemory() and Arena::allocate traps when it runs out.
   */
  static std::unique_ptr<MeshOpsWrapper>
  fromSolidName(const std::string &name) {
    last_mesh_op_result = MeshOpResult::OK;
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: fromSolidName unknown solid '%s' — ignored", name.c_str());
      last_mesh_op_result = MeshOpResult::UNKNOWN_NAME;
      return nullptr;
    }
    ToolingOpGuard guard;
    ensure_tooling_arenas();
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    const PolyMesh generated =
        entry->generate(tooling_scratch_a, tooling_scratch_b);
    if (hs_wasm::mesh_op_output_over_arena(
            generated.vertices.size(), generated.get_face_counts_size(),
            generated.get_faces_size(), 1, TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
            tooling_arena.get_offset(), tooling_arena.get_capacity())) {
      hs::log("WASM: fromSolidName '%s' does not fit the tooling arena (%zu of "
              "%zu bytes used) — ignored; call clearToolingMemory() to reclaim "
              "it, which invalidates every live mesh",
              name.c_str(), tooling_arena.get_offset(),
              tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return nullptr;
    }
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(generated, tooling_arena));
  }

  /**
   * @brief Returns the mesh vertices as a JS Float32Array.
   * @return Float32Array of flattened [x,y,z] triples, copied out of the mesh,
   *         or null if a clearToolingMemory() reclaimed this wrapper's storage;
   *         getLastResult() then reports STALE_WRAPPER.
   */
  val getVertices() const {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
    std::vector<float> data;
    data.reserve(mesh.vertices.size() * 3);
    for (const auto &v : mesh.vertices) {
      data.push_back(v.x);
      data.push_back(v.y);
      data.push_back(v.z);
    }
    return val::global("Float32Array")
        .new_(val(typed_memory_view(data.size(), data.data())));
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
  val getFaces() const {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
    size_t total = 0;
    for (size_t i = 0; i < mesh.get_face_counts_size(); ++i)
      total += mesh.get_face_counts_data()[i];
    HS_CHECK(total == mesh.get_faces_size(),
             "getFaces: face_counts sum disagrees with the flat index count");
    val out = val::object();
    out.set("indices", val::global("Uint16Array")
                           .new_(val(typed_memory_view(
                               mesh.get_faces_size(), mesh.get_faces_data()))));
    out.set("counts",
            val::global("Uint8Array")
                .new_(val(typed_memory_view(mesh.get_face_counts_size(),
                                            mesh.get_face_counts_data()))));
    return out;
  }

  /**
   * @brief Classifies faces by topology and returns the per-face codes.
   * @return JS Int32Array of one topology code per face, copied out of the
   *         mesh's now-populated topology buffer, or null when this wrapper's
   *         storage was reclaimed, the mesh is past
   *         MAX_MESH_CONNECTIVITY_ELEMENTS, or its topology block would not fit
   *         what is left of tooling_arena; getLastResult() names which. Null
   *         rather than an empty array so a caller can tell "no classification"
   *         from "no faces" with a plain truthiness test.
   * @details Same tooling-arena lifetime contract as getVertices(): the
   *          `topology` buffer lives in tooling_arena and is invalidated by the
   *          next mesh op / arena reset. The `.new_(Int32Array)(view)` form
   *          *copies* the typed_memory_view into a fresh JS array (it does not
   *          alias WASM memory), so the result is safe to hold across later
   *          calls — but if this is ever changed to return the view directly
   *          (as getPixels/bakeLut do), it MUST then be read before the next
   *          allocation, per that memory-view contract.
   */
  val classifyFaces() {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
    if (hs_wasm::tooling_mesh_over_ceiling(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), MAX_MESH_CONNECTIVITY_ELEMENTS)) {
      hs::log("WASM: classifyFaces mesh of %zu verts / %zu faces / %zu indices "
              "is past the %zu-element 16-bit connectivity range — ignored",
              mesh.vertices.size(), mesh.get_face_counts_size(),
              mesh.get_faces_size(), MAX_MESH_CONNECTIVITY_ELEMENTS);
      last_mesh_op_result = MeshOpResult::CONNECTIVITY_OVERFLOW;
      return val::null();
    }
    ensure_tooling_arenas();
    if (hs_wasm::mesh_op_output_over_arena(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), 1, TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
            tooling_arena.get_offset(), tooling_arena.get_capacity())) {
      hs::log("WASM: classifyFaces topology block does not fit the tooling "
              "arena (%zu of %zu bytes used) — ignored; call "
              "clearToolingMemory() to reclaim it, which invalidates every "
              "live mesh",
              tooling_arena.get_offset(), tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return val::null();
    }
    ToolingOpGuard guard;
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    MeshOps::classify_faces_by_topology(mesh, tooling_scratch_a,
                                        tooling_scratch_b, tooling_arena);
    // Copies (see the contract note above); does not alias WASM memory.
    return val::global("Int32Array")
        .new_(
            val(typed_memory_view(mesh.topology.size(), mesh.topology.data())));
  }

  // --- Conway/Goldberg operators -------------------------------------------

private:
  /**
   * @brief Highest vertex valence in this wrapper's mesh.
   * @return Faces meeting at the most-incident vertex, or 0 for an empty mesh.
   * @details Counts incidences in tooling_scratch_a (4 bytes per vertex, released
   *          before the operator runs), so the whole measurement is one pass over
   *          the flat index list on top of clearing the counters.
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
   *         MAX_MESH_FACE_DEGREE, or its output would not fit what is left of
   *         tooling_arena; getLastResult() names which.
   * @details Captures the shared operator boilerplate: reset both tooling scratch
   *          arenas, run the op into a fresh PolyMesh, finalize it into
   *          tooling_arena, and hand back a new wrapper. An input the operator
   *          would grow past the engine's 16-bit connectivity range is rejected
   *          here, at the untrusted JS boundary, rather than reaching
   *          build_half_edge_mesh's fail-fast trap partway through a composition;
   *          likewise an input holding a face or a vertex valence the operator
   *          would widen past narrow_face_count's 8-bit side count, and a result
   *          that would overrun tooling_arena, which accumulates one finalized
   *          mesh per chained wrapper until clearToolingMemory().
   */
  template <typename Op>
  std::unique_ptr<MeshOpsWrapper> apply(hs_wasm::MeshOpBounds bounds, Op &&op) {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return nullptr;
    if (hs_wasm::mesh_op_expansion_over_ceiling(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), bounds.elements,
            MAX_MESH_CONNECTIVITY_ELEMENTS)) {
      hs::log("WASM: MeshOps input mesh of %zu verts / %zu faces / %zu indices "
              "expands %zux, past the %zu-element 16-bit connectivity range — "
              "ignored",
              mesh.vertices.size(), mesh.get_face_counts_size(),
              mesh.get_faces_size(), bounds.elements,
              MAX_MESH_CONNECTIVITY_ELEMENTS);
      last_mesh_op_result = MeshOpResult::CONNECTIVITY_OVERFLOW;
      return nullptr;
    }
    ensure_tooling_arenas();
    // Before max_vertex_valence(), the first tooling-scratch use on this path.
    ToolingOpGuard guard;
    if (hs_wasm::mesh_op_output_over_arena(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), bounds.elements,
            TOOLING_ARENA_BYTES_PER_MESH_ELEMENT, tooling_arena.get_offset(),
            tooling_arena.get_capacity())) {
      hs::log("WASM: MeshOps result does not fit the tooling arena (%zu of %zu "
              "bytes used) — ignored; call clearToolingMemory() to reclaim it, "
              "which invalidates every live mesh",
              tooling_arena.get_offset(), tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return nullptr;
    }
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
   * @param arg Operator argument crossing the untrusted JS boundary.
   * @param op Operator name, for the rejection log message.
   * @return true if arg is finite; false (after logging) otherwise.
   * @details A non-finite fraction would flow straight into the geometry math and
   *          silently corrupt the mesh, so it is rejected (log + null) rather than
   *          producing NaN geometry.
   */
  bool finite_arg(float arg, const char *op) const {
    if (std::isfinite(arg))
      return true;
    hs::log("WASM: MeshOps::%s got a non-finite argument — ignored", op);
    last_mesh_op_result = MeshOpResult::NON_FINITE_ARG;
    return false;
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
    return apply({elements, degree, valence},                                  \
                 [](const PolyMesh &m, Arena &a, Arena &b) {                   \
                   return MeshOps::name(m, a, b);                              \
                 });                                                           \
  }
/**
 * @brief Defines a one-float-argument operator whose argument is a [0,1]
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Rejects a non-finite arg (finite_arg), then clamps the fraction to
 *          [0,1] (logging when it changed) so a direct/API caller passing a
 *          finite out-of-range value stays within the operator's documented
 *          domain and cannot trip truncate's or bevel's always-on HS_CHECK and
 *          abort the whole module.
 */
#define MESHOP_1U(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(float arg) {                            \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    if (hs_wasm::unit_fraction_out_of_range(arg))                              \
      hs::log("WASM: MeshOps::%s clamped %g to [0,1]", #name, arg);            \
    float t = hs_wasm::clamp_unit_fraction(arg);                               \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

/**
 * @brief Defines a one-float-argument operator whose argument is a [0,1)
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Like MESHOP_1U, but these operators assert `t < 1.0f`, so 1 is
 *          outside the domain and clamping to it would land on the trap instead
 *          of avoiding it.
 */
#define MESHOP_1H(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(float arg) {                            \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    if (hs_wasm::half_open_fraction_out_of_range(arg))                         \
      hs::log("WASM: MeshOps::%s clamped %g to [0,1)", #name, arg);            \
    float t = hs_wasm::clamp_half_open_fraction(arg);                          \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

  // The roster and its growth factors live in targets/wasm/mesh_op_bounds.h,
  // where a host test measures the real operators against them.
  MESHOP_LIST(MESHOP_0, MESHOP_1U, MESHOP_1H)

#undef MESHOP_0
#undef MESHOP_1U
#undef MESHOP_1H

  /**
   * Engine-enforced ceiling on relax smoothing passes: an independent
   * defense-in-depth limit for direct/API callers that bypass the editor's
   * 500-pass slider cap.
   */
  static constexpr int MAX_RELAX_ITERATIONS = 1000;

  /**
   * @brief Applies relax smoothing passes to the mesh.
   * @param iterations Number of smoothing passes; floored at 0 and clamped to
   *        MAX_RELAX_ITERATIONS.
   * @return Owning pointer to a new wrapper holding the relaxed mesh.
   * @details Explicit (not a MESHOP_* macro) because its int iteration count
   *          crosses the JS boundary unbounded: relax(1e9) would freeze the main
   *          thread for billions of passes, so the count is clamped rather than
   *          trusted.
   */
  std::unique_ptr<MeshOpsWrapper> relax(int iterations) {
    int clamped =
        hs_wasm::clamp_relax_iterations(iterations, MAX_RELAX_ITERATIONS);
    if (clamped != iterations)
      hs::log("WASM: MeshOps::relax clamped %d iterations to %d", iterations,
              clamped);
    iterations = clamped;
    return apply(hs_wasm::RELAX_BOUNDS,
                 [iterations](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::relax(m, a, b, iterations);
                 });
  }

  /**
   * Inclusive upper bound of the Hankin contact-angle domain: at pi/2 the
   * contact rays leave the edge perpendicular to it, and past that they tilt
   * back into the neighbouring face, mirroring an angle already in domain.
   */
  static constexpr float MAX_HANKIN_ANGLE = PI_F / 2.0f;

  /**
   * @brief Applies the Hankin interlace operator to the mesh.
   * @param radians Interlace angle in radians (the unit MeshOps::hankin
   *        expects), in the operator's [0, MAX_HANKIN_ANGLE] domain.
   * @return Owning pointer to a new wrapper holding the result, or null if the
   *         angle is non-finite or out of domain; getLastResult() names which.
   * @details Explicit (not a MESHOP_* macro) so the radians unit contract the JS
   *          caller relies on is carried here. The angle is rejected rather than
   *          clamped: MeshOps::hankin aliases an out-of-domain angle onto an
   *          in-domain pattern, so clamping would hand back geometry the caller
   *          did not ask for.
   */
  std::unique_ptr<MeshOpsWrapper> hankin(float radians) {
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
                   return MeshOps::hankin(m, a, b, radians);
                 });
  }

  /**
   * @brief Applies the chiral snub operator with explicit inset and twist.
   * @param t Inset factor of each face toward its centroid, clamped to [0, 1)
   *          (its documented domain) at the JS boundary.
   * @param twist Per-face rotation about the face normal, in radians (0 = none);
   *          unbounded, so only finiteness is checked.
   * @return Owning pointer to a new wrapper, or null if either arg is non-finite.
   * @details Explicit (not a MESHOP_* macro) because MeshOps::snub takes TWO
   *          float controls, which neither the zero-arg nor the one-float
   *          generator can express. Binding it via MESHOP_0 hardcodes the
   *          (0.5, 0.0) defaults and leaves both controls unreachable from JS;
   *          this 2-arg form exposes them to the solids editor.
   */
  std::unique_ptr<MeshOpsWrapper> snub(float t, float twist) {
    if (!finite_arg(t, "snub") || !finite_arg(twist, "snub"))
      return nullptr;
    if (hs_wasm::half_open_fraction_out_of_range(t))
      hs::log("WASM: MeshOps::snub clamped t=%g to [0,1)", t);
    float ct = hs_wasm::clamp_half_open_fraction(t);
    return apply(hs_wasm::SNUB_BOUNDS,
                 [ct, twist](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::snub(m, a, b, ct, twist);
                 });
  }
  /**
   * @brief Lists all available solids for the editor's solid picker.
   * @return JS array of {name, category} objects, one per registered solid.
   */
  static val getRegistry() {
    val registry = val::array();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      const auto &entry = Solids::get_entry(i);
      val item = val::object();
      item.set("name", val(entry.name));
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
   *         null for an unknown name or an entry without a recipe.
   * @details Pure table read: no arenas, no wrapper, no clearToolingMemory()
   *          pairing. Params cross in engine-native units (radians for hankin,
   *          raw t, relax iteration count), matching the MeshOps op bindings.
   *          An unknown name is rejected at the untrusted JS boundary
   *          (fromSolidName's precedent) rather than trapping; a recipe-less
   *          known entry is the normal not-morphable case and returns null
   *          without logging.
   */
  static val getRecipe(const std::string &name) {
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: getRecipe unknown solid '%s' — ignored", name.c_str());
      return val::null();
    }
    if (!entry->recipe)
      return val::null();
    const Solids::Recipe &recipe = *entry->recipe;
    // recipe.seed indexes simple_registry specifically (get_entry spans all
    // three registries and would misresolve it).
    val out = val::object();
    out.set("seed", val(Solids::simple_registry[recipe.seed].name));
    val ops = val::array();
    for (size_t i = 0; i < recipe.count; ++i) {
      const Solids::OpStep &step = recipe.steps[i];
      val item = val::object();
      item.set("op", val(op_name(step.op)));
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
   *         the largest counts and the solids that produce them.
   * @details Dev-only roster measurement for sizing MAX_VERTS-style constants;
   *          no UI consumer. Off by default; enable the HS_WASM_DEV_BINDINGS
   *          CMake option to compile + re-export it
   *          (`cmake --preset wasm-release -DHS_WASM_DEV_BINDINGS=ON`; see
   *          CMakeLists.txt). Measures each solid in the scratch arenas only —
   *          never tooling_arena, which backs live wrappers the JS side holds.
   */
  static val getMaxBounds() {
    int max_v = 0;
    int max_f = 0;
    int max_i = 0;
    const char *mv_name = "";
    const char *mf_name = "";
    const char *mi_name = "";

    ToolingOpGuard guard;
    ensure_tooling_arenas();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      // Measure in the scratch arenas only — never tooling_arena, which backs
      // live wrappers the JS side holds.
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

    val stats = val::object();
    stats.set("max_v", max_v);
    stats.set("v_name", val(mv_name));
    stats.set("max_f", max_f);
    stats.set("f_name", val(mf_name));
    stats.set("max_i", max_i);
    stats.set("i_name", val(mi_name));
    return stats;
  }
#endif

  /**
   * @brief Reports the engine arena metrics for the mesh tooling HUD.
   * @return JS object of {usage, high_water_mark, capacity} metrics per arena,
   *         in bytes.
   */
  static val getArenaMetrics() { return collect_arena_metrics(); }

  /**
   * @brief Reports why the most recent MeshOps call answered null.
   * @return OK when that call produced a result, otherwise its rejection reason.
   * @details Covers fromSolidName, getVertices, getFaces, classifyFaces and the
   *          operator methods. Read it immediately after the null; the next such
   *          call overwrites it.
   */
  static MeshOpResult getLastResult() { return last_mesh_op_result; }
};

/** @brief Registers the mesh editor bridge's enum and class with Embind. */
static void bind_mesh_ops() {
  enum_<MeshOpResult>("MeshOpResult")
      .value("OK", MeshOpResult::OK)
      .value("UNKNOWN_NAME", MeshOpResult::UNKNOWN_NAME)
      .value("CONNECTIVITY_OVERFLOW", MeshOpResult::CONNECTIVITY_OVERFLOW)
      .value("FACE_DEGREE_OVERFLOW", MeshOpResult::FACE_DEGREE_OVERFLOW)
      .value("ARENA_EXHAUSTED", MeshOpResult::ARENA_EXHAUSTED)
      .value("NON_FINITE_ARG", MeshOpResult::NON_FINITE_ARG)
      .value("ANGLE_OUT_OF_DOMAIN", MeshOpResult::ANGLE_OUT_OF_DOMAIN)
      .value("STALE_WRAPPER", MeshOpResult::STALE_WRAPPER);

  // No public .constructor<>(): all construction goes through fromSolidName so
  // JS cannot wrap an empty mesh past the operator boundary's wrapper_live().
  class_<MeshOpsWrapper>("MeshOps")
      .class_function("clearToolingMemory", &MeshOpsWrapper::clearToolingMemory)
      .class_function("getLastResult", &MeshOpsWrapper::getLastResult)
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
  // Bound from the same MESHOP_LIST that generates the wrapper methods, plus
  // MESHOP_IRREGULAR_LIST for the hand-written ops.
// Variadic so it can take both MESHOP_LIST's (name, expansion) entries and
// MESHOP_IRREGULAR_LIST's bare names.
#define MESHOP_BIND(name, ...) .function(#name, &MeshOpsWrapper::name)
          MESHOP_LIST(MESHOP_BIND, MESHOP_BIND, MESHOP_BIND)
              MESHOP_IRREGULAR_LIST(MESHOP_BIND);
#undef MESHOP_BIND
}
