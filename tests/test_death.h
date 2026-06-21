/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Death tests — the suite's coverage of the fail-fast philosophy the project
 * markets, spanning the memory/arena, math core, mesh-topology, registry,
 * container, spatial (KDTree), animation, canvas, scan, and plot seams (including non-finite-input
 * faults: NaN fed to normalize/slerp/make_rotation/make_basis must trap, not
 * propagate into geometry). An HS_CHECK violation is a deliberate __builtin_trap() that
 * aborts the whole process; the in-process HS_EXPECT_* harness cannot catch it.
 * So each trap is exercised in a CHILD process: the test binary re-exec's
 * itself with HS_DEATH_CASE=<name> (handled in main() before any module runs),
 * runs exactly one trap-triggering case, and the parent asserts the child died
 * by the *specific* trap status — clang lowers __builtin_trap() to an illegal
 * instruction (x86 ud2), so the child dies by SIGILL (POSIX) /
 * STATUS_ILLEGAL_INSTRUCTION (Windows). Asserting that exact status, not merely
 * "nonzero", means an unrelated child crash (a different signal, an ordinary
 * nonzero exit) can no longer masquerade as a passing death test.
 *
 * Cross-platform by design: the child is selected through an inherited env var
 * and spawned with std::system(), so no fork() (absent on Windows) is needed.
 * A control "spawn check" runs first; if the harness cannot re-exec itself
 * (unknown argv[0], a sandbox that blocks process creation), the death tests
 * are SKIPPED with a notice — EXCEPT under CI (the CI env var is set), where a
 * suite that cannot run is a hard FAILURE rather than a silent green skip, so a
 * regression in the fail-fast layer cannot slip through unobserved.
 */
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "tests/test_harness.h"

#include "core/3dmath.h"
#include "core/animation.h"
#include "core/canvas.h"
#include "core/geometry.h"
#include "core/filter.h"
#include "core/memory.h"
#include "core/mesh.h"
#include "core/plot.h"
#include "core/scan.h"
#include "core/solids.h"
#include "core/spatial.h"
#include "core/static_circular_buffer.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#endif

#if defined(_WIN32)
/**
 * @brief Local declaration of the Win32 SetErrorMode (no <windows.h>).
 * @param uMode Bitmask of error-mode flags; 0x0001|0x0002 =
 *              SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX.
 * @return The previous error-mode bitmask.
 * @details Declared locally to suppress the WER "stopped working" box for the
 *          children we intentionally crash. kernel32 is linked by default for a
 *          Windows-Clang console build.
 */
extern "C" __declspec(dllimport) unsigned int __stdcall
SetErrorMode(unsigned int uMode);
#endif

namespace hs_test {
namespace death {

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

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

/**
 * @brief Death case: arena over-allocation must trap.
 * @details Memory surface — requests more than the arena's capacity so
 *          allocate() fires HS_CHECK.
 */
inline void case_arena_oom() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(1024)); // > capacity -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x"); // keep the call live
}

/**
 * @brief Death case: rewinding the arena past its capacity must trap.
 * @details Memory surface — set_offset is the one seam that can break the
 *          offset <= capacity invariant the no-wrap bounds math in allocate()
 *          depends on, so an out-of-range rewind traps at the source.
 */
inline void case_arena_set_offset_overflow() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.set_offset(opaque<size_t>(sizeof(buf) + 1)); // > capacity -> HS_CHECK
}

/**
 * @brief Death case: ArenaVector fixed-capacity push_back overflow must trap.
 * @details Arena-container surface — a push_back past capacity fires HS_CHECK.
 */
inline void case_arena_vector_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2);
  v.push_back(1);
  v.push_back(2);
  v.push_back(opaque(3)); // exceeds capacity -> HS_CHECK
}

/**
 * @brief Death case: normalizing a degenerate (zero-length) vector must trap.
 * @details Math-core surface — length below epsilon fires the normalize guard.
 */
inline void case_normalize_zero() {
  Vector z{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  Vector n = z.normalized(); // length < eps -> HS_CHECK
  if (n.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: normalizing a NaN vector must trap.
 * @details Math-core surface — a NaN coordinate poisons the length to NaN, and
 *          `NaN >= epsilon` is false, so the normalize guard fires. The suite's
 *          NaN/Inf fault case: proves a non-finite producer is trapped at the
 *          math seam rather than silently propagating NaN into geometry.
 */
inline void case_normalize_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  Vector n = bad.normalized(); // length is NaN -> HS_CHECK fails
  if (n.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an out-of-range solids index must trap.
 * @details Lookup/registry surface — get_entry past NUM_ENTRIES fires HS_CHECK.
 */
inline void case_solids_index_oob() {
  const auto &e = Solids::get_entry(opaque<size_t>(Solids::NUM_ENTRIES));
  if (e.name == nullptr)
    std::printf("x");
}

/**
 * @brief Death case: looking up an unknown solid name must trap.
 * @details Registry-by-name surface — an unknown name has no valid fallback.
 */
inline void case_solids_unknown_name() {
  configure_arenas_default();
  PolyMesh m = Solids::get_by_name(persistent_arena, scratch_arena_a,
                                   scratch_arena_b, "definitely_not_a_solid");
  if (m.vertices.size() == 0x7fff)
    std::printf("x");
}

/**
 * @brief Death case: a StaticCircularBuffer index past the live count must trap.
 * @details Container surface — index >= count fires HS_CHECK.
 */
inline void case_circular_buffer_oob() {
  StaticCircularBuffer<int, 4> cb;
  cb.push_back(10);
  cb.push_back(20);
  int v = cb[opaque<size_t>(5)]; // index >= count -> HS_CHECK
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: front() on an empty StaticCircularBuffer must trap.
 * @details Container surface — the never-taken opaque(false) push keeps the
 *          optimizer from proving the buffer empty and folding the trap at
 *          compile time; is_empty() fires HS_CHECK.
 */
inline void case_circular_buffer_front_empty() {
  StaticCircularBuffer<int, 4> cb;
  if (opaque(false))
    cb.push_back(1);
  int v = cb.front(); // is_empty() -> HS_CHECK("front() on empty ...")
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: ArenaVector::append_bulk past its fixed capacity must trap.
 * @details Memory surface — a distinct seam from element-at-a-time push_back;
 *          the bulk memcpy path has its own size_+count guard.
 */
inline void case_arena_vector_append_bulk_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2); // exact capacity 2
  int src[4] = {1, 2, 3, 4};
  v.append_bulk(src, opaque<size_t>(4)); // 0 + 4 > 2 -> HS_CHECK
  if (v.size() == 0x7fff)
    std::printf("x");
}

/**
 * @brief Death case: requesting more KDTree neighbors than MAX_K must trap.
 * @details Spatial surface — k beyond the MAX_K-sized result/heap buffers makes
 *          nearest() trap rather than silently capping the result and masking
 *          the caller's sizing mistake.
 */
inline void case_spatial_knn_over_max() {
  static uint8_t buf[512];
  Arena a(buf, sizeof(buf));
  Vector pts[2] = {Vector(1.0f, 0.0f, 0.0f), Vector(0.0f, 1.0f, 0.0f)};
  KDTree tree(a, std::span<const Vector>(pts, 2));
  // Tree is non-empty and k > 0, so the k <= MAX_K guard is reached.
  auto r = tree.nearest(Vector(1.0f, 0.0f, 0.0f),
                        opaque<size_t>(KDTree::MAX_K + 1)); // k > MAX_K -> HS_CHECK
  if (r.size() == static_cast<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: an over-subscribed arena partition must trap.
 * @details Config surface — each request alone fits but the sum exceeds
 *          GLOBAL_ARENA_SIZE, so configure_arenas fires HS_CHECK.
 */
inline void case_arena_oversubscribed() {
  // Each request alone fits, but the sum exceeds GLOBAL_ARENA_SIZE -> HS_CHECK.
  configure_arenas(opaque(GLOBAL_ARENA_SIZE), opaque<size_t>(1024),
                   opaque<size_t>(1024));
}

/**
 * @brief Death case: relocating a retained (pinned) add_get() handle must trap.
 * @details Animation surface — step()'s compaction routes every relocation
 *          through TimelineEvent::move_into, which traps when the event was
 *          handed out via add_get(pin=true), converting the dangling-handle
 *          hazard into a fail-fast crash instead of silent corruption.
 */
inline void case_timeline_handled_relocation() {
  TimelineEvent src;
  src.handled = opaque(true); // as if handed out by add_get(pin=true)
  TimelineEvent dst;
  src.move_into(dst); // HS_CHECK(!handled) -> trap
}

/**
 * @brief Death case: a pinned (handled) animation that COMPLETES must trap.
 * @details Animation surface — the symmetric companion to
 *          case_timeline_handled_relocation, which guards the relocation path
 *          (move_into). A pinned-but-finite animation that finishes as the
 *          *last* event needs no relocation, so move_into never runs; step()'s
 *          completion branch would otherwise e.destroy() it and dangle the
 *          caller's retained pointer silently. The pin contract is
 *          pinned => infinite, so a pinned animation that naturally completes is
 *          misuse; the completion branch's HS_CHECK traps it. (A deliberate
 *          cancel() is exempt — see is_canceled() — so this case completes
 *          naturally rather than canceling.)
 */
inline void case_timeline_handled_completion() {
  // A do-nothing 8x8 effect whose fresh buffer_free() is true, so the Canvas
  // ctor does not spin. The animation never renders through it.
  struct DeathEffect : public Effect {
    DeathEffect() : Effect(8, 8) {}
    void draw_frame() override {}
    bool show_bg() const override { return false; }
  };
  static DeathEffect fx;
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  // pin=true marks the event handled; a 1-frame Transition is finite and the
  // sole event, so step() routes it through the completion/destroy branch.
  tl.add_get(0, Animation::Transition(v, 1.0f, 1, ease_mid), /*pin=*/true);
  tl.step(canvas); // t=1: done() && !repeats() && !canceled, keep=false -> trap
}

/**
 * @brief Death case: a second simultaneously-live Timeline must trap.
 * @details Animation surface — every Timeline shares the single global event
 *          array, so a second live instance would silently stomp the first's
 *          events; the construction guard traps instead. The real app holds
 *          exactly one (the old effect is destroyed before the next is built).
 */
inline void case_timeline_double_construct() {
  Timeline a;
  Timeline b; // second live ctor -> HS_CHECK(!global_timeline_live) -> trap
  if (global_timeline_num_events == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: narrowing an index past the int16 topology range must trap.
 * @details Mesh-topology surface — both conway.h and hankin.h route every output
 *          vertex/face-index narrowing through this shared MeshOps guard, so a
 *          future MAX_VERTS bump traps at the bench instead of silently wrapping
 *          an index and corrupting topology.
 */
inline void case_mesh_narrow_index() {
  size_t over = static_cast<size_t>(INT16_MAX) + 1;
  uint16_t i = MeshOps::narrow_index(opaque(over)); // > INT16_MAX -> HS_CHECK
  if (i == 0xBEEF)
    std::printf("x");
}

/**
 * @brief Death case: a NaN endpoint fed to slerp must trap.
 * @details Math-core surface — the NaN poisons interpolation through both
 *          branches into the final strict normalized(), which traps rather than
 *          emitting a NaN direction into geometry. Proves the non-finite input
 *          is caught at the slerp seam, not just at bare normalize().
 */
inline void case_slerp_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  Vector dst{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  Vector v = slerp(bad, dst, opaque(0.5f)); // NaN -> normalized() -> HS_CHECK
  if (v.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a NaN source must trap.
 * @details Math-core surface — the d-based parallel/antiparallel guards are
 *          NaN-false, so it falls through to cross(from,to).normalized(), which
 *          traps on the NaN-poisoned axis.
 */
inline void case_make_rotation_vectors_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector from{nan, opaque(0.0f), opaque(0.0f)};
  Vector to{opaque(0.0f), opaque(0.0f), opaque(1.0f)};
  Quaternion q = make_rotation(from, to); // NaN axis -> normalized() -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(axis, theta) with a NaN angle must trap.
 * @details Math-core surface — cos/sin of a NaN poison the quaternion, and its
 *          normalized() traps on the NaN magnitude.
 */
inline void case_make_rotation_angle_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector axis{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  Quaternion q = make_rotation(axis, nan); // NaN quat -> normalized() -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a NaN normal must trap.
 * @details Geometry surface — rotate(normal,.).normalized() is the first strict
 *          normalize in the basis construction and traps on the NaN-poisoned
 *          vector rather than returning a garbage frame.
 */
inline void case_make_basis_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector normal{nan, opaque(0.0f), opaque(0.0f)};
  Basis b = make_basis(Quaternion(), normal); // NaN -> normalized() -> HS_CHECK
  if (b.u.x == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: make_rotation(from, to) with a non-unit source must trap.
 * @details Math-core surface — the d-based parallel/antiparallel branches assume
 *          |from| = |to| = 1, so a finite but non-unit input must trap at the
 *          unit-vector guard rather than silently skewing the rotation angle.
 */
inline void case_make_rotation_nonunit() {
  Vector from{opaque(2.0f), opaque(0.0f), opaque(0.0f)}; // |from| = 2
  Vector to{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  Quaternion q = make_rotation(from, to); // |from| != 1 -> HS_CHECK
  if (q.r == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a live-source Driver built with a null speed pointer must trap.
 * @details Animation surface — the guard traps rather than dereferencing the
 *          null pointer in the member-init list.
 */
inline void case_driver_null_speed_src() {
  static float mutant = 0.0f;
  Animation::Driver d(mutant, opaque<const float *>(nullptr), 1.0f); // -> HS_CHECK
  (void)d;
  if (mutant == 42.0f)
    std::printf("x");
}

/**
 * @brief Concrete Effect for the canvas/scan death cases.
 * @details Fixed 32x16 so its trig LUTs match the well-exercised host config;
 *          exposes registerParam (via reg) and set_clip.
 */
struct DeathEffect : public Effect {
  /**
   * @brief Constructs the effect at the fixed 32x16 resolution.
   */
  DeathEffect() : Effect(32, 16) {}
  /**
   * @brief Draws one frame (no-op; the death cases never render).
   */
  void draw_frame() override {}
  /**
   * @brief Reports whether the effect paints a background.
   * @return Always false; no background is drawn.
   */
  bool show_bg() const override { return false; }
  /**
   * @brief Registers a parameter over the unit range, exposing registerParam.
   * @param n Parameter name.
   * @param p Pointer to the backing float storage.
   */
  void reg(const char *n, float *p) { registerParam(n, p, 0.0f, 1.0f); }
};

/**
 * @brief Death case: overflowing the fixed 32-slot ParamList must trap.
 * @details Canvas surface — registerParam traps rather than silently dropping a
 *          registration, which would desync the GUI and, on WASM, break the
 *          no-realloc memory-view invariant.
 */
inline void case_register_param_overflow() {
  DeathEffect fx;
  static float slot = 0.0f;
  for (int i = 0; i < opaque(64); ++i) // exceeds capacity 32 -> HS_CHECK
    fx.reg("p", &slot);
}

/**
 * @brief Death case: a scan clip whose x_end exceeds W must trap.
 * @details Scan surface — such a clip would index the trig LUTs out of bounds;
 *          Scan::Shader::draw enforces the LUT-domain invariant once per draw
 *          (not per pixel) before the loop runs.
 */
inline void case_scan_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  fx.set_clip(0, H, 0, opaque(W + 64)); // x_end > W -> LUT-domain HS_CHECK
  Canvas c(fx);
  Scan::Shader::draw<W, H, 1>(
      c, [](const Vector &) { return Color4(Pixel(0, 0, 0), 1.0f); });
}

/**
 * @brief Death case: a face vertex index past the edge-dedup bitset must trap.
 * @details Plot surface — a vertex index beyond the TriangularBitset<128>
 *          capacity makes the face-walk draw() overload trap on the cold
 *          per-edge setup path instead of silently dropping the edge, which
 *          would leave a wireframe with missing lines and mask the sizing bug.
 */
inline void case_plot_mesh_vertex_over_capacity() {
  constexpr int W = 32, H = 16;
  /**
   * @brief Minimal duck-typed mesh: one 2-gon face whose second index (130)
   *        exceeds the bitset capacity.
   * @details The trap fires before any vertex or pipeline access, so the vertex
   *          store only needs to satisfy the interface.
   */
  struct MockMesh {
    /**
     * @brief Stand-in vertex store satisfying the mesh interface.
     */
    struct Verts {
      /**
       * @brief Returns a fixed vertex for any index.
       * @param Unused vertex index.
       * @return A constant Vector{0,1,0}.
       */
      Vector operator[](size_t) const { return Vector{0.0f, 1.0f, 0.0f}; }
      /**
       * @brief Reports the vertex count.
       * @return Always 1.
       */
      size_t size() const { return 1; }
    } vertices;
    uint8_t fc[1];  /**< Face-counts data: a single 2-gon face. */
    uint16_t fi[2]; /**< Face-index data; second entry is over-capacity. */
    /**
     * @brief Builds the mock mesh with one over-capacity 2-gon face.
     * @details Stores the over-capacity index at runtime so the optimizer can't
     *          prove the trap at compile time and reshape the case (see opaque).
     */
    MockMesh() : fc{2}, fi{0, opaque<uint16_t>(130)} {}
    /**
     * @brief Returns the face-counts array.
     * @return Pointer to the face-counts data.
     */
    const uint8_t *get_face_counts_data() const { return fc; }
    /**
     * @brief Returns the number of faces.
     * @return Always 1.
     */
    size_t get_face_counts_size() const { return 1; }
    /**
     * @brief Returns the flat face-index array.
     * @return Pointer to the face-index data.
     */
    const uint16_t *get_faces_data() const { return fi; }
  } mesh;
  DeathEffect fx;
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::Mesh::draw<W, H>(pipe, c, mesh,
                         [](const Vector &, Fragment &) {}); // index 130 -> trap
}

/**
 * @brief Death case: extract_edges with an over-capacity vertex index must trap.
 * @details Plot surface — the precomputed-edge path traps on the same cold setup
 *          path as the face-walk draw() overload, rather than silently filtering
 *          the edge out (which would produce an edge list with missing lines and
 *          mask the sizing bug).
 */
inline void case_plot_extract_edges_vertex_over_capacity() {
  /**
   * @brief Minimal duck-typed mesh with one over-capacity 2-gon face.
   * @details Same shape as the draw() case (second index 130 > 128).
   */
  struct MockMesh {
    /**
     * @brief Stand-in vertex store satisfying the mesh interface.
     */
    struct Verts {
      /**
       * @brief Returns a fixed vertex for any index.
       * @param Unused vertex index.
       * @return A constant Vector{0,1,0}.
       */
      Vector operator[](size_t) const { return Vector{0.0f, 1.0f, 0.0f}; }
      /**
       * @brief Reports the vertex count.
       * @return Always 1.
       */
      size_t size() const { return 1; }
    } vertices;
    uint8_t fc[1];  /**< Face-counts data: a single 2-gon face. */
    uint16_t fi[2]; /**< Face-index data; second entry is over-capacity. */
    /**
     * @brief Builds the mock mesh with one over-capacity 2-gon face.
     */
    MockMesh() : fc{2}, fi{0, opaque<uint16_t>(130)} {}
    /**
     * @brief Returns the face-counts array.
     * @return Pointer to the face-counts data.
     */
    const uint8_t *get_face_counts_data() const { return fc; }
    /**
     * @brief Returns the number of faces.
     * @return Always 1.
     */
    size_t get_face_counts_size() const { return 1; }
    /**
     * @brief Returns the flat face-index array.
     * @return Pointer to the face-index data.
     */
    const uint16_t *get_faces_data() const { return fi; }
  } mesh;
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(scratch_arena_a, 8);
  Plot::Mesh::extract_edges(mesh, edges); // index 130 -> trap
}

/**
 * @brief Death case: a feedback downsample that doesn't divide the resolution must trap.
 * @details Filter surface — Pixel::Feedback::flush traps rather than silently
 *          turning the whole feedback effect into a no-op; a cold
 *          authoring/config error the project routes to HS_CHECK (enabled_
 *          remains the supported way to switch feedback off). The trap fires
 *          before any_pixel_lit / scratch allocation, so no buffers needed.
 */
inline void case_feedback_downsample_indivisible() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  Canvas c(fx);
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.downsample = opaque(5); // 32 % 5 != 0 -> HS_CHECK
  Filter::Pixel::Feedback<W, H> fb(style);
  fb.flush(
      c,
      ScreenTrailFn(
          [](float, float, float) { return Color4(Pixel(0, 0, 0), 0.0f); }),
      1.0f, [](float, float, const ::Pixel &, float, float) {});
}

/**
 * @brief Death case: Path::append_segment with zero samples must trap.
 * @details Animation surface — the t / samples term divides by zero
 *          (easing(0/0) = NaN) and the loop would silently append a garbage
 *          point; the samples >= 1 guard traps the authoring error on the cold
 *          path-construction seam instead.
 */
inline void case_path_append_zero_samples() {
  Path<32> path;
  path.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f,
                      opaque(0.0f),
                      [](float t) { return t; }); // samples < 1 -> HS_CHECK
}

/**
 * @brief A named death case selected by HS_DEATH_CASE in the child process.
 */
struct Case {
  const char *name; /**< Case selector matched against HS_DEATH_CASE. */
  void (*fn)();     /**< The trap-triggering case body to run. */
};

/**
 * @brief Returns the full death-case table.
 * @param n Out-param set to the number of cases in the table.
 * @return Pointer to the static case array.
 * @details Single source of truth shared by the child dispatcher and the
 *          parent's per-case spawn loop.
 */
inline const Case *all_cases(int &n) {
  static const Case cases[] = {
      {"arena_oom", case_arena_oom},
      {"arena_set_offset_overflow", case_arena_set_offset_overflow},
      {"arena_vector_overflow", case_arena_vector_overflow},
      {"normalize_zero", case_normalize_zero},
      {"normalize_nan", case_normalize_nan},
      {"solids_index_oob", case_solids_index_oob},
      {"solids_unknown_name", case_solids_unknown_name},
      {"circular_buffer_oob", case_circular_buffer_oob},
      {"circular_buffer_front_empty", case_circular_buffer_front_empty},
      {"arena_vector_append_bulk_overflow",
       case_arena_vector_append_bulk_overflow},
      {"spatial_knn_over_max", case_spatial_knn_over_max},
      {"arena_oversubscribed", case_arena_oversubscribed},
      {"timeline_handled_relocation", case_timeline_handled_relocation},
      {"timeline_handled_completion", case_timeline_handled_completion},
      {"timeline_double_construct", case_timeline_double_construct},
      {"mesh_narrow_index", case_mesh_narrow_index},
      {"slerp_nan", case_slerp_nan},
      {"make_rotation_vectors_nan", case_make_rotation_vectors_nan},
      {"make_rotation_angle_nan", case_make_rotation_angle_nan},
      {"make_rotation_nonunit", case_make_rotation_nonunit},
      {"make_basis_nan", case_make_basis_nan},
      {"driver_null_speed_src", case_driver_null_speed_src},
      {"path_append_zero_samples", case_path_append_zero_samples},
      {"register_param_overflow", case_register_param_overflow},
      {"scan_clip_out_of_bounds", case_scan_clip_out_of_bounds},
      {"plot_mesh_vertex_over_capacity", case_plot_mesh_vertex_over_capacity},
      {"plot_extract_edges_vertex_over_capacity",
       case_plot_extract_edges_vertex_over_capacity},
      {"feedback_downsample_indivisible",
       case_feedback_downsample_indivisible},
  };
  n = static_cast<int>(sizeof(cases) / sizeof(cases[0]));
  return cases;
}

/**
 * @brief Child entry point: runs exactly one named death case, then returns.
 * @param name Case selector; an unknown name (e.g. the "__spawn_check__"
 *             control) simply returns, so the child exits 0.
 * @details Called from main() when HS_DEATH_CASE is set. The case is expected to
 *          trap before returning; returning means it did NOT trap, so the child
 *          exits 0 and the parent flags it.
 */
inline void run_child_case(const char *name) {
#if defined(_WIN32)
  SetErrorMode(0x0001u | 0x0002u);
#endif
  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i)
    if (std::strcmp(cs[i].name, name) == 0) {
      cs[i].fn();
      return;
    }
}

/**
 * @brief Sets HS_DEATH_CASE in this process's environment.
 * @param name Case selector to publish; children inherit it via system().
 */
inline void set_case_env(const char *name) {
#if defined(_WIN32)
  _putenv_s("HS_DEATH_CASE", name);
#else
  setenv("HS_DEATH_CASE", name, 1);
#endif
}

/**
 * @brief Spawns the test binary as a child running the given death case.
 * @param name Case selector passed to the child via HS_DEATH_CASE.
 * @return The raw std::system() status of the child.
 * @details Child stdout/stderr are discarded; only the exit code matters.
 */
inline int spawn_child(const char *name) {
  set_case_env(name);
  std::string cmd;
#if defined(_WIN32)
  // cmd.exe /c strips one outer quote pair, so wrap the whole command again.
  cmd = "\"\"";
  cmd += self_exe();
  cmd += "\" >NUL 2>&1\"";
#else
  cmd = "\"";
  cmd += self_exe();
  cmd += "\" >/dev/null 2>&1";
#endif
  return std::system(cmd.c_str());
}

#if defined(_WIN32)
/**
 * @brief The Windows EXCEPTION_ILLEGAL_INSTRUCTION process exit code.
 * @details An unhandled trap sets it as the process exit code, which
 *          std::system() returns through cmd.exe.
 */
inline constexpr int kTrapStatus = static_cast<int>(0xC000001D);
#endif

/**
 * @brief Tests whether a std::system() status means the child died by the trap.
 * @param rc The raw std::system() return value to interpret.
 * @return True iff the child died by the specific illegal-instruction trap.
 * @details clang lowers __builtin_trap() to an illegal instruction, so a fired
 *          HS_CHECK kills the child with SIGILL (POSIX) /
 *          STATUS_ILLEGAL_INSTRUCTION (Windows). Requiring that exact status,
 *          rather than any nonzero exit, prevents an unrelated crash or an
 *          ordinary nonzero return from being misread as a passing death test.
 */
inline bool child_trapped(int rc) {
#if defined(_WIN32)
  return rc == kTrapStatus;
#else
  // Normally /bin/sh exec-optimizes the single child command, so SIGILL
  // propagates directly and WIFSIGNALED holds. A shell that instead forks and
  // waits relays the death as an ordinary exit with status 128+SIGILL; accept
  // that shape too so the death tests don't all fail under such a shell.
  if (rc == -1)
    return false;
  if (WIFSIGNALED(rc) && WTERMSIG(rc) == SIGILL)
    return true;
  return WIFEXITED(rc) && WEXITSTATUS(rc) == 128 + SIGILL;
#endif
}

/**
 * @brief Tests whether the child exited cleanly (exit code 0).
 * @param rc The raw std::system() return value to interpret.
 * @return True iff the child exited normally with status 0.
 * @details Used by the control spawn check.
 */
inline bool child_exited_clean(int rc) {
#if defined(_WIN32)
  return rc == 0;
#else
  return rc != -1 && WIFEXITED(rc) && WEXITSTATUS(rc) == 0;
#endif
}

/**
 * @brief Reports whether the suite is running under CI.
 * @return True iff the CI environment variable is set and non-empty.
 * @details GitHub Actions (and most CI providers) set CI=true. Under CI a death
 *          suite that cannot run must FAIL loudly, not skip silently.
 */
inline bool in_ci() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *ci = std::getenv("CI");
#pragma clang diagnostic pop
  return ci && ci[0] != '\0';
}

/**
 * @brief Reports that the death suite could not run.
 * @param why Human-readable reason the suite is unrunnable.
 * @param rc The associated child return code, for diagnostics.
 * @details Loud (counts as a failure) under CI, quiet skip otherwise.
 */
inline void report_unrunnable(const char *why, int rc) {
  if (in_ci()) {
    std::printf("  [FAIL] death tests: %s (rc=%d, CI=on)\n", why, rc);
    HS_EXPECT_TRUE(false && "death suite must run under CI");
  } else {
    std::printf("  [skip] death tests: %s (rc=%d)\n", why, rc);
  }
}

/**
 * @brief Parent entry point for the death module.
 * @return The module's failure count.
 * @details Spawn-checks the harness, then runs every case in a child and asserts
 *          each died by the exact trap status.
 */
inline int run_death_tests() {
  auto scope = hs_test::begin_module("death");

  if (!self_exe() || self_exe()[0] == '\0') {
    report_unrunnable("no argv[0] to re-exec", 0);
    return hs_test::end_module(scope);
  }

  // Control: a child given an unknown case must exit cleanly. If it doesn't,
  // this harness can't reliably spawn itself here (e.g. a sandbox) — skip
  // (or FAIL under CI) rather than emit false results for every case.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    report_unrunnable("cannot re-exec self", control);
    set_case_env("");
    return hs_test::end_module(scope);
  }

  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i) {
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc);
    HS_EXPECT_TRUE(trapped);
    std::printf("  [%s] trap fires: %-26s (child rc=%d)\n",
                trapped ? "ok" : "FAIL", cs[i].name, rc);
  }

  set_case_env(""); // leave the env clean for anything that runs after us
  return hs_test::end_module(scope);
}

} // namespace death
} // namespace hs_test
