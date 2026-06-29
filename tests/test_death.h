/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Death tests for the fail-fast (HS_CHECK / __builtin_trap) seams.
 *
 * An HS_CHECK violation traps and aborts the whole process, so the in-process
 * HS_EXPECT_* harness cannot catch it. Each trap is exercised in a CHILD
 * process: the test binary re-exec's itself with HS_DEATH_CASE=<name> (handled
 * in main() before any module runs), runs exactly one trap-triggering case, and
 * the parent asserts the child died by the *specific* trap status — clang lowers
 * __builtin_trap() to an illegal instruction (x86 ud2), so the child dies by
 * SIGILL (POSIX) / STATUS_ILLEGAL_INSTRUCTION (Windows).
 *
 * The child is selected through an inherited env var and spawned shell-free —
 * fork()+execv() on POSIX, _spawnv() on Windows — so no shell can mangle the
 * re-exec path. A control "spawn check" runs first; if the harness cannot
 * re-exec itself, the death tests are SKIPPED — EXCEPT under CI (the CI env var
 * is set), where a suite that cannot run is a hard FAILURE rather than a silent
 * green skip.
 */
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <string>

#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include "core/3dmath.h"
#include "core/animation.h"
#include "core/canvas.h"
#include "core/color.h"
#include "core/geometry.h"
#include "core/filter.h"
#include "core/memory.h"
#include "core/mesh.h"
#include "core/plot.h"
#include "core/scan.h"
#include "core/solids.h"
#include "core/spatial.h"
#include "core/static_circular_buffer.h"
#include "core/transformers.h"
#include "hardware/hd107s_frame.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <fcntl.h>    // open / O_WRONLY for the /dev/null redirect
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#include <unistd.h>   // fork / execv / dup2 / close / _exit — shell-free spawn
#else
#include <fcntl.h>   // _O_WRONLY for the NUL redirect
#include <io.h>      // _dup / _dup2 / _sopen_s / _close
#include <process.h> // _spawnv / _P_WAIT — shell-free child spawn
#include <share.h>   // _SH_DENYNO
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
 * @brief Death case: non-LIFO ScratchScope teardown must trap.
 * @details The scratch-arena sharing contract between Pixel::Feedback::flush and
 *          Plot::rasterize is safe because scratch_arena_a is a LIFO bump
 *          allocator — but only while scopes are torn down in stack order.
 *          ~ScratchScope enforces that: an outer scope rewinding while an inner
 *          one is still live leaves the arena offset below the inner's saved
 *          mark, and the inner's destructor HS_CHECKs offset >= saved_offset.
 *          Here the outer scope is destroyed first (std::optional::reset),
 *          rewinding to 0; destroying the inner then sees offset 0 < its saved
 *          mark and traps.
 */
inline void case_scratch_scope_non_lifo() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  std::optional<ScratchScope> outer;
  outer.emplace(a);                  // saves offset 0
  a.allocate(opaque<size_t>(8));     // offset -> 8
  std::optional<ScratchScope> inner;
  inner.emplace(a);                  // saves offset 8
  outer.reset();                     // non-LIFO: rewinds offset to 0
  inner.reset();                     // offset 0 < saved 8 -> HS_CHECK
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
  configure_arenas(opaque(GLOBAL_ARENA_SIZE), opaque<size_t>(1024),
                   opaque<size_t>(1024));
}

/**
 * @brief Trivial Cloneable whose clone() allocates from the destination arena,
 *        so a Persist restore measurably grows the persistent arena.
 */
struct PersistProbe {
  uint8_t *storage = nullptr; /**< Stand-in for arena-backed object state. */
  /**
   * @brief Clones by allocating fresh storage from @p arena.
   * @param src Source probe (unused beyond the Cloneable interface).
   * @param dst Destination probe receiving freshly allocated storage.
   * @param arena Arena the clone allocates from.
   */
  static void clone(const PersistProbe &src, PersistProbe &dst, Arena &arena) {
    (void)src;
    dst.storage = static_cast<uint8_t *>(arena.allocate(opaque<size_t>(8)));
  }
};

/**
 * @brief Death case: a Persist scope that forgets persistent_arena.reset() must trap.
 * @details Memory surface — without the rewind, ~Persist's restore clones the
 *          backup *after* the still-live object instead of over it, pushing the
 *          persistent offset past the construction watermark; the post-restore
 *          HS_CHECK fires.
 */
inline void case_persist_forgot_reset() {
  static uint8_t pbuf[256];
  static uint8_t sbuf[256];
  Arena persistent(pbuf, sizeof(pbuf));
  Arena scratch(sbuf, sizeof(sbuf));
  PersistProbe target;
  PersistProbe::clone(target, target, persistent); // the live object in persistent
  {
    Persist<PersistProbe> p(target, scratch, persistent);
    // A correct scope rewinds here (persistent.reset()); omitting it makes the
    // restore append past the watermark.
  } // ~Persist restore -> offset past watermark -> HS_CHECK
  if (target.storage == reinterpret_cast<uint8_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a swapped (unordered) TriangularBitset pair must trap.
 * @details Memory-safety surface — index() requires small < large < MAX_V; a
 *          swapped pair would alias the wrong bit and an out-of-range one would
 *          write adjacent memory, so the HS_CHECK traps the misuse on the cold
 *          edge-dedup setup path.
 */
inline void case_triangular_bitset_unordered_pair() {
  TriangularBitset<128> bits;
  bool hit = bits.test(opaque(5), opaque(3)); // small > large -> HS_CHECK
  if (hit)
    std::printf("x");
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
    bool strobe_columns() const override { return false; }
  };
  static DeathEffect fx;
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  // pin=true marks the event handled; a 1-frame Transition is finite and the
  // sole event, so step() routes it through the completion/destroy branch.
  tl.add_get(0, Animation::Transition(v, 1.0f, 1, ease_linear), /*pin=*/true);
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
 * @brief Death case: an active noise_transform fed a non-finite direction must
 *        trap, not propagate NaN/Inf into the rendered geometry.
 * @details Transformer surface — noise_transform's active path ends in
 *          (v + distortion).normalized(), whose zero/non-finite-length guard
 *          traps. A non-finite input direction is a logic bug (every caller feeds
 *          unit sphere directions), so it must fail fast here rather than emit a
 *          NaN dot somewhere downstream. The zero-amplitude short-circuit is the
 *          legitimate no-op and is covered in-process by test_transformers.h.
 */
inline void case_noise_transform_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  NoiseParams p;
  p.amplitude = opaque(0.5f); // active path (skips the zero-amplitude no-op)
  p.scale = opaque(4.0f);
  p.time = opaque(1.0f);
  Vector v{nan, opaque(0.0f), opaque(0.0f)};
  Vector r = noise_transform(v, p); // NaN -> normalized() -> HS_CHECK
  if (r.x == 42.0f)
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
  bool strobe_columns() const override { return false; }
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
 * @brief Minimal duck-typed mesh: one 2-gon face whose second index (130)
 *        exceeds the TriangularBitset<128> capacity. Shared by both the
 *        face-walk draw() and the extract_edges over-capacity death cases so the
 *        mock interface is defined once, not kept in sync across two copies.
 * @details The trap fires before any vertex or pipeline access, so the vertex
 *          store only needs to satisfy the interface.
 */
struct OverCapacityMockMesh {
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
  OverCapacityMockMesh() : fc{2}, fi{0, opaque<uint16_t>(130)} {}
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
  /**
   * @brief Returns the flat face-index array length.
   * @return Always 2.
   */
  size_t get_faces_size() const { return 2; }
};

/**
 * @brief Death case: a face vertex index past the edge-dedup bitset must trap.
 * @details Plot surface — a vertex index beyond the TriangularBitset<128>
 *          capacity makes the face-walk draw() overload trap on the cold
 *          per-edge setup path instead of silently dropping the edge, which
 *          would leave a wireframe with missing lines and mask the sizing bug.
 */
inline void case_plot_mesh_vertex_over_capacity() {
  constexpr int W = 32, H = 16;
  OverCapacityMockMesh mesh;
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
  OverCapacityMockMesh mesh;
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
 * @details Animation surface — a zero sample count divides by zero in the
 *          t / samples term (easing(0/0) = NaN) and the loop would silently
 *          append a garbage point; the samples >= 1 guard traps the authoring
 *          error on the cold path-construction seam instead.
 */
inline void case_path_append_zero_samples() {
  Path<32> path;
  path.append_segment([](float s) { return Vector(s, 0.0f, 0.0f); }, 1.0f,
                      opaque(0),
                      [](float t) { return t; }); // samples < 1 -> HS_CHECK
}

/**
 * @brief Death case: a Gradient stop position outside [0,1] must trap.
 * @details Color surface — a stop position indexes entries[256] via
 *          static_cast<int>(pos * 255); pos > 1 (or < 0) is an out-of-bounds
 *          table write. The constructor traps the authoring error always-on at
 *          the cold literal-construction seam rather than corrupting memory.
 */
inline void case_gradient_stop_out_of_range() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)},
                {1.5f, CPixel(255u, 255u, 255u)}}; // pos > 1 -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: descending (unsorted) Gradient stops must trap.
 * @details Color surface — segments are only filled when end > start, so a
 *          transposed/unsorted pair would silently degenerate to wrong output.
 *          The constructor requires ascending positions and traps otherwise.
 */
inline void case_gradient_stops_unsorted() {
  Gradient grad{{0.6f, CPixel(0u, 0u, 0u)},
                {0.3f, CPixel(255u, 255u, 255u)}}; // descending -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: a RandomTimer with min > max must trap.
 * @details Animation surface — reset() draws hs::rand_int(min, max + 1), a
 *          half-open range that is empty/inverted when min > max, giving an
 *          implementation-defined garbage delay. The constructor traps the
 *          inverted (or negative) range at the cold authoring seam.
 */
inline void case_random_timer_inverted_range() {
  Animation::RandomTimer timer(opaque(5), opaque(2),
                               [](Canvas &) {}); // min > max -> HS_CHECK
  (void)timer;
}

/**
 * @brief Death case: an out-of-range load() count must trap.
 * @details Hardware wire-format surface — HD107SFrame::load formerly clamped a
 *          count > N (and no-op'd a negative count), masking a caller sizing
 *          bug at this cold bind seam; it now HS_CHECKs count in [0, N].
 */
inline void case_hd107s_load_count_over_range() {
  static HD107SFrame<40> frame;
  static CRGB src[40] = {};
  frame.load(src, opaque(40 + 8)); // count > N -> HS_CHECK
  if (frame.data()[0] == 0x42)
    std::printf("x");
}

/**
 * @brief Death case: calling an empty (default-constructed) Fn must trap.
 * @details Concepts surface — hs::inplace_function routes an empty-state call
 *          through ipf_empty_ops::invoke, which fail-fast traps via check_fail
 *          rather than dereferencing the empty buffer (std::function would throw
 *          bad_function_call; the engine builds without exceptions). The
 *          never-taken opaque(false) assignment keeps the optimizer from proving
 *          the function empty and folding the trap at compile time. The non-trap
 *          value semantics (copy/move/empty operator bool) are covered in-process
 *          by tests/test_concepts.h.
 */
inline void case_empty_fn_call() {
  Fn<int(int), 16> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  if (v == 42)
    std::printf("x");
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
      {"scratch_scope_non_lifo", case_scratch_scope_non_lifo},
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
      {"persist_forgot_reset", case_persist_forgot_reset},
      {"triangular_bitset_unordered_pair",
       case_triangular_bitset_unordered_pair},
      {"timeline_handled_relocation", case_timeline_handled_relocation},
      {"timeline_handled_completion", case_timeline_handled_completion},
      {"timeline_double_construct", case_timeline_double_construct},
      {"mesh_narrow_index", case_mesh_narrow_index},
      {"slerp_nan", case_slerp_nan},
      {"make_rotation_vectors_nan", case_make_rotation_vectors_nan},
      {"make_rotation_angle_nan", case_make_rotation_angle_nan},
      {"make_rotation_nonunit", case_make_rotation_nonunit},
      {"make_basis_nan", case_make_basis_nan},
      {"noise_transform_nan", case_noise_transform_nan},
      {"driver_null_speed_src", case_driver_null_speed_src},
      {"path_append_zero_samples", case_path_append_zero_samples},
      {"register_param_overflow", case_register_param_overflow},
      {"scan_clip_out_of_bounds", case_scan_clip_out_of_bounds},
      {"plot_mesh_vertex_over_capacity", case_plot_mesh_vertex_over_capacity},
      {"plot_extract_edges_vertex_over_capacity",
       case_plot_extract_edges_vertex_over_capacity},
      {"feedback_downsample_indivisible",
       case_feedback_downsample_indivisible},
      {"gradient_stop_out_of_range", case_gradient_stop_out_of_range},
      {"gradient_stops_unsorted", case_gradient_stops_unsorted},
      {"random_timer_inverted_range", case_random_timer_inverted_range},
      {"hd107s_load_count_over_range", case_hd107s_load_count_over_range},
      {"empty_fn_call", case_empty_fn_call},
  };
  n = static_cast<int>(sizeof(cases) / sizeof(cases[0]));
  return cases;
}

/**
 * @brief Dedicated always-trapping case used only to probe the trap-relay shape.
 * @details Not part of all_cases(): run_child_case() dispatches it directly. It
 *          traps through the same HS_CHECK path as every real case, so its relay
 *          shape matches theirs, but it can never regress to not-trapping the way
 *          a real case might — so shape detection never rests on a real case.
 */
inline constexpr const char *kShapeProbeCase = "__shape_probe__";

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
  if (std::strcmp(name, kShapeProbeCase) == 0) {
    HS_CHECK(false, "death-harness trap-shape probe"); // always traps
    return;
  }
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
 * @param name Case selector to publish; the spawned child inherits it through
 *             the environment.
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
 * @return The child's raw status (_spawnv() on Windows, fork+execv wait status
 *         on POSIX). -1 on a spawn failure.
 * @details Child stdout/stderr are discarded; only the exit code matters.
 */
inline int spawn_child(const char *name) {
  set_case_env(name);
#if defined(_WIN32)
  // Shell-free spawn: hand argv straight to the CRT so no cmd.exe parsing can
  // mangle a self_exe() path containing &, %, ^, or quotes. Child stdout/stderr
  // go to NUL by redirecting fds 1/2 across the synchronous _P_WAIT spawn (the
  // child inherits the CRT fd table), then the descriptors are restored.
  std::fflush(stdout);
  std::fflush(stderr);
  int saved_out = _dup(1);
  int saved_err = _dup(2);
  int devnull = -1;
  // Redirect to NUL only when both originals were saved: a failed _dup leaves no
  // way to restore, so skipping the redirect keeps the parent's streams intact
  // (a muted parent would silence reporting for every remaining case).
  if (saved_out >= 0 && saved_err >= 0) {
    _sopen_s(&devnull, "NUL", _O_WRONLY, _SH_DENYNO, 0);
    if (devnull >= 0) {
      _dup2(devnull, 1);
      _dup2(devnull, 2);
    }
  }
  const char *argv[] = {self_exe(), nullptr};
  intptr_t rc = _spawnv(_P_WAIT, self_exe(), argv);
  // devnull is only open when both saves succeeded, so the restore is reached
  // only with valid descriptors.
  if (devnull >= 0) {
    _dup2(saved_out, 1);
    _dup2(saved_err, 2);
    _close(devnull);
  }
  if (saved_out >= 0)
    _close(saved_out);
  if (saved_err >= 0)
    _close(saved_err);
  return static_cast<int>(rc);
#else
  // Shell-free spawn: fork and execv the binary directly so no /bin/sh parsing
  // can mangle a self_exe() path containing a quote or shell metacharacter. The
  // child sends stdout/stderr to /dev/null and execs; the parent waits and
  // returns the raw wait status that classify_trap() decodes.
  std::fflush(stdout);
  std::fflush(stderr);
  const char *exe = self_exe();
  pid_t pid = fork();
  if (pid < 0)
    return -1;
  if (pid == 0) {
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull >= 0) {
      dup2(devnull, 1);
      dup2(devnull, 2);
      if (devnull > 2)
        close(devnull);
    }
    const char *argv[] = {exe, nullptr};
    execv(exe, const_cast<char *const *>(argv));
    _exit(127); // exec failed — never returns to the harness
  }
  int status = 0;
  if (waitpid(pid, &status, 0) < 0)
    return -1;
  return status;
#endif
}

#if defined(_WIN32)
/**
 * @brief The Windows EXCEPTION_ILLEGAL_INSTRUCTION process exit code.
 * @details An unhandled trap sets it as the process exit code, which _spawnv()
 *          returns directly to the parent.
 */
inline constexpr int kTrapStatus = static_cast<int>(0xC000001D);
#endif

/**
 * @brief How a child's illegal-instruction trap reaches the parent.
 * @details clang lowers __builtin_trap() to an illegal instruction, so a fired
 *          HS_CHECK kills the child with SIGILL (POSIX) /
 *          STATUS_ILLEGAL_INSTRUCTION (Windows). The harness spawns the child
 *          shell-free, so the trap normally arrives as that signal directly
 *          (Signal). Exit128 is kept defensively for any relay that instead
 *          surfaces the death as an ordinary exit with status 128+SIGILL: that
 *          is indistinguishable from a raw `exit(128+SIGILL)` at the wait-status
 *          level, so accepting BOTH unconditionally would let a child that
 *          genuinely exit(132)s read as a passing death test. The harness
 *          therefore probes which shape occurs once (run_death_tests), with a
 *          dedicated always-trapping sentinel, and then requires exactly that
 *          shape per case.
 */
enum class TrapShape { None, Signal, Exit128 };

/**
 * @brief Classifies a child wait status into the trap relay shape, if any.
 * @param rc The raw spawn_child() return value to interpret.
 * @return The TrapShape the status represents, or TrapShape::None.
 */
inline TrapShape classify_trap(int rc) {
#if defined(_WIN32)
  return rc == kTrapStatus ? TrapShape::Signal : TrapShape::None;
#else
  if (rc == -1)
    return TrapShape::None;
  if (WIFSIGNALED(rc) && WTERMSIG(rc) == SIGILL)
    return TrapShape::Signal;
  if (WIFEXITED(rc) && WEXITSTATUS(rc) == 128 + SIGILL)
    return TrapShape::Exit128;
  return TrapShape::None;
#endif
}

/**
 * @brief Tests whether a child died by the trap in the probed relay shape.
 * @param rc The raw spawn_child() return value to interpret.
 * @param expected The TrapShape the harness probed (never None).
 * @return True iff the child died by exactly that illegal-instruction relay.
 * @details Requiring the single probed shape — not "either signal OR 128+sig" —
 *          keeps the guarantee tight: a child that exit(128+SIGILL)s under a
 *          direct-relay environment no longer counts as a trap.
 */
inline bool child_trapped(int rc, TrapShape expected) {
  return expected != TrapShape::None && classify_trap(rc) == expected;
}

/**
 * @brief Tests whether the child exited cleanly (exit code 0).
 * @param rc The raw spawn_child() return value to interpret.
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
    // Count a skip — never a pass — so a green local run cannot be mistaken for
    // trap coverage; the banner is unmistakable and CI is the hard gate.
    ++hs_test::stats().skipped;
    std::printf("  [SKIPPED] death tests: %s (rc=%d) — 0 trap cases executed; "
                "CI is the hard gate\n",
                why, rc);
  }
}

/**
 * @brief Parent entry point for the death module.
 * @return The module's failure count.
 * @details Spawn-checks the harness, then runs every case in a child and asserts
 *          each died by the exact trap status.
 */
inline int run_death_tests() {
  hs_test::ModuleFixture fixture("death");

  if (!self_exe() || self_exe()[0] == '\0') {
    report_unrunnable("no argv[0] to re-exec", 0);
    return fixture.result();
  }

  // Control: a child given an unknown case must exit cleanly. If it doesn't,
  // this harness can't reliably spawn itself here (e.g. a sandbox) — skip
  // (or FAIL under CI) rather than emit false results for every case.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    report_unrunnable("cannot re-exec self", control);
    set_case_env("");
    return fixture.result();
  }

  int n;
  const Case *cs = all_cases(n);

  // Probe how a trap is relayed (direct SIGILL vs an exit 128+SIGILL) with a
  // dedicated always-trapping sentinel rather than a real case. A real case that
  // regressed to not trapping would otherwise corrupt shape detection and skip
  // the whole suite, instead of failing just that case in the loop below. The
  // sentinel traps through the same HS_CHECK path, so its shape matches the cases.
  TrapShape shape = classify_trap(spawn_child(kShapeProbeCase));
  if (shape == TrapShape::None) {
    report_unrunnable("trap-shape sentinel did not trap; cannot classify trap status",
                      0);
    set_case_env("");
    return fixture.result();
  }

  for (int i = 0; i < n; ++i) {
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc, shape);
    HS_EXPECT_TRUE(trapped);
    std::printf("  [%s] trap fires: %-26s (child rc=%d)\n",
                trapped ? "ok" : "FAIL", cs[i].name, rc);
  }

  set_case_env(""); // leave the env clean for anything that runs after us
  return fixture.result();
}

} // namespace death
} // namespace hs_test
