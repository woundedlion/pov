/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
 *
 * Dying by SIGILL alone proves only that SOMETHING trapped: under
 * -fsanitize-trap=undefined any UB in a case body lowers to the same illegal
 * instruction, so a case that stopped reaching its guard but tripped UB would
 * still look green. Each case therefore also pins the guard it must fire.
 * check_fail() logs "HS_CHECK failed: <file>:<line>: (<cond>) <msg>" and flushes
 * before trapping, so the child's stdout/stderr is captured to a file and the
 * parent requires the breadcrumb of that exact guard.
 *
 * The run closes with a coverage line: how many of the engine's HS_CHECK sites
 * a case actually pins, against every site in the tree. Both numbers are
 * derived — the denominator from a configure-time census of the sources
 * (tests/count_guard_sites.cmake), the numerator from the case table — so the
 * ratio cannot quietly disagree with either. It is a report, not a gate.
 */
#pragma once

#include <bit>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <string>

#include "death_guard_sites.h" // generated HS_CHECK census; see tests/CMakeLists.txt
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/test_shader_workbench.h" // ShaderWorkbenchWhiteBox, for the effect-side traps

#include "core/math/3dmath.h"
#include "core/animation/animation.h"
#include "core/animation/carousel.h"
#include "core/render/canvas.h"
#include "core/color/color.h"
#include "core/color/noise_hue_palette.h"
#include "core/math/geometry.h"
#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/control/registry.h"
#include "core/engine/memory.h"
#include "core/render/led.h"
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
#include "core/containers/triangular_bitset.h"
#include "core/containers/static_circular_buffer.h"
#include "core/animation/transformer.h"
#include "hardware/dma_led_controller.h"
#include "hardware/pov_sync.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <fcntl.h>    // open / O_WRONLY for the /dev/null redirect
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#include <unistd.h>   // fork / execv / dup2 / close / _exit — shell-free spawn
#else
#include <fcntl.h>   // _O_WRONLY / _O_CREAT / _O_TRUNC for the capture redirect
#include <io.h>      // _dup / _dup2 / _sopen_s / _close
#include <process.h> // _spawnv / _P_WAIT / _getpid — shell-free child spawn
#include <share.h>   // _SH_DENYNO
#include <sys/stat.h> // _S_IREAD / _S_IWRITE for the created capture file
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

/** @brief Death case: make() must preserve Arena's OOM trap contract. */
inline void case_arena_make_oom() {
  static uint8_t buf[sizeof(uint32_t) - 1];
  Arena a(buf, sizeof(buf));
  auto *p = a.make<uint32_t>(opaque<uint32_t>(7));
  if (p == reinterpret_cast<uint32_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a zero-size arena allocation must trap.
 * @details Memory surface — a zero-size request returns a bump pointer that
 *          reserves nothing and aliases the next allocation's address, so it is
 *          rejected as misuse rather than handed back as ownable storage.
 */
inline void case_arena_zero_size_alloc() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(0)); // size == 0 -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: a non-power-of-two allocation alignment must trap.
 * @details Memory surface — allocate()'s padding math is a modulo against the
 *          requested alignment, which only yields an aligned address for a
 *          power of two.
 */
inline void case_arena_bad_alignment() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(8), opaque<size_t>(3)); // -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: shrinking capacity below the live offset must trap.
 * @details Memory surface — set_capacity moves only the boundary, so a new
 *          capacity under the live offset would strand already-allocated
 *          content outside the arena instead of freeing it.
 */
inline void case_arena_set_capacity_below_offset() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.allocate(opaque<size_t>(32), 1);
  a.set_capacity(opaque<size_t>(16)); // < offset 32 -> HS_CHECK
}

/**
 * @brief Death case: growing capacity past the backing buffer must trap.
 * @details Memory surface — allocate() bounds-checks against the capacity, so a
 *          capacity beyond the buffer's real end would authorize allocations
 *          past it with no further guard.
 */
inline void case_arena_set_capacity_above_extent() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.set_capacity(opaque<size_t>(128)); // > extent 64 -> HS_CHECK
}

/**
 * @brief Death case: rebinding to a capacity past the buffer extent must trap.
 * @details Memory surface — the extent is the ceiling every later set_capacity()
 *          grow is bounded by, so a rebind that starts above it would authorize
 *          allocations past the backing buffer from the first call on.
 */
inline void case_arena_rebind_capacity_over_extent() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.rebind(buf, opaque<size_t>(128), opaque<size_t>(64)); // -> HS_CHECK
}

/**
 * @brief Death case: a mid-run resplit with live scratch content must trap.
 * @details Config surface — resplit_arenas rebases both scratch arenas, and a
 *          ScratchScope saved at offset 0 restores to 0 either way, so live
 *          scratch content would be silently rebased onto the new split.
 */
inline void case_resplit_scratch_not_empty() {
  configure_arenas_default();
  scratch_arena_a.allocate(opaque<size_t>(16));
  resplit_arenas(opaque(DEFAULT_PERSISTENT_SIZE),
                 opaque(DEFAULT_SCRATCH_A_SIZE),
                 opaque(DEFAULT_SCRATCH_B_SIZE)); // scratch live -> HS_CHECK
}

/**
 * @brief Death case: moving the arena offset forward must trap.
 * @details Memory surface — set_offset only ever rewinds. A forward move stays
 *          inside capacity yet hands back bytes already reclaimed, so the guard
 *          is monotone decrease, not a capacity bound.
 */
inline void case_arena_set_offset_forward() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  a.allocate(opaque<size_t>(8), opaque<size_t>(1));
  a.set_offset(opaque<size_t>(0));
  a.set_offset(opaque<size_t>(4)); // forward, still under capacity -> HS_CHECK
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
  outer.emplace(a);              // saves offset 0
  a.allocate(opaque<size_t>(8)); // offset -> 8
  std::optional<ScratchScope> inner;
  inner.emplace(a); // saves offset 8
  outer.reset();    // non-LIFO: rewinds offset to 0
  inner.reset();    // offset 0 < saved 8 -> HS_CHECK
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
 * @brief Death case: ArenaVector fixed-capacity emplace_back overflow must trap.
 * @details Arena-container surface — the in-place construction path carries its
 *          own capacity guard, distinct from push_back's copy path.
 */
inline void case_arena_vector_emplace_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2);
  v.emplace_back(1);
  v.emplace_back(2);
  v.emplace_back(opaque(3)); // exceeds capacity -> HS_CHECK
}

/**
 * @brief Death case: generate() with a scratch arena as its target must trap.
 * @details Generator surface — the depth-0 reset and the ScratchScope rewind
 *          would destroy output written into either engine scratch arena, so an
 *          aliasing target is rejected before the callback runs.
 */
inline void case_generate_target_is_scratch() {
  configure_arenas_default();
  int r = hs::generate(scratch_arena_a, [](Arena &, Arena &, Arena &) {
    return 0;
  }); // target aliases scratch_arena_a -> HS_CHECK
  if (r == 42)
    std::printf("x");
}

/**
 * @brief Recursive generator body: each level calls generate() once more.
 * @param target Arena forwarded to the next generate().
 * @param remaining Levels of nesting still to open; stops at 0.
 * @return Always 0.
 */
inline int nested_generate(Arena &target, Arena &, Arena &, int remaining) {
  if (remaining <= 0)
    return 0;
  return hs::generate(target, nested_generate, remaining - 1);
}

/**
 * @brief Death case: nesting generate() past MAX_GENERATE_DEPTH must trap.
 * @details Generator surface — every level stacks two ScratchScopes on a fixed
 *          scratch budget, so runaway reentrancy is capped at the wrapper rather
 *          than left to exhaust the arenas. The outermost call opens depth 1, so
 *          MAX_GENERATE_DEPTH further levels reach depth MAX_GENERATE_DEPTH + 1.
 */
inline void case_generate_recursion_too_deep() {
  configure_arenas_default();
  int r = hs::generate(persistent_arena, nested_generate,
                       opaque(hs::MAX_GENERATE_DEPTH));
  if (r == 42)
    std::printf("x");
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
 *          the bulk memcpy path has its own remaining-capacity guard.
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
  auto r =
      tree.nearest(Vector(1.0f, 0.0f, 0.0f),
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
  PersistProbe::clone(target, target,
                      persistent); // the live object in persistent
  {
    Persist<PersistProbe> p(target, scratch, persistent);
    // A correct scope rewinds here (persistent.reset()); omitting it makes the
    // restore append past the watermark.
  } // ~Persist restore -> offset past watermark -> HS_CHECK
  if (target.storage == reinterpret_cast<uint8_t *>(0x1))
    std::printf("x");
}

/**
 * @brief Cloneable payload that allocates nothing, so only the distinct-arena
 *        guard can trap a same-arena Persist.
 */
struct FlatProbe {
  int value = 0; /**< Whole payload; clone() copies it without an arena. */
  /**
   * @brief Clones by plain copy, leaving the arena offset untouched.
   * @param src Source probe.
   * @param dst Destination probe.
   * @param arena Unused; the payload needs no storage.
   */
  static void clone(const FlatProbe &src, FlatProbe &dst, Arena &arena) {
    (void)arena;
    dst.value = src.value;
  }
};

/**
 * @brief Death case: a Persist naming one arena for both roles must trap.
 * @details Memory surface — ~Persist's watermark restore assumes the backup
 *          outlives the rewind of the arena it restores into, which a single
 *          arena cannot provide. The payload allocates nothing, so the
 *          post-restore watermark check cannot fire and the distinct-arena
 *          guard is the only reachable trap.
 */
inline void case_persist_same_arena() {
  static uint8_t pbuf[256];
  Arena persistent(pbuf, sizeof(pbuf));
  FlatProbe target;
  target.value = opaque(7);
  Persist<FlatProbe> p(target, persistent, persistent); // -> HS_CHECK
  if (target.value == 42)
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
 *          handed out via add_get(Pin::PINNED), converting the dangling-handle
 *          hazard into a fail-fast crash instead of silent corruption.
 */
inline void case_timeline_handled_relocation() {
  TimelineEvent src;
  src.handled = opaque(true); // as if handed out by add_get(Pin::PINNED)
  TimelineEvent dst;
  src.move_into(dst); // HS_CHECK(!handled) -> trap
}

/**
 * @brief Death case: relocating into a slot that still owns an animation must
 *        trap.
 * @details Animation surface — move_into overwrites dst.manager/dst.iface, so a
 *          live destination would lose its animation's destructor. step()'s
 *          compaction only ever targets slots it has already vacated; the trap
 *          pins that invariant for every relocation path.
 */
inline void case_timeline_move_into_live_destination() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 10, ease_linear));
  tl.add(0, Animation::Transition(v, 1.0f, 10, ease_linear));
  global_timeline_events[opaque(0)].move_into(global_timeline_events[1]);
}

/**
 * @brief Death case: a negative timeline delay must trap.
 */
inline void case_timeline_negative_delay() {
  Timeline tl;
  float value = 0.0f;
  tl.add(opaque(-1), Animation::Transition(value, 1.0f, 1, ease_linear));
}

/**
 * @brief Death case: a timeline start past UINT32_MAX must trap.
 */
inline void case_timeline_start_overflow() {
  Timeline tl;
  global_timeline_t = opaque<uint32_t>(UINT32_MAX - 1);
  float value = 0.0f;
  tl.add(opaque(2), Animation::Transition(value, 1.0f, 1, ease_linear));
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
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  // add_get(Pin::PINNED) rejects a finite non-repeating animation up front (see
  // case_timeline_pinned_finite_animation), so the event is marked handled
  // directly to reach step()'s completion branch. A 1-frame Transition is finite
  // and the sole event, so step() routes it through completion/destroy.
  tl.add(0, Animation::Transition(v, 1.0f, 1, ease_linear));
  global_timeline_events[0].handled = opaque(true);
  tl.step(canvas); // t=1: done() && !repeats() && !canceled, keep=false -> trap
}

/**
 * @brief Death case: pinning a finite, non-repeating animation must trap.
 * @details Animation surface — add_get(Pin::PINNED) promises the caller a pointer
 *          valid across frames, which only holds for an animation that never
 *          completes on its own. The up-front check rejects the misuse at the
 *          add site instead of leaving it to step()'s completion guard, which
 *          fires only once the animation actually finishes.
 */
inline void case_timeline_pinned_finite_animation() {
  Timeline tl;
  float v = 0.0f;
  tl.add_get(0, Animation::Transition(v, 1.0f, opaque(1), ease_linear),
             Timeline::Pin::PINNED);
}

/**
 * @brief Death case: dropping a pinned add on a full timeline must trap.
 * @details Animation surface — the capacity guard returns nullptr, but an
 *          add_get(Pin::PINNED) caller retains that pointer across frames and no
 *          call site null-checks it. The guard traps on the pinned case so a
 *          full timeline fails at the add instead of at the first use of the
 *          stored handle.
 */
inline void case_timeline_pinned_add_on_full_timeline() {
  Timeline tl;
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, ease_linear));
  tl.add_get(0,
             Animation::PeriodicTimer(
                 1, [](Canvas &) {}, /*repeat=*/true),
             opaque(Timeline::Pin::PINNED));
}

/**
 * @brief Death case: a pinned one-shot timer must trap when it fires.
 * @details Animation surface — a one-shot RandomTimer/PeriodicTimer ends itself
 *          on its single trigger. Ending via finish() (not cancel()) keeps
 *          is_canceled() false, so the destroy of a pinned timer hits step()'s
 *          completion guard instead of slipping through its cancellation
 *          exemption and dangling the retained pointer silently.
 */
inline void case_timeline_pinned_one_shot_timer() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  tl.add(0, Animation::PeriodicTimer(1, [](Canvas &) {}, /*repeat=*/false));
  global_timeline_events[0].handled = opaque(true);
  tl.step(canvas); // t=1: fires, finish() -> done() && !canceled -> trap
}

/**
 * @brief Death case: clear()ing a pinned (handled) event must trap.
 * @details Animation surface — the third teardown path, alongside
 *          case_timeline_handled_relocation (move_into) and
 *          case_timeline_handled_completion (step's destroy branch). The public
 *          clear() would otherwise free an event whose animation pointer the
 *          caller still holds. ~Timeline reaches the same events through the
 *          unguarded reset_storage(), which is safe because no retained handle
 *          spans the instance boundary.
 */
inline void case_timeline_clear_pinned() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, ease_linear));
  global_timeline_events[0].handled = opaque(true);
  tl.clear(); // HS_CHECK(!handled) -> trap
}

/**
 * @brief Death case: clear()ing from a completion callback must trap.
 * @details Animation surface — step() runs post_callback() and only afterwards
 *          destroys the event, so a clear() inside that callback would free the
 *          callable whose frame is still executing. The trap sits at the top of
 *          clear(), ahead of destroy_events().
 */
inline void case_timeline_clear_during_step() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, ease_linear).then([&tl]() {
    tl.clear();
  }));
  tl.step(canvas); // t=1: completes -> callback -> clear() while stepping
}

/** @brief Death case: finite parameter animations reject the -1 sentinel. */
inline void case_finite_param_perpetual_duration() {
  float value = 0.0f;
  Animation::Transition transition(value, 1.0f, opaque(-1), ease_linear);
  if (transition.done())
    std::printf("x");
}

/** @brief Death case: a Transition target must be finite. */
inline void case_transition_nonfinite_target() {
  float value = 0.0f;
  Animation::Transition transition(
      value, opaque(std::numeric_limits<float>::quiet_NaN()), 1, ease_linear);
  if (transition.done())
    std::printf("x");
}

/** @brief Adds an event from a clear hook, violating the hook contract. */
inline void add_event_from_clear_hook(void *ctx) {
  static float value = 0.0f;
  static_cast<Timeline *>(ctx)->add(
      0, Animation::Transition(value, 1.0f, 1, ease_linear));
}

/** @brief Death case: clear hooks must not mutate timeline event storage. */
inline void case_timeline_clear_hook_adds_event() {
  Timeline tl;
  tl.add_clear_hook(&tl, add_event_from_clear_hook);
  tl.clear();
}

/** @brief Death case: a segue must target the already-flipped front slot. */
inline void case_mesh_carousel_unflipped_slot() {
  Timeline tl;
  MeshCarousel<> carousel;
  carousel.schedule_segue(tl, 1, [](Canvas &, float) {}, 4, 1);
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

/** @brief Death case: medial rejects an input used as both outputs. */
inline void case_medial_aliases_input() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::medial(mesh, mesh, mesh.vertices, target, temp);
}

/**
 * @brief Death case: needle rejects one arena passed as both target and temp.
 * @details needle is dual∘kis with the two arenas swapped between the legs, so a
 *          single arena has the second leg reading the block the first is
 *          overwriting.
 */
inline void case_needle_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::needle(mesh, arena, arena);
}

/** @brief Death case: zip rejects one arena passed as both target and temp. */
inline void case_zip_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::zip(mesh, arena, arena);
}

/** @brief Death case: gyro rejects one arena passed as both target and temp. */
inline void case_gyro_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::gyro(mesh, arena, arena);
}

/** @brief Death case: bevel rejects one arena passed as both target and temp. */
inline void case_bevel_aliased_arenas() {
  static uint8_t buf[64];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  MeshOps::bevel(mesh, arena, arena);
}

/**
 * @brief Death case: MeshOps::transform rejects a self-aliased destination.
 * @details set_borrowed() drops the source's owned topology before it is read,
 *          so a self-transform would report F faces over reclaimed arena
 *          bytes.
 */
inline void case_mesh_transform_aliases_source() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  MeshState mesh;
  MeshOps::transform(mesh, mesh, arena);
}

/** @brief Death case: chamfer rejects the face-collapse endpoint. */
inline void case_chamfer_collapsed_endpoint() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::chamfer(mesh, target, temp, opaque(1.0f));
}

/** @brief Death case: snub rejects the face-collapse endpoint. */
inline void case_snub_collapsed_endpoint() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::snub(mesh, target, temp, opaque(1.0f));
}

/** @brief Death case: a Conway morph operator rejects an empty mesh. */
inline void case_conway_empty_mesh() {
  static uint8_t target_buf[64];
  static uint8_t temp_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Death case: a Conway morph operator rejects an open one-face mesh. */
inline void case_conway_degenerate_mesh() {
  static uint8_t source_buf[1024];
  static uint8_t target_buf[4096];
  static uint8_t temp_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  mesh.vertices.bind(source, 3);
  mesh.vertices.push_back(Vector(1, 0, 0));
  mesh.vertices.push_back(Vector(0, 1, 0));
  mesh.vertices.push_back(Vector(0, 0, 1));
  mesh.face_counts.bind(source, 1);
  mesh.face_counts.push_back(3);
  mesh.faces.bind(source, 3);
  mesh.faces.push_back(0);
  mesh.faces.push_back(1);
  mesh.faces.push_back(2);
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Death case: a Conway morph operator traps on target exhaustion. */
inline void case_conway_target_exhausted() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[16];
  static uint8_t temp_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena temp(temp_buf, sizeof(temp_buf));
  PolyMesh mesh;
  build_solid<Solids::Tetrahedron>(mesh, source);
  MeshOps::truncate(mesh, target, temp, opaque(0.25f));
}

/** @brief Vertex-bit payload words a tetrahedron bake needs (3 per vertex). */
inline constexpr size_t TETRAHEDRON_BAKE_WORDS = 3 * 4;

/**
 * @brief Builds a tetrahedron plus a relax bake that matches it exactly.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param bits Payload storage, TETRAHEDRON_BAKE_WORDS words, filled with the
 *        mesh's own vertex bits; must outlive the relax_baked() call.
 * @return Bake relax_baked() accepts until the caller perturbs one field.
 * @details relax_baked does not check source positions, so the source mesh's
 *          own vertices are a legal payload.
 */
inline MeshOps::RelaxBake
build_matching_relax_bake(PolyMesh &mesh, Arena &arena, uint32_t *bits) {
  build_solid<Solids::Tetrahedron>(mesh, arena);
  uint32_t output_hash = MeshOps::FNV1A_BASIS;
  for (size_t i = 0; i < mesh.vertices.size(); ++i) {
    const Vector &v = mesh.vertices[i];
    bits[3 * i] = std::bit_cast<uint32_t>(v.x);
    bits[3 * i + 1] = std::bit_cast<uint32_t>(v.y);
    bits[3 * i + 2] = std::bit_cast<uint32_t>(v.z);
    for (size_t k = 0; k < 3; ++k)
      output_hash = MeshOps::fnv1a_step(output_hash, bits[3 * i + k]);
  }
  MeshOps::RelaxBake bake{};
  bake.name = "death_tetrahedron";
  bake.vertex_bits = bits;
  bake.vertex_count = static_cast<uint16_t>(mesh.vertices.size());
  bake.face_count = static_cast<uint16_t>(mesh.get_face_counts_size());
  bake.index_count = static_cast<uint16_t>(mesh.get_faces_size());
  bake.iterations = 0;
  bake.topology_hash = MeshOps::relax_topology_hash(mesh);
  bake.output_hash = output_hash;
  return bake;
}

/**
 * @brief Death case: relax_baked rejects a bake whose vertex count differs from
 *        the source mesh.
 * @details Baked-payload surface — the dimension check is what stops a payload
 *          baked against different geometry from being read past its end.
 */
inline void case_relax_baked_dimension_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.vertex_count = opaque<uint16_t>(bake.vertex_count + 1);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: relax_baked rejects a bake whose topology hash differs.
 * @details Baked-payload surface — dimensions alone do not pin connectivity, so
 *          this check is what stops a bake replaying onto a mesh of the same
 *          size but different face wiring.
 */
inline void case_relax_baked_topology_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.topology_hash = opaque(bake.topology_hash ^ 1u);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: relax_baked rejects a payload whose re-hash differs from
 *        the bake's output hash.
 * @details Baked-payload surface — the only check covering the vertex words
 *          themselves, so a corrupt or truncated flash payload stops here
 *          rather than shipping as geometry.
 */
inline void case_relax_baked_output_hash_mismatch() {
  static uint8_t source_buf[4096];
  static uint8_t target_buf[4096];
  Arena source(source_buf, sizeof(source_buf));
  Arena target(target_buf, sizeof(target_buf));
  uint32_t bits[TETRAHEDRON_BAKE_WORDS];
  PolyMesh mesh;
  MeshOps::RelaxBake bake = build_matching_relax_bake(mesh, source, bits);
  bake.output_hash = opaque(bake.output_hash ^ 1u);
  PolyMesh out = MeshOps::relax_baked(mesh, target, bake);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Builds a one-face PolyMesh with independently sized count and index
 * data.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param side_count Value stored in the face-count array.
 * @param num_indices Number of entries stored in the flat face-index array.
 */
inline void build_mismatched_polymesh(PolyMesh &mesh, Arena &arena,
                                      uint8_t side_count, size_t num_indices) {
  mesh.vertices.bind(arena, 4);
  for (size_t i = 0; i < 4; ++i)
    mesh.vertices.push_back(Vector{});
  mesh.face_counts.bind(arena, 1);
  mesh.face_counts.push_back(opaque(side_count));
  mesh.faces.bind(arena, num_indices);
  for (size_t i = 0; i < num_indices; ++i)
    mesh.faces.push_back(static_cast<uint16_t>(i));
}

/**
 * @brief Death case: HalfEdgeMesh rejects counts shorter than its flat index
 * list.
 */
inline void case_half_edge_face_counts_short() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  build_mismatched_polymesh(mesh, arena, 3, opaque<size_t>(4));
  HalfEdgeMesh half_edges(arena, mesh);
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: HalfEdgeMesh rejects counts longer than its flat index
 * list.
 */
inline void case_half_edge_face_counts_long() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh mesh;
  build_mismatched_polymesh(mesh, arena, 4, opaque<size_t>(3));
  HalfEdgeMesh half_edges(arena, mesh);
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: MeshOps::compile rejects counts shorter than its index
 * list.
 */
inline void case_mesh_compile_face_counts_short() {
  static uint8_t src_buf[1024];
  static uint8_t dst_buf[1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh mesh;
  MeshState compiled;
  build_mismatched_polymesh(mesh, src_arena, 3, opaque<size_t>(4));
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: MeshOps::compile rejects counts longer than its index
 * list.
 */
inline void case_mesh_compile_face_counts_long() {
  static uint8_t src_buf[1024];
  static uint8_t dst_buf[1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh mesh;
  MeshState compiled;
  build_mismatched_polymesh(mesh, src_arena, 4, opaque<size_t>(3));
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: MeshOps::compile rejects a face whose span ends past the
 *        16-bit index range even though its start offset still fits.
 */
inline void case_mesh_compile_face_span_over_16bit() {
  static uint8_t src_buf[256 * 1024];
  static uint8_t dst_buf[256 * 1024];
  static uint8_t scratch_buf[1024];
  Arena src_arena(src_buf, sizeof(src_buf));
  Arena dst_arena(dst_buf, sizeof(dst_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  // The last triangle starts at 65535 and ends at 65538.
  constexpr size_t FACES = 21846;
  PolyMesh mesh;
  mesh.vertices.bind(src_arena, 3);
  for (int i = 0; i < 3; ++i)
    mesh.vertices.push_back(Vector{});
  mesh.face_counts.bind(src_arena, FACES);
  mesh.faces.bind(src_arena, FACES * 3);
  for (size_t f = 0; f < FACES; ++f) {
    mesh.face_counts.push_back(opaque<uint8_t>(3));
    for (uint16_t k = 0; k < 3; ++k)
      mesh.faces.push_back(opaque(k));
  }
  MeshState compiled;
  MeshOps::compile(mesh, compiled, dst_arena, scratch);
}

/**
 * @brief Death case: update_hankin rejects a retained topology that no
 *        classification of this pattern's output ever produced.
 * @details The topology array survives an angle re-solve on purpose, so a mesh
 *          pointed at a new pattern would otherwise carry class ids that no
 *          longer match the faces written into it.
 */
inline void case_update_hankin_stale_topology() {
  static uint8_t geom_buf[64 * 1024];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  CompiledHankin compiled;
  MeshOps::compile_hankin(cube, compiled, geom, scratch);

  MeshState mesh;
  mesh.topology.bind(geom, compiled.face_counts.size());
  for (size_t i = 0; i < compiled.face_counts.size(); ++i)
    mesh.topology.push_back(opaque<uint16_t>(0));
  MeshOps::update_hankin(compiled, mesh, geom, opaque(0.0f));
  if (mesh.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: update_hankin rejects retained topology from a pattern
 *        whose census matches but whose connectivity does not.
 * @details Cube- and octahedron-seeded patterns agree on face, index and vertex
 *          counts, so only the topology key separates them.
 */
inline void case_update_hankin_dual_seed_topology() {
  static uint8_t geom_buf[192 * 1024];
  static uint8_t scratch_buf[128 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);
  CompiledHankin cube_pattern;
  MeshOps::compile_hankin(cube, cube_pattern, geom, scratch);

  MeshState mesh;
  MeshOps::update_hankin(cube_pattern, mesh, geom, opaque(0.0f));
  MeshOps::classify_faces_by_topology(mesh, scratch, scratch, geom);

  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, geom);
  CompiledHankin octa_pattern;
  MeshOps::compile_hankin(octa, octa_pattern, geom, scratch);
  MeshOps::update_hankin(octa_pattern, mesh, geom, opaque(0.0f));
  if (mesh.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: CompiledHankin::clone rejects a self-aliased destination.
 * @details Each vector is rebound from the arena before the copy, so a
 *          self-clone memcpy's a block onto itself from a stale source pointer.
 */
inline void case_hankin_clone_aliases_dst() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  CompiledHankin compiled;
  CompiledHankin::clone(compiled, compiled, arena);
}

/**
 * @brief Death case: a face-offsets span with the wrong length must trap.
 * @details Mesh-borrow surface — the accessors index offsets by face, so an
 *          offsets array that is not one entry per face would read past its end
 *          on the solid scan path; set_borrowed rejects it at the install site.
 */
inline void case_mesh_state_set_borrowed_offsets_count_mismatch() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 1);
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 3);
  for (uint16_t i = 0; i < 3; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 2);
  offsets.push_back(opaque<uint16_t>(0));
  offsets.push_back(opaque<uint16_t>(3));
  MeshState m;
  // 2 offsets for 1 face -> HS_CHECK
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: face offsets that do not span the flat faces list must
 *        trap.
 * @details Mesh-borrow surface — the last offset plus that face's count must
 *          reach the end of the flat list, or a walk of the final face reads
 *          short of the data the view claims to cover.
 */
inline void case_mesh_state_set_borrowed_offsets_short_span() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 1);
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 4);
  for (uint16_t i = 0; i < 4; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 1);
  offsets.push_back(opaque<uint16_t>(0));
  MeshState m;
  // 0 + 3 != 4 -> HS_CHECK
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: face offsets that are not the counts' prefix sum must trap.
 * @details Mesh-borrow surface — the count and span checks pass on the endpoints
 *          alone, so an interior offset off the prefix sum would walk one face
 *          over another's indices. The audit walk catches it.
 */
inline void case_mesh_state_set_borrowed_offsets_not_prefix_sum() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  ArenaVector<uint8_t> counts(arena, 2);
  counts.push_back(opaque<uint8_t>(3));
  counts.push_back(opaque<uint8_t>(3));
  ArenaVector<uint16_t> faces(arena, 6);
  for (uint16_t i = 0; i < 6; ++i)
    faces.push_back(opaque(i));
  ArenaVector<uint16_t> offsets(arena, 2);
  offsets.push_back(opaque<uint16_t>(1)); // prefix sum starts at 0
  offsets.push_back(opaque<uint16_t>(3));
  MeshState m;
  m.set_borrowed(ArenaSpan<uint8_t>(counts), ArenaSpan<uint16_t>(faces),
                 ArenaSpan<uint16_t>(offsets));
  if (m.num_faces() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Builds a PolyMesh from an explicit face-count and flat index list.
 * @param mesh Mesh to populate.
 * @param arena Arena backing the mesh arrays.
 * @param num_verts Vertex count; positions are all zero (never read).
 * @param counts Per-face side counts.
 * @param num_faces Number of entries in @p counts.
 * @param indices Flat per-face vertex index list.
 * @param num_indices Number of entries in @p indices.
 */
inline void build_polymesh(PolyMesh &mesh, Arena &arena, size_t num_verts,
                           const uint8_t *counts, size_t num_faces,
                           const uint16_t *indices, size_t num_indices) {
  mesh.vertices.bind(arena, num_verts);
  for (size_t i = 0; i < num_verts; ++i)
    mesh.vertices.push_back(Vector{});
  mesh.face_counts.bind(arena, num_faces);
  for (size_t i = 0; i < num_faces; ++i)
    mesh.face_counts.push_back(opaque(counts[i]));
  mesh.faces.bind(arena, num_indices);
  for (size_t i = 0; i < num_indices; ++i)
    mesh.faces.push_back(opaque(indices[i]));
}

/**
 * @brief Death case: a zero-side face must trap while building half-edges.
 * @details Mesh-topology surface — a zero-count face emits no half-edges yet
 *          still claims a face slot, whose half_edge entry would then point at
 *          the next face's loop. The trailing triangle keeps the flat index
 *          list non-empty so the pairing scratch is a real allocation.
 */
inline void case_half_edge_zero_side_face() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {0, 3};
  const uint16_t indices[] = {0, 1, 2};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 3, counts, 2, indices, 3);
  HalfEdgeMesh half_edges(arena, mesh); // face 0 has zero sides -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: >2 half-edges on one undirected edge must trap.
 * @details Mesh-topology surface — three faces share edge (0,1), so the pairing
 *          pass sees a run of three where a 2-manifold allows at most two.
 *          Pairing the first two would leave the third silently unpaired.
 */
inline void case_half_edge_non_manifold_edge() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 1, 3, 0, 1, 4};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 5, counts, 3, indices, 9);
  HalfEdgeMesh half_edges(arena, mesh); // 3 half-edges on (0,1) -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: two faces wound the same way around a shared edge must
 *        trap.
 * @details Mesh-topology surface — both triangles traverse edge (0,1) in the
 *          same direction, so the undirected pairing key matches and they would
 *          otherwise pair into a mesh that passes require_closed_manifold while
 *          every vertex_orbit walk through the pair runs backwards.
 */
inline void case_half_edge_inconsistent_winding() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 1, 3};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 4, counts, 2, indices, 6);
  HalfEdgeMesh half_edges(arena, mesh); // both edges run 0->1 -> HS_CHECK
  if (half_edges.faces.size() == opaque<size_t>(99))
    std::printf("x");
}

/**
 * @brief Death case: a face side count past uint8_t must trap.
 * @details Mesh-topology surface — every operator narrows its output valence
 *          through this shared guard, so a high-valence orbit traps instead of
 *          wrapping the uint8_t face_counts entry.
 */
inline void case_mesh_narrow_face_count() {
  uint8_t c = MeshOps::narrow_face_count(opaque(UINT8_MAX + 1)); // -> HS_CHECK
  if (c == 0xEE)
    std::printf("x");
}

/**
 * @brief Death case: an open mesh must trap the closed-manifold requirement.
 * @details Mesh-topology surface — operators size their output pools from
 *          E = I/2, so a lone triangle's three unpaired half-edges are rejected
 *          up front instead of overrunning a pool far from the cause.
 */
inline void case_mesh_require_closed_manifold() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3};
  const uint16_t indices[] = {0, 1, 2};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 3, counts, 1, indices, 3);
  HalfEdgeMesh half_edges(arena, mesh);
  // unpaired -> HS_CHECK
  MeshOps::require_closed_manifold(half_edges, arena, "death");
}

/**
 * @brief Death case: a bowtie vertex must trap the closed-manifold requirement.
 * @details Mesh-topology surface — two tetrahedra joined at vertex 0 are closed
 *          and edge-manifold, so only the fan pass catches them; the orbit
 *          scaffolding would otherwise emit one face from the first fan and
 *          silently drop the second.
 */
inline void case_mesh_require_vertex_manifold() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  const uint8_t counts[] = {3, 3, 3, 3, 3, 3, 3, 3};
  const uint16_t indices[] = {0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2,
                              0, 4, 5, 0, 5, 6, 0, 6, 4, 4, 6, 5};
  PolyMesh mesh;
  build_polymesh(mesh, arena, 7, counts, 8, indices, 24);
  HalfEdgeMesh half_edges(arena, mesh);
  // split fan at vertex 0 -> HS_CHECK
  MeshOps::require_closed_manifold(half_edges, arena, "death");
}

/**
 * @brief Death case: a half-edge mesh whose faces have different side counts
 *        must trap even when the census matches.
 * @details Mesh-topology surface — the reuse overloads walk face loops from the
 *          half-edge mesh while sizing and indexing from the source mesh, so a
 *          {3,4} pairing against a {4,3} source shares (V,F,I) yet emits every
 *          face from the wrong span.
 */
inline void case_mesh_require_matching_face_sides() {
  static uint8_t buf[2048];
  Arena arena(buf, sizeof(buf));
  const uint8_t he_counts[] = {3, 4};
  const uint8_t mesh_counts[] = {4, 3};
  const uint16_t indices[] = {0, 1, 2, 3, 4, 5, 6};
  PolyMesh he_source;
  build_polymesh(he_source, arena, 7, he_counts, 2, indices, 7);
  HalfEdgeMesh half_edges(arena, he_source);
  PolyMesh mesh;
  build_polymesh(mesh, arena, 7, mesh_counts, 2, indices, 7);
  // face 1 starts at 3 in he_source, 4 in mesh -> HS_CHECK
  MeshOps::require_matching_half_edges(half_edges, mesh, "death");
}

/**
 * @brief Death case: a HANKIN step with no contact angle must trap.
 * @details Recipe-replay surface — the zero default collapses every star point
 *          onto its corner, so an authored step that forgot its angle replays
 *          as a flat tiling instead of failing.
 */
inline void case_apply_step_hankin_no_angle() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::HANKIN, opaque(0.0f)}};
  PolyMesh mesh = Solids::build_steps(opaque<uint8_t>(1), steps, 1, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a BEVEL step with no depth must trap.
 * @details Recipe-replay surface — the composite lowers to ambo, truncate(t),
 *          so a zero default is the depthless truncate the lowered replay
 *          already traps on; the authored replay must not diverge from it.
 */
inline void case_apply_step_bevel_no_depth() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::BEVEL, opaque(0.0f)}};
  const Solids::Recipe recipe = {opaque<uint8_t>(1), steps, 1};
  PolyMesh mesh = Solids::build_recipe(recipe, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
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
  Quaternion q =
      make_rotation(from, to); // NaN axis -> normalized() -> HS_CHECK
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
  Quaternion q =
      make_rotation(axis, nan); // NaN quat -> normalized() -> HS_CHECK
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
  Animation::NoiseParams p;
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
  Animation::Driver d(mutant, opaque<const float *>(nullptr),
                      1.0f); // -> HS_CHECK
  (void)d;
  if (mutant == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a second TransformerPool::init_storage() must trap.
 * @details Transformer surface — a re-init would hand the pool a second block
 *          while spawned animations still hold Params references into the first,
 *          and would silently double-charge the persistent arena.
 */
inline void case_transformer_pool_init_storage_twice() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.init_storage(persistent_arena); // entities already set -> HS_CHECK
}

/**
 * @brief Death case: spawning before init_storage() must trap.
 * @details Transformer surface — the slot scan indexes the entity block, so a
 *          spawn on an un-initialized pool would dereference null instead of
 *          reporting the missed init() wiring.
 */
inline void case_transformer_pool_spawn_before_init() {
  Timeline tl;
  RippleTransformer<2> rt(tl);
  Animation::Ripple *p = rt.spawn(0, Vector(0, 1, 0), 0.2f, 4); // -> HS_CHECK
  if (p == reinterpret_cast<Animation::Ripple *>(0x1))
    std::printf("x");
}

/**
 * @brief Death case: an out-of-range active index must trap.
 * @details Transformer surface — active_params() indexes the compact active list,
 *          which is shorter than CAPACITY, so an index taken from the slot domain
 *          (or from a stale count) would read a dead slot as if it were live.
 */
inline void case_transformer_pool_active_index_oob() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  const Animation::RippleParams &p = rt.active_params(opaque(0)); // -> HS_CHECK
  if (p.amplitude == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: reclaimed storage landing at a new address must trap.
 * @details Transformer surface — spawned animations hold Params references into
 *          the slots, so the post-reset replay must re-land the blocks exactly
 *          where init_storage() put them. Here the arena is NOT reset first, so
 *          the replay appends past the originals and every live reference would
 *          be left pointing at abandoned bytes.
 */
inline void case_transformer_pool_reclaim_storage_moved() {
  configure_arenas_default();
  Timeline tl;
  RippleTransformer<2> rt(tl);
  rt.init_storage(persistent_arena);
  rt.reclaim_storage(persistent_arena); // blocks land elsewhere -> HS_CHECK
}

/**
 * @brief Death case: a pool outliving its Timeline must trap.
 * @details Transformer surface — the destructor reaches back into the timeline
 *          to drop the pool's clear hook, and the spawned completion callbacks
 *          reach back into the pool, so the two lifetimes are ordered. An owner
 *          that declares them the other way gets a dead reference here rather
 *          than at some later step().
 */
inline void case_transformer_pool_outlives_timeline() {
  configure_arenas_default();
  alignas(Timeline) static uint8_t tl_storage[sizeof(Timeline)];
  Timeline *tl = new (tl_storage) Timeline();
  RippleTransformer<2> rt(*tl);
  rt.init_storage(persistent_arena);
  tl->~Timeline();
  // ~RippleTransformer at scope exit -> HS_CHECK
}

/**
 * @brief Concrete Effect for the canvas death cases.
 * @details Defaults to 32x16 and exposes register_param via reg.
 */
struct DeathEffect : public Effect {
  /**
   * @brief Constructs the effect at the requested resolution.
   * @param width Canvas width in pixels.
   * @param height Canvas height in pixels.
   */
  DeathEffect(int width = 32, int height = 16) : Effect(width, height) {}
  /**
   * @brief Draws one frame (no-op; the death cases never render).
   */
  void draw_frame() override {}
  /**
   * @brief Registers a parameter over the unit range, exposing register_param.
   * @param n Parameter name.
   * @param p Pointer to the backing float storage.
   */
  void reg(const char *n, float *p) { register_param(n, p, 0.0f, 1.0f); }
  /**
   * @brief Registers an integer parameter, exposing register_int_param.
   * @param n Parameter name.
   * @param p Pointer to the backing uint8_t storage.
   * @param min Minimum value, inclusive.
   * @param max Maximum value, inclusive.
   */
  void reg_int(const char *n, uint8_t *p, int min, int max) {
    register_int_param(n, p, min, max);
  }
};

/**
 * @brief Death case: a second simultaneously-live Effect must trap.
 * @details Canvas surface — the structural twin of case_timeline_double_construct.
 *          Every Effect aliases the same two static framebuffers and double-buffer
 *          indices, so a second live instance would scribble over the first's
 *          frames; the construction guard traps instead. The real app builds the
 *          next effect only after destroying the outgoing one.
 */
inline void case_effect_double_construct() {
  DeathEffect a;
  DeathEffect b; // second live Effect ctor -> HS_CHECK(!s_alive) -> trap
  if (opaque(a.strobe_columns()))
    std::printf("x");
}

/** @brief Death case: a zero Effect width must trap. */
inline void case_effect_width_zero() { DeathEffect fx(opaque(0), 16); }

/** @brief Death case: a zero Effect height must trap. */
inline void case_effect_height_zero() { DeathEffect fx(32, opaque(0)); }

/** @brief Death case: an Effect width above MAX_W must trap. */
inline void case_effect_width_over_max() {
  DeathEffect fx(opaque(MAX_W + 1), 16);
}

/** @brief Death case: an Effect height above MAX_H must trap. */
inline void case_effect_height_over_max() {
  DeathEffect fx(32, opaque(MAX_H + 1));
}

/**
 * @brief Initializes a particle system with the requested maximum lifetime.
 * @param max_life Maximum lifetime passed to ParticleSystem::init.
 */
inline void init_particle_system_with_lifetime(float max_life) {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, max_life);
}

/** @brief Death case: a zero particle lifetime must trap. */
inline void case_particle_lifetime_zero() {
  init_particle_system_with_lifetime(opaque(0.0f));
}

/** @brief Death case: a NaN particle lifetime must trap. */
inline void case_particle_lifetime_nan() {
  init_particle_system_with_lifetime(
      opaque(std::numeric_limits<float>::quiet_NaN()));
}

/** @brief Death case: a particle lifetime above uint16_t must trap. */
inline void case_particle_lifetime_over_max() {
  init_particle_system_with_lifetime(opaque(65536.0f));
}

/** @brief Pipeline stub for the particle-render lifetime death case. */
struct DeathPlotPipeline {
  void plot(Canvas &, const Vector &, const Pixel &, float, float) {}
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/** @brief Death case: a nonempty render with zero max_life must trap. */
inline void case_particle_render_zero_lifetime() {
  configure_arenas_default();
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, 1.0f);
  ps.spawn(Vector(1, 0, 0), Vector(), 0);
  ps.max_life = 0;

  DeathEffect fx;
  Canvas canvas(fx);
  DeathPlotPipeline pipeline;
  Plot::ParticleSystem::draw<32, 16>(pipeline, canvas, ps,
                                     [](const Vector &, Fragment &) {});
}

/**
 * @brief Death case: a second simultaneously-live correction guard must trap.
 * @details LED surface — NoColorCorrection and NoTempCorrection share one
 *          liveness flag and set the global FastLED correction/temperature, so a second
 *          live guard of either type would leave the wrong baseline on the earlier
 *          guard's exit; the construction guard traps instead.
 */
inline void case_correction_guard_double_construct() {
  NoColorCorrection a;
  NoColorCorrection b; // second live guard -> liveness HS_CHECK -> trap
  if (correction_guard_live() == opaque(true))
    std::printf("x");
}

/**
 * @brief Death case: a live NoColorCorrection plus a NoTempCorrection must trap.
 * @details LED surface — the two guard types share the one liveness flag, so a
 *          second live guard of the OTHER type is as unsafe as a same-type
 *          double-construct; this is the case the shared "either type" contract
 *          exists to guarantee. The construction guard traps on either.
 */
inline void case_correction_guard_cross_type() {
  NoColorCorrection a;
  NoTempCorrection b; // second live guard of a different type -> trap
  if (correction_guard_live() == opaque(true))
    std::printf("x");
}

/**
 * @brief Death case: overflowing the fixed ParamList must trap.
 * @details Canvas surface — register_param traps rather than silently dropping a
 *          registration, which would desync the GUI and, on WASM, break the
 *          no-realloc memory-view invariant.
 */
inline void case_register_param_overflow() {
  DeathEffect fx;
  static float slot = 0.0f;
  // Distinct names, so the capacity guard fires ahead of the duplicate guard.
  constexpr int CAPACITY = static_cast<int>(Effect::ParamList::FIXED_CAPACITY);
  static char names[CAPACITY + 1][8];
  for (int i = 0; i < opaque(CAPACITY + 1); ++i) {
    std::snprintf(names[i], sizeof(names[i]), "p%d", i);
    fx.reg(names[i], &slot);
  }
}

/**
 * @brief Death case: an integer param bound the target cannot store must trap.
 * @details Canvas surface — a value write narrows through
 *          static_cast<Integer>(float), which is undefined once the registered
 *          range leaves the storage type.
 */
inline void case_register_int_param_range() {
  DeathEffect fx;
  static uint8_t slot = 0;
  fx.reg_int("count", &slot, 0, opaque(256));
}

/**
 * @brief Death case: set_clip rejects x_end beyond the canvas width.
 */
inline void case_set_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  fx.set_clip(0, H, 0, opaque(W + 1));
}

/** @brief Death case: clip state cannot change after a frame begins. */
inline void case_set_clip_mid_frame() {
  DeathEffect fx;
  Canvas canvas(fx);
  fx.set_clip(0, fx.height(), 0, fx.width());
}

/** @brief Death case: the publication envelope rejects values above one. */
inline void case_output_envelope_out_of_range() {
  DeathEffect fx;
  fx.set_output_envelope(opaque(1.01f));
}

/**
 * @brief Death case: set_clip_x rejects x_end beyond the canvas width.
 */
inline void case_set_clip_x_out_of_bounds() {
  constexpr int W = 32;
  DeathEffect fx;
  fx.set_clip_x(0, opaque(W + 1));
}

/**
 * @brief Death case: an arc start outside [0, w) must trap.
 * @details Clip surface — arcs_overlap wraps the seam-relative offset with one
 *          conditional add instead of a modulo, which only lands in range while
 *          both starts are already reduced; a start outside the cylinder would
 *          silently report the wrong overlap. Both lengths are positive and
 *          under w so the early-out branches do not preempt the guard.
 */
inline void case_arcs_overlap_start_out_of_range() {
  bool hit = ClipRegion::arcs_overlap(opaque(-1), opaque(2), opaque(0),
                                      opaque(2), opaque(8));
  if (hit)
    std::printf("x");
}

/**
 * @brief Death case: the shader rejects a clip beyond its LUT width.
 * @details Direct construction exercises the downstream guard independently
 *          of the Effect setters.
 */
inline void case_scan_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = opaque(W + 1);
  cr.y_end = H;
  cr.margin = 0;
  cr.w = W;
  cr.h = H;
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: the shader rejects a render band past the canvas rows.
 * @details The rows the guard admits subscript the canvas, so a band reaching
 *          past H must be rejected even though the phi LUT holds those rows.
 */
inline void case_scan_clip_rows_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = W;
  cr.y_end = opaque(H + 1);
  cr.margin = 0;
  cr.w = W;
  cr.h = opaque(H + 1);
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later Face claimed.
 * @details The second build retargets the first face's spans, so the first no
 *          longer describes its own geometry.
 */
inline void case_face_scratch_retargeted() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  Basis basis = make_basis(Quaternion(), Vector(0, 1, 0));
  Vector verts[6];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 6; ++i) {
    float a = (2.0f * PI_F * i) / 6.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
  }
  for (int i = 0; i < 3; ++i) {
    idx_a[i] = static_cast<uint16_t>(i);
    idx_b[i] = static_cast<uint16_t>(i + 3);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face first(std::span<const Vector>(verts, 6),
                  std::span<const uint16_t>(idx_a, 3), scratch, HV, H);
  SDF::Face second(std::span<const Vector>(verts, 6),
                   std::span<const uint16_t>(idx_b, 3), scratch, HV, H);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/**
 * @brief Death case: a scan rejects a canvas that is not its <W, H>.
 * @details Direct construction exercises the guard independently of the draw
 *          primitives that call it.
 */
inline void case_scan_canvas_dim_mismatch() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, opaque(H + 1));
  Canvas c(fx);
  Scan::check_canvas_dims<W, H>(c);
}

/**
 * @brief Death case: a scan rejects a direct-raster sink prepared for no
 *        canvas.
 * @details The sink writes through its cached framebuffer base, which the
 *          double-buffered canvas leaves pointing at the buffer being scanned
 *          out until prepare() runs.
 */
inline void case_scan_pipeline_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  Scan::check_pipeline_prepared(sink, c);
}

/**
 * @brief Death case: erasing a direct-raster sink prepared for no canvas.
 * @details PipelineRef drops prepared_for(), so a draw taking the erased handle
 *          cannot check the cached base; the erasure checks it instead.
 */
inline void case_pipeline_ref_erase_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  PipelineRef erased(sink, c);
  (void)erased;
}

/**
 * @brief Death case: a face index past the vertex pool must trap.
 * @details SDF::Face reads the vertex span with operator[], which only asserts,
 *          so a stale index domain would read arbitrary memory as a Vector on
 *          device. Scan::Mesh validates the flat index array once per mesh.
 */
inline void case_scan_mesh_face_index_out_of_range() {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[1024];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  mesh.vertices.bind(geom, 3);
  for (int i = 0; i < 3; ++i)
    mesh.vertices.push_back(Vector(0.0f, 1.0f, 0.0f));
  mesh.face_counts.bind(geom, 1);
  mesh.face_counts.push_back(opaque<uint8_t>(3));
  mesh.face_offsets.bind(geom, 1);
  mesh.face_offsets.push_back(opaque<uint16_t>(0));
  mesh.faces.bind(geom, 3);
  mesh.faces.push_back(opaque<uint16_t>(0));
  mesh.faces.push_back(opaque<uint16_t>(1));
  mesh.faces.push_back(opaque<uint16_t>(3)); // only 3 vertices -> HS_CHECK

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Scan::Mesh::draw<W, H>(
      pipe, c, mesh, [](const Vector &, Fragment &) {}, scratch);
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
  Plot::Mesh::draw<W, H>(
      pipe, c, mesh, [](const Vector &, Fragment &) {}); // index 130 -> trap
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
 *          authoring/config error the project routes to HS_CHECK (enabled
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
  fb.flush(c, 1.0f);
}

/**
 * @brief Death case: retuning a screen trail to a non-positive lifetime must trap.
 * @details Filter surface — Screen::Trails::set_lifetime carries the
 *          constructor's bound, so a slider that reaches zero traps here rather
 *          than dividing by it in flush()'s fade progress.
 */
inline void case_screen_trails_set_lifetime_nonpositive() {
  Filter::Screen::Trails<8> trails(4);
  trails.set_lifetime(opaque(0));
}

/**
 * @brief Death case: a ring index past the last ring must trap.
 * @details SphericalFieldLayout surface — the chain walk saturates at row H-1
 *          while the offset keeps accumulating, so an out-of-range index would
 *          otherwise hand back a Ring pointing past the sample array.
 */
inline void case_spherical_field_ring_index_oob() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4);
  const auto ring = layout.ring(opaque(layout.ring_count()));
  if (ring.offset == 42)
    std::printf("x");
}

/**
 * @brief Death case: populating past the last ring must trap.
 * @details next_ring() saturates at the last ring, so an overrunning band would
 *          re-populate it and leave the caller believing it wrote fresh rings.
 */
inline void case_spherical_field_populate_ring_end_oob() {
  constexpr hs::SphericalFieldLayout<32, 16, 0> layout(4);
  static float values[layout.sample_count()];
  hs::SphericalField<float, 32, 16, 0> field(values, layout);
  field.populate(0, opaque(layout.ring_count()),
                 [](const Vector &v, const auto &) { return v.y; });
  if (values[0] == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a harmonic mode past the float factorial range must trap.
 * @details (l + |m|)! overflows to infinity, so the factorial ratio collapses
 *          to 0 and every sample of the mode comes back black.
 */
inline void case_spherical_harmonic_normalization_overflow() {
  const float n = SHMath::normalization(opaque(20), 20);
  if (n == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an order past the degree must trap.
 * @details reduced_legendre() has no term to recur on for |m| > l and returns
 *          0, so every sample of the mode comes back black.
 */
inline void case_spherical_harmonic_order_over_degree() {
  const float n = SHMath::normalization(opaque(2), 3);
  if (n == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a negative flat harmonic index must trap.
 * @details sqrtf of a negative argument is NaN, and the cast of a NaN to int
 *          is undefined, so the decoded level would be arbitrary.
 */
inline void case_spherical_harmonic_decode_negative_index() {
  auto [l, m] = SHMath::decode_lm(opaque(-1));
  if (l == 42 && m == 42)
    std::printf("x");
}

/**
 * @brief Death case: an infill band past the rendered domain must trap.
 * @details A south_infill wider than H puts every row at full longitude
 *          resolution, multiplying sample_count() by the spacing; the arena
 *          would then overflow at an unrelated call site.
 */
inline void case_spherical_field_infill_over_domain() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, opaque(17));
  if (layout.sample_count() == 42)
    std::printf("x");
}

/**
 * @brief Death case: a negative feedback fade must trap in sync_hue.
 * @details Style is a public aggregate, so nothing but a slider bound keeps fade
 *          non-negative. logf of a negative yields NaN, the hue matrix carries it
 *          into every feedback pixel, and float_to_pixel16 clamps NaN to 65535 —
 *          a white buffer with no other symptom. The guard also catches a NaN
 *          fade, which compares false against zero.
 */
inline void case_feedback_negative_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(-0.25f); // logf(negative) -> HS_CHECK
  style.sync_hue();
}

/**
 * @brief Death case: an infinite feedback fade must trap in sync_hue.
 * @details A comparison against zero admits +INFINITY, whose -logf is -inf; the
 *          turn-to-cos/sin reduction of an infinite angle is NaN, which spreads
 *          through all nine hue matrix entries and whitens every feedback pixel.
 */
inline void case_feedback_infinite_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(std::numeric_limits<float>::infinity());
  style.sync_hue();
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
 * @brief Death case: AlphaFalloffShade requires a non-null callback.
 */
inline void case_alpha_falloff_null() {
  auto fn = opaque<AlphaFalloffShade::FalloffFunction>(nullptr);
  AlphaFalloffShade shade(fn);
  (void)shade;
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null palette source.
 */
inline void case_noise_hue_palette_null_source() {
  static Pixel hue_rotation_lut[1];
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), hue_rotation_lut,
               hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-rotation LUT.
 */
inline void case_noise_hue_palette_null_rotation_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const Pixel *>(nullptr), hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-noise LUT.
 */
inline void case_noise_hue_palette_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static Pixel hue_rotation_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, hue_rotation_lut, opaque<const int8_t *>(nullptr));
}

/**
 * @brief Death case: rebaking an endpoint-aliasing blend result must trap.
 * @details Color surface — bake_palette_blend's w <= 0 fast path hands @c dst
 *          the @c from endpoint's LUT storage rather than baking a copy, so a
 *          rebake through @c dst would silently rewrite the endpoint every
 *          other consumer samples.
 */
inline void case_baked_palette_rebake_aliased() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPalette from, to, dst;
  from.bake(a, src);
  to.bake(a, src);
  bake_palette_blend(dst, a, from, to, opaque(0.0f)); // dst aliases from
  dst.rebake(src);                                    // -> HS_CHECK
}

/**
 * @brief Death case: cloning a BakedPalette from itself must trap.
 * @details Color surface — clone_from allocates fresh storage into this handle
 *          before reading @c src, so a self-clone memcpys uninitialized arena
 *          onto itself and leaves the LUT filled with garbage.
 */
inline void case_baked_palette_clone_from_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPalette lut;
  lut.bake(a, src);
  BakedPalette &self = *opaque(&lut);
  lut.clone_from(self, a); // -> HS_CHECK
}

/**
 * @brief Death case: blending a BakedPalette with itself as an endpoint must
 *        trap.
 * @details Color surface — bake_blend reallocates this handle before walking
 *          the endpoints, so an endpoint that is the output reads the fresh
 *          uninitialized arena instead of the baked LUT.
 */
inline void case_baked_palette_bake_blend_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPalette from, dst;
  from.bake(a, src);
  dst.bake(a, src);
  BakedPalette &self = *opaque(&dst);
  dst.bake_blend(a, from, self, opaque(0.5f)); // -> HS_CHECK
}

/**
 * @brief Death case: a Gradient built from an empty stop list must trap.
 * @details Color surface — with no stops the 256-entry LUT keeps its
 *          value-initialized state, so every lookup returns black. The
 *          constructor requires at least one stop rather than yielding a
 *          silently all-black palette.
 */
inline void case_gradient_no_stops() {
  Gradient grad({}); // empty stop list -> HS_CHECK
  (void)grad;
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
 * @brief A wedged transport: never completes, and traps on the watchdog consult.
 * @details Stands in for a TeensySPIDMA whose completion ISR never fires. Its
 *          checkStaleTransfer() plays the role the real driver's watchdog does on
 *          the overrun-drop path — trapping a permanently in-flight channel.
 */
struct WedgedStrip {
  /**
   * @brief Matches the transport ctor contract; the clock is unused.
   */
  explicit WedgedStrip(uint32_t) {}
  /**
   * @brief No-op init.
   */
  void init() {}
  /**
   * @brief Always reports the channel busy, so every submit hits the overrun path.
   * @return Always false.
   */
  bool isComplete() const { return false; }
  /**
   * @brief Traps: a wedged channel surfaced from the overrun-drop path.
   */
  void checkStaleTransfer() { HS_CHECK(false, "DMA channel wedged"); }
  /**
   * @brief Unreachable here (the overrun path returns before transmitting).
   */
  void transmitAsync(const uint8_t *, size_t) {}
};

/**
 * @brief Death case: submitFrame() consults the watchdog on overrun, so a wedged
 *        channel traps rather than dropping frames forever.
 * @details Controller surface — pins the ordering that the overrun-drop branch
 *          calls checkStaleTransfer() before bumping the counter and returning
 *          false. The transfer-stale predicate itself is covered in-process by
 *          test_dma_core.h; this proves the controller wires the drop path into
 *          the watchdog so a permanently in-flight transport fails fast.
 */
inline void case_dma_controller_wedged_overrun() {
  static DMALEDController<8, WedgedStrip> ctl;
  bool ok =
      ctl.submitFrame(opaque(false)); // busy -> checkStaleTransfer -> trap
  if (ok)
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
 *          by tests/test_concepts.h. Host/WASM only: the device Fn backend
 *          returns a zero-initialized R instead of trapping (row 9 of
 *          docs/ledgers/device_host_divergence_ledger.md).
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
 * @brief Death case: registering two effects under one name must trap.
 * @details Registry surface — the name keys the factory lookup and the
 *          HS_EFFECT_LIST anti-drift oracle, so a duplicate (an effect header
 *          pulled into two translation units) would leave a shadowed entry that
 *          can never be selected. The append-time guard traps at static-init.
 */
inline void case_effect_registry_duplicate_name() {
  EffectRegistration reg{};
  reg.name = "DeathDuplicate";
  EffectRegistry::add(reg);
  EffectRegistry::add(reg); // name already present -> HS_CHECK
}

/**
 * @brief Death case: a Flywheel period of zero must trap at construction.
 * @details POV-sync surface — position() divides the int32 elapsed window by the
 *          period, so a zero divides by zero and an over-large one voids the
 *          signed-safe coast window; the constructor rejects both before the
 *          driver ever schedules a column.
 */
inline void case_flywheel_period_zero() {
  pov::sync::Config cfg;
  cfg.cycles_per_half_rev = opaque<uint32_t>(0);
  pov::sync::Flywheel fw(cfg); // period 0 -> HS_CHECK
  (void)fw;
}

/**
 * @brief Death case: scheduling a boundary burst onto a busy wire must trap.
 * @details Sync surface — the emitter drives one shared wire, so a second
 *          burst scheduled while pulses are still due would splice two symbols
 *          into one count and decode as a third.
 */
inline void case_symbol_emitter_overlapping_burst() {
  const pov::sync::Config cfg = pov::sync::phantasm_config(
      opaque<uint32_t>(600000000u), opaque<uint32_t>(480u), opaque(288),
      opaque(4));
  pov::sync::SymbolEmitter em;
  const bool first = em.schedule_boundary(
      pov::sync::Symbol::ZERO, opaque<uint32_t>(0), opaque<uint32_t>(0), cfg);
  const bool second = // pulses still due -> HS_CHECK
      em.schedule_boundary(pov::sync::Symbol::HALF, opaque<uint32_t>(0),
                           opaque<uint32_t>(0), cfg);
  if (first && second)
    std::printf("x");
}

/**
 * @brief Death case: a virtual height of one row must trap in the phi mapping.
 * @details Geometry surface — the row-to-angle scale divides by (h_virt - 1),
 *          so a single-row canvas would map every row to a non-finite phi.
 */
inline void case_y_to_phi_degenerate_height() {
  float phi = y_to_phi_virtual(opaque(0.0f), opaque(1)); // divisor 0 -> trap
  if (phi == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: reading an orientation frame past the history must trap.
 * @details Geometry surface — the motion-blur history is a fixed array whose
 *          live prefix is num_frames long, so an index past it would read a
 *          stale or never-written quaternion instead of failing.
 */
inline void case_orientation_frame_index_oob() {
  Orientation<> orientation; // constructed with one frame
  const Quaternion &q = orientation.get(opaque(3));
  if (q.r == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a non-unit quaternion must trap.
 * @details Geometry surface — the rotation assumes a unit quaternion, so a
 *          finite but over-long one would scale and shear the frame rather than
 *          rotate it; the guard fires before the axes are built.
 */
inline void case_make_basis_nonunit_quaternion() {
  Quaternion q(opaque(2.0f), opaque(0.0f), opaque(0.0f), opaque(0.0f));
  Vector normal{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  Basis b = make_basis(q, normal); // |q| = 2 -> HS_CHECK
  if (b.u.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a polygon with fewer than three sides must trap.
 * @details SDF surface — the sector fold divides a full turn by the side count,
 *          so a 2-gon has no interior for the distance to be measured against.
 */
inline void case_sdf_polygon_side_count() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  SDF::PlanarPolygon poly(b, opaque(0.5f), opaque(2), opaque(0.0f));
  if (poly.apothem == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: an angular repeat around a non-unit axis must trap.
 * @details SDF surface — the sector fold rotates the query point about the
 *          axis, so a non-unit one scales every folded copy off the sphere.
 */
inline void case_sdf_angular_repeat_nonunit_axis() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(1.0f), opaque(0.1f));
  Vector axis{opaque(0.0f), opaque(2.0f), opaque(0.0f)};
  SDF::AngularRepeat<SDF::Ring> rep(ring, opaque(4), axis); // non-unit -> trap
  if (rep.sector == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a knot ring with no cells must trap.
 * @details SDF surface — the per-pixel cell index divides the azimuth by
 *          2π/n, so n == 0 wraps to knots[-1] on every probe.
 */
inline void case_sdf_distorted_ring_zero_knots() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  const float knots[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(0),
                          opaque(0.0f), pf); // no knot cells -> trap
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a twist warp around a zero-radius torus must trap.
 * @details SDF warp surface — the Lipschitz bound scales by 2/R, so a zero
 *          major radius hands the rasterizer a non-finite step bound.
 */
inline void case_sdf_twist_zero_major_radius() {
  SDF::Warp::Twist tw(opaque(2), opaque(0.1f), opaque(0.0f)); // R = 0 -> trap
  if (tw.two_over_r == opaque(42.0f))
    std::printf("x");
}

/** @brief Draw callback for the OpLeg construction death cases; never runs. */
inline void death_opleg_draw(Canvas &, MeshState &,
                             const Animation::OpLeg::Shading &) {}

/**
 * @brief A palette handoff complete enough to clear the OpLeg handoff guard.
 * @return A handoff naming a default bank and a one-face departed palette.
 * @details Lets a case reach a guard that sits behind the handoff check; the
 *          LUTs are never sampled, since every such case traps in the
 *          constructor before the first frame.
 */
inline Animation::OpLeg::PaletteHandoff death_opleg_handoff() {
  static const BakedPaletteBank bank;
  static const uint8_t face_palette[1] = {0};
  return {.bank = &bank,
          .prev_face_palette = face_palette,
          .prev_faces = 1,
          .prev_face_centroid = nullptr,
          .correspondence = Animation::OpLeg::FaceCorrespondence::GEOMETRIC};
}

/**
 * @brief Death case: choosing an edge from a node outside the graph must trap.
 * @details ConwayGraph surface — no EDGES row touches such a node, so the
 *          weighted pick would have nothing to divide by.
 */
inline void case_pick_next_edge_unknown_node() {
  uint8_t visits[ConwayGraph::NUM_NODES] = {};
  const int e = ConwayGraph::pick_next_edge(opaque<int>(ConwayGraph::NUM_NODES),
                                            -1, 0, visits, 0u);
  if (e == opaque<int>(-42))
    std::printf("x");
}

/** @brief A well-formed non-settling graph edge for the OpLeg death cases. */
inline constexpr ConwayGraph::EdgeSpec death_opleg_edge{
    .from_node = 0,
    .to_node = 0,
    .seed_solid = 0,
    .op = ConwayGraph::MorphOp::TRUNCATE,
    .t_from = 0.0f,
    .t_to = 0.4f,
    .twist_from = 0.0f,
    .twist_to = 0.0f,
    .settle = false,
    .reseed = ConwayGraph::Reseed::NONE,
    .bridge = false};

/**
 * @brief Death case: an edge-sweep leg without a graph edge must trap.
 * @details OpLeg surface — the constructor reads the edge's operator and settle
 *          flag on its first line, so a null edge is a null dereference.
 */
inline void case_opleg_edge_sweep_no_edge() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed, Animation::OpLeg::EdgeSweepSpec{}, arena,
                       death_opleg_draw, handoff); // null edge -> HS_CHECK
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a leg with a non-positive sweep length must trap.
 * @details OpLeg surface — the per-frame sweep parameter divides by the frame
 *          count, and a zero-frame leg would also complete before drawing.
 */
inline void case_opleg_zero_sweep_frames() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                       .sweep_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a leg built without a palette handoff must trap.
 * @details OpLeg surface — the departed node's per-face palette keys every
 *          blend ramp the leg bakes, so an absent bank leaves the whole
 *          crossfade unresolvable rather than merely uncolored.
 */
inline void case_opleg_incomplete_palette_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff; // no bank, no per-face palette
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                       .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: settle frames that contradict the edge must trap.
 * @details OpLeg surface — the edge's settle flag decides whether the leg
 *          computes a relaxed endpoint at all, so a settle window on a
 *          non-settling edge would slerp toward vertices nothing produced.
 */
inline void case_opleg_edge_settle_mismatch() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::EdgeSweepSpec{
                           .edge = &death_opleg_edge,
                           .reverse = false,
                           .sweep_frames = opaque(1),
                           .settle_frames = opaque(1)}, // edge does not settle
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: an edge-sweep leg without a palette handoff must trap.
 * @details OpLeg surface — the edge-sweep constructor carries its own handoff
 *          guard, distinct from the param-sweep one.
 */
inline void case_opleg_edge_sweep_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::EdgeSweepSpec{.edge = &death_opleg_edge,
                                      .reverse = false,
                                      .sweep_frames = opaque(1),
                                      .settle_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a hankin leg without a palette handoff must trap.
 * @details OpLeg surface — the hankin constructor carries its own handoff
 *          guard.
 */
inline void case_opleg_hankin_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.1f),
                                        .theta_end = opaque(0.5f),
                                        .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a hankin leg sweeping to a smaller angle must trap.
 * @details OpLeg surface — the leg sweeps the slerp fraction outward from the
 *          collapsed corner, which is monotone only while the arrival angle is
 *          the larger of the two.
 */
inline void case_opleg_hankin_backward_theta() {
  static uint8_t buf[8192];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  const Animation::OpLeg::PaletteHandoff handoff = death_opleg_handoff();
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.9f),
                                        .theta_end = opaque(0.3f),
                                        .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a relax leg with neither a bake nor iterations must trap.
 * @details OpLeg surface — the leg needs a relaxed endpoint to slerp to, which
 *          is either the shipped bake or the result of live iterations; with
 *          neither it would slerp the seed onto itself for its whole duration.
 */
inline void case_opleg_relax_no_iterations() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::RelaxSpec{.iterations = opaque(0),
                                                   .bake = nullptr,
                                                   .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a relax leg without a palette handoff must trap.
 * @details OpLeg surface — the relax constructor carries its own handoff guard.
 */
inline void case_opleg_relax_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::RelaxSpec{.iterations = opaque(1),
                                                   .bake = nullptr,
                                                   .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a medial leg without a palette handoff must trap.
 * @details OpLeg surface — the medial constructor carries its own handoff
 *          guard.
 */
inline void case_opleg_medial_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::MedialSpec{.sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a reconcile leg without endpoints must trap.
 * @details OpLeg surface — the leg slerps every seed vertex to an authored
 *          position, so an absent endpoint array is the whole leg's target.
 */
inline void case_opleg_reconcile_no_endpoints() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::ReconcileSpec{
                           .to_positions = nullptr, .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a reconcile leg without a palette handoff must trap.
 * @details OpLeg surface — the reconcile constructor carries its own handoff
 *          guard.
 */
inline void case_opleg_reconcile_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  static const Vector endpoints[1] = {Vector(0, 0, 1)};
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ReconcileSpec{.to_positions = endpoints,
                                      .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a gated-swap leg with no gate window must trap.
 * @details OpLeg surface — the leg runs 2*gate_frames + 1 frames around the
 *          swap, so a zero gate leaves the swap frame with no approach or
 *          departure to blend across.
 */
inline void case_opleg_gated_swap_zero_gate_frames() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a gated-swap leg without a palette handoff must trap.
 * @details OpLeg surface — the gated-swap constructor carries its own handoff
 *          guard.
 */
inline void case_opleg_gated_swap_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a shading lookup past the leg's face table must trap.
 * @details OpLeg surface — the ramp index is read straight from the per-face
 *          table, so an out-of-range face would shade through whatever follows
 *          it instead of failing.
 */
inline void case_opleg_shading_face_out_of_range() {
  static BakedPalette ramps[1];
  static const uint8_t face_ramp[1] = {0};
  const Animation::OpLeg::Shading shading{
      .ramps = ramps, .face_ramp = face_ramp, .faces = opaque<size_t>(1)};
  if (&shading.ramp_for(opaque<size_t>(3)) == ramps)
    std::printf("x");
}

/**
 * @brief Death case: a Motion over an empty path must trap on the first step.
 * @details Animation surface — an unfilled Path samples the origin, and the
 *          angle between two zero vectors is a NaN that reaches the Orientation
 *          as a NaN quaternion rather than a visible fault.
 */
inline void case_motion_empty_path_origin_sample() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  Orientation<4> orientation;
  static Path<32> path; // never appended -> get_point returns the origin
  Animation::Motion<W, 4> motion(orientation, path, opaque(10));
  motion.step(c);
}

/**
 * @brief Death case: a negative equator sample count must trap.
 * @details Spherical-field surface — the count sizes every ring's longitude
 *          walk, so a negative one underflows the per-ring sample allocation.
 */
inline void case_spherical_field_negative_equator_samples() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, 0, opaque(-1));
  if (layout.sample_count() == 42)
    std::printf("x");
}

/**
 * @brief Death case: a spherical polygon wider than a hemisphere must trap.
 * @details SDF surface — beyond the hemisphere the cap fold changes sign, so
 *          the shape must be built inverted about its antipode instead.
 */
inline void case_sdf_spherical_polygon_radius_over_hemisphere() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  SDF::SphericalPolygon poly(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (poly.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a flower wider than a hemisphere must trap.
 * @details SDF surface — the petal cap bound is taken about the antipode, so a
 *          radius past the hemisphere inverts the band it derives.
 */
inline void case_sdf_flower_radius_over_hemisphere() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a zero-radius flower must trap.
 * @details SDF surface — the petal parameter divides by the circumradius, so a
 *          zero radius hands every probe a non-finite distance.
 */
inline void case_sdf_flower_zero_radius() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(0.0f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT for a degenerate polygon must trap.
 * @details SDF class-LUT surface — fewer than three vertices leaves no closed
 *          boundary for the crossing test, so every sample would read as
 *          outside.
 */
inline void case_sdf_class_lut_too_few_vertices() {
  static const float poly_xy[4] = {-0.5f, -0.5f, 0.5f, -0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(2), opaque(4), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT on a single-cell grid must trap.
 * @details SDF class-LUT surface — the cell step divides by (n - 1), so a
 *          resolution below 2 makes the whole domain non-finite.
 */
inline void case_sdf_class_lut_grid_too_small() {
  static const float poly_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(3), opaque(1), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: binding a class LUT at a vertex offset outside the face
 *        must trap.
 * @details SDF class-LUT surface — the offset indexes the canonical polygon
 *          cyclically, so an out-of-range one correlates the face against
 *          storage past the shape.
 */
inline void case_sdf_bind_class_lut_offset_out_of_range() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  Basis basis = make_basis(Quaternion(), Vector(0, 1, 0));
  Vector verts[3];
  uint16_t idx[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts, 3),
                 std::span<const uint16_t>(idx, 3), scratch, HV, H);
  static const float canon_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  SDF::ClassLut lut;
  if (face.bind_class_lut(&lut, canon_xy, opaque(7), false))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring built with a null shift callback must
 *        trap.
 * @details SDF ring surface — the callback is invoked per azimuth on every
 *          probe, so a null one faults deep inside the rasterizer instead.
 */
inline void case_sdf_distorted_ring_null_shift() {
  const Basis b{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
  ScalarFn shift; // default-constructed -> empty
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), shift, opaque(0.1f),
                          opaque(0.0f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: fixing a ShaderWorkbench preset view to nothing must trap.
 * @details Effect surface — an empty view leaves preset_count_for_view() at
 *          zero, so preset selection would index an empty roster.
 */
inline void case_shader_workbench_empty_preset_view() {
  using WB = shader_workbench_tests::ShaderWorkbenchWhiteBox;
  WB::SB sb;
  WB::set_fixed_preset_view(sb, std::span<const uint8_t>());
}

/**
 * @brief Death case: a preset view naming a preset past the roster must trap.
 * @details Effect surface — the view indirects into PRESETS, so an out-of-range
 *          entry reads a config off the end of the table.
 */
inline void case_shader_workbench_preset_view_index_out_of_range() {
  using WB = shader_workbench_tests::ShaderWorkbenchWhiteBox;
  static const uint8_t indices[1] = {200};
  WB::SB sb;
  WB::set_fixed_preset_view(sb, std::span<const uint8_t>(indices, 1));
}

/**
 * @brief Death case: resolving a preset past the view must trap.
 * @details Effect surface — the lookup indexes PRESETS through the view, so an
 *          out-of-range index hands the pipeline a config read off the table.
 */
inline void case_shader_workbench_preset_for_view_out_of_range() {
  using WB = shader_workbench_tests::ShaderWorkbenchWhiteBox;
  WB::SB sb;
  if (reinterpret_cast<uintptr_t>(
          &WB::preset_for_view(sb, opaque<size_t>(200))) == 0x1)
    std::printf("x");
}

/**
 * @brief Death case: an operator table whose entry decreases carrier family
 *        rank must trap.
 * @details Interpreter surface — compile() only matches adjacent carriers, so
 *          a rank-decreasing entry would run a chain the type system forbids.
 */
inline void case_chain_table_rank_decreases() {
  static Pullback::Interp::ChainProgram program;
  static Pullback::Interp::OperatorDescriptor descriptor{};
  descriptor.input = Pullback::Interp::CarrierId::COLOR;
  descriptor.output = Pullback::Interp::CarrierId::SPHERE;
  alignas(std::max_align_t) static uint8_t block_a[64];
  alignas(std::max_align_t) static uint8_t block_b[64];
  program.bind_storage(
      block_a, block_b, opaque<size_t>(sizeof(block_a)),
      std::span<const Pullback::Interp::OperatorDescriptor>(&descriptor, 1));
}

/**
 * @brief A named death case selected by HS_DEATH_CASE in the child process.
 */
struct Case {
  const char *name;       /**< Case selector matched against HS_DEATH_CASE. */
  void (*fn)();           /**< The trap-triggering case body to run. */
  const char *guard_file; /**< Basename check_fail() logs for the guard that
                               must fire; nothing else counts as a pass. */
  const char *guard_text; /**< Expected "(condition) message" tail of that
                               guard's breadcrumb line. */
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
      {"arena_oom", case_arena_oom, "memory.h",
       "(false) Arena::allocate: out of memory"},
      {"arena_make_oom", case_arena_make_oom, "memory.h",
       "(false) Arena::allocate: out of memory"},
      {"arena_zero_size_alloc", case_arena_zero_size_alloc, "memory.h",
       "(size > 0) Arena::allocate: zero-size request"},
      {"arena_bad_alignment", case_arena_bad_alignment, "memory.h",
       "(align != 0 && (align & (align - 1)) == 0) Arena::allocate: "
       "alignment "},
      {"arena_set_capacity_below_offset", case_arena_set_capacity_below_offset,
       "memory.h",
       "(offset <= new_capacity) Arena::set_capacity below the live offset "
       "would strand content"},
      {"arena_set_capacity_above_extent", case_arena_set_capacity_above_extent,
       "memory.h",
       "(new_capacity <= extent) Arena::set_capacity past the backing buffer "
       "would hand out bytes the arena does not own"},
      {"arena_rebind_capacity_over_extent",
       case_arena_rebind_capacity_over_extent, "memory.h",
       "(new_capacity <= buffer_extent) Arena::rebind capacity exceeds its "
       "backing buffer"},
      {"resplit_scratch_not_empty", case_resplit_scratch_not_empty,
       "memory.cpp",
       "(scratch_arena_a.get_offset() == 0 && scratch_arena_b.get_offset() "
       "== 0) resplit_arenas: both scratch arenas must be empty"},
      {"arena_set_offset_forward", case_arena_set_offset_forward, "memory.h",
       "(new_offset <= offset) Arena::set_offset: "},
      {"scratch_scope_non_lifo", case_scratch_scope_non_lifo, "memory.h",
       "(arena.get_offset() >= saved_offset) "},
      {"arena_vector_overflow", case_arena_vector_overflow, "memory.h",
       "(element_count < element_capacity) ArenaVector push_back exact "
       "capacity exceeded!"},
      {"arena_vector_emplace_overflow", case_arena_vector_emplace_overflow,
       "memory.h",
       "(element_count < element_capacity) ArenaVector emplace_back exact "
       "capacity exceeded!"},
      {"generate_target_is_scratch", case_generate_target_is_scratch,
       "memory.h",
       "(&target != &scratch_arena_a && &target != &scratch_arena_b) "
       "generate: target must not alias an engine scratch arena"},
      {"generate_recursion_too_deep", case_generate_recursion_too_deep,
       "memory.h",
       "(depth <= MAX_GENERATE_DEPTH) generate: recursion too deep"},
      {"normalize_zero", case_normalize_zero, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"normalize_nan", case_normalize_nan, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"solids_index_oob", case_solids_index_oob, "solids.h",
       "(index < static_cast<size_t>(NUM_ENTRIES)) Solids::get_entry: index "
       "out of range"},
      {"solids_unknown_name", case_solids_unknown_name, "solids.h",
       "(entry) Solids::get_by_name: unknown solid name"},
      {"circular_buffer_oob", case_circular_buffer_oob,
       "static_circular_buffer.h",
       "(index < count) StaticCircularBuffer index out of range"},
      {"circular_buffer_front_empty", case_circular_buffer_front_empty,
       "static_circular_buffer.h",
       "(!is_empty()) front() on empty StaticCircularBuffer"},
      {"arena_vector_append_bulk_overflow",
       case_arena_vector_append_bulk_overflow, "memory.h",
       "(count <= element_capacity - element_count) ArenaVector bulk append "
       "exceeds capacity!"},
      {"spatial_knn_over_max", case_spatial_knn_over_max, "kd_tree.h",
       "(k <= static_cast<size_t>(MAX_K)) KDTree::nearest k exceeds MAX_K"},
      {"arena_oversubscribed", case_arena_oversubscribed, "memory.cpp",
       "(false) "},
      {"persist_forgot_reset", case_persist_forgot_reset, "memory.h",
       "(persistent.get_offset() <= persistent_offset_at_ctor) Persist: "
       "restore grew the persistent arena past its construction watermark — "
       "the caller did not rewind/reset it during the scope, so the restore "
       "appended a duplicate instead of reconstructing"},
      {"persist_same_arena", case_persist_same_arena, "memory.h",
       "(&scratch_arena != &restore_arena) Persist: scratch and persistent "
       "must be distinct arenas — the dtor's watermark restore assumes the "
       "backup lives in a different arena than the one it restores into"},
      {"triangular_bitset_unordered_pair",
       case_triangular_bitset_unordered_pair, "triangular_bitset.h",
       "(small >= 0 && small < large && large < MAX_V) "
       "TriangularBitset::index: pair "},
      {"timeline_handled_relocation", case_timeline_handled_relocation,
       "timeline.h", "(!handled) "},
      {"timeline_move_into_live_destination",
       case_timeline_move_into_live_destination, "timeline.h",
       "(!dst.manager) move_into would leak the destination's live animation"},
      {"timeline_negative_delay", case_timeline_negative_delay, "timeline.h",
       "(in_frames >= 0) Timeline delay must be non-negative"},
      {"finite_param_perpetual_duration", case_finite_param_perpetual_duration,
       "params.h",
       "(duration >= 0) finite parameter animation duration must be >= 0"},
      {"transition_nonfinite_target", case_transition_nonfinite_target,
       "params.h", "(std::isfinite(to)) Transition target must be finite"},
      {"timeline_start_overflow", case_timeline_start_overflow, "timeline.h",
       "(delay <= UINT32_MAX - global_timeline_t) Timeline start frame "
       "overflow"},
      {"timeline_handled_completion", case_timeline_handled_completion,
       "timeline.h", "(!e.handled || anim->is_canceled()) "},
      {"timeline_pinned_finite_animation",
       case_timeline_pinned_finite_animation, "timeline.h",
       "(!animation.is_finite() || animation.repeats()) pinned animation "
       "must be infinite or repeating"},
      {"timeline_pinned_add_on_full_timeline",
       case_timeline_pinned_add_on_full_timeline, "timeline.h",
       "(pin == Pin::UNPINNED) Timeline full, dropped a pinned animation"},
      {"timeline_pinned_one_shot_timer", case_timeline_pinned_one_shot_timer,
       "timeline.h", "(!e.handled || anim->is_canceled()) "},
      {"timeline_clear_pinned", case_timeline_clear_pinned, "timeline.h",
       "(!global_timeline_events[i].handled) clear() would destroy a pinned "
       "animation"},
      {"timeline_clear_during_step", case_timeline_clear_during_step,
       "timeline.h",
       "(!stepping) clear() from inside step() would destroy the animation "
       "whose callback is running"},
      {"timeline_clear_hook_adds_event", case_timeline_clear_hook_adds_event,
       "timeline.h",
       "(global_timeline_num_events == event_count) clear hook added or "
       "removed timeline events"},
      {"mesh_carousel_unflipped_slot", case_mesh_carousel_unflipped_slot,
       "carousel.h",
       "(slot == front) MeshCarousel segue scheduled before incoming slot "
       "flip"},
      {"timeline_double_construct", case_timeline_double_construct,
       "timeline.h", "(!global_timeline_live) "},
      {"transformer_pool_init_storage_twice",
       case_transformer_pool_init_storage_twice, "transformer.h",
       "(!entities) TransformerPool: init_storage() called twice"},
      {"transformer_pool_spawn_before_init",
       case_transformer_pool_spawn_before_init, "transformer.h",
       "(entities) TransformerPool: call init_storage() before spawn"},
      {"transformer_pool_active_index_oob",
       case_transformer_pool_active_index_oob, "transformer.h",
       "(k >= 0 && k < active_slot_count) TransformerPool: active index out of "
       "range"},
      {"transformer_pool_reclaim_storage_moved",
       case_transformer_pool_reclaim_storage_moved, "transformer.h",
       "(e == entities && s == active_slots) TransformerPool: reclaimed "
       "storage moved"},
      {"transformer_pool_outlives_timeline",
       case_transformer_pool_outlives_timeline, "transformer.h",
       "(global_timeline_live) TransformerPool outlived its Timeline: declare "
       "the Timeline before the pools that schedule on it"},
      {"effect_double_construct", case_effect_double_construct, "canvas.h",
       "(!s_alive) Effect: a second Effect was constructed while one is "
       "still alive; buffer_a/buffer_b are shared static storage (one live "
       "Effect only)"},
      {"effect_width_zero", case_effect_width_zero, "canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 0 x "
       "16 are outside 1..288 x 1..144"},
      {"effect_height_zero", case_effect_height_zero, "canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 32 x "
       "0 are outside 1..288 x 1..144"},
      {"effect_width_over_max", case_effect_width_over_max, "canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 289 x "
       "16 are outside 1..288 x 1..144"},
      {"effect_height_over_max", case_effect_height_over_max, "canvas.h",
       "(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H) Effect dimensions 32 x "
       "145 are outside 1..288 x 1..144"},
      {"particle_lifetime_zero", case_particle_lifetime_zero, "sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"particle_lifetime_nan", case_particle_lifetime_nan, "sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"particle_lifetime_over_max", case_particle_lifetime_over_max,
       "sprites.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem max_life must be finite and in [1, 65535]"},
      {"particle_render_zero_lifetime", case_particle_render_zero_lifetime,
       "particles.h",
       "(std::isfinite(max_life) && max_life >= 1.0f && max_life <= "
       "65535.0f) ParticleSystem render max_life must be finite and in [1, "
       "65535]"},
      {"correction_guard_double_construct",
       case_correction_guard_double_construct, "led.h",
       "(!correction_guard_live()) at most one correction guard may be "
       "live at a time (see contract above)"},
      {"correction_guard_cross_type", case_correction_guard_cross_type, "led.h",
       "(!correction_guard_live()) at most one correction guard may be "
       "live at a time (see contract above)"},
      {"mesh_narrow_index", case_mesh_narrow_index, "mesh.h",
       "(i <= static_cast<size_t>(INT16_MAX)) mesh index exceeds int16_t "
       "topology range (oversized mesh?)"},
      {"medial_aliases_input", case_medial_aliases_input, "conway.h",
       "(&mesh != &out_a) medial input mesh must not alias output mesh"},
      {"needle_aliased_arenas", case_needle_aliased_arenas, "conway.h",
       "(&target != &temp) needle: target and temp must differ"},
      {"zip_aliased_arenas", case_zip_aliased_arenas, "conway.h",
       "(&target != &temp) zip: target and temp must differ"},
      {"gyro_aliased_arenas", case_gyro_aliased_arenas, "conway.h",
       "(&target != &temp) gyro: target and temp must differ"},
      {"bevel_aliased_arenas", case_bevel_aliased_arenas, "conway.h",
       "(&target != &temp) bevel: target and temp must differ"},
      {"mesh_transform_aliases_source", case_mesh_transform_aliases_source,
       "conway.h",
       "(&mesh != &transformed) MeshOps::transform source mesh must not alias "
       "the destination"},
      {"chamfer_collapsed_endpoint", case_chamfer_collapsed_endpoint,
       "conway.h", "(t >= 0.0f && t < 1.0f) chamfer: t out of [0,1)"},
      {"snub_collapsed_endpoint", case_snub_collapsed_endpoint, "conway.h",
       "(t >= 0.0f && t < 1.0f) snub: t out of [0,1)"},
      {"conway_empty_mesh", case_conway_empty_mesh, "mesh.h",
       "(total_indices > 0) half-edge mesh requires at least one face index"},
      {"conway_degenerate_mesh", case_conway_degenerate_mesh, "mesh.h",
       "(he_mesh.half_edges[i].pair != HE_NONE) MeshOps::truncate requires a "
       "closed manifold (unpaired half-edge)"},
      {"conway_target_exhausted", case_conway_target_exhausted, "memory.h",
       "(false) "},
      {"relax_baked_dimension_mismatch", case_relax_baked_dimension_mismatch,
       "conway.h",
       "(V == bake.vertex_count && F == bake.face_count && I == "
       "bake.index_count) relax_baked: source dimensions differ"},
      {"relax_baked_topology_mismatch", case_relax_baked_topology_mismatch,
       "conway.h",
       "(relax_topology_hash(mesh) == bake.topology_hash) relax_baked: source "
       "topology differs"},
      {"relax_baked_output_hash_mismatch",
       case_relax_baked_output_hash_mismatch, "conway.h",
       "(output_hash == bake.output_hash) relax_baked: output hash differs"},
      {"half_edge_face_counts_short", case_half_edge_face_counts_short,
       "mesh.h",
       "(counted_indices == total_indices) mesh face counts do not span flat "
       "index length"},
      {"half_edge_face_counts_long", case_half_edge_face_counts_long, "mesh.h",
       "(count <= total_indices - counted_indices) mesh face counts exceed "
       "flat index length"},
      {"mesh_compile_face_counts_short", case_mesh_compile_face_counts_short,
       "mesh.h",
       "(counted_indices == total_indices) mesh face counts do not span flat "
       "index length"},
      {"mesh_compile_face_counts_long", case_mesh_compile_face_counts_long,
       "mesh.h",
       "(count <= total_indices - counted_indices) mesh face counts exceed "
       "flat index length"},
      {"mesh_compile_face_span_over_16bit",
       case_mesh_compile_face_span_over_16bit, "mesh.h",
       "(current_offset + count <= UINT16_MAX) mesh face_offsets exceeds "
       "16-bit index range"},
      {"update_hankin_stale_topology", case_update_hankin_stale_topology,
       "hankin.h",
       "(out_mesh.topology.size() == 0 || out_mesh.topology_key == "
       "compiled.topology_key) update_hankin: reused out_mesh carries "
       "a topology from a different compiled pattern (clear it first)"},
      {"update_hankin_dual_seed_topology",
       case_update_hankin_dual_seed_topology, "hankin.h",
       "(out_mesh.topology.size() == 0 || out_mesh.topology_key == "
       "compiled.topology_key) update_hankin: reused out_mesh carries "
       "a topology from a different compiled pattern (clear it first)"},
      {"hankin_clone_aliases_dst", case_hankin_clone_aliases_dst, "hankin.h",
       "(&src != &dst) CompiledHankin::clone src must not alias dst"},
      {"mesh_state_set_borrowed_offsets_count_mismatch",
       case_mesh_state_set_borrowed_offsets_count_mismatch, "mesh_state.h",
       "(face_offsets_span.size() == face_counts_span.size()) "
       "MeshState::set_borrowed: one face offset per face count required"},
      {"mesh_state_set_borrowed_offsets_short_span",
       case_mesh_state_set_borrowed_offsets_short_span, "mesh_state.h",
       "(static_cast<size_t>(face_offsets_span[last]) + "
       "face_counts_span[last] == faces_span.size()) MeshState::set_borrowed: "
       "face offsets do not span faces"},
      {"mesh_state_set_borrowed_offsets_not_prefix_sum",
       case_mesh_state_set_borrowed_offsets_not_prefix_sum, "mesh_state.h",
       "(offsets_are_prefix_sum(face_counts_span, face_offsets_span)) "
       "MeshState::set_borrowed: face offsets are not the prefix sum of the "
       "face counts"},
      {"half_edge_zero_side_face", case_half_edge_zero_side_face, "mesh.h",
       "(count > 0) half-edge mesh face has zero sides"},
      {"half_edge_non_manifold_edge", case_half_edge_non_manifold_edge,
       "mesh.h", "(j - i <= 2) non-manifold edge: >2 half-edges share an edge"},
      {"half_edge_inconsistent_winding", case_half_edge_inconsistent_winding,
       "mesh.h",
       "(out.half_edges[a].vertex != out.half_edges[b].vertex) half-edge "
       "mesh faces are inconsistently wound"},
      {"mesh_narrow_face_count", case_mesh_narrow_face_count, "mesh.h",
       "(count >= 0 && count <= UINT8_MAX) mesh face side count exceeds "
       "uint8_t range"},
      {"mesh_require_closed_manifold", case_mesh_require_closed_manifold,
       "mesh.h",
       "(he_mesh.half_edges[i].pair != HE_NONE) MeshOps::death requires a "
       "closed manifold (unpaired half-edge)"},
      {"mesh_require_vertex_manifold", case_mesh_require_vertex_manifold,
       "mesh.h",
       "(walked == fan_size[origin]) MeshOps::death requires a vertex "
       "manifold (split vertex fan)"},
      {"mesh_require_matching_face_sides",
       case_mesh_require_matching_face_sides, "mesh.h",
       "(static_cast<size_t>(he_mesh.faces[fi].half_edge) == face_offset) "
       "MeshOps::death: half-edge mesh face sides differ from the source mesh"},
      {"apply_step_hankin_no_angle", case_apply_step_hankin_no_angle,
       "recipe.h",
       "(step.param > 0.0f) apply_step: HANKIN step has no contact angle"},
      {"apply_step_bevel_no_depth", case_apply_step_bevel_no_depth, "recipe.h",
       "(step.param > 0.0f) apply_step: BEVEL step has no depth"},
      {"slerp_nan", case_slerp_nan, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"make_rotation_vectors_nan", case_make_rotation_vectors_nan, "3dmath.h",
       "(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ && "
       "std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ) "
       "make_rotation(from, to): inputs must be unit vectors"},
      {"make_rotation_angle_nan", case_make_rotation_angle_nan, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"make_rotation_nonunit", case_make_rotation_nonunit, "3dmath.h",
       "(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ && "
       "std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ) "
       "make_rotation(from, to): inputs must be unit vectors"},
      {"make_basis_nan", case_make_basis_nan, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"noise_transform_nan", case_noise_transform_nan, "3dmath.h",
       "(m2 >= math::EPS_NORMALIZE_SQ) "},
      {"driver_null_speed_src", case_driver_null_speed_src, "params.h",
       "(speed_src != nullptr) Driver: live speed_src is null"},
      {"path_append_zero_samples", case_path_append_zero_samples, "motion.h",
       "(samples >= 1) "},
      {"motion_empty_path_origin_sample", case_motion_empty_path_origin_sample,
       "motion.h",
       "(dot(current_v, current_v) >= math::EPS_LEN_SQ && dot(target_v, "
       "target_v) >= math::EPS_LEN_SQ) Motion: path sampled at the origin "
       "(empty or origin-crossing path)"},
      {"alpha_falloff_null", case_alpha_falloff_null, "composition.h",
       "(fn != nullptr) AlphaFalloffShade: falloff function must not be null"},
      {"noise_hue_palette_null_source", case_noise_hue_palette_null_source,
       "noise_hue_palette.h",
       "(source != nullptr) NoiseHuePalette bound to null source"},
      {"noise_hue_palette_null_rotation_lut",
       case_noise_hue_palette_null_rotation_lut, "noise_hue_palette.h",
       "(hue_rotation_lut != nullptr) NoiseHuePalette bound to null "
       "hue-rotation LUT"},
      {"noise_hue_palette_null_noise_lut",
       case_noise_hue_palette_null_noise_lut, "noise_hue_palette.h",
       "(hue_noise_lut != nullptr) NoiseHuePalette bound to null hue-noise "
       "LUT"},
      {"baked_palette_rebake_aliased", case_baked_palette_rebake_aliased,
       "composition.h",
       "(!aliased) BakedPalette::rebake through an aliasing handle"},
      {"baked_palette_clone_from_self", case_baked_palette_clone_from_self,
       "composition.h", "(&src != this) BakedPalette::clone_from from itself"},
      {"baked_palette_bake_blend_self", case_baked_palette_bake_blend_self,
       "composition.h",
       "(&from != this && &to != this) BakedPalette::bake_blend endpoint is "
       "the output"},
      {"register_param_overflow", case_register_param_overflow, "canvas.h",
       "(parameters.count < parameters.capacity()) register_param: "
       "exceeded ParamList capacity"},
      {"register_int_param_range", case_register_int_param_range, "canvas.h",
       "(range_fits) register_int_param: [min,max] must fit the target "
       "integer type"},
      {"set_clip_out_of_bounds", case_set_clip_out_of_bounds, "canvas.h",
       "(y0 >= 0 && y0 <= y1 && y1 <= clip_region.h && x0 >= 0 && x0 <= x1 "
       "&& x1 <= clip_region.w) set_clip band must be non-inverted and "
       "within canvas bounds"},
      {"set_clip_mid_frame", case_set_clip_mid_frame, "canvas.h",
       "(!canvas_active) clip cannot change while a frame is active"},
      {"output_envelope_out_of_range", case_output_envelope_out_of_range,
       "canvas.h",
       "(std::isfinite(value) && value >= 0.0f && value <= 1.0f) output "
       "envelope must be finite and in [0,1]"},
      {"set_clip_x_out_of_bounds", case_set_clip_x_out_of_bounds, "canvas.h",
       "(x0 >= 0 && x0 <= x1 && x1 <= clip_region.w) set_clip_x band must be "
       "non-inverted and within canvas width"},
      {"arcs_overlap_start_out_of_range", case_arcs_overlap_start_out_of_range,
       "clip.h", "(s1 >= 0 && s1 < w && s2 >= 0 && s2 < w) "},
      {"scan_clip_out_of_bounds", case_scan_clip_out_of_bounds, "shader.h",
       "(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 && "
       "cr.render_y_end() <= H) "},
      {"scan_clip_rows_out_of_bounds", case_scan_clip_rows_out_of_bounds,
       "shader.h",
       "(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 && "
       "cr.render_y_end() <= H) "},
      {"face_scratch_retargeted", case_face_scratch_retargeted, "face.h",
       "(!scratch_owner || scratch_owner->claim_seq == scratch_claim) "
       "SDF::Face scanned after a later Face claimed its scratch buffer"},
      {"scan_canvas_dim_mismatch", case_scan_canvas_dim_mismatch, "raster.h",
       "(canvas.width() == W && canvas.height() == H) "},
      {"scan_pipeline_not_prepared", case_scan_pipeline_not_prepared,
       "raster.h",
       "(pipeline.prepared_for(canvas)) direct raster pipeline not prepared "
       "for this canvas"},
      {"pipeline_ref_erase_not_prepared", case_pipeline_ref_erase_not_prepared,
       "concepts.h",
       "(t.prepared_for(cv)) direct raster pipeline not prepared for this "
       "canvas"},
      {"scan_mesh_face_index_out_of_range",
       case_scan_mesh_face_index_out_of_range, "mesh.h",
       "(static_cast<size_t>(faces[k]) < num_verts) mesh face index exceeds "
       "the vertex pool"},
      {"plot_mesh_vertex_over_capacity", case_plot_mesh_vertex_over_capacity,
       "mesh.h", "(large < DEDUP_CAPACITY) "},
      {"plot_extract_edges_vertex_over_capacity",
       case_plot_extract_edges_vertex_over_capacity, "mesh.h",
       "(large < DEDUP_CAPACITY) "},
      {"feedback_downsample_indivisible", case_feedback_downsample_indivisible,
       "pixel_feedback.h",
       "(downsample > 0 && W % downsample == 0) feedback downsample 5 must "
       "be > 0 and divide width 32"},
      {"screen_trails_set_lifetime_nonpositive",
       case_screen_trails_set_lifetime_nonpositive, "screen_trails.h",
       "(new_lifetime > 0) Screen::Trails: lifetime 0 must be positive"},
      {"spherical_field_ring_index_oob", case_spherical_field_ring_index_oob,
       "spherical_field.h",
       "(y < H - 1) SphericalFieldLayout: ring index 5 out of range"},
      {"spherical_field_populate_ring_end_oob",
       case_spherical_field_populate_ring_end_oob, "spherical_field.h",
       "(ring_end < layout.ring_count()) SphericalField::populate: ring_end 5 "
       "past the last ring 4"},
      {"spherical_harmonic_normalization_overflow",
       case_spherical_harmonic_normalization_overflow, "spherical_harmonics.h",
       "(l + abs_m <= MAX_FACTORIAL_ARGUMENT) spherical harmonic "
       "normalization: l + |m| = 40 collapses the float factorial ratio"},
      {"spherical_harmonic_order_over_degree",
       case_spherical_harmonic_order_over_degree, "spherical_harmonics.h",
       "(l >= 0 && abs_m <= l) spherical harmonic normalization: order 3 is "
       "outside [-2, 2]"},
      {"spherical_harmonic_decode_negative_index",
       case_spherical_harmonic_decode_negative_index, "spherical_harmonics.h",
       "(idx >= 0) decode_lm: flat index -1 is negative"},
      {"spherical_field_infill_over_domain",
       case_spherical_field_infill_over_domain, "spherical_field.h",
       "(north_infill >= 0 && south_infill >= 0 && north_infill + south_infill "
       "<= H) SphericalFieldLayout: infills 0 + 17 must be non-negative and "
       "fit within H = 16"},
      {"spherical_field_negative_equator_samples",
       case_spherical_field_negative_equator_samples, "spherical_field.h",
       "(equator_samples >= 0) SphericalFieldLayout: equator_samples -1 must "
       "be >= 0"},
      {"feedback_negative_fade", case_feedback_negative_fade, "styles.h",
       "(std::isfinite(fade) && fade >= 0.0f) Feedback::Style::fade must be "
       "finite and >= 0"},
      {"feedback_infinite_fade", case_feedback_infinite_fade, "styles.h",
       "(std::isfinite(fade) && fade >= 0.0f) Feedback::Style::fade must be "
       "finite and >= 0"},
      {"gradient_no_stops", case_gradient_no_stops, "color.h",
       "(points.size() > 0) Gradient requires at least one stop"},
      {"gradient_stop_out_of_range", case_gradient_stop_out_of_range, "color.h",
       "(stop.first >= 0.0f && stop.first <= 1.0f) Gradient stop position "
       "out of [0,1]"},
      {"gradient_stops_unsorted", case_gradient_stops_unsorted, "color.h",
       "(stop.first >= prev_check) Gradient stops must be sorted ascending"},
      {"random_timer_inverted_range", case_random_timer_inverted_range,
       "timers.h", "(min >= 0 && min <= max) "},
      {"dma_controller_wedged_overrun", case_dma_controller_wedged_overrun,
       "test_death.h", "(false) DMA channel wedged"},
      {"empty_fn_call", case_empty_fn_call, "inplace_function.h",
       "(vtable != empty) empty hs::inplace_function called"},
      {"effect_registry_duplicate_name", case_effect_registry_duplicate_name,
       "registry.h",
       "(existing.name != reg.name) effect header included by more than one "
       "translation unit: effects/DeathDuplicate.h"},
      {"flywheel_period_zero", case_flywheel_period_zero, "pov_sync_flywheel.h",
       "(p > 0 && p <= static_cast<uint32_t>(INT32_MAX) / MIN_SAFE_HALF_REVS) "
       "Flywheel: cycles_per_half_rev outside the range position()'s int32 "
       "elapsed window holds for MIN_SAFE_HALF_REVS of coast"},
      {"symbol_emitter_overlapping_burst",
       case_symbol_emitter_overlapping_burst, "pov_sync_emitter.h",
       "(pulses_left == 0 && queue_pos >= queue_len) "
       "SymbolEmitter::schedule_boundary: wire busy"},
      {"y_to_phi_degenerate_height", case_y_to_phi_degenerate_height,
       "geometry.h", "(h_virt > 1) y_to_phi_virtual: h_virt must be > 1"},
      {"orientation_frame_index_oob", case_orientation_frame_index_oob,
       "geometry.h", "(i >= 0 && i < num_frames) "},
      {"make_basis_nonunit_quaternion", case_make_basis_nonunit_quaternion,
       "geometry.h",
       "(std::abs(orientation.squared_magnitude() - 1.0f) < "
       "math::EPS_UNIT_QUAT_SQ) "},
      {"sdf_polygon_side_count", case_sdf_polygon_side_count, "shapes.h",
       "(sides >= 3) "},
      {"sdf_spherical_polygon_radius_over_hemisphere",
       case_sdf_spherical_polygon_radius_over_hemisphere, "shapes.h",
       "(radius <= 1.0f) "},
      {"sdf_flower_radius_over_hemisphere",
       case_sdf_flower_radius_over_hemisphere, "shapes.h", "(radius <= 1.0f) "},
      {"sdf_flower_zero_radius", case_sdf_flower_zero_radius, "shapes.h",
       "(radius > 0.0f) "},
      {"sdf_class_lut_too_few_vertices", case_sdf_class_lut_too_few_vertices,
       "face.h",
       "(count >= 3) build_canonical_distance_lut requires at least 3 polygon "
       "vertices"},
      {"sdf_class_lut_grid_too_small", case_sdf_class_lut_grid_too_small,
       "face.h",
       "(n >= 2) build_canonical_distance_lut requires a grid resolution of at "
       "least 2"},
      {"sdf_bind_class_lut_offset_out_of_range",
       case_sdf_bind_class_lut_offset_out_of_range, "face.h",
       "(vert_offset >= 0 && vert_offset < count) bind_class_lut: vertex "
       "offset outside the face"},
      {"sdf_distorted_ring_null_shift", case_sdf_distorted_ring_null_shift,
       "rings.h", "(sf) DistortedRing: shift_fn must be non-null"},
      {"shader_workbench_empty_preset_view",
       case_shader_workbench_empty_preset_view, "ShaderWorkbench.h",
       "(!source_indices.empty()) set_fixed_preset_view: empty preset view"},
      {"shader_workbench_preset_view_index_out_of_range",
       case_shader_workbench_preset_view_index_out_of_range,
       "ShaderWorkbench.h",
       "(index < PRESETS.size()) set_fixed_preset_view: preset index out of "
       "range"},
      {"chain_table_rank_decreases", case_chain_table_rank_decreases,
       "interpreter.h",
       "(entry.input <= entry.output) ChainProgram::bind_storage: operator "
       "family rank decreases"},
      {"shader_workbench_preset_for_view_out_of_range",
       case_shader_workbench_preset_for_view_out_of_range, "ShaderWorkbench.h",
       "(index < preset_count_for_view()) preset_for_view: index out of "
       "range"},
      {"sdf_angular_repeat_nonunit_axis", case_sdf_angular_repeat_nonunit_axis,
       "csg.h", "(fabsf(ax.length() - 1.0f) < 1e-3f) "},
      {"sdf_distorted_ring_zero_knots", case_sdf_distorted_ring_zero_knots,
       "rings.h", "(kn != nullptr && n >= 1) "},
      {"sdf_twist_zero_major_radius", case_sdf_twist_zero_major_radius,
       "volume.h", "(R > 0.0f) "},
      {"pick_next_edge_unknown_node", case_pick_next_edge_unknown_node,
       "conway_graph.h", "(n > 0) pick_next_edge: node outside the graph"},
      {"opleg_edge_sweep_no_edge", case_opleg_edge_sweep_no_edge, "opleg.h",
       "(spec.edge) OpLeg: edge sweep carries no graph edge"},
      {"opleg_zero_sweep_frames", case_opleg_zero_sweep_frames, "opleg.h",
       "(spec.sweep_frames >= 1) OpLeg needs a positive sweep length"},
      {"opleg_incomplete_palette_handoff",
       case_opleg_incomplete_palette_handoff, "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: param sweep leg has an incomplete palette handoff"},
      {"opleg_edge_settle_mismatch", case_opleg_edge_settle_mismatch, "opleg.h",
       "(spec.settle_frames >= 0 && edge.settle == (spec.settle_frames > 0)) "
       "OpLeg: settle frames disagree with the edge"},
      {"opleg_edge_sweep_incomplete_handoff",
       case_opleg_edge_sweep_incomplete_handoff, "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: edge sweep leg has an incomplete palette handoff"},
      {"opleg_hankin_incomplete_handoff", case_opleg_hankin_incomplete_handoff,
       "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: hankin sweep leg has an incomplete palette handoff"},
      {"opleg_hankin_backward_theta", case_opleg_hankin_backward_theta,
       "opleg.h",
       "(spec.theta_start <= spec.theta_end) OpLeg: hankin leg sweeps back to "
       "a smaller contact angle"},
      {"opleg_relax_no_iterations", case_opleg_relax_no_iterations, "opleg.h",
       "(spec.bake || spec.iterations >= 1) OpLeg: relax leg needs a positive "
       "iteration count"},
      {"opleg_relax_incomplete_handoff", case_opleg_relax_incomplete_handoff,
       "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: relax leg has an incomplete palette handoff"},
      {"opleg_medial_incomplete_handoff", case_opleg_medial_incomplete_handoff,
       "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: medial leg has an incomplete palette handoff"},
      {"opleg_reconcile_no_endpoints", case_opleg_reconcile_no_endpoints,
       "opleg.h",
       "(spec.to_positions) OpLeg: reconcile leg carries no "
       "endpoints"},
      {"opleg_reconcile_incomplete_handoff",
       case_opleg_reconcile_incomplete_handoff, "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: reconcile leg has an incomplete palette handoff"},
      {"opleg_gated_swap_zero_gate_frames",
       case_opleg_gated_swap_zero_gate_frames, "opleg.h",
       "(spec.gate_frames >= 1) OpLeg needs a positive gate length"},
      {"opleg_gated_swap_incomplete_handoff",
       case_opleg_gated_swap_incomplete_handoff, "opleg.h",
       "(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0) "
       "OpLeg: gated swap leg has an incomplete palette handoff"},
      {"opleg_shading_face_out_of_range", case_opleg_shading_face_out_of_range,
       "opleg.h", "(face < faces) OpLeg::Shading: ramp face out of range"},
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
inline constexpr const char *SHAPE_PROBE_CASE = "__shape_probe__";

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
  if (std::strcmp(name, SHAPE_PROBE_CASE) == 0) {
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
 * @details Keeps the LAST CHILD_OUTPUT_CAP-1 bytes: check_fail() flushes its
 *          breadcrumb immediately before trapping, so it is the final text a
 *          trapping child writes.
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
 * @param file Basename check_fail() logs for the guard's source file.
 * @param text Expected "(condition) message" tail of the breadcrumb line.
 * @return True iff some captured line reads
 *         "HS_CHECK failed: <file>:<line>: <text>...".
 * @details The line number is skipped rather than pinned, so editing above a
 *          guard does not churn the table while the guard's identity stays
 *          exact.
 */
inline bool breadcrumb_names_guard(const char *out, const char *file,
                                   const char *text) {
  char prefix[128];
  std::snprintf(prefix, sizeof(prefix), "HS_CHECK failed: %s:", file);
  const size_t prefix_len = std::strlen(prefix);
  const size_t text_len = std::strlen(text);
  for (const char *p = std::strstr(out, prefix); p;
       p = std::strstr(p + 1, prefix)) {
    const char *q = p + prefix_len;
    while (*q >= '0' && *q <= '9')
      ++q;
    if (q[0] != ':' || q[1] != ' ')
      continue;
    if (std::strncmp(q + 2, text, text_len) == 0)
      return true;
  }
  return false;
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
 * @details Child stdout/stderr are redirected to child_capture_path() and the
 *          tail is loaded into child_output(), so the caller can require the
 *          HS_CHECK breadcrumb of the guard the case is supposed to fire.
 */
inline int spawn_child(const char *name) {
  set_case_env(name);
  const char *capture_path = child_capture_path();
  // Drop the previous spawn's capture up front: a spawn that fails before the
  // redirect takes effect must leave an EMPTY capture, never the last child's
  // breadcrumb, which would read as a pass for this case.
  child_output()[0] = '\0';
  std::remove(capture_path);
#if defined(_WIN32)
  // Shell-free spawn: hand argv straight to the CRT so no cmd.exe parsing can
  // mangle a self_exe() path containing &, %, ^, or quotes. Child stdout/stderr
  // reach the capture file by redirecting fds 1/2 across the synchronous
  // _P_WAIT spawn (the child inherits the CRT fd table), then the descriptors
  // are restored.
  std::fflush(stdout);
  std::fflush(stderr);
  int saved_out = _dup(1);
  int saved_err = _dup(2);
  int capture = -1;
  // Redirect only when both originals were saved: a failed _dup leaves no way
  // to restore, so skipping the redirect keeps the parent's streams intact
  // (a muted parent would silence reporting for every remaining case).
  if (saved_out >= 0 && saved_err >= 0) {
    _sopen_s(&capture, capture_path, _O_WRONLY | _O_CREAT | _O_TRUNC,
             _SH_DENYNO, _S_IREAD | _S_IWRITE);
    if (capture >= 0) {
      _dup2(capture, 1);
      _dup2(capture, 2);
    }
  }
  const char *argv[] = {self_exe(), nullptr};
  intptr_t rc = _spawnv(_P_WAIT, self_exe(), argv);
  // capture is only open when both saves succeeded, so the restore is reached
  // only with valid descriptors.
  if (capture >= 0) {
    _dup2(saved_out, 1);
    _dup2(saved_err, 2);
    _close(capture);
  }
  if (saved_out >= 0)
    _close(saved_out);
  if (saved_err >= 0)
    _close(saved_err);
  load_child_output();
  return static_cast<int>(rc);
#else
  // Shell-free spawn: fork and execv the binary directly so no /bin/sh parsing
  // can mangle a self_exe() path containing a quote or shell metacharacter. The
  // child sends stdout/stderr to the capture file and execs; the parent waits
  // and returns the raw wait status that classify_trap() decodes.
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
  if (waitpid(pid, &status, 0) < 0)
    return -1;
  load_child_output();
  return status;
#endif
}

#if defined(_WIN32)
/**
 * @brief The Windows EXCEPTION_ILLEGAL_INSTRUCTION process exit code.
 * @details An unhandled trap sets it as the process exit code, which _spawnv()
 *          returns directly to the parent.
 */
inline constexpr int TRAP_STATUS = static_cast<int>(0xC000001D);
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
  return rc == TRAP_STATUS ? TrapShape::Signal : TrapShape::None;
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
    // trap coverage; the banner is unmistakable and CI is the hard gate. It is
    // a whole-suite skip so the module's zero-assertion guard does not fire.
    hs_test::skip_suite();
    std::printf("  [SKIPPED] death tests: %s (rc=%d) — 0 trap cases executed; "
                "CI is the hard gate\n",
                why, rc);
  }
}

/**
 * @brief Counts the distinct guard sites in @p file the case table pins.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param file Source-file basename to count pins for.
 * @return Distinct pins naming that file, never above its real site count.
 * @details Distinct by (file, condition text), which is what
 *          breadcrumb_names_guard() compares: guards whose breadcrumbs read
 *          identically (the same condition and message repeated in several
 *          constructors) are one covered site, because no case can prove which
 *          of them fired.
 */
inline int pinned_guards_in(const Case *cs, int n, const char *file) {
  int pinned = 0;
  for (int i = 0; i < n; ++i) {
    if (std::strcmp(cs[i].guard_file, file) != 0)
      continue;
    bool duplicate = false;
    for (int j = 0; j < i && !duplicate; ++j)
      duplicate = std::strcmp(cs[j].guard_file, file) == 0 &&
                  std::strcmp(cs[j].guard_text, cs[i].guard_text) == 0;
    if (!duplicate)
      ++pinned;
  }
  return pinned;
}

/**
 * @brief Ratchet floor on the HS_CHECK sites the case table pins.
 * @details Raise it after adding cases; lower it only alongside a deliberate
 *          removal of the engine guards those cases target.
 */
constexpr int MIN_COVERED_GUARD_SITES = 138;

/** @brief One file's approved count of guard sites no case pins. */
struct GuardGapAllowance {
  const char *file; /**< Source-file basename, as the census names it. */
  int gap;          /**< Unpinned HS_CHECK sites tolerated in that file. */
};

/**
 * @brief Per-file allowances for the sites the case table leaves unpinned.
 * @details MIN_COVERED_GUARD_SITES is a whole-suite total, so a subsystem can
 *          land a dozen guards with no case at all and still clear it. Every
 *          census file's gap is gated against its row here, and a file with no
 *          row must be fully pinned — so new guards red the death module unless
 *          the same commit either pins them or writes the wider gap down. Each
 *          row is exact in both directions: a row that over-approves reds the
 *          module too, so pinning a guard or deleting one forces the row down
 *          in the same commit instead of leaving an allowance nothing spends.
 */
inline constexpr GuardGapAllowance GUARD_GAP_ALLOW[] = {
    {"animation.h", 3},
    {"carousel.h", 3},
    {"motion.h", 6},
    {"opleg.h", 42},
    {"params.h", 11},
    {"segue.h", 1},
    {"sprites.h", 10},
    {"timeline.h", 9},
    {"timers.h", 1},
    {"color.h", 1},
    {"composition.h", 33},
    {"generative_palette.h", 4},
    {"palette_cycler.h", 8},
    {"choreography.h", 1},
    {"memory.cpp", 1},
    {"memory.h", 3},
    {"reaction_graph.h", 2},
    {"static_circular_buffer.h", 4},
    {"transformer.h", 4},
    {"3dmath.h", 5},
    {"geometry.h", 16},
    {"lenses.h", 2},
    {"spherical_field.h", 2},
    {"waves.h", 1},
    {"conway.h", 31},
    {"conway_graph.h", 1},
    {"hankin.h", 9},
    {"mesh.h", 21},
    {"mesh_classes.h", 6},
    {"mesh_state.h", 3},
    {"recipe.h", 13},
    {"solid_generators.h", 5},
    {"solids.h", 2},
    {"kd_tree.h", 5},
    {"canvas.h", 29},
    {"common.h", 4},
    {"csg.h", 2},
    {"cull.h", 1},
    {"face.h", 3},
    {"interpreter.h", 13},
    {"led.h", 3},
    {"operator_model.h", 1},
    {"operators.h", 1},
    {"pixel_feedback.h", 7},
    {"raster.h", 9},
    {"screen_trails.h", 2},
    {"shader.h", 2},
    {"shading.h", 1},
    {"shapes.h", 32},
    {"volume.h", 7},
    {"world_trails.h", 3},
    {"Fishbowl.h", 1},
    {"Comets.h", 2},
    {"DisplacementField.h", 1},
    {"DreamBalls.h", 11},
    {"Dynamo.h", 1},
    {"GnomonicStars.h", 1},
    {"HankinSolids.h", 13},
    {"IslamicStars.h", 28},
    {"MeshFeedback.h", 1},
    {"MindSplatter.h", 3},
    {"MobiusRings.h", 1},
    {"ReactionDiffusionBase.h", 2},
    {"RingShower.h", 1},
    {"ShaderWorkbench.h", 14},
    {"ShaderChain.h", 1},
    {"ShapeShifter.h", 2},
    {"dma_led.h", 4},
    {"pov_segmented.h", 9},
    {"pov_single.h", 9},
    {"Profile.ino", 3},
    {"phantasm_target.h", 2},
    {"engine_bindings.h", 5},
    {"mesh_ops_bindings.h", 2},
};

/**
 * @brief Looks up a file's approved unpinned-site count.
 * @param file Source-file basename from the census.
 * @return The approved gap, or 0 for a file with no allowance row.
 */
inline int allowed_guard_gap(const char *file) {
  for (const GuardGapAllowance &a : GUARD_GAP_ALLOW)
    if (std::strcmp(a.file, file) == 0)
      return a.gap;
  return 0;
}

/**
 * @brief Prints what fraction of the engine's fail-fast surface is pinned.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @details Both sides are derived, never written down: the denominator is the
 *          generated HS_CHECK census (death_guard_sites.h) and the numerator is
 *          the case table itself, so neither can drift from what it measures.
 *          Cases pinning a file the census does not know — the harness's own
 *          trap stand-ins, and any mistyped basename — count in neither and are
 *          reported separately rather than silently dropped. The pinned count is
 *          gated against MIN_COVERED_GUARD_SITES and each file's unpinned
 *          remainder against GUARD_GAP_ALLOW; the ratio itself is reported but
 *          not gated, since new engine guards move the denominator without
 *          weakening any case.
 */
inline void report_guard_coverage(const Case *cs, int n) {
  int covered = 0;
  int off_census = 0;
  int unapproved_gaps = 0;
  int stale_allowances = 0;
  constexpr int GAPS = 5;
  const GuardSiteCount *worst[GAPS] = {};
  int worst_gap[GAPS] = {};
  for (const GuardSiteCount &f : GUARD_SITE_COUNTS) {
    int pinned = pinned_guards_in(cs, n, f.file);
    if (pinned > f.sites)
      pinned = f.sites;
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
    off_census += in_census ? 0 : 1;
  }
  std::printf("  guard coverage: %d/%d HS_CHECK sites pinned by a case (%d%%), "
              "%d case(s) outside the census\n",
              covered, GUARD_SITE_TOTAL,
              GUARD_SITE_TOTAL ? covered * 100 / GUARD_SITE_TOTAL : 0,
              off_census);
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
  HS_EXPECT_GE(covered, MIN_COVERED_GUARD_SITES);
  HS_EXPECT_EQ(unapproved_gaps, 0);
  HS_EXPECT_EQ(stale_allowances, 0);
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

  // Exact roster size, so a silently dropped case fails here rather than
  // hiding under slack. Update when adding or removing cases.
  constexpr int DEATH_CASE_COUNT = 178;
  HS_EXPECT_EQ(n, DEATH_CASE_COUNT);

  // Probe how a trap is relayed (direct SIGILL vs an exit 128+SIGILL) with a
  // dedicated always-trapping sentinel rather than a real case. A real case that
  // regressed to not trapping would otherwise corrupt shape detection and skip
  // the whole suite, instead of failing just that case in the loop below. The
  // sentinel traps through the same HS_CHECK path, so its shape matches the cases.
  TrapShape shape = classify_trap(spawn_child(SHAPE_PROBE_CASE));
  if (shape == TrapShape::None) {
    report_unrunnable(
        "trap-shape sentinel did not trap; cannot classify trap status", 0);
    set_case_env("");
    return fixture.result();
  }

  // The same sentinel proves the capture channel: without the child's
  // breadcrumb every per-case guard check below would fail identically and
  // point at the cases instead of at the broken redirect.
  if (!breadcrumb_names_guard(child_output(), "test_death.h",
                              "(false) death-harness trap-shape probe")) {
    report_unrunnable("cannot capture child output; which guard fired is "
                      "unverifiable",
                      0);
    set_case_env("");
    return fixture.result();
  }

  for (int i = 0; i < n; ++i) {
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc, shape);
    // Dying is not enough: the child must die at THIS case's guard. Any other
    // trap — UB lowered to the same illegal instruction, or a guard the case
    // hits on its way to the one it targets — fails here.
    bool at_guard = breadcrumb_names_guard(child_output(), cs[i].guard_file,
                                           cs[i].guard_text);
    HS_EXPECT_TRUE(trapped);
    HS_EXPECT_TRUE(at_guard);
    std::printf("  [%s] trap fires: %-26s (child rc=%d)\n",
                trapped && at_guard ? "ok" : "FAIL", cs[i].name, rc);
    if (trapped && !at_guard)
      std::printf("      expected %s: %s\n      child logged: %s\n",
                  cs[i].guard_file, cs[i].guard_text, child_output());
  }

  std::remove(child_capture_path());
  set_case_env(""); // leave the env clean for anything that runs after us
  report_guard_coverage(cs, n);
  return fixture.result();
}

} // namespace death_tests
} // namespace hs_test
