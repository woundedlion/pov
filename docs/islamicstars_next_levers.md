# Brief: IslamicStars — remaining optimization levers

You are picking up performance work on `c:\work\Holosphere` (Holosphere /
Phantasm POV renderer, C++20, Teensy 4.0 @ 600 MHz Cortex-M7). Assume none of
the prior sessions' context. Read this whole brief before touching anything.

The effect is **23 green / 0 yellow / 1 red**. One solid
(`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`, F=1082) misses the 62.5 ms
display window. Everything below is sized against that solid.

---

## 0. Ground rules (non-negotiable)

**Never write to `C:\work\Holosphere` directly.** It is a shared main tree with
concurrent peer sessions committing under the *same git identity* as you. Work
in your own worktree:

```
cd /c/work/Holosphere && BASE=$(git rev-parse refs/heads/master)
git worktree add ../hs-wt-<yourtask> -b work/<yourtask> $BASE
```

Always `git -C C:/work/hs-wt-<yourtask> ...`; never `cd` into a tree and rely on
cwd. Land only by rebasing onto live `refs/heads/master` and doing an **attached
`git merge --ff-only`** in the main tree, asserting
`git merge-base --is-ancestor` first. `master` is append-only. If the merge
refuses over main-tree WIP, **stop and surface it** — never stash or discard a
peer's work. Peers leave WIP in the main tree routinely; leave it alone.

**The device is a single shared Teensy.** Never run `pio run -t upload` or
`profile_capture.py` directly — they bypass the host-global lock, and an upload
issued while a peer holds the port reports SUCCESS *without flashing*, so you can
capture a clean, plausible log of the peer's firmware. The only supported path is
`tools/profile_one.sh`, which takes the lock. `HS_DEVICE_WAIT=<sec>` to queue.
`bash tools/device_lock.sh status` to check.

`tools/profile_one.sh` always builds **`/c/work/Holosphere`** (hardcoded `cd`),
whatever tree you invoke it from. **You cannot profile your worktree — land
first.** For a change that is bit-exact by construction this is safe; for
anything else, prove correctness before landing, then profile, then revert if the
device disagrees (that is exactly what happened to the interleaving lever below).

**Style:** invoke the `code-style` skill first and obey it. Terse factual
comments only — no narration, no justifying correct code, no history, no
finding/ticket references. `core/render/sdf.h` and `core/color/color.h` carry
pre-existing whole-file clang-format drift against local clang-format v22:
hand-format your own lines and commit with `HS_SKIP_FORMAT=1`; **never** run
`clang-format -i`. No `Co-Authored-By` line.

**Gates after every commit:**
- `export EMSDK=C:/work/emsdk; cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF;
  cmake --build --preset tests -j 8; ctest --preset tests` → **51/51**
- `pio run -e phantasm` → `[teensy-gate] phantasm: PASS`, and report the RAM1
  `code` figure **and the per-commit delta**.

**ITCM budget: 190,424 B of 196,608, headroom 6,184 B.** This has vetoed changes
before. If a high-value change is blocked by the ceiling, the owner has
pre-authorized trading away low-value `HS_O3` regions. Measured de-O3 yields
(per-file, by deleting that file's `HS_O3_FN` tokens):

| file | frees | impact |
|---|--:|---|
| `core/color/color.h` | 6,048 B | high — MeshFeedback is colour-bound |
| `core/engine/styles.h` | 2,032 B | low — data-only by convention |
| `core/render/filter.h` | 880 B | medium |
| `core/math/3dmath.h` | 656 B | medium |
| `core/vendor/FastNoiseLite.h` | 464 B | low |
| `core/math/geometry.h` | 48 B | negligible |
| `core/render/shading.h`, `core/render/plot.h` | 0 B | attribute changes no codegen |

~4 KB is available from `styles.h` + `filter.h` + `3dmath.h` + `FastNoiseLite.h`
without touching the colour path. Verify each on device before keeping it.

---

## 1. Where things stand — measured, tip `5b1cbdb7`

Capture: `bash tools/profile_one.sh IslamicStars profile 210 16 -D
HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4`, read from the worst
solid's 8 clean-hold windows. Full report:
`docs/profiles/shipping/profile_islamicstars_teensy_2026-07-19.md`.

```
frame                    124.38 ms   (8 fps: render 67.09, peak 75.1)
  is_timeline_step        67.09 ms
      is_mesh_scan        61.91 ms
        scan_mesh_raster  48.78 ms   79% of scan
          raster_scan     47.83 ms
            raster_shade  13.67 ms     <- shade stage, TASK D
              filter_blend 3.66 ms
          (probe = raster_scan - raster_shade = 34.16 ms, 51% of render)
        scan_face_setup   12.99 ms   21% of scan
          face_bounds      4.68 ms
          face_project     3.20 ms
          face_azimuth     1.39 ms
          face_thetas      1.17 ms
          face_phi_extent  0.49 ms
          face_pole/edges/sectors  0.81 ms
      is_mesh_transform    4.39 ms   <- never examined
      is_face_offsets      0.74 ms
  is_buffer_wait          57.28 ms   (round-up idle, not work)

ISR (concurrent, absorbed by every foreground scope):
  isr_wake        10.35 ms/frame   8.32% CPU   <- TASK C
  isr_pack         8.24 ms/frame   6.62% CPU   <- TASK C
  isr_dma_submit   0.28 ms/frame   0.22% CPU
  total           18.87 ms/frame  15.16% CPU
```

**The budget that matters.** One 62.5 ms window at 15.16% ISR leaves **~53.0 ms**
of foreground render. The worst shape needs **67.1 → 53.0 ms, a 1.27× speedup**.
Note this can be attacked from *either* end: cutting render cost, or cutting ISR
load (which raises the foreground budget). Nobody has tried the second.

Static profile of the hot symbols in the shipping `phantasm` image:

| symbol | insns | vmrs | vdiv | vsqrt | vcvt | vldr[sp |
|---|--:|--:|--:|--:|--:|--:|
| `SDF::Face::Face` | 1263 | **59** | **14** | **9** | 3 | 0 |
| `Scan::rasterize_face` | 997 | 29 | 4 | 3 | 11 | 1 |
| `Scan::Mesh::draw` | 396 | 5 | 3 | 3 | 3 | 0 |
| `POVSegmented::flywheel_isr` | 232 | 0 | 0 | 0 | 0 | 0 |
| `HD107SFrame<72>::packPixel` | 53 | 0 | 0 | 0 | 0 | 0 |
| `DMALEDController::submitFrame` | 77 | 0 | 0 | 0 | 0 | 0 |

Regenerate with:
```
OD=~/.platformio/packages/toolchain-gccarmnoneeabi-teensy/bin/arm-none-eabi-objdump.exe
"$OD" -d .pio/build/phantasm/firmware.elf > ph.dis
```

---

## 2. Methodology — what actually predicts device here

Read this section twice. Prior investigations burned many hours acting on
numbers that did not transfer.

- **Static instruction/`vmrs` counts over-predict, weighted by trip count.** A
  −21.8% instruction / −31% `vmrs` cut to the `Face` ctor bought **−0.43 ms**,
  while removing *two FP ops per edge* from the probe loop bought **−1.68 ms**.
  The ctor runs 1082×/frame; the probe loop runs ~43,500×. Always multiply a
  static delta by its trip count before believing it.
- **Host `-O2` does not predict device `-Os`+regions.** Host CSEs sqrt work GCC
  does not, so a redundancy lever can look dead on host and be alive on device
  (and vice versa). Never report a speedup from host timing.
- **Run-to-run noise on `scan_mesh_raster` is ~0.3 ms (0.6%).** Do not trust a
  single-capture delta smaller than that. Attribute changes by *scope* — changes
  in different `HS_PROFILE` scopes A/B cleanly in one capture.
- **Read per-shape numbers only from clean-hold windows** (`scan_mesh_raster`
  calls == 16 × F, F = `face_counts.size()`). Windows straddling a transition
  draw two meshes and mis-state everything. `parse_profile.py <log> presets
  --gate scan_mesh_raster` does this for you.
- **Always `python tools/parse_profile.py <log> validate` before trusting a
  capture, and check the log's mtime** — `profile_one.sh`'s verify has been seen
  to pass on a stale log when the serial capture flaked.
- **A spill is not always register pressure.** See TASK A.

---

## 3. TASK A — register pressure

**Status: the pixel loop is solved; the ctor is clean; the remaining question is
elsewhere.** Do not redo the solved parts.

What was found and fixed: `rasterize_face` was inlined into the enormous
`Scan::Mesh::draw`, and the allocator spilled the per-face constants (`center`,
`basis_u`, `basis_w`) so the per-pixel projection **reloaded them from the stack
every pixel** — 22 FP stack reloads in the pixel loop. Forcing `rasterize_face`
out of line (`5b1cbdb7`, `__attribute__((noinline)) inline`) dropped that to 1 at
**+0 bytes**, worth **−0.80 ms**. It is bit-exact by construction: no FP
expression crosses the call boundary.

**The transferable lesson: before concluding a hot loop is register-bound, check
whether it is welded into a mega-function.** The slack may already exist and be
wasted. The diagnostic is `vldr sN, [sp` density inside the loop body.

Where that leaves you:

- `SDF::Face::Face` is now its own symbol with **0 stack reloads** — it is *not*
  register-bound. Its cost is FPU stalls, not spills (TASK B).
- `rasterize_face` has **1** stack reload. Effectively clean.
- **Unexamined:** whether any *other* hot path in this effect is still welded
  into a caller and silently spilling. Sweep `vldr sN, [sp` density per symbol
  across the image and rank by known hotness. `Filter::Pixel::Feedback::flush`
  (1840 insns, 19 sp-reloads) and `Scan::Volume::draw` (760 insns, 19
  sp-reloads) are the two visibly spilling symbols in the image, but neither is
  on IslamicStars' path — they belong to MeshFeedback and Raymarch. If you find a
  spilling symbol that *is* on this effect's path, the de-inline lever above is
  the first thing to try, and it is nearly free.

---

## 4. TASK B — FPU pipeline stalls in the render path

Cortex-M7 stall costs: `vmrs APSR_nzcv, fpscr` ~5 cyc (FPU→core flag transfer,
not pipelined); `vdiv`/`vsqrt` ~14 cyc and **not pipelined** — they block the FPU
pipe; `vcvt` ~2; core↔FPU `vmov` ~2.

### B1 — `SDF::Face::Face`: 59 `vmrs`, 14 `vdiv`, 9 `vsqrt` in 1263 instructions

Runs **1082×/frame**, costs **12.99 ms/frame (21% of scan)**. The largest single
concentration of stall content on this effect's path. A prior round already took
it 1718 → ~1300 instructions and 96 → 59 `vmrs` for −1.61 ms, so the cheap
`__builtin_fminf`/`fmaxf` and modulo-elimination passes are **done** — there are
zero runtime `% count` sites left in `sdf.h`. What remains is the divide/sqrt
density: 14 `vdiv` + 9 `vsqrt` ≈ 320 cycles of non-pipelined stall per face, ≈
0.58 ms/frame at 1082 faces.

These come from repeated `normalize()`/`normalized()` calls and `1.0f/x` in
`setup_frame_and_polygon` and `compute_full_bounds`. **Sharing one reciprocal
across normalizations is NOT bit-exact and is forbidden** unless you can prove
otherwise (you almost certainly cannot — a prior attempt measured 47.3% of
5.71 M planes differing in bits). Report the cost as characterized rather than
forcing it. The honest question to answer is whether any of these divides are
*redundant* (same value computed twice) rather than *inherent* — a redundancy is
removable bit-exactly, a shared reciprocal is not.

### B2 — the per-probe `vdiv` and `vsqrt` (highest trip count on the effect)

Per probe, on ~43,500 probes/frame:
- `float inv_cos = 1.0f / cos_angle;` — one `vdiv`, **100% of probes**.
- `sqrtf(d)` at the tail of `plane_dist_exact`/`plane_dist_sector` — one
  `vsqrt`, **~89% of probes** (every non-convex probe).

At ~14 cyc each non-pipelined, that is ~28 cyc of the ~460 cyc probe — **~6%**.

**The `vsqrt` has a real bit-exact escape and nobody has tried it.** The caller
rejects on `d >= pixel_width` where `res.dist = raw - thickness` and, for the
98.9% of faces with `linear_dist`, `raw = ±sqrt(dsq)`. For an *outside* probe the
reject test `sqrt(dsq) - thickness >= pixel_width` is equivalent to
`dsq >= (pixel_width + thickness)²`, and the threshold is loop-invariant. So the
square root is only needed for probes that **survive** the reject — about
**42%** of probes (18,100 shade candidates of ~43,500). Skipping it on the other
58% saves ~8 cyc/probe ≈ 1.5–2% of the probe stage.

This requires a reject-first API change (`distance()` currently computes
everything then hands `dist` to the caller). It is bit-exact for survivors —
they still take the same `sqrtf` — but you must handle: the *inside* case (sign
negative, never rejects), the `!linear_dist` faces (1.1%, where `raw =
fast_atan2(plane_dist, 1)` — `fast_atan2` is monotonic in its first argument, so
the same threshold trick works but the constant must be transformed correctly, or
just fall back to the current path for those faces), and `thickness` being 0 in
the shipping config. **Prove it with the framebuffer oracle before believing it.**

### B3 — the remaining 29 `vmrs` in `rasterize_face`

Disassemble and attribute each. The exact-walk crossing test was already reduced
to two strict compares per edge; the sector search was moved to integer keys
(`angle_key`), removing ~5 per sector probe. What is left is likely the AA/alpha
path and the reject compares. Note **a `vcmp`+`vmrs` feeding control flow cannot
be removed by `vsel`** — `vsel` reads the same APSR flags. The only real escapes
are (a) integer-domain comparison via an order-preserving key (the transform in
`angle_key` is the correct one — `(u & 0x80000000) ? (0u - u) : (u + 0x80000000)`
maps ±0.0 to the same key and survives device FTZ), or (b) `__builtin_fminf`/
`fmaxf` where the result is a value rather than a branch.

### B4 — 11 `vcvt` in `rasterize_face`

Typically `static_cast<int>(floorf(x))` idioms or fixed-point packing. Cheap
individually (~2 cyc) but worth a look given the trip count.

---

## 5. TASK C — ISR prep and processing (the biggest unexamined lever)

**18.87 ms/frame, 15.16% of CPU, and nobody has ever optimised it.** Cutting it
raises the foreground render budget directly: at 15.16% the budget is 53.0 ms and
the worst shape needs 1.27×; halve the ISR load and the budget becomes ~57.8 ms
and the requirement drops to ~1.16×.

### C1 — `isr_pack` is load-bound, not FP-bound: 8.24 ms/frame

The critical measured fact: **the ISR path contains zero floating-point stalls**
— `flywheel_isr` (232 insns), `HD107SFrame<72>::packPixel` (53 insns) and
`submitFrame` (77 insns) have **0 `vmrs`, 0 `vdiv`, 0 `vcvt`** between them. Do
not go looking for FPU stalls here; it is a **memory** problem.

`isr_pack` fires ~287×/frame (once per column) at **28.74 µs average = ~17,240
cycles per call**. With 72 LEDs per column that is roughly **240 cycles per
pixel** against a 53-instruction `packPixel` — an effective **CPI of ~4.5**.
That gap is load latency, and the likely cause is in `core/engine/memory.cpp`:
the framebuffers (`buffer_a`/`buffer_b`) **carry `DMAMEM`, so they live in OCRAM,
not the zero-wait DTCM** the arenas use. The DMA engine cannot reach DTCM, so
they must be there — but the *pack* reads them with the CPU.

Leads, in order of suspicion:

1. **Access pattern.** Confirm what stride `packPixel` walks the framebuffer in.
   A column-major read across a row-major buffer is a worst case for the OCRAM
   burst/cache behaviour. If the pack reads a column while the framebuffer is
   laid out by row, each pixel is a separate cache line touch. Reordering the
   framebuffer, or packing in the buffer's natural order, could be a large win.
2. **Cache state.** Check whether the framebuffer region is cacheable and whether
   the D-cache is being invalidated/cleaned around the DMA handoff. An
   unnecessary invalidate per column would explain a lot; so would a
   non-cacheable mapping. Look at how the buffers are declared and whether
   `arm_dcache_*` calls bracket the pack.
3. **Wider loads.** 21 `ldr` in a 53-instruction `packPixel` — if it is loading
   byte/halfword components separately, `ldrd`/`ldm` or a 32-bit load with shifts
   could cut the load count.
4. **Do less.** Is the pack converting a format it need not? `Pixel16` → HD107S
   frame conversion may be foldable into how the renderer writes pixels.

**Measure before optimising.** `HS_ISR_PROFILE` scopes already exist
(`core/engine/platform.h`, used in `pov_segmented.h`); add finer ones inside the
pack if needed — but note ISR-side profiling must use `HS_ISR_PROFILE`/
`hs::IsrCycleStats`, **not** `HS_PROFILE` (the main-loop registry is non-atomic).

### C2 — `isr_wake`: 10.35 ms/frame, and the average is skewed

`n≈2290/frame`, avg 4.52 µs, but **min 0.63 µs and max 199 µs**. The mean is
dragged by outliers, so the typical wake is far cheaper than 4.52 µs and the
total is dominated by a tail. **Characterise the distribution before optimising
the body** — if the cost is a few long preemption/contention events rather than
per-call work, the fix is scheduling, not the 232 instructions in
`flywheel_isr`. Look for what blocks or preempts during those spikes (the pack
ISR, the DMA completion, or foreground code disabling interrupts).

### C3 — the honest caveat

ISR cost scales with display cadence, not with effect cost: at 8 fps a frame is
125 ms and sees ~2× the ISR fires of a 62.5 ms frame. So per-*window* ISR load is
roughly constant, and the percentage figures move as render time moves. Quote
ISR cost in **µs/window or CPU%**, and state which frame length you measured at.

---

## 6. TASK D — the shade stage: 13.67 ms/frame, never analysed

`raster_shade` is **20% of render** at 18,100 events/frame (1.75 per quadrant
pixel — genuine overdraw from overlapping faces and AA fringes), **453 cyc per
event**: ~121 cyc in `filter_blend` and **~332 cyc in the fragment shader**.

The shader is `shade_mesh_topology` (see `effects/IslamicStars.h`) and it has
**never been examined for stalls**. 332 cyc × 18,100 = ~10 ms/frame. Start by
finding its symbol in the disassembly (it may be inlined into `rasterize_face`,
which carries 11 `vcvt` and 29 `vmrs`), and apply the TASK B patterns. Given the
trip count is the highest of any scope on this effect, a per-event saving here
multiplies harder than anywhere else — this is the same reasoning that made the
2-FP-op probe-loop change worth −1.68 ms.

Also unexamined: **`is_mesh_transform`, 4.39 ms/frame**, transforming 3240
vertices per frame. Straightforward vector work; check it for the same
divide/normalise density as the `Face` ctor.

---

## 7. TASK E — the remaining red is a transition cost, not steady state

The worst solid holds **61.91 ms** of scan against the 62.5 ms window — its
steady render is 67.1 ms, 4.6 ms over — but its **bucket peak is 75.1 ms**.
Buckets charge each preset the transition that follows it, where **two meshes are
drawn in one frame**.

So a large part of what keeps this effect red is the shape transition, not the
held shape. The pattern that fixed the analogous problem in MeshFeedback is a
**drain transition** (`5f975feb`: emit → drain-to-black → next preset), which
took that effect's peak 180.8 → 87.7 ms. Whether an equivalent is acceptable here
is an *art* decision — it changes how shapes hand over — so **surface it to the
owner rather than implementing it unilaterally.**

If steady-state work continues instead, note that the remaining margin is small:
4.6 ms of 67.1, i.e. ~7%.

---

## 8. Measured dead — do not re-attempt

Each of these cost real time to close. They are closed.

- **Interleaving two columns through the exact walk** (`814655b6`, reverted in
  `f0916471`). Paired adjacent columns so each edge record loaded once (4
  loads/pixel vs 7) and the FPU had two independent chains. Bit-exact
  (framebuffers byte-identical, 24 solids × 64 orientations, negative control
  confirmed live on 16/24). **Regressed on device**: `is_mesh_scan` 61.9 → 63.2,
  a second shape slipped red, +4,720 B ITCM. The paired-path overhead exceeded
  the latency hidden on the **narrow spans (median 9–15 px)** that dominate this
  solid. The MeshFeedback interleave (`dedd3dde`) works because its composite is
  one long uniform chain over wide runs; this is not that shape. **Do not retry
  pairing on the walks.** If you ever retry it, it must be on a scope with long
  uniform runs — the shade stage, not the probe.
- **Per-row scan spans (L1).** Implemented in full and measured dead: covering
  the AA fringe costs ~7 columns of padding per span end against rows ~20 columns
  wide, so at zero pixel leakage probes go **up 0–3.6%**. The "57–59% saving"
  figure measured lit pixels and ignored the fringe the scan must cover. Split
  intervals lose too. See `docs/single_pass_mesh_raster_spec.md`.
- **Linear scan replacing the sector binary search.** Not bit-exact — K2 faces
  are only *weakly* monotone, so counting keys ≤ q does not reproduce a binary
  search on a non-sorted array (4,277 index divergences in 22.4 M probes) — and
  ~2.1× the dynamic work at mean face count 16.77.
- **`inv_sector_span` reciprocal multiply.** Worth ~3% of the sector stage; costs
  **1 divergent pixel in 112 M probes**. Out under a bit-exact bar.
- **Scratch-buffer placement.** `global_arena_block` links at `0x20000ce0`, size
  `0x4a800`, inside DTCM (`0x20000000`); OCRAM starts at `0x20200000`. Every
  `EdgePacked` load is already zero-wait. (Note this is the *arena*; the
  framebuffers are a different story — see TASK C.)
- **Widening the `linear_dist` threshold**, **sinking `Vector p` below the reject
  test**, **templating the scan on LUT presence**, and **any scheme that reduces
  probe count by tightening per-face scan spans** (four separate attempts).
- **The polar AA pad (`1/sin φ`).** The under-coverage is real (805 px/frame) but
  invisible — the mesh tiles the sphere, so a clipped fringe pixel is painted by
  the neighbouring face; **0 true holes across 192 draws**. The fix costs ~+25%
  probes. Documented in `docs/probe_path_open_items.md`.

**Two specific traps:**

- **A sign-bit test on a float difference is only safe for STRICT comparisons.**
  For `a > b` it is correct (equality gives `+0.0`, sign clear, "not greater" —
  which agrees). For a non-strict `a <= b` it inverts on an exact hit. And
  subnormal differences can flush to zero on device where no host A/B catches it.
  The clean solution is the order-preserving integer key in `angle_key`.
- **`hs::clamp` already emits `vminnm`/`vmaxnm`.** Do not "fix" it.

---

## 9. Correctness bar and the proof harnesses

**Output must be bit-exact.** Prove it, do not assert it. Anything that changes a
float result needs a much stronger argument than a speedup. Three oracles exist
and are worth reusing — all are ad-hoc scratch tools, not committed, and are
rebuilt against two checkouts and `cmp`'d:

- **Framebuffer equality** — renders every Islamic-registry solid through
  `Scan::Mesh::draw` over N orientations and FNV-hashes the RGB16 framebuffer.
  This is the right oracle for anything touching the raster or probe path. Build
  it from `tests/aa_audit_main.cpp`'s `render()` pattern.
- **`Face` construction** — placement-news each face into a poisoned buffer and
  hashes every output field plus all 12 scratch spans. The right oracle for
  anything touching the ctor. A prior round ran 9,646,000 faces byte-identical.
- **Per-probe dump** — dumps `(sector index, returned float)` per probe. A prior
  round ran 112,068,276 probes byte-identical.

**Always include a negative control.** A matching hash is worthless if the code
path was never executed. Perturb your new path by an epsilon and confirm the hash
*moves*, and report how many solids it moved (a prior pairing proof moved 16/24 —
the other 8 legitimately never take that path). One committed harness exists:
`tests/aa_audit_main.cpp` (`-DHS_AA_AUDIT`, zero-cost when the macro is absent).

---

## 10. Reporting standard

For each lever: what you measured and the exact command, what you changed, the
before/after disassembly evidence (instructions and `vmrs`/`vdiv`/`vsqrt`/`vcvt`/
`vldr[sp` for the affected symbol), the bit-exactness proof with sample counts
and the `cmp` result, `ctest`, and phantasm RAM1 + headroom + per-commit delta.

**State explicitly if a lead came back empty.** A well-measured negative is a good
outcome here — most of this problem's dead ends were only closed because someone
reported one honestly, and section 8 exists entirely because of that. Do not
report a speedup from host timing alone. If the device disagrees with your static
prediction, say so and revert; that has already happened twice and both reverts
were the right call.
