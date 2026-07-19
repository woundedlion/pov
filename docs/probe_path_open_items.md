# Brief: IslamicStars probe path — two open items

You are picking up work on `c:\work\Holosphere` (Holosphere / Phantasm POV
renderer, C++20, Teensy 4.0 @ 600 MHz). Assume none of the prior session's
context. Read this whole brief before touching anything.

---

## Status

**Task A — open.** The walks were re-measured on device after the three earlier
fixes: buckets 469.6 cyc/probe, exact walk 387.4 cyc/event = **64.6 cyc/edge**,
so the gap to the ~20-25 target is real. Two leads are now closed:

- **Scratch-buffer placement (was lead 1): DEAD.** `global_arena_block` links at
  `0x20000ce0`, size `0x4a800`, inside DTCM (`0x20000000`); OCRAM starts at
  `0x20200000`. `EdgePacked` loads are already zero-wait.
- **Sector-walk modulo: FIXED.** `packed_edges[(s + k + count) % count]` lowered
  to `sdiv` + `mls` feeding the address the five edge loads depend on. Replaced
  by one wrap correction; sector stage 434.5 -> 415.0 cyc/event (-4.5%), every
  other stage flat, worst solid's scan 70.37 -> 69.84 ms.

Still open: the exact walk's parity test carries **3 `vmrs` per edge** (`ylo <=
py`, `py < yhi`, `px < isx`), up to 15 of its 64.6 cycles. A branchless form
still needs APSR flags — `vsel` reads them too — and the arithmetic dodge is the
sign-bit trick ruled out below. `SDF::Face::Face` (96 `vmrs`, 17 `vdiv`, 14
`vsqrt`, runs 1082x/frame) is unexamined and is the largest unclaimed target.

**Task B — closed, measured, not fixed.** The `1/sin(phi)` under-coverage is
real: 805.5 missed fringe px/frame registry-mean, 1,223.9 on the worst solid,
confined to rows 0-35 and 108-143, max alpha 0.50 (mean 0.084). It is also
invisible, for a structural reason: **the mesh tiles the sphere**, so a fringe
pixel clipped off one face is painted opaque by its neighbour — **0 true holes
across 192 draws under both a flat and a per-face-primary shader**. The
conservative pad costs **+25.6% probes whole-canvas, +50.8% on the polar rows**
(~+16.6 ms/frame on a polar segment at 520 cyc/probe), against an effect already
missing cadence on 3 shapes. Not worth it. The audit harness is landed
(`-DHS_AA_AUDIT`, `tests/aa_audit_main.cpp`); the pad change is not.

Residual and independent of the pad: four solids still miss fringe after a
corrected pad (`truncatedIcosidodecahedron_truncate50d_ambo_dual` 486.5/frame,
`truncatedIcosahedron_ambo_relax100_hk54_needle` 199.2,
`snubDodecahedron_truncate5d_ambo_dual` 158.2,
`truncatedIcosahedron_truncate50d_ambo_dual` 141.0). That defect is in the
azimuth-interval construction, not the pad. Same zero-holes verdict.

**Task C — audit regenerated; one fix landed.** The snapshot table below is not
comparable to a fresh run (different classifier basis, not image drift).
`gamut_channel_exit`, the snapshot's top-ranked row, is **unreachable in
MeshFeedback**: it is called only from `gamut_clip_analytic`, taken only when the
LUT is disarmed, and `init_gamut_lut`'s one caller repo-wide is
`MeshFeedback.h:129`. `gamut_bracket_refine` was the confirmed-hot row and went
48 -> 16 `vmrs`.

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
`git merge-base --is-ancestor` first. `master` is append-only and a
`reference-transaction` hook enforces it. If the merge refuses over main-tree
WIP, **stop and surface it** — never stash or discard a peer's work.

**The device is a single shared Teensy.** Never run `pio run -t upload` or
`profile_capture.py` directly — they bypass the host-global lock and an upload
issued while a peer holds the port reports SUCCESS *without flashing*, so you
can capture a clean, plausible log of the peer's firmware. The only supported
path is `tools/profile_one.sh`, which takes the lock. `HS_DEVICE_WAIT=<sec>` to
queue. `bash tools/device_lock.sh status` to check.

`tools/profile_one.sh` always builds **`/c/work/Holosphere`** (hardcoded `cd`),
whatever tree you invoke it from. **You cannot profile your worktree.** Land
first, or hand the run back to the coordinator.

**Style:** invoke the `code-style` skill first and obey it. Terse factual
comments only — no narration, no justifying correct code, no history in
comments, no finding/ticket references. `core/render/sdf.h` has pre-existing
whole-file clang-format drift against local clang-format v22: hand-format your
own lines and commit with `HS_SKIP_FORMAT=1`; **never** run `clang-format -i`.
No `Co-Authored-By` line.

**Gates after every commit:**
- `export EMSDK=C:/work/emsdk; cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF;
  cmake --build --preset tests -j 8; ctest --preset tests` → **51/51**
- `pio run -e phantasm` → `[teensy-gate] phantasm: PASS`, and **report the RAM1
  `code` figure and headroom**. Headroom is ~6 KB of a 196,608 B ceiling; this
  is a real constraint that has vetoed changes before.

---

## 1. Background — how the current numbers were obtained

IslamicStars rasterizes a polyhedral mesh as per-face SDFs
(`Scan::Mesh::draw` → `rasterize_face` in `core/render/scan.h`, `SDF::Face` in
`core/render/sdf.h`). Its heaviest solids miss the 62.5 ms display window.

Four separate schemes to reduce the **probe count** (pixels evaluated per face)
were built and all failed. The lever turned out to be **per-probe cost**, which
had never been measured on hardware because the scan counters were not compiled
into the profile image. Two instrumentation flags now exist — use them, do not
rebuild this capability:

- **`-D HS_SCAN_METRICS`** — per-probe counters (`pixels_tested`, `culled`,
  path split convex/sector/exact-walk/lut). Read with
  `python tools/parse_profile.py <log> metrics`.
- **`-D HS_PROBE_BREAKDOWN`** — per-probe *stage* cycle buckets (point,
  project, convex, sector, exact, pack, alpha). Read with
  `python tools/parse_profile.py <log> probe`. It self-measures the CYCCNT read
  cost (~1.06 cyc, negligible) and subtracts it — the `net` column is what you
  want. Both are **zero-cost when the flag is absent** (verified byte-identical
  `teensy_size`).

Capture line (cycler; the epoch stretch and trans-speed knobs change how long a
shape *holds*, never its per-frame cost):

```
bash tools/profile_one.sh IslamicStars profile 210 16 \
  -D HS_PROBE_BREAKDOWN -D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4
```

Then **always** `python tools/parse_profile.py <log> validate` before trusting a
capture, and check the log's **mtime** — `profile_one.sh`'s verify has been seen
to pass on a stale log when the serial capture flaked.

**Read per-shape numbers only from clean-hold windows**: a window is clean for a
shape when `scan_mesh_raster` calls == 16 × F, where F = `face_counts.size()`.
Windows straddling a shape transition draw two meshes and mis-state everything.

### Current measured state (worst solid: `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`, F=1082)

Pre-optimisation baseline, from `HS_PROBE_BREAKDOWN`:

| stage | events/frame | cyc/event | cyc/probe |
|---|--:|--:|--:|
| exact walk | 19,180 | 463.6 | 204.3 |
| sector walk | 18,337 | 465.4 | 196.1 |
| project | 43,516 | 29.5 | 29.5 |
| point | 43,516 | 28.4 | 28.4 |
| convex | 4,960 | 216.5 | 24.7 |
| pack | 43,516 | 18.0 | 18.0 |
| alpha | 19,038 | 35.9 | 15.7 |
| **buckets** | | | **516.7** |
| **measured** | | | **583.6** |

(The 66.9 residual is `Vector p` construction, the mask check, the reject
compare and loop overhead — deliberately outside the buckets.)

Path mix: convex 11.4% / sector 42.1% / exact 44.1% / lut 0%. Diagnosis was
**stalls, not bloat**: 246 instructions/probe at **2.33 CPI**, ~55–60% stall,
concentrated in ~19 `vmrs APSR_nzcv, fpscr` per probe (the Cortex-M7 FPU→core
flag transfer, ~5 cycles each, not pipelined).

Three fixes then landed (`d59384d6`, `ab8bad4e`, `6cab193d`), taking the probe
to **520.5 cyc (−10.8%)**, verified bit-exact over 57.9 M probes:
`__builtin_fminf` so GCC emits `vminnm`; a min/max bracket replacing two
`vcmpe`+`vmrs` bool materialisations per edge; and deleting a stale `noinline`
on the three `plane_dist_*` routines (its ~40%-spill-cliff comment did not
reproduce — phantasm actually *shrank* 80 B).

Roster went 18 green / 6 red → **21 green / 0 yellow / 3 red**, worst peak
89.4 → 80.8 ms. Some of that is concurrent peer work (a color/gamut rewrite and
hankin effect changes), so do not treat the whole delta as attributable.

---

## TASK A — the edge walk is still ~77 cycles per edge

Both walks cost ~465 cyc/event for ~5–6 edges (the exact walk covers all 6 on
the degree-6 faces that dominate this solid; the sector walk covers ±2 = 5).
That is **~77 cycles per edge** — pre-fix; post-fix it should be nearer ~68,
but **nobody has re-measured**.

Per edge the arithmetic is a point-to-segment distance in the face's 2D
gnomonic plane, from `FaceScratchBuffer::EdgePacked {vx, vy, ex, ey,
inv_len_sq, inv_ej, next_vy, pad}`:

```
t   = clamp(((px-vx)*ex + (py-vy)*ey) * inv_len_sq, 0, 1)
dx  = px - vx - t*ex;  dy = py - vy - t*ey
dsq = dx*dx + dy*dy
d   = min(d, dsq)                       // now vminnm
plus a crossing/parity test on vy / next_vy
```

~15 flops and ~7 loads. On an M7 with a pipelined FPU that should be ~20–25
cycles. **There is still a ~3× gap. Find it.**

### Step A1 — re-measure first (mandatory)

Every number above predates the three landed fixes. Re-run the breakdown on
current master and report the post-fix per-stage table. If the walks are
already near 25 cyc/edge, stop and say so — do not optimise a number you have
not re-confirmed. **This investigation has been burned repeatedly by acting on
stale or host-simulated figures.**

### Step A2 — leads, in my order of suspicion

1. **Where does `FaceScratchBuffer` actually live?** It is allocated from a
   scratch arena in `Scan::Mesh::draw`. If that arena is in OCRAM/DMAMEM rather
   than DTCM, every `EdgePacked` load is a slow off-chip access and 6 loads per
   edge with immediate `vfma` dependencies would explain the entire gap. **Check
   the linker map / section placement.** This is my strongest lead and it is
   cheap to test. See `docs/` for the DTCM bank ledger and note that
   `DMAMEM` has previously been silently dropped on template statics in this
   codebase — a variable can end up somewhere you did not intend.
2. **Load-use stalls.** 6 `vldr` per edge feeding `vfma` immediately. Try
   restructuring `EdgePacked` (it is already 32 B / 8 floats), software
   pipelining the next edge's loads, or hoisting `px`/`py`-invariant terms.
3. **Remaining `vmrs`.** The three landed fixes cut ~5.6 per probe of ~19.
   Disassemble and count what is left; each is ~5 cycles.
4. **Register pressure.** A prior attempt to hoist a whole-`Face` local copy
   made things *worse* (sp-loads 33→49, +208 B) — the loop is already
   register-bound with a 588 B frame and 53 spill stores. Do not simply retry
   that; if you attack pressure, do it by *shrinking* live state.
5. **`clamp` and the parity test** — check what they emit now.

### Step A3 — method

Disassemble with the PlatformIO `arm-none-eabi-objdump` against the `profile`
ELF, and report **instructions/probe and `vmrs`/probe** for every variant. That
static metric predicted the last round's device result accurately, so use it to
triage before spending device runs.

**Correctness bar:** the geometry is contract-tested and the output must be
**bit-exact**. Prove it, do not assert it: build a host sweep that renders the
real registry solids and compares framebuffers (or dumps `dist`/`raw_dist`/
`size` per probe) before and after. The prior round did 57.9 M probes,
byte-identical. Anything that changes a float result needs a much stronger
argument than a speedup.

**Do not re-attempt (measured dead):** widening the `linear_dist` threshold to
avoid the atan (it is 3 cyc/probe amortised — `linear_dist = size < 0.2f`
already covers 98.9% of faces, and GCC folds the call because `x` is the
literal `1.0f`); sinking `Vector p` below the reject test (a wash); templating
the scan on LUT presence to kill the dead `lut_data` test (duplicates the scan
body against ~6 KB of headroom); and **any** scheme that reduces probe *count*
by tightening per-face scan spans — four of those have failed, see
`docs/single_pass_mesh_raster_spec.md` for why the family is closed.

**A sign-bit XOR on `(vy - py)` for the parity test is WRONG** — do not
"optimise" it that way. `x - x` yields `+0.0`, whose sign bit is *clear* while
`vy > py` is false, so the test inverts on an exact hit; and a subnormal
difference can flush to zero on device where no host A/B would catch it.

---

## TASK B — latent polar AA under-coverage (correctness)

**The bug.** `SDF::Face::get_horizontal_intervals` (core/render/sdf.h) pads its
emitted azimuth interval by a **constant**:

```
float pad = 1.25f * (2.0f * PI_F / W);   // radians of theta
```

`clip_rejects` uses the same constant. But the renderer paints every pixel with
`d < pixel_width`, where `d` is an (approximately angular) **distance**. At
polar angle φ, an angular distance `R` subtends `Δθ = R / sin(φ)` of azimuth.
So covering a one-pixel fringe requires `1.25 * (2π/W) / sin(φ)` — and the pad
is short by exactly `1/sin(φ)` wherever the face is away from the equator.

**Evidence it is real.** A per-row span implementation (since reverted) painted
**1–19 more pixels/frame** than the constant-pad extent path — i.e. the extent
path is clipping AA fringe that the geometry says should be painted. Separately,
a probe measured the true azimuth reach of the fringe as **18.98 columns per
side at row 3**, 4.80 at row 12, 2.48 at row 24, 1.76 at row 36 — the 1.25
figure only holds near the equator.

**Impact is small**: ~0.1% of shade events, concentrated in the polar caps,
losing the outermost AA column. Nobody has confirmed it is *visible*.

### What to do

1. **Quantify before fixing.** Render the real solids natively and diff the
   framebuffer against a correctly-padded build. Where are the missing pixels,
   how many, what alpha did they carry, and is the result visible (e.g. edge
   hardening or a seam near the poles)? There is a framebuffer-dump harness
   pattern in the repo — render to PNG rather than theorising.
2. **Then price the fix.** The obvious correction — pad by
   `1.25 * (2π/W) / sin(φ)` — is not free: near the pole it widens every face's
   scan span substantially, and at **~520 cycles per probe** extra probes are
   expensive on an effect that already has 3 shapes missing cadence. Note
   `get_horizontal_intervals` currently **ignores its row argument** (the
   parameter is unnamed) and `rasterize_face` evaluates it once at `y_lo`, so a
   per-row pad is not directly available; the cheap conservative form is the
   worst-case `1/sin(φ)` over the face's `[y_min, y_max]` band.
3. **Report the trade and recommend** — do not land a correctness fix with a
   material frame-time cost without surfacing the number. If it is visually
   undetectable and costs several ms, the right answer may be to leave it and
   document it.

Whatever you conclude, **add a regression test** if you fix it: the repo's
convention is that every new test is wired into the runner
(`tests/test_sdf.h` / `tests/test_mesh_raster.h`, registered in the
corresponding `run_*_tests()`; the modules are already in `_hs_test_modules`).
`tests/test_sdf.h` already has `expect_face_cull_covers_fringe`, which
brute-forces "every pixel with `distance < pixel_width` must be visited" — but
only over synthesized *convex* polygons near the equator, which is precisely why
this slipped through. A polar case belongs there.

---

## Reporting

For each task: what you measured (with the command), what you changed, the
disassembly evidence, the bit-exactness proof, ctest, and phantasm RAM1 +
headroom. State explicitly if a lead came back empty — a well-measured negative
is a good outcome here and several of this problem's dead ends were only closed
because someone reported one honestly. Do not report a speedup from host timing
alone; host `-O2` has repeatedly failed to predict device `-Os`+regions on this
code.

---

## TASK C — repo-wide FPU stall audit

The `vmrs` finding in Task A is **not** IslamicStars-specific. A static audit of
the shipping `phantasm` image found **1,467 `vmrs`, 560 `vdiv`, 548 `vcvt`, 583
core<->FPU `vmov`, 181 `vsqrt`** across 262 functions.

### The audit script

Regenerate this after every change rather than working from the snapshot below:

```bash
OD=~/.platformio/packages/toolchain-gccarmnoneeabi-teensy/bin/arm-none-eabi-objdump.exe
"$OD" -d .pio/build/phantasm/firmware.elf > ph.dis
# then, per function, count: vmrs (*5 cyc), vdiv/vsqrt (*14, NOT pipelined),
# vcvt (*2), core<->FPU vmov (*2), and rank by the weighted sum.
```

### Snapshot, ranked by static stall content x known hotness

| function | insns | vmrs | vdiv | sqrt | vcvt | why it matters |
|---|--:|--:|--:|--:|--:|---|
| `SDF::Face::Face` | 1718 | **96** | 19 | 14 | 3 | runs **1082x/frame** on the worst solid; nobody has looked at the ctor |
| `DisplacementField::draw_rings` | 2576 | 85 | 20 | 7 | 37 | green, only ~3 ms margin |
| `gamut_channel_exit` | 855 | **102** | 21 | 1 | 0 | **highest vmrs density in the image**; color path, every blending effect |
| `Scan::Volume::draw` | 760 | 51 | 12 | 15 | 3 | Raymarch — 100 % spill, worst effect on the roster |
| `Plot::rasterize` (5 variants) | ~1000 | 48/33/27/25 | ~13 | 6 | — | MindSplatter, DreamBalls, Comets |
| `HopfFibration::render_trails` | 1366 | 64 | 10 | 6 | 10 | green, thin margin |
| `Filter::Pixel::Feedback::flush` | 794 | 19 | 7 | 2 | 6 | **MeshFeedback's dominant scope** (red) |
| `BZReactionDiffusion::shade_pixel` | 602 | 15 | 2 | 0 | **39** | its GS twin is the worst effect on the roster |
| `slerp` / `make_basis` / `quaternion_from_basis` / `Motion::path_frame` | ~150-200 | 4-9 | **11-15** | 4-5 | 0 | brutal divide density for their size |

### Three patterns, in value order

1. **`vmrs` clusters.** `if (a < b) b = a` on floats emits `vcmp` + `vmrs
   APSR_nzcv, fpscr` + a predicated `vmov`; `__builtin_fminf`/`fmaxf` emit
   `vminnm`/`vmaxnm` instead. This is the fix already proven in `d59384d6` — it
   is pure missed optimisation, GCC picks `vminnm` for `hs::clamp` two lines
   away. Bit-exact for non-NaN inputs; confirm no NaN reaches the site and say
   why in the report, not in a comment.
2. **Divide clusters in small math helpers.** `slerp` carries 11 `vdiv` + 4
   `vsqrt` in 201 instructions — mostly stall. Same shape in `make_basis`,
   `quaternion_from_basis`, `path_frame`. These look like repeated
   normalisations that could share one reciprocal, or use an `rsqrt`
   approximation where the precision budget allows. `vdiv`/`vsqrt` are ~14
   cycles and **not pipelined** — they block the FPU pipe.
3. **`vcvt` clusters** in per-pixel code — `BZReactionDiffusion::shade_pixel`
   (39) and the `Feedback::flush` lambdas (29 each), typically from
   `static_cast<int>(floorf(x))` idioms or fixed-point packing.

### The trap

**These are STATIC counts. A `vmrs` in a cold function costs nothing.** The
table is already filtered against the profile reports' dominant scopes, but each
row still needs its execution frequency confirmed before anyone spends effort on
it. Acting on an unverified static or host-side number is exactly what cost the
prior investigation six hours and four reverted schemes. Confirm hotness from
`docs/profiles/shipping/` or a `HS_PROFILE` scope first, then optimise.

Per-effect reports and the ranked roster live in `docs/profiles/` (gitignored —
regenerate with the teensy-profile skill if absent).
