# Agent brief — IslamicStars render chain to peak < 60 ms on every shape

**Status: OPEN brief.** Drive the IslamicStars render chain to peak < 60 ms on
every shape, in both the shipping selective-`-O3` image and the global-`-O3`
image. The branch layout described under "Where the work lives" predates the
opchain_morph landing and no longer tracks the tree.

## Objective

Optimize **IslamicStars** (the opchain_morph version) so that a full 24-shape
on-device cycle profiles with **peak frame render < 60 ms on every shape, 0
spilled frames**, in **both** the shipping selective-O3 image and the global-O3
image. This is the gate to re-landing the reverted opchain_morph feature.

Colour convention (owner, binary): **0 spills = green, any spill = red, no
yellow.** "Green with margin" is the goal here: peak render **< 60 ms** (a
comfortable margin under the 62.5 ms / 16 fps frame budget), holds *and*
transitions.

## Where the work lives

- **Branch:** `opchain/dissolve` (`7fd8155d`), which is master `009d3536` + the
  full opchain_morph feature (77 commits) + both ITCM fixes + the Dissolve seed
  intro. Do **not** touch master; master was deliberately rewound to the green
  pre-opchain tip. Work on `opchain/dissolve` or a child branch off it.
- **Effect:** `effects/IslamicStars.h` (spawn scheduler, `spawn_shape`,
  `start_build_leg`, `schedule_dual_bridge` / `schedule_dt_macro` /
  `schedule_dtd_macro` / `schedule_macro_truncate`, `draw_shape`, `draw_sprite`).
- **Morph machinery:** `core/animation/opleg.h` (`OpLeg`: leg constructors,
  `run_op`, `arrival_mesh`, `hankin_at`/`relax_at`/`medial_at`, `slerp_vertices`,
  `finish_frame`, cohort/coarsening).
- **Mesh operators:** `core/mesh/` (`MeshOps::dual/kis/ambo/truncate/…/relax`,
  `solids.h`, `recipe.h`, `relax_bakes_generated.h`).
- **Profiling:** the `teensy-profile` skill. Both boards (COM3/COM4) are on the
  bench.

## My insights on the problem (root-caused this session)

The profile (COM3, 210 s, w16, `TRANS_SPEED=4`, `EPOCH_REVS=1920`, full cycle)
of the landed feature showed **peak 150.08 ms, 19/3296 frames spilled**, 6 red
presets. Diagnosis:

**There are two distinct cost regimes, and they compound on a few frames.**

1. **Heavy intermediates (steady build-leg frames).** The smooth morph grows the
   mesh far past its endpoints: the needle balloons from **F=32 to F=720 faces**
   (`Built Shape: …needle (V=362, E=1080, F=720, I=2160)`). Drawing that
   intermediate every frame is ~**47 ms** (`scan_mesh_raster` ≈ 720 calls/frame,
   `is_mesh_scan` ≈ 99 % of the timeline). IslamicStars is **rasterizer-bound**,
   so this is per-face scan cost, not transform. 47 ms is green but only ~15 ms
   under budget — thin.

2. **In-frame Conway rebuilds (the spikes).** The morph advances by applying
   Conway operators (`hk_conway_sweep`) over the build window. On the ~one-in-six
   frame where a new operator (or a **live relax** — see below) is applied, the
   topology rebuild runs **synchronously inside the render frame**:
   `render(47) + rebuild(~100) = ~150 ms`. These spikes come in **pairs ~6
   frames apart** (per stage transition) and are only **~19 frames across the
   whole 3296-frame cycle (0.58 %)** — but each is a visible ~7 fps hitch, and
   each such preset is red.

**Key mechanical facts that constrain the fix:**

- **You cannot fix a spike by cheapening the render alone** — the render (47 ms)
  and the rebuild (~100 ms) are *both* on that frame; killing the spike means
  addressing the **rebuild's placement** (amortize / pre-bake / move off-frame)
  or eliminating in-frame rebuilds entirely.
- **The idle slack is on the wrong frames.** `is_buffer_wait` (display-sync
  idle) is 27–48 % on the *cheap* frames and ~0 on the spike frames. There is
  ~15 ms of slack on each of the ~5 cheap frames between spikes (~75 ms total) —
  *nearly* enough to hide one ~100 ms rebuild if you could spread it, but the
  47 ms steady draw means each cheap frame only has ~15 ms of headroom, so
  amortization alone likely needs a **companion steady-draw trim** to fully
  clear 60 ms.
- **Suspect live relax in the morph.** The heaviest reds are relax-bearing
  (`…relax100…needle`, `…bevel5_relax_hk77`, `…hk35…relax_hk42`,
  `bevel2_relax_gyro`). Master bakes converged relax (`fac4b98a`) and Blocker 1
  made the branch replay that bake **for recipe steps that mirror a baked
  generator** — but the morph's **intermediate** relax steps inside OpLeg legs
  may still relax **live** (many iterations) on-frame. Confirm whether
  `hk_conway_sweep` on the spike frames includes live `MeshOps::relax`; if so,
  baking / capping those intermediate relaxes could remove most of the spike
  cheaply. **Check this first — it may be the single biggest lever.**

**The 6 red presets (all same mechanism), worst first:**
`…ambo_relax100_hk54_needle` 150.1 (8/619) · `…truncate50d_ambo_dual` 143.9 ·
`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` 130.4 · `…bevel5_relax_hk77`
111.1 · `snubDodecahedron_truncate5d_ambo_dual` 109.1 ·
`dodecahedron_bevel2_relax_gyro` 108.1. The O3 twin has the same split (peak
138.3). Every one of the 24 **clean-holds** already sustains 16 fps — the reds
are purely transition/rebuild frames.

## Decompose the work

- **Problem A — kill the in-frame rebuild spikes** (turns 6 reds green). Ranked
  candidate approaches:
  1. **Eliminate live intermediate relax** — bake or cap the relax iterations
     the morph legs run on-frame (verify with insight above). Cheapest if it
     applies.
  2. **Amortize the Conway rebuild** across the cheap slerp frames using their
     buffer-wait slack, so no frame does render + full rebuild. No memory cost,
     preserves the visual exactly. Conway ops aren't naturally incremental, so
     this may mean building the *next* stage's mesh a slice at a time on the
     preceding frames (double-buffer the stage mesh), then swapping.
  3. **Pre-bake the stage meshes at spawn** into a ring; the build window becomes
     pure slerp+draw. Deterministic kill, but costs spawn-frame time (watch the
     **epoch-commit watchdog** — an init that overruns the K-rev window traps the
     board) and persistent arena (budget is tight, ~90/106 KB — see
     `project_dtcm_bank_ledger` / IslamicStars arena survey tests).
  4. **Cap the heaviest intermediates' face count** (needle/dual-bridge 720 →
     ~360). Simplest, but changes the morph's look on those shapes — get owner
     sign-off before shipping a visual change.

- **Problem B — steady-state headroom** (buys margin so A's amortization can
  land under 60 ms, and protects the ~47–62 ms holds). IslamicStars is
  rasterizer-bound, so this is per-face scan cost:
  - **Port any raster optimizations the opchain branch is missing.** The
    pre-opchain master IslamicStars carried landed raster wins the opchain fork
    may predate: **sector-walk locate (K1+K2)**, **dead face-edge-normal
    removal**, **cached face probe flags**, **destination-alias plot**. Diff the
    opchain `draw_shape`/`Scan::Mesh::draw`/`SDF::Face` path against master and
    graft what's absent — likely the cheapest steady-state wins. (Memories:
    `project_islamicstars_sector_walk`, `project_islamicstars_o3_gap_located`,
    `project_islamicstars_rasterizer_bound`, `project_islamicstars_probe_bound`.)
  - Revisit the ripple/amplitude and column-cull levers
    (`project_islamicstars_perf_ripple_fold`, `project_plot_column_cull_ledger`).

## Method / rig (this is a measure-driven perf task)

- **Measure on device, every lever.** Host/WASM timing lies for IslamicStars
  (rasterizer-bound; GCC CSEs differently) — trust device pixel counts and
  per-frame `r=` values, not wall clock or host runs. Use the `teensy-profile`
  skill; success metric = `parse_profile.py <log> buckets` showing **every preset
  peak < 60 ms and 0 spills**, both configs, full cycle wrapped.
- The profiler already pinpoints the spikes: the `r=` per-frame lines > 62500,
  and the windows carrying `hk_conway_sweep`. Instrument the morph legs with
  `HS_PROFILE` sub-scopes (e.g. separate `is_conway_rebuild` from
  `is_mesh_scan`) if you need to attribute the ~100 ms precisely.
- **Bisect before optimizing:** confirm whether each red's spike is live-relax,
  operator rebuild, or two-mesh overlap, per preset — they may not all share one
  fix.
- **Guard determinism:** sim and device must stay bit-identical (wrap guards,
  `hs::shuffle`, no `<random>`, no `Date.now`/`Math.random` equivalents). Any
  bake must verify host==device (see the `HS_RELAX_BAKE_VERIFY` ctest pattern).
- **Don't regress correctness:** the per-leg palette crossfade and the smooth
  morph visual are owner-validated — keep them. Any face-count cap or bake is a
  visual change needing sign-off.

## Guardrails / constraints

- **No heap — arena only** on the render path. Fail-fast asserts (always-on
  trap), no OS patching, fix at the app layer.
- **ITCM budget:** the branch has ~1 KB phantasm margin (with Dissolve). New hot
  code must fit; route spawn-time code to flash via `phantasm.ld` (the OpLeg /
  half-edge-sort routing is the established pattern) rather than growing ITCM.
  The pre-commit hook builds all 3 pio images — the granule cliff surfaces at
  commit time.
- **Persistent arena** budget is tight (~90/106 KB); pre-baking stages must fit
  and is gated by the IslamicStars arena/roster tests.
- **Landing:** isolated worktree + FF-only (the `implement` skill). **Do not
  land to master until the full 24-shape cycle profiles green (< 60 ms, 0
  spills) in both configs** — owner rule. The reference-transaction hook blocks
  non-FF master moves.

## Definition of done

1. Full 24-shape cycle, **both** configs (shipping + O3): **peak render < 60 ms
   on every preset, 0 spilled frames** (`parse_profile.py buckets`, wrapped
   cycle validated).
2. Native suite green; all 3 pio images green (ITCM fits with margin).
3. Morph visually smooth and colours unchanged (or visual changes owner-approved).
4. Then re-land `opchain/dissolve` (+ Dissolve) onto master FF-only, reinstall
   the wasm simulator, and write both profile reports + update the 3 profile
   READMEs under the binary colour convention.

## Fast start

```sh
# reproduce the baseline profile on the branch (pin the tree!)
HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=<opchain/dissolve worktree> \
  bash tools/profile_one.sh IslamicStars profile    210 16 \
  "-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4"
python tools/parse_profile.py build/prof/islamicstars_ship.log buckets   # reds first
# then inspect the spike frames:
grep -E '^f [0-9]+ .* r=[0-9]+' build/prof/islamicstars_ship.log \
  | awk '{split($4,a,"=");if(a[2]+0>62500)print}'
```
