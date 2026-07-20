# Brief: generalized hot-path LUT compression

You are continuing a line of work that just landed on `c:\work\Holosphere`
(Holosphere / Phantasm POV renderer, C++20, Teensy 4.0 @ 600 MHz Cortex-M7).
One large lookup table has been compressed and moved out of the L1 cache; your
job is to find the *next* table worth the same treatment — **if any**. Read this
whole brief before touching anything. The most important section is §2: most
LUTs are NOT candidates, and the failure mode is compressing one that didn't need
it.

---

## 0. Ground rules (non-negotiable)

**Never write to `C:\work\Holosphere` directly.** Shared main tree, concurrent
peer sessions committing under the *same* git identity. Work in your own
worktree, land by rebasing onto live `refs/heads/master` and an **attached
`git merge --ff-only`** in the main tree, asserting
`git merge-base --is-ancestor` first. `master` is append-only. If a merge refuses
over main-tree WIP, **stop and surface it** — never stash or discard a peer's
work. A peer may land a build-breaking WIP on master (it happened this session: a
Teensyduino/gcc-upgrade WIP that overflowed phantasm) — if the tree won't build
phantasm, stop and surface, don't build on top of it. New files need
`git add` before `git commit <path>` (a bare `git commit <path>` silently
matches nothing for an untracked file).

**The device is a single shared Teensy.** Never run `pio run -t upload` or
`profile_capture.py` directly. The only supported path is `tools/profile_one.sh`,
which takes the host-global lock. `HS_DEVICE_WAIT=<sec>` to queue.
`profile_one.sh` builds **`/c/work/Holosphere`** (hardcoded), so **you cannot
profile a worktree — land first** (for a bit-exact change this is safe). Toggle a
diagnostic via `-D` flags forwarded to the build.

**Style:** invoke the `code-style` skill and obey it. `core/render/sdf.h`,
`core/color/color.h`, and generated tables carry whole-file clang-format drift
against local clang-format v22 — hand-format your own lines, wrap generated
tables in `// clang-format off/on`, and commit with `HS_SKIP_FORMAT=1`; **never**
run `clang-format -i` on an existing file. No `Co-Authored-By` line.

**Gates after every commit:**
- `export EMSDK=C:/work/emsdk; cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF;
  cmake --build --preset tests -j 8; ctest --preset tests` → **51/51**
- `pio run -e phantasm` → `[teensy-gate] phantasm: PASS`, and report RAM1 `code`,
  RAM1 `variables` (DTCM), FLASH `data`, **and the per-commit delta of each**.

---

## 1. What was proven — the template (commit range `f9a5b70e`..`59fed997`)

`linear_to_srgb_lut` (the display-path sRGB encode) was a **65,536-entry / 64 KB**
byte table in flash, hit 3×/pixel by the POV pack ISR. Measured on device: it
**thrashed the 32 KB L1 D-cache** against the framebuffer read — costing ~13 k
cyc/col on `isr_pack` (of ~17 k) **and** bleeding ~7 ms into every color-rich
effect's render via cache contention. Split of the pack cost:

| pack config (worst shape, `isr_pack` cyc/col) | value |
|---|--:|
| full (64 KB flash LUT) | 17,102 |
| 8 KB split-decode in flash | 8,491 (61% of thrash gone — still evicted) |
| 8 KB split-decode in DTCM | 4,842 (87%) |
| **2-region 1.5 KB split-decode in DTCM (shipped)** | **5,830 (80%)** |
| no-LUT floor | 3,047 |

The fix that shipped: a **bit-exact split-decode** that shrinks the table 64 KB →
1.5 KB and lives in **DTCM** (zero-wait, bypasses L1, so nothing evicts it). It
took the worst IslamicStars shape's *steady* render red → green, freed 64 KB
flash, cost +1.5 KB DTCM / +240 B ITCM. See `core/color/srgb_decode.h`,
`core/color/srgb_decode_lut.h`, `scripts/generate_srgb_decode.cpp`, and the
`unit_srgb_decode` test in `tests/test_color.h`.

**The technique (reuse it):**
1. **Split-decode.** A monotone `input → output` map compresses to a table of
   *base + step-threshold* per bucket, decoded as `base + (frac >= step)` — a
   single branchless compare — **iff every bucket holds ≤1 output step**. Pick the
   bucket width so that holds (the sRGB curve needed 16-wide). Verify the ≤1-step
   property in the generator; it is what makes the decode branchless.
2. **Two-region** when one uniform bucket width is too big: a fine region where
   the curve is steep + a coarse region where it's flat, selected by one branch.
   Got sRGB to 1.5 KB (from 8 KB uniform). Costs ~+1–2 cyc/lookup for the region
   select.
3. **Placement.** If it thrashes L1, it must go in **DTCM** (a non-`const` static,
   filled once at static-init from a `const` flash source; a `const` table routes
   to flash via `phantasm.ld` and stays cacheable). `always_inline` the decode
   (the `-Os` build otherwise emits a *call* per lookup — that alone was ~2 k
   cyc/col of overhead).

---

## 2. The candidate criterion — READ THIS BEFORE COMPRESSING ANYTHING

**Compression helps only a table that is BOTH large enough to thrash L1 AND hot.**
The sRGB LUT was an extreme case (64 KB ≫ 32 KB L1). Compressing a table that
already fits L1 **adds decode compute for zero cache benefit** — a net loss. Do
not compress on principle.

A candidate must clear **both** bars:
- **Large:** the table's *working set* competes with the 32 KB L1 (i.e. hundreds
  of KB, or tens of KB accessed scattered so it can't stay resident). A 256-entry
  / 512 B table is resident forever — skip it.
- **Hot:** read per-pixel / per-probe / per-column, so its cost multiplies by a
  large trip count. A table touched once per frame is irrelevant.

**Measure the thrash before compressing.** The `9bd8713b` diagnostic pattern:
flag-gate a variant that nulls out the table (replace its lookup with a cheap
stand-in) and compare `isr_pack` / the relevant scope full-vs-nulled on device.
If nulling the table barely moves the scope, it wasn't thrashing — stop.

---

## 3. The actual candidates (surveyed; sizes are real, verdicts are hypotheses)

| table | where | size | hot? | first read |
|---|---|--:|---|---|
| `linear_to_srgb_lut` | color_luts.h | 64 KB | per-pixel (pack) | **DONE** — the template |
| `GAMUT_LUT` | gamut_lut.h | **512 KB flash** (512×256×2, uint16) | gamut-clip effects | **Most promising, but already downsampled at runtime** — `init_gamut_lut` downsamples the flash master into an arena LUT (`g_gamut_lut`), so the *runtime* footprint is far smaller and depends on the downsample factor. Investigate the runtime working set and whether the gamut-clip path (`color.h:840`) thrashes; the flash master is cold (init only). |
| trig / "pixel→vector" | `TrigLUT` (geometry.h, plot.h) | `sin_theta`/`cos_theta[W=288]`, `sin_phi[H_VIRT]` ≈ 1–2 KB each | per-pixel in the scan (`Vector p(sp*cos_theta[x], cp, sp*sin_theta[x])`) | Hot but **small** — likely resident, likely NOT a thrash candidate. This is the "pixel→vector" path. **Measure before assuming** — but expect it fits L1 and compressing it only adds compute. |
| `srgb_to_linear_lut` | color_luts.h | 512 B (256×uint16) | per-pixel where used | **Skip** — fits L1 trivially. |
| SDF congruence-class LUT | `build_canonical_distance_lut`, sdf.h/scan.h | per-shape, staged | per-probe | Deforming-effect-hostile (see `project_congruence_class_lut_landed` memory) — facility only; do not re-wire. Likely not a compression target. |
| FastNoiseLite gradient/randvec | vendor/FastNoiseLite.h | 256–512 entries | per-noise-sample | Small; skip unless a noise effect profiles LUT-bound. |

**Honest read:** after sRGB, the field thins out fast. `GAMUT_LUT` is the only
other genuinely large table, and it is *already* mitigated by downsampling — the
open question is whether its runtime form still thrashes on the effects that clip
to gamut (profile one of those, not IslamicStars). Everything else is small
enough that the criterion in §2 probably rules it out. **A well-measured "nothing
else is worth it" is a valid and likely outcome — report it rather than forcing a
compression.**

---

## 4. The hard constraint: DTCM is nearly full

After the sRGB decode, **DTCM has only ~1.2 KB free** (10 banks = 327,680 B;
314,176 B variables + 12,288 B stack). The arena owns the rest and is sized for
GSReactionDiffusion (~298 KB, ~7 KB global margin) — so you **cannot** shave the
arena without GS optimization (its `static_assert` refuses it). So:

- A new DTCM-resident table must fit ~1.2 KB, **or** you must first free DTCM
  (shrink another table, or GS work — out of scope unless the win justifies it).
- **Before spending scarce DTCM, test whether a small-enough compressed table
  stays resident in *cacheable flash*.** The sRGB 8 KB flash table still thrashed
  (61%), but a ~1.5 KB table might survive L1 under render pressure — if it does,
  it needs **no DTCM at all**. Measure flash-resident vs DTCM (the `48cfcb0c`
  three-way did exactly this) before deciding.

ITCM headroom: ~5.8 KB (`code` 190,760 / 196,608). The stack floor (12,288 B) is a
safety reserve — do not shrink it.

---

## 5. Methodology — the confounds that burned this session

Every "quick diagnostic" this session introduced a confound. Do not repeat them:
- **`const`/single-value read confound.** A diagnostic that reads the *same*
  input every iteration trivializes the table's working set — it measures
  "table removed", not "table resident". Use *real, varied* inputs.
- **Access-pattern confound.** Changing *which* entries are read changes the
  table's working set, not just its footprint. You cannot cleanly isolate
  "footprint" from "workload" with a lookup-remap; only a real layout change does.
- **Inlining confound.** Under `-Os` the decode is *not* inlined unless
  `always_inline`d — a lazy-copy guard or an out-of-line function adds a call per
  lookup (~2 k cyc/col) and poisons the measurement. **Disassemble and confirm 0
  calls** (`objdump -d ... | grep 'bl.*<your_decode>'`).
- **Stale flash.** `profile_one.sh`'s first upload after a peer session can report
  SUCCESS *without flashing*, handing you the previous firmware. **Gate every
  run:** a FULL baseline must reproduce its known value (sRGB pack ~17 k cyc/col)
  before you trust the experimental capture. Build a `pk()` check into the script.
- **Epoch wrap.** A 210 s capture can wrap the profiler epoch (frame N→1),
  tripping `parse_profile.py validate` (`INVALID`) even though the pre-wrap data
  is fine — but the wrap also corrupts `buckets` tier counts. Read per-shape
  numbers from **clean-hold windows** (`windows --scope`), not the wrapped
  `buckets`, and prefer a shorter capture or trust a VALID run.
- **Trust pixel counts / scope cycles, not cross-session wall-clock.** Run-to-run
  `scan_mesh_raster` noise is ~0.3 ms. Attribute by *scope* in one capture.

Standard capture: `bash tools/profile_one.sh <Effect> profile 210 16 -D
HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4`, then `python
tools/parse_profile.py <log> validate` and `... presets --gate scan_mesh_raster`.
`isr_pack` is an `HS_ISR_PROFILE` scope (atomic, ISR-safe) — read its
`min/avg/max` cyc/col from the worst-shape window. Pick the *effect* whose hot
path uses the table under test (gamut → a gamut-clipping effect, not IslamicStars).

---

## 6. Proof bar — bit-exact, always

Output must be **bit-exact** vs the original table. The generator (a) asserts the
≤1-step property, (b) verifies `decode(v) == table[v]` for **every** input, then
(c) emits the compressed table. Add a CI unit test that re-checks the exhaustive
equivalence (mirror `test_linear_to_srgb8_decode_matches_lut` in
`tests/test_color.h`) and wire it into the runner (memory:
`feedback_new_tests_into_ci`). **Include a negative control** — perturb the decode
and confirm the test fails (a `>=`→`>` flip moved 236 entries; a matching hash on
an unexercised path is worthless). The framebuffer oracle is the wrong tool for a
display-transform LUT (the transform is downstream of the framebuffer the oracle
hashes) — the exhaustive `decode==table` test *is* the complete proof for a pure
lookup substitution.

---

## 7. Where to start

1. **Profile a gamut-clipping effect** with a diagnostic that nulls the runtime
   gamut LUT, and see whether `g_gamut_lut` thrashes at all. If it doesn't, the
   generalized-compression thesis is mostly spent after sRGB — say so.
2. If it does, characterize the runtime working set (downsample factor) and decide
   split-decode vs. a coarser downsample, weighed against the ~1.2 KB DTCM budget
   and the flash-resident test in §4.
3. Only then touch trig/`pixel→vector` — and expect §2 to rule it out.

Report each candidate as: measured thrash (or absence), the compression + placement
if warranted, the bit-exactness proof with sample counts, `ctest`, and phantasm
RAM1/DTCM/FLASH deltas. **A rigorously-measured negative is the expected and
valuable outcome for most of these.**
