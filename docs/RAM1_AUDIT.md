# RAM1 Byte-Level Audit — Teensy 4.0 FlexRAM

**Build analyzed:** `.pio/build/phantasm/firmware.elf` — the **real** phantasm config: `-Os` + newlib-nano (`libc_nano`/`libstdc++_nano`, via `tools/teensy_nano.py`) + the custom flash-routing linker script `tools/phantasm.ld`. teensy@5.0.0 / Teensyduino 1.59.0 / arm-gcc 11.3.1, all effects via `HS_EFFECT_LIST`.
**Toolchain:** `arm-none-eabi-{size,nm,readelf}` + `teensy_size` (authoritative region split).

> ⚠️ **This build overflows RAM1 by 27,680 bytes — a regression of exactly one 32 KiB FlexRAM bank from the prior +5,088 B baseline.** The committed `holosphere` env (RingSpin only) *fits* with ~31 KB headroom; the `phantasm` env (every effect compiled in) is the stressing case and is what this audit catalogs.
>
> **What changed (byte-exact):** at the +5k baseline (budget note: RAM1 used 519,200) ITCM code fit *exactly* in 5 × 32,768 = **163,840 B (5 banks)**, DTCM vars 355,360, stack +5,088. ITCM code has since grown **10,088 B to 173,928**, spilling into a **6th bank**. Banks are all-or-nothing, so that crossing transferred a whole 32,768 B from DTCM to ITCM: `5,088 − 32,768 = −27,680`. DTCM data is unchanged (355,360 then and now) — **this is a pure code-size regression, not a data regression.** The fix is therefore entirely on the ITCM side: claw code back under 163,840 (see §5 R1) and the 6th bank — and the overflow — disappears.

---

## 1. What RAM1 is

On the IMXRT1062 (Teensy 4.0), **RAM1 = the 512 KiB FlexRAM** (`524,288 B`), partitioned at boot into 32 KiB banks assigned to either **ITCM** (instruction TCM — fast code) or **DTCM** (data TCM — `.data`/`.bss`/stack). RAM2 (OCRAM, `0x2020_0000`, holding `.bss.dma` = the two framebuffers) is a *separate* 512 KiB region and is **out of scope** here.

| RAM1 sub-region | Holds (sections) | Address |
|---|---|---|
| **ITCM** | `.text.itcm`, `.fini` | `0x0000_0000` |
| **DTCM** | `.data`, `.bss`, stack (grows down from top) | `0x2000_0000` |

`teensy_size` ground truth:

```
RAM1: variables:355360  code:173928  padding:22680   free for local variables:-27680
```

| Component | Bytes | % of 524,288 | Notes |
|---|---:|---:|---|
| DTCM `variables` (`.data`+`.bss`) | 355,360 | 67.8 % | dominated by one symbol — see §3 |
| ITCM `code` (`.text.itcm`+`.fini`) | 173,928 | 33.2 % | rounds up to **6 × 32 KiB banks** |
| ITCM bank **padding** | 22,680 | 4.3 % | dead space inside the 6th bank |
| Stack (`free for local variables`) | **−27,680** | — | **overflow** |
| **Total demanded** | **551,968** | 105.3 % | over the 524,288 wall by 27,680 |

### The bank arithmetic (the central lever)

ITCM is allocated in whole 32,768-byte banks; DTCM gets whatever banks are left.

```
ITCM code 173,928  → ceil(173928 / 32768) = 6 banks = 196,608 allocated
                     padding = 196,608 − 173,928 = 22,680  (wasted RAM1)
DTCM available = 524,288 − 196,608 = 327,680  (10 banks)
DTCM used      = 355,360
stack free     = 327,680 − 355,360 = −27,680   ← overflow
```

**Crossing one bank boundary is worth 32,768 bytes of DTCM.** If ITCM code drops to **≤ 163,840** (5 banks) — a cut of just **10,088 B** — DTCM available jumps to 360,448 and stack free becomes **+5,088 B**. That single structural move converts the overflow into a (tight) fit. Pushing the stack to the budgeted 16 KiB floor needs that bank drop **plus** ~11 KB of DTCM trims (§3 / §5).

---

## 2. ITCM catalog — 173,928 B of code in fast RAM

Teensy places **all** code in ITCM by default; only functions explicitly tagged `FLASHMEM`/`PROGMEM` go to flash (which has **1.63 MB free**). So ITCM is "every function not opted out." The padding (22,680 B) is unreclaimable *except* by crossing the bank boundary above.

### 2a. Top ITCM consumers (per-frame HOT — keep in ITCM)

These run per-pixel / per-frame; moving them to flash **would** cost performance — leave them.

| Bytes | Symbol (demangled) |
|---:|---|
| 3,312 / 3,224 / 3,052 / 2,824 / 2,240 | `Plot::rasterize<288,144,…>` (5 pipeline specializations) |
| 1,816 | `Filter::Pixel::Feedback<288,144>::flush` |
| 1,200 ×8 | `Scan::scan_region<288,144,…>` (SDF Flower/Star/Ring/SphericalPolygon/PlanarPolygon/DistortedRing/Face) |
| 1,172 | `Voronoi<288,144>::draw_frame` |
| 1,152 | `usb_isr` (ISR — must stay) |
| 1,082 | `ShapeShifter<288,144>::dispatchScan` |
| 872 / 796 | `Scan::scan_region` / `Scan::Volume::draw` (torus/twist) |
| 764 | `FastNoiseLite::SingleOpenSimplex2<float>` (per-pixel noise) |
| 708 ×2 | `SDF::Face::distance` / `SDF::SphericalPolygon` ctor |
| 660 | `pov::sync::SyncBoard::tick` |
| 652 | `DistortedRing::drawFn`, `Liquid2D::draw_frame` shader |
| 644 | `Animation::ParticleSystem::step_particle` |

### 2b. Setup-only (COLD) code currently in ITCM — relocate to flash, zero perf cost

These run **once at effect construction / mesh build**, never per-frame, yet occupy fast RAM. **Several are already tagged `FLASHMEM` in source but landed in ITCM anyway** (§5, finding R1).

| Bytes | Symbol | Source tag | In flash? |
|---:|---|---|---|
| 1,528 | `MeshOps::classify_faces_by_topology` `.constprop.0` | — | **no (ITCM)** |
| 1,360 | `MeshOps::snub<PolyMesh>` | `FLASHMEM` (conway.h:1087) | **no (ITCM)** ⚠ |
| 948 | `MeshOps::truncate<PolyMesh>` | `FLASHMEM` (conway.h:659) | **no (ITCM)** ⚠ |
| 930 | `MeshOps::compile_hankin<PolyMesh>` | — | **no (ITCM)** |
| 868 | `Solids::SolidBuilder::expand` | `FLASHMEM` (solids.h) | **no (ITCM)** ⚠ |
| 864 | `Solids::SolidBuilder::relax` | `FLASHMEM` (solids.h) | **no (ITCM)** ⚠ |
| 696 | `MeshOps::ambo<PolyMesh>` | `FLASHMEM` (conway.h) | **no (ITCM)** ⚠ |
| 652 ×2 | `HalfEdgeMesh` ctor (C1+C2) | — | no (ITCM) |
| 446 | `std::__adjust_heap` helpers (sort inside classify_faces) | — | no (ITCM) |
| 192 | `compile_hankin` lambda | — | no (ITCM) |
| 868 (`MeshOps::truncate` already counted)| … | | |
| ~430 | `PolyMesh` ctor, `MeshOps::normalize/narrow_*/hash_combine/fmix32`, `ArenaVector<…>::bind/emplace/append` `.part.0` clones | — | no (ITCM) |

**Setup-class subtotal in ITCM ≈ 27.5 KB** (grep-bucketed; the construction-only core is ~9–10 KB even excluding anything ambiguous). Relocating ~10 KB of it clears the 6th bank.

> The `Solids::Archimedean::snubCube`/`snubDodecahedron` factories *did* honor `FLASHMEM` (they sit in flash section `[2]`), and the linker emitted `…_veneer` trampolines for `SolidBuilder::{expand,relax,snub}` — proof the call crosses the flash↔ITCM boundary, i.e. the builders are *meant* to be cold/flash-resident.

---

## 3. DTCM catalog — 355,360 B (every byte)

DTCM is **96.6 % one object.** `.bss` = 352,608, `.data` = 2,752.

| Bytes | Symbol | Kind | Reclaimable? |
|---:|---|---|---|
| **343,040** | `global_arena_block` (`_ZL18…`) | effect scratch arena | **§5 R4** — right-size (the only large lever) |
| 2,500 | `hs::random()::gen` = `std::mt19937` | RNG state | **§5 R2** — swap to small PRNG (−2,470) |
| 1,280 | `POVSegmented<288,4,480>::ledController_` | LED driver state | required |
| 1,152 | `TrigLUT<288,144>::sin_theta` | runtime LUT | **§5 R3** → flash |
| 1,152 | `TrigLUT<288,144>::cos_theta` | runtime LUT | **§5 R3** → flash |
| 704 | `_VectorsRam` (`.data`) | relocated IRQ vectors | core — keep |
| 640 | `endpoint_queue_head` (`.data`) | USB EP queue | core — keep |
| 588 | `TrigLUT::sin_phi` | runtime LUT | **§5 R3** → flash |
| 588 | `TrigLUT::cos_phi` | runtime LUT | **§5 R3** → flash |
| 588 | `PhiLUT<144>::data` | runtime LUT | **§5 R3** → flash |
| 316 | `POVSegmented::sync_` | POV sync state | required |
| 312 | `__sf` (`.bss`) | newlib stdio `FILE` ×3 | **§5 R6** if no stdio |
| 256 | `rx_transfer` | USB DMA desc | core — keep |
| 128 ×4 | `isr_table_gpio1..4` (`.data`) | GPIO ISR dispatch | core — keep |
| 128 | `tx_transfer` | USB DMA desc | core — keep |
| 96 | `funct_table` | USB | core — keep |
| 80 | `SPI` (`.data`) | SPI object | **§5 R7** if SPI unused |
| 76 | `_impure_data` (`.data`) | newlib reent | **§5 R6** with `__sf` |
| 40 | `microsoft_os_compatible_id_desc` (`.data`) | WebUSB/MS-OS desc | **§5 R7** if unused |
| 32 ×3 | `hs_ctr_*` profiling counters (`scan_face_setup`, `filter_blend`, `mesh_raster`) | `HS_PROFILE` instrumentation | **§5 R8** — drop in release |
| 32 ×2 | `endpoint0_transfer_{data,ack}` | USB | core — keep |
| 28 | `HardwareSerialIMXRT::s_serials_with_serial_events` | core | keep |
| 22 ×2 | `usb_string_serial_number(_default)` | USB | keep |
| 18 ×2 | `microsoft_os_string_desc`, `device_descriptor` | USB | partly trimmable (R7) |
| 16 ×4 | `Serial`, `FastLED`, `scratch_arena_{a,b}`, `persistent_arena` | objects/handles | keep |
| 16 ×2 | `rx_index`, `rx_count` | USB | keep |
| <16 each | `_impure_ptr`, `__brkval`, `F_CPU_ACTUAL`, `F_BUS_ACTUAL`, `IntervalTimer::nvic_priorites`, `rx_list`, … | core scalars | keep |

`.data` (initialized — costs RAM1 **and** a flash copy of the init image): **2,146 B**, all of it Teensy core/USB infrastructure (`_VectorsRam`, `endpoint_queue_head`, `isr_table_gpio*`, `SPI`, descriptors). Nothing application-owned — no action.

**Everything in DTCM that is *not* the arena totals only 12,320 B.** Of that, ~7.4 KB is genuinely reclaimable (RNG 2,470 + LUTs 4,068 + stdio ~390 + profiling ~96 + MS-OS descriptors ~58). The arena is the structural consumer.

---

## 4. Region boundary sanity (out of scope, for reference)

| Section | Region | Bytes |
|---|---|---:|
| `.text.headers/.code/.progmem/.csf` + `.ARM.exidx` | FLASH (`0x6000_0000`) | 211,260 code / 182,700 data — 1.63 MB free |
| `.bss.dma` (`buffer_a`+`buffer_b` framebuffers) | **RAM2/OCRAM** | 518,272 (≈6 KB free) |

RAM2 is full but separate; do **not** move RAM1 data there — there is no room.

---

## 5. Reclamation plan — prioritized, performance-preserving

Ordered by impact. Every item is **perf-neutral or perf-positive** unless flagged.

### P0 — Structural (fixes the overflow), **zero runtime cost**

**R1 — ✅ Repatriate `FLASHMEM` mesh/solid construction code that leaked into ITCM.** *(DONE: ITCM code 173,928 → 159,976 B — the 6th 32 KiB bank is released, free-for-locals −27,680 → +5,088 B, overflow resolved.)*
The Conway operators and `SolidBuilder` are tagged `FLASHMEM` in source but `readelf` shows them in section `[4] .text.itcm`. Root cause: under `-O3`, GCC's IPA passes emit `.constprop`/`.isra`/`.part` **clones** and template/`WEAK` COMDAT instantiations that **drop the `section(".flashmem")` attribute**. Confirmed leaked symbols (all setup-only): `MeshOps::{snub 1360, truncate 948, ambo 696, compile_hankin 930+192, classify_faces 1528+446}`, `SolidBuilder::{expand 868, relax 864}`, `HalfEdgeMesh` ctor 652×2.
**Fix landed:** an `HS_COLD` macro (`FLASHMEM __attribute__((noinline, noclone))` on GCC; no-op off-device) applied to the `static` free-function Conway/Hankin operators (conway.h, hankin.h) and `classify_faces_by_topology` / `compile` (mesh.h). The measured win comes from `noinline`/`noclone` collapsing the duplicate inlined + `.constprop` copies these operators stamped across the ~40 solid factories — the `static` template instantiations themselves stay in ITCM (the `.flashmem` attribute is not honored on them) but now exist as single copies. COMDAT members (`SolidBuilder` methods, `HalfEdgeMesh` ctor, `build_from_flat`) were left unmarked: a section attribute on a COMDAT function is a section-type conflict.
Remedies (verify each with `readelf -sW | grep … → Ndx 3`):
- Add `__attribute__((noinline, noclone))` alongside `FLASHMEM` on these definitions so IPA can't spawn attribute-less clones.
- Or move the heavy template bodies to an explicitly-instantiated `.cpp` compiled with the section attribute.
- Or wrap each in a non-template, out-of-line `FLASHMEM` shim that the hot path calls once.
Target: ITCM ≤ 163,840 (cut ≥ 10,088). Moving the list above (~9–10 KB) plus `PolyMesh` ctor / `MeshOps::normalize`/`narrow_*` / `ArenaVector::*.part.0` (~0.4 KB) clears it.

**Net after R1:** stack free −27,680 → **+5,088**. Combine with one DTCM item below (e.g. R2+R3 = ~6.5 KB) to reach a comfortable stack margin; with all of R2+R3+R6 you reach ~+12 KB, near the 16 KiB budgeted floor.

### P1 — DTCM data, perf-neutral or **faster**

**R2 — ✅ Replace `std::mt19937` with a small PRNG.** *(DONE: DTCM variables 355,360 → 352,864 B, −2,496 B; stack +5,088 → +7,584 B. Added `hs::Pcg32` (PCG XSH-RR 64/32, 16 B state) behind `hs::random()`; device + host use the identical seeded type so parity holds; all 35 native tests pass.)* *(−2,470 B DTCM, and faster.)*
`hs::random()` is a process-wide `mt19937(1337)` whose 624-word state is the 2,500-byte DTCM symbol. A `pcg32` or `xoshiro128**` has 16–32 B of state and a shorter critical path. ⚠ **Parity caveat:** `platform.h` documents this generator as the *single* determinism source shared by sim and device; swap **both** sides to the identical new PRNG in one change and re-tune any sequence-sensitive effect. Net: −~2,470 B RAM1, +throughput.

**R3 — ✅ (cos_theta-only, in-DTCM) Recover `cos_theta` from `sin_theta`.** *(DONE: DTCM 352,864 → 352,000 B, −864 B; stack +7,584 → +8,448 B. The full flash move was REJECTED — moving per-pixel LUT reads to flash risks L1-D-cache eviction by streaming framebuffer writes, a hot-path cost only measurable on hardware. Instead `cos_theta[x]` is read as `sin_theta[x + W/4]` from a sin table extended by W/4 entries: stays in zero-wait DTCM, no modulo, constant-offset load → zero perf cost. All 35 native tests pass.)* *(Original proposal: −4,068 B DTCM via constexpr flash move, marginal perf.)*
`TrigLUT<288,144>::{sin,cos}_theta` (1,152 ×2), `{sin,cos}_phi` (588 ×2) and `PhiLUT<144>::data` (588) are runtime-filled `static std::array` in `.bss`. Make them `static constexpr` (a `constexpr` sine polynomial — `std::sinf` isn't `constexpr`) so they land in `.rodata`/FLASH (1.6 MB free). ⚠ Reads are per-pixel (`filter.h:1006`), but the tables are tiny and stay in the flash cache — cost is marginal. Float values shift slightly vs runtime `sinf`; acceptable under the project's "parity = integer wrap guards, not floats" contract. *Cheaper half-measure:* derive `cos` from `sin` via a quarter-table index offset → −1,740 B with zero precision change and no flash move.

### P2 — The arena (the real consumer)

**R4 — Right-size `global_arena_block` (343,040 B = 96.6 % of DTCM).** *(every byte trimmed → a byte of stack, zero perf cost up to the true high-water.)*
This is the dominant lever and the budget's tuning target (`min 327,680 / max 360,448`). The repo already stack-paints to find the worst-effect stack peak (`tests/stack_measure.cpp`); apply the same high-water instrumentation to the arena across **all** effects, then set the arena to `peak + margin`. The fail-fast trap (project policy: arena overflow must crash, never mask) makes an aggressive cut safe — a regression traps loudly rather than corrupting. Even a 5 % trim (~17 KB) more than covers the residual overflow after R1.

### P3 — Small / marginal wins (exhaustive)

**R5 — Eliminate the dead base-object constructor.** *(up to −3,652 B ITCM + −708 B.)*
`SDF::Face` ctor is emitted twice (C1 complete + C2 base, **3,652 B each**, both in ITCM — and Face is constructed per-frame at `scan.h:802`, so this is hot fast-RAM). `SDF::SphericalPolygon` ctor likewise (708 ×2). If `Face`/`SphericalPolygon` are never used as base classes, the C2 (base) variant is dead. Verify it is unreferenced, then ensure `-ffunction-sections` + `-Wl,--gc-sections` removes it (arm-gcc ld has no `--icf`, so identical-code-folding won't merge them automatically — GC is the lever). **Do not** move the surviving ctor to flash (hot path).

**R6 — Drop newlib stdio if unused.** *(−~390 B DTCM: `__sf` 312 + `_impure_data` 76 + `_impure_ptr` 4.)*
These are the `FILE` table for `stdin/out/err`. If the device only uses `Serial.print*` (not `printf`/`fwrite`), nothing references them and they can be excluded (nano specs + avoid pulling `_impure_ptr`). Verify no `<cstdio>` stream use on the device path.

**R7 — Strip unused USB descriptors / peripherals.** *(−~140 B DTCM.)*
`microsoft_os_compatible_id_desc` (40) + `microsoft_os_string_desc` (18) are WebUSB/MS-OS-2.0 descriptors — drop if no WinUSB binding is needed. `SPI` object (80) is dead if the LED path is FlexIO/DMA-only (confirm no `SPI.begin`).

**R8 — Compile out `HS_PROFILE` counters in release.** *(−~96 B DTCM + the ITCM of the profiling IIFEs.)*
The `hs_ctr_scan_face_setup`/`filter_blend`/`mesh_raster` counters (32 B ×3) and their `HS_PROFILE` scopes are instrumentation. Gate them behind a release flag so they (and the cycle-counter code they pull into ITCM) vanish from the shipping image.

**R9 — ITCM bank padding (22,680 B) is not directly addressable.** It is intra-bank slack; it only converts to usable DTCM when code crosses the 163,840 boundary (R1). Listed for completeness — there is no "free it in place" option.

---

## 6. Bottom line

| Action | RAM1 effect | Perf |
|---|---:|---|
| R1 repatriate leaked FLASHMEM (drop 1 ITCM bank) | **+32,768 B to DTCM** (−27,680 overflow → +5,088 stack) | none |
| R2 mt19937 → pcg32 | +~2,470 B | faster |
| R3 trig LUTs → flash | +~4,068 B (or +1,740 half-measure) | marginal |
| R4 arena right-size | +N KB (largest lever) | none ≤ high-water |
| R5 dead C2 ctor GC | +up to ~4,360 B ITCM | none |
| R6 stdio drop | +~390 B | none |
| R7 USB/SPI strip | +~140 B | none |
| R8 profiling off | +~96 B + ITCM | faster |

**R1 alone resolves the overflow** by reclaiming a full FlexRAM bank at zero performance cost — it is the highest-leverage fix and should land first. R2–R8 then restore a healthy stack margin without touching the per-frame hot paths. The arena (R4) is the only place large amounts of further DTCM can come from, and it is safe to trim under the project's fail-fast arena policy.
