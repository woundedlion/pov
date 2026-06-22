# Design Spec: Teensy 4 CI Gate — Build, Warning Hygiene, Image Size & Memory Layout

**Status:** Phase-0 scaffold + toolchain spike done. Landed & self-tested (report-only):
`platformio.ini` + the `tools/teensy_pre.py` sketch-placement hook, the post-build gate
(`tools/teensy_gate.py` + `tools/teensy_gate_extra.py`), budgets (`tools/teensy_budgets.json`), the
warning ratchet (`tools/teensy_warnings.py` + baseline), the host self-tests
(`tools/teensy_gate_tests/`, green), the `teensy-gate-tests` CI job, the `just teensy-size` recipe,
and the Teensy-4.0 doc reconciliation (§14.5).

**Phase-0 spike results (2026-06) — see §16:** the toolchain question is resolved favourably
(branch 1, exact bench parity): `platform = teensy@5.0.0` natively delivers Teensyduino 1.59.0 +
arm-none-eabi-gcc 11.3.1, `platformio==6.1.19`, `FastLED@3.10.3`. Two config corrections the real
build forced: (a) the Teensy core default `-std=gnu++17` must be `build_unflags`'d so the engine's
C++20 compiles; (b) PlatformIO's sketch discovery globs `$PROJECT_SRC_DIR/*.ino` and ignores
`build_src_filter`, so the spec's filter-only plan (§5/§6) can't place the `.ino` — a pre-script
selects it instead. The build now compiles `core/*.cpp` + the converted sketch cleanly.

**The gate earned its keep, then went green (2026-06).** The first headless device build surfaced
**real device-only breakage** invisible to the WASM/native CI; all of it is now FIXED (§16):
Arduino `TWO_PI` macro collisions (renamed to `kTwoPi` in `effects/{Flyby,Liquid2D,Raymarch}.h`),
the wrong FastLED pin (3.10.3 → the bench's **3.4.0**, which moved types into `namespace fl`),
`effects_legacy.h` `Pixel16`→`CRGB` explicit-conversion bit-rot, an ambiguous `pov` vs the `namespace
pov` (sketch var → `g_pov`), and a `color.h` include collision with FastLED's own `color.h`
(`hardware/` now uses `"core/color.h"`). **Result:** **Holosphere builds GREEN and PASSES the gate
end-to-end** on the real ELF (calibrated budgets + the real `_ZL18global_arena_block` mangled arena
symbol). **Phantasm compiles + links but OVERFLOWS RAM1 by ~243 KB** (all effects → ~381 KB ITCM +
~380 KB DTCM ≫ 512 KB) — the gate correctly fails it; shrinking it (-Os / FLASHMEM / fewer effects /
smaller arena) is owner-scope (§4.1 relief valve). The warning baseline is captured (1 first-party
warning — 55 `HS_CHECK` sites dedupe to one). Remaining: enable the `teensy-size` CI job
(`vars.TEENSY_GATE_ENABLED`) once Phantasm fits.
**Scope:** Add an automated gate that compiles the Teensy 4 firmware images, then validates
build success, warning hygiene, image size, and memory-region layout against checked-in budgets.
**Out of scope:** On-device flashing/running, hardware-in-the-loop tests, the WASM/native-test
paths (already covered by `ci.yml`).

---

## 1. Problem statement

The project has three build targets sharing one engine: the **WASM** simulator, the **native
unit suite**, and the **Teensy firmware** (`targets/Holosphere/Holosphere.ino`,
`targets/Phantasm/Phantasm.ino`). The first two are covered exhaustively by
`.github/workflows/ci.yml` (sharded Clang tests, ASan/UBSan, Windows, WASM build + runtime smoke
+ install provenance). **The firmware path is built nowhere automated** — it compiles only on a
developer machine through Visual Micro (VMicro).

Consequences of the gap:

- **Device-only code is uncompiled in CI.** Everything behind `#ifdef ARDUINO` / `CORE_TEENSY`
  — the eDMA register programming in `hardware/dma_led.h`, `pov_single.h`, `pov_segmented.h`,
  `IntervalTimer` ISR wiring, `DMAMEM`/`FLASHMEM` placement — can break without any signal until
  a human opens VMicro. `docs/device_host_divergence_ledger.md` row 4 records this directly:
  *"the bare DMA/register writes are device-only by nature and not host-reachable."*
- **No size or memory-layout signal.** The two framebuffers (`MAX_W*MAX_H` = 288×144 ×
  `sizeof(Pixel16)` = 6 B = **243 KiB each**; 486 KiB together) live in `DMAMEM` (OCRAM/RAM2,
  512 KiB total), the **335 KB arena lives in DTCM/RAM1** (`global_arena_block`, no `DMAMEM` — it
  is the largest single RAM1 object, [memory.cpp:43-53](../core/memory.cpp)), and the 90 KB
  reaction-graph table lives in flash. Both buffers are fixed at `MAX_W*MAX_H` for **both**
  targets regardless of virtual resolution, so OCRAM is genuinely tight (486 KiB of 512 KiB). A
  change that pushes any region over budget — or silently relocates a large buffer from flash to
  RAM, or grows the arena and squeezes the DTCM stack — is invisible until it either fails to link
  on-device or crashes at runtime where *"there is no debugger attached and no console to read"*
  (README §2).
- **No warning hygiene baseline.** arm-none-eabi-gcc with `-Wall -Wextra` surfaces warnings the
  Clang/Emscripten builds do not. Today nothing watches them.

This spec designs a **compile + link + measure** gate (no on-device execution) that closes the
compile/link/size blind spot, run headlessly in CI and reproducible locally.

---

## 2. Goals & non-goals

### Goals
1. Build both Teensy 4 firmware images headlessly, deterministically, version-pinned.
2. Fail CI on any **compile or link error** in device-only code paths.
3. Enforce **warning hygiene** (no *new* warnings, via a committed baseline ratchet — §7.2).
4. Enforce **image size** budgets per memory region (FLASH, RAM1/ITCM+DTCM, RAM2/OCRAM).
5. Validate **memory layout** invariants beyond raw totals — e.g. the framebuffers stay in
   OCRAM, the reaction-graph table stays in flash, the DTCM stack headroom is preserved.
6. Be **reproducible locally** with one command, so a developer who uses VMicro day-to-day can
   pre-flight the gate before pushing.
7. **Coexist with VMicro** — VMicro remains the edit/flash/debug tool; this gate never owns that
   workflow and the two must not fight over the working tree.

### Non-goals
- Replacing VMicro as the local development or upload path.
- Running effects on real hardware (no runner has a Teensy; that stays manual — README §11,
  memory `project_no_device_build_in_ci`).
- Re-validating sim≠device numeric divergence — that is the host device-value tests' job
  (`h_offset_renorm_check`, `fastmath_clamp_check`, the divergence ledger). This gate is about
  *does it build for the device and fit*, not *does it compute the same pixels*.

---

## 3. Why PlatformIO

The firmware today builds via VMicro, a Visual Studio GUI front-end over the Arduino/Teensyduino
toolchain. VMicro is excellent for interactive edit/flash/debug but is a poor CI citizen: it is
Windows + Visual Studio bound, GUI-oriented, and its project state (`__vm/`, `*.vcxproj`) is
machine-local and already `.gitignore`d.

Candidate headless build drivers:

| Option | Verdict |
|---|---|
| **PlatformIO** (`platform = teensy`, `framework = arduino`) | **Chosen.** Pure CLI, cross-platform, Linux-runner friendly, declarative `platformio.ini`, exact version pinning of platform/framework/toolchain/libs, built-in per-region size report, post-build hook API (`extra_scripts`) for custom size/layout parsing, first-class GitHub Actions support, trivial local reproduction (`pio run`). |
| **arduino-cli + teensy** | Viable but heavier to pin reproducibly; Teensy support via the third-party board manager URL is less ergonomic to cache; weaker hooks for size parsing. |
| **VMicro headless / MSBuild** | Ties CI to Windows + Visual Studio + a licensed/GUI tool; fragile and slow on runners. Rejected. |
| **Bare Makefile over Teensyduino** | Maximum control, maximum maintenance; we'd reinvent dependency/version management PlatformIO already does. Rejected. |

PlatformIO uses the **same underlying toolchain** Teensyduino/VMicro use (`arm-none-eabi-gcc`,
the Teensy core, `teensy_size`), so a green PlatformIO build is strong evidence the VMicro build
is also healthy — without claiming bit-identical artifacts.

---

## 4. Coexistence with VMicro

The two systems live side by side, each owning a different job:

| Concern | VMicro (local) | PlatformIO (CI + local pre-flight) |
|---|---|---|
| Edit / IntelliSense / flash / debug | ✅ primary | ❌ |
| Headless build for CI | ❌ | ✅ |
| Size / layout / warning gate | ❌ | ✅ |
| Project state on disk | `__vm/`, `*.vcxproj.user` (gitignored) | `.pio/` (add to `.gitignore`) |

Design rules to keep them from interfering:

- **Single source of truth for sources.** Both drivers compile the *same* `.ino` + `core/` +
  `hardware/` files. Neither generates or rewrites source. `platformio.ini` references the
  existing tree; it does not fork a copy.
- **No shared generated dirs.** Add `/.pio/` to `.gitignore` alongside the existing `__vm/`
  entry. PlatformIO never writes into `__vm/` and VMicro never writes into `.pio/`.
- **Include-path parity.** The committed VMicro project puts **`core/`, `effects/`, and
  `hardware/`** on the include path (vcxproj `IncludePath`, §4.1) — note `effects/`, which README
  §11's "`../../core;../../hardware`" omits. `platformio.ini` reproduces all three via
  `build_flags = -I core -I effects -I hardware`, so a file that compiles under one resolves the
  same headers under the other.
- **Library parity.** Both depend on `FastLED` (+ `SPI`, bundled with the core). PlatformIO
  pins `FastLED` in `lib_deps`; the pinned version should match what the developer has installed
  for VMicro to avoid "builds in CI, not locally" surprises (documented, not enforced).
- **Build-option parity (most consequential — flash size hinges on it).** VMicro bakes in Teensy
  "menu" selections that PlatformIO otherwise sets from its *own* defaults; the CI size numbers are
  only comparable to the flashed image if these match. **These have now been captured from the
  committed VMicro project** — see §4.1; `platformio.ini` must pin every row of that table.

### 4.1 Captured VMicro configuration (source of truth for parity)

Read from the committed Holosphere VMicro project — `targets/Holosphere/Holosphere.vcxproj`
(`UserProperties`), `__vm/Compile.vmps.xml` / `Configuration.Release.vmps.xml`, and
`__vm/.Holosphere.vsarduino.h`. This is the authoritative set the size gate must reproduce:

| Option | VMicro value | Evidence | PlatformIO encoding |
|---|---|---|---|
| Board | **Teensy 4.0** (`teensy40`, `Platform=teensy4`) | `.vsarduino.h` line 9; vmps `Platform` | `board = teensy40` |
| Core / Teensyduino | **1.59.0** (`TEENSYDUINO=159`) | vcxproj include path `.../avr/1.59.0/...`; vmps | pinned via `platform = teensy@<ver mapping to TD 1.59>` |
| GCC toolchain | **arm-none-eabi-gcc 11.3.1** | vcxproj `RemoteCppCompileToolExe` `.../teensy-compile/11.3.1/...` | comes with the pinned platform — verify it resolves to 11.3.1 |
| **Optimization** | **`-O3` ("Fastest", no LTO)** — owner-confirmed (`o3std`) | vcxproj `custom_teensy40_opt="o3std"` (the explicit override; last-build state still showed the `o2std`/`-O2` default) | `build_unflags = <PIO/Teensy default -O*>` **+** `build_flags = -O3`; no LTO |
| CPU speed | **600 MHz** (`F_CPU=600000000`) — board default | vcxproj `PreprocessorDefinitions` `F_CPU=600000000` | `board_build.f_cpu = 600000000L` |
| USB type | **Serial** (`USB_SERIAL`) — board default, no override | vcxproj `PreprocessorDefinitions` `USB_SERIAL` | `-D USB_SERIAL` (Teensy default; pin it explicitly) |
| Keyboard layout | **US English** (`LAYOUT_US_ENGLISH`) — default | vcxproj `PreprocessorDefinitions` `LAYOUT_US_ENGLISH` | `-D LAYOUT_US_ENGLISH` |
| CPU/FPU flags | `-mcpu=cortex-m7 -mfloat-abi=hard -mfpu=fpv5-d16 -mthumb` (`v7e-m+dp/hard`) | vmps; include path `thumb/v7e-m+dp/hard` | set by `board = teensy40` (verify) |
| C++ dialect | **`gnu++20`** (C: `gnu11`) | vcxproj `CppLanguageStandard=gnu++20`, `CLanguageStandard=gnu11` | Teensy core default is `gnu++17` — may need `build_unflags=-std=gnu++17` + `build_flags=-std=gnu++20` |
| Exceptions / RTTI | **OFF** — `-fno-exceptions -fno-rtti` | vmps default **and** explicit user flag | confirm parity (matters for the `std::nothrow` OOM path in `Phantasm.ino`) |
| Extra user C++ flags | `-fno-threadsafe-statics -fno-unwind-tables -fno-asynchronous-unwind-tables -fno-use-cxa-atexit` | vcxproj `VM_ADDITIONAL_COMPILER_CPP_FLAGS` | add verbatim to `build_flags` — they shrink code and must match |
| Warning suppressions | **`-Wno-psabi -Wno-deprecated -Wno-attributes`** | vcxproj `VM_ADDITIONAL_COMPILER_CPP_FLAGS` | add to `build_flags`; feeds §7.2 (these are already-tolerated warnings) |
| Include dirs | **`core/`, `effects/`, `hardware/`** (+ sketch dir + repo root via `-I …/../../`) | vcxproj `IncludePath`; user flag `-I{build.source.path}/../../` | `-I core -I effects -I hardware` — **`effects/` is parity; the repo root (via `src_dir=.`) is the real requirement**, see §6 |
| Core/lib optimization | `OptimiseLibs=True`, `OptimiseCore=True` | vmps `Compile` attrs | PlatformIO optimizes core/libs with the env flags by default |
| Linker script / FlexRAM split | Teensy core's `.ld` (ITCM↔DTCM FlexRAM bank config) | from `board = teensy40` + core 1.59 | **automatic** (same core ⇒ same `.ld`); see note below — load-bearing for the DTCM/RAM1 budget |

**FlexRAM / linker-script parity (the mechanism behind the RAM1 budget matching).** The Teensy 4
RAM1 (512 KiB) is split between ITCM (code) and DTCM (data + stack) in **32 KiB FlexRAM banks** by
the Teensy core's linker script. Because both VMicro and PlatformIO build with the *same* core
(1.59.0), they get the *same* `.ld` and the *same* split rule — so DTCM size (hence arena + stack
headroom, §7.4 #1/#4) matches the flashed image without any extra pinning. Two consequences worth
naming: it is parity we get for free *only as long as the core pin holds* (a core bump can change
the `.ld`), and the split is **quantized to 32 KiB** — when ITCM crosses a bank boundary, DTCM-free
jumps by a whole 32 KiB step. That's fine for the gate (it still fails when `dtcm_free_min_bytes`
is breached), but explains why that floor can move in 32 KiB increments rather than smoothly.

**Two corrections this capture forces on earlier drafts:**

1. **Optimization is `-O3` ("Fastest", no LTO) — owner-confirmed.** The committed override
   `custom_teensy40_opt="o3std"` (`-O3`) is authoritative; the `o2std`/`-O2` seen in the last-build
   vmps state was the stale default, and the earlier `-Os` recollection did not hold (no project
   file selects `osstd`). `-O3` is the **largest** image of the three (`-Os < -O2 < -O3`), so:
   - Calibrate `flash_max_bytes` (§8) against the **-O3** build, and
   - CI must `build_unflags` PlatformIO's Teensy default `-O*` and force `-O3` (no LTO) — the PIO
     default is unlikely to be `-O3`, so without this the gate would measure a *smaller* image than
     the bench flashes.
   - **Relief valve:** if Phase-0 calibration shows the `-O3` image doesn't fit (flash or OCRAM),
     dropping the bench to `-Os` is the obvious shrink — but that's an owner build-config change,
     mirrored in `platformio.ini`, not something the gate decides. For now the gate targets `-O3`.
2. **The repo root must be on the include path; `effects/` is kept for VMicro parity.** Correction
   to earlier drafts: the effect headers are *not* included by bare name — `core/effects.h` pulls
   them PATH-PREFIXED (`#include "effects/BZReactionDiffusion.h"`), and those headers in turn include
   `core/...` and `effects/...` siblings the same way. Those path-prefixed includes resolve from the
   **repo root**, which PlatformIO puts on the include path automatically via `src_dir = .` — *that*
   is the load-bearing requirement, not `-I effects`. `-I core` / `-I hardware` are still needed for
   the BARE `effects.h` / `effects_legacy.h` / `pov_single.h` includes in the `.ino`. `-I effects` is
   retained to mirror VMicro's `../../effects` (vcxproj line 67) and guard a future bare-name include,
   but the current tree does not require it (§6). README §11's include-dir hint is updated to add
   `../../effects` for VMicro parity regardless.

> **Scope of the capture.** Only **Holosphere** has a committed VMicro/VS project; **Phantasm has
> none** to read. The table above is therefore Holosphere's actual config; Phantasm is assumed to
> use the **same** board/optimization/USB/toolchain (same author, same bench) and only adds
> `-D USE_DMA_LEDS`. Phase 0 must confirm that assumption with the owner before locking Phantasm's
> flash budget — it is the one remaining unread option set.

---

## 5. Source layout mapping

PlatformIO's default convention is `src/` with `.cpp`. This repo keeps per-target `.ino` entry
points under `targets/`. We map without moving files using per-environment `src_dir`/`src_filter`
(or `build_src_dir`):

Both targets are **Teensy 4.0** (project-owner confirmed; the 2 MB flash figure follows from
that — see the board note below). They differ in output path and design resolution, not board.
**The gate builds whatever each `.ino` actually commits** — Holosphere's design resolution is
96×20 (README §1) but `Holosphere.ino` is currently set to 288×144 top to bottom (a WIP state);
Phantasm is 288×144 and selects the DMA HD107S path (`#define USE_DMA_LEDS`). The size budgets key
on the committed image, not the design figure.

> **Phantasm is one image shared across four boards.** `Phantasm.ino` reads a hardware ID at boot
> ([lines 11-23](../targets/Phantasm/Phantasm.ino)) and the *same* binary runs on all four
> ID-selected Teensys (a 4-node segmented display, 72 LEDs each). For this compile/size gate that
> is a single image to build and budget; the "288×144" is the assembled virtual display, not one
> board's output.

> **Static footprints are *mostly* shared, but RAM2 is not identical.** The framebuffers
> ([memory.cpp:139-141](../core/memory.cpp)) and the arena are sized by compile-time constants
> (`MAX_W=288`, `MAX_H=144`, `GLOBAL_ARENA_SIZE=335 KB`), **not** by virtual resolution — so both
> allocate the same two 243 KiB buffers in OCRAM and the same 335 KB arena in DTCM. But Phantasm's
> `USE_DMA_LEDS` path adds an OCRAM consumer Holosphere lacks: the double-buffered
> `DMAMEM DMALEDController` eDMA TX frame buffers ([pov_segmented.h:802-810](../hardware/pov_segmented.h)),
> where Holosphere's FastLED path uses a non-`DMAMEM` `CRGB leds_[]`. So the targets diverge in
> **flash** (different effects + USB/driver code) **and** in **RAM2/OCRAM** (Phantasm's DMA TX
> buffers) — which, given OCRAM's structural tightness (§8), means **Phantasm needs its own RAM2
> ceiling**, not a shared one. They do share the RAM1/DTCM (arena-dominated) footprint.

```
targets/
├── Holosphere/Holosphere.ino   ← env:holosphere  (board = teensy40; committed 288×144, design 96×20)
├── Phantasm/Phantasm.ino       ← env:phantasm    (board = teensy40, 288×144, USE_DMA_LEDS, 4 boards)
└── wasm/…                       ← untouched (CMake/Emscripten)
core/        memory.cpp, reaction_graph.cpp, *.h   ← compiled into both (lib or extra source)
hardware/    *.h                                    ← header-only, on include path
```

Key decisions:

- **Keep `.ino`, don't convert to `.cpp`.** PlatformIO compiles `.ino` (it runs the same
  prototype-injection preprocessing Arduino does). Converting to `.cpp` would diverge from the
  VMicro source of truth and risk auto-prototype differences. The `.ino` stays canonical.
- **`core/memory.cpp` and `core/reaction_graph.cpp` are real translation units** (the WASM build
  lists them explicitly in `CMakeLists.txt`). They must be added to the firmware build's source
  set — via `build_src_filter` including `core/*.cpp`, or by treating `core/` as a private
  library. The rest of `core/` and all of `hardware/` are header-only and reached via `-I`.
- **`effects_legacy.h` is in scope here even though it's excluded from code review.**
  `Holosphere.ino` includes it, so the firmware image must compile it. (The review-scope
  exclusion in `prompts/analysis.txt` is about *quality grading*, not *buildability*.)
- **Same board, different flags.** Both envs set `board = teensy40`. Phantasm adds
  `build_flags += -D USE_DMA_LEDS` to mirror the `#define USE_DMA_LEDS` the `.ino` already sets
  (keep the `#define` in the `.ino` as the source of truth; the env flag is only needed if a
  target ever relies on it being set before the first include). The two builds share the same
  board and RAM1 footprint but differ in flash *and* in RAM2 (Phantasm's DMA TX buffers, above) —
  so the budgets (§8) share RAM1 but give Phantasm its own flash and RAM2 ceilings.

> **Board note — confirmed Teensy 4.0, committed docs are stale (must reconcile, not "later").**
> The board pin is the single load-bearing assumption of the flash gate: T4.0 flash is **2 MB**,
> T4.1 is **8 MB**, so a wrong board makes `flash_max_bytes` meaningless and could false-fail a
> valid image on day one. The project owner has **confirmed both targets are Teensy 4.0** (there is
> no 4.1), so the gate pins `board = teensy40` and budgets flash against 2 MB. But `Phantasm.ino`'s
> header ([line 7](../targets/Phantasm/Phantasm.ino)) and README §1 still say *4× Teensy 4.1* —
> those committed lines are now **confirmed stale** and reconciling them to 4.0 is a **prerequisite
> of Phase 0/1**, not optional tidy-up, precisely because the budget depends on the board being
> unambiguous in the tree. Separately, `Holosphere.ino` is committed at 288×144 (header *and* body)
> while its design resolution is 96×20 — the gate budgets the committed image, and the resolution
> reconciliation is a parallel cleanup (§14 decision 5).

---

## 6. `platformio.ini` design (illustrative, not final)

```ini
[platformio]
src_dir = .                           # repo root — sources live in core/ and targets/, no src/
# Keep PlatformIO's tree out of the repo root; gitignored.
# Sources are referenced in place — no files move.

[env]
platform = teensy@<pinned>            # pin to the build giving Teensy core 1.59.0 / arm-gcc 11.3.1 (§4.1)
framework = arduino
board_build.f_cpu = 600000000L        # 600 MHz (§4.1)
build_cache_dir = .pio/build_cache    # object cache (mirrors ci.yml's ccache) — gitignored, CI-cached (§10)
lib_ldf_mode = chain                  # constrain the Library Dependency Finder; confirm scope in Phase-0 spike (§6)
build_unflags = -Os                   # drop PlatformIO/Teensy default so the -O3 below wins (§4.1)
build_flags =
    -O3                               # owner-confirmed "Fastest" (o3std), no LTO (§4.1)
    -I core
    -I effects                        # VMicro parity; repo-root (src_dir=.) resolves effects/*.h (§4.1)
    -I hardware
    -D USB_SERIAL                     # USB type (§4.1)
    -D LAYOUT_US_ENGLISH              # keyboard layout (§4.1)
    -std=gnu++20                      # match VMicro CppLanguageStandard (Teensy default is gnu++17)
    ; VMicro user flags, copied verbatim so code size matches the bench (§4.1):
    -fno-exceptions -fno-rtti -fno-threadsafe-statics
    -fno-unwind-tables -fno-asynchronous-unwind-tables -fno-use-cxa-atexit
    ; VMicro warning suppressions — these are already-tolerated; keep them so the §7.2
    ; baseline isn't polluted by warnings the bench build never shows:
    -Wno-psabi -Wno-deprecated -Wno-attributes
    -Wall -Wextra                     # warning hygiene (baseline ratchet, not -Werror — §7.2)
lib_deps =
    fastled/FastLED@<pinned>
# src_dir is the repo root, so the filter must START by excluding everything and
# then add back ONLY the wanted TUs — otherwise targets/wasm/, build*/ , .pio/,
# obj/, out/ etc. get swept in. Both targets compile the two real core TUs:
build_src_filter =
    -<*>
    +<core/memory.cpp>
    +<core/reaction_graph.cpp>

[env:holosphere]
board = teensy40                       # 96×20
build_src_filter = ${env.build_src_filter} +<targets/Holosphere/Holosphere.ino>

[env:phantasm]
board = teensy40                       # 288×144, DMA HD107S path
build_flags = ${env.build_flags} -D USE_DMA_LEDS
build_src_filter = ${env.build_src_filter} +<targets/Phantasm/Phantasm.ino>
```

> **Load-bearing mechanic — de-risk in Phase 0.** This is the one part most likely to misbehave on
> the first try, so it is **not** "just an implementation detail." `build_src_filter` resolves
> relative to a single `src_dir`; the repo has sources in both `core/` and `targets/<X>/` and no
> `src/`, which forces `src_dir = .` (repo root). At repo root the filter must default-exclude
> (`-<*>`) and add back only the exact wanted TUs, or it will pull in `targets/wasm/wasm.cpp`,
> the CMake `build*/` trees, `.pio/`, `obj/`, `out/`, and FastLED sources. Phase 0 (§13) must
> include a spike that confirms each env compiles *exactly* its `.ino` + the two `core/*.cpp` TUs
> and nothing else (verify via the build's compiled-source list), and pins down the precise
> exclude set. The spike must also: (a) confirm the option **name** for the pinned PlatformIO Core
> version — the key was `src_filter` in older PIO and `build_src_filter` in newer releases, use
> whichever the pin accepts; and (b) confirm the **Library Dependency Finder** (`lib_ldf_mode`)
> does not wander the repo root (`targets/wasm/`, `build*/`, `.pio/`) — with `src_dir = .` the LDF
> can scan far more than intended, a perf and surprise-dependency risk adjacent to the filter
> hazard. The design contract otherwise stands: **two named environments, one real device board,
> the existing sources in place, pinned everything.**

### Version pinning (reproducibility contract)
Mirror `ci.yml`'s discipline (it pins clang-18, emsdk 5.0.0 by commit SHA, ubuntu-24.04). Pin:
- `platform = teensy@…` — **pin by git URL + commit SHA, not a version tag.** `ci.yml` pins emsdk
  by commit SHA precisely because a tag is mutable/re-pointable; a `teensy@X.Y.Z` tag is the weaker
  bar the spec elsewhere criticizes. Use `platform = https://github.com/platformio/platform-teensy.git#<sha>`
  (or a vendored fork) so the toolchain can't shift under a retagged release.
- `framework-arduinoteensy` package version (via the platform pin),
- `fastled/FastLED@X.Y.Z` (pin to a tag at minimum; SHA if you want full parity with the bar above),
- the PlatformIO Core version itself (install a fixed `platformio==X.Y` in the CI step),
- the runner image (`ubuntu-24.04`, not `-latest`).

A toolchain bump is then a deliberate, reviewable one-line change — and because the size budgets
are toolchain-sensitive (a new GCC changes codegen size), the pin is what makes the size gate
meaningful frame-over-frame.

> **⚠ Phase-0 gating spike (do this FIRST — it's go/no-go, not a verification line item).** The
> whole "green PlatformIO ≈ healthy VMicro" argument assumes PlatformIO can actually deliver the
> **bench toolchain: Teensyduino 1.59.0 + arm-gcc 11.3.1** (§4.1). The upstream `platform-teensy`
> package has historically *lagged* Teensyduino releases, so there may be **no published platform
> version that yields exactly TD 1.59 / gcc 11.3.1.** Resolve this before anything downstream,
> because budgets and the warning baseline are toolchain-sensitive. Decision tree with a written
> fallback:
> 1. **A pinnable platform commit yields TD 1.59 / gcc 11.3.1** → pin it (above); done.
> 2. **It doesn't** → override the toolchain via `platform_packages` pointing at the
>    Teensyduino-provided `toolchain-gccarmnoneeabi` + framework (or vendor a platform fork) to
>    force 11.3.1.
> 3. **Neither is practical** → accept a *different but pinned* toolchain, and **recalibrate**
>    budgets + warning baseline against it — explicitly noting that CI then approximates (not
>    mirrors) the bench compiler. This is acceptable for a size/link gate but must be stated, not
>    silent.
>
> Everything else in this spec is downstream of this spike: if none of the three branches is
> acceptable, the *driver* (PlatformIO) is reconsidered, not just a number.

---

## 7. Validation dimensions

### 7.1 Build success
`pio run -e holosphere -e phantasm`. Any compile or link failure fails the job. This alone closes
the largest part of the gap: device-only `#ifdef ARDUINO` code now has a compiler pointed at it
on every push.

### 7.2 Warning hygiene — baseline ratchet *(decided)*
Policy: **baseline ratchet**, not hard `-Werror`. The tree is not yet known to be warning-clean
under arm-none-eabi-gcc, and a `-Werror` flip would either block landing the build gate or force a
large up-front cleanup. The ratchet lets the build gate land now while still preventing *new*
warnings.

- Compile with `-Wall -Wextra` (defer `-Wshadow`/`-Wconversion` — likely noisy; revisit once the
  baseline is clean).
- **Capture a committed baseline as an unordered, deduplicated set.** Store a fingerprint of the
  current first-party warning set (e.g. `tools/teensy_warning_baseline.txt`: normalized `warning:`
  lines with volatile bits — absolute paths, line numbers, column — stripped so the fingerprint is
  stable across unrelated edits). **The comparison must be set-based (sorted, de-duplicated), not a
  line-ordered diff:** PlatformIO builds in parallel (`-j`), so warning *emission order* is
  nondeterministic run-to-run — an ordered diff would flap green/red on identical inputs. The gate
  compares the build's warning *set* against the baseline *set* and **fails only on members not in
  the baseline**. A real fix that removes a baseline warning is fine; the gate flags *additions*.
- **Library noise excluded.** FastLED and the Teensy core emit their own warnings. The plan is
  `-isystem` for vendored/core includes so their warnings never enter the set — **but demoting the
  Teensy core to `-isystem` under PlatformIO is non-trivial (the core is added by the framework,
  not by us), so make "confirm the core/FastLED can actually be `-isystem`'d" a Phase-0 spike**, not
  an assumed capability. The first-party path filter (keep only `core/`, `effects/`, `hardware/`,
  `targets/` warnings) is an independent backstop, so the design is robust even if the `-isystem`
  demotion proves awkward.
- **The warning gate must run on a *clean* (cache-disabled) build — it cannot share the cached size
  build.** The object cache (§10, `build_cache_dir`) suppresses compiler invocation on a hit, so a
  cached TU emits **no** warnings: a warm build produces a *smaller* warning set than a cold one,
  and a new warning introduced in a **header** is invisible whenever the dependent TU is served from
  cache. Set-based comparison fixes emission *order*, not this. Resolution: the **size** build stays
  cached (objects don't change between warm runs), but the **warning** capture runs a separate
  build with the cache disabled (e.g. a dedicated env or `PLATFORMIO_BUILD_CACHE_DIR=` empty / a
  forced clean) so every first-party TU re-emits its warnings every run. Decouple the two; don't let
  the cache that makes the size job fast silently blind the warning ratchet. (If a full clean build
  is too slow, the documented fallback is scoping the ratchet to changed TUs — but then the
  header-introduced-warning blind spot must be stated, not hidden.)
- **`effects_legacy.h` warnings will be in the baseline.** It is review-excluded but in *build*
  scope (`Holosphere.ino` includes it, §5), so whatever it warns under arm-gcc lands in the
  committed first-party baseline. That's correct (the firmware does compile it); just don't be
  surprised to see review-excluded code in `teensy_warning_baseline.txt`.
- **Updating the baseline is a reviewed act.** A `--update-baseline` script regenerates the file;
  it lands in the same PR as the change that legitimately alters the warning set, so the diff is
  visible in review (same discipline as the LUT/reaction-graph provenance gates).
- **Phasing:** ship report-only for a PR or two to capture a stable baseline, then enforce
  (Phase 2, §13).

### 7.3 Image size
Teensy's post-build `teensy_size` (run automatically by the Teensy platform) prints per-region
usage, e.g.:

```
teensy_size: FLASH: code:NNNNN, data:NNNN, headers:NNNN   free for files: NNNNNN
teensy_size: RAM1: variables:NNNNN, code:NNNNN, padding:NNNNN   free for local variables: NNNNN
teensy_size: RAM2: variables:NNNNN   free for malloc/new: NNNNNN
```

- **RAM1** = ITCM (fast code) + DTCM (fast `.bss`/`.data` + stack). **The 335 KB arena
  (`global_arena_block`) lives here** — it is a plain static array with no `DMAMEM` qualifier
  ([memory.cpp:43-53](../core/memory.cpp)) and is the dominant RAM1 occupant. The stack also grows
  down within DTCM.
- **RAM2** = OCRAM: `DMAMEM` variables — the two 243 KiB framebuffers and the timeline event
  buffer — plus the heap (`malloc`/`new`: the `POVSegmented`/`POVDisplay` driver objects). The
  arena is **not** here, and is not heap-allocated.
- **FLASH** = code + `const`/`PROGMEM`/`FLASHMEM` data (incl. the 90 KB reaction-graph table) +
  headers.

The gate parses these and enforces a **budget per region per target** (§8). Thresholds carry
headroom — failing *before* 100% so there's room to flash and to grow the stack, not at the
linker's hard wall.

### 7.4 Memory layout
Raw totals miss structural regressions. Validate **invariants** by inspecting section/symbol
placement from the ELF. **Classify each symbol by its load address against the Teensy 4 memory
map, not by an `nm` type letter** — DTCM `.bss` and OCRAM `.dmabuffers` are *both* NOBITS, so `nm`
reports both as `b`/`B` and cannot tell a DTCM buffer from an OCRAM one (it would miss the exact
DMAMEM→DTCM regression invariant #2 targets). Use `arm-none-eabi-readelf -s`/`-S` (symbol →
section), or each symbol's address bucketed by region — ITCM `0x0000_0000`, FLASH `0x6000_0000`,
DTCM `0x2000_0000`, OCRAM `0x2020_0000` — cross-checked against the `.map`. `size -A` still gives
per-section totals; `nm` alone is insufficient for the placement checks below.

**Match symbols by their mangled / linkage names — the parser must not assume source spellings.**
`Effect::buffer_a` is a class static, so it appears mangled (`_ZN6Effect8buffer_aE`);
`global_arena_block` is `static` (internal linkage) and may show as a local symbol. The classifier
must either demangle (`arm-none-eabi-c++filt`, or `readelf`'s demangle) or match the mangled form
directly, and the §9.1 golden fixtures must contain **real `readelf -s` output with the actual
mangled strings** so the tests exercise the real matching, not a tidied alias. This is an easy spot
to be silently wrong (a name that never matches → the invariant never fires → false-green).

1. **Arena in DTCM/RAM1, and its size ≈ 335 KB (assert the magnitude, not just the section).**
   `global_arena_block` must stay in DTCM (`.bss`, no `DMAMEM`) **and** measure ~`GLOBAL_ARENA_SIZE`
   = 335 KB ([memory.h:33](../core/memory.h)). Pinning the *magnitude* matters: under
   `HS_TEST_BUILD` the same constant is **8 MB** ([memory.h:31](../core/memory.h)), so if that
   test-only macro ever leaked into the firmware build the arena would silently balloon 24× — a
   "still in DTCM" check passes, a "~335 KB ± tolerance" check catches it. It is the **largest RAM1
   object and the most likely RAM1-budget mover** — a change to `GLOBAL_ARENA_SIZE`, or accidentally
   tagging the block `DMAMEM` (which would shove 335 KB into already-tight OCRAM), shows up here too.
2. **Framebuffers in OCRAM, not DTCM.** `buffer_a`/`buffer_b` (the `DMAMEM Pixel[MAX_W*MAX_H]`
   arrays) must land in the OCRAM/`.dmabuffers` section. If a refactor drops the `DMAMEM`
   qualifier they'd silently move into DTCM and blow the fast-RAM/stack budget — a layout bug a
   size *total* might not catch but an *attribution* check will.
3. **Reaction-graph table in flash, not RAM.** The 90 KB `reaction_graph.cpp` table must be in a
   flash/`.rodata` section; a missing `const` qualifier would pull it into RAM and is an instant
   90 KB regression. Key the check on the **section** (flash/`.rodata`), *not* on a literal
   `PROGMEM` marker: on the Teensy 4's memory-mapped flash `PROGMEM` is essentially decorative —
   it's the `const` that keeps the table out of RAM ([reaction_graph.cpp:4](../core/reaction_graph.cpp)).
4. **DTCM stack headroom (link-time availability, not runtime depth).** Assert RAM1 "free for local
   variables" (= DTCM minus the arena and other `.bss`/`.data`) stays above a floor. This measures
   *how much room the linker left for the stack*, **not** actual runtime stack depth — it is not
   stack-overflow protection, it just guarantees the static allocations didn't crowd the stack space
   out. The stack grows down into this region at runtime. Mirrors the WASM build's stack-size
   vigilance (`-sSTACK_SIZE`, the debug 64 KB override note in `CMakeLists.txt`).
5. **(Optional) named-symbol size caps** for the largest known objects, so an unexpected balloon
   in a specific buffer is attributed, not just summed.

---

## 8. Budget / threshold model

A checked-in, reviewable budget file (e.g. `tools/teensy_budgets.json`) — analogous in spirit to
the committed LUT/reaction-graph provenance the existing CI guards:

Both targets are Teensy 4.0, so both share the same hard region sizes (2 MiB FLASH, 512 KiB RAM1
ITCM+DTCM, 512 KiB RAM2 OCRAM). The framebuffers and arena are sized by compile-time constants
identically for both (§5), so the two targets share their **RAM1/DTCM** footprint. They diverge in
**FLASH** (different effects + DMA vs FastLED path) and, *in principle*, in **RAM2/OCRAM** —
Phantasm's `USE_DMA_LEDS` path adds the `DMAMEM DMALEDController` eDMA TX buffers Holosphere lacks
(§5).

> **Right-sizing: the OCRAM *layout invariant* is the real value; a *separate numeric* RAM2 ceiling
> is a Phase-0 call, not an a-priori requirement.** The DMA TX delta is the double-buffered frame
> for ~72 LEDs/segment — order of a few hundred bytes against OCRAM's ~26 KiB margin. So the thing
> that actually protects OCRAM is the **`dma_tx_buffers_section: OCRAM` layout check** (a dropped
> `DMAMEM` there is a correctness bug regardless of bytes), not a second hand-maintained number.
> **Default to a single shared OCRAM ceiling + the layout invariant**, and only split into a
> separate, tighter Phantasm `ram2_max_bytes` if Phase-0's recorded as-built delta proves it earns
> the extra budget to maintain. The JSON below shows the *split* form for completeness; collapsing
> `ram2_max_bytes` to one shared value is the expected outcome.

```jsonc
{
  "holosphere": {
    "board": "teensy40",
    "flash_max_bytes":  1900000,   // < 2 MiB T4.0 flash, with headroom
    "ram1_max_bytes":   480000,    // DTCM dominated by the ~335 KB arena (+ stack, .bss) (of 512 KiB)
    "ram2_max_bytes":   510000,    // OCRAM: 2×243 KiB buffers = 497,664 B ALONE — see tightness note
    "dtcm_free_min_bytes": 16384,  // stack headroom floor (DTCM after arena/.bss)
    "layout": {
      "arena_section": "DTCM",
      "framebuffers_section": "OCRAM",
      "reaction_graph_section": "FLASH"
    }
  },
  "phantasm": {
    "board": "teensy40",
    "flash_max_bytes":  1900000,   // calibrate separately — different effect set + DMA/USB code
    "ram1_max_bytes":   480000,    // shared with Holosphere (same arena/buffers)
    "ram2_max_bytes":   "<tighter>", // OCRAM: framebuffers + the DMALEDController TX buffers — its own cap
    "dtcm_free_min_bytes": 16384,
    "layout": {
      "arena_section": "DTCM",
      "framebuffers_section": "OCRAM",
      "dma_tx_buffers_section": "OCRAM",   // DMALEDController must stay DMA-reachable (OCRAM), see §5
      "reaction_graph_section": "FLASH"
    }
  }
}
```

Two callouts the calibration in Phase 0 must respect:

- **The arena dominates RAM1.** `ram1_max_bytes` and `dtcm_free_min_bytes` are governed by
  `GLOBAL_ARENA_SIZE` (335 KB) plus the stack; that constant is the most likely thing to move the
  RAM1 budget. The layout check (§7.4 #1) pins it in DTCM *and* near 335 KB. Note `dtcm_free_min_bytes`
  can step by a whole **32 KiB** when ITCM crosses a FlexRAM bank boundary (§4.1 FlexRAM note) — a
  jump in the floor, not smooth drift; the gate still fires correctly, it's just lumpy.
- **OCRAM headroom is structurally small, not a free parameter.** The two static framebuffers are
  **497,664 B of OCRAM's 524,288 B** — leaving only ~**26 KiB** for the timeline event buffer
  (`global_timeline_events`, `TIMELINE_MAX_EVENTS=64` — small, ~a few KiB) and all `new`/`malloc`
  (the driver objects; Phantasm's `DMALEDController` TX buffers also sit here). The known consumers
  are small, so the **current image very likely fits today** — but the margin is thin enough that
  the gate's value is real. A placeholder like `ram2_max_bytes: 500000` is *barely above the buffers
  alone* — that is not "headroom." **Phase 0 must record the as-built OCRAM-free number** so the
  framing is explicit: if Phantasm is already near the wall, the gate's first act is *surfacing
  existing tightness*, not just *preventing future* growth. Either way the buffers are fixed at
  `MAX_W*MAX_H`, so the only give is the heap.

**Enforcement: absolute ceilings only** *(decided)*. The gate hard-fails if any region exceeds its
cap — catching "won't fit / no room to flash" and the layout invariants (§7.4). A per-build
**regression-delta** check (fail on growth beyond a tolerance vs a committed size baseline) was
*considered and deferred*: it adds a second baseline to maintain and would make routine effect
work churn the baseline file. Ceilings with deliberate headroom give the needed protection without
that overhead; the delta gate can be revisited later if size creep proves to be a real problem.

**Raising a ceiling is a reviewed, one-line edit** to `tools/teensy_budgets.json`, landed in the
*same* PR as the change that needs it — the symmetric escape hatch to the warning baseline's
`--update-baseline` (§7.2). A budget bump is therefore visible in review (you can see the headroom
shrink) rather than silent; provide a small `--show`/`--update` helper that re-reads the current
build and prints the would-be values to make the edit mechanical.

Headroom numbers above are placeholders — calibrate from the first real build's `teensy_size`
output before enabling enforcement.

---

## 9. Measurement mechanism

A small **PlatformIO `extra_scripts` post-build hook** (Python) is the cleanest integration point:

1. After link, locate the ELF (`$BUILD_DIR/firmware.elf`).
2. Run `teensy_size` (or `arm-none-eabi-size -A`) and parse region totals.
3. Run `arm-none-eabi-readelf -s -S` (and/or parse the `.map`) and classify the framebuffer,
   arena, DMA-TX, and reaction-graph symbols **by load address against the Teensy 4 memory map**
   (§7.4) — DTCM/OCRAM are both NOBITS so `nm` type letters can't separate them. This is what makes
   the placement invariants (arena→DTCM, framebuffers/DMA-TX→OCRAM, table→FLASH) actually
   detectable.
4. Compare against `tools/teensy_budgets.json`; emit a clear pass/fail report.
5. **Fail the build on violation.** Attach the hook via
   `env.AddPostAction("$BUILD_DIR/${PROGNAME}.elf", check_fn)` and have `check_fn` call
   `sys.exit(1)` (or `raise`) on any violation — a post-action that merely prints does **not** fail
   `pio run`; the non-zero exit / exception is what propagates. Emit GitHub Actions annotations
   (`::error::…`) first so violations surface inline on the PR, matching `ci.yml` convention.

Keeping the parser in a versioned script (not inline YAML) means it is unit-testable and reused
verbatim by the local pre-flight command (§11). The toolchain binaries (`teensy_size`,
`arm-none-eabi-readelf`, `-size`) ship inside the pinned PlatformIO Teensy package, so no extra
install.

### 9.1 The gate must gate itself — negative/golden tests *(highest-value deliverable)*

A size/layout check that can never fail is **worse than none**: it renders green while covering
nothing, and a subtle address-bucket bug yields permanent false-green. So the parser and *every*
layout invariant ship with tests that prove each one **fails when it should** — this is a required
deliverable of Phase 0/1, not an optional nicety:

- **Golden fixtures.** Commit captured `arm-none-eabi-size -A` / `readelf -s -S` output from a known
  build under `tools/teensy_gate_tests/`, and assert the parser extracts the expected region totals
  and symbol→section classifications from them (no toolchain needed to run the parser tests).
- **Deliberately-broken cases — one per invariant, each must turn the check red:**
  - framebuffer fixture with `DMAMEM` dropped → symbol lands in DTCM → **framebuffers→OCRAM check
    must fail**;
  - reaction-graph fixture with `const` dropped → symbol in RAM → **table→FLASH check must fail**;
  - arena fixture at the 8 MB `HS_TEST_BUILD` size → **arena ~335 KB magnitude check must fail**
    (§7.4 #1);
  - an over-cap totals fixture → **each region ceiling must fail**.
- **Address-bucketing unit tests** for the Teensy 4 memory-map classifier (ITCM/FLASH/DTCM/OCRAM
  boundaries), since that logic is the load-bearing replacement for `nm` and the easiest place for a
  silent off-by-one to hide.

These run as ordinary host Python tests (fast, no ARM toolchain), so they can live in the existing
CI test job or a tiny dedicated step — independent of the slow PlatformIO build.

---

## 10. CI integration

A new job in `.github/workflows/ci.yml` (it already fans out tests/sanitizers/wasm/provenance;
this slots in alongside):

A **single job builds both envs** (no matrix) — they share one Teensy toolchain, so a matrix would
download it twice for marginal isolation. `pio run` builds the envs sequentially and reports per-env
size, and the gate script fails the step on the first env that violates a budget; this is the same
invocation the local recipe uses (§11):

```yaml
  teensy-size:
    name: Teensy 4 firmware (build + size + layout)
    runs-on: ubuntu-24.04
    # Skip the heavy toolchain restore on changes that can't affect the firmware image
    # (docs, the WASM/daydream side, the native tests). Scope to what the image is built
    # from. NOTE: path filters only apply to pull_request/push events; the required-check
    # below must tolerate "skipped" as success on out-of-scope PRs (see Notes).
    # (Apply via the workflow-level on.<event>.paths, or a paths-filter step.)
    #   paths: [ core/**, effects/**, hardware/**, targets/**, platformio.ini,
    #            tools/teensy_*, .github/workflows/ci.yml ]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5        # PlatformIO is a Python package
        with: { python-version: '3.x' }
      - name: Cache PlatformIO toolchain
        uses: actions/cache@v4
        with:
          path: ~/.platformio
          # Target-agnostic (toolchain doesn't vary by env, so both share one entry),
          # but DO include everything that changes ~/.platformio contents: the pinned
          # platform version, the FastLED pin, and the PlatformIO Core version — else a
          # dependency bump silently restores a stale cache.
          key: pio-teensy-<platform-ver>-<fastled-ver>-<pio-core-ver>
          # Optional: hashFiles('platformio.ini') to auto-rotate on any pin change.
      - name: Cache build objects               # mirrors ci.yml's ccache discipline
        uses: actions/cache@v4
        with:
          # build_cache_dir from platformio.ini — without this the 7,685-line
          # reaction_graph.cpp + FastLED recompile from scratch on every run, ×2 envs.
          path: .pio/build_cache
          # github.sha makes each commit write a fresh entry, restoring from the newest
          # prefix match — standard, but it does churn the cache (Actions evicts at the
          # 10 GB repo budget, LRU). Fine for object caches; just don't be surprised.
          key: pio-objs-<platform-ver>-${{ hashFiles('platformio.ini') }}-${{ github.sha }}
          restore-keys: pio-objs-<platform-ver>-
      - name: Install pinned PlatformIO
        run: pip install 'platformio==<pinned>'
      - name: Build + size + layout gate (both targets)
        run: pio run -e holosphere -e phantasm   # extra_script enforces per-env budgets, sets exit code
      - name: Upload firmware artifacts            # optional: .hex/.elf for inspection
        uses: actions/upload-artifact@v4
        with: { name: firmware, path: .pio/build/*/firmware.* }
```

Notes:
- **No hardware, no upload** — build + link + measure only. Runs on a stock Linux runner.
- **Required status check, flipped at Phase 1.** "Fail CI" only blocks merges if `teensy-size` is a
  *required* status check in branch protection. The phasing maps directly: Phase 0 runs it
  report-only (not required), Phase 1 makes it required once budgets are calibrated. When using
  `paths:` scoping, make the required check tolerate a *skipped* run as success (a path-filtered job
  that doesn't run must not block out-of-scope PRs) — e.g. a always-passing companion job, or
  `dorny/paths-filter` gating the heavy step rather than the whole job.
- **Two caches, like `ci.yml`.** `~/.platformio` (toolchain/packages) *and* `build_cache_dir`
  (objects). Without the object cache the big TU + FastLED rebuild cold every run; with it, only
  changed TUs recompile — the same reason `ci.yml` runs ccache.
- **One job, both envs.** Cold-cache cost is one toolchain download, not two. If per-target
  *independent* failure isolation later proves worth the extra runner, fall back to a matrix — but
  keep the cache key target-agnostic (above) so the warm-cache case still shares one entry.
- **Cache `~/.platformio`** (the toolchain + framework packages) keyed on the pinned platform
  version (target-agnostic), mirroring the ccache/emsdk caching already in `ci.yml`. First run
  downloads the Teensy toolchain (~hundreds of MB); subsequent runs restore it.
- This job is **additive** — it does not touch the existing tests/wasm jobs, and PlatformIO needs
  no `EMSDK`, so it composes cleanly with the memory note `project_no_device_build_in_ci` (which
  this gate partially supersedes for *compile/link/size*, while on-device *execution* stays
  manual).

---

## 11. Local parity (for a VMicro developer)

A developer who flashes via VMicro should be able to run the identical gate before pushing. Add a
`justfile` recipe (the repo already routes tasks through `just`):

```
# Build the Teensy firmware images and enforce size/layout budgets (CI parity).
teensy-size:
    pio run -e holosphere -e phantasm
```

`pio` is the only new prerequisite (`pip install platformio`); the Teensy toolchain auto-installs
on first `pio run`. This gives a one-command local pre-flight without disturbing the VMicro project.
The contract is **same pass/fail under the headroom'd ceilings**, *not* identical bytes: the dev
runs PIO on Windows while CI runs `ubuntu-24.04`, and arm-gcc — though largely host-deterministic —
can shift a few bytes via `__FILE__`/path-length or any build-timestamp injection. Same honesty as
the VMicro≠PIO note (§4): the gate asserts *fit*, not *bit-identity*. Document the recipe in README
§11 alongside the existing firmware build instructions.

---

## 12. Failure modes & reporting

| Failure | Signal | Example cause |
|---|---|---|
| Compile error in device path | `pio run` non-zero, compiler diagnostic | broke `#ifdef ARDUINO` code, missing header on `-I` path |
| Link error | linker diagnostic | undefined symbol only referenced on device, region overflow at hard wall |
| New warning | `::error::` from warning gate | baseline-ratchet trip (warning not in the committed baseline) |
| Region over budget | `::error::` from size parser, region + bytes + cap | added a large effect, dropped `FLASHMEM` |
| Layout invariant violated | `::error::` naming the symbol + wrong section | dropped `DMAMEM`/`const`, table moved to RAM |

All violations use the `::error::` annotation convention already established across `ci.yml` so
they render inline on the PR.

---

## 13. Rollout plan (phased)

1. **Phase 0 — spikes + land `platformio.ini` + the post-build script, report-only.** Build both
   targets in CI; print sizes/warnings; **do not fail** on thresholds yet (not a required check
   yet). Phase 0 must resolve the prerequisites that make the numbers meaningful:
   - **(0) FIRST: the toolchain go/no-go spike (§6).** Can PlatformIO deliver TD 1.59 / arm-gcc
     11.3.1, or which fallback (platform_packages override / fork / accept-and-recalibrate)? Settle
     this before any calibration — everything below is toolchain-sensitive, and a no-go could change
     the driver, not just a number.
   - (a) the `build_src_filter` spike (§6) confirming each env compiles exactly its `.ino` + the two
     `core/*.cpp`, **and** the `lib_ldf_mode` scope check (§6);
   - (b) **build-option parity** (§4.1) — already captured (optimization `-O3`, `f_cpu`, USB, layout,
     gnu++20, the `-fno-*`/`-Wno-*` flags); confirm they reproduce VMicro's image, and confirm the
     **`-isystem` demotion** of FastLED/the Teensy core actually works (§7.2);
   - (c) reconcile the stale "Teensy 4.1" references to 4.0 so the board/flash budget is unambiguous
     in the tree (§5 board note);
   - (d) **record the as-built per-region numbers** — especially OCRAM-free (§8) — so the
     regression-prevention-vs-existing-tightness framing is explicit;
   - (e) write the **negative/golden parser tests** (§9.1) so every invariant is proven to fail.

   Then calibrate real budgets (incl. Phantasm's separate RAM2 cap) and capture the initial warning
   baseline (as an unordered set, §7.2). (One or two PRs.)
2. **Phase 1 — enforce build + size (ceilings) + layout, and make the check required.** Set the
   calibrated ceilings and the layout invariants to hard-fail; with the negative tests (§9.1) green,
   flip `teensy-size` to a **required status check** (§10). Add the `just teensy-size` recipe and
   README docs (including the "Teensy 4.0 for both targets" fix, §14 decision 5).
3. **Phase 2 — enforce the warning baseline ratchet.** With the committed baseline in place
   (Phase 0) and `-isystem` excluding FastLED/core, fail on any first-party warning not in the
   baseline set.

Each phase is independently revertable and adds signal without destabilizing the existing green
CI. (A size regression-delta gate was considered and deferred — see §8 — so there is no Phase for
it.)

### 13.1 Acceptance criteria — "the gate is done when…"

1. Both envs (`holosphere`, `phantasm`) build **green** under the pinned toolchain, compiling
   exactly their `.ino` + the two `core/*.cpp` (filter spike confirmed, §6).
2. The build options reproduce the VMicro image: `-O3`, 600 MHz, `USB_SERIAL`, `gnu++20`, and the
   captured `-fno-*`/`-Wno-*` flags (§4.1).
3. All layout invariants — arena→DTCM ≈335 KB, framebuffers/DMA-TX→OCRAM, table→FLASH — are
   **demonstrably fail-then-pass** via the negative/golden fixtures (§9.1).
4. Per-region budgets are calibrated from a real `-O3` build, with the as-built OCRAM-free recorded
   (§8); Phantasm carries its own FLASH and RAM2 ceilings.
5. The warning baseline is captured (unordered set) and the ratchet fails on a synthetic new
   warning but not on reordering (§7.2).
6. `teensy-size` is a **required status check** that tolerates a path-skipped run as success (§10),
   and `just teensy-size` reproduces the CI *pass/fail* locally (§11).
7. **Wall-clock budget met.** Warm-cache run (toolchain + objects restored) completes under a stated
   target (e.g. **< ~3 min**); cold run (toolchain download + full `-O3` rebuild of the 7,685-line
   `reaction_graph.cpp` + FastLED × 2 envs) under a looser one. If warm runs blow the target, the
   object cache (§10) isn't engaging — investigate before flipping the check required.

---

## 14. Resolved decisions

All five prior open questions have been decided; the spec body reflects them.

1. **Targets — both.** Gate `Holosphere` *and* `Phantasm` (each currently committed at 288×144) —
   built in a single job (`pio run -e holosphere -e phantasm`), not a matrix, since both share one
   toolchain (§10). Both are **Teensy 4.0** (owner-confirmed). They share the RAM1/arena footprint
   and differ mainly in **flash** (different effects + DMA/USB code); the small OCRAM delta
   (Phantasm's DMA TX buffers) is guarded primarily by the `dma_tx_buffers_section: OCRAM` **layout
   invariant**, with a *separate* numeric RAM2 ceiling only if Phase-0's as-built delta warrants it
   (default: one shared OCRAM ceiling). Per-target **flash** ceilings stand. (§5, §8, §10)
2. **Warning policy — baseline ratchet.** Not hard `-Werror`. Capture a committed first-party
   warning baseline and fail only on new warnings, with `-isystem` excluding FastLED/the Teensy
   core. (§7.2, Phase 2 in §13)
3. **Size enforcement — absolute ceilings only.** Per-region caps with headroom; no
   regression-delta gate (considered and deferred to avoid a second baseline to churn). (§8)
4. **Workflow placement — inside `ci.yml`.** A new single `teensy-size` job (both envs in one
   `pio run`) alongside the existing tests/wasm/provenance jobs; shared triggers/concurrency, no
   separate workflow file. (§10)
5. **Fix the docs — agreed, and a Phase-0 prerequisite for the board pin.** The owner confirms
   **Teensy 4.0** for both (no 4.1). The stale "4× Teensy 4.1" lines in `Phantasm.ino:7` and
   README §1 must be reconciled to 4.0 **before** the flash budget is locked — the board fixes the
   flash ceiling (2 MB vs 8 MB), so an unambiguous tree is a prerequisite, not downstream tidy-up
   (§5 board note, Phase 0 in §13). The `Holosphere.ino` 288×144-vs-96×20 resolution mismatch is a
   parallel cleanup; the gate budgets the committed image regardless. The board lines **have now
   been reconciled to 4.0** (`Phantasm.ino:7`, README §1 and §11); the `Holosphere.ino` resolution
   mismatch remains the parallel cleanup.

---

## 15. Risks & mitigations

- **PlatformIO may not provide the bench toolchain (the one risk that can change the driver, not
  just a number).** `platform-teensy` has historically lagged Teensyduino, so TD 1.59 / arm-gcc
  11.3.1 may not be cleanly pinnable. Mitigate with the **first** Phase-0 spike (§6) and its written
  three-branch fallback (pin / `platform_packages` override / accept-and-recalibrate); a hard no-go
  reopens the tool choice (§3). Everything else is downstream of this.
- **Object cache can blind the warning ratchet.** A cached TU emits no warnings, so a warm build's
  warning set shrinks and header-introduced warnings hide. Mitigate by running the warning capture on
  a cache-disabled build, decoupled from the cached size build (§7.2).
- **Source-filter mechanic may leak unintended TUs (highest *build-config* risk).** With `src_dir = .`,
  `build_src_filter` must default-exclude and re-add only the wanted TUs or it sweeps in
  `targets/wasm/`, the CMake `build*/` trees, `.pio/`, etc. Mitigate with the Phase-0 spike (§6)
  that asserts each env compiles exactly its `.ino` + the two `core/*.cpp` and nothing else.
- **OCRAM is structurally tight (~26 KiB free after the two fixed framebuffers).** Not a tuning
  knob — the buffers are fixed at `MAX_W*MAX_H`. The known consumers (timeline events; Phantasm's
  small `DMALEDController` TX buffers) are sub-KiB, so the image very likely fits today; the value is
  the **layout invariant** keeping those buffers in OCRAM (`dma_tx_buffers_section: OCRAM`) plus the
  recorded as-built margin (§8). A future OCRAM consumer, or a bump to `MAX_W`/`MAX_H`, is the real
  risk the `ram2` ceiling + layout check must catch (§7.4 #2, §8).
- **Build-option drift between VMicro and PlatformIO defeats the size numbers.** If CI compiles at
  a different optimization level / USB type / `f_cpu` than the menu options used to flash, "fits in
  CI" ≠ "fits on the bench." Mitigate by pinning those options in `platformio.ini` to the captured
  VMicro selections (§4, Phase-0 deliverable).
- **Toolchain download weight / flakiness.** Mitigate with `~/.platformio` caching keyed on the
  pinned platform + FastLED + PIO-Core versions, target-agnostic so one entry serves both envs (§10).
- **arm-gcc surfaces warnings Clang/Emscripten don't.** Expected and *valuable* — that's coverage
  the other builds can't give. Phase the warning gate (§13) so it doesn't block landing the build
  gate.
- **Budgets become stale / bit-rot.** Keep the ceilings in a single reviewed JSON
  (`tools/teensy_budgets.json`), calibrated from real `teensy_size` output, not guesses; the
  warning baseline has its own reviewed `--update-baseline` path (§7.2).
- **VMicro/PlatformIO drift** (a file builds in one, not the other). Mitigate with strict
  include-path + library parity (§4) and the local `just teensy-size` pre-flight (§11). Note: the
  two will never be bit-identical artifacts; the contract is *both compile/link cleanly*, not
  *same bytes*.
- **`.ino` auto-prototype differences** between Arduino's and PlatformIO's preprocessors. Low risk
  (both implement the same Arduino convention); keeping `.ino` canonical rather than converting to
  `.cpp` avoids introducing a third behavior.
- **This gate does not run effects on-device.** It is a compile/link/size gate, not a behavioral
  one. On-device execution and the sim≠device numeric forks remain covered by manual flashing and
  the host device-value tests respectively (divergence ledger). Stated as a non-goal (§2) so the
  gate isn't mistaken for hardware validation.
```

---

## 16. Phase-0 spike log (2026-06)

The go/no-go toolchain spike (§6) was run with a real PlatformIO build. Outcome and the concrete
values now pinned:

**Toolchain — branch 1 (exact bench parity), no override needed.** `platform = teensy@5.0.0`
natively installs, for Teensy 4, the *same* PJRC toolchain Teensyduino/VMicro use:
`framework-arduinoteensy @ 1.159.0` (Teensyduino 1.59.0) and `toolchain-gccarmnoneeabi-teensy @
1.110301.0` (arm-none-eabi-gcc 11.3.1). The `toolchain-gccarmnoneeabi @ ~1.80201.0` (gcc 8.2.1) in
`platform.json` is for the older Teensy 3.x/LC boards, not Teensy 4 — an earlier draft misread it.
Pinned: `platform = teensy@5.0.0`, `platformio == 6.1.19`, `fastled/FastLED @ 3.10.3`. CI thus
*mirrors* the bench compiler. (Pin remains a version tag; a future SHA pin is optional hardening.)

**Two config corrections the build forced (now in `platformio.ini` / `tools/teensy_pre.py`):**
1. **`-std=gnu++17` must be unflagged.** The Teensy core appends `-std=gnu++17` after our flags, so
   `-std=gnu++20` lost and the engine's C++20 `concept` (memory.h) failed. Fixed via
   `build_unflags = ... -std=gnu++17`.
2. **Sketch discovery ignores `build_src_filter`.** PlatformIO finds the sketch only by globbing
   `$PROJECT_SRC_DIR/*.ino` (top level; `pioino.FindInoNodes`). With `src_dir = .` that finds no
   `.ino`, so `setup`/`loop` don't link. The spec's filter-only plan (§5/§6) cannot place the
   sketch. Fix: keep `src_dir = .` (so `core/*.cpp` build as project sources with full LDF include
   paths — FastLED, framework SPI) and override `FindInoNodes` per-env in `tools/teensy_pre.py`;
   `build_src_filter` then adds the converted `targets/<X>/<X>.ino.cpp`. Verified: `core/memory.cpp`,
   `core/reaction_graph.cpp`, and the converted sketch all compile.

**Headline finding (now resolved) — the gate surfaced real device-only bit-rot the WASM/native CI
cannot see.** Each was fixed; the device build is the first thing to compile this code under
`Arduino.h` + arm-gcc:
- **Arduino `TWO_PI` macro collision** — `effects/{Flyby,Liquid2D,Raymarch}.h` declared a local
  `constexpr float TWO_PI`, colliding with `wiring.h`'s `#define TWO_PI`. **Fixed:** renamed → `kTwoPi`.
- **Wrong FastLED pin** — I first pinned 3.10.3; the bench/Teensyduino 1.59 bundles **FastLED 3.4.0**,
  and 3.7+ moved `CHSVPalette16`/`HUE_RED`/`CEveryNMillis` and `CRGB`/`CHSV` into `namespace fl`.
  **Fixed:** pinned `FastLED@3.4.0` (§4.1 parity) — cleared a whole error class at once.
- **`effects_legacy.h` `Pixel16`→`CRGB`** (lines 25, 42, 100–112, 517, 743): `operator CRGB()` was made
  `explicit`, so the legacy effects no longer convert implicitly. **Fixed:** explicit `static_cast<CRGB>`.
- **Ambiguous `pov`** — the sketch var collided with the hardware `namespace pov` (pov_segment_map.h).
  **Fixed:** sketch var → `g_pov` in both `.ino`s.
- **`color.h` include collision** — FastLED ships its own `src/color.h`; the device `-I` order made
  `hardware/{dma_led,hd107s_frame}.h`'s bare `"color.h"` resolve to FastLED's (no `Pixel16`/LUTs),
  while host builds got the engine's. **Fixed:** those two now `#include "core/color.h"`.

**Result — Holosphere builds GREEN and the gate PASSES end-to-end on the real ELF.** Real calibrated
budgets (FLASH 154,620 / RAM1 459,072 free 65,216 / RAM2 518,272 free 6,016) and the **real** mangled
arena symbol `_ZL18global_arena_block` (the source spelling would never have matched — §7.4 vindicated).
The layout invariants are confirmed on real ELFs: arena→DTCM (343,040 B), framebuffers→OCRAM, and (on
a Phantasm ELF) reaction_graph→FLASH (`0x60003dc0`, 92,160 B). Parser quirks the real output exposed,
now handled: `teensy_size` prints to **stderr**, and `readelf` prints large sizes in **hex**. The
warning baseline is captured: **1** first-party warning (55 `HS_CHECK` call sites dedupe to one —
the set-based ratchet working as designed; and the object cache was observed blinding it, §7.2).

**Phantasm OVERFLOWS RAM1 by ~243 KB** (all effects via `HS_EFFECT_LIST` → ~381 KB ITCM + ~380 KB
DTCM). Its ELF was captured for symbol analysis by temporarily neutralizing `teensy_size`'s overflow
exit (reverted). RAM2 is ~518 KB / 524 KB for **both** targets — the structural OCRAM tightness the
spec predicted (§8). Shrinking Phantasm (-Os / FLASHMEM / fewer effects / smaller arena) is owner-scope.
