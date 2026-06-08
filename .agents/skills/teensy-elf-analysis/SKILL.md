---
name: Teensy ELF Analysis
description: Analyze Teensy 4.x ARM ELF binaries for code size, ITCM usage, and symbol bloat
---

# Teensy ELF Analysis

## ELF Location

The Teensy ELF is built by Visual Micro to:
```
C:\Users\gabe\AppData\Local\Temp\VMBuilds\Holosphere\teensy40\Release\Holosphere.ino.elf
```

## ARM Toolchain

The ARM toolchain is installed via Arduino/Teensyduino at:
```
C:\Users\gabe\AppData\Local\Arduino15\packages\teensy\tools\teensy-compile\11.3.1\arm\bin\
```

Key tools (prefix all with the path above):
- `arm-none-eabi-nm.exe` — symbol table (names, sizes, sections)
- `arm-none-eabi-objdump.exe` — disassembly and section analysis
- `arm-none-eabi-size.exe` — section size summary
- `arm-none-eabi-strings.exe` — embedded string literals
- `arm-none-eabi-readelf.exe` — ELF header and section details

## Common Analysis Commands

All commands use `cmd /c` prefix for PowerShell compatibility.

### Section sizes overview
```
cmd /c "arm-none-eabi-size.exe -A Holosphere.ino.elf"
```
Key sections:
- `.text.itcm` — ITCM (instruction tightly-coupled memory, 128KB max on Teensy 4.0)
- `.text.progmem` — PROGMEM (flash, no size pressure)
- `.text.code` — general code
- `.data` — initialized data in RAM2
- `.bss` — zero-initialized data in RAM2
- `.bss.dma` — DMA-accessible memory

### Largest symbols (sorted by size)
```
cmd /c "arm-none-eabi-nm.exe -S --size-sort Holosphere.ino.elf"
```
Output format: `address size type name`
- `T` = .text (code), `t` = static/local code
- `W` = weak symbol
- `B`/`b` = .bss, `D`/`d` = .data

### Find specific symbols
```
cmd /c "arm-none-eabi-nm.exe -S --size-sort Holosphere.ino.elf 2>&1 | findstr /i keyword"
```

### Find string literals in the binary
```
cmd /c "arm-none-eabi-strings.exe Holosphere.ino.elf 2>&1 | findstr /i keyword"
```

### Disassemble a specific section
```
cmd /c "arm-none-eabi-objdump.exe -d -j .text.itcm Holosphere.ino.elf"
```

## Teensy 4.0 Memory Map

| Region | Size | Section | Purpose |
|--------|------|---------|---------|
| ITCM | 128 KB max (shared with DTCM) | `.text.itcm` | Fast instruction memory, all non-FLASHMEM code |
| DTCM | 512 KB - ITCM | `.data`, `.bss` | Fast data (stack, globals) |
| RAM2 | 512 KB | `.bss.dma` | DMA-accessible, `DMAMEM` variables |
| Flash | 2 MB | `.text.progmem` | `PROGMEM` data, `FLASHMEM` code |

ITCM + DTCM share 512 KB of RAM1. ITCM is allocated in 32KB granules. Reducing ITCM below a granule boundary frees 32KB for DTCM.

## Common Bloat Sources

1. **printf/fprintf** (~1.4KB) — pulled in by `assert()` via newlib's `__assert_func`. Fix: define `NDEBUG` or stub `__assert_func`.
2. **C++ demangler** (~15KB) — pulled in by `std::function` → `__cxa_throw` → `__verbose_terminate_handler`. Fix: stub `__verbose_terminate_handler` (already done in `memory.cpp`).
3. **Template instantiation bloat** — each `Effect<W,H>` instantiation duplicates the entire pipeline. Use `FLASHMEM` for cold-path code.
4. **Math library** — `atan2f`, `acosf`, `powf` in hot loops. Replace with fast approximations where visual difference is imperceptible.
5. **Double-precision literals** — `0.5` instead of `0.5f` promotes operations to double, pulling in soft-float routines.

## Annotation Attributes

- `FLASHMEM` — places function in flash (.text.progmem), freeing ITCM
- `PROGMEM` — places data in flash
- `DMAMEM` — places data in DMA-accessible RAM2
- `FASTRUN` — explicitly places function in ITCM (default for most code)
