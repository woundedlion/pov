# Teensy gate fixtures & self-tests

These prove the size/layout gate (`tools/teensy_gate.py`) **fails when it should**
— a check that can never fail is worse than none (docs/teensy_ci_gate_spec.md §9.1).
They are pure host Python (`unittest`, no ARM toolchain, no PlatformIO).

```
python -m unittest discover -s tools/teensy_gate_tests
```

## Fixtures

| File | Role |
|---|---|
| `good_teensy_size.txt` | A passing `teensy_size` region report. |
| `good_readelf_syms.txt` | `readelf -s` with the real **mangled** symbol names (`_ZN6Effect8buffer_aE`, `_ZN13ReactionGraph9neighborsE`) and the internal-linkage `global_arena_block` (a LOCAL symbol). |
| `good_readelf_secs.txt` | `readelf -S` section headers (ndx → name/addr). |
| `good_size_a.txt` | `size -A` for the address-bucketing test. |
| `broken_framebuffer_dtcm_syms.txt` | `DMAMEM` dropped → `Effect::buffer_a/b` land in DTCM → framebuffers→OCRAM must fail. |
| `broken_reaction_graph_ram_syms.txt` | `const` dropped → `ReactionGraph::neighbors` in DTCM → table→FLASH must fail. |
| `broken_arena_8mb_syms.txt` | `global_arena_block` at the 8 MB `HS_TEST_BUILD` size → arena magnitude must fail. |
| `broken_missing_symbol_syms.txt` | a framebuffer renamed → "symbol not found" hard-fails (never a silent skip). |
| `broken_over_cap_teensy_size.txt` | FLASH + RAM2 over cap, DTCM headroom under floor → every region check fails. |

## ⚠ These are synthetic-but-realistic, not hardware-captured (yet)

The addresses, sizes and mangled names are hand-built to match the Teensy 4
memory map and the real source symbols, but they were **not** captured from an
actual `-O3` firmware link (no Teensy toolchain in this environment). A Phase-0
deliverable is to **replace the `good_*` fixtures with truly-captured output**
from a real build (`arm-none-eabi-size -A`, `readelf -s -S`) so the parser tests
exercise the exact toolchain formatting, then re-derive the broken variants from
that capture (spec §9.1, §13).
