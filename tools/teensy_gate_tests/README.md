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
| `broken_negative_free_teensy_size.txt` | Negative RAM1 stack headroom remains parseable and fails the configured floor. |
| `real/verbose_build_log.txt` | Verbatim `pio run -v` compiler invocations — first-party then third-party, Windows then Linux-CI — for the warning ratchet's capture-evidence guard. The Windows home directory is rewritten to `C:\Users\dev`; nothing else is edited. |
| `real/cold_env_section.txt` | Verbatim `holosphere` section of a `pio run -v` with `.pio/build_cache` deleted: PlatformIO's banner (the `build_src_filter` the expectation is derived from) plus all three first-party compiles. |
| `real/warm_env_section.txt` | The same section from the next run, reusing that cache: two `Retrieved … from cache` lines in place of the core compiles. The ratchet must fail on it. |
| `real/phantasm_teensy_size.txt` | Verbatim `teensy_size` for `[env:phantasm]` — the calibrated region totals the end-to-end `evaluate()` test runs against. |
| `real/phantasm_readelf_syms.txt` | `readelf -sW` rows for the layout symbols `tools/teensy_budgets.json` names, excerpted verbatim from the full table (the whole dump is ~600 KB). |
| `real/phantasm_readelf_secs.txt` | Verbatim `readelf -SW`; also pins `.ARM.exidx` to flash, where `tools/phantasm.ld` routes it. |
| `real/phantasm_size_a.txt` | Verbatim `size -A` for the fallback bucketing path. |

## ⚠ These are synthetic-but-realistic, not hardware-captured (yet)

The addresses, sizes and mangled names are hand-built to match the Teensy 4
memory map and the real source symbols, but they were **not** captured from an
actual `-O3` firmware link (no Teensy toolchain in this environment). A Phase-0
deliverable is to **replace the `good_*` fixtures with truly-captured output**
from a real build (`arm-none-eabi-size -A`, `readelf -s -S`) so the parser tests
exercise the exact toolchain formatting, then re-derive the broken variants from
that capture (spec §9.1, §13).
