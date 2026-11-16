"""Teensy 4 address regions and FlexRAM geometry."""

# ---------------------------------------------------------------------------
# Teensy 4 (i.MX RT1062) memory map, half-open [lo, hi). Classify a symbol by
# its VMA, not by section name or nm letter.
#   ITCM  fast code              (part of RAM1)
#   DTCM  fast data + stack      (part of RAM1)
#   OCRAM DMAMEM + heap          (RAM2)
#   FLASH code + const/PROGMEM
# RAM1 in the budget == ITCM ∪ DTCM.
# ---------------------------------------------------------------------------
MEMORY_MAP: tuple[tuple[str, int, int], ...] = (
    ("ITCM", 0x00000000, 0x00080000),   # 512 KiB window; ITCM is carved from RAM1
    ("DTCM", 0x20000000, 0x20080000),   # 512 KiB window
    ("OCRAM", 0x20200000, 0x20280000),  # 512 KiB RAM2
    ("FLASH", 0x60000000, 0x60200000),  # 2 MiB T4.0 flash
)

#: Program-image capacity of FLASH: the linker scripts' `LENGTH = 1984K`.
FLASH_PROGRAM_BYTES = 0x1F0000

#: FlexRAM allocation granule. The 512 KiB of RAM1 is split ITCM/DTCM in whole
#: banks of this size, so ITCM code always occupies a bank-rounded footprint.
FLEXRAM_BANK_BYTES = 0x8000
