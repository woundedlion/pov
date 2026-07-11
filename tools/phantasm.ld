/* ============================================================================
 * Custom Teensy 4.0 linker script for the Phantasm SIZE build — a verbatim copy
 * of framework-arduinoteensy/cores/teensy4/imxrt1062.ld with TWO routing changes
 * that move RAM1 (FlexRAM) content to FLASH. Both rely on -ffunction-sections /
 * -fdata-sections (on by default), so each function/const is its own input
 * section and can be matched individually:
 *
 *   #5  cold code -> FLASH: library .text (libm/libc_nano/libstdc++_nano/libgcc)
 *       and the per-effect construct_effect<> factories run once at spawn / are
 *       cold, so they leave ITCM. Routed in .text.code BEFORE .text.itcm's
 *       catch-all *(.text*), so the linker assigns them to flash first.
 *
 *   #3  const rodata -> FLASH: the stock script sends *(.rodata*) to .data (DTCM,
 *       RAM1). Teensy 4 flash is memory-mapped + cached, so const tables/strings/
 *       vtables read fine in place (and skip the boot-time copy to RAM). Routed in
 *       .text.progmem; the .data rodata line is removed so nothing double-claims.
 *
 * Everything else (FlexRAM bank math, heap, stack, csf) is byte-identical to the
 * stock script. Selected per-env via board_build.ldscript (phantasm only).
 * ============================================================================ */
MEMORY
{
	ITCM (rwx):  ORIGIN = 0x00000000, LENGTH = 512K
	DTCM (rwx):  ORIGIN = 0x20000000, LENGTH = 512K
	RAM (rwx):   ORIGIN = 0x20200000, LENGTH = 512K
	FLASH (rwx): ORIGIN = 0x60000000, LENGTH = 1984K
}

ENTRY(ImageVectorTable)

SECTIONS
{
	.text.headers : {
		KEEP(*(.flashconfig))
		FILL(0xFF)
		. = ORIGIN(FLASH) + 0x1000;
		KEEP(*(.ivt))
		KEEP(*(.bootdata))
		. = ALIGN(1024);
	} > FLASH

	.text.code : {
		KEEP(*(.startup))
		*(.flashmem*)
		/* #5: cold library code + per-effect factories run from flash, not ITCM. */
		*libm.a:*(.text .text.*)
		*libc_nano.a:*(.text .text.*)
		*libstdc++_nano.a:*(.text .text.*)
		*libg_nano.a:*(.text .text.*)
		*libgcc.a:*(.text .text.*)
		*(.text.*construct_effect*)
		/* Effect spawn-time setup (init()/setup_presets()) runs once per effect
		 * (~every 120 s), never per frame — route to flash too. */
		*(.text._ZN*4initEv)
		*(.text._ZN*13setup_presetsEv)
		/* HS_COLD_MEMBER (attribute cold) sections. */
		*(.text.unlikely*)
		/* Annotation-unreachable cold code: compiler-generated ctors
		 * (Animation graph nodes, PolyMesh), sort instantiations from
		 * classify_faces_impl, setup_presets/compile_hankin inner lambdas
		 * (preset/geometry builders only -- stored draw callbacks must not
		 * be defined inside them), and framework one-time setup. */
		*(.text._ZN9Animation*C1E*)
		*(.text._ZN9Animation*C2E*)
		*(.text._ZN8PolyMeshC1E*)
		*(.text._ZN8PolyMeshC2E*)
		*(.text.*classify_faces_impl*)
		*(.text._ZZN*13setup_presetsEv*)
		*(.text.*compile_hankin*)
		*(.text.set_arm_clock)
		*(.text._ZN8SPIClass5beginEv)
		*(.text._ZN13IntervalTimer11beginCycles*)
		*(.text._ZN10DMAChannel5begin*)
		. = ALIGN(4);
		KEEP(*(.init))
		__preinit_array_start = .;
		KEEP (*(.preinit_array))
		__preinit_array_end = .;
		__init_array_start = .;
		KEEP (*(.init_array))
		__init_array_end = .;
		. = ALIGN(4);
	} > FLASH

	.text.progmem : {
		*(.progmem*)
		/* #3: route all const rodata to flash instead of DTCM (RAM1). */
		*(SORT_BY_ALIGNMENT(SORT_BY_NAME(.rodata*)))
		. = ALIGN(4);
	} > FLASH

	.text.itcm : {
		. = . + 32; /* MPU to trap NULL pointer deref */
		*(.fastrun)
		*(.text*)
		. = ALIGN(16);
	} > ITCM  AT> FLASH

	/* Unwind tables to FLASH: the only .ARM.exidx left comes from the prebuilt
	 * libs whose .text we just moved to flash (engine/framework/FastLED build
	 * -fno-unwind-tables). Keeping exidx in ITCM while its .text is in flash
	 * overflows the R_ARM_PREL31 (±1 GB) relocation, so it must follow .text. */
	.ARM.exidx : {
		__exidx_start = .;
		*(.ARM.exidx* .ARM.extab.text* .gnu.linkonce.armexidx.*)
		__exidx_end = .;
	} > FLASH

	.data : {
    		*(.endpoint_queue)
		*(SORT_BY_ALIGNMENT(SORT_BY_NAME(.data*)))
   		 KEEP(*(.vectorsram))
	} > DTCM  AT> FLASH

	.bss ALIGN(4) : {
		*(SORT_BY_ALIGNMENT(SORT_BY_NAME(.bss*)))
		*(SORT_BY_ALIGNMENT(SORT_BY_NAME(COMMON)))
		. = ALIGN(32);
		. = . + 32; /* MPU to trap stack overflow */
	} > DTCM

	.bss.dma (NOLOAD) : {
		*(.hab_log)
		*(.dmabuffers)
		. = ALIGN(32);
	} > RAM

	.text.csf : {
		FILL(0xFF)
		. = ALIGN(1024);
		KEEP(*(.csf))
		__text_csf_end = .;
	} > FLASH

	_stext = ADDR(.text.itcm);
	_etext = ADDR(.text.itcm) + SIZEOF(.text.itcm) + SIZEOF(.ARM.exidx);
	_stextload = LOADADDR(.text.itcm);

	_sdata = ADDR(.data);
	_edata = ADDR(.data) + SIZEOF(.data);
	_sdataload = LOADADDR(.data);

	_sbss = ADDR(.bss);
	_ebss = ADDR(.bss) + SIZEOF(.bss);

	_heap_start = ADDR(.bss.dma) + SIZEOF(.bss.dma);
	_heap_end = ORIGIN(RAM) + LENGTH(RAM);

	_itcm_block_count = (SIZEOF(.text.itcm) + SIZEOF(.ARM.exidx) + 0x7FFF) >> 15;
	_flexram_bank_config = 0xAAAAAAAA | ((1 << (_itcm_block_count * 2)) - 1);
	_estack = ORIGIN(DTCM) + ((16 - _itcm_block_count) << 15);

	_flashimagelen = __text_csf_end - ORIGIN(FLASH);
	_teensy_model_identifier = 0x24;

	.debug_info     0 : { *(.debug_info) }
	.debug_abbrev   0 : { *(.debug_abbrev) }
	.debug_line     0 : { *(.debug_line) }
	.debug_frame    0 : { *(.debug_frame) }
	.debug_str      0 : { *(.debug_str) }
	.debug_loc      0 : { *(.debug_loc) }

}
