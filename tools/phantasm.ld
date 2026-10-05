/* ============================================================================
 * Custom Teensy 4.0 linker script for the Phantasm environments.
 * FORKED FROM: framework-arduinoteensy 1.162.0
 *
 * FLASH routing covers cold libraries and effect factories, const rodata,
 * spawn-time init/setup_presets, cold-attribute and hot-leaf sections, cold
 * constructors/destructors and mesh builders, and framework one-time setup.
 * The hot flywheel libgcc divmod also uses cached flash, retaining ITCM capacity
 * at the cost of cache-dependent latency.
 * The Teensy core USB/serial stack and one-time setup objects also use flash.
 * .fini is placed explicitly. Function/data sections allow individual matches
 * before the ITCM catch-all; _etext and _itcm_block_count exclude .ARM.exidx.
 *
 * tests/check_linker_patches.cmake pins the fork version and selected library,
 * factory and rodata routes. A framework update requires comparing the full
 * script against the new stock imxrt1062.ld and preserving all routing rules.
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
		/* Startup and topology/control sections. */
		*(.text.startup*)
		*(.text._ZN6Solids12SolidBuilder*)
		*(.text._ZN7MeshOps5cloneI8PolyMeshE*)
		*(.text._ZN8TimelineC*Ev*)
		*(.text._ZN8TimelineD*Ev*)
		*(.text._ZN8Timeline14destroy_eventsEv*)
		/* Libraries (including hot libgcc divmod) and factories use flash. */
		*libm.a:*(.text .text.*)
		*libc_nano.a:*(.text .text.*)
		*libstdc++_nano.a:*(.text .text.*)
		*libg_nano.a:*(.text .text.*)
		*libgcc.a:*(.text .text.*)
		*(.text.*construct_effect*)
		/* Effect spawn-time setup: init()/setup_presets(). */
		*(.text._ZN*4initEv)
		*(.text._ZN*13setup_presetsEv)
		/* HS_COLD_MEMBER (attribute cold) sections. */
		*(.text.unlikely*)
		/* Selectively optimized hot leaves intentionally executed from cache. */
		*(.text.hot*)
		/* Sections without source-level cold annotations: generated ctors
		 * (Animation graph nodes, PolyMesh), sort instantiations from
		 * classify_faces_impl, setup_presets/compile_hankin inner lambdas
		 * (preset/geometry builders only -- stored draw callbacks must not
		 * be defined inside them), and framework one-time setup. */
		*(.text._ZN9Animation*C1E*)
		*(.text._ZN9Animation*C2E*)
		*(.text._ZN8PolyMeshC1E*)
		*(.text._ZN8PolyMeshC2E*)
		/* Destructor sections for Animation graph nodes, <288,144> effects,
		 * Persist/ReactionDiffusionBase, and PolyMesh use flash. */
		*(.text._ZN9Animation*D0Ev*)
		*(.text._ZN9Animation*D1Ev*)
		*(.text._ZN9Animation*D2Ev*)
		*(.text._ZN8PolyMeshD0Ev*)
		*(.text._ZN8PolyMeshD1Ev*)
		*(.text._ZN8PolyMeshD2Ev*)
		*(.text._ZN7PersistI*D0Ev*)
		*(.text._ZN7PersistI*D1Ev*)
		*(.text._ZN7PersistI*D2Ev*)
		*(.text._ZN*ILi288ELi144EED0Ev*)
		*(.text._ZN*ILi288ELi144EED1Ev*)
		*(.text._ZN*ILi288ELi144EED2Ev*)
		*(.text.*classify_faces_impl*)
		*(.text.*KDTree5build*)
		*(.text._ZZN*13setup_presetsEv*)
		*(.text.*compile_hankin*)
		*(.text._ZN14CompiledHankin*)
		*(.text._ZN8Timeline7add_getI*)
		*(.text._ZZN8Timeline7add_getI*)
		/* OpLeg scheduling and mesh build sections: timeline construction,
		 * mesh move assignment. MeshState move
		 * assignment also runs during per-frame compilation. */
		*(.text._ZN8Timeline3addIN9Animation*)
		*(.text._ZNO9Animation13AnimationBaseINS_5OpLeg*)
		*(.text._ZN9MeshStateaSEOS_*)
		*(.text._ZN11ArenaVectorIiEaSEOS0_*)
		/* Half-edge topology construction: the edge-record sort family and
		 * mesh normalize also serve per-frame truncate/expand sweeps. */
		*(.text.*HalfEdgePairRecord*)
		*(.text._ZN7MeshOpsL9normalizeI*)
		/* inplace_function invokers for stored Canvas& and ParticleSystem&
		 * callbacks can run for multiple events or emitters per frame.
		 * The per-pixel path uses FunctionRef, a distinct type. */
		*(.text._ZZN6teensy23inplace_function_detail6vtableIvJR6Canvas*)
		*(.text._ZZN6teensy23inplace_function_detail6vtableIvJRN9Animation14ParticleSystem*)
		/* Teensy core USB/serial stack, one-time pin/PWM/RTC/temperature
		 * setup, and the SPI library (the LED DMA driver drives LPSPI4
		 * registers directly after SPI.begin()). */
		*libFrameworkArduino.a:usb.c.o(.text .text.*)
		*libFrameworkArduino.a:usb_serial.c.o(.text .text.*)
		*libFrameworkArduino.a:usb_desc.c.o(.text .text.*)
		*libFrameworkArduino.a:usb_inst.cpp.o(.text .text.*)
		*libFrameworkArduino.a:pwm.c.o(.text .text.*)
		*libFrameworkArduino.a:digital.c.o(.text .text.*)
		*libFrameworkArduino.a:tempmon.c.o(.text .text.*)
		*libFrameworkArduino.a:rtc.c.o(.text .text.*)
		*libFrameworkArduino.a:Print.cpp.o(.text .text.*)
		*libFrameworkArduino.a:nonstd.c.o(.text .text.*)
		*libFrameworkArduino.a:startup.c.o(.text.unused_interrupt_vector)
		*libFrameworkArduino.a:interrupt.c.o(.text.attachInterrupt)
		*libSPI.a:*(.text .text.*)
		/* Spawn, preset-change and palette-change templates: solid builders,
		 * vertex reconciliation and its sort, relax-bake hashing, palette
		 * recipe validation, timeline storage reset, preset switching,
		 * MindSplatter warp scheduling, mesh edge extraction, and the
		 * composed-effect projection-walk constructor with its two-vector
		 * make_rotation. */
		*(.text._ZN6Solids11to_polymesh*)
		*(.text._ZN6Solids9get_entry*)
		*(.text._ZN7MeshOps17relax_source_hash*)
		*(.text._ZZN7MeshOpsL18reconcile_vertices*)
		*(.text._ZSt*MeshOpsL18reconcile_vertices*)
		*(.text._ZN17GenerativePalette15validate_finite*)
		*(.text._ZN17GenerativePalette10axis_state*)
		*(.text._ZN17GenerativePalette6finite*)
		*(.text._ZNK17GenerativePalette16morph_compatible*)
		*(.text._ZN14PaletteRecipes7harmony*)
		*(.text._ZN8Timeline13reset_storage*)
		*(.text._ZN10PresetHost13change_preset*)
		*(.text._ZN12MindSplatterILi288ELi144EE13schedule_warp*)
		*(.text._ZN12MindSplatterILi288ELi144EE12perform_warp*)
		*(.text._ZN4Plot4Mesh13extract_edges*)
		*(.text._ZN4Plot4Mesh20for_each_unique_edge*)
		*(.text._ZN11ArenaVectorIN4math6Snorm3EE*)
		*(.text._ZN4math13make_rotationERKNS_6VectorES2_)
		*(.text._ZN8Pullback19ProjectionWalkStateILb1EEC*)
		*(.text.set_arm_clock)
		*(.text._ZN8SPIClass5beginEv)
		*(.text._ZN13IntervalTimer11beginCycles*)
		*(.text._ZN10DMAChannel5begin*)
		. = ALIGN(4);
		KEEP(*(.init))
		/* Placed explicitly: .fini matches no other output section, so the
		 * linker would orphan it into ITCM past .text.itcm, outside the boot
		 * copy's [_stext, _etext) range. */
		KEEP(*(.fini))
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
		/* Route all const rodata to flash instead of DTCM (RAM1). */
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
	_etext = ADDR(.text.itcm) + SIZEOF(.text.itcm);
	_stextload = LOADADDR(.text.itcm);

	_sdata = ADDR(.data);
	_edata = ADDR(.data) + SIZEOF(.data);
	_sdataload = LOADADDR(.data);

	_sbss = ADDR(.bss);
	_ebss = ADDR(.bss) + SIZEOF(.bss);

	_heap_start = ADDR(.bss.dma) + SIZEOF(.bss.dma);
	_heap_end = ORIGIN(RAM) + LENGTH(RAM);

	_itcm_block_count = (SIZEOF(.text.itcm) + 0x7FFF) >> 15;
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
