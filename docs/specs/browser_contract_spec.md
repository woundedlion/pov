# Browser engine contracts

*Status: IMPLEMENTED. Holosphere owns the browser ABI and installed contracts.*

Embind emits the bound method signatures and numeric enum members during the
WASM link. JavaScript payload types are registered in
`targets/wasm/ts_types.h`; authored object schemas and loader hooks live in
`targets/wasm/engine_types.d.ts`. The install derives browser aliases from the
emitted declarations. Class-returning methods admit null; constructors retain
their emitted non-null return types.

`scripts/export_engine_contract.mjs` loads the built engine and exports its
simulation cadence, enum numbers, and bound instance method names. Browser
frame timing and test fakes consume that output. Display labels and test-specific
optional calls remain authored by their consumers.

`hardware/pov_segment_map.h` owns row-band ordering and partition arithmetic.
The device uses integer coordinates; the standalone browser geometry build uses
double coordinates to retain JavaScript's admitted integer range. Its generated
JavaScript initializes a fixed 64 KB heap and a 1 KB stack, without loading
another WASM binary. Segment calculations remain synchronous after module
initialization and run when configuring a partition; pixel rendering introduces
no geometry calls per frame.

The install and verified engine bundle carry declarations, object schemas,
runtime contracts, and the geometry module together with the engine. These
assets are ignored install outputs in Daydream. The bundle checksum manifest
covers them. Worker warming includes the geometry module; the site manifest
publishes both generated runtime modules.
