# Chain snapshots

**Status: IMPLEMENTED.** This contract covers the native chain host, typed
operator state codecs, WASM authoring capability, and daydream's worker,
persistence and geometry reconstruction. Shader documents describe authoring
intent; snapshots describe the next executable frame.

`ShaderChain` is the sole simulator shader host. Firmware excludes the chain
interpreter through `HS_ENABLE_CHAIN_INTERPRETER`. Document imports, links and
snapshots accept their current formats only; retired names and formats are rejected.

## Wire format

`ShaderChainBindings.getSnapshot()` returns an owned plain object:

```text
{
  schemaVersion: 2,
  chain: [{instance, operator}],
  parameters: [{name, value}],
  runtime: [{instance, kind, state}],
  paletteBank: {chroma, hues: [u32, u32, u32], cycles: [clock, clock, clock]},
  animationsPaused: boolean
}
clock = {frame: u32, nextSequence: u32, fadeActive: boolean, displayDirty: boolean}
```

Names are `instance.field-id`; values are finite binary32 numbers. Topology
values are integral catalog indices. Instance and operator IDs obey the catalog
budgets. Program and parameter order are stable. Parameter names are unique;
omitted names receive catalog defaults. A supplied runtime array has exactly one
entry for every stateful instance and none for stateless instances. Omitting
runtime initializes fresh state. Omitting the palette bank initializes the three
generated harmonies with the restored colorizer's chroma.

| Kind | State properties | Operators |
|---|---|---|
| `spatial-walk-v2` | `noiseSeed` i32, `walkTime` u32, `position`, `direction`, `wander`, `angularVelocity`, `spinPhase` | Rotate and all Project operators |
| `source-clock-v1` | `primary`, `secondary`, `angle` | Grid, TwinWave, Rings, Spiral, Fractal, Tessellation |
| `noise-clock-v1` | `phase`, `noiseSeed` i32 | Curl/Direct displacement, VectorNoise/CurlFlow warp, projected/spherical noise samples |
| `phase-clock-v1` | `phase` | WaveShear, Vortex, MirrorTile, PolarChart |
| `ripple-clock-v1` | `phase` | Ripple displacement |
| `affine-clock-v1` | `phase`, `rotation` | Affine warp |
| `color-clock-v1` | `oscillationPhase`, `hueNoisePhase`, `hueNoiseSeed` i32 | Generated-palette colorizers |
| `spherical-rings-v2` | `walk` spatial state, `phase` | Spherical Rings |

Vectors are `[x,y,z]`; quaternions are `[real,x,y,z]`. Spatial vectors and
quaternions are finite unit values; position and direction are perpendicular.
Loop phases lie in `[0,1)`, source and spin angles in `(-2*pi,2*pi)`, and ripple
phase is finite and nonnegative. Lattice, fixed lenses, Transfer and Coverage
operators are stateless. Fresh walks start at `[0,1,0]` with direction `[0,0,-1]`,
identity orientations, zero clocks and angular velocity.

Palette hues and sequence clocks preserve both endpoints of each generated fade,
including inactive harmonies. Recipes and baked tables are reconstructed;
pointer values, allocator addresses and derived cache generations are excluded.
All three hue words equal `(nextSequence-1)*159` modulo uint32. Fade frame and
sequence bounds must pass the palette cycler's validation.

## Atomic restoration

`restoreSnapshot()` returns `Module.ChainSnapshotRestoreResult`: `APPLIED`,
`NOT_SHADER_CHAIN`, `UNSUPPORTED_VERSION`, `INVALID_LENGTH`, `INVALID_VALUE` or
`INVALID_CHAIN`. Compare enum values explicitly. Every refusal preserves the
program, parameter definitions, schema generation, clocks, noise, palettes and
pause state.

The host compiles into its inactive arena with state migration disabled, then
validates parameter ranges, topology indices, operator predicates, edge-distance
dependencies and complete typed runtime. It flips the active arena only after
all checks pass. Palette admission precedes compilation; restoration after the
flip cannot fail. Cache inputs are invalidated and rebuilt from restored state.

The WASM capability clones caller data before decoding and verifies its owning
effect incarnation and schema generation again before commit. Effect replacement,
resolution change, geometry reconstruction and deletion invalidate old handles.
Geometry reconstruction captures this same full snapshot, including hidden
palette clocks and noise seeds, before rebuilding the effect. Worker messages
carry owned snapshots and never raw runtime bytes.

## Current operators and verification

The catalog publishes one current identity per operator family. Superseded IDs
are rejected. `sample.grid.v3` and `sample.twin-wave.v3` accept drift from zero
through two. `warp.affine.v3` accepts periods from 1/64 through 100. The current
Peirce, Bonne and Airocean operators expose coordinate scale from 1/4 through
four; Peirce also exposes diamond, square, horizontal and vertical layouts and
scroll from minus one through one.

`tests/data/chain_capture_fixtures.jsonl` fixes the canonical capture programs and
before-frame metadata. `tools/gen_chain_capture_fixtures.py` emits typed snapshots
for the native and WASM producers. Both producers evaluate the chain and shared
exact projection/color kernels for approximation oracles. Native tests also compare compiled composed programs to document-built chains,
reject invalid snapshots atomically and replay future frames after a palette
handoff. The WASM smoke test checks capability ownership, malformed input,
roundtrip frame identity and reentrant replacement.
