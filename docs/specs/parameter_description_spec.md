# Parameter descriptions

**Status: IMPLEMENTED typed registration and shared field descriptions.**
Runtime descriptors retain their existing ABI.

## Typed registration

`core/control/param_spec.h` defines `ParamSpec<T>`. `ParamHost` accepts a name,
a live `T*` target, and this construction-time description through its protected
`register_param` overload. Supported targets are `float`, `bool`, and integer
or enum types of at most 32 bits. The description owns no runtime storage.

Every description specifies bounds, animation, readonly and preset-export
policies, optional display/export label arrays, and an initial-value policy.
Names, labels, export literals, target storage and external descriptor storage
must outlive their registered use.

Float bounds use floats. Integer and enum bounds use signed 64-bit construction
time values, so validation sees the authored bounds before conversion to the
float-based runtime interface. Registration checks storage representability and
exact float representation before assembling a descriptor. Boolean targets
have the range `[0,1]` and no option labels.

`ParamSpec<T>::enumerated` constructs the contiguous option range
`[0, option_count - 1]`. An explicitly constructed dropdown description must
provide those same bounds. Typed enums require option labels; integer and float
targets can optionally use them. Display labels and optional export literals
are indexed in the same order and cannot contain null entries. Export literals
require a display-label array. Array lengths are the caller's responsibility.

## Initial values and writes

Registration does not mutate the target. `REQUIRE_IN_RANGE` is the default;
it rejects a starting value outside the published range, including non-finite
float values. Published float bounds must be finite.

`PRESERVE_REQUESTED_FLOAT` is available only with the parameter GUI bridge and
only for ordinary float targets. It admits a finite requested value outside the
published range without changing that value. It does not admit out-of-range
enums or integer targets, and it does not weaken later write normalization.

Readonly, animation and preset policies are independent for every target type.
Readonly writes return `READONLY` without modifying the value or engaging the
animation pause. An accepted animated write engages the pause. Integer and
option writes round and clamp; ordinary float writes clamp; bool writes retain
the existing threshold behavior. Preset exclusion affects exports independently
of edit and animation policy.

## Compatibility and ownership

The named registration helpers delegate to the typed overload. All typed
validation and runtime `ParamDef` assembly happen in that overload. The runtime
descriptor layout, target representations and read/write paths remain unchanged.

One successful registration appends one descriptor in call order and advances
the parameter-schema generation once. Parameter writes do not change that
generation. The existing reset, external-storage and metadata-mutation contracts
continue to apply.

Behavioral tests live in `tests/test_canvas.h`; invalid ranges, labels and
initial-value policies are covered in `tests/test_death.h`.

## Field descriptions

`core/control/fields.h` defines `Control::Field<Owner, Value>`. Each description
names a typed member pointer, a machine identifier, an optional slider name,
and its `ParamSpec<Value>`. Heterogeneous tuples can describe floats, integers,
enums and booleans together. `Control::FieldGroup` describes a nested aggregate
without changing either aggregate's layout. The control layer depends on no
renderer types.

The same descriptions drive ordered registration, constexpr range validation
and interpolation. A null slider name excludes registration while retaining
validation and interpolation. Validation rejects non-finite floats, unordered
or non-finite bounds, and integral bounds outside the target representation.
An explicit validation bound may differ from the registered slider bound;
ShapeShifter validates authored Count values against its authored capacity,
then clamps adopted Count to the canvas drawing limit.

Continuous domains support endpoint-exact linear interpolation, geometric
positive interpolation, shortest periodic arcs in radians or turns, and snap.
`RAW_LINEAR` retains the arithmetic of existing handwritten transitions.
Integral, enum and bool fields switch at the midpoint by default, or hold the
start value until completion under `SNAP`. Descriptions can independently
exclude validation and interpolation. Interpolation writes only described,
included fields into the destination, preserving telemetry and untabled state.

Choreographed effects provide `parameter_fields()` or declare `Params::FIELDS`. Their base supplies
registration, default range validation, and a default blend hook. Effects retain
custom cross-field admission, adoption and blend hooks: HyperLattice checks
supported pattern/view combinations and refreshes configuration metadata;
Comets validates its path function; feedback noise bindings remain effect-owned.
MindSplatter's live particle count is readonly and excluded from validation and
interpolation. The eight authored effects preserve their parameter layouts,
names, order, bounds, labels, animation flags, presets and transition arithmetic.

Pullback's scalar arrays remain adapters over these descriptions. Their existing
field IDs, curves, catalog representation, renderer gates and topology metadata
remain in `core/render/pullback/fields.h`. Renderer gates decide registration;
shared control operations implement scalar validation and interpolation. The
operator catalog and generated pattern bytes remain unchanged.

`tests/test_param_marshal.h` covers typed/nested registration order, invalid
description bounds, discrete transitions, curve domains, excluded state, and
invalid snapshot ranges across all eight authored effects. Existing pullback,
choreography, preset, effect and capture tests pin their behavior.
