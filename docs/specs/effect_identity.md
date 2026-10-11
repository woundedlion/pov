# Effect identity contract

Status: IMPLEMENTED.

`HS_EFFECT_IDENTITY_LIST` in [`core/platform/effect_identity.h`](../../core/platform/effect_identity.h)
is the authoritative implementation-class, canonical ID, display alias, legacy ID
and immutable RNG-key roster. These fields are declared separately.
Each effect exposes the same typed `hs::EffectIdentity` as `EFFECT_IDENTITY`,
and its `EFFECT_ID` is a `constexpr std::string_view` equal to `machine_id`.
The identity is independent of render resolution.

Machine IDs use lowercase alphanumeric words separated by single hyphens,
begin with a lowercase letter and remain stable across class renames.
`display_name` retains the class/browser name used by current UI state and
the keys of `getEffectSizes()` and `getEffectPresetCounts()`.
`legacy_id` retains the original persisted ID. Factory lookup and `setEffect()`
accept the canonical ID, display/browser name and legacy ID with exact,
case-sensitive matching. Empty keys and prefixes do not match. Different
effects may not share any accepted lookup key; registration checks enforce
this namespace before allocation.

`rng_seed_key` is immutable and retains the original ID bytes. Installers
seed the stream through the typed identity or `effect_seed_key<E>()`, then
apply the unchanged `stable_effect_seed()` hash. Changing a machine ID or
lookup alias does not rekey the effect. Device, profile and WASM installers
use the same key; an unregistered profile type retains its declared ID or
class-name fallback.

`HolosphereEngine::getEffectIdentities()` exports an array of
`{id, display_name, legacy_ids}` objects generated from registrations.
`legacy_ids` contains distinct accepted aliases other than the canonical ID;
RNG keys are internal. Consumers can resolve existing saved values without
maintaining a second ID map. Existing class/browser map keys remain valid.

Composed shader documents use the canonical ID as `effect_id`.
[`scripts/generate_composed_presets.mjs`](../../scripts/generate_composed_presets.mjs)
reads the central roster and emits the wrapper's typed identity declaration.
Generation rejects missing identities, duplicate IDs and document/identity
mismatches.

[`tests/test_effect_factory.h`](../../tests/test_effect_factory.h) pins original
key and seed-hash snapshots, the identity type, aliases at each factory
resolution, installation stream restarts and the Phantasm seed table.
