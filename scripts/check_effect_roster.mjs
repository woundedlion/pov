// Cross-check the standard HS_EFFECT_LIST roster against the REGISTER_EFFECT
// calls under effects/. The standalone WASM-only Shader workbench is excluded.
// The WASM startup check only
// compares the registry SIZE to HS_EFFECT_COUNT, so an effect whose header is
// never #included from core/engine/effects.h (no #include, no X() row) is absent from
// BOTH the registry and the count — they still agree and nothing fails, while
// the native smoke suite (driven by the X-macro list) silently never runs it.
// Scanning the headers on disk catches that: a REGISTER_EFFECT with no X() row,
// or an X() row with no REGISTER_EFFECT, fails here.
import { loadEffectRoster, loadRegisteredEffects } from './effect_roster.mjs';

const roster = new Set(await loadEffectRoster());
const registered = new Set(await loadRegisteredEffects());

const unlisted = [...registered].filter(e => !roster.has(e)); // REGISTER_EFFECT, no X() row
const unregistered = [...roster].filter(e => !registered.has(e)); // X() row, no REGISTER_EFFECT

if (unlisted.length || unregistered.length) {
  if (unlisted.length)
    console.error(`::error::effects/**/*.h register standard effects absent from HS_EFFECT_LIST (add an #include + X() row in core/engine/effects.h): ${unlisted.join(', ')}`);
  if (unregistered.length)
    console.error(`::error::HS_EFFECT_LIST names effects with no REGISTER_EFFECT under effects/: ${unregistered.join(', ')}`);
  process.exitCode = 1;
} else {
  console.log(`Effect roster matches registrations: all ${roster.size} effects listed and registered.`);
}
