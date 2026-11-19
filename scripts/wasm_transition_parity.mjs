import { PLANAR_SOURCE_OPERATOR, declarationFromCatalogField, evaluateTransition } from './shader_workbench.mjs';

const CAMERA_WALK_SEED = 1337;
const DOCUMENT_CHANNEL_TOLERANCE = 2;

/** Largest absolute RGB16 channel-code difference between two renders. */
export function peakChannelDifference(compiled, interpreter) {
  if (compiled.length === 0 || compiled.length !== interpreter.length)
    throw new Error('Compared framebuffers must have the same nonzero length');
  let peak = 0;
  for (let index = 0; index < compiled.length; ++index)
    peak = Math.max(peak, Math.abs(compiled[index] - interpreter[index]));
  return peak;
}

function transitionWrite({ id: name, storage, domain }, value) {
  if (storage !== 'enum8') return { name, value };
  if (typeof value === 'string') return { name, value: domain.values.indexOf(value) };
  return { name, value: {
    from: domain.values.indexOf(value.from), to: domain.values.indexOf(value.to), mix: value.mix,
  } };
}

const transitionWrites = (descriptor, values) =>
  descriptor.parameters.map(parameter => transitionWrite(parameter, values[parameter.id]));

function numericSourceMorph(descriptor, { from, to, ...morph }) {
  const endpoint = values => Object.fromEntries(descriptor.parameters
    .filter(({ id }) => Object.hasOwn(values, id))
    .map(parameter => {
      const { name, value } = transitionWrite(parameter, values[parameter.id]);
      return [name, value];
    }));
  return { ...morph, from: endpoint(from), to: endpoint(to) };
}

function checkLatticeBinding(bindings, require) {
  const admits = (program, from, to) => bindings.admitsPlanarSourceTransition(program, from, to);
  const source = { instance: 'sample', operator: PLANAR_SOURCE_OPERATOR };
  const affine = { instance: 'affine', operator: 'warp.affine.v4' };
  const grid = { 'sample.source': 0, 'sample.grid-pattern-freq': 2 };
  const bands = { 'sample.source': 2, 'sample.bands-pattern-freq': 2, 'sample.bands-angle': 0 };
  const spiral = { 'sample.source': 3 };
  require(admits([source], grid, spiral), 'An affine-free morph was refused');
  for (const warp of [{ 'affine.translation-x': 1 }, { 'affine.path-length': 1 }]) {
    require(admits([affine, source], { ...warp, ...grid }, { ...warp, ...bands }),
      'An affine morph with matching lattices was refused');
    require(!admits([affine, source], { ...warp, ...grid }, { ...warp, ...spiral }),
      'An affine morph with different lattices was admitted');
  }
}

/** Checks GridWarp preset identity and renders an edited ShaderChain transition. */
export function verifyWasmPresetTransition(compiledModule, interpreterModule, compiled, catalog) {
  const { document, preset_bank_digest: bankDigest } = compiled;
  const { descriptor: fixed, effect_id: effectId, preset_bank: authored } = document;
  const [origin] = authored.presets;
  const from = origin.preset_id;
  const to = 'preview-spiral';
  const sample = fixed.chain.find(slot => slot.operator === 'sample.grid.v4');
  const prefix = sample.label + '.';
  const sourceId = prefix + 'source';
  const mappingId = 'colorize.palette-mapping';
  const planar = catalog.operators.find(op => op.id === PLANAR_SOURCE_OPERATOR);
  const parameters = planar.params.filter(field => field.id === 'source' ||
    field.id.startsWith('grid-') || field.id.startsWith('spiral-'))
    .map(field => declarationFromCatalogField(sample.label, field, planar.id));
  const replaced = new Set(parameters.filter(({ id }) => id.startsWith(prefix + 'grid-'))
    .map(({ id }) => prefix + id.slice(prefix.length + 5)));
  const descriptor = { ...fixed,
    chain: fixed.chain.map(slot => slot === sample ? { ...slot, operator: planar.id } : slot),
    parameters: [...fixed.parameters.filter(({ id }) => !replaced.has(id)), ...parameters],
  };
  const values = Object.fromEntries(Object.entries(origin.values).filter(([id]) => !replaced.has(id)));
  for (const { id, default: value } of parameters) {
    const legacyId = prefix + id.slice(prefix.length + 5);
    values[id] = id.startsWith(prefix + 'grid-') ? origin.values[legacyId] : value;
  }
  values[sourceId] = 'grid';
  const target = { ...origin, preset_id: to, values: {
    ...values, [sourceId]: 'spiral', [mappingId]: 'bell',
  } };
  const edge = { from, to, duration: 8, easing: 'EASE_IN_OUT_SIN', path_policy: 'parallel' };
  const bank = { ...authored, presets: [{ ...origin, values }, target], edges: [edge] };
  const active = new compiledModule.HolosphereEngine();
  const interpreter = new interpreterModule.HolosphereEngine();
  let bindings;
  const require = (test, message) => {
    if (!test) throw new Error(message);
  };
  try {
    for (const [engine, module, effect] of [
      [active, compiledModule, effectId], [interpreter, interpreterModule, 'ShaderChain'],
    ]) {
      const result = engine.setResolution(96, 20);
      require(result === module.ResolutionSetResult.RESIZED || result === module.ResolutionSetResult.ALREADY_ACTIVE,
        'Could not configure the parity resolution');
      require(engine.setEffect(effect) === module.EffectSetResult.INSTALLED, `Could not install ${effect}`);
    }
    require(active.getPresetBankDigestForEffect(effectId) === bankDigest &&
      active.getPresetBankDigest() === bankDigest, 'Wasm preset bank differs from the document');
    require(JSON.stringify(active.getPresetIds()) === JSON.stringify(['grid-warp']) &&
      authored.presets.length === 1 && authored.edges.length === 0,
      'GridWarp must ship with one preset and no edges');
    require(!active.startPresetTransition(0, 1), 'GridWarp accepted an extra preset');

    bindings = interpreter.getShaderChainBindings();
    require(bindings, 'ShaderChain bindings unavailable');
    checkLatticeBinding(interpreterModule.ShaderChainBindings, require);
    const program = descriptor.chain.map(({ label: instance, operator }) => ({ instance, operator }));
    require(bindings.setShaderChain(program).code === 'APPLIED', 'Could not install the authored chain');
    require(bindings.setShaderChainParameters(transitionWrites(descriptor, values)) ===
      interpreterModule.ParamSetResult.APPLIED, 'Could not apply the chain origin');

    const snapshot = bindings.getSnapshot();
    const camera = snapshot.runtime.find(entry => entry.kind === 'spatial-walk-v2' && entry.instance === 'camera');
    require(camera, 'Chain camera clock unavailable');
    camera.state.noiseSeed = CAMERA_WALK_SEED;
    snapshot.animationsPaused = false;
    snapshot.paletteBank.chroma = values['colorize.palette-chroma'];
    require(bindings.restoreSnapshot(snapshot) === interpreterModule.ChainSnapshotRestoreResult.APPLIED,
      'Could not align the chain camera clock');
    interpreter.setAnimationsPaused(false);
    active.drawFrame();
    interpreter.drawFrame();
    const originDifference = peakChannelDifference(active.getPixels(), interpreter.getPixels());
    for (let frame = 1; frame <= edge.duration; ++frame) {
      const sample = evaluateTransition(descriptor, bank, from, to, frame, catalog);
      require(sample.sourceMorph, `Chain source endpoints unavailable at frame ${frame}`);
      require(bindings.setShaderChainSourceMorphParameters(transitionWrites(descriptor, sample.values),
        numericSourceMorph(descriptor, sample.sourceMorph)) ===
        interpreterModule.ParamSetResult.APPLIED, `Chain refused transition frame ${frame}`);
      require(!interpreter.getAnimationsPaused(), `Transition write paused animations at frame ${frame}`);
      interpreter.drawFrame();
    }
    interpreter.setTemporalSampling(false);
    const midpoint = evaluateTransition(descriptor, bank, from, to, Math.floor(edge.duration / 2), catalog).values;
    const originRender = (parameter, value) => {
      require(bindings.restoreSnapshot(snapshot) === interpreterModule.ChainSnapshotRestoreResult.APPLIED,
        'Could not restore the origin clock');
      require(bindings.setShaderChainTransitionParameters([transitionWrite(parameter, value)]) ===
        interpreterModule.ParamSetResult.APPLIED, 'Chain refused an origin test write');
      interpreter.drawFrame();
      return interpreter.getPixels().slice();
    };
    for (const id of [sourceId, mappingId]) {
      const parameter = descriptor.parameters.find(parameter => parameter.id === id);
      const held = originRender(parameter, values[id]);
      const blended = originRender(parameter, midpoint[id]);
      const incoming = originRender(parameter, target.values[id]);
      require(peakChannelDifference(blended, held) > DOCUMENT_CHANNEL_TOLERANCE &&
        peakChannelDifference(blended, incoming) > DOCUMENT_CHANNEL_TOLERANCE,
        `${id} blend rendered as a held endpoint`);
    }
    return { frames: edge.duration, originDifference };
  } finally {
    bindings?.delete();
    interpreter.delete();
    active.delete();
  }
}
