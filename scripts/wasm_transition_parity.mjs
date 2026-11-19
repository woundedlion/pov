import { PLANAR_SOURCE_OPERATOR, evaluateTransition } from './shader_workbench.mjs';

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

/** Drives an authored GridWarp edge and measures the two runtimes frame by frame. */
export function verifyWasmPresetTransition(compiledModule, interpreterModule, compiled, catalog) {
  const { document, preset_bank_digest: bankDigest } = compiled;
  const { descriptor, effect_id: effectId, preset_bank: bank } = document;
  const [from, to] = bank.choreography.generated_order;
  const edge = bank.edges.find(item => item.from === from && item.to === to);
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
    require(JSON.stringify(active.getPresetIds()) === JSON.stringify(bank.choreography.generated_order),
      'Wasm preset IDs differ from generated order');
    require(!active.startPresetTransition(0, 2), 'Compiled effect accepted a non-successor edge');
    require(edge && edge.duration > 1, 'The pilot needs a timed departure');
    const cameraIndex = active.getParameterDefinitions().findIndex(({ name }) => name === 'Camera Wander');
    require(cameraIndex >= 0, 'Compiled camera control unavailable');

    bindings = interpreter.getShaderChainBindings();
    require(bindings, 'ShaderChain bindings unavailable');
    checkLatticeBinding(interpreterModule.ShaderChainBindings, require);
    const program = descriptor.chain.map(({ label: instance, operator }) => ({ instance, operator }));
    require(bindings.setShaderChain(program).code === 'APPLIED', 'Could not install the authored chain');
    const { values } = bank.presets.find(item => item.preset_id === from);
    const source = descriptor.chain.find(slot => slot.operator === PLANAR_SOURCE_OPERATOR);
    const target = bank.presets.find(item => item.preset_id === to);
    const sourceId = source && `${source.label}.source`;
    require(sourceId && values[sourceId] !== target.values[sourceId],
      'The pilot departure must cross source kinds');
    const mappingId = Object.keys(values).find(id => id.endsWith('.palette-mapping'));
    require(mappingId && values[mappingId] !== target.values[mappingId],
      'The pilot departure must cross palette mappings');
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
    require(active.startPresetTransition(0, 1), 'Could not arm the compiled departure');

    let peak = originDifference, peakFrame = 0;
    let firstMismatch = originDifference > DOCUMENT_CHANNEL_TOLERANCE ? 0 : null;
    let peakCameraDifference = 0;
    for (let frame = 1; frame <= edge.duration; ++frame) {
      const sample = evaluateTransition(descriptor, bank, from, to, frame, catalog);
      require(bindings.setShaderChainTransitionParameters(transitionWrites(descriptor, sample.values)) ===
        interpreterModule.ParamSetResult.APPLIED, `Chain refused transition frame ${frame}`);
      require(!interpreter.getAnimationsPaused(), `Transition write paused animations at frame ${frame}`);
      active.drawFrame();
      interpreter.drawFrame();
      peakCameraDifference = Math.max(peakCameraDifference,
        Math.abs(active.getParamValues()[cameraIndex] - sample.values['camera.wander']));
      const difference = peakChannelDifference(active.getPixels(), interpreter.getPixels());
      if (difference > peak) {
        peak = difference;
        peakFrame = frame;
      }
      if (firstMismatch === null && difference > DOCUMENT_CHANNEL_TOLERANCE) firstMismatch = frame;
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
    return { frames: edge.duration, peak, peakFrame, firstMismatch, peakCameraDifference,
      tolerance: DOCUMENT_CHANNEL_TOLERANCE };
  } finally {
    bindings?.delete();
    interpreter.delete();
    active.delete();
  }
}
