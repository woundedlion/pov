import { writeFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';

const [jsPath, resolution, outputPath] = process.argv.slice(2);
if (!jsPath || !resolution || !outputPath) {
  throw new Error('usage: pullback_capture_wasm.mjs module.js WxH output.bin');
}
const match = /^(96)x(20)$|^(288)x(144)$/.exec(resolution);
if (!match) throw new Error(`unsupported pullback capture resolution: ${resolution}`);
const [width, height] = resolution.split('x').map(Number);
const { default: createHolosphereModule } = await import(pathToFileURL(jsPath));
const Module = await createHolosphereModule({
  print: () => {},
  printErr: (message) => console.error(`[wasm:err] ${message}`),
});
const engine = new Module.HolosphereEngine();
const output = [];
const u16 = (value) => output.push(value & 0xff, (value >>> 8) & 0xff);
const u32 = (value) => { u16(value); u16(value >>> 16); };
for (const byte of Buffer.from('HSPB')) output.push(byte);
u16(1); u16(width); u16(height); u16(12); u16(4); u32(48);

const resolutionResult = engine.setResolution(width, height);
if (resolutionResult !== Module.ResolutionSetResult.RESIZED &&
    resolutionResult !== Module.ResolutionSetResult.ALREADY_ACTIVE) {
  throw new Error(`setResolution(${resolution}) failed`);
}
const fractions = [0.25, 0.5, 0.75];
for (let preset = 0; preset < 12; preset++) {
  for (let caseIndex = 0; caseIndex < 4; caseIndex++) {
    if (engine.setEffect('ShaderBall') !== Module.EffectSetResult.INSTALLED) {
      throw new Error('setEffect(ShaderBall) failed');
    }
    if (!engine.selectPreset(preset)) throw new Error(`selectPreset(${preset}) failed`);
    engine.setAnimationsPaused(true);
    if (caseIndex !== 0) {
      const parameters = engine.getParameterDefinitions().filter(
        (parameter) => !parameter.readonly && parameter.options === undefined);
      for (let index = 0; index < parameters.length; index++) {
        const parameter = parameters[index];
        let value;
        if (typeof parameter.requestedValue === 'boolean') {
          value = caseIndex === 1 ? 0 : caseIndex === 2 ? 1 :
            fractions[index % fractions.length] > 0.5 ? 1 : 0;
        } else {
          const fraction = caseIndex === 1 ? 0 : caseIndex === 2 ? 1 :
            fractions[index % fractions.length];
          value = parameter.min + fraction * (parameter.max - parameter.min);
        }
        if (engine.setParameter(parameter.name, value) !== Module.ParamSetResult.APPLIED) {
          throw new Error(`setParameter(${parameter.name}) failed`);
        }
      }
    }
    engine.drawFrame();
    const pixels = engine.getPixels();
    if (pixels.length !== width * height * 3) throw new Error('pixel span mismatch');
    u16(preset); u16(caseIndex);
    for (const channel of pixels) u16(channel);
  }
}
engine.delete();
await writeFile(outputPath, Buffer.from(output));
