import { exitAfterStderr } from './exit.mjs';
import { writeFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';

const [modulePath, outputPath] = process.argv.slice(2);
if (!modulePath || !outputPath) {
  console.error('usage: export_engine_catalog.mjs <module.js> <catalog.json>');
  await exitAfterStderr(2);
}

try {
  const { default: createModule } = await import(pathToFileURL(modulePath));
  const module = await createModule();
  const catalog = module.HolosphereEngine.getShaderChainCatalog();
  const parsed = JSON.parse(catalog);
  if (!parsed || !Array.isArray(parsed.operators) || !parsed.operators.length
      || !Array.isArray(parsed.carriers) || !parsed.carriers.length || !parsed.budgets)
    throw new Error('engine returned an incomplete operator catalog');
  await writeFile(outputPath, `${catalog}\n`);
} catch (error) {
  console.error(`catalog export failed: ${error.message}`);
  process.exitCode = 1;
}
