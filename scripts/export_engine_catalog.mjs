import { writeFile } from 'node:fs/promises';
import { pathToFileURL } from 'node:url';

const [modulePath, outputPath] = process.argv.slice(2);
if (!modulePath || !outputPath) {
  console.error('usage: export_engine_catalog.mjs <module.js> <catalog.json>');
  process.exit(2);
}

const { default: createModule } = await import(pathToFileURL(modulePath));
const module = await createModule();
const catalog = module.HolosphereEngine.getShaderChainCatalog();
JSON.parse(catalog);
await writeFile(outputPath, `${catalog}\n`);
