import { readFile } from 'node:fs/promises';
import {
  classifyExport,
  compileShaderDocument,
  parseShaderDocument,
} from './shader_workbench.mjs';

const usage = () => {
  console.error(
    'Usage:\n' +
    '  node scripts/shader_workbench_cli.mjs check <document.shader.json>\n' +
    '  node scripts/shader_workbench_cli.mjs descriptor <document.shader.json>\n' +
    '  node scripts/shader_workbench_cli.mjs classify <document.shader.json> <registry.json> <capability-profile>',
  );
};

const [command, documentPath, registryPath, capabilityProfile] = process.argv.slice(2);
if (!command || !documentPath) {
  usage();
  process.exitCode = 2;
} else {
  const source = await readFile(documentPath, 'utf8');
  const compiled = compileShaderDocument(source);
  if (command === 'check') {
    console.log(JSON.stringify({
      status: compiled.status,
      descriptor_digest: compiled.descriptor_digest,
      preset_bank_digest: compiled.preset_bank_digest,
      diagnostics: compiled.diagnostics,
    }, null, 2));
    if (compiled.status !== 'VALID') process.exitCode = 1;
  } else if (command === 'descriptor') {
    if (compiled.status === 'INVALID') {
      console.error(JSON.stringify(compiled.diagnostics, null, 2));
      process.exitCode = 1;
    } else {
      console.log(compiled.descriptor_json);
    }
  } else if (command === 'classify' && registryPath && capabilityProfile) {
    const registry = parseShaderDocument(await readFile(registryPath, 'utf8'));
    const result = classifyExport(compiled, registry, capabilityProfile);
    console.log(JSON.stringify(result, null, 2));
    if (result.kind === 'REJECTED') process.exitCode = 1;
  } else {
    usage();
    process.exitCode = 2;
  }
}
