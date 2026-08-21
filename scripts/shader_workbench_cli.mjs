import { readFile } from 'node:fs/promises';
import {
  ShaderDocumentError,
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

/** An invocation or input fault reported as one line, not as a stack trace. */
class CliError extends Error {}

const readText = async (path) => {
  try {
    return await readFile(path, 'utf8');
  } catch (cause) {
    throw new CliError(`cannot read ${path}: ${cause.message}`);
  }
};

// Root fields every shader document carries, v1 six-role and v2 chain alike.
// Other JSON under patterns/ — shaderball_migration.json, for one — would
// otherwise compile far enough to report its first stray key as a schema
// error, which reads as a broken document rather than as the wrong file.
const DOCUMENT_ROOT = ['descriptor', 'preset_bank'];

const readDocument = async (path) => {
  const document = parseShaderDocument(await readText(path));
  if (typeof document !== 'object' || document === null || Array.isArray(document) ||
      DOCUMENT_ROOT.some((field) => !(field in document)))
    throw new CliError(
      `${path} is not a shader document: a document's root carries ` +
      `${DOCUMENT_ROOT.join(' and ')}.`);
  return document;
};

const run = async () => {
  const [command, documentPath, registryPath, capabilityProfile] = process.argv.slice(2);
  if (!command || !documentPath) {
    usage();
    process.exitCode = 2;
    return;
  }
  const catalogPath = new URL('./engine_catalog.json', import.meta.url);
  const catalog = JSON.parse(await readText(catalogPath));
  const compiled = compileShaderDocument(await readDocument(documentPath), { catalog });
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
    const registry = parseShaderDocument(await readText(registryPath));
    const result = classifyExport(compiled, registry, capabilityProfile);
    console.log(JSON.stringify(result, null, 2));
    if (result.kind === 'REJECTED') process.exitCode = 1;
  } else {
    usage();
    process.exitCode = 2;
  }
};

try {
  await run();
} catch (error) {
  if (error instanceof CliError) console.error(error.message);
  else if (error instanceof ShaderDocumentError)
    console.error(`${error.code} at ${error.path}: ${error.message}`);
  else throw error;
  process.exitCode = 2;
}
