// Discovery and compilation of the committed patterns/ shader documents.
import { readdir, readFile } from 'node:fs/promises';
import { compileShaderDocument } from './shader_workbench.mjs';

const PATTERNS_DIR = new URL('../patterns/', import.meta.url);
const CATALOG_URL = new URL('./engine_catalog.json', import.meta.url);

/**
 * @returns {Promise<Object>} The wasm32 operator catalog every committed
 *   document is compiled against.
 */
export async function loadOperatorCatalog() {
  return JSON.parse(await readFile(CATALOG_URL, 'utf8'));
}

/**
 * Compiles every committed pattern document, in sorted filename order.
 *
 * Non-VALID compiles are returned, not thrown.
 *
 * @param {Object} catalog Operator catalog to compile against.
 * @returns {Promise<{name: string, source: string, compiled: Object}[]>} One
 *   entry per patterns/*.shader.json, with the source LF-normalized.
 */
export async function compilePatternDocuments(catalog) {
  const names = (await readdir(PATTERNS_DIR))
    .filter((name) => name.endsWith('.shader.json')).sort();
  const documents = [];
  for (const name of names) {
    const source = (await readFile(new URL(name, PATTERNS_DIR), 'utf8'))
      .replaceAll('\r\n', '\n');
    documents.push({
      name,
      source,
      compiled: compileShaderDocument(source, { catalog }),
    });
  }
  return documents;
}
