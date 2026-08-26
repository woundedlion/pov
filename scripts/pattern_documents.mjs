// Discovery and compilation of the committed patterns/ shader documents, shared
// by the WASM smoke gate, the promoted-digest test and the generator's
// canonicality check.
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
 * Compiles every committed pattern document, in directory order.
 *
 * The caller decides what a non-VALID compile means: the smoke gate throws, the
 * digest test asserts, the generator collects it as noncanonical.
 *
 * @param {Object} catalog Operator catalog to compile against.
 * @returns {Promise<{name: string, source: string, compiled: Object}[]>} One
 *   entry per patterns/*.shader.json, with the source LF-normalized.
 */
export async function compilePatternDocuments(catalog) {
  const names = (await readdir(PATTERNS_DIR))
    .filter((name) => name.endsWith('.shader.json'));
  const documents = [];
  for (const name of names) {
    // The committed blobs are LF; core.autocrlf hands a Windows checkout CRLF.
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
