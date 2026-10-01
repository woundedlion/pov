import { realpathSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

export function isMain(url) {
  return Boolean(process.argv[1])
    && realpathSync(fileURLToPath(url)) === realpathSync(process.argv[1]);
}

/** Drain pending diagnostics before terminating a command. */
export async function exitAfterStderr(code) {
  await new Promise((resolve) => process.stderr.write('', resolve));
  process.exit(code);
}
