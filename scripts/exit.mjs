/** Drain pending diagnostics before terminating a command. */
export async function exitAfterStderr(code) {
  await new Promise((resolve) => process.stderr.write('', resolve));
  process.exit(code);
}
