export function moduleRosterFailures(tracked, loaded, exemptions) {
  const sources = new Set(tracked.filter((file) =>
    /\.(?:mjs|js)$/.test(file) && !/\.test\.(?:mjs|js)$/.test(file)));
  const failures = [];
  if (sources.size === 0) failures.push('no first-party modules selected');
  for (const [file, reason] of Object.entries(exemptions)) {
    if (typeof reason !== 'string' || reason.trim().length === 0)
      failures.push(`${file}: exemption needs a written reason`);
    if (!sources.has(file)) failures.push(`${file}: exemption is not a tracked source`);
    if (loaded.has(file)) failures.push(`${file}: loaded module no longer needs exemption`);
  }
  for (const file of sources) {
    if (!loaded.has(file) && !Object.hasOwn(exemptions, file))
      failures.push(`${file}: module was never loaded by a test`);
  }
  return failures;
}
