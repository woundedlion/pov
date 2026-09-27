import { test } from 'node:test';
import assert from 'node:assert/strict';
import { descendToHonoredResolution, loadEffectForCapture } from './screenshot_resolution.mjs';

// Returns a load() stub that reports `offers[resolution]` as the selected
// effect, recording the order it was called in.
function stubLoader(offers, calls = []) {
  return async (effect, resolution) => {
    calls.push(resolution);
    return Object.prototype.hasOwnProperty.call(offers, resolution)
      ? offers[resolution] : 'FallbackEffect';
  };
}

test('stops at the first resolution that offers the effect', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingSpin', ['big', 'mid', 'small'],
    stubLoader({ big: 'RingSpin', mid: 'RingSpin' }, calls));
  assert.deepEqual(result, { resolution: 'big', honored: true });
  assert.deepEqual(calls, ['big']);
});

test('descends past resolutions the app falls back on', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingShower', ['big', 'mid', 'small'],
    stubLoader({ small: 'RingShower' }, calls));
  assert.deepEqual(result, { resolution: 'small', honored: true });
  assert.deepEqual(calls, ['big', 'mid', 'small']);
});

test('reports no honored resolution when every load falls back', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'Ghost', ['big', 'small'], stubLoader({}, calls));
  assert.deepEqual(result, { resolution: 'small', honored: false });
  assert.deepEqual(calls, ['big', 'small']);
});

test('an empty resolution list honors nothing and loads nothing', async () => {
  const calls = [];
  const result = await descendToHonoredResolution(
    'RingSpin', [], stubLoader({}, calls));
  assert.deepEqual(result, { resolution: null, honored: false });
  assert.deepEqual(calls, []);
});

test('a page reporting no effect param is never honored', async () => {
  const result = await descendToHonoredResolution(
    'RingSpin', ['big', 'small'], stubLoader({ big: null, small: undefined }));
  assert.equal(result.honored, false);
});

test('a name the requested effect merely prefixes is not honored', async () => {
  const result = await descendToHonoredResolution(
    'Ring', ['big'], stubLoader({ big: 'RingSpin' }));
  assert.equal(result.honored, false);
});

test('capture descent retries a selection that changes during the capture offset', async () => {
  const visits = [];
  let resolution;
  let settled;
  const page = {
    async goto(url) {
      resolution = new URL(url).searchParams.get('resolution');
      visits.push(resolution);
      settled = false;
    },
    async waitForSelector() {},
    async waitForFunction() {},
    async evaluate() { return resolution === 'big' && settled ? 'Fallback' : 'Dynamo'; },
    async waitForTimeout() { settled = true; },
  };
  const result = await descendToHonoredResolution('Dynamo', ['big', 'small'],
    (effect, size) => loadEffectForCapture(page, 'http://fixture/', effect, size, 1));
  assert.deepEqual(result, { resolution: 'small', honored: true });
  assert.deepEqual(visits, ['big', 'small']);
});

test('capture descent rejects a hydrated fallback without waiting for capture', async () => {
  const page = {
    async goto() {},
    async waitForSelector() {},
    async waitForFunction() {},
    async evaluate() { return null; },
    async waitForTimeout() { assert.fail('unselected effects must not wait for capture'); },
  };
  assert.equal(await loadEffectForCapture(page, 'http://fixture/', 'Dynamo', 'big', 1), null);
});

if (process.env.HS_BROWSER_TESTS === '1') {
  test('browser capture descends using hydrated selection while the URL still names the request', async () => {
    const { chromium } = await import('playwright');
    const browser = await chromium.launch({ headless: true });
    try {
      const page = await browser.newPage();
      const visits = [];
      await page.route('http://fixture/**', async (route) => {
        const resolution = new URL(route.request().url()).searchParams.get('resolution');
        visits.push(resolution);
        const selected = resolution === 'big' ? 'Fallback' : 'Dynamo';
        await route.fulfill({ contentType: 'text/html', body: `
          <canvas id="canvas"></canvas><div id="loading-overlay"></div>
          <button class="effect-button" aria-selected="false" data-effect="${selected}"></button>
          <script>
            setTimeout(() => {
              document.querySelector('button').setAttribute('aria-selected', 'true');
              document.getElementById('loading-overlay').remove();
            }, 20);
          </script>` });
      });
      const result = await descendToHonoredResolution('Dynamo', ['big', 'small'],
        (effect, size) => loadEffectForCapture(page, 'http://fixture/', effect, size, 1));
      assert.deepEqual(result, { resolution: 'small', honored: true });
      assert.deepEqual(visits, ['big', 'small']);
      assert.equal(new URL(page.url()).searchParams.get('effect'), 'Dynamo');
    } finally {
      await browser.close();
    }
  });
}
