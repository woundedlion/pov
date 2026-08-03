import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import { deflateSync } from 'node:zlib';
import { REPO_ROOT } from './effect_roster.mjs';
import { inspectPng } from './png_probe.mjs';
import { GALLERY_HEIGHT, GALLERY_WIDTH } from './screenshot_capture_config.mjs';

const SIGNATURE = Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);

const CRC_TABLE = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
  CRC_TABLE[n] = c >>> 0;
}
function crc32(bytes) {
  let c = 0xffffffff;
  for (let i = 0; i < bytes.length; i++)
    c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

function chunk(type, data) {
  const head = Buffer.alloc(8);
  head.writeUInt32BE(data.length, 0);
  head.write(type, 4, 'latin1');
  const crc = Buffer.alloc(4);
  crc.writeUInt32BE(crc32(Buffer.concat([head.subarray(4), data])), 0);
  return Buffer.concat([head, data, crc]);
}

// Minimal 8-bit RGBA image with every scanline filtered as None.
function makePng(width, height) {
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(width, 0);
  ihdr.writeUInt32BE(height, 4);
  ihdr[8] = 8;
  ihdr[9] = 6;
  const raw = Buffer.alloc(height * (1 + width * 4));
  return Buffer.concat([
    SIGNATURE,
    chunk('IHDR', ihdr),
    chunk('IDAT', deflateSync(raw)),
    chunk('IEND', Buffer.alloc(0)),
  ]);
}

test('a well-formed PNG reports its declared dimensions', () => {
  assert.deepEqual(inspectPng(makePng(7, 3)), { width: 7, height: 3 });
});

test('an empty file is rejected', () => {
  assert.throws(() => inspectPng(Buffer.alloc(0)), /bad signature/);
});

test('a text file named .png is rejected', () => {
  assert.throws(() => inspectPng(Buffer.from('not an image, just prose\n')),
    /bad signature/);
});

test('a truncated datastream is rejected', () => {
  const png = makePng(7, 3);
  assert.throws(() => inspectPng(png.subarray(0, png.length - 20)),
    /truncated|no IEND/);
});

test('a flipped pixel byte is caught by the chunk CRC', () => {
  const png = makePng(7, 3);
  const idat = png.indexOf('IDAT', 0, 'latin1');
  png[idat + 6] ^= 0xff;
  assert.throws(() => inspectPng(png), /fails its CRC/);
});

test('a header-only file with no pixel data is rejected', () => {
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(4, 0);
  ihdr.writeUInt32BE(4, 4);
  ihdr[8] = 8;
  ihdr[9] = 6;
  const png = Buffer.concat([
    SIGNATURE, chunk('IHDR', ihdr), chunk('IEND', Buffer.alloc(0)),
  ]);
  assert.throws(() => inspectPng(png), /no IDAT/);
});

test('pixel data shorter than the declared raster is rejected', () => {
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(7, 0);
  ihdr.writeUInt32BE(3, 4);
  ihdr[8] = 8;
  ihdr[9] = 6;
  const png = Buffer.concat([
    SIGNATURE,
    chunk('IHDR', ihdr),
    chunk('IDAT', deflateSync(Buffer.alloc(4))),
    chunk('IEND', Buffer.alloc(0)),
  ]);
  assert.throws(() => inspectPng(png), /inflated pixel data is 4 bytes/);
});

test('appended trailing bytes are rejected', () => {
  const png = Buffer.concat([makePng(7, 3), Buffer.from('junk')]);
  assert.throws(() => inspectPng(png), /trailing bytes after IEND/);
});

test('a zero dimension in the header is rejected', () => {
  assert.throws(() => inspectPng(makePng(0, 3)), /zero dimension/);
});

test('every committed gallery PNG decodes at the stored size', async () => {
  for (const effect of ['IslamicStars', 'Thrusters', 'Voronoi']) {
    const file = join(REPO_ROOT, 'docs', 'screenshots', `${effect}.png`);
    assert.deepEqual(inspectPng(await readFile(file)),
      { width: GALLERY_WIDTH, height: GALLERY_HEIGHT }, effect);
  }
});
