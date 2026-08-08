// Dependency-free PNG validator for the docs/screenshots gallery gate. Walks
// the chunk stream, checks every CRC, and inflates the IDAT payload, so a
// zero-byte, truncated, mislabelled, or bit-rotted file cannot pass as an image.
import { inflateSync } from 'node:zlib';

const SIGNATURE = Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
// PNG color type -> legal bit depths and samples per pixel (1 and 5 do not exist).
const COLOR_TYPES = new Map([
  [0, { depths: [1, 2, 4, 8, 16], channels: 1 }],
  [2, { depths: [8, 16], channels: 3 }],
  [3, { depths: [1, 2, 4, 8], channels: 1 }],
  [4, { depths: [8, 16], channels: 2 }],
  [6, { depths: [8, 16], channels: 4 }],
]);

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

function readHeader(data) {
  if (data.length !== 13) throw new Error('IHDR is not 13 bytes');
  const width = data.readUInt32BE(0);
  const height = data.readUInt32BE(4);
  const depth = data[8];
  const colorType = data[9];
  const spec = COLOR_TYPES.get(colorType);
  if (width === 0 || height === 0) throw new Error('IHDR declares a zero dimension');
  if (spec === undefined) throw new Error(`IHDR color type ${colorType} is not a PNG color type`);
  if (!spec.depths.includes(depth))
    throw new Error(`IHDR bit depth ${depth} is not legal for color type ${colorType}`);
  return { width, height, depth, colorType, channels: spec.channels, interlace: data[12] };
}

/**
 * Validates a PNG datastream and returns its pixel dimensions.
 * @param {Buffer} bytes Whole file contents.
 * @returns {{width: number, height: number}} Dimensions declared by IHDR.
 * @throws {Error} If the signature, chunk layout, a CRC, or the pixel data is malformed.
 */
export function inspectPng(bytes) {
  if (bytes.length < SIGNATURE.length || !bytes.subarray(0, SIGNATURE.length).equals(SIGNATURE))
    throw new Error('not a PNG (bad signature)');

  let offset = SIGNATURE.length;
  let header = null;
  let ended = false;
  let palette = false;
  const pixelChunks = [];
  while (!ended) {
    if (offset === bytes.length) throw new Error('no IEND chunk');
    if (offset + 12 > bytes.length) throw new Error('truncated chunk header');
    const length = bytes.readUInt32BE(offset);
    if (offset + 12 + length > bytes.length) throw new Error('truncated chunk data');
    const type = bytes.toString('latin1', offset + 4, offset + 8);
    const data = bytes.subarray(offset + 8, offset + 8 + length);
    if (crc32(bytes.subarray(offset + 4, offset + 8 + length))
        !== bytes.readUInt32BE(offset + 8 + length))
      throw new Error(`chunk ${type} fails its CRC`);
    if (header === null && type !== 'IHDR') throw new Error('first chunk is not IHDR');
    if (type === 'IHDR') {
      if (header !== null) throw new Error('duplicate IHDR chunk');
      header = readHeader(data);
    } else if (type === 'PLTE') {
      palette = true;
    } else if (type === 'IDAT') {
      pixelChunks.push(data);
    } else if (type === 'IEND') {
      ended = true;
    }
    offset += 12 + length;
  }
  if (offset !== bytes.length) throw new Error('trailing bytes after IEND');
  if (pixelChunks.length === 0) throw new Error('no IDAT pixel data');
  if (header.colorType === 3 && !palette) throw new Error('palette image has no PLTE chunk');

  let raw;
  try {
    raw = inflateSync(Buffer.concat(pixelChunks));
  } catch {
    throw new Error('IDAT pixel data does not inflate');
  }
  // Adam7 splits the raster into seven passes with their own row padding; the
  // gallery is never interlaced, so only the flat layout is size-checked.
  if (header.interlace === 0) {
    const stride = Math.ceil(header.width * header.channels * header.depth / 8);
    const expected = header.height * (stride + 1);
    if (raw.length !== expected)
      throw new Error(`inflated pixel data is ${raw.length} bytes, expected ${expected}`);
  }
  return { width: header.width, height: header.height };
}
