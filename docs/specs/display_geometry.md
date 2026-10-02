# Display geometry

**Status: IMPLEMENTED.**

The display samples a complete mathematical sphere at the physical LED-center
latitudes. The north and south shaft caps contain no display rows. Row zero and
the last row are ordinary latitude rings when their centers are away from the
poles.

## Profiles and calibration

Firmware defaults to the physical profile (`HS_DISPLAY_PROFILE=1`).
`HS_DISPLAY_PROFILE=0` selects an ideal pole-to-pole grid. The native regression
suite selects the ideal profile explicitly; its physical geometry executable
exercises the cropped profile separately.

The provisional physical calibration is 2% of the north-to-south arc at each
end: first LED center 3.6 degrees, last LED center 176.4 degrees. These are
assumptions pending measurement. Set `HS_DISPLAY_NORTH_FRACTION` and
`HS_DISPLAY_SOUTH_FRACTION` to the measured LED-center polar angles divided by
180 degrees. The defaults are `0.02f` and `0.98f`; unequal caps are supported.
Firmware accepts these settings as compiler definitions. Daydream defaults to
full coverage and exposes separate Top cap (%) and Bottom cap (%) global
controls, from 0 to 25 percent each. Setting both to 2 previews the provisional
physical calibration. A zero cap places that endpoint at the mathematical pole.

For H rows, north angle N and south angle S, the pitch is `(S-N)/(H-1)`.
Forward mapping is `phi = N + row*pitch`; inverse mapping is
`row = (phi-N)/pitch`. Calibration angles stay fixed when display resolution
changes. Neither conversion clamps. A direction inside a cap maps outside the
physical row range.

`math::DisplayGeometry` supplies the resolution-specific mapping and
`math::LatitudeGeometry` carries runtime row bounds. Firmware uses compile-time
constants; WASM enables `HS_RUNTIME_DISPLAY_GEOMETRY` for live calibration. Pixel conversion, lookup
tables, angular bounds, rasterizer sampling and field geometry use this mapping.
The old explicit virtual-grid overloads and `HS_TEST_H_OFFSET` remain solely
for legacy callers and regression coverage.

## Visibility and poles

The framebuffer stores exactly H physical rows. Physical antialias splats
discard missing-row contributions without redistributing their brightness.
Footprints overlapping an endpoint LED can still light it. Geometry wholly
inside the missing region produces no pixels.

Crossing a framebuffer edge does not reflect a sample. Pole reflection occurs
only beyond polar angles zero and pi, and shifts longitude by half a turn.
Fractional sampling reflects before choosing its interpolation footprint.
Integer-only sampling rejects reflected locations that do not coincide with
the row lattice. Endpoint rings collapse longitude only when the selected
profile actually places them at a mathematical pole.

Spherical display fields and pixel feedback contain visible samples only;
missing samples contribute zero. World-space simulations remain independent
of this display crop. Dynamo's strand covers the full sphere using fractional
display coordinates, including hidden cap nodes.

## Simulator contract

The WASM module exports initial `DISPLAY_PROFILE`, `DISPLAY_NORTH_PHI` and
`DISPLAY_SOUTH_PHI` values. Daydream applies its global cap settings through
`setDisplayCaps(topPercent, bottomPercent)` and reads back the accepted angles
through `getDisplayNorthPhi()` and `getDisplaySouthPhi()`. Both controls default
to zero. Changes made while the module loads are applied before the first frame.

A geometry change refreshes the engine lookup tables and recreates the active
effect to invalidate geometry-dependent caches. Effect configuration, preset,
parameters, pause state and clipping are retained. Ordinary effects restart
animation and trail history. `ShaderChain` restores its complete executable
[snapshot](chain_snapshot_spec.md), including operator clocks, noise seeds and
generated palettes, after rebuilding geometry-dependent caches.
Invalid requests and unchanged values leave the running effect intact.
The controls persist across effect and resolution changes, propagate to segment
workers, and refresh the displayed LED mesh even while playback is paused.

## Validation

The physical regression executable checks endpoint-ring longitude separation,
forward/inverse mapping, asymmetric calibration, cap clipping, actual-pole
reflection, field sampling, feedback and angular scan bounds. Existing ideal
and legacy-offset tests cover their explicit compatibility profiles. Daydream
tests check matching endpoint placement, startup hydration and cache identity.
