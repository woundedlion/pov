# Display geometry

The display samples a complete mathematical sphere at the physical LED-center
latitudes. The north and south shaft caps contain no display rows. Row zero and
the last row are ordinary latitude rings when their centers are away from the
poles.

## Profiles and calibration

`HS_DISPLAY_PROFILE=1` selects the physical profile on firmware and WASM.
`HS_DISPLAY_PROFILE=0` selects an ideal pole-to-pole grid. The native regression
suite selects the ideal profile explicitly; its physical geometry executable
exercises the cropped profile separately.

The provisional physical calibration is 2% of the north-to-south arc at each
end: first LED center 3.6 degrees, last LED center 176.4 degrees. These are
assumptions pending measurement. Set `HS_DISPLAY_NORTH_FRACTION` and
`HS_DISPLAY_SOUTH_FRACTION` to the measured LED-center polar angles divided by
180 degrees. The defaults are `0.02f` and `0.98f`; unequal caps are supported.
The WASM CMake configuration exposes these three settings as cache variables.
Firmware accepts the same names as compiler definitions.

For H rows, north angle N and south angle S, the pitch is `(S-N)/(H-1)`.
Forward mapping is `phi = N + row*pitch`; inverse mapping is
`row = (phi-N)/pitch`. Calibration angles stay fixed when display resolution
changes. Neither conversion clamps. A direction inside a cap maps outside the
physical row range.

`math::DisplayGeometry` is the compile-time mapping and
`math::LatitudeGeometry` is its runtime counterpart. Pixel conversion, lookup
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

The WASM module exports `DISPLAY_PROFILE`, `DISPLAY_NORTH_PHI` and
`DISPLAY_SOUTH_PHI`. Daydream reads them before constructing LED geometry and
uses both angles in its matrix-cache identity. It does not select a different
geometry independently of the engine. Rebuild the module to switch between
ideal and physical profiles or change calibration.

## Validation

The physical regression executable checks endpoint-ring longitude separation,
forward/inverse mapping, asymmetric calibration, cap clipping, actual-pole
reflection, field sampling, feedback and angular scan bounds. Existing ideal
and legacy-offset tests cover their explicit compatibility profiles. Daydream
tests check matching endpoint placement, startup hydration and cache identity.
