# Pullback temporal sampling

Status: IMPLEMENTED for GridSpace.

## Timing contract

`Pullback::Spec::TEMPORAL_SAMPLING` selects column-time rendering. The default
is false. GridSpace enables it.

For an even canvas width `W`, absolute column `x` samples the interval from the
previous animation endpoint to the current endpoint at
`(x % (W / 2)) / float(W / 2)`. Opposing arms therefore share a time, and a
clipped renderer uses the same phase as a full-frame renderer. The interval is
one effect frame; the renderer introduces up to one frame of animation latency.
On the first draw, both endpoints are the current state.

The next sweep begins at the current endpoint. Temporal progression across the
scan boundary is consequently one column interval, rather than a frame jump.
This describes the rotating display; a flat snapshot still contains columns
from different times.

## State and resources

The composed runtime advances animation once per frame. It then prepares a
pipeline once per visible column and shades its rows using that pipeline.
Geometry is evaluated once per pixel. Quaternions use spherical interpolation,
parameters use their declared curves, and periodic clocks take the shortest arc.

The previous endpoint owns copies of its palette and hue tables. Current tables
are baked once per frame. At intermediate times, the colorizer evaluates both
resource endpoints at the column's interpolated geometry and controls, then
blends the resulting colors. It does not rebuild noise tables per column.
Effects without the opt-in allocate no history and retain frame-wide shading.

## Opt-in requirements

Clock advances must remain below half their period, and wrapping a phase must
preserve the spatial field. GridSpace's whole-cell affine winding meets this
condition. Continuous parameter curves are needed for continuous changes;
discrete controls and held Mobius coefficients remain discontinuous. A changing
periodic-ripple period needs a separate wrap-continuity review before opt-in.

## Verification

The composed-effect tests check opposing-arm timing, forward and reverse clock
rollover, quaternion interpolation, retained resource ownership, and actual
GridSpace pixels under a right-half clip. Color-policy tests check both endpoints
and intermediate palette, hue, and alpha interpolation. Firmware size/layout
and native arena/stack gates cover the added storage and prepared state.
