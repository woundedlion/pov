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

The simulator's `HolosphereEngine.setTemporalSampling` selects this behavior
at runtime. With POV fidelity disabled, every column uses the current endpoint.
With fidelity enabled, columns use their emission times and the simulator applies
its visual response model. History advances in either mode, so toggling does not
leave an old endpoint behind. Each new engine defaults to temporal sampling;
the simulator sets its selected mode before the first frame. Physical firmware
always enables temporal sampling for opted-in effects.

## State and resources

The composed runtime advances animation once per frame. It then prepares a
pipeline once per visible column and shades its rows using that pipeline.
Geometry is evaluated once per pixel. Quaternions use spherical interpolation,
parameters use their declared curves, and periodic clocks take the shortest arc.
Quaternion arc setup and periodic clock deltas are prepared once per frame.
Unchanged parameter blocks share an interior sample of their declared curves.

`Spec::FACTORED_COLUMNS` separately enables a prepared meridian basis for
pipelines beginning with rotation and projection. GridSpace enables it for both
temporal and instantaneous rendering. Two transformed basis vectors replace
the two quaternion rotations at each pixel. The remaining stages still run per
pixel; custom effect shading overrides retain their original path. This changes
floating-point evaluation order. Tests bound GridSpace's authored preset to
16 RGB16 codes maximum error and 0.05 mean error against direct rotation.

The previous endpoint owns copies of its palette and hue tables. Current tables
are refreshed at frame boundaries. At intermediate times, the colorizer evaluates
both resource endpoints at the column's interpolated geometry and controls, then
blends the resulting colors. It does not rebuild noise tables per column.
Effects without the opt-in allocate no history and retain frame-wide shading.

Temporal effects using noise-driven hue automatically keep float copies of the
quantized hue and noise tables. Conversion occurs when a table is baked; history
owns the corresponding previous copies. The float copies preserve the stored
channel and noise values and avoid converting each pixel's lookup taps. The
original tables remain available to ordinary endpoint and fallback sampling.
Constant palette alpha is detected once per frame; varying alpha still samples
the palette. These resources use the existing persistent arena budget.

The temporal colorizer shares mapping, opacity, cube addressing and palette
value addressing between endpoints. Column preparation selects separate noise
paths for unit gain and brightness bounded to `[0, 1]`. These paths keep the Q16
sampling coordinates, combine floating-point RGB samples, and round once at
the output. Gains outside that interval retain per-endpoint rounding and
clamping. The color-policy tests allow at most two RGB16 codes of difference
from independently rounded endpoint sampling. Alpha retains endpoint
interpolation; compiler reassociation can change its last floating-point bits.
Endpoint fractions retain the ordinary single-endpoint colorizer.

Noise tables are compared once per frame. Equal tables share one noise sample
and its hue coordinates across the two palette lookups. Distinct tables retain
independent noise samples. The history continues to own its tables in both
cases. Equal constant endpoint alpha is also shared by the column colorizer.

## Opt-in requirements

Clock advances must remain below half their period, and wrapping a phase must
preserve the spatial field. GridSpace's whole-cell affine winding meets this
condition. Continuous parameter curves are needed for continuous changes;
discrete controls and held Mobius coefficients remain discontinuous. A changing
periodic-ripple period needs a separate wrap-continuity review before opt-in.

## Verification

The composed-effect tests check opposing-arm timing, forward and reverse clock
rollover, quaternion interpolation, retained resource ownership, and actual
GridSpace pixels under a clip crossing an opposing-arm boundary. They also check
instantaneous rendering and resuming temporal sampling after disabling it.
Color-policy tests check both endpoints
and intermediate palette, hue, and alpha interpolation. Firmware size/layout
and native arena/stack gates cover the added storage and prepared state.
