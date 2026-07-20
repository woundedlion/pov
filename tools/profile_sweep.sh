#!/bin/bash
# profile_sweep.sh <group: g1_ship..g4_ship, g1_o3..g4_o3>
# Phantasm-roster profiling sweep, one group per invocation (see the
# teensy-profile skill for the per-effect duration/knob rationale).
# Covers the 22 effects in the 288x144 Phantasm playlist. Dynamo and Thrusters
# are Holosphere 96x20-only (HS_PHANTASM_EFFECT_LIST excludes them), so they are
# not profiled here.
# Sequential profile_one.sh calls per group. A failed capture is recorded and
# the group carries on: the board can drop USB mid-capture, and aborting the
# group turns one dropped effect into every later one missing (a single hang
# once cost five O3 captures). The group exits non-zero listing what failed, so
# a re-run only needs those effects.
# A capture must fit inside one epoch: crossing a boundary re-inits the effect
# mid-run, and an init that overruns the K-revolution commit window traps the
# board. The RD sims want 130 s of data, past the 120 s default, so they stretch
# the epoch to 150 s rather than lose a regime to a shorter capture.
set -u
P="$(dirname "$0")/profile_one.sh"
FAILED=()

# run <Effect> <env> <seconds> <window> [extra flags]
run() {
  if ! bash "$P" "$@"; then
    echo ">>> FAILED: $1 [$2] — continuing with the rest of the group" >&2
    FAILED+=("$1/$2")
  fi
}
case "$1" in
g1_ship)
  run BZReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  run ChaoticStrings profile 70 32
  run DisplacementField profile 70 32
  run GnomonicStars profile 70 32
  run GSReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  ;;
g2_ship)
  run HopfFibration profile 70 32
  run MobiusGrid profile 70 32
  run PetalFlow profile 70 32
  run Raymarch profile 70 32
  run RingShower profile 70 32
  run RingSpin profile 70 32
  run Voronoi profile 70 32
  ;;
g3_ship)
  run Liquid2D profile 70 16
  run ShapeShifter profile 70 16
  run MindSplatter profile 110 16
  run HankinSolids profile 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"
  run SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"
  ;;
g4_ship)
  run Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run MeshFeedback profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run Flyby profile 310 16 "-D HS_PROFILE_EPOCH_REVS=2560"
  run IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  run DreamBalls profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  ;;
g1_o3)
  run BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  run ChaoticStrings profile_o3 70 32
  run DisplacementField profile_o3 70 32
  run GnomonicStars profile_o3 70 32
  run GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  ;;
g2_o3)
  run HopfFibration profile_o3 70 32
  run MobiusGrid profile_o3 70 32
  run PetalFlow profile_o3 70 32
  run Raymarch profile_o3 70 32
  run RingShower profile_o3 70 32
  run RingSpin profile_o3 70 32
  run Voronoi profile_o3 70 32
  ;;
g3_o3)
  run Liquid2D profile_o3 70 16
  run ShapeShifter profile_o3 70 16
  run MindSplatter profile_o3 110 16
  run HankinSolids profile_o3 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"
  run SphericalHarmonics profile_o3 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"
  ;;
g4_o3)
  run Comets profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run MeshFeedback profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run Flyby profile_o3 310 16 "-D HS_PROFILE_EPOCH_REVS=2560"
  run IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  run DreamBalls profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  ;;
*) echo "unknown group $1"; exit 1;;
esac
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "GROUP $1 INCOMPLETE — ${#FAILED[@]} failed: ${FAILED[*]}"
  exit 1
fi
echo "GROUP $1 DONE"
