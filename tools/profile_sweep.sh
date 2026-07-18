#!/bin/bash
# profile_sweep.sh <group: g1_ship..g4_ship, g1_o3..g4_o3>
# Phantasm-roster profiling sweep, one group per invocation (see the
# teensy-profile skill for the per-effect duration/knob rationale).
# Covers the 22 effects in the 288x144 Phantasm playlist. Dynamo and Thrusters
# are Holosphere 96x20-only (HS_PHANTASM_EFFECT_LIST excludes them), so they are
# not profiled here.
# Sequential profile_one.sh calls per group; stops on first failure.
# A capture must fit inside one epoch: crossing a boundary re-inits the effect
# mid-run, and an init that overruns the K-revolution commit window traps the
# board. The RD sims want 130 s of data, past the 120 s default, so they stretch
# the epoch to 150 s rather than lose a regime to a shorter capture.
set -e
P="$(dirname "$0")/profile_one.sh"
case "$1" in
g1_ship)
  bash "$P" BZReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  bash "$P" ChaoticStrings profile 70 32
  bash "$P" DisplacementField profile 70 32
  bash "$P" GnomonicStars profile 70 32
  bash "$P" GSReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  ;;
g2_ship)
  bash "$P" HopfFibration profile 70 32
  bash "$P" MobiusGrid profile 70 32
  bash "$P" PetalFlow profile 70 32
  bash "$P" Raymarch profile 70 32
  bash "$P" RingShower profile 70 32
  bash "$P" RingSpin profile 70 32
  bash "$P" Voronoi profile 70 32
  ;;
g3_ship)
  bash "$P" Liquid2D profile 70 16
  bash "$P" ShapeShifter profile 70 16
  bash "$P" MindSplatter profile 110 16
  bash "$P" HankinSolids profile 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"
  bash "$P" SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"
  ;;
g4_ship)
  bash "$P" Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  bash "$P" MeshFeedback profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  bash "$P" Flyby profile 310 16 "-D HS_PROFILE_EPOCH_REVS=2560"
  bash "$P" IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  ;;
g1_o3)
  bash "$P" BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  bash "$P" ChaoticStrings profile_o3 70 32
  bash "$P" DisplacementField profile_o3 70 32
  bash "$P" GnomonicStars profile_o3 70 32
  bash "$P" GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"
  ;;
g2_o3)
  bash "$P" HopfFibration profile_o3 70 32
  bash "$P" MobiusGrid profile_o3 70 32
  bash "$P" PetalFlow profile_o3 70 32
  bash "$P" Raymarch profile_o3 70 32
  bash "$P" RingShower profile_o3 70 32
  bash "$P" RingSpin profile_o3 70 32
  bash "$P" Voronoi profile_o3 70 32
  ;;
g3_o3)
  bash "$P" Liquid2D profile_o3 70 16
  bash "$P" ShapeShifter profile_o3 70 16
  bash "$P" MindSplatter profile_o3 110 16
  bash "$P" HankinSolids profile_o3 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"
  bash "$P" SphericalHarmonics profile_o3 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"
  ;;
g4_o3)
  bash "$P" Comets profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  bash "$P" MeshFeedback profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  bash "$P" Flyby profile_o3 310 16 "-D HS_PROFILE_EPOCH_REVS=2560"
  bash "$P" IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  bash "$P" DreamBalls profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  ;;
*) echo "unknown group $1"; exit 1;;
esac
echo "GROUP $1 DONE"
