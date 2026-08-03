#!/bin/bash
# profile_sweep.sh <group: g1_ship..g4_ship | check>
# Phantasm-roster profiling sweep, one group per invocation (see the
# teensy-profile skill for the per-effect duration/knob rationale).
# Covers the 288x144 Phantasm playlist. Dynamo and Thrusters are Holosphere
# 96x20-only (HS_PHANTASM_EFFECT_LIST excludes them), so they are not profiled
# here.
# The per-effect duration/window/epoch knobs are hand-tuned, so the groups below
# cannot be generated from the roster -- but their union is cross-checked
# against HS_PHANTASM_EFFECT_LIST on every invocation, and `check` runs that
# alone (no board needed). A newly added Phantasm effect otherwise goes
# unprofiled with nothing reporting it.
# Sequential profile_one.sh calls per group. A failed capture is recorded and
# the group carries on: the board can drop USB mid-capture, and aborting the
# group turns one dropped effect into every later one missing (a single hang
# once cost five captures). The group exits non-zero listing what failed, so
# a re-run only needs those effects.
# A capture must fit inside one epoch: crossing a boundary re-inits the effect
# mid-run, and an init that overruns the K-revolution commit window traps the
# board. The RD sims want 130 s of data, past the 120 s default, so they stretch
# the epoch to 150 s rather than lose a regime to a shorter capture.
set -u
P="$(dirname "$0")/profile_one.sh"
EFFECTS_H="$(dirname "$0")/../core/engine/effects.h"
FAILED=()
GROUP=${1:-}

# The firmware playlist, in HS_PHANTASM_EFFECT_LIST order. The macro body runs
# from the #define to the first line without a trailing backslash.
playlist_roster() {
  tr -d '\r' < "$EFFECTS_H" \
    | sed -n '/^#define HS_PHANTASM_EFFECT_LIST(X)/,/[^\\]$/p' \
    | sed -nE 's/^[[:space:]]*X\(([A-Za-z0-9_]+)\).*/\1/p'
}

# Every effect this script profiles, read off its own `run` lines.
swept_roster() {
  tr -d '\r' < "$0" | sed -nE 's/^[[:space:]]+run[[:space:]]+([A-Za-z0-9_]+).*/\1/p'
}

check_roster() {
  local playlist swept missing extra
  playlist="$(playlist_roster | sort)"
  swept="$(swept_roster | sort -u)"
  if [ -z "$playlist" ]; then
    echo "profile_sweep: parsed no effects from $EFFECTS_H — the roster scan broke" >&2
    return 1
  fi
  missing="$(comm -23 <(printf '%s\n' "$playlist") <(printf '%s\n' "$swept"))"
  extra="$(comm -13 <(printf '%s\n' "$playlist") <(printf '%s\n' "$swept"))"
  if [ -n "$missing" ] || [ -n "$extra" ]; then
    [ -n "$missing" ] && echo "profile_sweep: in HS_PHANTASM_EFFECT_LIST but in no group: $missing" >&2
    [ -n "$extra" ] && echo "profile_sweep: swept but not in HS_PHANTASM_EFFECT_LIST: $extra" >&2
    echo "profile_sweep: add the effect to a group below (with its duration/window" >&2
    echo "and epoch knobs), or drop it if the playlist no longer carries it." >&2
    return 1
  fi
  echo "profile_sweep: $(printf '%s\n' "$playlist" | wc -l | tr -d ' ') effect(s) match HS_PHANTASM_EFFECT_LIST"
}

# run <Effect> <env> <seconds> <window> [extra flags]
run() {
  if ! bash "$P" "$@"; then
    echo ">>> FAILED: $1 [$2] — continuing with the rest of the group" >&2
    FAILED+=("$1/$2")
  fi
}

# Ahead of the dispatch: a divergence must surface before a 20-minute capture,
# not after it.
check_roster || exit 1

case "$GROUP" in
check)
  exit 0
  ;;
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
  run MeshFeedback profile 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"
  run Flyby profile 310 16 "-D HS_PROFILE_EPOCH_REVS=2560"
  run IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  run DreamBalls profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  ;;
*) echo "unknown group $GROUP"; exit 1;;
esac
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "GROUP $GROUP INCOMPLETE — ${#FAILED[@]} failed: ${FAILED[*]}"
  exit 1
fi
echo "GROUP $GROUP DONE"
