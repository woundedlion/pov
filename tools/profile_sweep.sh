#!/bin/bash
# profile_sweep.sh <group: g1_ship..g6_ship | check>
# Phantasm-playlist profiling sweep: one group of sequential profile_one.sh
# calls per invocation. A failed capture is recorded, the group carries on, and
# the run exits non-zero listing the failures. The groups' union is checked
# against HS_PHANTASM_EFFECT_LIST first; `check` runs only that (no board).
# A capture plus at least 20 s of attach slack must fit inside one epoch
# (HS_PROFILE_EPOCH_REVS): crossing a boundary re-inits the effect mid-run, and
# an init that overruns the commit window traps the board.
set -uo pipefail
P="$(dirname "$0")/profile_one.sh"
PLAYLIST_H="$(dirname "$0")/../targets/Phantasm/phantasm_playlist.h"
FAILED=()
GROUP=${1:-}

# The firmware playlist, in HS_PHANTASM_EFFECT_LIST order. X-row spelling
# mirrored by tools/docs_check.py.
playlist_roster() {
  tr -d '\r' < "$PLAYLIST_H" \
    | sed -n '/^#define HS_PHANTASM_EFFECT_LIST(X)/,/[^\\]$/p' \
    | sed -nE 's/^[[:space:]]*X\([[:space:]]*([A-Za-z0-9_]+)[[:space:]]*,.*/\1/p'
}

# Enumerate the reachable groups without flashing a board.
swept_roster() {
  local group
  for group in g{1..6}_ship; do
    HS_PROFILE_ROSTER_DUMP=1 bash "$0" "$group" || return 1
  done
}

check_roster() {
  local playlist swept missing extra
  if ! playlist="$(playlist_roster | sort)" || ! swept="$(swept_roster | sort)"; then
    echo "profile_sweep: cannot read and sort the effect rosters" >&2
    return 1
  fi
  if [ -z "$playlist" ]; then
    echo "profile_sweep: parsed no effects from $PLAYLIST_H — the roster scan broke" >&2
    return 1
  fi
  if ! missing="$(comm -23 <(printf '%s\n' "$playlist") <(printf '%s\n' "$swept"))" ||
      ! extra="$(comm -13 <(printf '%s\n' "$playlist") <(printf '%s\n' "$swept"))"; then
    echo "profile_sweep: cannot compare the effect rosters" >&2
    return 1
  fi
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
  if [ "${HS_PROFILE_ROSTER_DUMP:-0}" = 1 ]; then
    printf '%s\n' "$1"
    return 0
  fi
  if ! bash "$P" "$@"; then
    echo ">>> FAILED: $1 [$2] — continuing with the rest of the group" >&2
    FAILED+=("$1/$2")
  fi
}

# Before any capture, so a roster divergence fails fast.
if [ "${HS_PROFILE_ROSTER_DUMP:-0}" != 1 ]; then
  check_roster || exit 1
fi

case "$GROUP" in
check)
  exit 0
  ;;
g1_ship)
  run Galaxies profile 60 16
  run BZReactionDiffusion profile 60 32
  run Fishbowl profile 60 32
  run DisplacementField profile 121 32
  run GnomonicStars profile 30 32
  run GSReactionDiffusion profile 30 32
  ;;
g2_ship)
  run HopfFibration profile 60 32
  run MobiusGrid profile 105 16
  run PetalFlow profile 30 32
  run Raymarch profile 30 32
  run RingSpin profile 60 32
  run Voronoi profile 30 32
  ;;
g3_ship)
  run ShapeShifter profile 269 16
  run MindSplatter profile 80 16
  run HankinSolids profile 170 16
  run SphericalHarmonics profile 96 16
  ;;
g4_ship)
  run Comets profile 120 16
  run MeshFeedback profile 181 16
  run IslamicStars profile 306 16
  run DreamBalls profile 200 16
  ;;
g5_ship)
  run AlienBrain profile 240 16
  run KaleidoscopeHexSoft profile 60 32
  run AlienOcean profile 30 32
  run AlienCore profile 30 32
  run KaleidoscopeMandala profile 105 32
  run GridSpace profile 60 32
  run HyperLattice profile 300 16
  run ChromaticLichen profile 30 32
  run MermaidSkin profile 30 32
  run AshCloud profile 60 32
  ;;
g6_ship)
  run KaleidoscopePentBright profile 60 32
  run KaleidoscopeHexOil profile 105 16
  run KaleidoscopeSmooth profile 240 16
  run KaleidoscopeHexBright profile 105 32
  run KaleidoscopeFlowers profile 173 16
  run CosmicEyeball profile 60 32
  ;;
*) echo "unknown group $GROUP"; exit 1;;
esac
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "GROUP $GROUP INCOMPLETE — ${#FAILED[@]} failed: ${FAILED[*]}"
  exit 1
fi
[ "${HS_PROFILE_ROSTER_DUMP:-0}" = 1 ] || echo "GROUP $GROUP DONE"
