#!/bin/bash
# profile_sweep.sh <group: g1_ship..g6_ship | check>
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
# board. The profile image's epoch is one hour (targets/Profile/Profile.ino);
# HS_PROFILE_EPOCH_REVS replaces it (hardware/pov_segmented.h), so an override
# only ever shortens it. Attach latency must fit between the capture duration
# and the epoch boundary (at least 20 s of slack for the entries below).
# parse_profile validate rejects captures that cross an epoch reset.
set -uo pipefail
P="$(dirname "$0")/profile_one.sh"
PLAYLIST_H="$(dirname "$0")/../targets/Phantasm/phantasm_playlist.h"
FAILED=()
GROUP=${1:-}

# The firmware playlist, in HS_PHANTASM_EFFECT_LIST order. The macro body runs
# from the #define to the first line without a trailing backslash.
# Shared X-row spelling with tools/docs_check.py and scripts/effect_roster.mjs:
# whitespace inside the parens is tolerated so a reformat to `X( Foo , 12 )`
# cannot silently shrink one roster reader's count. Anchoring at the line start
# is what excludes a commented-out row.
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

# Ahead of the dispatch: a divergence must surface before a 20-minute capture,
# not after it.
if [ "${HS_PROFILE_ROSTER_DUMP:-0}" != 1 ]; then
  check_roster || exit 1
fi

case "$GROUP" in
check)
  exit 0
  ;;
g1_ship)
  run BZReactionDiffusion profile 130 32
  run Fishbowl profile 70 32
  run DisplacementField profile 70 32
  run GnomonicStars profile 70 32
  run GSReactionDiffusion profile 130 32
  ;;
g2_ship)
  run HopfFibration profile 70 32
  run MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  run PetalFlow profile 70 32
  run Raymarch profile 70 32
  run RingShower profile 70 32
  run RingSpin profile 70 32
  run Voronoi profile 70 32
  ;;
g3_ship)
  run ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  run MindSplatter profile 110 16
  run HankinSolids profile 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"
  run SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"
  ;;
g4_ship)
  run Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run MeshFeedback profile 420 16
  run IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"
  run DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"
  ;;
g5_ship)
  run AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"
  run KaleidoscopeHexSoft profile 70 32
  run AlienOcean profile 70 32
  run AlienCore profile 70 32
  run KaleidoscopeMandala profile 150 32
  run GridSpace profile 70 32
  run HyperLattice profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  run LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"
  run ChromaticLichen profile 70 32
  run MermaidSkin profile 70 32
  run AshCloud profile 70 32
  ;;
g6_ship)
  run KaleidoscopePentBright profile 70 32
  run KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"
  run KaleidoscopeStainedGlass profile 70 32
  run KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run KaleidoscopeHexBright profile 150 32
  run KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"
  run CosmicEyeball profile 70 32
  ;;
*) echo "unknown group $GROUP"; exit 1;;
esac
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "GROUP $GROUP INCOMPLETE — ${#FAILED[@]} failed: ${FAILED[*]}"
  exit 1
fi
[ "${HS_PROFILE_ROSTER_DUMP:-0}" = 1 ] || echo "GROUP $GROUP DONE"
