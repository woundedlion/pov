#!/bin/bash
# profile_one.sh <Effect> <env:profile|profile_o3> <seconds> <window> [extra flags...]
# Builds+flashes the profile image for one effect and captures its serial dump
# to build/prof/<effect>_<tag>.log. Verifies the capture header, and for a
# cycling effect that a preset marker appears (guards against a stale build
# silently flashing old code — see the validation notes in the skill). On a
# marker/header mismatch it wipes the env build dir and retries once.
#
# Takes a per-board device lock (tools/device_lock.sh) around the whole
# build+flash+capture, so concurrent agents run on different boards, and queue
# instead of clobbering when every board is busy. HS_DEVICE_WAIT=<s> to queue
# rather than fail fast.
#
# With several Teensys attached, hs_device_acquire enumerates them, claims the
# first free one, and exports HS_TEENSY_PORT for it, which pins both the flash
# and the capture to that board. Never leave the board to auto-search when more
# than one is attached: the teensy-gui loader refuses to choose ("Found 2 Teensy
# boards") yet still exits SUCCESS, leaving the previous image running, and the
# capture's first-VID-match can read the other board — a plausible log of the
# wrong firmware. HS_TEENSY_PORT=<COMn> set by hand pins the run to one board
# and skips the search (use it to profile a specific board).
#
# HS_PROFILE_DEEP=1 additionally enables the HS_PROFILE_DEEP sub-scopes (the
# per-pixel/per-cell/per-face counters in shared render code, off by default
# because they tax every effect's numbers). A deep capture writes its own
# _deep.log so it cannot overwrite the roster log the standard reports cite.
#
# HS_PROFILE_TREE=<path> builds that checkout instead of the main one, so a
# branch can be profiled before it lands. Everything this script writes is
# relative to the tree (build/prof, .pio), so a worktree keeps its own logs and
# object dir and cannot mix results with the main tree's. Without it the default
# below is built no matter which directory the script is invoked from — running
# it inside a worktree would otherwise profile master's code under the branch's
# name, and nothing in the capture would say so.
set -eo pipefail
. "$(dirname "$0")/device_lock.sh"
EFFECT=$1; ENV=$2; SECONDS_ARG=$3; WINDOW=$4; shift 4
EXTRA="$*"
case "$ENV" in
  profile) TAG=ship;;
  profile_o3) TAG=o3;;
  *) echo "unknown profile environment: $ENV" >&2; exit 1;;
esac
LOWER=$(echo "$EFFECT" | tr '[:upper:]' '[:lower:]')
DEEP=""
DEEP_SUFFIX=""
if [ -n "$HS_PROFILE_DEEP" ] && [ "$HS_PROFILE_DEEP" != "0" ]; then
  DEEP="-D HS_PROFILE_DEEP_ENABLE"
  DEEP_SUFFIX="_deep"
fi
OUT=build/prof/${LOWER}_${TAG}${DEEP_SUFFIX}.log
cd "${HS_PROFILE_TREE:-/c/work/Holosphere}"
mkdir -p build/prof
export PLATFORMIO_BUILD_FLAGS="-D HS_PROFILE_TARGET=$EFFECT -D HS_PROFILE_WINDOW=$WINDOW $DEEP $EXTRA"

# Cyclers emit a per-advance marker; a capture of one must contain it.
CYCLERS="Liquid2D ShapeShifter MindSplatter DreamBalls Comets Flyby MeshFeedback HankinSolids SphericalHarmonics IslamicStars"
MARKER=""
case " $CYCLERS " in *" $EFFECT "*)
  case "$EFFECT" in
    ShapeShifter) MARKER="Shape:";;
    SphericalHarmonics) MARKER="Mode:";;
    IslamicStars) MARKER="Spawning Shape:";;
    HankinSolids) MARKER="Loading shape:";;
    *) MARKER="Preset:";;
  esac ;;
esac

TEENSY_TOOLS=${HS_TEENSY_TOOLS:-$HOME/.platformio/packages/tool-teensy}

# teensy_post_compile identifies a board by its USB location, not its COM name,
# and rejects a bare -portlabel. Resolve the location from the loader's own
# listing so a replug (which renumbers the hub port) needs no config change.
flash() {
  if [ -z "$HS_TEENSY_PORT" ]; then
    pio run -e "$ENV" -t upload 2>&1 | tail -2
    return
  fi
  pio run -e "$ENV" 2>&1 | tail -2
  local line loc label
  line=$("$TEENSY_TOOLS/teensy_ports.exe" -L | grep " $HS_TEENSY_PORT " || true)
  [ -n "$line" ] || { echo "no Teensy at $HS_TEENSY_PORT"; return 1; }
  loc=$(echo "$line" | awk '{print $1}')
  label=$(echo "$line" | awk '{print $2" "$3" "$4}')
  "$TEENSY_TOOLS/teensy_post_compile.exe" -file=firmware \
    -path="$(cygpath -w "$PWD/.pio/build/$ENV")" \
    -tools="$(cygpath -w "$TEENSY_TOOLS")" -board=TEENSY40 -reboot \
    "-port=$loc" "-portlabel=$label" -portprotocol=Teensy
}

capture() {
  echo "=== $EFFECT [$ENV] board=${HS_TEENSY_PORT:-auto} window=$WINDOW seconds=$SECONDS_ARG deep=${DEEP:-off} extra='$EXTRA'"
  flash
  # Let the capture's stderr through: it dies on a device trap (USB drops) and
  # on a port already held by a peer, and those look identical from the exit
  # code alone. Under set -e this aborts the run, so without the message the
  # failure reaches the caller with no reason attached.
  if ! python tools/profile_capture.py --seconds "$SECONDS_ARG" --out "$OUT" >/dev/null; then
    echo "CAPTURE FAILED (device trap, or port held by a peer?): $OUT" >&2
    return 1
  fi
}

verify() {
  grep -q "=== profile $EFFECT " "$OUT" || { echo "BAD/NO HEADER in $OUT"; return 1; }
  local bad ok configs matching_configs
  bad=$(grep -c "=== profile " "$OUT" || true)
  ok=$(grep -c "=== profile $EFFECT " "$OUT" || true)
  [ "$bad" = "$ok" ] || { echo "WINDOW NAME MISMATCH ($ok/$bad) — contention?"; return 1; }
  configs=$(grep -c "^profile harness:" "$OUT" || true)
  matching_configs=$(grep -c "^profile harness: effect=$EFFECT config=$TAG " "$OUT" || true)
  [ "$configs" -gt 0 ] || { echo "NO CONFIG TAG in $OUT"; return 1; }
  [ "$configs" = "$matching_configs" ] || {
    echo "CONFIG TAG MISMATCH ($matching_configs/$configs expected $TAG)"; return 1;
  }
  if [ -n "$MARKER" ]; then
    grep -q "$MARKER" "$OUT" || { echo "NO '$MARKER' MARKER — stale build?"; return 1; }
  fi
  # A flash that did not take leaves the previous image running, and the name
  # checks above pass whenever it happens to be the same effect (an -Os log
  # then publishes as the -O3 twin). The frame counter is the tell: a freshly
  # flashed board starts at 1, and the capture attaches within the 30 s
  # connect window, so a first frame past ~600 means no reboot happened.
  local first
  first=$(grep -m1 -oE "^f [0-9]+" "$OUT" | awk '{print $2}')
  [ -n "$first" ] || { echo "NO FRAME LINES in $OUT"; return 1; }
  [ "$first" -lt 600 ] || { echo "STALE IMAGE: first frame $first — the upload did not take"; return 1; }
  return 0
}

# ETA covers a clean rebuild + the capture + one retry: overshooting only
# delays a stale-break (safe), undershooting invites a peer to evict a live
# capture (not), and a crashed holder is reaped by the PID check regardless.
hs_device_acquire "$EFFECT" "$ENV" $((SECONDS_ARG * 2 + 900)) || exit 1
trap hs_device_release EXIT INT TERM

capture
if ! verify; then
  echo ">>> verify failed; wiping .pio/build/$ENV and retrying clean"
  rm -rf ".pio/build/$ENV"
  capture
  verify || { echo "FAILED after clean rebuild: $OUT"; exit 1; }
fi
NWIN=$(grep -c "=== profile $EFFECT " "$OUT")
NMARK=$([ -n "$MARKER" ] && grep -c "$MARKER" "$OUT" || echo "n/a")
echo "OK $OUT: $NWIN windows, markers($MARKER)=$NMARK"
