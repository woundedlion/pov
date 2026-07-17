#!/bin/bash
# profile_one.sh <Effect> <env:profile|profile_o3> <seconds> <window> [extra flags...]
# Builds+flashes the profile image for one effect and captures its serial dump
# to build/prof/<effect>_<tag>.log. Verifies the capture header, and for a
# cycling effect that a preset marker appears (guards against a stale build
# silently flashing old code — see the validation notes in the skill). On a
# marker/header mismatch it wipes the env build dir and retries once.
#
# Takes the host-global device lock (tools/device_lock.sh) around the whole
# build+flash+capture, so concurrent agents queue instead of clobbering each
# other. HS_DEVICE_WAIT=<s> to queue rather than fail fast when it is busy.
set -e
. "$(dirname "$0")/device_lock.sh"
EFFECT=$1; ENV=$2; SECONDS_ARG=$3; WINDOW=$4; shift 4
EXTRA="$*"
TAG=$([ "$ENV" = "profile_o3" ] && echo o3 || echo ship)
LOWER=$(echo "$EFFECT" | tr '[:upper:]' '[:lower:]')
OUT=build/prof/${LOWER}_${TAG}.log
cd /c/work/Holosphere
mkdir -p build/prof
export PLATFORMIO_BUILD_FLAGS="-D HS_PROFILE_TARGET=$EFFECT -D HS_PROFILE_WINDOW=$WINDOW $EXTRA"

# Cyclers emit a per-advance marker; a capture of one must contain it.
CYCLERS="Liquid2D ShapeShifter MindSplatter DreamBalls Comets Flyby MeshFeedback HankinSolids SphericalHarmonics IslamicStars"
MARKER=""
case " $CYCLERS " in *" $EFFECT "*)
  case "$EFFECT" in
    ShapeShifter|MeshFeedback) MARKER="Shape:";;
    SphericalHarmonics) MARKER="Mode:";;
    IslamicStars) MARKER="Spawning Shape:";;
    HankinSolids) MARKER="Loading shape:";;
    *) MARKER="Preset:";;
  esac ;;
esac

capture() {
  echo "=== $EFFECT [$ENV] window=$WINDOW seconds=$SECONDS_ARG extra='$EXTRA'"
  pio run -e "$ENV" -t upload 2>&1 | tail -2
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
  local bad ok
  bad=$(grep -c "=== profile " "$OUT" || true)
  ok=$(grep -c "=== profile $EFFECT " "$OUT" || true)
  [ "$bad" = "$ok" ] || { echo "WINDOW NAME MISMATCH ($ok/$bad) — contention?"; return 1; }
  if [ -n "$MARKER" ]; then
    grep -q "$MARKER" "$OUT" || { echo "NO '$MARKER' MARKER — stale build?"; return 1; }
  fi
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
