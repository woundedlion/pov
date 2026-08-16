#!/bin/bash
# profile_one.sh <Effect-or-ID> <env:profile|profile_o3> <seconds> <window> [extra flags...]
# Builds+flashes the profile image for one effect and captures its serial dump
# to build/prof/<effect>_<tag>.log. Verifies the capture header, and for a
# cycling effect that a preset marker appears (guards against a stale build
# silently flashing old code — see the validation notes in the skill). On a
# marker/header mismatch it wipes the env build dir and retries once.
#
# Host: Windows + Git Bash. The flash and the ELF attestation shell out to the
# PlatformIO loader's .exe tools through cygpath, and device_lock.sh enumerates
# boards by COM name; on any other host the run fails at the flash step and the
# per-board lock degrades to one bench-wide lock.
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
# object dir and cannot mix results with the main tree's. Without it the main
# checkout is built no matter which directory the script is invoked from —
# profiling a worktree's code under master's name (or the reverse) would leave
# nothing in the capture to say so.
# HS_PROFILE_MINDSPLATTER=counts|stalls builds a dedicated MindSplatter
# instrumentation image and writes a suffixed log. Count images also enable the
# generic Plot counters; neither image is valid for timing comparisons.
set -eo pipefail
. "$(dirname "$0")/device_lock.sh"
# Without this, a short/split argument list makes `shift 4` fail and set -e
# aborts with no message.
[ $# -ge 4 ] || {
  echo "usage: profile_one.sh <Effect-or-ID> <profile|profile_o3> <seconds> <window> [extra flags...]" >&2
  exit 1
}
EFFECT=$1; ENV=$2; SECONDS_ARG=$3; WINDOW=$4; shift 4
case "$EFFECT" in
  signal-weave) EFFECT=SignalWeave;;
  kaleido-wave) EFFECT=KaleidoWave;;
  alien-ocean) EFFECT=AlienOcean;;
  glitch-grid) EFFECT=GlitchGrid;;
  facet-wave) EFFECT=FacetWave;;
  contour-lattice) EFFECT=ContourLattice;;
  curl-lattice) EFFECT=CurlLattice;;
  prism-lattice) EFFECT=PrismLattice;;
  vector-facets) EFFECT=VectorFacets;;
  facet-grid) EFFECT=FacetGrid;;
  hex-wave) EFFECT=HexWave;;
  equator-grid) EFFECT=EquatorGrid;;
  cosmic-eyeball) EFFECT=CosmicEyeball;;
  mobius-grid) EFFECT=MobiusGrid;;
esac
for n in "$SECONDS_ARG" "$WINDOW"; do
  case $n in
    ''|*[!0-9]*)
      echo "seconds and window must be integers: '$SECONDS_ARG' '$WINDOW'" >&2
      exit 1;;
  esac
done
EXTRA="$*"
PROFILE_PRESET=$(printf '%s\n' "$EXTRA" |
  sed -nE 's/.*-D[[:space:]]*HS_PROFILE_PRESET=([0-9]+).*/\1/p')
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
REPLAY_SUFFIX=""
case " $EXTRA " in
  *" HS_MINDSPLATTER_REPLAY_AB "*|*"HS_MINDSPLATTER_REPLAY_AB="*)
    REPLAY_SUFFIX="_replay_ab";;
  *" HS_MINDSPLATTER_REPLAY "*|*"HS_MINDSPLATTER_REPLAY="*)
    REPLAY_SUFFIX="_replay";;
esac
MSP_FLAGS=""
MSP_SUFFIX=""
MSP_MARKER=""
case "${HS_PROFILE_MINDSPLATTER:-}" in
  "") ;;
  counts)
    [ "$EFFECT" = MindSplatter ] || {
      echo "HS_PROFILE_MINDSPLATTER requires effect MindSplatter" >&2; exit 1;
    }
    MSP_FLAGS="-D HS_PROFILE_MINDSPLATTER_COUNTS -D HS_PLOT_COUNTS"
    MSP_SUFFIX="_msp_counts"
    MSP_MARKER="msp counts particles:"
    ;;
  stalls)
    [ "$EFFECT" = MindSplatter ] || {
      echo "HS_PROFILE_MINDSPLATTER requires effect MindSplatter" >&2; exit 1;
    }
    MSP_FLAGS="-D HS_PROFILE_MINDSPLATTER_STALLS"
    MSP_SUFFIX="_msp_stalls"
    MSP_MARKER="msp stall: stage=history_vertex"
    ;;
  *)
    echo "HS_PROFILE_MINDSPLATTER must be counts or stalls" >&2
    exit 1
    ;;
esac
MODE_SUFFIX="${REPLAY_SUFFIX}${MSP_SUFFIX}"
OUT=${HS_PROFILE_OUT:-build/prof/${LOWER}_${TAG}${DEEP_SUFFIX}${MODE_SUFFIX}.log}
# Every per-run artifact hangs off the log's own stem, so HS_PROFILE_OUT moves
# the whole set. Spelling out the default naming here instead would let an
# overridden run (profile_islamic_big.sh) write its ELFs, maps, build logs and
# envdumps over the standard run's, and leave the standard run's .provenance
# naming the other configuration's ELF.
STEM=${OUT%.log}
PROVENANCE_OUT=$STEM.provenance
PROFILE_BUILD_LOG=${STEM}_build.log
PHANTASM_BUILD_LOG=${STEM}_phantasm_build.log
PROFILE_ENVDUMP=${STEM}_envdump.txt
PHANTASM_ENVDUMP=${STEM}_phantasm_envdump.txt
PROFILE_ELF=.pio/build/$ENV/firmware.elf
PROFILE_MAP=.pio/build/$ENV/firmware.map
ATTEST_DIR=$(dirname "$STEM")/attest/$(basename "$STEM")
PHANTASM_ELF=$ATTEST_DIR/phantasm.elf
PHANTASM_MAP=$ATTEST_DIR/phantasm.map
ARM_READELF=${HS_ARM_READELF:-$HOME/.platformio/packages/toolchain-gccarmnoneeabi-teensy/bin/arm-none-eabi-readelf.exe}
ATTESTED_COMPILER=

# The main checkout: git's common dir is the main repo's .git from every
# worktree, so its parent is that checkout wherever this script is invoked from.
main_tree() {
  local common
  common=$(git -C "$(dirname "$0")" rev-parse --path-format=absolute \
             --git-common-dir) || return 1
  dirname "$common"
}
TREE=${HS_PROFILE_TREE:-$(main_tree)} ||
  { echo "cannot resolve the main checkout from $(dirname "$0")" >&2; exit 1; }
cd "$TREE" || { echo "no such tree: $TREE" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
export PLATFORMIO_BUILD_FLAGS="-D HS_PROFILE_TARGET=$EFFECT -D HS_PROFILE_WINDOW=$WINDOW $DEEP $MSP_FLAGS $EXTRA"

# Cyclers emit a per-advance marker; a capture of one must contain it.
CYCLERS="SignalWeave CurlLattice FacetGrid EquatorGrid MobiusGrid ShapeShifter MindSplatter DreamBalls Comets MeshFeedback HankinSolids SphericalHarmonics IslamicStars"
MARKER=""
case " $CYCLERS " in *" $EFFECT "*)
  case "$EFFECT" in
    SphericalHarmonics) MARKER="Mode:";;
    IslamicStars) MARKER="Spawning Shape:";;
    HankinSolids) MARKER="Loading shape:";;
    *) MARKER="Preset:";;
  esac ;;
esac
if [ -n "$PROFILE_PRESET" ]; then
  MARKER="Profile preset:"
fi
if [ -n "$REPLAY_SUFFIX" ]; then
  MARKER="replay corpus:"
fi

TEENSY_TOOLS=${HS_TEENSY_TOOLS:-$HOME/.platformio/packages/tool-teensy}

build_image() {
  local env=$1 log=$2
  if ! pio run -e "$env" >"$log" 2>&1; then
    tail -30 "$log"
    return 1
  fi
  tail -2 "$log"
}

elf_compiler() {
  "$ARM_READELF" -p .comment "$1" 2>/dev/null |
    sed -n 's/^  \[[^]]*\]  //p' | sort -u
}

elf_abi() {
  "$ARM_READELF" -A "$1" 2>/dev/null |
    sed -n \
      -e '/Tag_CPU_name:/p' \
      -e '/Tag_CPU_arch:/p' \
      -e '/Tag_FP_arch:/p' \
      -e '/Tag_ABI_HardFP_use:/p' \
      -e '/Tag_ABI_VFP_args:/p' |
    sed 's/^[[:space:]]*//'
}

package_fingerprint() {
  sed -n \
    -e '/^PLATFORM:/p' \
    -e '/^ - framework-arduinoteensy /p' \
    -e '/^ - toolchain-gccarmnoneeabi-teensy /p' "$1"
}

file_sha256() {
  sha256sum "$1" | awk '{print $1}'
}

assert_matching_toolchains() {
  local profile_compiler phantasm_compiler profile_abi phantasm_abi
  local profile_packages phantasm_packages
  profile_compiler=$(elf_compiler "$PROFILE_ELF") ||
    { echo "PROFILE PROVENANCE FAIL: cannot inspect $PROFILE_ELF"; return 1; }
  phantasm_compiler=$(elf_compiler "$PHANTASM_ELF") ||
    { echo "PROFILE PROVENANCE FAIL: cannot inspect $PHANTASM_ELF"; return 1; }
  [ -n "$profile_compiler" ] ||
    { echo "PROFILE PROVENANCE FAIL: no compiler in $PROFILE_ELF"; return 1; }
  [ "$profile_compiler" = "$phantasm_compiler" ] || {
    echo "PROFILE PROVENANCE FAIL: compiler mismatch"
    echo "  profile:  $profile_compiler"
    echo "  phantasm: $phantasm_compiler"
    return 1
  }

  profile_abi=$(elf_abi "$PROFILE_ELF") ||
    { echo "PROFILE PROVENANCE FAIL: cannot read ARM ABI from $PROFILE_ELF"; return 1; }
  phantasm_abi=$(elf_abi "$PHANTASM_ELF") ||
    { echo "PROFILE PROVENANCE FAIL: cannot read ARM ABI from $PHANTASM_ELF"; return 1; }
  [ -n "$profile_abi" ] ||
    { echo "PROFILE PROVENANCE FAIL: no ARM ABI in $PROFILE_ELF"; return 1; }
  [ "$profile_abi" = "$phantasm_abi" ] ||
    { echo "PROFILE PROVENANCE FAIL: ARM ABI mismatch"; return 1; }

  profile_packages=$(package_fingerprint "$PROFILE_BUILD_LOG")
  phantasm_packages=$(package_fingerprint "$PHANTASM_BUILD_LOG")
  [ -n "$profile_packages" ] ||
    { echo "PROFILE PROVENANCE FAIL: no package fingerprint"; return 1; }
  [ "$profile_packages" = "$phantasm_packages" ] || {
    echo "PROFILE PROVENANCE FAIL: PlatformIO package mismatch"
    diff -u <(echo "$phantasm_packages") <(echo "$profile_packages") || true
    return 1
  }
  ATTESTED_COMPILER=$profile_compiler
}

build_and_attest() {
  mkdir -p "$ATTEST_DIR"
  PLATFORMIO_BUILD_FLAGS='' pio run -e phantasm -t envdump >"$PHANTASM_ENVDUMP"
  PLATFORMIO_BUILD_FLAGS='' build_image phantasm "$PHANTASM_BUILD_LOG"
  cp .pio/build/phantasm/firmware.elf "$PHANTASM_ELF"
  cp .pio/build/phantasm/firmware.map "$PHANTASM_MAP"
  pio run -e "$ENV" -t envdump >"$PROFILE_ENVDUMP"
  build_image "$ENV" "$PROFILE_BUILD_LOG"
  assert_matching_toolchains

  local source_sha profile_sha phantasm_sha artifact_dir
  local profile_flags_sha phantasm_flags_sha
  source_sha=$(git rev-parse HEAD)
  profile_sha=$(file_sha256 "$PROFILE_ELF")
  phantasm_sha=$(file_sha256 "$PHANTASM_ELF")
  profile_flags_sha=$(file_sha256 "$PROFILE_ENVDUMP")
  phantasm_flags_sha=$(file_sha256 "$PHANTASM_ENVDUMP")
  artifact_dir=build/prof/artifacts/${LOWER}_${TAG}_${source_sha:0:12}_${profile_sha:0:12}
  mkdir -p "$artifact_dir"
  cp "$PROFILE_ELF" "$artifact_dir/profile.elf"
  cp "$PROFILE_MAP" "$artifact_dir/profile.map"
  cp "$PHANTASM_ELF" "$artifact_dir/phantasm.elf"
  cp "$PHANTASM_MAP" "$artifact_dir/phantasm.map"
  cp "$PROFILE_BUILD_LOG" "$artifact_dir/profile_build.log"
  cp "$PHANTASM_BUILD_LOG" "$artifact_dir/phantasm_build.log"
  cp "$PROFILE_ENVDUMP" "$artifact_dir/profile_envdump.txt"
  cp "$PHANTASM_ENVDUMP" "$artifact_dir/phantasm_envdump.txt"
  git status --porcelain >"$artifact_dir/source_status.txt"
  git diff --binary >"$artifact_dir/source.diff"

  {
    echo "profile provenance: version=1"
    echo "profile provenance: source_sha=$source_sha"
    echo "profile provenance: compiler=$ATTESTED_COMPILER"
    echo "profile provenance: profile_elf_sha256=$profile_sha"
    echo "profile provenance: phantasm_elf_sha256=$phantasm_sha"
    echo "profile provenance: profile_envdump_sha256=$profile_flags_sha"
    echo "profile provenance: phantasm_envdump_sha256=$phantasm_flags_sha"
    echo "profile provenance: artifact_profile_elf=$artifact_dir/profile.elf"
    echo "profile provenance: artifact_phantasm_elf=$artifact_dir/phantasm.elf"
    echo "profile provenance: artifact_dir=$artifact_dir"
  } >"$PROVENANCE_OUT"
}

# teensy_post_compile identifies a board by its USB location, not its COM name.
flash() {
  [ -n "$HS_TEENSY_PORT" ] ||
    { echo "device lock did not pin a Teensy port"; return 1; }
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
  build_and_attest
  flash
  # Let the capture's stderr through: it dies on a device trap (USB drops) and
  # on a port already held by a peer, and those look identical from the exit
  # code alone. Under set -e this aborts the run, so without the message the
  # failure reaches the caller with no reason attached.
  if ! python tools/profile_capture.py --seconds "$SECONDS_ARG" --out "$OUT" >/dev/null; then
    echo "CAPTURE FAILED (device trap, or port held by a peer?): $OUT" >&2
    return 1
  fi
  cat "$PROVENANCE_OUT" >>"$OUT"
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
  local provenance profile_sha artifact_profile_elf actual_profile_sha
  provenance=$(grep -c "^profile provenance: version=1$" "$OUT" || true)
  [ "$provenance" = 1 ] ||
    { echo "NO/INVALID PROFILE PROVENANCE ($provenance)"; return 1; }
  profile_sha=$(sed -n 's/^profile provenance: profile_elf_sha256=//p' "$OUT")
  artifact_profile_elf=$(sed -n 's/^profile provenance: artifact_profile_elf=//p' "$OUT")
  [ -n "$profile_sha" ] && [ -f "$artifact_profile_elf" ] ||
    { echo "MISSING PROFILE ARTIFACT"; return 1; }
  actual_profile_sha=$(file_sha256 "$artifact_profile_elf")
  [ "$profile_sha" = "$actual_profile_sha" ] ||
    { echo "PROFILE ARTIFACT HASH MISMATCH"; return 1; }
  if [ -n "$PROFILE_PRESET" ]; then
    grep -q "^Profile preset: $PROFILE_PRESET/" "$OUT" || {
      echo "PROFILE PRESET MISMATCH (expected $PROFILE_PRESET)"; return 1;
    }
  fi
  if [ -n "$MARKER" ]; then
    grep -q "$MARKER" "$OUT" || { echo "NO '$MARKER' MARKER — stale build?"; return 1; }
  fi
  if [ -n "$MSP_MARKER" ]; then
    grep -q "$MSP_MARKER" "$OUT" || {
      echo "NO '$MSP_MARKER' INSTRUMENTATION; stale build?"; return 1;
    }
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
