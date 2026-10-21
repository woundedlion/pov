#!/bin/bash
# profile_one.sh [--dry-run|--build-only] [--elf <entry>/firmware.elf]
#                <Effect-or-ID> <env:profile|profile_o3> <seconds> <window> [extra flags...]
# Builds+flashes the profile image for one effect and captures its serial dump
# to a log under build/prof/ by default, then verifies the capture header and, for a
# cycling effect, a preset marker; on a mismatch it rebuilds the profile image
# clean and retries once. The shipping phantasm image must pass its size/layout
# gates and supplies the ELF for compiler/ABI attestation. Each image is reused
# from tools/fw_cache.sh when a committed checkout already built it with the
# same flags and toolchain, and is stored there after a build. Host: Windows +
# Git Bash.
#
# Options:
#   --dry-run     print each image's cache key and hit or miss; no lock, build or board
#   --build-only  build, attest and cache both images under the checkout lock; no board
#   --elf <path>  flash this cached profile image; refused unless its entry's key
#                 matches this checkout, environment and flags
#
# Env knobs:
#   HS_DEVICE_WAIT=<s>      queue for a free board, and let the pinned board drop
#                           off the port list for up to HS_DEVICE_ABSENT_GRACE
#                           (tools/device_lock.sh)
#   HS_TEENSY_PORT=<COMn>   pin one board
#   HS_ARM_READELF          readelf for ELF attestation
#   HS_PROFILE_DEEP=1       enable the HS_PROFILE_DEEP sub-scopes; writes _deep.log
#   HS_PIO_JOBS=<n>         `pio run -j` for the builds (tools/fw_cache.sh)
#   HS_PROFILE_TREE=<path>  checkout to build (default: the one holding this script)
#   HS_PROFILE_OUT=<path>   capture log path; all per-run artifacts share its stem
#   HS_PROFILE_EXPECT_SHAPE=<name> require "Spawning Shape: <name>" in the capture
#                           (used by profile_islamic_big.sh for IslamicStars)
#   HS_PROFILE_MINDSPLATTER=counts|stalls  MindSplatter instrumentation image with
#                           a suffixed log; not valid for timing comparisons
#   HS_PROFILE_PC=1         PC-sampler image (-D HS_PROFILE_PC) with a _pc log;
#                           dumps PC_DUMP_MARGIN_S before the capture ends and
#                           writes tools/pc_profile.py's report to <stem>_report.txt
set -euo pipefail
# shellcheck source-path=SCRIPTDIR source=device_lock.sh
. "$(dirname "$0")/device_lock.sh"
# shellcheck source-path=SCRIPTDIR source=teensy_flash.sh
. "$(dirname "$0")/teensy_flash.sh"
# shellcheck source-path=SCRIPTDIR source=fw_cache.sh
. "$(dirname "$0")/fw_cache.sh"
usage() {
  echo "usage: profile_one.sh [--dry-run|--build-only] [--elf <path>] <Effect-or-ID> <profile|profile_o3> <seconds> <window> [extra flags...]" >&2
  exit 1
}
MODE=capture
PREBUILT=""
while [ $# -gt 0 ]; do
  case $1 in
    --dry-run) MODE=dry-run; shift;;
    --build-only) MODE=build-only; shift;;
    --elf)
      [ $# -ge 2 ] || usage
      [ -f "$2" ] || { echo "no such ELF: $2" >&2; exit 1; }
      PREBUILT=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
      shift 2;;
    *) break;;
  esac
done
# Under set -e, `shift 4` on a short list would abort with no message.
[ $# -ge 4 ] || usage
EFFECT=$1; ENV=$2; SECONDS_ARG=$3; WINDOW=$4; shift 4
# The checkout containing the invoked script, including a linked worktree.
profile_tree() {
  git -C "$(dirname "$0")" rev-parse --path-format=absolute --show-toplevel
}
TREE=${HS_PROFILE_TREE:-$(profile_tree)} ||
  { echo "cannot resolve the profile checkout from $(dirname "$0")" >&2; exit 1; }

effect_class() {
  awk -v wanted="$1" '
    /^[[:space:]]*X\(/ {
      name=$0; sub(/^[[:space:]]*X\(/, "", name); sub(/,.*/, "", name)
      id=$0; sub(/^[^"]*"/, "", id); sub(/".*/, "", id)
      if (wanted == id || wanted == name) { print name; found=1; exit }
    }
    END { if (!found) exit 1 }
  ' "$TREE/core/platform/effect_identity.h"
}
REQUESTED_EFFECT=$EFFECT
EFFECT=$(effect_class "$REQUESTED_EFFECT") || {
  echo "unknown profile effect: $REQUESTED_EFFECT" >&2; exit 1;
}
for n in "$SECONDS_ARG" "$WINDOW"; do
  case $n in
    ''|*[!0-9]*)
      echo "seconds and window must be integers: '$SECONDS_ARG' '$WINDOW'" >&2
      exit 1;;
  esac
done
EXTRA="$*"
# Seconds the PC dump leaves between itself and the end of the capture.
PC_DUMP_MARGIN_S=8
# Least seconds the capture waits for the flashed board's serial port.
CAPTURE_CONNECT_S=30
PC_FLAGS=""
PC_SUFFIX=""
PC_MARKER=""
if [ -n "${HS_PROFILE_PC:-}" ] && [ "${HS_PROFILE_PC:-}" != "0" ]; then
  [ "$SECONDS_ARG" -gt $((PC_DUMP_MARGIN_S + 2)) ] || {
    echo "HS_PROFILE_PC needs more than $((PC_DUMP_MARGIN_S + 2)) seconds" >&2
    exit 1
  }
  PC_FLAGS="-D HS_PROFILE_PC -D HS_PROFILE_PC_MS=$(((SECONDS_ARG - PC_DUMP_MARGIN_S) * 1000))"
  PC_SUFFIX="_pc"
  PC_MARKER="pc profile end"
fi
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
if [ -n "${HS_PROFILE_DEEP:-}" ] && [ "${HS_PROFILE_DEEP:-}" != "0" ]; then
  DEEP="-D HS_PROFILE_DEEP_ENABLE"
  DEEP_SUFFIX="_deep"
fi
REPLAY_SUFFIX=$(printf '%s\n' "$EXTRA" |
  sed -nE 's/.*-D[[:space:]]*HS_MINDSPLATTER_REPLAY(_AB)?([=[:space:]]|$).*/_replay\1/p' |
  tr '[:upper:]' '[:lower:]')
MSP_FLAGS=""
MSP_SUFFIX=""
MSP_MARKER=""
case "${HS_PROFILE_MINDSPLATTER:-}" in
  "") ;;
  counts)
    [ "$EFFECT" = MindSplatter ] || {
      echo "HS_PROFILE_MINDSPLATTER requires effect MindSplatter" >&2; exit 1;
    }
    MSP_FLAGS="-D HS_PROFILE_PLOT_RENDER_COUNTS -D HS_PLOT_COUNTS"
    MSP_SUFFIX="_msp_counts"
    MSP_MARKER="plot render counts particles:"
    ;;
  stalls)
    [ "$EFFECT" = MindSplatter ] || {
      echo "HS_PROFILE_MINDSPLATTER requires effect MindSplatter" >&2; exit 1;
    }
    MSP_FLAGS="-D HS_PROFILE_PLOT_STALLS"
    MSP_SUFFIX="_msp_stalls"
    MSP_MARKER="plot stall: stage=history_vertex"
    ;;
  *)
    echo "HS_PROFILE_MINDSPLATTER must be counts or stalls" >&2
    exit 1
    ;;
esac
MODE_SUFFIX="${REPLAY_SUFFIX}${MSP_SUFFIX}${PC_SUFFIX}"
OUT=$(hs_native_path "${HS_PROFILE_OUT:-build/prof/${LOWER}_${TAG}${DEEP_SUFFIX}${MODE_SUFFIX}.log}")
# Every per-run artifact shares the log's stem, and the provenance lines name
# them in a form the Windows Python report tools open.
STEM=${OUT%.log}
PROVENANCE_OUT=$STEM.provenance
PROFILE_BUILD_LOG=${STEM}_build.log
PHANTASM_BUILD_LOG=${STEM}_phantasm_build.log
PROFILE_ENVDUMP=${STEM}_envdump.txt
PHANTASM_ENVDUMP=${STEM}_phantasm_envdump.txt
PROFILE_DIR=.pio/build/$ENV
PROFILE_ELF=$PROFILE_DIR/firmware.elf
PROFILE_MAP=$PROFILE_DIR/firmware.map
ARTIFACT_BASE=$(dirname "$STEM")/artifacts/$(basename "$STEM")
ATTEST_DIR=$(dirname "$STEM")/attest/$(basename "$STEM")
PHANTASM_ELF=$ATTEST_DIR/phantasm.elf
PHANTASM_MAP=$ATTEST_DIR/phantasm.map
ARM_READELF=${HS_ARM_READELF:-${PLATFORMIO_CORE_DIR:-$HOME/.platformio}/packages/toolchain-gccarmnoneeabi-teensy/bin/arm-none-eabi-readelf.exe}
ATTESTED_COMPILER=

cd "$TREE" || { echo "no such tree: $TREE" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"
export PLATFORMIO_BUILD_FLAGS="-D HS_PROFILE_TARGET=$EFFECT -D HS_PROFILE_WINDOW=$WINDOW $DEEP $MSP_FLAGS $PC_FLAGS $EXTRA"

# Cyclers emit a per-advance marker; a capture of one must contain it.
CYCLERS="AlienBrain KaleidoscopeSmooth KaleidoscopeFlowers KaleidoscopeMandala KaleidoscopeHexOil KaleidoscopeHexBright HyperLattice MobiusGrid ShapeShifter MindSplatter DreamBalls Comets MeshFeedback HankinSolids SphericalHarmonics IslamicStars DisplacementField"
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


build_image() {
  local env=$1 log=$2
  if ! hs_fw_pio "$env" >"$log" 2>&1; then
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

envdump() {
  hs_fw_pio "$1" -t envdump | "$_HS_LOCK_PYTHON" tools/profile_envdump.py >"$2"
}

# obtain_image <env> <envdump> <build-log>: reuses the cached image, or builds
# and caches it; sets IMAGE_DIR and IMAGE_SOURCE.
obtain_image() {
  local env=$1 dump=$2 log=$3 rc=0
  if [ "$env" = "$ENV" ] && [ -n "$FW_REBUILD" ]; then
    hs_fw_cache_key "$env" || rc=$?
  elif hs_fw_cache_lookup "$env"; then
    cp "$HS_FW_ENTRY/build.log" "$log"
    if [ -f "$HS_FW_ENTRY/envdump.txt" ]; then
      cp "$HS_FW_ENTRY/envdump.txt" "$dump"
    else
      envdump "$env" "$dump"
    fi
    IMAGE_DIR=$HS_FW_ENTRY
    IMAGE_SOURCE=cache:$HS_FW_KEY
    echo "fw_cache: $env reused $HS_FW_ENTRY"
    return 0
  else
    rc=$?
  fi
  [ "$rc" != 2 ] || return 1
  envdump "$env" "$dump"
  build_image "$env" "$log"
  IMAGE_DIR=$(hs_fw_build_dir)/$env
  IMAGE_SOURCE=built
  if [ "$rc" != 3 ] && hs_fw_cache_store "$env" "$IMAGE_DIR" "$log" "$dump" "$FW_REBUILD"; then
    IMAGE_DIR=$HS_FW_ENTRY
    IMAGE_SOURCE=built:$HS_FW_KEY
  fi
}

# Adopts the --elf cache entry; sets IMAGE_DIR and IMAGE_SOURCE.
prebuilt_image() {
  local entry
  entry=$(dirname "$PREBUILT")
  if [ "$(basename "$PREBUILT")" != firmware.elf ] || ! hs_fw_cache_key "$ENV" ||
     ! hs_fw_cache_check_entry "$entry" "$HS_FW_KEY"; then
    echo "PROFILE PROVENANCE FAIL: $PREBUILT is not a cached $ENV image of this committed checkout with flags '$PLATFORMIO_BUILD_FLAGS'"
    return 1
  fi
  cp "$entry/build.log" "$PROFILE_BUILD_LOG"
  if [ -f "$entry/envdump.txt" ]; then
    cp "$entry/envdump.txt" "$PROFILE_ENVDUMP"
  else
    envdump "$ENV" "$PROFILE_ENVDUMP"
  fi
  IMAGE_DIR=$entry
  IMAGE_SOURCE=prebuilt:$HS_FW_KEY
}

build_and_attest() {
  local phantasm_source
  mkdir -p "$ATTEST_DIR"
  _hs_resolve_python
  PLATFORMIO_BUILD_FLAGS='' obtain_image phantasm "$PHANTASM_ENVDUMP" "$PHANTASM_BUILD_LOG"
  phantasm_source=$IMAGE_SOURCE
  cp "$IMAGE_DIR/firmware.elf" "$PHANTASM_ELF"
  cp "$IMAGE_DIR/firmware.map" "$PHANTASM_MAP"
  if [ -n "$PREBUILT" ]; then
    prebuilt_image
  else
    obtain_image "$ENV" "$PROFILE_ENVDUMP" "$PROFILE_BUILD_LOG"
  fi
  PROFILE_DIR=$IMAGE_DIR
  PROFILE_ELF=$PROFILE_DIR/firmware.elf
  PROFILE_MAP=$PROFILE_DIR/firmware.map
  assert_matching_toolchains

  local source_sha profile_sha phantasm_sha artifact_dir
  local profile_flags_sha phantasm_flags_sha
  source_sha=$(git rev-parse HEAD)
  profile_sha=$(file_sha256 "$PROFILE_ELF")
  phantasm_sha=$(file_sha256 "$PHANTASM_ELF")
  profile_flags_sha=$(file_sha256 "$PROFILE_ENVDUMP")
  phantasm_flags_sha=$(file_sha256 "$PHANTASM_ENVDUMP")
  artifact_dir=${ARTIFACT_BASE}_${source_sha:0:12}_${profile_sha:0:12}
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
  git diff --binary HEAD >"$artifact_dir/source.diff"

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
    echo "profile provenance: profile_image=$IMAGE_SOURCE"
    echo "profile provenance: phantasm_image=$phantasm_source"
  } >"$PROVENANCE_OUT"
}


capture() {
  echo "=== $EFFECT [$ENV] board=${HS_TEENSY_PORT:-auto} window=$WINDOW seconds=$SECONDS_ARG deep=${DEEP:-off} extra='$EXTRA'"
  build_and_attest
  hs_teensy_flash "$PROFILE_DIR"
  local connect
  connect=$(hs_device_absent_budget)
  [ "$connect" -ge "$CAPTURE_CONNECT_S" ] || connect=$CAPTURE_CONNECT_S
  # stderr passes through: a device trap and a peer-held port share an exit code.
  if ! "$_HS_LOCK_PYTHON" tools/profile_capture.py --seconds "$SECONDS_ARG" \
      --connect-timeout "$connect" --out "$OUT" >/dev/null; then
    echo "CAPTURE FAILED (device trap, or port held by a peer?): $OUT" >&2
    return 1
  fi
  cat "$PROVENANCE_OUT" >>"$OUT"
}

RETRY_CACHE=""
FW_REBUILD=""

# A --elf image is reflashed as is; otherwise the profile image is rebuilt
# clean and replaces its cache entry.
prepare_retry() {
  [ -z "$PREBUILT" ] || return 0
  rm -rf "$(hs_fw_build_dir)/${ENV:?}"
  RETRY_CACHE=$(mktemp -d "${TMPDIR:-/tmp}/hs-profile-cache.XXXXXX")
  export PLATFORMIO_BUILD_CACHE_DIR="$RETRY_CACHE"
  FW_REBUILD=1
}

cleanup() {
  hs_device_release
  hs_tree_release
  [ -z "$RETRY_CACHE" ] || rm -rf "$RETRY_CACHE"
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
  if [ -n "${HS_PROFILE_EXPECT_SHAPE:-}" ]; then
    grep -Fq "Spawning Shape: $HS_PROFILE_EXPECT_SHAPE" "$OUT" || {
      echo "PROFILE SHAPE MISMATCH (expected $HS_PROFILE_EXPECT_SHAPE)"; return 1;
    }
  fi
  if [ -n "$MSP_MARKER" ]; then
    grep -q "$MSP_MARKER" "$OUT" || {
      echo "NO '$MSP_MARKER' INSTRUMENTATION; stale build?"; return 1;
    }
  fi
  if [ -n "${PC_MARKER:-}" ]; then
    grep -qx "$PC_MARKER" "$OUT" || {
      echo "NO '$PC_MARKER' LINE; PC dump missing or cut short"; return 1;
    }
  fi
  # A flash that did not take leaves the previous image running; a freshly
  # flashed board starts at frame 1.
  local first
  first=$(grep -m1 -oE "^f [0-9]+" "$OUT" | awk '{print $2}')
  [ -n "$first" ] || { echo "NO FRAME LINES in $OUT"; return 1; }
  [ "$first" -lt 600 ] || { echo "STALE IMAGE: first frame $first — the upload did not take"; return 1; }
  return 0
}

# Each image's key and cache state; touches no lock, build or board.
dry_run() {
  local env state
  for env in phantasm "$ENV"; do
    state=0
    if [ "$env" = phantasm ]; then
      PLATFORMIO_BUILD_FLAGS='' hs_fw_cache_lookup "$env" || state=$?
    else
      hs_fw_cache_lookup "$env" || state=$?
    fi
    case $state in
      0) echo "dry-run: $env hit $HS_FW_KEY $HS_FW_ENTRY";;
      1) echo "dry-run: $env miss $HS_FW_KEY";;
      3) echo "dry-run: $env uncached (uncommitted changes)";;
      *) return 1;;
    esac
  done
  if [ -n "$PREBUILT" ]; then
    if [ "$(basename "$PREBUILT")" = firmware.elf ] && hs_fw_cache_key "$ENV" &&
       hs_fw_cache_check_entry "$(dirname "$PREBUILT")" "$HS_FW_KEY"; then
      echo "dry-run: --elf accepted $PREBUILT"
    else
      echo "dry-run: --elf refused $PREBUILT"
      return 1
    fi
  fi
  echo "dry-run: flags=$PLATFORMIO_BUILD_FLAGS"
}

if [ "$MODE" = dry-run ]; then
  dry_run
  exit
fi

# Armed before the claim (release is a no-op without a token); INT/TERM exit
# and leave the release to the EXIT trap.
hs_exit_traps cleanup

# ETA covers a clean rebuild, capture, and one retry.
acquire_tree_lock "$((SECONDS_ARG * 2 + 900))" || exit $?

if [ "$MODE" = build-only ]; then
  build_and_attest
  sed -n 's/^profile provenance: \(artifact_profile_elf\|profile_image\)=/\1 /p' "$PROVENANCE_OUT"
  exit 0
fi

hs_device_acquire "$EFFECT" "$ENV" $((SECONDS_ARG * 2 + 900)) || exit $?

capture
if ! verify; then
  echo ">>> verify failed; wiping the $ENV build and retrying with an isolated cache"
  prepare_retry
  capture
  verify || { echo "FAILED after clean rebuild: $OUT"; exit 1; }
fi
NWIN=$(grep -c "=== profile $EFFECT " "$OUT")
# grep -c prints 0 and exits 1 on no match; a `||` fallback would append.
NMARK="n/a"
if [ -n "$MARKER" ]; then
  NMARK=$(grep -c "$MARKER" "$OUT" || true)
fi
echo "OK $OUT: $NWIN windows, markers($MARKER)=$NMARK"
if [ -n "$PC_SUFFIX" ]; then
  PC_REPORT=${STEM}_report.txt
  "$_HS_LOCK_PYTHON" tools/pc_profile.py "$OUT" >"$PC_REPORT" ||
    { echo "PC REPORT FAILED: $PC_REPORT"; exit 1; }
  echo "PC report: $PC_REPORT"
  sed -n '3,30p' "$PC_REPORT"
fi
