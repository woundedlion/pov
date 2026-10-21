#!/bin/bash
# Cache of built Teensy images keyed by the committed source tree, the
# PlatformIO environment, the build-flag variables and the installed PlatformIO
# toolchain. A checkout with uncommitted or untracked changes is never cached.
# Entries are shared by every worktree of the repository and hold
# firmware.elf/.hex/.map, build.log, envdump.txt (when captured) and meta.
#
# usage: fw_cache.sh build <env>...   reuse each image or `pio run -e <env>`
#                                     under the checkout lock and store it
#        fw_cache.sh key <env>        print the key; 3 = uncommitted changes
#        fw_cache.sh get <env>        print the entry; 1 = miss, 3 = uncommitted
#        fw_cache.sh tree             print the clean checkout's tree hash
#        fw_cache.sh put <env> <tree> <build-log>
#                                     store <env> from a build of
#                                     clean tree <tree> that just finished
#        fw_cache.sh prune            keep the HS_FW_CACHE_KEEP newest entries
#
# Env knobs:
#   HS_FW_CACHE_DIR   cache root (default <git-common-dir>/fw-cache)
#   HS_FW_CACHE_KEEP  entries kept after a store (default 50)
#   HS_FW_CACHE=0     never reuse or store
#   HS_PIO_JOBS       `pio run -j` for builds
#   HS_FW_BUILD_DIRS_KEEP  flag-set build directories kept (default 8)

# Environment variables PlatformIO folds into a build.
HS_FW_CACHE_VARS="PLATFORMIO_BUILD_FLAGS PLATFORMIO_BUILD_UNFLAGS PLATFORMIO_BUILD_SRC_FLAGS PLATFORMIO_SRC_BUILD_FLAGS PLATFORMIO_BUILD_SRC_FILTER PLATFORMIO_EXTRA_SCRIPTS"
HS_FW_TOOLCHAIN=""
HS_FW_KEY=""
HS_FW_MATERIAL=""
HS_FW_ENTRY=""
HS_FW_IMAGE=""

hs_fw_cache_root() {
  if [ -n "${HS_FW_CACHE_DIR:-}" ]; then
    echo "$HS_FW_CACHE_DIR"
  else
    echo "$(git rev-parse --path-format=absolute --git-common-dir)/fw-cache"
  fi
}

# The checkout's tree hash; rc 3 when anything is uncommitted or untracked.
hs_fw_cache_tree() {
  local status
  status=$(git status --porcelain) || return 2
  if [ -n "$status" ]; then
    echo "fw_cache: the checkout has uncommitted changes; image not cached" >&2
    return 3
  fi
  git rev-parse 'HEAD^{tree}'
}

hs_fw_cache_toolchain() {
  local core=${PLATFORMIO_CORE_DIR:-$HOME/.platformio} manifest version
  pio --version || return 2
  for manifest in "$core"/platforms/teensy*/platform.json \
      "$core"/packages/*teensy*/package.json; do
    [ -f "$manifest" ] || continue
    version=$(sed -n 's/^  "version": *"\([^"]*\)".*/\1/p' "$manifest" | head -n 1)
    printf '%s %s\n' "${manifest#"$core"/}" "$version"
  done
}

# One VAR=value line per build-flag variable, whitespace-normalized.
hs_fw_cache_flags() {
  local var value
  for var in $HS_FW_CACHE_VARS; do
    value=$(printf '%s' "${!var:-}" | tr -s ' \t\r\n' ' ' | sed 's/^ //; s/ $//')
    printf '%s=%s\n' "$var" "$value"
  done
}

# Sets HS_FW_KEY and HS_FW_MATERIAL for <env> under the current environment.
hs_fw_cache_key() {
  local env=$1 tree
  tree=$(hs_fw_cache_tree) || return
  if [ -z "$HS_FW_TOOLCHAIN" ]; then
    HS_FW_TOOLCHAIN=$(hs_fw_cache_toolchain) ||
      { echo "fw_cache: cannot fingerprint the PlatformIO toolchain" >&2; return 2; }
  fi
  HS_FW_MATERIAL=$(
    printf 'fw-cache 1\ntree %s\nenv %s\n' "$tree" "$env"
    hs_fw_cache_flags
    printf '%s\n' "$HS_FW_TOOLCHAIN"
  )
  HS_FW_KEY=$(printf '%s\n' "$HS_FW_MATERIAL" | sha256sum | cut -c1-32)
}

# PlatformIO wipes an environment's build directory when its flags change, so
# each flag set builds under its own directory: .pio/build without flags,
# .pio/flags/<hash> otherwise. PLATFORMIO_BUILD_DIR overrides.
hs_fw_build_dir() {
  local flags
  if [ -n "${PLATFORMIO_BUILD_DIR:-}" ]; then
    echo "$PLATFORMIO_BUILD_DIR"
    return
  fi
  flags=$(hs_fw_cache_flags)
  if ! grep -q '=.' <<<"$flags"; then
    echo .pio/build
  else
    echo ".pio/flags/$(printf '%s\n' "$flags" | sha256sum | cut -c1-12)"
  fi
}

# hs_fw_pio <env> [pio run args...]: `pio run` in the flag set's build
# directory, with -j $HS_PIO_JOBS when set. Keeps the HS_FW_BUILD_DIRS_KEEP
# (default 8) most recently used flag-set directories.
hs_fw_pio() {
  local env=$1 dir stale
  shift
  dir=$(hs_fw_build_dir)
  case $dir in
    .pio/flags/*)
      mkdir -p "$dir" && touch "$dir"
      # shellcheck disable=SC2012
      ls -1t .pio/flags | tail -n +$((${HS_FW_BUILD_DIRS_KEEP:-8} + 1)) |
        while IFS= read -r stale; do rm -rf ".pio/flags/$stale"; done
      ;;
  esac
  PLATFORMIO_BUILD_DIR=$dir pio run -e "$env" ${HS_PIO_JOBS:+-j "$HS_PIO_JOBS"} "$@"
}

# rc 0 when <entry> is an intact image stored under <key>.
hs_fw_cache_check_entry() {
  local entry=$1 key=$2 recorded actual
  [ -f "$entry/meta" ] && [ -f "$entry/firmware.hex" ] && [ -f "$entry/build.log" ] || return 1
  [ "$(sed -n 's/^key //p' "$entry/meta")" = "$key" ] || return 1
  recorded=$(sed -n 's/^elf_sha256 //p' "$entry/meta")
  actual=$(sha256sum "$entry/firmware.elf" 2>/dev/null | cut -d' ' -f1) || return 1
  [ -n "$recorded" ] && [ "$recorded" = "$actual" ]
}

# Sets HS_FW_KEY; on a hit sets HS_FW_ENTRY. rc 0 hit, 1 miss, 3 uncommitted,
# 2 error.
hs_fw_cache_lookup() {
  local env=$1 entry
  HS_FW_ENTRY=""
  hs_fw_cache_key "$env" || return
  [ "${HS_FW_CACHE:-1}" != 0 ] || return 1
  entry=$(hs_fw_cache_root)/$HS_FW_KEY
  if hs_fw_cache_check_entry "$entry" "$HS_FW_KEY"; then
    touch "$entry"
    HS_FW_ENTRY=$entry
    return 0
  fi
  return 1
}

hs_fw_cache_prune() {
  local root keep stale
  root=$(hs_fw_cache_root)
  keep=${HS_FW_CACHE_KEEP:-50}
  [ -d "$root" ] || return 0
  find "$root" -mindepth 1 -maxdepth 1 -name '.tmp.*' -mmin +1440 -exec rm -rf {} + || :
  # shellcheck disable=SC2012
  ls -1t "$root" | tail -n +$((keep + 1)) | while IFS= read -r stale; do
    rm -rf "${root:?}/$stale"
  done
}

# hs_fw_cache_store <env> <image-dir> <build-log> [envdump] [replace]
# Stores the image of the key HS_FW_MATERIAL describes, provided the checkout
# still has that tree; sets HS_FW_ENTRY.
hs_fw_cache_store() {
  local env=$1 dir=$2 log=$3 envdump=${4:-} replace=${5:-} root tmp tree entry
  [ "${HS_FW_CACHE:-1}" != 0 ] || return 1
  [ -n "$HS_FW_KEY" ] && [ "$(sed -n 's/^env //p' <<<"$HS_FW_MATERIAL")" = "$env" ] ||
    { echo "fw_cache: no key for $env" >&2; return 2; }
  tree=$(hs_fw_cache_tree) || return
  [ "$tree" = "$(sed -n 's/^tree //p' <<<"$HS_FW_MATERIAL")" ] ||
    { echo "fw_cache: the checkout changed during the $env build; not cached" >&2; return 1; }
  root=$(hs_fw_cache_root)
  mkdir -p "$root" || return 2
  tmp=$(mktemp -d "$root/.tmp.XXXXXX") || return 2
  if ! cp "$dir/firmware.elf" "$dir/firmware.hex" "$tmp/" ||
     ! cp "$log" "$tmp/build.log" ||
     { [ -f "$dir/firmware.map" ] && ! cp "$dir/firmware.map" "$tmp/"; } ||
     { [ -n "$envdump" ] && ! cp "$envdump" "$tmp/envdump.txt"; }; then
    rm -rf "$tmp"
    return 2
  fi
  {
    printf 'key %s\n%s\n' "$HS_FW_KEY" "$HS_FW_MATERIAL"
    printf 'commit %s\n' "$(git rev-parse HEAD)"
    printf 'elf_sha256 %s\n' "$(sha256sum "$tmp/firmware.elf" | cut -d' ' -f1)"
  } >"$tmp/meta"
  entry=$root/$HS_FW_KEY
  [ -z "$replace" ] || rm -rf "$entry"
  if ! mv -T "$tmp" "$entry" 2>/dev/null; then
    rm -rf "$tmp"
    hs_fw_cache_check_entry "$entry" "$HS_FW_KEY" || return 2
  fi
  HS_FW_ENTRY=$entry
  hs_fw_cache_prune
}

# Reuses or builds <env> and sets HS_FW_IMAGE to the directory to flash.
hs_fw_cache_build() {
  local env=$1 rc=0 log
  if hs_fw_cache_lookup "$env"; then
    HS_FW_IMAGE=$HS_FW_ENTRY
    tail -n 15 "$HS_FW_ENTRY/build.log"
    echo "fw_cache: $env reused $HS_FW_ENTRY"
    return 0
  else
    rc=$?
  fi
  [ "$rc" != 2 ] || return 2
  log=$(mktemp "${TMPDIR:-/tmp}/hs-fw-build.XXXXXX") || return 2
  if ! (set -o pipefail; hs_fw_pio "$env" 2>&1 | tee "$log"); then
    rm -f "$log"
    return 1
  fi
  HS_FW_IMAGE=$(hs_fw_build_dir)/$env
  if [ "$rc" = 1 ] && hs_fw_cache_store "$env" "$HS_FW_IMAGE" "$log"; then
    HS_FW_IMAGE=$HS_FW_ENTRY
    echo "fw_cache: $env stored $HS_FW_ENTRY"
  fi
  rm -f "$log"
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  set -uo pipefail
  usage() {
    echo "usage: $0 build <env>... | key <env> | get <env> | tree | put <env> <tree> <build-log> | prune" >&2
    exit 2
  }
  cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)" || exit 2
  case "${1:-}" in
    build)
      shift
      [ "$#" -gt 0 ] || usage
      # shellcheck source-path=SCRIPTDIR source=device_lock.sh
      . "$(dirname "$0")/device_lock.sh"
      TREE=$PWD
      hs_exit_traps
      acquire_tree_lock || exit "$?"
      status=0
      for env in "$@"; do
        hs_fw_cache_build "$env" || status=1
      done
      exit "$status"
      ;;
    key)
      [ "$#" = 2 ] || usage
      hs_fw_cache_key "$2" || exit
      echo "$HS_FW_KEY"
      ;;
    get)
      [ "$#" = 2 ] || usage
      hs_fw_cache_lookup "$2" || exit
      echo "$HS_FW_ENTRY"
      ;;
    tree)
      [ "$#" = 1 ] || usage
      hs_fw_cache_tree
      ;;
    put)
      [ "$#" = 4 ] || usage
      hs_fw_cache_key "$2" || exit
      [ "$(sed -n 's/^tree //p' <<<"$HS_FW_MATERIAL")" = "$3" ] ||
        { echo "fw_cache: the checkout is no longer at tree $3; not cached" >&2; exit 1; }
      hs_fw_cache_store "$2" "$(hs_fw_build_dir)/$2" "$4" || exit
      echo "$HS_FW_ENTRY"
      ;;
    prune)
      [ "$#" = 1 ] || usage
      hs_fw_cache_prune
      ;;
    *) usage ;;
  esac
fi
