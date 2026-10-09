#!/bin/bash
# Advisory host-global, PER-BOARD lock for the bench Teensys. Source it, then
# wrap any flash+capture in hs_device_acquire / hs_device_release.
#
# Windows + Git Bash. Boards are enumerated (COMn names) by the PlatformIO
# loader's teensy_ports.exe; an explicit HS_TEENSY_PORT skips enumeration when
# the loader is absent. Acquire claims the first free board and exports
# HS_TEENSY_PORT so flash and capture pin the same device. Claims live outside
# worktrees at "$HS_DEVICE_LOCK-<COMn>.d".
#
# Env knobs:
#   HS_DEVICE_LOCK   lock path base (default $TMPDIR/holosphere-teensy-device)
#   HS_DEVICE_WAIT   seconds to wait for a busy device (default 0 = fail fast)
#   HS_DEVICE_FORCE  1 = break someone else's live lock
#   HS_DEVICE_STALE_GRACE  seconds an incomplete claim, or a dead-owner claim
#                    past its deadline, survives before it is stale (default 120)
#   HS_SESSION      owner label recorded in the claim
#   HS_PYTHON       Python interpreter for host lock operations
#   HS_TEENSY_PORT   pin to one board (COMn) instead of searching for a free one
#   HS_TEENSY_TOOLS  tool-teensy dir holding teensy_ports.exe (board enumeration)

HS_DEVICE_STALE_GRACE=${HS_DEVICE_STALE_GRACE:-120}
_hs_lock_base() {
  echo "${HS_DEVICE_LOCK:-${TMPDIR:-${TMP:-/tmp}}/holosphere-teensy-device}"
}
# <port> names the lock directory for one board.
_hs_lock_dir() {
  local base; base=$(_hs_lock_base)
  base=${base%.d}
  echo "$base-$1.d"
}
_hs_now() { date +%s; }

# Attached Teensys, one COM name per line, in teensy_ports.exe's enumeration
# order (the authority the flash resolves against). A caller-set HS_TEENSY_PORT
# is checked against it. rc 1 = the pin is not attached; rc 2 = enumeration
# failed.
hs_device_ports() {
  local tools=${HS_TEENSY_TOOLS:-${PLATFORMIO_CORE_DIR:-$HOME/.platformio}/packages/tool-teensy}
  local attached="" listing rc=0 enumerated=0
  if [ -x "$tools/teensy_ports.exe" ]; then
    enumerated=1
    listing=$("$tools/teensy_ports.exe" -L 2>/dev/null) || rc=$?
    if [ "$rc" -ne 0 ]; then
      echo "device: $tools/teensy_ports.exe -L failed (rc $rc)." >&2
      echo "Board enumeration is the lock's authority; without it a flash can" >&2
      echo "silently pick another session's board or not happen at all." >&2
      echo "Unplug/replug the boards, or unset HS_TEENSY_TOOLS if it is wrong." >&2
      return 2
    fi
    attached=$(printf '%s\n' "$listing" | awk '$2 ~ /^COM[0-9]+$/ {print $2}')
  elif [ -z "${HS_TEENSY_PORT:-}" ]; then
    echo "device: $tools/teensy_ports.exe not found (set HS_TEENSY_TOOLS)" >&2
    return 2
  fi
  if [ -n "${HS_TEENSY_PORT:-}" ]; then
    if [ "$enumerated" -eq 1 ] &&
       ! printf '%s\n' "$attached" | grep -qxF "$HS_TEENSY_PORT"; then
      echo "device: HS_TEENSY_PORT=$HS_TEENSY_PORT is not attached." >&2
      echo "The loader enumerates: $(printf '%s' "$attached" | tr '\n' ' ')" >&2
      echo "Replug that board, or unset HS_TEENSY_PORT to claim a free one." >&2
      return 1
    fi
    echo "$HS_TEENSY_PORT"
    return 0
  fi
  [ -n "$attached" ] && printf '%s\n' "$attached"
  return 0
}

# Our claim token: only the holder may release, so an evicted owner's late
# release cannot unlock the new holder.
_HS_TOKEN=""
_HS_LOCK_DIR=""
# The board this process holds; also exported as HS_TEENSY_PORT on acquire.
HS_DEVICE_PORT=""
# Set when acquire exported that pin itself, so release can drop it; a pin
# the caller set outlives the claim.
_HS_PORT_EXPORTED=""

_hs_lock_field() {  # <dir> <field>
  sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1
}

_hs_lock_born() {  # <dir> — claim mtime in epoch seconds, empty if unknown
  local born
  born=$(stat -c %Y "$1" 2>/dev/null) ||
    born=$(stat -f %m "$1" 2>/dev/null) || return 0
  printf '%s\n' "$born"
}

_hs_holder_desc() {  # <dir>
  local d=$1
  echo "$(_hs_lock_field "$d" port) held by session $(_hs_lock_field "$d" session) (pid $(_hs_lock_field "$d" pid))"
  echo "  effect=$(_hs_lock_field "$d" effect) env=$(_hs_lock_field "$d" env)" \
       "since $(_hs_lock_field "$d" started_h)"
  echo "  expected free by $(_hs_lock_field "$d" deadline_h)"
}

# Live holders retain ownership until release or an explicit forced eviction.
_hs_lock_is_stale() {  # <dir>
  local d=$1 now deadline pid started born
  now=$(_hs_now)
  pid=$(_hs_lock_field "$d" pid)
  if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then return 1; fi
  deadline=$(_hs_lock_field "$d" deadline)
  started=$(_hs_lock_field "$d" started)
  if [ -z "$deadline" ]; then
    born=$(_hs_lock_born "$d")
    [ -z "$born" ] && return 1                      # cannot date it: leave it alone
    [ "$now" -gt $((born + HS_DEVICE_STALE_GRACE)) ] && return 0
    return 1                                        # mid-write, not abandoned
  fi
  [ "$now" -gt $((deadline + HS_DEVICE_STALE_GRACE)) ] && return 0
  if [ -n "$pid" ] && [ -n "$started" ] && [ "$now" -gt $((started + 60)) ]; then
    return 0
  fi
  return 1
}

_HS_LOCK_HELPER="${BASH_SOURCE[0]//\\//}"
_HS_LOCK_HELPER="$(cd "$(dirname "$_HS_LOCK_HELPER")" && pwd)/device_lock_guard.py"

_hs_resolve_python() {
  [ -n "${_HS_LOCK_PYTHON:-}" ] && return 0
  local candidate
  for candidate in "${HS_PYTHON:-}" python3 python; do
    if "$candidate" --version >/dev/null 2>&1; then
      _HS_LOCK_PYTHON="$candidate"
      return 0
    fi
  done
  echo "device: no working Python for the lock guard (set HS_PYTHON)" >&2
  return 2
}

# _hs_break_lock <dir> <token> — evicts only the claim identified by token.
_hs_break_lock() {
  _hs_resolve_python || return 2
  "$_HS_LOCK_PYTHON" "$_HS_LOCK_HELPER" break "$1" "$2"
}

_hs_break_stale() {
  if [ -n "$2" ]; then
    _hs_break_lock "$1" "$2"
  else
    _hs_resolve_python || return 2
    "$_HS_LOCK_PYTHON" "$_HS_LOCK_HELPER" break-empty "$1" \
      "$(($(_hs_now) - HS_DEVICE_STALE_GRACE))"
  fi
}

# _hs_try_claim <dir> <port> <effect> <env> <eta> — claim-or-fail; success
# pins this shell's HS_TEENSY_PORT to the board won.
_hs_try_claim() {
  local d=$1 port=$2 effect=$3 env=$4 eta=$5
  local token info now result
  _hs_resolve_python || return 2
  token="$$-$(_hs_now)-$RANDOM"; now=$(_hs_now)
  info=$(
    echo "token=$token"
    echo "session=${HS_SESSION:-local}"
    echo "pid=$$"
    echo "host=$(hostname)"
    echo "port=${port}"
    echo "effect=$effect"
    echo "env=$env"
    echo "started=$now"
    echo "started_h=$(date '+%H:%M:%S')"
    echo "deadline=$((now + eta))"
    echo "deadline_h=$(date -d "@$((now + eta))" '+%H:%M:%S' 2>/dev/null || echo '?')"
  )
  if printf '%s\n' "$info" | "$_HS_LOCK_PYTHON" "$_HS_LOCK_HELPER" claim "$d"; then
    :
  else
    result=$?
    [ "$result" -eq 1 ] && return 1
    echo "device: lock guard could not run (exit $result)" >&2
    return 2
  fi
  if [ "$(_hs_lock_field "$d" token)" != "$token" ]; then
    echo "device: cannot record the claim in $d/info — leaving ${port} unclaimed" >&2
    _hs_break_lock "$d" "$token" || :
    return 1
  fi
  # Shell state only once the claim is on disk; a failed claim leaves no pin.
  _HS_TOKEN="$token"
  _HS_LOCK_DIR="$d"
  HS_DEVICE_PORT="$port"
  [ -n "${HS_TEENSY_PORT:-}" ] || _HS_PORT_EXPORTED=1
  export HS_TEENSY_PORT="$port"
  return 0
}

# hs_device_acquire <effect> <env> <eta_seconds>
# Claims the first free attached board and exports HS_TEENSY_PORT for it.
# Blocks per HS_DEVICE_WAIT, else fails fast (rc 1) naming every board's holder.
hs_device_acquire() {
  local effect=$1 env=$2 eta=$3
  local waited=0 wait_for=${HS_DEVICE_WAIT:-0} forced=0 p port d
  _hs_resolve_python || return 2
  local root; root=$(dirname "$(_hs_lock_base)")
  if [ ! -d "$root" ]; then
    echo "device: lock root $root does not exist, so no claim can be recorded." >&2
    echo "Create it, or point HS_DEVICE_LOCK at a base whose parent exists." >&2
    return 2
  fi
  while :; do
    # Re-enumerated every round (a board can be replugged while we wait); rc 2
    # (enumeration broke) must stay distinct from 1 (pin not attached).
    local ports; ports=$(hs_device_ports) || return $?
    if [ -z "$ports" ]; then
      if [ "$wait_for" -gt 0 ] && [ "$waited" -lt "$wait_for" ]; then
        [ "$waited" -ne 0 ] || echo "device: waiting for a Teensy to enumerate" >&2
        sleep 5; waited=$((waited + 5)); continue
      fi
      echo "device: no Teensy is enumerated; attach/replug the board" >&2
      return 1
    fi
    for p in $ports; do
      port=$p
      d=$(_hs_lock_dir "$port")
      if _hs_try_claim "$d" "$port" "$effect" "$env" "$eta"; then
        echo "device: using ${port} (lock $d)" >&2
        return 0
      else
        [ "$?" -eq 1 ] || return 2
      fi
    done
    # Stale locks are broken only once every board is busy.
    for p in $ports; do
      port=$p
      d=$(_hs_lock_dir "$port")
      local token; token=$(_hs_lock_field "$d" token)
      if _hs_lock_is_stale "$d"; then
        # Read while the claim still exists; the break destroys its info.
        local desc; desc=$(_hs_holder_desc "$d")
        if _hs_break_stale "$d" "$token"; then
          echo "device lock is stale (holder gone or past its ETA) — breaking it" >&2
          echo "$desc" >&2
          _hs_try_claim "$d" "$port" "$effect" "$env" "$eta" && {
            echo "device: using ${port} (lock $d)" >&2; return 0; }
        fi
      fi
    done
    if [ "${HS_DEVICE_FORCE:-0}" = "1" ] && [ "$forced" -eq 0 ]; then
      forced=1
      # ports is one board per line; force takes the first.
      p=${ports%%$'\n'*}
      port=$p
      d=$(_hs_lock_dir "$port")
      if _hs_try_claim "$d" "$port" "$effect" "$env" "$eta"; then
        return 0
      fi
      echo "HS_DEVICE_FORCE=1 — breaking a LIVE device lock" >&2
      _hs_holder_desc "$d" >&2
      # A changed token belongs to a peer; report busy.
      if _hs_break_stale "$d" "$(_hs_lock_field "$d" token)" &&
          _hs_try_claim "$d" "$port" "$effect" "$env" "$eta"; then
        return 0
      fi
    fi
    # hs_device_status returns 1 here; `|| :` spares a caller's `set -e`.
    if [ "$wait_for" -gt 0 ] && [ "$waited" -lt "$wait_for" ]; then
      if [ "$waited" = 0 ]; then
        echo "ALL DEVICES BUSY — waiting up to ${wait_for}s" >&2
        hs_device_status >&2 || :
      fi
      sleep 5; waited=$((waited + 5)); continue
    fi
    echo "ALL DEVICES BUSY:" >&2
    hs_device_status >&2 || :
    echo "Every attached Teensy is in use. Wait, set HS_DEVICE_WAIT=<s>, attach" >&2
    echo "another board, or coordinate with those sessions. Do NOT flash: an" >&2
    echo "upload now would corrupt a capture and may silently not flash yours." >&2
    return 1
  done
}

# Releases only our own claim, and drops the pin acquire exported.
hs_device_release() {
  local d=$_HS_LOCK_DIR
  [ -n "$_HS_TOKEN" ] && [ -n "$d" ] || return 0
  _hs_break_lock "$d" "$_HS_TOKEN" || :   # a peer holds it now: leave it alone
  _HS_TOKEN=""
  _HS_LOCK_DIR=""
  # Read by whoever sourced this file, not here.
  # shellcheck disable=SC2034
  HS_DEVICE_PORT=""
  if [ -n "$_HS_PORT_EXPORTED" ]; then
    unset HS_TEENSY_PORT
    _HS_PORT_EXPORTED=""
  fi
}

# Reports every attached board. rc 0 if one is claimable, 1 if none is available,
# and 2 if enumeration fails or the configured pin is unattached.
hs_device_status() {
  local ports p port d rc=1
  ports=$(hs_device_ports) || return 2
  if [ -z "$ports" ]; then
    echo "no Teensy is enumerated"
    return 1
  fi
  for p in $ports; do
    port=$p
    d=$(_hs_lock_dir "$port")
    if [ ! -d "$d" ]; then
      echo "${port} free ($d)"; rc=0
    elif _hs_lock_is_stale "$d"; then
      echo "${port} lock STALE (breakable): $(_hs_holder_desc "$d")"; rc=0
    else
      echo "${port} BUSY: $(_hs_holder_desc "$d")"
    fi
  done
  return $rc
}

TREE_LOCK=""
TREE_TOKEN=""
acquire_tree_lock() {
  local now info old_token attempt result
  _hs_resolve_python || return 2
  TREE_LOCK="$TREE/.profile-lock"
  now=$(_hs_now)
  TREE_TOKEN="$$-$now-$RANDOM"
  info=$(printf 'token=%s\npid=%s\nstarted=%s\ndeadline=%s\n' \
    "$TREE_TOKEN" "$$" "$now" "$((now + ${1:-900}))")
  for attempt in 1 2; do
    if printf '%s\n' "$info" | "$_HS_LOCK_PYTHON" "$_HS_LOCK_HELPER" claim "$TREE_LOCK"; then
      return 0
    else
      result=$?
      if [ "$result" -ne 1 ]; then
        echo "build checkout lock guard could not run (exit $result): $TREE_LOCK" >&2
        return 2
      fi
    fi
    old_token=$(_hs_lock_field "$TREE_LOCK" token)
    if [ "$attempt" -eq 1 ] && _hs_lock_is_stale "$TREE_LOCK"; then
      _hs_break_stale "$TREE_LOCK" "$old_token" || :
    else
      break
    fi
  done
  echo "build checkout is already claimed: $TREE_LOCK" >&2
  return 1
}

# `tree <command...>` runs a command under the checkout build lock.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  case "${1:-status}" in
    status) hs_device_status;;
    ports) hs_device_ports;;
    tree)
      shift
      [ "$#" -gt 0 ] || { echo "usage: $0 tree <command...>" >&2; exit 2; }
      TREE=$(git -C "$(dirname "$0")" rev-parse --show-toplevel) || exit 2
      cd "$TREE" || exit 2
      trap '[ -z "$TREE_TOKEN" ] || _hs_break_lock "$TREE_LOCK" "$TREE_TOKEN" || :' EXIT
      trap 'exit 130' INT
      trap 'exit 143' TERM
      acquire_tree_lock || exit "$?"
      "$@"
      ;;
    *) echo "usage: $0 status|ports|tree <command...>   (acquire/release are for sourcing)" >&2; exit 2;;
  esac
fi
