#!/bin/bash
# Advisory host-global, PER-BOARD lock for the bench Teensys. Source it, then
# wrap any flash+capture in hs_device_acquire / hs_device_release
# (profile_one.sh does).
#
# The lock is host-global, NOT repo-local: a device is one physical board on
# one COM port, while concurrent sessions each work from their own worktree
# with their own build/. A lock under build/ would hand every worktree its own
# lock and every holder a green light. Path: "$HS_DEVICE_LOCK-<COMn>.d", base
# $HS_DEVICE_LOCK defaulting to $TMPDIR/holosphere-teensy-device — a fixed
# location outside every worktree.
#
# One lock per attached board, not one for the bench: with several Teensys
# plugged in, sessions run in parallel on different boards instead of queueing
# on one. hs_device_acquire enumerates the attached boards and claims the first
# whose lock is free, then exports HS_TEENSY_PORT so the flash and the capture
# both pin to *that* board — the loader's auto-search refuses to choose between
# two boards (prints "Found 2 Teensy boards", flashes nothing, exits SUCCESS)
# and the capture's first-VID-match would read whichever enumerated first.
# HS_TEENSY_PORT set by the caller pins the run to one board and disables the
# search.
#
# Why a lock at all: a peer's `pio run -t upload` mid-capture corrupts the log,
# and an upload issued while the port is held reports SUCCESS *without
# flashing*, so the loser can be either side and neither is told. The window
# spans the whole build+flash+capture, so the lock is taken before the build
# and held through the capture.
#
# Held state lives in one directory created by mkdir (atomic on Windows and
# POSIX alike; a plain -f test would race). `info` inside it names the holder.
#
# Env knobs:
#   HS_DEVICE_LOCK   override the lock path base (per-board suffix still added)
#   HS_DEVICE_WAIT   seconds to wait for a busy device (default 0 = fail fast)
#   HS_DEVICE_FORCE  1 = break someone else's lock (see the warning below)
#   HS_TEENSY_PORT   pin to one board (COMn) instead of searching for a free one
#   HS_TEENSY_TOOLS  tool-teensy dir holding teensy_ports.exe (board enumeration)

HS_DEVICE_STALE_GRACE=${HS_DEVICE_STALE_GRACE:-120}
_hs_lock_base() {
  echo "${HS_DEVICE_LOCK:-${TMPDIR:-${TMP:-/tmp}}/holosphere-teensy-device}"
}
# <port> — the lock dir for one board. Portless (single-board host with no
# loader to enumerate with) keeps the historical bench-wide path.
_hs_lock_dir() {
  local base; base=$(_hs_lock_base)
  base=${base%.d}
  [ -n "$1" ] && echo "$base-$1.d" || echo "$base.d"
}
_hs_now() { date +%s; }

# Attached Teensys, one COM name per line, in the loader's own enumeration
# order. teensy_ports.exe is the authority, not pyserial: the loader has been
# seen listing a board pyserial no longer enumerated, and it is what the flash
# resolves its USB location against. Empty output = enumerate-less host; the
# caller falls back to the portless lock and the loader's auto-search.
hs_device_ports() {
  if [ -n "$HS_TEENSY_PORT" ]; then echo "$HS_TEENSY_PORT"; return 0; fi
  local tools=${HS_TEENSY_TOOLS:-$HOME/.platformio/packages/tool-teensy}
  [ -x "$tools/teensy_ports.exe" ] || return 0
  "$tools/teensy_ports.exe" -L 2>/dev/null | awk '$2 ~ /^COM[0-9]+$/ {print $2}'
}

# Our claim token: only the holder may release, so a stale-break followed by a
# late release from the evicted owner cannot unlock the new holder's device.
_HS_TOKEN=""
_HS_LOCK_DIR=""
# The board this process holds; also exported as HS_TEENSY_PORT on acquire.
HS_DEVICE_PORT=""

_hs_lock_field() {  # <dir> <field>
  sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1
}

_hs_lock_born() {  # <dir> — claim mtime in epoch seconds, empty if unknown
  stat -c %Y "$1" 2>/dev/null
}

_hs_holder_desc() {  # <dir>
  local d=$1
  echo "$(_hs_lock_field "$d" port) held by session $(_hs_lock_field "$d" session) (pid $(_hs_lock_field "$d" pid))"
  echo "  effect=$(_hs_lock_field "$d" effect) env=$(_hs_lock_field "$d" env)" \
       "since $(_hs_lock_field "$d" started_h)"
  echo "  expected free by $(_hs_lock_field "$d" deadline_h)"
}

# A live holder mid-capture must never be evicted, so staleness is deliberately
# conservative: the deadline is the holder's own ETA plus a grace, and a dead
# PID only counts once the claim is old enough that a just-started peer racing
# us cannot be mistaken for a corpse. An unreadable claim gets the same benefit
# of the doubt: acquire mkdirs the lock before it writes info, so a peer that
# just won the mkdir reads back blank for a moment, and calling that stale
# hands two sessions the device at once.
_hs_lock_is_stale() {  # <dir>
  local d=$1 now deadline pid started born
  now=$(_hs_now)
  deadline=$(_hs_lock_field "$d" deadline)
  started=$(_hs_lock_field "$d" started)
  if [ -z "$deadline" ]; then
    born=$(_hs_lock_born "$d")
    [ -z "$born" ] && return 1                      # cannot date it: leave it alone
    [ "$now" -gt $((born + HS_DEVICE_STALE_GRACE)) ] && return 0
    return 1                                        # mid-write, not abandoned
  fi
  [ "$now" -gt $((deadline + HS_DEVICE_STALE_GRACE)) ] && return 0
  pid=$(_hs_lock_field "$d" pid)
  if [ -n "$pid" ] && [ -n "$started" ] && [ "$now" -gt $((started + 60)) ]; then
    kill -0 "$pid" 2>/dev/null || return 0          # holder died mid-run
  fi
  return 1
}

# _hs_try_claim <dir> <port> <effect> <env> <eta> — mkdir-or-fail, then record
# the claim. Success also pins this shell's HS_TEENSY_PORT to the board won, so
# the flash and the capture can never drift onto a peer's device.
_hs_try_claim() {
  local d=$1 port=$2 effect=$3 env=$4 eta=$5
  mkdir "$d" 2>/dev/null || return 1
  _HS_TOKEN="$$-$(_hs_now)-$RANDOM"
  _HS_LOCK_DIR="$d"
  HS_DEVICE_PORT="$port"
  [ -n "$port" ] && export HS_TEENSY_PORT="$port"
  local now; now=$(_hs_now)
  # Written before we hand out the lock so a peer never reads a half-claim.
  {
    echo "token=$_HS_TOKEN"
    echo "session=${CLAUDE_SESSION_ID:-${HS_SESSION:-local}}"
    echo "pid=$$"
    echo "host=$(hostname)"
    echo "port=${port:-auto}"
    echo "effect=$effect"
    echo "env=$env"
    echo "started=$now"
    echo "started_h=$(date '+%H:%M:%S')"
    echo "deadline=$((now + eta))"
    echo "deadline_h=$(date -d "@$((now + eta))" '+%H:%M:%S' 2>/dev/null || echo '?')"
  } >"$d/info"
  return 0
}

# hs_device_acquire <effect> <env> <eta_seconds>
# Claims the first free attached board and exports HS_TEENSY_PORT for it.
# Blocks per HS_DEVICE_WAIT, else fails fast (rc 1) naming every board's holder.
hs_device_acquire() {
  local effect=$1 env=$2 eta=$3
  local waited=0 wait_for=${HS_DEVICE_WAIT:-0} p port d
  while :; do
    # Re-enumerated every round: a board can be plugged in (or replugged onto a
    # new COM name) while we wait, and that board is a free one.
    local ports; ports=$(hs_device_ports)
    [ -n "$ports" ] || ports="-"       # "-" = portless: no loader to ask
    for p in $ports; do
      port=$([ "$p" = "-" ] && echo "" || echo "$p")
      d=$(_hs_lock_dir "$port")
      if _hs_try_claim "$d" "$port" "$effect" "$env" "$eta"; then
        echo "device: using ${port:-auto-search} (lock $d)" >&2
        return 0
      fi
    done
    # Only once every board is busy: breaking a claim is a last resort, so a
    # stale lock on board A must never be preferred over a free board B.
    for p in $ports; do
      port=$([ "$p" = "-" ] && echo "" || echo "$p")
      d=$(_hs_lock_dir "$port")
      if _hs_lock_is_stale "$d"; then
        echo "device lock is stale (holder gone or past its ETA) — breaking it" >&2
        _hs_holder_desc "$d" >&2
        rm -rf "$d"
        _hs_try_claim "$d" "$port" "$effect" "$env" "$eta" && {
          echo "device: using ${port:-auto-search} (lock $d)" >&2; return 0; }
      fi
    done
    if [ "${HS_DEVICE_FORCE:-0}" = "1" ]; then
      port=$([ "${ports%% *}" = "-" ] && echo "" || echo "${ports%% *}")
      d=$(_hs_lock_dir "$port")
      echo "HS_DEVICE_FORCE=1 — breaking a LIVE device lock" >&2
      _hs_holder_desc "$d" >&2
      rm -rf "$d"
      continue
    fi
    if [ "$wait_for" -gt 0 ] && [ "$waited" -lt "$wait_for" ]; then
      [ "$waited" = 0 ] && { echo "ALL DEVICES BUSY — waiting up to ${wait_for}s" >&2
                             hs_device_status >&2; }
      sleep 5; waited=$((waited + 5)); continue
    fi
    echo "ALL DEVICES BUSY:" >&2
    hs_device_status >&2
    echo "Every attached Teensy is in use. Wait, set HS_DEVICE_WAIT=<s>, attach" >&2
    echo "another board, or coordinate with those sessions. Do NOT flash: an" >&2
    echo "upload now would corrupt a capture and may silently not flash yours." >&2
    return 1
  done
}

# Releasing only our own claim keeps an evicted holder's teardown from
# unlocking the device out from under whoever legitimately took it next.
hs_device_release() {
  local d=$_HS_LOCK_DIR
  [ -n "$_HS_TOKEN" ] && [ -n "$d" ] || return 0
  if [ "$(_hs_lock_field "$d" token)" = "$_HS_TOKEN" ]; then
    rm -rf "$d"
  fi
  _HS_TOKEN=""; _HS_LOCK_DIR=""; HS_DEVICE_PORT=""
}

# Reports every attached board. rc 0 if at least one is claimable.
hs_device_status() {
  local ports p port d free=1
  ports=$(hs_device_ports)
  [ -n "$ports" ] || ports="-"
  for p in $ports; do
    port=$([ "$p" = "-" ] && echo "" || echo "$p")
    d=$(_hs_lock_dir "$port")
    if [ ! -d "$d" ]; then
      echo "${port:-auto} free ($d)"; free=0
    elif _hs_lock_is_stale "$d"; then
      echo "${port:-auto} lock STALE (breakable): $(_hs_holder_desc "$d")"; free=0
    else
      echo "${port:-auto} BUSY: $(_hs_holder_desc "$d")"
    fi
  done
  return $free
}

# `bash tools/device_lock.sh status` / `ports` for a quick check from any shell.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  case "${1:-status}" in
    status) hs_device_status;;
    ports) hs_device_ports;;
    *) echo "usage: $0 status|ports   (acquire/release are for sourcing)" >&2; exit 2;;
  esac
fi
