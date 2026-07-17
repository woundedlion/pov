#!/bin/bash
# Advisory host-global lock for the bench Teensy. Source it, then wrap any
# flash+capture in hs_device_acquire / hs_device_release (profile_one.sh does).
#
# The lock is host-global, NOT repo-local: the device is one physical board on
# one COM port, while concurrent sessions each work from their own worktree
# with their own build/. A lock under build/ would hand every worktree its own
# lock and every holder a green light. Path: $HS_DEVICE_LOCK, else
# $TMPDIR/holosphere-teensy-device.d — a fixed location outside every worktree.
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
#   HS_DEVICE_LOCK   override the lock path
#   HS_DEVICE_WAIT   seconds to wait for a busy device (default 0 = fail fast)
#   HS_DEVICE_FORCE  1 = break someone else's lock (see the warning below)

HS_DEVICE_STALE_GRACE=${HS_DEVICE_STALE_GRACE:-120}
_hs_lock_dir() {
  echo "${HS_DEVICE_LOCK:-${TMPDIR:-${TMP:-/tmp}}/holosphere-teensy-device.d}"
}
_hs_now() { date +%s; }

# Our claim token: only the holder may release, so a stale-break followed by a
# late release from the evicted owner cannot unlock the new holder's device.
_HS_TOKEN=""

_hs_lock_field() {  # <dir> <field>
  sed -n "s/^$2=//p" "$1/info" 2>/dev/null | head -1
}

_hs_lock_born() {  # <dir> — claim mtime in epoch seconds, empty if unknown
  stat -c %Y "$1" 2>/dev/null
}

_hs_holder_desc() {  # <dir>
  local d=$1
  echo "held by session $(_hs_lock_field "$d" session) (pid $(_hs_lock_field "$d" pid))"
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

# hs_device_acquire <effect> <env> <eta_seconds>
# Blocks per HS_DEVICE_WAIT, else fails fast (rc 1) naming the holder.
hs_device_acquire() {
  local effect=$1 env=$2 eta=$3
  local d waited=0 wait_for=${HS_DEVICE_WAIT:-0}
  d=$(_hs_lock_dir)
  while :; do
    if mkdir "$d" 2>/dev/null; then
      _HS_TOKEN="$$-$(_hs_now)-$RANDOM"
      local now; now=$(_hs_now)
      # Written before we hand out the lock so a peer never reads a half-claim.
      {
        echo "token=$_HS_TOKEN"
        echo "session=${CLAUDE_SESSION_ID:-${HS_SESSION:-local}}"
        echo "pid=$$"
        echo "host=$(hostname)"
        echo "effect=$effect"
        echo "env=$env"
        echo "started=$now"
        echo "started_h=$(date '+%H:%M:%S')"
        echo "deadline=$((now + eta))"
        echo "deadline_h=$(date -d "@$((now + eta))" '+%H:%M:%S' 2>/dev/null || echo '?')"
      } >"$d/info"
      return 0
    fi
    if _hs_lock_is_stale "$d"; then
      echo "device lock is stale (holder gone or past its ETA) — breaking it" >&2
      _hs_holder_desc "$d" >&2
      rm -rf "$d"
      continue
    fi
    if [ "${HS_DEVICE_FORCE:-0}" = "1" ]; then
      echo "HS_DEVICE_FORCE=1 — breaking a LIVE device lock" >&2
      _hs_holder_desc "$d" >&2
      rm -rf "$d"
      continue
    fi
    if [ "$wait_for" -gt 0 ] && [ "$waited" -lt "$wait_for" ]; then
      [ "$waited" = 0 ] && { echo "DEVICE BUSY — waiting up to ${wait_for}s" >&2
                             _hs_holder_desc "$d" >&2; }
      sleep 5; waited=$((waited + 5)); continue
    fi
    echo "DEVICE BUSY: $(_hs_holder_desc "$d")" >&2
    echo "Another agent is using the Teensy. Wait for it, set HS_DEVICE_WAIT=<s>," >&2
    echo "or coordinate with that session. Do NOT flash: an upload now would" >&2
    echo "corrupt its capture and may silently not flash yours." >&2
    return 1
  done
}

# Releasing only our own claim keeps an evicted holder's teardown from
# unlocking the device out from under whoever legitimately took it next.
hs_device_release() {
  local d; d=$(_hs_lock_dir)
  [ -n "$_HS_TOKEN" ] || return 0
  if [ "$(_hs_lock_field "$d" token)" = "$_HS_TOKEN" ]; then
    rm -rf "$d"
  fi
  _HS_TOKEN=""
}

hs_device_status() {
  local d; d=$(_hs_lock_dir)
  if [ ! -d "$d" ]; then echo "device free ($d)"; return 0; fi
  _hs_lock_is_stale "$d" && { echo "device lock STALE (breakable):"; _hs_holder_desc "$d"; return 0; }
  echo "device BUSY:"; _hs_holder_desc "$d"; return 1
}

# `bash tools/device_lock.sh status` for a quick check from any shell.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  case "${1:-status}" in
    status) hs_device_status;;
    *) echo "usage: $0 status   (acquire/release are for sourcing)" >&2; exit 2;;
  esac
fi
