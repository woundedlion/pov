#!/bin/bash
# Source after device_lock.sh; flash the board claimed by hs_device_acquire.

_HS_FLASH_HELPER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/teensy_upload.py"

# hs_teensy_flash <image-dir>: uploads <image-dir>/firmware.hex.
hs_teensy_flash() {
  local image=$1
  local teensy_tools=${HS_TEENSY_TOOLS:-${PLATFORMIO_CORE_DIR:-$HOME/.platformio}/packages/tool-teensy}
  [ -n "${_HS_TOKEN:-}" ] || { echo "flash requires a device lock" >&2; return 1; }
  [ -n "${HS_TEENSY_PORT:-}" ] ||
    { echo "device lock did not pin a Teensy port" >&2; return 1; }
  local backend=${HS_UPLOAD_BACKEND:-auto} resolved resolve_status
  case "$backend" in
    auto)
      if [ -n "${HS_TYCMD:-}" ]; then
        backend=headless
      else
        _hs_resolve_python || return 2
        if resolved=$("$_HS_LOCK_PYTHON" "$_HS_FLASH_HELPER" --resolve 2>&1); then
          backend=headless
        else
          resolve_status=$?
          if [ "$resolve_status" -eq 2 ]; then
            backend=vendor
          else
            printf '%s\n' "$resolved" >&2
            return "$resolve_status"
          fi
        fi
      fi
      ;;
    headless|vendor) ;;
    *) echo "unknown HS_UPLOAD_BACKEND=$backend" >&2; return 1;;
  esac
  if [ "$backend" = headless ]; then
    [ -n "${_HS_LOCK_DIR:-}" ] ||
      { echo "device lock directory is missing" >&2; return 1; }
    _hs_resolve_python || return 2
    "$_HS_LOCK_PYTHON" "$_HS_FLASH_HELPER" \
      --port "$HS_TEENSY_PORT" --lock-dir "$(hs_native_path "$_HS_LOCK_DIR")" \
      --lock-token "$_HS_TOKEN" --absent-seconds "$(hs_device_absent_budget)" \
      "$(hs_native_path "$image/firmware.hex")"
    return $?
  fi
  local line loc label waited=0 budget
  budget=$(hs_device_absent_budget)
  while line=$("$teensy_tools/teensy_ports.exe" -L | grep " $HS_TEENSY_PORT " || true)
        [ -z "$line" ]; do
    [ "$waited" -lt "$budget" ] || { echo "no Teensy at $HS_TEENSY_PORT" >&2; return 1; }
    [ "$waited" -ne 0 ] ||
      echo "device: $HS_TEENSY_PORT is off the port list; waiting up to ${budget}s" >&2
    sleep 5; waited=$((waited + 5))
  done
  loc=$(echo "$line" | awk '{print $1}')
  label=$(echo "$line" | awk '{print $2" "$3" "$4}')
  "$teensy_tools/teensy_post_compile.exe" -file=firmware \
    -path="$(cygpath -aw "$image")" \
    -tools="$(cygpath -w "$teensy_tools")" -board=TEENSY40 -reboot \
    "-port=$loc" "-portlabel=$label" -portprotocol=Teensy
}
