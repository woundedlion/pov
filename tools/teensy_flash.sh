#!/bin/bash
# Source after device_lock.sh; flash the board claimed by hs_device_acquire.

# hs_teensy_flash <image-dir>: uploads <image-dir>/firmware.hex.
hs_teensy_flash() {
  local image=$1
  local teensy_tools=${HS_TEENSY_TOOLS:-${PLATFORMIO_CORE_DIR:-$HOME/.platformio}/packages/tool-teensy}
  [ -n "${_HS_TOKEN:-}" ] || { echo "flash requires a device lock" >&2; return 1; }
  [ -n "${HS_TEENSY_PORT:-}" ] ||
    { echo "device lock did not pin a Teensy port" >&2; return 1; }
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
