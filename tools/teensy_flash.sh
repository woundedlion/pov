#!/bin/bash
# Source after device_lock.sh; flash the board claimed by hs_device_acquire.

hs_teensy_flash() {
  local env=$1
  local teensy_tools=${HS_TEENSY_TOOLS:-$HOME/.platformio/packages/tool-teensy}
  [ -n "${_HS_TOKEN:-}" ] || { echo "flash requires a device lock" >&2; return 1; }
  [ -n "${HS_TEENSY_PORT:-}" ] ||
    { echo "device lock did not pin a Teensy port" >&2; return 1; }
  local line loc label
  line=$("$teensy_tools/teensy_ports.exe" -L | grep " $HS_TEENSY_PORT " || true)
  [ -n "$line" ] || { echo "no Teensy at $HS_TEENSY_PORT" >&2; return 1; }
  loc=$(echo "$line" | awk '{print $1}')
  label=$(echo "$line" | awk '{print $2" "$3" "$4}')
  "$teensy_tools/teensy_post_compile.exe" -file=firmware \
    -path="$(cygpath -w "$PWD/.pio/build/$env")" \
    -tools="$(cygpath -w "$teensy_tools")" -board=TEENSY40 -reboot \
    "-port=$loc" "-portlabel=$label" -portprotocol=Teensy
}
