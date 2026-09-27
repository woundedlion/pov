# HyperLattice repeatability capture, 2026-09-27

Both captures use clean source `6a7d7160e180e55da8cdd85f0d9d0bfee2cdb478`,
the two authored presets, COM3, 100 seconds, and 16-frame counter windows.
The supported `tools/profile_one.sh` wrapper held the board lock throughout
each build, flash, and capture. Both untouched logs pass
`python tools/parse_profile.py <capture> validate`, including both presets and
the wrap, with no epoch reset. Shipping captured at 01:02 local time; O3 at 01:05.

- `ship_cycle.txt` and `o3_cycle.txt` retain the raw bytes and setup frame.
- Each `.provenance` records source, compiler, ELF hashes, environment hashes,
  and the local artifact directory containing the retained ELF binaries.
- `*_envdump.json` preserves each dump as UTF-8 text with original line endings.
  Encoding `content` as UTF-8 reproduces the bytes whose SHA-256 is recorded
  in both the JSON and corresponding provenance. The build text files are
  readable logs with trailing whitespace removed.
- Each `*_summary.json` records all-live telemetry, window scopes, selected
  held/transition windows, ISR statistics, and the raw capture hash.
- `comparison.json` intersects current and previous frame numbers and compares
  exactly frames 2–1577. The previous shipping capture lasted 120 seconds;
  comparing its full aggregate against this 100-second capture changes phase
  weighting and would give a misleading apparent speedup.

Runtime statistics exclude only setup frame 1 and retain all subsequent rows,
including nine trailing frames after the last complete window. Counter/ISR
statistics exclude the first entire window and use frames 17–1568. Clean-hold
windows exclude transition marker frames through 240 subsequent frames.
All 1,576 live frames in each configuration stay below 62.5 ms.

No shipping code changed since the previous captured source `8441a47efd` in
`core/`, `effects/`, `hardware/`, `targets/`, or `platformio.ini`. These captures
precede the new experimental-preset work and do not measure its performance.
Both wrapper runs also build the ordinary `phantasm` image successfully;
its size is RAM1 code 195,544 B, variables 314,784 B, padding 1,064 B, free
12,896 B. The report-only update does not change engine or driver headers.
