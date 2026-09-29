# Sphere-only shell flight capture evidence

Four 30-second fixed-preset captures on COM4, window 16. Raw serial logs include per-frame telemetry and appended provenance. Summaries exclude setup frame 1, while scope/ISR summaries use complete later windows. Prior matching captures remain in [the 2026-09-28 evidence](../hyperlattice_shell_flight_2026-09-28/README.md).

Final ELF/map bundles are preserved outside the removable worker checkout at C:/work/Holosphere/.git/implement-sessions/shells-profile-artifacts. Its archive_manifest.json maps original worker paths to preserved directories and records SHA-256 hashes. Eight ELFs and eight maps were verified after copying; ELF hashes match unchanged provenance. Original worker paths in provenance may disappear after checkout cleanup.

| Preset/config | Raw serial | Validation | Runtime summary | Provenance |
|---|---|---|---|---|
| 5 ship | [log](shells_5_ship.txt) | [validation](shells_5_ship_validate.txt) | [summary](shells_5_ship_summary.json) | [provenance](shells_5_ship.provenance) |
| 5 o3 | [log](shells_5_o3.txt) | [validation](shells_5_o3_validate.txt) | [summary](shells_5_o3_summary.json) | [provenance](shells_5_o3.provenance) |
| 6 ship | [log](shells_6_ship.txt) | [validation](shells_6_ship_validate.txt) | [summary](shells_6_ship_summary.json) | [provenance](shells_6_ship.provenance) |
| 6 o3 | [log](shells_6_o3.txt) | [validation](shells_6_o3_validate.txt) | [summary](shells_6_o3_summary.json) | [provenance](shells_6_o3.provenance) |
