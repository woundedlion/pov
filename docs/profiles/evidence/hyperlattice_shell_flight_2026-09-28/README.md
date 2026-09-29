# HyperLattice shell-flight capture evidence

Four fixed-preset captures on COM4, two authored shell presets and two optimization configurations. Raw serial logs include frame telemetry and appended build provenance. JSON preserves exact build log and flag text with original file hashes.

On the capture host, the four final ELF/map bundles are archived at C:/work/Holosphere/.git/implement-sessions/hyper-five-artifacts. Its archive_manifest.json maps the original worker paths to their preserved locations and records SHA-256 hashes. All eight ELFs and eight maps were verified after copying; ELF hashes also match the unchanged provenance records. The original worker paths in those records are historical and may no longer exist after checkout cleanup.

| Preset/config | Raw serial | Validation | Runtime summary | Provenance |
|---|---|---|---|---|
| 5 ship | [log](hyperlattice_shells_ship.txt) | [validation](hyperlattice_shells_ship_validate.txt) | [summary](hyperlattice_shells_ship_summary.json) | [provenance](hyperlattice_shells_ship.provenance) |
| 5 o3 | [log](hyperlattice_shells_o3.txt) | [validation](hyperlattice_shells_o3_validate.txt) | [summary](hyperlattice_shells_o3_summary.json) | [provenance](hyperlattice_shells_o3.provenance) |
| 6 ship | [log](hyperlattice_shells_close_ship.txt) | [validation](hyperlattice_shells_close_ship_validate.txt) | [summary](hyperlattice_shells_close_ship_summary.json) | [provenance](hyperlattice_shells_close_ship.provenance) |
| 6 o3 | [log](hyperlattice_shells_close_o3.txt) | [validation](hyperlattice_shells_close_o3_validate.txt) | [summary](hyperlattice_shells_close_o3_summary.json) | [provenance](hyperlattice_shells_close_o3.provenance) |
