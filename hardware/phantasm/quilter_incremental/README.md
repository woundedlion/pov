# Historical incremental Quilter routing input

This project is the protected snapshot that was used to complete these two
connections:

- U1 pin 13 to `MASTER_EN`
- U1 pin 11 to R_PD pin 2 on `SYNC_PULLDOWN`

The completed board is now `../phantasm.kicad_pcb`. Do not use this historical
snapshot for fabrication.

In Quilter:

1. Confirm that every component is shown as pre-placed and that there are pins
   left to route. If it reports zero pins to route, stop and recheck the files.
2. Preserve the input four-layer stackup and the 4 mil / 4 mil fabrication
   rules.
3. Select **Preserve copper on internal layers**.
4. Keep the pre-routed traces enabled. Quilter treats them as fixed routing.
5. Review the two bypass-capacitor assignments, but do not allow component
   placement changes; the validated C1 placements are already inside the
   outline.

For a new full placement run, use the project in `../unplaced/` and follow the
current instructions in `../README.md`.
