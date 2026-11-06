# Rev 1.2 routed candidate

This is **unvalidated Quilter output, not ordered**. The original Quilter job ID
and export date are unavailable. [candidate.json](candidate.json) records the
known provenance and hashes of the complete project, including libraries.
Changes to these inputs require updating the hashes and provenance; passing
that check does not establish fabrication readiness.

The committed candidate has `copper_finish "None"` and inner-zone minimum
width/thermal settings below the fabrication floors. Compare its placement and
connector keepouts against the current generated input before fabrication.
The candidate also lacks the `CUT TEENSY VIN-VUSB` silkscreen legend and the
per-pin J2/J3 marks present on the generated rev 1.2 board.
Do not upload this directory to a board house. Use the
[validation and release workflow](../../README.md#validation-and-release).
