# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Capture outcome: **failed — arena out of memory during the shipping traversal**.
Replaces `profile_hankinsolids_teensy_2026-10-22.md`.
The preceding capture is not a current-source timing result.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | Planned 170000 ms capture with 16-frame windows; default shipping traversal; no ordered-cycle override |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HankinSolids profile 170 16` |

Profile ELF SHA-256: `95feb44cb69587d436a637f8d62f67272ffc2286d07e6603aab5540df08f3342`. The adjacent provenance, envdump,
build logs and ELF artifacts remain in `build/prof/`. The profile epoch is
3600000 ms, beyond the planned capture plus 20000 ms attach slack.

## Frame cadence

**New aggregate peak, mean and spilled-frame count are unavailable.**
The wrapper reports a serial drop after 6360 lines and removes the incomplete
raw capture. The final retained frame line is frame 2305, with render
13.338 ms and wall 53.585 ms. The subsequent arena trap prevents completion
of the show segment. That single retained frame is not a peak or mean for
the failed run and does not establish a zero-spill count.

The planned duration is the current Phantasm show segment. Actual elapsed
capture time was not retained. A display window is 62.5 ms. The run fails
to complete its show; frame-budget compliance cannot be established.

## Phase-by-phase readout

The retained exception tail in `build/prof/remaining_groups.log` gives this
allocation diagnostic:

```text
[OOM] Arena @2003d4e0: req 576, offset 24256, pad 0 / cap 24576
HS_CHECK failed: memory.cpp:61: (false) Arena::allocate: out of memory
```

The tail retains only the last counter dump and frame line. It is insufficient
to reconstruct setup, preset coverage, a full phase table or an aggregate.
No complete cycle or cycle wrap is claimed.

### Per-preset table

Unavailable: the incomplete raw log is removed on serial failure. Preset
coverage and complete-window denominators cannot be recovered from the tail.

### Per-pixel figures

Blended-pixel counts are not available from the retained failure tail.

## Column-ISR / DMA marshaling cost

The final retained snapshot has inclusive `isr_wake` 3.68%, `isr_pack` 2.12%
and `isr_dma_submit` 0.21%. Its separately timestamped interval is 1006.323 ms.
Pack and submit are nested in wake and are not additive. This is a snapshot
before the failure, not a capture-wide CPU share or budget guarantee.

## Summary ranking

No current-source peak, mean or spill ranking is available. The index marks
the capture failure rather than carrying the preceding report’s timing forward.
The default shipping traversal reaches an arena-allocation trap.

README cells: peak 🔴 unavailable (trap), spilled 🔴 unavailable (trap).

## Caveats

- The capture is incomplete and supplies no full-show timing aggregate.
- The retained frame and ISR snapshot apply only to their stated scope.
- `parse_profile.py` windows, presets, buckets and validate cannot run without
  the deleted raw log; no validator success is claimed.
- COM4 re-enumerates and accepts the following SphericalHarmonics upload.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`,
`HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures
through the device lock. Failure diagnostics and source attestation are retained
locally; raw captures and images are not checked into the archive.
