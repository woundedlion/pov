# Effects excluded from the Phantasm image

Composed-effect rows captured before 2026-10-17 predate per-column temporal sampling; [GridSpace's matched A/B](../shipping/profile_gridspace_teensy_2026-10-17.md) measures about +11% mean live-render time for that path.

JewelMelt, KaleidoscopeMelt, KaleidoscopeStainedGlass and RingShower remain available in the full roster and are omitted from Phantasm by product curation. These captures preserve their earlier segmented-device measurements.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [JewelMelt selective o3](profile_jewelmelt_selective_o3_teensy_2026-10-14.md) ● | `fx_shader_draw` | 🟢 46.35 | 🟢 0/1087 (0.00%) | 2026-10-14 20:16 |
| [KaleidoscopeMelt selective o3](profile_kaleidoscopemelt_selective_o3_teensy_2026-10-11.md) ● | `fx_shader_draw` | 🟢 46.26 | 🟢 0/1087 (0.00%) | 2026-10-11 11:38 |
| [KaleidoscopeStainedGlass selective o3](profile_kaleidoscopestainedglass_selective_o3_teensy_2026-10-06.md) | `fx_shader_draw` | 🟢 43.00 | 🟢 0/1087 (0.00%) | 2026-10-06 18:50 |
| [KaleidoscopeStainedGlass global o3](profile_kaleidoscopestainedglass_global_o3_teensy_2026-08-26.md) | `fx_shader_draw` | 🟢 46.99 | 🟢 0/1088 (0%) | 2026-08-26 02:51 |
| [RingShower selective o3](profile_ringshower_selective_o3_teensy_2026-10-06.md) | `rsh_draw_rings` | 🟢 4.32 | 🟢 0/1087 (0.00%) | 2026-10-06 18:20 |
| [RingShower global o3](profile_ringshower_global_o3_teensy_2026-08-26.md) | `rsh_draw_rings` | 🟢 3.86 | 🟢 0/1088 (0%) | 2026-08-26 01:36 |

Low-resolution-only (96×20) effects:

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [MobiusRings Hole A/B](profile_mobiusrings_teensy_2026-10-13.md) | `mr_draw_grid` | 🟢 12.30 | 🟢 0/1088 (0.00%) | 2026-10-13 |
