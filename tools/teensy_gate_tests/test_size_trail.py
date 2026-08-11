#!/usr/bin/env python3
"""Host tests for the per-commit size trail (tools/teensy_size_trail.py).

Two things carry the whole feature and both are pure functions, so both are
tested here rather than against a firmware binary:
  * the ELF32 section-header parser — the reason the recorder needs no ARM
    toolchain. Exercised against a SYNTHETIC ELF assembled in-process (no
    committed firmware image) plus the malformed variants that must raise
    rather than silently report zero bytes.
  * the regression classifier — deltas are per-environment, a first row is
    never a regression, and a shrink is never one either.

Run:  python -m unittest discover -s tools/teensy_gate_tests
"""

import os
import struct
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import teensy_size_trail as tst   # noqa: E402

_EHDR_SIZE = 52
_SHDR_SIZE = 40


def make_elf(sections, *, ei_class=tst._ELFCLASS32, ei_data=tst._ELFDATA2LSB,
             magic=tst._ELF_MAGIC, shstrndx=None, shnum=None):
    """Assemble a minimal ELF32 image carrying `sections` = {name: sh_size}."""
    names = list(sections) + [".shstrtab"]
    strtab = bytearray(b"\0")
    offsets = {}
    for name in names:
        offsets[name] = len(strtab)
        strtab += name.encode() + b"\0"

    str_off = _EHDR_SIZE
    shoff = str_off + len(strtab)
    shoff += -shoff % 4

    headers = [b"\0" * _SHDR_SIZE]                       # SHN_UNDEF
    for name, size in sections.items():
        headers.append(struct.pack(tst._SHDR32, offsets[name], 1, 0x3, 0, 0,
                                   size, 0, 0, 4, 0))
    headers.append(struct.pack(tst._SHDR32, offsets[".shstrtab"], 3, 0, 0,
                               str_off, len(strtab), 0, 0, 1, 0))

    out = bytearray(magic + bytes([ei_class, ei_data, 1, 0]) + b"\0" * 8)
    out += struct.pack(tst._EHDR32, 2, 40, 1, 0, 0, shoff, 0, _EHDR_SIZE, 0, 0,
                       _SHDR_SIZE,
                       len(headers) if shnum is None else shnum,
                       len(headers) - 1 if shstrndx is None else shstrndx)
    assert len(out) == _EHDR_SIZE
    out += strtab
    out += b"\0" * (-len(out) % 4)
    for header in headers:
        out += header
    return bytes(out)


#: Section sizes shaped like a real phantasm link.
FIRMWARE = {
    ".text.headers": 5120,
    ".text.code": 141404,
    ".text.progmem": 1205872,
    ".text.itcm": 195696,
    ".data": 3776,
    ".bss": 309984,
    ".bss.dma": 520064,
    ".comment": 69,
}


class ElfParser(unittest.TestCase):
    def test_reads_every_section_size(self):
        sizes = tst.parse_elf_section_sizes(make_elf(FIRMWARE))
        self.assertEqual(sizes[".text.itcm"], 195696)
        self.assertEqual(sizes[".bss.dma"], 520064)
        self.assertEqual(sizes[".comment"], 69)

    def test_projects_sections_onto_tracked_regions(self):
        regions = tst.regions_from_sections(
            tst.parse_elf_section_sizes(make_elf(FIRMWARE)))
        self.assertEqual(regions["itcm"], 195696)
        self.assertEqual(regions["code"], 141404)
        self.assertEqual(regions["progmem"], 1205872)
        self.assertEqual(regions["data"], 3776)
        self.assertEqual(regions["bss"], 309984)
        self.assertEqual(regions["dma"], 520064)
        # ram1 is derived: ITCM code + DTCM variables, no FlexRAM padding.
        self.assertEqual(regions["ram1"], 195696 + 3776 + 309984)
        self.assertEqual(set(regions), set(tst.REGIONS))

    def test_absent_section_reads_zero(self):
        regions = tst.regions_from_sections(
            tst.parse_elf_section_sizes(make_elf({".text.itcm": 8})))
        self.assertEqual(regions["itcm"], 8)
        self.assertEqual(regions["dma"], 0)
        self.assertEqual(regions["ram1"], 8)

    def test_rejects_non_elf(self):
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(b"MZ\x90\x00" + b"\0" * 512)

    def test_rejects_elf64(self):
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(make_elf(FIRMWARE, ei_class=2))

    def test_rejects_big_endian(self):
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(make_elf(FIRMWARE, ei_data=2))

    def test_rejects_truncated_image(self):
        blob = make_elf(FIRMWARE)
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(blob[:_EHDR_SIZE + 8])
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(blob[:16])

    def test_rejects_out_of_range_string_table_index(self):
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(make_elf(FIRMWARE, shstrndx=99))

    def test_rejects_missing_section_table(self):
        with self.assertRaises(tst.ElfFormatError):
            tst.parse_elf_section_sizes(make_elf(FIRMWARE, shnum=0))


class Collect(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        self.warnings = []

    def _env(self, name, blob):
        (self.root / name).mkdir(parents=True)
        (self.root / name / tst.ELF_NAME).write_bytes(blob)

    def test_skips_environments_that_did_not_build(self):
        self._env("phantasm", make_elf(FIRMWARE))
        found = tst.collect(self.root, ("phantasm", "holosphere"),
                            warn=self.warnings.append)
        self.assertEqual(set(found), {"phantasm"})
        self.assertEqual(len(self.warnings), 1)
        self.assertIn("holosphere", self.warnings[0])

    def test_corrupt_elf_is_skipped_not_fatal(self):
        self._env("phantasm", make_elf(FIRMWARE))
        self._env("holosphere", b"not an elf at all")
        found = tst.collect(self.root, ("phantasm", "holosphere"),
                            warn=self.warnings.append)
        self.assertEqual(set(found), {"phantasm"})
        self.assertIn("cannot read", self.warnings[0])

    def test_no_elf_at_all_yields_nothing(self):
        self.assertEqual(tst.collect(self.root, ("phantasm",),
                                     warn=self.warnings.append), {})


SHA = "d" * 40


class Backfill(unittest.TestCase):
    """`backfill` attributes sizes to the commit that produced them. `pio run`'s
    status cannot: it is non-zero for an over-budget commit that linked fine. The
    ELF mtime is the discriminator — a link restamps it, a failed compile leaves
    the previous commit's ELF untouched."""

    def _backfill(self, *, relink):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            worktree = root / "wt"
            worktree.mkdir()
            elf = root / "build" / "holosphere" / tst.ELF_NAME
            elf.parent.mkdir(parents=True)
            elf.write_bytes(make_elf(FIRMWARE))
            stale = os.stat(elf).st_mtime_ns - 10**9
            os.utime(elf, ns=(stale, stale))
            trail = root / "trail.tsv"

            def fake_pio(cmd, **kwargs):
                if relink:
                    # A link rewrites the ELF byte-identically when every object
                    # came from the cache; only its mtime moves.
                    elf.write_bytes(make_elf(FIRMWARE))
                return subprocess.CompletedProcess(cmd, 0)

            args = types.SimpleNamespace(
                worktree=str(worktree), trail=str(trail), env=["holosphere"],
                build_dir=str(root / "build"), rev_range="a..b")
            with mock.patch.object(tst, "_git", return_value=SHA), \
                 mock.patch.object(
                     tst, "head_stamp",
                     return_value=(SHA, "2026-08-03T00:00:00-07:00", "s")), \
                 mock.patch.object(tst.subprocess, "run", side_effect=fake_pio):
                rc = tst.cmd_backfill(args)
            return rc, tst.read_trail(trail)

    def test_unrelinked_env_is_not_recorded_with_the_previous_elf(self):
        rc, rows = self._backfill(relink=False)
        self.assertEqual(rc, 0)
        self.assertEqual(rows, [])

    def test_identical_relink_is_still_recorded(self):
        rc, rows = self._backfill(relink=True)
        self.assertEqual(rc, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].sha, SHA)
        self.assertEqual(rows[0].sizes,
                         tst.regions_from_sections(FIRMWARE))


#: Absolute on both POSIX and Windows: a drive-less path would be resolved
#: against the test runner's cwd drive and stop matching.
_FAKE_COMMON = Path(Path(__file__).resolve().anchor) / "repo" / ".git"
_FAKE_GITDIR = _FAKE_COMMON / "worktrees" / "wt"


class DefaultPaths(unittest.TestCase):
    """The trail is shared by every worktree; the pending record is not.

    The pending record holds one commit in flight: sharing it would let two
    worktrees committing at once stamp each other's sizes onto the wrong sha.
    """

    @staticmethod
    def _rev_parse(args, cwd=None):
        return {"--git-common-dir": str(_FAKE_COMMON),
                "--absolute-git-dir": str(_FAKE_GITDIR)}[args[1]]

    def test_trail_is_on_the_shared_common_dir(self):
        with mock.patch.object(tst, "_git", self._rev_parse):
            self.assertEqual(tst.default_trail(_FAKE_COMMON.parent),
                             _FAKE_COMMON / tst.TRAIL_NAME)

    def test_pending_is_per_worktree(self):
        with mock.patch.object(tst, "_git", self._rev_parse):
            self.assertEqual(tst.default_pending(_FAKE_COMMON.parent),
                             _FAKE_GITDIR / tst.PENDING_NAME)


def _row(sha, env, subject="s", date="2026-08-03T00:00:00-07:00", **sizes):
    full = {r: 0 for r in tst.REGIONS}
    full.update(sizes)
    return tst.TrailRow(sha=sha, date=date, env=env, sizes=full,
                        subject=subject)


class TrailStorage(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "trail.tsv"

    def test_round_trips_rows_and_writes_one_header(self):
        tst.append_rows(self.path, [_row("a" * 40, "phantasm", itcm=100)])
        tst.append_rows(self.path, [_row("b" * 40, "phantasm", itcm=200)])
        text = self.path.read_text(encoding="utf-8")
        self.assertEqual(text.count("#sha"), 1)
        rows = tst.read_trail(self.path)
        self.assertEqual([r.sha[0] for r in rows], ["a", "b"])
        self.assertEqual(rows[1].sizes["itcm"], 200)

    def test_subject_whitespace_cannot_break_a_row(self):
        tst.append_rows(self.path,
                        [_row("c" * 40, "phantasm", subject="one\ttwo\nthree")])
        rows = tst.read_trail(self.path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].subject, "one two three")

    def test_malformed_rows_are_skipped(self):
        tst.append_rows(self.path, [_row("d" * 40, "phantasm", itcm=1)])
        with self.path.open("a", encoding="utf-8", newline="\n") as fh:
            fh.write("garbage\n")
            fh.write("\t".join(["e" * 40, "date", "phantasm"]
                               + ["nope"] * len(tst.REGIONS) + ["s"]) + "\n")
        self.assertEqual([r.sha[0] for r in tst.read_trail(self.path)], ["d"])

    def test_missing_trail_reads_empty(self):
        self.assertEqual(tst.read_trail(self.path / "nope"), [])


class RegressionClassifier(unittest.TestCase):
    def test_deltas_compare_within_an_environment(self):
        rows = [_row("1" * 40, "phantasm", itcm=100),
                _row("1" * 40, "holosphere", itcm=900),
                _row("2" * 40, "phantasm", itcm=140),
                _row("2" * 40, "holosphere", itcm=900)]
        deltas = tst.compute_deltas(rows)
        self.assertTrue(deltas[0].first and deltas[1].first)
        self.assertEqual(deltas[2].changes["itcm"], 40)
        self.assertEqual(deltas[3].changes["itcm"], 0)

    def test_backfilled_older_rows_are_ordered_by_date_not_append_order(self):
        # `backfill` appends commits older than the hook-written rows; in append
        # order the +40 growth would be reported as a -40 shrink at the parent.
        rows = [_row("2" * 40, "phantasm", date="2026-08-03T00:00:00-07:00",
                     itcm=140),
                _row("1" * 40, "phantasm", date="2026-08-01T00:00:00-07:00",
                     itcm=100)]
        deltas = tst.compute_deltas(rows)
        self.assertEqual([d.row.sha[0] for d in deltas], ["1", "2"])
        self.assertEqual(deltas[1].changes["itcm"], 40)
        self.assertEqual([(r.region, r.delta)
                          for r in tst.find_regressions(rows, region="itcm")],
                         [("itcm", 40)])

    def test_rows_across_a_utc_offset_change_order_by_instant(self):
        # A DST fall-back: 01:30-07:00 (08:30Z) precedes 01:15-08:00 (09:15Z),
        # but the raw strings sort the other way and blame the shrink.
        rows = [_row("1" * 40, "phantasm", date="2026-11-01T01:30:00-07:00",
                     itcm=100),
                _row("2" * 40, "phantasm", date="2026-11-01T01:15:00-08:00",
                     itcm=140, subject="grew")]
        deltas = tst.compute_deltas(rows)
        self.assertEqual([d.row.sha[0] for d in deltas], ["1", "2"])
        self.assertEqual(deltas[1].changes["itcm"], 40)
        self.assertEqual([(r.sha[0], r.delta)
                          for r in tst.find_regressions(rows, region="itcm")],
                         [("2", 40)])

    def test_unparseable_date_sorts_oldest_without_raising(self):
        rows = [_row("1" * 40, "phantasm", date="not-a-date", itcm=100),
                _row("2" * 40, "phantasm", date="2026-08-01T00:00:00-07:00",
                     itcm=140)]
        deltas = tst.compute_deltas(rows)
        self.assertEqual([d.row.sha[0] for d in deltas], ["1", "2"])
        self.assertEqual(deltas[1].changes["itcm"], 40)

    def test_first_row_of_an_environment_is_never_a_regression(self):
        found = tst.find_regressions([_row("1" * 40, "phantasm", itcm=195696)])
        self.assertEqual(found, [])

    def test_shrink_is_not_a_regression(self):
        rows = [_row("1" * 40, "phantasm", itcm=200),
                _row("2" * 40, "phantasm", itcm=100)]
        self.assertEqual(tst.find_regressions(rows), [])

    def test_growth_reported_largest_first_across_regions(self):
        rows = [_row("1" * 40, "phantasm", itcm=100, bss=100, ram1=200),
                _row("2" * 40, "phantasm", itcm=110, bss=400, ram1=510,
                     subject="fat commit")]
        found = tst.find_regressions(rows)
        self.assertEqual([(f.region, f.delta) for f in found],
                         [("ram1", 310), ("bss", 300), ("itcm", 10)])
        self.assertEqual(found[0].sha, "2" * 40)
        self.assertEqual(found[0].subject, "fat commit")
        self.assertEqual(found[0].size, 510)

    def test_min_delta_filters_noise(self):
        rows = [_row("1" * 40, "phantasm", itcm=100),
                _row("2" * 40, "phantasm", itcm=104)]
        self.assertEqual(tst.find_regressions(rows, min_delta=1)[0].delta, 4)
        self.assertEqual(tst.find_regressions(rows, min_delta=8), [])

    def test_region_filter_reports_only_that_region(self):
        rows = [_row("1" * 40, "phantasm", itcm=100, bss=100),
                _row("2" * 40, "phantasm", itcm=110, bss=999)]
        found = tst.find_regressions(rows, region="itcm")
        self.assertEqual([(f.region, f.delta) for f in found], [("itcm", 10)])

    def test_regression_attributes_the_growing_commit(self):
        rows = [_row("1" * 40, "phantasm", itcm=100, subject="innocent"),
                _row("2" * 40, "phantasm", itcm=100, subject="also innocent"),
                _row("3" * 40, "phantasm", itcm=180, subject="guilty")]
        found = tst.find_regressions(rows)
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0].subject, "guilty")


class Rendering(unittest.TestCase):
    def test_show_marks_deltas_and_honours_last(self):
        rows = [_row("1" * 40, "phantasm", itcm=100),
                _row("2" * 40, "phantasm", itcm=140),
                _row("3" * 40, "phantasm", itcm=140)]
        text = tst.render_show(rows, last=2)
        self.assertNotIn("1" * 8, text)
        self.assertIn("+40", text)

    def test_empty_renders_are_not_errors(self):
        self.assertIn("empty", tst.render_show([]))
        self.assertIn("no growth", tst.render_regressions([]))

    def test_negative_last_keeps_every_row_instead_of_dropping_the_oldest(self):
        rows = [_row("1" * 40, "phantasm", itcm=100),
                _row("2" * 40, "phantasm", itcm=140)]
        self.assertIn("1" * 8, tst.render_show(rows, last=-1))

    def test_show_rejects_a_last_below_one(self):
        args = tst.build_parser().parse_args(["show", "--last", "-5"])
        self.assertEqual(tst.cmd_show(args), 2)


if __name__ == "__main__":
    unittest.main()
