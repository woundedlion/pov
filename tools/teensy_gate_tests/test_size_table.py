#!/usr/bin/env python3
"""Host tests for the `just teensy-size` summary table (tools/teensy_size_table.py).

Covers the wrapper's own logic — attributing teensy_size lines to the right env
in a multi-env `pio run` log and rendering them side-by-side. The teensy_size
line PARSING is teensy_gate.parse_teensy_size, already covered by
test_teensy_gate.py.

Run:  python -m unittest discover -s tools/teensy_gate_tests
"""

import contextlib
import io
import sys
import unittest
from unittest import mock
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import teensy_size_table as tst   # noqa: E402

FIX = Path(__file__).resolve().parent / "fixtures"


def _env_chunk(env, teensy_size_fixture=None, extra=()):
    lines = [f"Processing {env} (board: teensy40; framework: arduino; "
             f"platform: teensy@5.0.0)",
             "Building in release mode"]
    if teensy_size_fixture:
        lines += (FIX / teensy_size_fixture).read_text(
            encoding="utf-8").splitlines()
    lines += list(extra)
    return lines


class SplitByEnv(unittest.TestCase):
    def test_attributes_size_lines_to_the_emitting_env(self):
        log = (["PlatformIO Core 6.1.19"]
               + _env_chunk("holosphere", "good_teensy_size.txt")
               + _env_chunk("phantasm", "broken_over_cap_teensy_size.txt"))
        order, sizes = tst.collect_sizes(log)
        self.assertEqual(order, ["holosphere", "phantasm"])
        # good fixture: FLASH code 158,788 / broken fixture: 1,958,788
        self.assertEqual(sizes["holosphere"]["flash"]["components"]["code"],
                         158788)
        self.assertEqual(sizes["phantasm"]["flash"]["components"]["code"],
                         1958788)

    def test_preamble_before_first_env_is_dropped(self):
        # A teensy_size-shaped line before any "Processing" banner must not be
        # attributed to the first env (there is nothing to attribute it to).
        stray = "teensy_size:   RAM2: variables:1   free for malloc/new: 2"
        log = [stray] + _env_chunk("holosphere", "good_teensy_size.txt")
        _, sizes = tst.collect_sizes(log)
        self.assertEqual(sizes["holosphere"]["ram2"]["components"]["variables"],
                         497920)

    def test_env_without_size_output_still_listed(self):
        log = (_env_chunk("holosphere", "good_teensy_size.txt")
               + _env_chunk("phantasm8", extra=["*** [link] Error 1"]))
        order, sizes = tst.collect_sizes(log)
        self.assertEqual(order, ["holosphere", "phantasm8"])
        self.assertEqual(sizes["phantasm8"], {})


class RenderTable(unittest.TestCase):
    def test_renders_all_regions_and_dashes_for_failed_env(self):
        log = (_env_chunk("holosphere", "good_teensy_size.txt")
               + _env_chunk("phantasm8"))
        order, sizes = tst.collect_sizes(log)
        table = tst.render_table(order, sizes)
        header = table.splitlines()[0]
        self.assertIn("holosphere", header)
        self.assertIn("phantasm8", header)
        rows = {line.split("  ")[0].strip(): line
                for line in table.splitlines()[2:]}
        self.assertEqual(set(rows), {label for _, _, label in tst._ROWS})
        # good_teensy_size.txt: RAM1 code 62,240; gate total = 351280+62240+30496.
        self.assertIn("62,240", rows["RAM1 code (ITCM)"])
        self.assertIn("444,016", rows["RAM1 used (gate total)"])
        self.assertIn("26,368", rows["RAM2 free (malloc/new)"])
        # phantasm8 produced no size output -> '-' cells, not a dropped column.
        self.assertTrue(rows["FLASH code"].rstrip().endswith("-"))


class Main(unittest.TestCase):
    class _FakePio:
        def __init__(self, lines, rc):
            self.stdout = io.StringIO("".join(f"{line}\n" for line in lines))
            self.rc = rc

        def wait(self):
            return self.rc

    def _run(self, lines, rc):
        for target, name, value in (
                (tst.shutil, "which", "pio"),
                (tst.subprocess, "Popen", self._FakePio(lines, rc))):
            patcher = mock.patch.object(target, name, return_value=value)
            patcher.start()
            self.addCleanup(patcher.stop)
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()):
            with contextlib.redirect_stderr(err):
                return tst.main(["holosphere"]), err.getvalue()

    def test_malformed_size_line_still_propagates_the_pio_exit_code(self):
        log = _env_chunk("holosphere", extra=[
            "teensy_size:   FLASH: unavailable   free for files: 1883136"])
        rc, err = self._run(log, 3)
        self.assertEqual(rc, 3)
        self.assertIn("unparsable teensy_size output", err)

    def test_a_clean_run_returns_the_pio_exit_code(self):
        log = _env_chunk("holosphere", "good_teensy_size.txt")
        self.assertEqual(self._run(log, 0)[0], 0)

if __name__ == "__main__":
    unittest.main()
