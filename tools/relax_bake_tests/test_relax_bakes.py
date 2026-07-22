#!/usr/bin/env python3
"""Host tests for deterministic relax-bake asset hashing."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import relax_bakes  # noqa: E402


class TextHashing(unittest.TestCase):
    def test_line_endings_have_one_identity(self):
        variants = (b"first\nsecond\n", b"first\r\nsecond\r\n",
                    b"first\rsecond\r")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.txt"
            hashes = []
            for contents in variants:
                path.write_bytes(contents)
                hashes.append(relax_bakes.sha256_text(path))
        self.assertEqual(len(set(hashes)), 1)

    def test_freshness_accepts_crlf_checkout(self):
        source_root = relax_bakes.ROOT
        relative_paths = (
            Path("core/mesh/relax_bakes_generated.h"),
            Path("core/mesh/relax_bakes_manifest.json"),
            *(Path(path) for path in relax_bakes.INPUT_PATHS),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for relative in relative_paths:
                source = source_root / relative
                contents = source.read_bytes()
                contents = contents.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(contents.replace(b"\n", b"\r\n"))

            with (
                mock.patch.object(relax_bakes, "ROOT", root),
                mock.patch.object(
                    relax_bakes, "ASSET",
                    root / "core/mesh/relax_bakes_generated.h"),
                mock.patch.object(
                    relax_bakes, "MANIFEST",
                    root / "core/mesh/relax_bakes_manifest.json"),
            ):
                relax_bakes.check(None)


if __name__ == "__main__":
    unittest.main()
