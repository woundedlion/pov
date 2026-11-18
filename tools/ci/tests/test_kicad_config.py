"""Tests for seeding KiCad library tables from the selected installation."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import kicad_config


class KiCadConfigTests(unittest.TestCase):
    def configure(self, cli, config):
        with (patch.object(kicad_config, "generation_blockers", return_value=[]),
              patch.object(kicad_config, "kicad_cli", return_value=str(cli)),
              patch.object(sys, "argv", ["kicad_config.py", str(config)])):
            return kicad_config.main()

    def test_copies_both_library_tables_for_each_install_layout(self):
        tables = {"fp-lib-table": b"footprints", "sym-lib-table": b"symbols"}
        for cli, templates in (
            ("KiCad.app/Contents/MacOS/kicad-cli", "KiCad.app/Contents/SharedSupport/template"),
            ("bin/kicad-cli", "share/kicad/template"),
        ):
            with self.subTest(cli=cli), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stock = root / templates
                stock.mkdir(parents=True)
                for name, content in tables.items():
                    (stock / name).write_bytes(content)
                config = root / "config"
                self.assertEqual(self.configure(root / cli, config), 0)
                version = ".".join(kicad_config.KICAD_VERSION.split(".")[:2])
                self.assertEqual({name: (config / version / name).read_bytes() for name in tables},
                                 tables)

    def test_missing_table_names_the_searched_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stock = root / "share/kicad/template"
            stock.mkdir(parents=True)
            (stock / "fp-lib-table").touch()
            config = root / "config"
            errors = io.StringIO()
            with contextlib.redirect_stderr(errors):
                self.assertEqual(self.configure(root / "bin/kicad-cli", config), 1)
            self.assertIn(f"sym-lib-table not found in {stock}", errors.getvalue())
            self.assertFalse(config.exists())
