#!/usr/bin/env python3
"""Seed an isolated KiCad configuration with the stock library tables.

usage: kicad_config.py <config-home>   the KICAD_CONFIG_HOME the PCB tests use
"""

from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "hardware/phantasm/gen"))
from kicad_environment import generation_blockers
from kicad_common import KICAD_VERSION, kicad_cli


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    blockers = generation_blockers()
    if blockers:
        print("; ".join(blockers), file=sys.stderr)
        return 1
    config = Path(sys.argv[1]) / ".".join(KICAD_VERSION.split(".")[:2])
    cli = kicad_cli()
    bin_dir = Path(shutil.which(cli) or cli).resolve().parent
    prefix = bin_dir.parent
    bundled = bin_dir.name == "MacOS" and prefix.name == "Contents" and prefix.parent.suffix == ".app"
    templates = prefix / ("SharedSupport/template" if bundled else "share/kicad/template")
    tables = [templates / name for name in ("fp-lib-table", "sym-lib-table")]
    for table in tables:
        if not table.is_file():
            print(f"KiCad library table {table.name} not found in {templates}", file=sys.stderr)
            return 1
    config.mkdir(parents=True, exist_ok=True)
    for table in tables:
        shutil.copyfile(table, config / table.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
