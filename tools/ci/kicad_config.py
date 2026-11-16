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
    templates = Path(kicad_cli()).resolve().parents[1] / "share/kicad/template"
    config.mkdir(parents=True, exist_ok=True)
    for table in ("fp-lib-table", "sym-lib-table"):
        shutil.copyfile(templates / table, config / table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
