"""Validate candidates and publish source-bound fabrication releases."""

import argparse
import contextlib
from decimal import Decimal
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import board_metadata
import connectivity
import fab
import pcb
import sexp
from kicad_common import F, is_copper_pour, net_name

PROJECT = "phantasm"
TRANSIENT_DIRS = {".history", "out", ".git", "releases"}
SOURCE_SUFFIXES = {".kicad_pcb", ".kicad_pro", ".kicad_sch", ".kicad_sym",
                   ".kicad_mod", ".kicad_dru", ".kicad_wks", ".wrl", ".step", ".stp"}
TABLES = {"fp-lib-table", "sym-lib-table"}
CANDIDATE = Path(__file__).resolve().parents[1] / "1.2" / "routed"


def sources(project):
    required = {f"{PROJECT}{suffix}" for suffix in
                (".kicad_pcb", ".kicad_pro", ".kicad_sch")}
    paths = sorted(path for path in project.rglob("*")
                   if path.is_file() and not TRANSIENT_DIRS.intersection(path.relative_to(project).parts) and
                   (path.suffix in SOURCE_SUFFIXES or path.name in TABLES))
    if not required <= {path.name for path in paths if path.parent == project}:
        raise ValueError("project requires matching phantasm PCB, schematic and project files")
    if any(path.is_symlink() for path in paths):
        raise ValueError("project sources must not be symbolic links")
    return {path.relative_to(project).as_posix(): fab.sha256_file(path) for path in paths}


def check_candidate(project):
    status = json.loads((project / "candidate.json").read_text(encoding="utf-8"))
    if (status.get("revision") != project.parent.name or status.get("status") != "unvalidated"
            or status.get("ordered") is not False or status.get("producer") != "Quilter"
            or not status.get("provenance")):
        raise ValueError("candidate must declare revision, unvalidated status, Quilter provenance and ordered=false")
    if status.get("source_sha256") != sources(project):
        raise ValueError("candidate source hashes changed; update provenance and source_sha256")
    return status


@contextlib.contextmanager
def project_context(project, output, revision):
    values = {"PROJ": str(project), "OUT": str(output), "JLC": str(output / "jlc"),
              "PCB": str(project / f"{PROJECT}.kicad_pcb"),
              "SCH": str(project / f"{PROJECT}.kicad_sch"),
              "PRO": str(project / f"{PROJECT}.kicad_pro")}
    if revision == "1.2":
        values["KNOWN_PARITY_WARNING_COUNTS"] = pcb.GENERATED_LAND_DEVIATIONS.copy()
        values["OPTIONAL_PARITY_WARNING_MAX"] = pcb.GENERATED_OPTIONAL_DEVIATIONS.copy()
    previous = {key: getattr(fab, key) for key in values}
    try:
        for key, value in values.items():
            setattr(fab, key, value)
        yield
    finally:
        for key, value in previous.items():
            setattr(fab, key, value)


def fabrication_options(project, revision):
    metadata = board_metadata.load_board(project / f"{PROJECT}.kicad_pcb")
    if revision == "1.2":
        if metadata.copper_finish != "ENIG":
            raise ValueError(f"copper finish is {metadata.copper_finish!r}; rev 1.2 requires ENIG")
        expected = (("F.Cu", Decimal("0.035")), ("In1.Cu", Decimal("0.0152")),
                    ("In2.Cu", Decimal("0.0152")), ("B.Cu", Decimal("0.035")))
        if metadata.thickness_mm != Decimal("1.6") or metadata.copper_stackup != expected:
            raise ValueError("rev 1.2 requires the 1.6 mm four-layer generator copper stackup")
    board = fab.read_board(project / f"{PROJECT}.kicad_pcb")
    stackup = F(F(board, "setup")[0], "stackup")[0]
    dielectrics = [[str(layer[1]), str(sexp.val(layer, "type")[0]),
                    str(sexp.val(layer, "thickness")[0]),
                    str(sexp.val(layer, "material")[0])]
                   for layer in F(stackup, "layer")
                   if sexp.val(layer, "type") in (["core"], ["prepreg"])]
    if revision == "1.2":
        if dielectrics != [["dielectric 1", "prepreg", "0.2104", "FR4"],
                           ["dielectric 2", "core", "1.065", "FR4"],
                           ["dielectric 3", "prepreg", "0.2104", "FR4"]]:
            raise ValueError("rev 1.2 requires the generator FR4 dielectric stackup")
        ground_layers = {layer for zone in F(board, "zone")
                         if is_copper_pour(zone) and str(sexp.val(zone, "net_name", [net_name(zone)])[0]).lstrip("/") == "GND"
                         for layer in connectivity.zone_layers(zone, metadata.copper_layers)}
        if not {"In1.Cu", "In2.Cu"} <= ground_layers:
            raise ValueError("rev 1.2 requires GND pours on both inner copper layers")
    return {"copper_finish": metadata.copper_finish,
            "dielectric_stackup_mm": dielectrics,
            "thickness_mm": str(metadata.thickness_mm),
            "copper_stackup_mm": [[name, str(value)] for name, value in metadata.copper_stackup],
            "assembly_side": fab.ASSEMBLY_SIDE}


def preflight(project, revision, check_copper=True):
    board_path = project / f"{PROJECT}.kicad_pcb"
    board = fab.read_board(board_path)
    checks = {
        "fabrication_options": lambda: fabrication_options(project, revision),
        "mechanical_width": lambda: board_metadata.validate_mechanical_width(
            board_metadata.load_board(board_path)),
        "project_rules": lambda: fab.validate_project_rules(project / f"{PROJECT}.kicad_pro"),
        "zones": lambda: fab.validate_zone_geometry(board_path, board=board),
        "vias": lambda: fab.validate_via_geometry(board_path, min_vias=1, board=board),
        "solder_mask": lambda: fab.validate_solder_mask(board_path, board=board),
        "plot_origin": lambda: fab.validate_plot_origin(board_path, board=board),
    }
    failures = []
    for name, check in checks.items():
        try:
            check()
        except ValueError as error:
            failures.append(f"{name}: {error}")
    if failures:
        raise ValueError("\n".join(failures))
    if not check_copper:
        return
    broken = connectivity.opens(board)
    if broken:
        raise ValueError(connectivity.report(broken))


def snapshot(project, destination):
    hashes = sources(project)
    destination.mkdir()
    for name in hashes:
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(project / name, target)
    if sources(destination) != hashes or sources(project) != hashes:
        raise ValueError("sources changed while snapshotting; retry")
    return hashes


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def inventory(directory):
    return {path.relative_to(directory).as_posix(): fab.sha256_file(path)
            for path in sorted(directory.rglob("*")) if path.is_file()
            and path != directory / "release.json"}


def verify(directory):
    record = json.loads((directory / "release.json").read_text(encoding="utf-8"))
    if record.get("status") != "validated" or record.get("files") != inventory(directory):
        raise ValueError("release contents do not match validated manifest")
    if record.get("source_sha256") != sources(directory / "source"):
        raise ValueError("release source hashes do not match")
    if record.get("evaluated_source_sha256") != sources(directory / "evaluated"):
        raise ValueError("evaluated source hashes do not match")
    fab.verify_package(str(directory / "output" / "jlc"),
                       str(directory / "output" / "jlc" / fab.SUMS_FILE))
    return fab.sha256_file(directory / "release.json")


def build(project, revision, destination=None):
    if destination is not None and (destination == project or project in destination.parents):
        raise ValueError("release output must be outside the source project")
    if destination is not None and destination.exists():
        raise ValueError("release destination already exists; choose a new directory")
    parent = destination.parent if destination else None
    if parent:
        parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="phantasm-release-", dir=parent))
    try:
        source = stage / "source"
        hashes = snapshot(project, source)
        preflight(source, revision, check_copper=False)
        evaluated = stage / "evaluated"
        snapshot(source, evaluated)
        output = stage / "output"
        output.mkdir()
        subprocess.run([fab.kicad_cli(), "pcb", "drc", "--refill-zones", "--save-board",
                        "--format", "json", "-o", str(output / "refill-drc.json"),
                        str(evaluated / f"{PROJECT}.kicad_pcb")], check=True)
        evaluated_hashes = sources(evaluated)
        preflight(evaluated, revision)
        with project_context(evaluated, output, revision):
            fab.fabricate(revision=revision, validate_only=destination is None, refill_zones=False)
        if sources(evaluated) != evaluated_hashes or sources(source) != hashes or sources(project) != hashes:
            raise ValueError("sources changed during validation/export; release refused")
        if destination is None:
            shutil.rmtree(stage)
            print("Validation passed; no release created.")
            return
        version = subprocess.check_output([fab.kicad_cli(), "--version"], text=True).strip()
        write_json(stage / "release.json", {
            "status": "validated", "revision": revision, "ordered": False,
            "source_sha256": hashes, "evaluated_source_sha256": evaluated_hashes,
            "fabrication_options": fabrication_options(evaluated, revision),
            "tools": {"kicad_cli": version, "python": sys.version,
                      "generator_sha256": {path.name: fab.sha256_file(path)
                                           for path in sorted(Path(__file__).parent.glob("*.py"))}},
            "validation": ["mechanical_width", "project_rules", "zones", "vias", "solder_mask", "plot_origin",
                           *(["finish_stackup"] if revision == "1.2" else []),
                           "connectivity", "ERC", "DRC", "parity",
                           "netlist", "assembly", "export_content"],
            "files": inventory(stage),
        })
        verify(stage)
        stage.rename(destination)
    except BaseException:
        print(f"Reports kept in {stage}", file=sys.stderr)
        raise
    print(f"Release: {destination}\nRelease SHA256: {verify(destination)}")


def record_order(package, expected_hash, details, output):
    if output.resolve() == package.resolve() or package.resolve() in output.resolve().parents:
        raise ValueError("order record must be outside the release package")
    actual = verify(package)
    if actual != expected_hash:
        raise ValueError("release SHA256 differs from the reviewed package")
    release = json.loads((package / "release.json").read_text(encoding="utf-8"))
    for field in ("vendor", "order_id", "ordered_at"):
        if not isinstance(details.get(field), str) or not details[field].strip():
            raise ValueError(f"order details require {field}")
    for field in ("assembly_preview_reviewed", "mechanical_access_reviewed"):
        if details.get(field) is not True:
            raise ValueError(f"order details require {field}=true")
    if details.get("fabrication_options") != release["fabrication_options"]:
        raise ValueError("actual order options differ from the validated fabrication options")
    with output.open("x", encoding="utf-8") as stream:
        json.dump({"release_sha256": actual, "revision": release["revision"],
                   "order": details}, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "release"):
        command = commands.add_parser(name)
        command.add_argument("--revision", choices=("1.1", "1.2"), required=True)
        command.add_argument("--project", type=Path, required=True)
        if name == "release":
            command.add_argument("--output", type=Path, required=True)
    command = commands.add_parser("verify")
    command.add_argument("package", type=Path)
    command = commands.add_parser("record-order")
    command.add_argument("--package", type=Path, required=True)
    command.add_argument("--release-sha256", required=True)
    command.add_argument("--details", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    commands.add_parser("check-candidates")
    args = parser.parse_args(argv)
    try:
        if args.command == "check-candidates":
            projects = sorted(CANDIDATE.parents[1].glob("*/routed"))
            if not projects or list(CANDIDATE.parents[1].glob("*/final")):
                raise ValueError("expected routed candidates with provenance; final directories are forbidden")
            for project in projects:
                check_candidate(project)
            print("Candidate status and provenance verified; fabrication readiness is not asserted.")
        elif args.command == "record-order":
            record_order(args.package.resolve(), args.release_sha256,
                         json.loads(args.details.read_text(encoding="utf-8")), args.output)
        elif args.command == "verify":
            print(f"Release SHA256: {verify(args.package.resolve())}")
        else:
            build(args.project.resolve(), args.revision,
                  args.output.resolve() if args.command == "release" else None)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
