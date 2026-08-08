#!/usr/bin/env python3
"""Per-commit Teensy 4 firmware size trail (ELF parser / recorder / classifier).

The pre-commit hook already links all three images to run the size gate
(tools/teensy_gate.py) and then discards the numbers. This module captures them
so a later ITCM/RAM1 regression can be attributed to a commit.

Like teensy_gate.py it holds NO PlatformIO dependency and shells out to no ARM
toolchain: section sizes are read straight out of the ELF32 section-header table
(stdlib `struct` only), so it runs as an ordinary host Python module and is
exercised by tools/teensy_gate_tests/test_size_trail.py against a synthetic ELF.

Storage (deliberately NOT a tracked file):
  * The trail lives at `$(git rev-parse --git-common-dir)/teensy-size-trail.tsv`.
    It is local-only — a tracked per-commit log would conflict on every branch,
    and the repo lands branches through a fast-forward-only integrator.
  * `--git-common-dir`, never `--git-dir`: inside a `git worktree` the latter
    points at `.git/worktrees/<name>`, which would fragment the trail into one
    file per worktree. Every worktree appends to the one trail.

Capture is two-phase because pre-commit does not know the commit sha yet:
  * `record`  — pre-commit, after the gate passes: parse the built ELFs into a
    pending JSON record.
  * `commit`  — post-commit: stamp the pending record with HEAD's sha, committer
    date and subject, append one row per environment, drop the pending file.

Neither phase may fail a commit. Every subcommand degrades to a warning on
stderr and a non-zero status the hooks swallow.
"""

from __future__ import annotations

import argparse
import json
import os
import struct
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# ELF32 section-header parsing
# ---------------------------------------------------------------------------
_ELF_MAGIC = b"\x7fELF"
_ELFCLASS32 = 1
_ELFDATA2LSB = 1
#: e_ident[16] then the ELF32 header fields the section table needs.
_EHDR32 = "<HHIIIIIHHHHHH"
_EHDR32_OFF = 16
_SHDR32 = "<IIIIIIIIII"


class ElfFormatError(ValueError):
    """The blob is not an ELF32 little-endian image this parser can read."""


def parse_elf_section_sizes(data: bytes) -> dict[str, int]:
    """Return {section name: sh_size} for an ELF32 little-endian image.

    Section SIZES are the figure of interest; addresses are the gate's business
    (teensy_gate.region_for_address), not the trail's.
    """
    if len(data) < _EHDR32_OFF + struct.calcsize(_EHDR32):
        raise ElfFormatError("file shorter than an ELF32 header")
    if data[:4] != _ELF_MAGIC:
        raise ElfFormatError("bad ELF magic")
    if data[4] != _ELFCLASS32:
        raise ElfFormatError(f"not ELFCLASS32 (EI_CLASS={data[4]})")
    if data[5] != _ELFDATA2LSB:
        raise ElfFormatError(f"not little-endian (EI_DATA={data[5]})")

    (_type, _machine, _version, _entry, _phoff, shoff, _flags, _ehsize,
     _phentsize, _phnum, shentsize, shnum, shstrndx) = struct.unpack_from(
        _EHDR32, data, _EHDR32_OFF)

    if shnum == 0 or shoff == 0:
        raise ElfFormatError("no section header table")
    if shentsize < struct.calcsize(_SHDR32):
        raise ElfFormatError(f"section header entry size {shentsize} too small")
    if shoff + shnum * shentsize > len(data):
        raise ElfFormatError("section header table runs past end of file")
    if shstrndx >= shnum:
        raise ElfFormatError(f"section-name string table index {shstrndx} "
                             f"out of range ({shnum} sections)")

    headers = [struct.unpack_from(_SHDR32, data, shoff + i * shentsize)
               for i in range(shnum)]
    str_off, str_size = headers[shstrndx][4], headers[shstrndx][5]
    if str_off + str_size > len(data):
        raise ElfFormatError("section-name string table runs past end of file")
    strtab = data[str_off:str_off + str_size]

    sizes: dict[str, int] = {}
    for name_off, _sh_type, _sh_flags, _addr, _off, size, *_rest in headers:
        if name_off >= len(strtab):
            raise ElfFormatError(f"section name offset {name_off} outside the "
                                 f"string table")
        end = strtab.find(b"\0", name_off)
        if end < 0:
            raise ElfFormatError("unterminated section name")
        name = strtab[name_off:end].decode("utf-8", "replace")
        if name:
            sizes[name] = size
    return sizes


# ---------------------------------------------------------------------------
# Section -> tracked region
# ---------------------------------------------------------------------------
#: Teensy 4 linker sections the trail follows, in report order.
_SECTION_REGION: tuple[tuple[str, str], ...] = (
    (".text.itcm", "itcm"),        # fast code in RAM1 — the binding figure
    (".text.code", "code"),        # flash-resident code
    (".text.progmem", "progmem"),  # flash const data
    (".data", "data"),             # DTCM initialized variables
    (".bss", "bss"),               # DTCM zero-init variables (the arena)
    (".bss.dma", "dma"),           # OCRAM DMAMEM (the framebuffers)
)

#: `ram1` is derived, not a section: ITCM code + DTCM variables. It excludes the
#: FlexRAM bank padding teensy_size reports, so it tracks growth, not the
#: gate's absolute RAM1 figure.
_RAM1_PARTS = ("itcm", "data", "bss")

#: Every column the trail carries a byte count for.
REGIONS: tuple[str, ...] = tuple(r for _, r in _SECTION_REGION) + ("ram1",)

#: Environments the pre-commit hook links, hence the default record set.
DEFAULT_ENVIRONMENTS: tuple[str, ...] = (
    "holosphere", "phantasm", "holosphere_dma")

ELF_NAME = "firmware.elf"


def regions_from_sections(sizes: dict[str, int]) -> dict[str, int]:
    """Project parsed section sizes onto the tracked regions."""
    out = {region: sizes.get(section, 0) for section, region in _SECTION_REGION}
    out["ram1"] = sum(out[part] for part in _RAM1_PARTS)
    return out


def regions_from_elf(path: str | Path) -> dict[str, int]:
    """Parse one firmware ELF into its tracked region sizes."""
    return regions_from_sections(parse_elf_section_sizes(Path(path).read_bytes()))


# ---------------------------------------------------------------------------
# Trail storage (TSV, one row per commit x environment)
# ---------------------------------------------------------------------------
TRAIL_NAME = "teensy-size-trail.tsv"
PENDING_NAME = "teensy-size-trail.pending.json"

_COLUMNS: tuple[str, ...] = ("sha", "date", "env") + REGIONS + ("subject",)
_HEADER = "#" + "\t".join(_COLUMNS)


@dataclass(frozen=True)
class TrailRow:
    sha: str
    date: str       # committer date, ISO-8601 with offset (%cI)
    env: str
    sizes: dict[str, int]
    subject: str

    def to_tsv(self) -> str:
        fields = [self.sha, self.date, self.env]
        fields += [str(self.sizes.get(r, 0)) for r in REGIONS]
        fields.append(_flatten(self.subject))
        return "\t".join(fields)


def _flatten(text: str) -> str:
    """Collapse anything that would break the one-row-per-line TSV."""
    return " ".join(text.split())


def parse_trail(text: str) -> list[TrailRow]:
    """Parse trail TSV text. Malformed rows are skipped, never fatal."""
    rows: list[TrailRow] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < len(_COLUMNS):
            continue
        try:
            sizes = {r: int(fields[3 + i]) for i, r in enumerate(REGIONS)}
        except ValueError:
            continue
        rows.append(TrailRow(sha=fields[0], date=fields[1], env=fields[2],
                             sizes=sizes, subject=fields[len(_COLUMNS) - 1]))
    return rows


def read_trail(path: str | Path) -> list[TrailRow]:
    p = Path(path)
    if not p.exists():
        return []
    return parse_trail(p.read_text(encoding="utf-8"))


def append_rows(path: str | Path, rows: list[TrailRow]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    new = not p.exists() or p.stat().st_size == 0
    with p.open("a", encoding="utf-8", newline="\n") as fh:
        if new:
            fh.write(_HEADER + "\n")
        for row in rows:
            fh.write(row.to_tsv() + "\n")


# ---------------------------------------------------------------------------
# Deltas + regression classifier
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Delta:
    """One row's byte change per region against the previous row of its env."""

    row: TrailRow
    changes: dict[str, int]   # region -> signed delta; empty for a first row

    @property
    def first(self) -> bool:
        return not self.changes


#: Sort key for an unparseable committer date: oldest, so one bad row cannot
#: displace the real ordering around it.
_UNDATED = datetime(1970, 1, 1, tzinfo=timezone.utc)


def instant(date: str) -> datetime:
    """The absolute instant of a trail row's ISO-8601 committer date."""
    try:
        parsed = datetime.fromisoformat(date)
    except ValueError:
        return _UNDATED
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def compute_deltas(rows: list[TrailRow]) -> list[Delta]:
    """Annotate trail rows with per-region deltas, comparing within an env.

    Rows are ordered by committer date, not append order: `backfill` appends
    commits older than everything the hooks already recorded, and comparing
    those in append order reports every delta backwards. Ordering on the raw
    `%cI` string would compare local wall clocks, so rows recorded under two
    UTC offsets (a DST shift, a travelling committer) can invert a delta pair
    and blame the commit that shrank; sort on the parsed instant instead. The
    sort is stable, so same-instant rows keep append order.
    """
    rows = sorted(rows, key=lambda r: instant(r.date))
    previous: dict[str, TrailRow] = {}
    out: list[Delta] = []
    for row in rows:
        prev = previous.get(row.env)
        changes = ({r: row.sizes.get(r, 0) - prev.sizes.get(r, 0)
                    for r in REGIONS} if prev else {})
        out.append(Delta(row=row, changes=changes))
        previous[row.env] = row
    return out


@dataclass(frozen=True)
class Regression:
    """One region that GREW at one commit."""

    sha: str
    date: str
    env: str
    region: str
    delta: int
    size: int
    subject: str


def find_regressions(rows: list[TrailRow], *, min_delta: int = 1,
                     region: str | None = None) -> list[Regression]:
    """Every (commit, env, region) whose size grew by >= min_delta, largest first.

    A first row for an env has no predecessor and so is never a regression:
    reporting its whole size as growth would bury every real one.
    """
    wanted = REGIONS if region is None else (region,)
    out: list[Regression] = []
    for delta in compute_deltas(rows):
        for name in wanted:
            grew = delta.changes.get(name, 0)
            if grew >= min_delta:
                out.append(Regression(
                    sha=delta.row.sha, date=delta.row.date, env=delta.row.env,
                    region=name, delta=grew, size=delta.row.sizes.get(name, 0),
                    subject=delta.row.subject))
    out.sort(key=lambda r: r.delta, reverse=True)
    return out


# ---------------------------------------------------------------------------
# Git plumbing
# ---------------------------------------------------------------------------
class GitError(RuntimeError):
    pass


def _git(args: list[str], cwd: str | Path | None = None) -> str:
    try:
        proc = subprocess.run(["git", *args], cwd=cwd, check=False,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
    except OSError as exc:
        raise GitError(f"cannot run git: {exc}") from exc
    if proc.returncode != 0:
        raise GitError(f"git {' '.join(args)}: {proc.stderr.strip()}")
    return proc.stdout.strip()


def git_common_dir(cwd: str | Path | None = None) -> Path:
    """Absolute path of the repo's COMMON git dir (shared by all worktrees)."""
    raw = Path(_git(["rev-parse", "--git-common-dir"], cwd))
    if raw.is_absolute():
        return raw
    return (Path(cwd or os.getcwd()) / raw).resolve()


def default_trail(cwd: str | Path | None = None) -> Path:
    return git_common_dir(cwd) / TRAIL_NAME


def default_pending(cwd: str | Path | None = None) -> Path:
    return git_common_dir(cwd) / PENDING_NAME


def head_stamp(cwd: str | Path | None = None, rev: str = "HEAD"
               ) -> tuple[str, str, str]:
    """(sha, committer date, subject) for a revision."""
    out = _git(["log", "-1", "--format=%H%n%cI%n%s", rev], cwd)
    sha, date, subject = (out.split("\n", 2) + ["", "", ""])[:3]
    return sha, date, subject


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------
def collect(build_dir: str | Path, envs: tuple[str, ...],
            warn=print) -> dict[str, dict[str, int]]:
    """Parse each environment's linked ELF under build_dir.

    An environment that did not build, or whose ELF will not parse, is warned
    about and skipped — a partial record beats no record.
    """
    root = Path(build_dir)
    found: dict[str, dict[str, int]] = {}
    for env in envs:
        elf = root / env / ELF_NAME
        if not elf.is_file():
            warn(f"no {elf} - skipping env '{env}'")
            continue
        try:
            found[env] = regions_from_elf(elf)
        except (ElfFormatError, OSError) as exc:
            warn(f"cannot read {elf} ({exc}) - skipping env '{env}'")
    return found


def _warn(message: str) -> None:
    print(f"[size-trail] {message}", file=sys.stderr)


def cmd_record(args) -> int:
    envs = tuple(args.env) if args.env else DEFAULT_ENVIRONMENTS
    found = collect(args.build_dir, envs, warn=_warn)
    if not found:
        _warn("no firmware ELF parsed - nothing recorded.")
        return 1
    out = Path(args.out) if args.out else default_pending()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"envs": found}, indent=1, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(f"[size-trail] recorded {', '.join(sorted(found))} -> {out}")
    return 0


def cmd_commit(args) -> int:
    pending = Path(args.pending) if args.pending else default_pending()
    if not pending.is_file():
        return 0  # no capture this commit (doc-only, gate skipped): not an error
    trail = Path(args.trail) if args.trail else default_trail()
    try:
        payload = json.loads(pending.read_text(encoding="utf-8"))
        envs = payload["envs"]
        sha, date, subject = head_stamp(args.repo, args.rev)
        rows = [TrailRow(sha=sha, date=date, env=env,
                         sizes={r: int(sizes.get(r, 0)) for r in REGIONS},
                         subject=subject)
                for env, sizes in sorted(envs.items())]
        append_rows(trail, rows)
    except (OSError, ValueError, KeyError, GitError) as exc:
        print(f"[size-trail] could not append to {trail} ({exc}).",
              file=sys.stderr)
        pending.unlink(missing_ok=True)
        return 1
    pending.unlink(missing_ok=True)
    summary = " ".join(f"{row.env} itcm={row.sizes['itcm']:,}" for row in rows)
    print(f"[size-trail] {sha[:8]} {summary}")
    return 0


def _fmt_delta(value: int) -> str:
    return f"{value:+,}" if value else "."


def render_show(rows: list[TrailRow], last: int | None = None) -> str:
    """Render the trail with per-row deltas, oldest first."""
    deltas = compute_deltas(rows)
    if last:
        deltas = deltas[-last:]
    if not deltas:
        return "(trail empty)"
    head = f"{'sha':<9} {'env':<15} " + " ".join(f"{r:>10}" for r in REGIONS)
    lines = [head, "-" * len(head)]
    for d in deltas:
        sizes = " ".join(f"{d.row.sizes.get(r, 0):>10,}" for r in REGIONS)
        lines.append(f"{d.row.sha[:8]:<9} {d.row.env:<15} {sizes}"
                     f"  {d.row.date[:10]} {d.row.subject}")
        if not d.first:
            marks = " ".join(f"{_fmt_delta(d.changes.get(r, 0)):>10}"
                             for r in REGIONS)
            lines.append(f"{'':<9} {'':<15} {marks}")
    return "\n".join(lines)


def cmd_show(args) -> int:
    trail = Path(args.trail) if args.trail else default_trail()
    rows = read_trail(trail)
    if not rows:
        print(f"[size-trail] {trail}: no rows yet.")
        return 0
    print(render_show(rows, args.last))
    return 0


def render_regressions(items: list[Regression]) -> str:
    if not items:
        return "(no growth recorded)"
    head = (f"{'delta':>10} {'now':>10}  {'region':<8} {'env':<15} "
            f"{'sha':<9} subject")
    lines = [head, "-" * len(head)]
    for r in items:
        lines.append(f"{r.delta:>+10,} {r.size:>10,}  {r.region:<8} "
                     f"{r.env:<15} {r.sha[:8]:<9} {r.subject}")
    return "\n".join(lines)


def cmd_regressions(args) -> int:
    trail = Path(args.trail) if args.trail else default_trail()
    if args.region and args.region not in REGIONS:
        print(f"[size-trail] unknown region '{args.region}' "
              f"(known: {', '.join(REGIONS)})", file=sys.stderr)
        return 2
    print(render_regressions(find_regressions(
        read_trail(trail), min_delta=args.min_delta, region=args.region)))
    return 0


def _elf_stamp(build_dir: Path, env: str) -> int | None:
    """mtime of an environment's linked ELF in ns, or None when it is absent."""
    try:
        return (build_dir / env / ELF_NAME).stat().st_mtime_ns
    except OSError:
        return None


def cmd_backfill(args) -> int:
    """Build each commit of a rev-range in a worktree and record its sizes.

    SLOW: one full three-image firmware link per commit (minutes each).

    `pio run` exits non-zero when a commit fails the size gate, but the ELF is
    LINKED BEFORE the gate runs, so pio's status cannot separate an over-budget
    commit (worth trailing) from one that never compiled. The ELF mtime can: a
    relink always restamps it — a full object-cache hit still links — while a
    failed compile never reaches the linker and leaves the previous commit's
    ELF in place. An environment whose ELF did not change is skipped, never
    recorded under this commit's sha.
    """
    worktree = Path(args.worktree)
    if not worktree.is_dir():
        print(f"[size-trail] no worktree at {worktree}", file=sys.stderr)
        return 2
    trail = Path(args.trail) if args.trail else default_trail()
    envs = tuple(args.env) if args.env else DEFAULT_ENVIRONMENTS
    build_dir = Path(args.build_dir) if args.build_dir else worktree / ".pio" / "build"

    try:
        revs = _git(["rev-list", "--reverse", args.rev_range], worktree).split()
    except GitError as exc:
        print(f"[size-trail] {exc}", file=sys.stderr)
        return 2
    have = {(r.sha, r.env) for r in read_trail(trail)}

    for rev in revs:
        sha, date, subject = head_stamp(worktree, rev)
        todo = [e for e in envs if (sha, e) not in have]
        if not todo:
            print(f"[size-trail] {sha[:8]} already trailed - skipping.")
            continue
        try:
            _git(["checkout", "--detach", "--force", rev], worktree)
        except GitError as exc:
            print(f"[size-trail] {sha[:8]}: checkout failed ({exc})",
                  file=sys.stderr)
            continue
        pio = ["pio", "run", "-d", str(worktree)]
        for env in todo:
            pio += ["-e", env]
        before = {env: _elf_stamp(build_dir, env) for env in todo}
        try:
            subprocess.run(pio, check=False)
        except OSError as exc:
            print(f"[size-trail] pio not runnable ({exc})", file=sys.stderr)
            return 2
        linked = [env for env in todo
                  if _elf_stamp(build_dir, env) != before[env]]
        for env in todo:
            if env not in linked:
                print(f"[size-trail] {sha[:8]}: '{env}' ELF not relinked - "
                      f"build failed; skipping.", file=sys.stderr)
        if not linked:
            continue
        found = collect(build_dir, tuple(linked),
                        warn=lambda m: print(f"[size-trail] {sha[:8]}: {m}",
                                             file=sys.stderr))
        if not found:
            continue
        append_rows(trail, [TrailRow(sha=sha, date=date, env=env, sizes=sizes,
                                     subject=subject)
                            for env, sizes in sorted(found.items())])
        print(f"[size-trail] {sha[:8]} recorded {', '.join(sorted(found))}")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Per-commit Teensy firmware size trail. The trail is local "
                    "only (git common dir), never committed.")
    sub = p.add_subparsers(dest="command", required=True)

    rec = sub.add_parser("record", help="parse built ELFs into a pending record")
    rec.add_argument("--build-dir", default=".pio/build")
    rec.add_argument("--out", help=f"default: <git-common-dir>/{PENDING_NAME}")
    rec.add_argument("--env", action="append",
                     help="environment to record (repeatable); default "
                          + ", ".join(DEFAULT_ENVIRONMENTS))
    rec.set_defaults(func=cmd_record)

    com = sub.add_parser("commit", help="stamp the pending record and append it")
    com.add_argument("--pending", help=f"default: <git-common-dir>/{PENDING_NAME}")
    com.add_argument("--trail", help=f"default: <git-common-dir>/{TRAIL_NAME}")
    com.add_argument("--repo", help="repo/worktree to read the commit from")
    com.add_argument("--rev", default="HEAD")
    com.set_defaults(func=cmd_commit)

    show = sub.add_parser("show", help="print the trail with per-row deltas")
    show.add_argument("--trail")
    show.add_argument("--last", type=int, help="only the last N rows")
    show.set_defaults(func=cmd_show)

    reg = sub.add_parser("regressions", help="commits where a region grew")
    reg.add_argument("--trail")
    reg.add_argument("--min-delta", type=int, default=1)
    reg.add_argument("--region", help="one of: " + ", ".join(REGIONS))
    reg.set_defaults(func=cmd_regressions)

    bf = sub.add_parser("backfill",
                        help="build each commit of a rev-range and record it "
                             "(SLOW: one firmware link per commit)")
    bf.add_argument("rev_range", metavar="rev-range")
    bf.add_argument("--worktree", required=True,
                    help="checkout to build in; it is hard-reset per commit")
    bf.add_argument("--build-dir", help="default: <worktree>/.pio/build")
    bf.add_argument("--trail")
    bf.add_argument("--env", action="append")
    bf.set_defaults(func=cmd_backfill)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except GitError as exc:
        print(f"[size-trail] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
