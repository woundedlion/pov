#!/usr/bin/env python3
# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Check literal death-pin uniqueness and cross-directory breadcrumb aliases,
and reject raw assert() in guard sources.

With --new-guards BASE TIP, instead fail when TIP adds an always-on guard whose
text BASE contains nowhere and no death case pins.
"""

import argparse
import ast
from pathlib import Path
import re
import subprocess

STRING = r'"(?:[^"\\]|\\.)*"'
CHAR = r"'(?:[^'\\\n]|\\.)+'"
TOKEN = re.compile(STRING + r"|//[^\n]*|/\*.*?\*/|[(),]|[^(),\"/]+|/", re.S)


def literals(text):
    """Decode adjacent C++ string literals."""
    return "".join(ast.literal_eval(s) for s in re.findall(STRING, text))


def pins(text):
    """Read the literal fields of the death-case table."""
    pattern = (r'\{\s*"([^"\n]+)"\s*,\s*(\w+)\s*,\s*'
               r'"([^"\n]+)"\s*,\s*((?:' + STRING + r'\s*)+)'
               r'(?:,\s*(?:true|false)\s*)?'
               r'(?:,\s*(?:' + STRING + r'|nullptr)\s*)?\}')
    rows = []
    for match in re.finditer(pattern, text):
        name, function, file, message = match.groups()
        if function != "case_" + name:
            raise ValueError(f"{name}: expected case_{name}, found {function}")
        rows.append((name, file, literals(message)))
    heads = re.findall(r'\{\s*"[^"\n]+"\s*,\s*\w+\s*,', text)
    if len(rows) != len(heads):
        raise ValueError(f"parsed {len(rows)} of {len(heads)} death-case rows")
    return rows


def guard_texts(text):
    """Read macro condition/message pairs with balanced parentheses."""
    for _, guard in guard_sites(text):
        yield guard


def guard_sites(text):
    """Yield (offset into text, condition/message pair) for each guard."""
    # Blank comments in place so match offsets still index the original text.
    text = re.sub(CHAR + r"|" + STRING + r"|//[^\n]*|/\*.*?\*/",
                  lambda m: m[0] if m[0][0] in "\"'" else re.sub(r"[^\n]", " ", m[0]),
                  text, flags=re.S)
    for match in re.finditer(r"HS_(?:AUDIT_)?CHECK\s*\(", text):
        depth = 1
        condition = []
        message = []
        field = condition
        for token in TOKEN.finditer(text, match.end()):
            value = token[0]
            if value == "(":
                depth += 1
            elif value == ")":
                depth -= 1
                if depth == 0:
                    break
            elif value == "," and depth == 1:
                if field is message:
                    break
                field = message
                continue
            field.append(value)
        yield match.start(), "(" + "".join(condition).strip() + ") " + literals("".join(message))
    direct = r"(?:hs::)?check_fail\s*\(\s*HS_CHECK_SITE\s*\(\s*(" + STRING + r")\s*\)\s*,\s*((?:" + STRING + r"\s*)+)"
    for match in re.finditer(direct, text):
        yield match.start(), "(" + literals(match[1]) + ") " + literals(match[2])


RAW_ASSERT = re.compile(r"(?<!\w)assert\s*\(|#\s*include\s*<(?:cassert|assert\.h)>")


def raw_assert_lines(text):
    """Line numbers of raw assert() calls and <cassert> includes outside comments and strings."""
    code = re.sub(CHAR + r"|" + STRING + r"|//[^\n]*|/\*.*?\*/", lambda m: re.sub(r"[^\n]", " ", m[0]), text,
                  flags=re.S)
    return [code.count("\n", 0, m.start()) + 1 for m in RAW_ASSERT.finditer(code)]


def compact(text):
    return re.sub(r"\s+", "", text)


def guard_matches(guard, pin):
    """Match a literal breadcrumb prefix against a guard's printf message."""
    guard, pin = compact(guard), compact(pin)
    depth = 0
    message_start = 0
    for index, token in enumerate(guard):
        if token == "(":
            depth += 1
        elif token == ")":
            depth -= 1
            if depth == 0:
                message_start = index + 1
                break
    formats = list(re.compile(r"%(?:[-+ #0]*)(?:\d+|\*)?(?:\.(?:\d+|\*))?(?:hh|ll|[hljztL])?[diuoxXfFeEgGaAcsp]").finditer(guard, message_start))
    if not formats:
        return guard.startswith(pin)
    if guard[:formats[0].start()].startswith(pin):
        return True
    pieces = []
    start = 0
    for field in formats:
        pieces.append(re.escape(guard[start:field.start()]))
        pieces.append(".*?")
        start = field.end()
    pieces.append(re.escape(guard[start:]))
    return re.match("".join(pieces), pin) is not None


GUARD_RENAME_PINS = (
    ('core/animation/opleg/sweep.h', '(vertex < V) OpLeg: dual face index out of range', '(vertex < v) OpLeg: dual face index out of range'),
    ('core/animation/opleg/sweep.h', '(out == V) OpLeg: dual corners incomplete', '(out == v) OpLeg: dual corners incomplete'),
    ('core/mesh/conway.h', '(static_cast<size_t>(faces[offset + k]) < V) kis: face vertex index out of range', '(static_cast<size_t>(faces[offset + k]) < v) kis: face vertex index out of range'),
    ('core/mesh/conway.h', '(V <= UINT16_MAX) reconcile_vertices: endpoint exceeds index capacity', '(v <= UINT16_MAX) reconcile_vertices: endpoint exceeds index capacity'),
    ('core/mesh/hankin.h', '((I / 2) + I <= MeshLimits::MAX_VERTICES) Hankin output vertex count exceeds int16_t index range (oversized source mesh?)', '((i_value / 2) + i_value <= MeshLimits::MAX_VERTICES) Hankin output vertex count exceeds int16_t index range (oversized source mesh?)'),
    ('core/mesh/hankin.h', '(count + 2 <= (int)(2 * I)) Hankin rosette winding overflow', '(count + 2 <= (int)(2 * i_value)) Hankin rosette winding overflow'),
    ('core/render/sdf/face_class_bake.h', '(mesh.get_topology_size() == F) build_mesh_class_bake requires classify_faces_by_topology first', '(mesh.get_topology_size() == f_value) build_mesh_class_bake requires classify_faces_by_topology first'),
    ('core/render/sdf/face_class_bake.h', '(mesh.get_face_offsets_size() == F) build_mesh_class_bake requires one face offset per face', '(mesh.get_face_offsets_size() == f_value) build_mesh_class_bake requires one face offset per face'),
    ('core/render/sdf/face_class_bake.h', '(static_cast<size_t>(idx[k]) < V) mesh face vertex index out of range', '(static_cast<size_t>(idx[k]) < v_value) mesh face vertex index out of range'),
    ('targets/wasm/engine_bindings.h', '(setEffect(NAME) == EffectSetResult::INSTALLED) Geometry rebuild must reinstall the current effect', '(setEffect(name) == EffectSetResult::INSTALLED) Geometry rebuild must reinstall the current effect'),
    ('targets/wasm/engine_bindings.h', '(current_effect->synchronizePreset(PRESET)) Geometry rebuild must restore the selected preset', '(current_effect->synchronizePreset(preset)) Geometry rebuild must restore the selected preset'),
)


SOURCE_SUFFIXES = (".h", ".cpp", ".ino")
WAIVER = re.compile(r"//[ \t]*no death case:[ \t]*\S")


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True, encoding="utf-8").stdout


def blobs(repo, commit, paths):
    """Read paths as of commit in one batch; missing paths are omitted."""
    if not paths:
        return {}
    request = "".join(f"{commit}:{path}\n" for path in paths).encode("utf-8")
    out = subprocess.run(["git", "-C", str(repo), "cat-file", "--batch"],
                         input=request, check=True, capture_output=True).stdout
    found = {}
    offset = 0
    for path in paths:
        end = out.index(b"\n", offset)
        header = out[offset:end].split()
        offset = end + 1
        if header[-1] == b"missing":
            continue
        size = int(header[2])
        found[path] = out[offset:offset + size].decode("utf-8")
        offset += size + 1
    return found


def waived(text, offset):
    """True when a `// no death case: <reason>` comment sits on the guard's line or the line above."""
    start = text.rfind("\n", 0, offset) + 1
    above = text.rfind("\n", 0, max(start - 1, 0)) + 1
    end = text.find("\n", offset)
    return WAIVER.search(text[above:end if end >= 0 else len(text)]) is not None


def renamed_guards(before, after, file):
    """Match explicit identifier-only census pins at unchanged guard positions."""
    old_guards = list(guard_texts(before))
    new_guards = list(guard_texts(after))
    if len(old_guards) != len(new_guards):
        return set()
    aliases = {(compact(old), compact(new)) for path, old, new in GUARD_RENAME_PINS if path == file}
    matched = set()
    token = re.compile(STRING + "|" + CHAR + r"|(?:0[xX][0-9a-fA-F]+|[0-9]+)(?:\.[0-9]*)?(?:[eEpP][+-]?[0-9]+)?[uUlLfF]*|[A-Za-z_][A-Za-z_0-9]*|[^\s]")
    identifier = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")
    for old, new in zip(old_guards, new_guards):
        if (compact(old), compact(new)) not in aliases:
            continue
        def parts(guard):
            depth = 0
            for match in TOKEN.finditer(guard):
                if match[0] == "(":
                    depth += 1
                elif match[0] == ")":
                    depth -= 1
                    if depth == 0:
                        return token.findall(guard[:match.end()]), guard[match.end():]
            return [], guard
        old_tokens, old_message = parts(old)
        new_tokens, new_message = parts(new)
        if old_message != new_message or len(old_tokens) != len(new_tokens):
            continue
        changes = [(a, b) for a, b in zip(old_tokens, new_tokens) if a != b]
        mapping = dict(changes)
        if (changes and all(identifier.fullmatch(a) and identifier.fullmatch(b) for a, b in changes)
                and len(set(mapping.values())) == len(mapping)
                and [mapping.get(a, a) for a in old_tokens] == new_tokens):
            matched.add(compact(new))
    return matched


def new_guards(repo, base, tip):
    """List "path:line: guard" for each guard tip adds over base that no case pins."""
    directories = blobs(repo, tip, ["tests/guard_directories.txt"])[
        "tests/guard_directories.txt"].split()
    changed = [path for path in git(repo, "diff", "--name-only", "--no-renames", "--diff-filter=AM",
                                    base, tip, "--", *directories).splitlines()
               if path.endswith(SOURCE_SUFFIXES)]
    if not changed:
        return []
    bearing = [hit.split(":", 1)[1] for hit in git(
        repo, "grep", "-l", "-E", r"HS_(AUDIT_)?CHECK|check_fail", base, "--",
        *directories).splitlines()]
    base_sources = blobs(repo, base, [p for p in bearing if p.endswith(SOURCE_SUFFIXES)])
    existing = {compact(guard)
                for text in base_sources.values()
                for guard in guard_texts(text)}
    case_pins = pins(blobs(repo, tip, ["tests/test_death.h"])["tests/test_death.h"])
    missing = []
    for path, text in blobs(repo, tip, changed).items():
        renames = renamed_guards(base_sources.get(path, ""), text, path)
        for offset, guard in guard_sites(text):
            if text.startswith("HS_AUDIT_CHECK", offset) or compact(guard) in existing or compact(guard) in renames or waived(text, offset):
                continue
            if any(file == path and guard_matches(guard, pin) for _, file, pin in case_pins):
                continue
            missing.append(f"{path}:{text.count(chr(10), 0, offset) + 1}: {guard}")
    return missing


def check_new_guards(repo, base, tip):
    missing = new_guards(repo, base, tip)
    for site in missing:
        print(f"new guard without a death case: {site}")
    if missing:
        raise SystemExit(f"{len(missing)} new guard(s) need a death case in tests/test_death.h "
                         "or a `// no death case: <reason>` comment")
    print("death pins: every new guard has a death case")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--new-guards", nargs=2, metavar=("BASE", "TIP"))
    args = parser.parse_args()
    if args.new_guards:
        check_new_guards(args.root, *args.new_guards)
        return
    root = args.root
    case_pins = pins((root / "tests/test_death.h").read_text(encoding="utf-8"))
    pinned_names = {Path(file).name for _, file, _ in case_pins}
    sources = {}
    asserts = []
    directories = (root / "tests/guard_directories.txt").read_text(encoding="utf-8").splitlines()
    for directory in directories:
        for path in (root / directory).rglob("*"):
            if not path.is_file() or not path.name.endswith(SOURCE_SUFFIXES):
                continue
            text = path.read_text(encoding="utf-8")
            name = path.relative_to(root).as_posix()
            asserts += [f"{name}:{line}" for line in raw_assert_lines(text)]
            if path.name in pinned_names:
                sources[name] = list(guard_texts(text))
    if asserts:
        raise SystemExit("raw assert vanishes under NDEBUG on device; use HS_CHECK, or HS_AUDIT_CHECK "
                         "on hot paths:\n  " + "\n  ".join(asserts))
    for name, file, text in case_pins:
        for other, candidates in sources.items():
            if other == file or Path(other).name != Path(file).name:
                continue
            if any(guard_matches(guard, text) for guard in candidates):
                raise SystemExit(f"{name}: breadcrumb for {file} also matches {other}")
    checked = 0
    for name, file, text in case_pins:
        matches = sum(guard_matches(guard, text) for guard in sources.get(file, []))
        if matches != 1:
            raise SystemExit(f"{name}: pin identifies {matches} guard sites in {file}, expected 1")
        checked += 1
    if checked == 0:
        raise SystemExit("no literal death pins found")
    print(f"death pins: {checked} literal cases each identify one guard site")


if __name__ == "__main__":
    main()
