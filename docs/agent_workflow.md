# Agent session workflow (non-negotiable)

Ground rules for any agent session that commits to this repository. The design
briefs and ledgers in `docs/` assume them.

---

## 1. Never write to the main tree

`C:\work\Holosphere` is a shared main tree, with concurrent peer sessions
committing under the *same git identity* as you. Work in your own worktree:

```
cd /c/work/Holosphere && BASE=$(git rev-parse refs/heads/master)
git worktree add ../hs-wt-<yourtask> -b work/<yourtask> $BASE
```

Always `git -C C:/work/hs-wt-<yourtask> ...`; never `cd` into a tree and rely on
cwd. Land only by rebasing onto live `refs/heads/master` and doing an **attached
`git merge --ff-only`** in the main tree, asserting
`git merge-base --is-ancestor` first. `master` is append-only and a
`reference-transaction` hook enforces it.

If a merge refuses over main-tree WIP, **stop and surface it** — never stash or
discard a peer's work. Peers leave WIP in the main tree routinely; leave it
alone. A peer may also land a build-breaking WIP on `master`: if the tree will
not build phantasm, stop and surface it rather than building on top of it.

New files need `git add` before `git commit <path>` — a bare `git commit <path>`
silently matches nothing for an untracked file.

---

## 2. The device is a single shared Teensy

Never run `pio run -t upload` or `profile_capture.py` directly. They bypass the
host-global lock, and an upload issued while a peer holds the port reports
SUCCESS *without flashing*, so you can capture a clean, plausible log of the
peer's firmware. The only supported path is `tools/profile_one.sh`, which takes
the lock. `HS_DEVICE_WAIT=<sec>` queues; `bash tools/device_lock.sh status`
checks it. Toggle a diagnostic through `-D` flags forwarded to the build.

`tools/profile_one.sh` always builds **`/c/work/Holosphere`** (hardcoded `cd`),
whatever tree you invoke it from. **You cannot profile a worktree — land
first.** For a change that is bit-exact by construction this is safe; for
anything else, prove correctness before landing, then profile, then revert if
the device disagrees.

---

## 3. Style

Invoke the `code-style` skill first and obey it. Terse factual comments only —
no narration, no justifying correct code, no history, no finding/ticket
references. `core/render/sdf.h`, `core/color/color.h`, and generated tables
carry whole-file clang-format drift against local clang-format v22: hand-format
your own lines and wrap generated tables in `// clang-format off/on`; **never**
run `clang-format -i` on an existing file. The pre-commit format gate is live
and blocks unconditionally — a missing binary, or one whose major is not 22,
fails the commit rather than waving it through, so never bypass it with
`HS_SKIP_FORMAT=1` or `--no-verify`. No `Co-Authored-By` line.

---

## 4. Gates after every commit

- `export EMSDK=C:/work/emsdk; cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF;
  cmake --build --preset tests -j 8; ctest --preset tests` → **every registered
  test passes, zero failures** (report the count you saw; it grows as tests
  land, so a fixed number is not the criterion)
- `pio run -e phantasm` → `[teensy-gate] phantasm: PASS`, and report RAM1
  `code`, RAM1 `variables` (DTCM), FLASH `data`, **and the per-commit delta of
  each**.

RAM1 `code` (ITCM) has a 196,608 B ceiling and only a few KB of headroom. It is
a real constraint that has vetoed changes before.
