# Contributing

This repository holds the Holosphere engine and firmware; the browser simulator
lives in the sibling [daydream](https://github.com/woundedlion/daydream)
repository and is built and installed from here. Read `README.md` §1–2 first —
it is the architecture reference, and its file map is gated against the tracked
tree. [§11 Building](README.md#11-building) has the commands for all three
targets: the Teensy firmware, the WASM module, and the native test suite.

## Licensing

`LICENSE` grants PolyForm Noncommercial 1.0.0 to everything outside `effects/`
and reserves all rights over `effects/`, with four named exceptions. Every
tracked C/C++ source carries the header its path is granted, and
`tools/license_check.py` gates that in CI. A new source file needs the header
before it lands; a new generated one needs its generator to emit it.

## Landing model

`master` is **fast-forward only**. `.githooks/reference-transaction` refuses any
non-fast-forward move of the ref, so a rewind needs a deliberate one-shot token
(the hook's header documents it). The workflow this enforces:

1. Branch in a worktree taken off the ref, never off a possibly-stale `HEAD`:
   `git worktree add <path> -b <branch> $(git rev-parse refs/heads/master)`.
2. Commit there. Never commit in someone else's tree, and never `cd` into a
   tree — pass `git -C <path>` so a stray working directory cannot commit in the
   wrong one.
3. Land by rebasing onto the live tip and merging fast-forward only:
   `git -C <main> merge --ff-only <branch>`. A refusal means a peer moved the ref
   or the main tree carries overlapping work; re-read the tip and rebase again.
   Never stash or discard another session's work to make a landing fit.

One logical change is one commit, with an imperative subject naming the
component (`scan: clamp the row index before the cast`). Commit messages carry
**no `Co-Authored-By` line**.

Agent sessions that commit here work under the additional ground rules in
`docs/agent_workflow.md`, which bind worktree discipline, the host-global lock
on the single shared Teensy, comment style, and the per-commit gates.

## Design specs

`docs/specs/` holds the design specifications — the pullback pipeline, its
stage families and preview interpreter, the shader workbench's chain schema and
editor, Phantasm's frame-sync protocol, and the Phantasm segment board. Each
is the source of truth for one contract that spans several files, so a change
that moves such a contract carries the spec update with it.

A change that stays inside one file needs no spec. A new subsystem other code
will be written against does: land the spec with the implementation, not after
it.

## Gates

Every gate below runs in `.github/workflows/ci.yml` behind one aggregate
`CI green` check. `.githooks/pre-commit` is a fast staged-file format/lint
prefilter; the protected branch's `CI green` status is the authoritative
correctness gate.

- **`.githooks/pre-commit`** — checks staged first-party C++ with clang-format
  and runs ruff/eslint on staged Python/JavaScript. A required tool missing for
  an applicable change fails the commit rather than skipping the check.
  Configuring the `tests` preset points `core.hooksPath` at `.githooks` for you.
- **clang-format is pinned to major 22.** A different major reflows unrelated
  code, so the hook fails rather than trusting an off-major verdict. Install the
  pin (`pip install clang-format==22.1.8`) or point `CLANG_FORMAT` at a
  `clang-format-22` binary. Every external tool version is single-sourced
  through `tools/build_pins.py`, whose `--check` fails a partial bump.
- **Native suite:** `cmake --preset tests && cmake --build --preset tests` then
  `ctest --preset tests --output-on-failure --no-tests=error`. Set
  `HS_EFFECTS_FULL=1` to reproduce the full-resolution master leg locally.
- **Node script suite:** `npm test` runs every `scripts/*.test.mjs` — the
  shader-workbench schema and digest contracts, the WASM smoke predicates, the
  engine bindings contract, the profile roster and the PNG probe. CI runs it as
  its own `scripts-unit-tests` job.
- **Lint:** the CI `lint` job has five legs — `ruff` over the Python tooling,
  `eslint` over the JavaScript, `shellcheck` over every tracked `*.sh` and
  `.githooks/*`, and a `just --evaluate` / `just --summary` parse of the
  `justfile`, plus the profiling-roster cross-check in `tools/profile_sweep.sh`.
  `just lint` runs the four lint/roster checks locally; the hook only checks
  staged Python and JavaScript, so CI remains authoritative.
- **Documentation:** `python tools/docs_check.py` validates fences, links,
  anchors and every backticked repo path, and the README's file map must list a
  new tracked path; `just docs-check` runs it locally.
  `python tools/docs_images.py` resolves every documented `<img>` against the
  tracked tree. It only reports; `--stage` copies the images into a built
  Doxygen tree and is the sole mode that writes.
- **License headers:** `python tools/license_check.py`. The staged-tree
  pre-commit check runs it; `just license-headers` also runs the checker's own
  unit tests.
- **Simulator:** in the daydream checkout, `npm ci` then `npm test`; its
  `pre-push` hook runs lint, typecheck, the import-map check and the JS suite,
  and refuses a push from a tree that cannot run them. The `daydream-consumer`
  job runs the same suite here: it installs the verified bundle the `wasm` job
  built over a daydream checkout pinned in `tools/build_pins.py`, so a change to
  anything this repository installs there fails before it is mirrored. daydream
  pins this repository the other way round in `holosphere_wasm.sha`, which makes
  the pair circular — land the daydream side first, then move the pin here.

`--no-verify` is the explicit emergency escape from the local prefilter. It
does not bypass protected-branch CI.

## Reporting a vulnerability

Report privately through GitHub's security advisories for this repository
("Security" → "Report a vulnerability") rather than opening a public issue.
