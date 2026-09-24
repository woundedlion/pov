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
`CI green` check. `.githooks/pre-commit` is a fast prefilter over the staged
tree — format, lint, documentation, build pins and license headers; the
protected branch's `CI green` status is the authoritative correctness gate.

- **`.githooks/pre-commit`** — rejects staged whitespace errors, checks staged
  first-party C++ with clang-format, and runs ruff/eslint on staged
  Python/JavaScript. It then runs `tools/docs_check.py` (without `--sync`),
  `tools/docs_images.py`, `tools/build_pins.py --check` and
  `tools/license_check.py` against an isolated checkout of the index, so the
  verdict is on what is being committed rather than on the working tree. A
  required tool missing for an applicable change fails the commit rather than
  skipping the check. Configuring the `tests` preset points `core.hooksPath` at
  `.githooks` for you.
- **Shell edits:** the pre-commit hook runs `shellcheck -x` on staged `.sh` files and `.githooks/` scripts. Install its pinned prerequisite with `pip install --require-hashes -r requirements/shellcheck.txt`.
- **clang-format is pinned to major 22.** A different major reflows unrelated
  code, so the hook fails rather than trusting an off-major verdict. Install the
  pin (`pip install clang-format==22.1.8`) or point `CLANG_FORMAT` at a
  `clang-format-22` binary. Every external tool version is single-sourced
  through `tools/build_pins.py`, whose `--check` fails a partial bump.
- **Native suite:** `cmake --preset tests && cmake --build --preset tests` then
  `ctest --preset tests --output-on-failure --no-tests=error`. Every CI leg
  drives `HS_SMOKE_FRAMES=120`; at the 8-frame default no preset transition
  arms. Set `HS_EFFECTS_FULL=1` to reproduce the full-resolution master leg
  locally.
- **Node script suite:** `npm test` runs every `scripts/*.test.mjs` — the
  shader-workbench schema and digest contracts, the WASM smoke predicates, the
  engine bindings contract, the profile roster and the PNG probe. CI runs it as
  its own `scripts-unit-tests` job.
- **Native variants and coverage:** `sanitizers`, `thread-sanitizer`,
  `optimized-tests` and `windows-tests` exercise distinct runtime and platform
  configurations. `code-coverage` enforces aggregate and directory coverage;
  `shard-coverage` checks test-module selection across CI legs.
- **Generated-source provenance:** `lut-provenance`, `reaction-graph-provenance`,
  `gamut-lut-provenance`, `srgb-decode-provenance` and `patterns-provenance`
  regenerate their artifacts and compare them with committed bytes. Change the
  generator or authored input and regenerate; hand-editing its output fails.
- **Firmware:** `teensy-gate-tests` tests the gate tooling, `teensy-size` enforces
  firmware memory budgets, and `teensy-warnings` checks every PlatformIO
  environment with the pinned toolchain. The local entry points are
  `just python-test`, `just teensy-size` and `just teensy-warnings`.
- **Published artifacts:** `wasm` builds and verifies the engine bundle;
  `screenshot-gallery` checks capture membership and images. `docs-doxygen`
  builds the API reference with warnings treated as errors. These complement
  the Markdown and image-reference checks below.
- **Host-Python tool suites:** `just python-test` and CI run
  `python tools/run_python_tests.py`. It discovers every tracked suite, rejects
  empty suites, and propagates failures. Install `requirements/numpy.txt` first;
  no ARM toolchain or KiCad is required.
- **Lint:** the CI `lint` job has seven legs — `tools/eol_gate.sh` first, so
  every later leg reads the line endings `.gitattributes` declares, then
  `ruff` over the Python tooling, `eslint` over the JavaScript, `shellcheck`
  over every tracked `*.sh` and `.githooks/*`, `actionlint` over
  `.github/workflows/*.yml` (which pipes every `run:` body through
  `shellcheck`, since no workflow is a `*.sh` file), a `just --evaluate` /
  `just --summary` parse of the `justfile`, and the profiling-roster
  cross-check in `tools/profile_sweep.sh`. `just lint` runs five of them
  locally (line endings, `ruff`, `eslint`, `shellcheck` and the roster
  check); the hook lints only staged Python and JavaScript, so CI remains
  authoritative.
- **Documentation:** the ci.yml docs-markdown job runs `tools/docs_check.py`
  without `--sync`: fences, links, anchors, every backticked repo path, the
  README's file map against the tracked tree and its effect counts against
  `HS_EFFECT_LIST`. `just docs-sync` regenerates the maps and counts first,
  then runs the same checker, so the repaired diff lands with the change.
  `python tools/docs_images.py` resolves every documented `<img>` against the
  tracked tree. It only reports; `--stage` copies the images into a built
  Doxygen tree and is the sole mode that writes.
- **License headers:** `python tools/license_check.py`. The staged-tree
  pre-commit check and `just license-headers` run it; `just python-test` runs
  the checker's unit tests.
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
