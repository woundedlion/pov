# Contributing

This repository holds the Holosphere engine and firmware; the browser simulator
lives in the sibling [daydream](https://github.com/woundedlion/daydream)
repository and is built and installed from here. Read `README.md` §3–4 first —
its file map is gated against the tracked tree, and §4 describes the architecture.
[§11 Building](README.md#11-building) has the commands for all three
targets: the Teensy firmware, the WASM module, and the native test suite.

## Licensing

`LICENSE` grants PolyForm Noncommercial 1.0.0 to everything outside `effects/`
except four named paths: `workbench/` and `core/engine/effects_legacy.h` are
reserved like `effects/`; `core/math/projections.h` and
`core/vendor/FastNoiseLite.h` are MIT. All rights over `effects/` are reserved. Every
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
component (`scan: clamp the row index before the cast`).

## Design specs

[The specifications index](docs/specs/README.md) lists the design contracts
that span several files. A change that moves such a contract carries the
corresponding spec update with it.

A change that stays inside one file needs no spec. A new subsystem other code
will be written against does: land the spec with the implementation, not after
it.

## Gates

Host tooling and Git hooks require Python 3.11 or newer. Set `HS_PYTHON` to
select a supported interpreter; CI uses the version pinned in `tools/build_pins.py`.

The engine CI gates run in `.github/workflows/ci.yml` behind the aggregate
`CI green` check. Local hooks and the sibling simulator checks below have
separate entry points. `.githooks/pre-commit` is a fast prefilter over the staged
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
  Changes to the profile archive, its validator, or its effect roster inputs
  also run the archive structure and coverage gate on that staged snapshot,
  using the pinned Node version.
- **Shell edits:** the pre-commit hook runs `shellcheck -x` on staged `.sh` files and `.githooks/` scripts. Install its pinned prerequisite with `pip install --require-hashes -r requirements/shellcheck.txt`.
- **clang-format is pinned to 22.1.8.** A different major reflows unrelated
  code, so the hook fails rather than trusting an off-major verdict. Install the
  pin (`pip install clang-format==22.1.8`) or point `CLANG_FORMAT` at a
  `clang-format-22` binary reporting that exact version. Every external tool version is single-sourced
  through `tools/build_pins.py`, whose `--check` fails a partial bump.
- **Guard tests:** Prove that a valid companion fixture is accepted. For the
  guarded fixture, assert the intended result or diagnostic. Deleting the
  guarded line must fail the test; an unrelated rejection cannot satisfy it.
- **Native suite:** `cmake --preset tests && cmake --build --preset tests` then
  `ctest --preset tests --output-on-failure --no-tests=error`. Every CI leg
  drives `HS_SMOKE_FRAMES=120`; at the 8-frame default no preset transition
  arms. Set `HS_EFFECTS_FULL=1` to reproduce the full-resolution leg (also run on pull requests)
  locally.
- **Node script suite:** `npm test` runs every `scripts/*.test.mjs` — the
  shader-workbench schema and digest contracts, the WASM smoke predicates, the
  engine bindings contract, the profile roster and the PNG probe. CI runs it as
  its own `scripts-unit-tests` job. Set `HS_BROWSER_TESTS=1` when running
  `node --test scripts/screenshot_resolution.test.mjs` to include the Chromium
  capture-descent test; first run `npm ci` and `npx playwright install chromium`.
- **Native variants and coverage:** `sanitizers`, `thread-sanitizer`,
  `optimized-tests` and `windows-tests` exercise distinct runtime and platform
  configurations. `optimized-tests` also runs the death harness with `NDEBUG`.
  `code-coverage` enforces aggregate and directory coverage.
- **Generated-source provenance:** `lut-provenance`, `reaction-graph-provenance`,
  `gamut-lut-provenance`, `srgb-decode-provenance` and `patterns-provenance`
  regenerate their artifacts and compare them with committed bytes. Change the
  generator or authored input and regenerate; hand-editing its output fails.
- **Firmware:** `teensy-size` enforces firmware memory and layout budgets, and
  `teensy-warnings` checks every PlatformIO environment with the pinned toolchain.
  The local entry points are `just teensy-size` and `just teensy-warnings`.
- **Published artifacts:** `wasm` builds and verifies the engine bundle;
  `screenshot-gallery` checks capture membership and images. `docs-doxygen`
  builds the API reference with warnings treated as errors. These complement
  the Markdown and image-reference checks below.
- **Host-Python tool suites:** `just python-test` and the CI `python-tests`
  job run `python tools/run_python_tests.py` across all tracked Python suites,
  including firmware gates, profiling tools, build checks and PCB generators.
  It rejects empty suites and propagates failures. The CI job also checks
  routed PCB metadata with `hardware/phantasm/gen/board_metadata.py --check`. Install `requirements/numpy.txt`, `requirements/ruff.txt` and
`requirements/shellcheck.txt` first (the same set the CI job installs);
  no ARM toolchain or KiCad is required.
- **PCB:** `pcb-tests` runs the KiCad-backed generation, DRC, parity and
  fabrication suite (`python tools/run_python_tests.py --suite hardware/phantasm/gen/tests`)
  against pinned KiCad 10.0.4; locally it needs that KiCad install.
- **Lint:** `just lint` checks `just --fmt --check`, working-tree whitespace,
  declared line endings, Python selection and ruff, JavaScript selection and
  eslint, tracked shell files with shellcheck, the profiling roster, and
  actionlint last. CI applies the
  corresponding checks; the hook lints staged Python, JavaScript, and shell files.
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
  `pre-push` hook runs lint, typecheck, the import-map and Tailwind stylesheet freshness
  checks, and three workflow helper tests. The full JavaScript suite is separate. daydream's deployment gate
  tests an immutable Holosphere/daydream pair before publishing; Holosphere CI
  does not run a daydream consumer job. See daydream's deployment documentation.

`--no-verify` is the explicit emergency escape from the local prefilter. It
does not bypass protected-branch CI.

## Reporting a vulnerability

Report privately through GitHub's security advisories for this repository
("Security" → "Report a vulnerability") rather than opening a public issue.
