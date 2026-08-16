# Contributing

This repository holds the Holosphere engine and firmware; the browser simulator
lives in the sibling [daydream](https://github.com/woundedlion/daydream)
repository and is built and installed from here. Read `README.md` §1–2 first —
it is the architecture reference, and its file map is gated against the tracked
tree.

## Licensing

`LICENSE` grants PolyForm Noncommercial 1.0.0 to everything outside `effects/`
and reserves all rights over `effects/`, with three named exceptions. Every
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

## Gates

Every gate below except the simulator suite runs in `.github/workflows/ci.yml`
behind one aggregate `CI green` check; the simulator suite runs in daydream's
own workflows, and no check here covers it. `.githooks/pre-commit` runs the
format, lint/docs, native-suite and Teensy gates locally on every commit; the
license-header check and the simulator suite have no local hook and are run by
hand.

- **`.githooks/pre-commit`** — four gates, each keyed on the staged paths: a
  `clang-format` check over the staged first-party sources, the Python /
  JavaScript / Markdown lint and docs checks, a build and run of the native
  suite, and the Teensy size/layout gate. Configuring the `tests` preset points
  `core.hooksPath` at `.githooks` for you.
- **clang-format is pinned to major 22.** A different major reflows unrelated
  code, so the hook fails rather than trusting an off-major verdict. Install the
  pin (`pip install clang-format==22.1.8`) or point `CLANG_FORMAT` at a
  `clang-format-22` binary. Every external tool version is single-sourced
  through `tools/build_pins.py`, whose `--check` fails a partial bump.
- **Native suite:** `cmake --preset tests && cmake --build --preset tests` then
  `ctest --preset tests --output-on-failure --no-tests=error`. `HS_EFFECTS_FULL=1`
  selects the full-resolution effect tier CI runs; a green hook alone is the
  QUICK tier.
- **Documentation:** `python tools/docs_check.py` validates fences, links,
  anchors and every backticked repo path, and the README's file map must list a
  new tracked path; the pre-commit hook and `just docs-check` both run it.
  `python tools/docs_images.py` resolves every documented `<img>` against the
  tracked tree. It only reports; `--stage` copies the images into a built
  Doxygen tree and is the sole mode that writes.
- **License headers:** `python tools/license_check.py`. No hook runs it; `just
  license-headers` does, alongside the checker's own unit tests — otherwise the
  first evidence of a missing header is a red `License headers match LICENSE`
  job.
- **Simulator:** in the daydream checkout, `npm ci` then `npm test`; its
  `pre-push` hook runs lint, typecheck, the import-map check and the JS suite,
  and refuses a push from a tree that cannot run them. This repository's CI
  never runs it, so a green `CI green` says nothing about the simulator.

`HS_SKIP_TESTS=1` stands the native suite, the lint checks and the Teensy gate
down for one commit; `HS_SKIP_LINT=1` and `HS_SKIP_TEENSY=1` stand down one
apiece. The format check runs ahead of all three and stands down only for
`HS_SKIP_FORMAT=1` or `--no-verify`, which disables the whole hook — an
unformatted commit reds CI for whoever pushes next.

## Reporting a vulnerability

Report privately through GitHub's security advisories for this repository
("Security" → "Report a vulnerability") rather than opening a public issue.
