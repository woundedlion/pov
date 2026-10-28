#!/bin/sh
# Format, lint and documentation checks that predict the fast, deterministic CI
# gates, run on the index (the pre-commit hook) or on a commit range.
#
# usage: commit_checks.sh                 the index against HEAD
#        commit_checks.sh <base>..<tip>   <tip>'s tree against its merge base
#                                         with <base>: the range's net change
set -eu

SCRIPT_DIR="$(CDPATH='' cd "$(dirname "$0")" && pwd)"
# shellcheck source=tools/find_python.sh
. "$SCRIPT_DIR/find_python.sh"

LABEL=pre-commit
SUBJECT=staged
DIFF_BASE=
RANGE_DIR=
STAGED_BLOB=
SNAPSHOT_ROOT=
cleanup() {
  rm -rf "$RANGE_DIR" "$SNAPSHOT_ROOT" "$STAGED_BLOB"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

usage() {
  echo "usage: commit_checks.sh [<base>..<tip>]" >&2
  exit 2
}

cd "$(git rev-parse --show-toplevel)"

case $# in
  0) ;;
  1)
    case $1 in
      *...*) usage ;;
      ?*..?*) ;;
      *) usage ;;
    esac
    LABEL=commit-checks
    SUBJECT="range $1"
    if ! base="$(git rev-parse --verify --quiet "${1%%..*}^{commit}")" ||
       ! tip="$(git rev-parse --verify --quiet "${1#*..}^{commit}")"; then
      echo "[$LABEL] cannot resolve $1 to two commits."
      exit 2
    fi
    if ! DIFF_BASE="$(git merge-base "$base" "$tip")"; then
      echo "[$LABEL] $1 has no merge base."
      exit 2
    fi
    RANGE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/hs-commit-checks.XXXXXX")"
    GIT_INDEX_FILE="$RANGE_DIR/index"
    export GIT_INDEX_FILE
    if ! git read-tree "$tip"; then
      echo "[$LABEL] cannot read the tree of $tip."
      exit 1
    fi
    ;;
  *) usage ;;
esac

CHANGED="$(git -c core.quotePath=false diff --cached --name-only --no-renames --diff-filter=ACMRD ${DIFF_BASE:+"$DIFF_BASE"})"
PRESENT="$(git -c core.quotePath=false diff --cached --name-only --diff-filter=ACMR ${DIFF_BASE:+"$DIFF_BASE"})"
if [ -z "$CHANGED" ]; then
  [ -z "$DIFF_BASE" ] || echo "[$LABEL] $SUBJECT changes no files."
  exit 0
fi

STAGED_BLOB="$(mktemp "${TMPDIR:-/tmp}/hs-staged.XXXXXX")"

# Lints read the staged blob; a failed read must fail the hook, not lint as
# empty input.
extract_staged() {
  if ! GIT_INDEX_FILE="$INDEX_PATH" git show ":$1" >"$STAGED_BLOB"; then
    echo "[$LABEL] cannot read staged $1 from the index."
    exit 1
  fi
}

if ! git -c core.whitespace=blank-at-eol,blank-at-eof diff --cached --check ${DIFF_BASE:+"$DIFF_BASE"}; then
  echo "[$LABEL] $SUBJECT whitespace errors must be fixed."
  exit 1
fi

PYTHON_BIN="$(hs_find_python)" || PYTHON_BIN=
if [ -z "$PYTHON_BIN" ]; then
  echo "[$LABEL] no working Python found; documentation checks need one (set HS_PYTHON)."
  exit 1
fi
if ! "$PYTHON_BIN" -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
  echo "[$LABEL] Python 3.11 or newer is required; set HS_PYTHON to a supported interpreter."
  exit 1
fi

SNAPSHOT_ROOT="$(mktemp -d "$PWD/.hs-pre-commit.XXXXXX")"
SNAPSHOT="$SNAPSHOT_ROOT/repo"
mkdir "$SNAPSHOT"
if ! git checkout-index --all --force --prefix="$SNAPSHOT/"; then
  echo "[$LABEL] cannot materialize staged snapshot."
  exit 1
fi
LOCAL_GIT_ENV="$(git rev-parse --local-env-vars)"
INDEX_PATH="$(git rev-parse --path-format=absolute --git-path index)"
# Fold a split index's shared index in while the source git dir resolves it.
SNAPSHOT_INDEX="$SNAPSHOT_ROOT/index"
cp "$INDEX_PATH" "$SNAPSHOT_INDEX"
if ! GIT_INDEX_FILE="$SNAPSHOT_INDEX" git -c core.splitIndex=false update-index --no-split-index; then
  echo "[$LABEL] cannot build a standalone snapshot index."
  exit 1
fi

(
  for name in $LOCAL_GIT_ENV; do unset "$name"; done
  git -C "$SNAPSHOT" init --quiet
  mv "$SNAPSHOT_INDEX" "$SNAPSHOT/.git/index"

  require_clang_format() {
    CF="${CLANG_FORMAT:-clang-format}"
    if ! command -v "$CF" >/dev/null 2>&1; then
      echo "[$LABEL] $CF not found; install the pinned formatter."
      exit 1
    fi
    # Install clang-format 22 through the pinned requirements.
    "$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check-tool clang-format
  }

  require_ruff() {
    if ! command -v ruff >/dev/null 2>&1; then
      echo "[$LABEL] ruff not found; pip install -r requirements/ruff.txt"
      exit 1
    fi
    "$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check-tool ruff
  }

  require_eslint() {
    if [ ! -x node_modules/.bin/eslint ]; then
      echo "[$LABEL] node_modules is missing; run npm ci."
      exit 1
    fi
  }

  FORMAT="$(printf '%s\n' "$PRESENT" \
    | grep -E '\.(h|hpp|cpp|cc|inl|ino)$' \
    | grep -vE '(^|/)core/vendor/FastNoiseLite\.h$|(^|/)core/color/color_luts\.h$|(^|/)core/color/gamut_lut\.h$|(^|/)core/color/srgb_decode_lut\.h$|(^|/)effects/generated/|(^|/)core/mesh/relax_bakes_generated\.h$|(^|/)core/spatial/reaction_graph\.cpp$|(^|/)tests/mindsplatter_replay_corpus\.h$' \
    || true)"

  if [ -n "$FORMAT" ]; then
    require_clang_format
    CF="${CLANG_FORMAT:-clang-format}"
    failed=
    while IFS= read -r file; do
      extract_staged "$file"
      "$CF" --dry-run --Werror --style=file --assume-filename="$SNAPSHOT/$file" \
        <"$STAGED_BLOB" >/dev/null 2>&1 \
        || failed="$failed $file"
    done <<EOF
$FORMAT
EOF
    if [ -n "$failed" ]; then
      echo "[$LABEL] clang-format failed:$failed"
      exit 1
    fi
  fi

  PYTHON_FILES="$(printf '%s\n' "$PRESENT" | sed -n '/\.py$/p')"
  if [ -n "$PYTHON_FILES" ]; then
    require_ruff
    while IFS= read -r file; do
      extract_staged "$file"
      ruff check --no-cache --config "$SNAPSHOT/ruff.toml" --stdin-filename "$file" - <"$STAGED_BLOB"
    done <<EOF
$PYTHON_FILES
EOF
  fi

  JAVASCRIPT_FILES="$(printf '%s\n' "$PRESENT" | sed -n -E '/\.(js|cjs|mjs)$/p')"
  if [ -n "$JAVASCRIPT_FILES" ]; then
    require_eslint
    while IFS= read -r file; do
      extract_staged "$file"
      node_modules/.bin/eslint --config "$SNAPSHOT/eslint.config.mjs" --stdin --stdin-filename "$file" <"$STAGED_BLOB"
    done <<EOF
$JAVASCRIPT_FILES
EOF
  fi

  if printf '%s\n' "$CHANGED" | grep -qx '.clang-format'; then
    require_clang_format
    (cd "$SNAPSHOT" && bash tools/clang_format_gate.sh)
  fi
  if printf '%s\n' "$CHANGED" | grep -qx 'ruff.toml'; then
    require_ruff
    (cd "$SNAPSHOT" && ruff check --no-cache .)
  fi
  if printf '%s\n' "$CHANGED" | grep -qx 'eslint.config.mjs'; then
    require_eslint
    ESLINT="$PWD/node_modules/.bin/eslint"
    (cd "$SNAPSHOT" && "$ESLINT" .)
  fi

  SHELL_FILES="$(printf '%s\n' "$PRESENT" | grep -E '\.sh$|^\.githooks/' || true)"
  if [ -n "$SHELL_FILES" ]; then
    if ! command -v shellcheck >/dev/null 2>&1; then
      echo "[$LABEL] shellcheck not found; pip install -r requirements/shellcheck.txt"
      exit 1
    fi
    "$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check-tool shellcheck
    while IFS= read -r file; do
      (cd "$SNAPSHOT" && shellcheck -x "$file")
    done <<EOF
$SHELL_FILES
EOF
  fi

  # No --sync: a repair would land in the snapshot, not the index.
  "$PYTHON_BIN" "$SNAPSHOT/tools/docs_check.py" --root "$SNAPSHOT" --skip-checkout daydream
  "$PYTHON_BIN" "$SNAPSHOT/tools/docs_images.py" --root "$SNAPSHOT"
  "$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check
  "$PYTHON_BIN" "$SNAPSHOT/tools/license_check.py" --root "$SNAPSHOT"
  "$PYTHON_BIN" "$SNAPSHOT/tools/include_paths.py" --root "$SNAPSHOT"

  if printf '%s\n' "$CHANGED" | grep -Eq '^docs/profiles/|^scripts/(check_profiles|effect_roster|exit)\.mjs$|^targets/(effects\.h|Phantasm/phantasm_playlist\.h)$|^\.githooks/pre-commit$|^tools/commit_checks\.sh$'; then
    if ! command -v node >/dev/null 2>&1; then
      echo "[$LABEL] node not found; profile archive checks need the pinned Node version."
      exit 1
    fi
    "$PYTHON_BIN" "$SNAPSHOT/tools/build_pins.py" --check-tool node
    node "$SNAPSHOT/scripts/check_profiles.mjs"
  fi
)

echo "[$LABEL] $SUBJECT formatting, lint, and documentation checks passed."
