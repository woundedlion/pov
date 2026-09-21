#!/usr/bin/env bash
# Assert every tracked file whose .gitattributes entry declares an `eol` value
# actually carries those line endings, in the index blob and in the working copy.
#
# A working copy that diverged from its eol=lf blob -- checked out before the
# attribute was added, or rewritten by a CRLF-emitting tool -- is invisible to
# `git status`, so the shell, clang-format and mirror gates then read different
# bytes locally than they do on a fresh checkout.
#
# usage: eol_gate.sh
set -eu

if [ "$#" -ne 0 ]; then
  echo "usage: $0" >&2
  exit 2
fi

tmp=$(mktemp)
trap 'rm -f -- "$tmp"' EXIT

git ls-files --eol -z > "$tmp"

declared=0
diverged=0

# Record layout: "i/<kind>  w/<kind>  attr/<attributes>" then a TAB then the
# path. `none` is a file with no line ending at all and `-text` is binary
# content; neither contradicts a declared eol.
while IFS= read -r -d '' record; do
  info=${record%%$'\t'*}
  path=${record#*$'\t'}

  case $info in
    *'eol=lf'*) want=lf ;;
    *'eol=crlf'*) want=crlf ;;
    *) continue ;;
  esac
  declared=$((declared + 1))

  index=${info#i/}
  index=${index%% *}
  work=${info#*w/}
  work=${work%% *}

  for probe in "index:$index" "worktree:$work"; do
    where=${probe%%:*}
    have=${probe#*:}
    case $have in
      "$want" | none | -text) continue ;;
    esac
    echo "::error file=$path::$where line endings are $have, .gitattributes declares eol=$want"
    diverged=$((diverged + 1))
  done
done < "$tmp"

if [ "$declared" -eq 0 ]; then
  echo "no tracked file declares an eol attribute -- the selection is broken" >&2
  exit 1
fi

if [ "$diverged" -ne 0 ]; then
  cat >&2 <<'REMEDY'

To restore a diverged working copy:   rm -f -- <path> && git checkout -- <path>
To restore a diverged index blob:     git add --renormalize -- <path>
REMEDY
  exit 1
fi

echo "eol clean: $declared files declare an eol attribute."
