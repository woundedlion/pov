# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Placement and inlining annotation contracts."""

from pathlib import Path
import re
import subprocess
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[1]
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))

from doc_coverage import CPP_LEXEMES  # noqa: E402


def inline_conflicts(source: str) -> list[int]:
    """Return lines combining cold placement and forced hot inlining."""
    source = re.sub(r"\\\r?\n", "", source)
    masked = CPP_LEXEMES.sub(lambda m: re.sub(r"[^\n]", " ", m.group()), source)
    directives = re.compile(r"(?m)^[ \t]*#[^\n]*")
    declarations = [(m.start(), m.group()) for m in directives.finditer(masked)
                    if re.match(r"[ \t]*#[ \t]*define\b", m.group())]
    code = directives.sub(lambda m: " " * len(m.group()), masked)
    declarations.extend((m.start(), m.group()) for m in re.finditer(r"[^;{}]+", code))
    conflicts = []
    for offset, declaration in declarations:
        tokens = list(re.finditer(r"\bHS_(?:FLASH|HOT)_INLINE\b", declaration))
        if {token.group() for token in tokens} == {"HS_FLASH_INLINE", "HS_HOT_INLINE"}:
            conflicts.append(source.count("\n", 0, offset + tokens[0].start()) + 1)
    return sorted(set(conflicts))


class InlineAttributeTests(unittest.TestCase):
    def test_rejects_stacked_attributes(self):
        for source in (
            "HS_FLASH_INLINE HS_HOT_INLINE static void f();",
            "HS_HOT_INLINE\nHS_FLASH_INLINE static void f() {}",
            "HS_FLASH_INLINE /* placement */ HS_HOT_INLINE void f();",
            "#define COMBINED HS_FLASH_INLINE \\\nHS_HOT_INLINE\n",
        ):
            with self.subTest(source=source):
                self.assertTrue(inline_conflicts(source))

    def test_accepts_distinct_declarations(self):
        self.assertEqual(inline_conflicts(
            "HS_FLASH_INLINE void setup(); HS_O3_FN HS_HOT_INLINE void shade() {}"), [])

    def test_ignores_comments_strings_and_preprocessor_conditions(self):
        self.assertEqual(inline_conflicts(
            '// HS_FLASH_INLINE HS_HOT_INLINE\n'
            '/* HS_FLASH_INLINE HS_HOT_INLINE */\n'
            'auto s = "HS_FLASH_INLINE HS_HOT_INLINE";\n'
            'auto raw = R"tag(HS_FLASH_INLINE HS_HOT_INLINE)tag";\n'
            '#if defined(HS_FLASH_INLINE) && defined(HS_HOT_INLINE)\n'
            '#endif\n'), [])

    def test_macro_definitions_are_separate(self):
        self.assertEqual(inline_conflicts(
            '#define HS_FLASH_INLINE __attribute__((cold))\n'
            '#define HS_HOT_INLINE __attribute__((always_inline))\n'), [])

    def test_first_party_declarations_keep_annotations_exclusive(self):
        paths = subprocess.check_output(
            ['git', '-C', str(ROOT), 'ls-files', '-z', '--', '*.h', '*.hpp',
             '*.cpp', '*.cc', '*.inl', '*.ino'], text=True, encoding='utf-8').split('\0')
        self.assertTrue(any(paths))
        for path in filter(None, paths):
            if path.startswith('core/vendor/'):
                continue
            with self.subTest(path=path):
                self.assertEqual(inline_conflicts((ROOT / path).read_text(encoding='utf-8')), [])


if __name__ == '__main__':
    unittest.main()
