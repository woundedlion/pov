import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import docs_check
import docs_sync
import engine_source_state


class EngineSourceState(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.git("init", "-q")
        self.git("config", "user.name", "Source state test")
        self.git("config", "user.email", "source-state@example.invalid")
        self.git("config", "core.autocrlf", "false")
        self.git("config", "core.hooksPath", ".no-hooks")
        self.original = "Authored explanation.\n\n<!-- docs-check: tree exhaustive -->\n```\n├── absent.txt  Stale entry\n├── old.txt  Existing description\n└── new.txt  New source\n```\n"
        self.write("README.md", self.original)
        self.write("old.txt", "source")
        self.write("new.txt", "added source")
        self.git("add", "README.md", "old.txt", "new.txt")
        self.git("commit", "-qm", "fixture")

    def git(self, *arguments):
        return subprocess.check_output(
            ["git", "-C", str(self.root), *arguments], stderr=subprocess.STDOUT)

    def write(self, name, text):
        (self.root / name).write_text(text, encoding="utf-8", newline="\n")

    def generate(self):
        _, entries = docs_check._tracked_entries(self.root)
        with contextlib.redirect_stdout(io.StringIO()):
            generated = docs_sync.sync_text(
                PurePosixPath("README.md"), self.original, entries, {}, {})
        self.assertIn("new.txt", generated)
        self.write("README.md", generated)
        return generated

    def test_clean_checkout_and_generated_documentation_are_clean(self):
        self.assertEqual(engine_source_state.changed_sources(self.root), [])
        self.generate()
        self.assertEqual(engine_source_state.changed_sources(self.root), [])
        self.git("add", "README.md")
        self.assertEqual(engine_source_state.changed_sources(self.root), [])

    def test_authored_prose_change_remains_dirty(self):
        generated = self.generate()
        self.write("README.md", generated.replace("Authored explanation", "Edited explanation"))
        self.assertEqual(engine_source_state.changed_sources(self.root), ["README.md"])

    def test_source_edits_remain_dirty_with_generated_docs(self):
        self.generate()
        self.write("old.txt", "changed source")
        self.assertIn("old.txt", engine_source_state.changed_sources(self.root))

    def test_deleted_and_new_documents_remain_dirty(self):
        (self.root / "README.md").unlink()
        self.assertIn("README.md", engine_source_state.changed_sources(self.root))
        self.write("README.md", self.original)
        self.write("new.md", "New prose\n")
        self.git("add", "new.md")
        self.assertIn("new.md", engine_source_state.changed_sources(self.root))
