import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import docs_images as di  # noqa: E402


class TestDocsImages(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name) / "repo"
        self.html = self.repo / "build" / "docs" / "html"
        self.html.mkdir(parents=True)
        self.addCleanup(self.tmp.cleanup)

    def write_page(self, body, name="index.html"):
        (self.html / name).write_text(body, encoding="utf-8")

    def write_asset(self, relative, data=b"png"):
        path = self.repo / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def test_repo_relative_image_is_staged_into_the_artifact(self):
        self.write_page('<img src="docs/screenshots/IslamicStars.png" width="640">')
        self.write_asset("docs/screenshots/IslamicStars.png", b"pixels")
        errors, staged = di.resolve(self.html, self.repo)
        self.assertEqual(errors, [])
        self.assertEqual(staged, 1)
        copied = self.html / "docs" / "screenshots" / "IslamicStars.png"
        self.assertEqual(copied.read_bytes(), b"pixels")

    def test_image_already_in_the_artifact_is_not_restaged(self):
        self.write_page('<img src="inherit_graph_0.png" alt="">')
        (self.html / "inherit_graph_0.png").write_bytes(b"generated")
        errors, staged = di.resolve(self.html, self.repo)
        self.assertEqual((errors, staged), ([], 0))

    def test_reference_missing_everywhere_is_reported(self):
        self.write_page('<img alt="gone" src="docs/screenshots/Gone.png">')
        errors, staged = di.resolve(self.html, self.repo)
        self.assertEqual(staged, 0)
        self.assertEqual(len(errors), 1)
        self.assertIn("docs/screenshots/Gone.png", errors[0])

    def test_absolute_and_data_sources_are_left_alone(self):
        self.write_page('<img src="https://example.com/a.png">'
                        '<img src="//example.com/b.png">'
                        '<IMG SRC="data:image/png;base64,AAAA">')
        self.assertEqual(di.resolve(self.html, self.repo), ([], 0))

    def test_reference_escaping_the_artifact_is_reported(self):
        self.write_page('<img src="../../../secrets/key.png">')
        self.write_asset("secrets/key.png")
        errors, _ = di.resolve(self.html, self.repo)
        self.assertEqual(len(errors), 1)
        self.assertIn("outside the artifact", errors[0])

    def test_query_and_fragment_are_stripped_before_resolving(self):
        self.write_page('<img src="docs/screenshots/A.png?v=2#top">')
        self.write_asset("docs/screenshots/A.png")
        self.assertEqual(di.resolve(self.html, self.repo), ([], 1))

    def test_percent_escapes_are_decoded(self):
        self.write_page('<img src="docs/screenshots/Ring%20Spin.png">')
        self.write_asset("docs/screenshots/Ring Spin.png")
        self.assertEqual(di.resolve(self.html, self.repo), ([], 1))

    def test_nested_pages_resolve_against_their_own_directory(self):
        (self.html / "sub").mkdir()
        self.write_page('<img src="../docs/screenshots/B.png">', name="sub/page.html")
        self.write_asset("docs/screenshots/B.png")
        self.assertEqual(di.resolve(self.html, self.repo), ([], 1))

    def test_empty_source_is_reported(self):
        self.write_page('<img src="">')
        errors, _ = di.resolve(self.html, self.repo)
        self.assertEqual(len(errors), 1)
        self.assertIn("no path component", errors[0])


if __name__ == "__main__":
    unittest.main()
