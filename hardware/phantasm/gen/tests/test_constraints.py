import json
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
PROJECT = GEN.parent / "phantasm.kicad_pro"
UNPLACED_PROJECT = GEN.parent / "unplaced" / "phantasm_unplaced.kicad_pro"
sys.path.insert(0, str(GEN))

from constraints import (DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS,  # noqa: E402
                         UNPLACED_DEFAULT_CLASS, UNPLACED_RULES)


def default_class(project):
    return next(item for item in project["net_settings"]["classes"]
                if item["name"] == "Default")


class CommittedProjectConstraintTests(unittest.TestCase):
    def test_project_satisfies_fabrication_constraints(self):
        project = json.loads(PROJECT.read_text(encoding="utf-8"))
        rules = project["board"]["design_settings"]["rules"]
        for field, minimum in RULE_MINIMUMS.items():
            with self.subTest(field=field):
                self.assertGreaterEqual(rules[field], minimum)

        default = default_class(project)
        for field, minimum in DEFAULT_CLASS_MINIMUMS.items():
            with self.subTest(field=f"Default.{field}"):
                self.assertGreaterEqual(default[field], minimum)


class UnplacedProjectConstraintTests(unittest.TestCase):
    """The unplaced Quilter project has no generator, so nothing but this pin
    would notice its constraints drifting to the routed board's floors."""

    def setUp(self):
        self.project = json.loads(UNPLACED_PROJECT.read_text(encoding="utf-8"))

    def test_rules_match_the_captured_values(self):
        rules = self.project["board"]["design_settings"]["rules"]
        for field, expected in UNPLACED_RULES.items():
            with self.subTest(field=field):
                self.assertEqual(rules[field], expected)

    def test_default_net_class_matches_the_captured_values(self):
        default = default_class(self.project)
        for field, expected in UNPLACED_DEFAULT_CLASS.items():
            with self.subTest(field=f"Default.{field}"):
                self.assertEqual(default[field], expected)

    def test_constraints_stay_wider_than_the_routed_board(self):
        rules = self.project["board"]["design_settings"]["rules"]
        for field, minimum in RULE_MINIMUMS.items():
            with self.subTest(field=field):
                self.assertGreaterEqual(rules[field], minimum)
        default = default_class(self.project)
        for field, minimum in DEFAULT_CLASS_MINIMUMS.items():
            with self.subTest(field=f"Default.{field}"):
                self.assertGreaterEqual(default[field], minimum)


if __name__ == "__main__":
    unittest.main()
