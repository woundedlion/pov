import json
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
UNPLACED_PROJECT = GEN.parent / "1.2" / "phantasm.kicad_pro"
sys.path.insert(0, str(GEN))

from constraints import (DEFAULT_CLASS_MINIMUMS, NEW_LAYOUT_RULES, RULE_MINIMUMS,  # noqa: E402
                         UNPLACED_DEFAULT_CLASS, UNPLACED_RULES, rule_checks,
                         rule_shortfalls)


def default_class(project):
    return next(item for item in project["net_settings"]["classes"]
                if item["name"] == "Default")


class UnplacedProjectConstraintTests(unittest.TestCase):
    """Captured unplaced values that the generator's constraint floors must not undercut."""

    def setUp(self):
        self.project = json.loads(UNPLACED_PROJECT.read_text(encoding="utf-8"))

    def test_rules_match_the_captured_values(self):
        rules = self.project["board"]["design_settings"]["rules"]
        for field, expected in UNPLACED_RULES.items():
            with self.subTest(field=field):
                self.assertEqual(rules[field], expected)

    def test_unplaced_project_satisfies_new_layout_floors(self):
        settings = self.project["board"]["design_settings"]
        for field, minimum in NEW_LAYOUT_RULES.items():
            with self.subTest(field=field):
                self.assertGreaterEqual(settings["rules"][field], minimum)
        self.assertEqual(settings["rule_severities"]["silk_over_copper"], "error")

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


class RuleCheckTests(unittest.TestCase):
    def project(self):
        return {
            "board": {"design_settings": {
                "rules": {**RULE_MINIMUMS, **NEW_LAYOUT_RULES},
                "rule_severities": {"silk_over_copper": "error"}}},
            "net_settings": {"classes": [dict(DEFAULT_CLASS_MINIMUMS, name="Default")]},
        }

    def test_checks_every_floor_and_the_silk_severity(self):
        checks = rule_checks(self.project(), RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS)
        self.assertEqual(
            set(checks),
            {*RULE_MINIMUMS, *NEW_LAYOUT_RULES, "rule_severities.silk_over_copper",
             *(f"Default.{field}" for field in DEFAULT_CLASS_MINIMUMS)})
        self.assertTrue(all(met for _, _, met in checks.values()))
        self.assertEqual(rule_shortfalls(self.project(), RULE_MINIMUMS,
                                         DEFAULT_CLASS_MINIMUMS), {})

    def test_shortfalls_are_the_unmet_checks(self):
        project = self.project()
        project["board"]["design_settings"]["rule_severities"]["silk_over_copper"] = "warning"
        project["net_settings"]["classes"][0]["via_drill"] = 0
        self.assertEqual(
            rule_shortfalls(project, RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS),
            {"rule_severities.silk_over_copper": ("warning", "error"),
             "Default.via_drill": (0, DEFAULT_CLASS_MINIMUMS["via_drill"])})


if __name__ == "__main__":
    unittest.main()
