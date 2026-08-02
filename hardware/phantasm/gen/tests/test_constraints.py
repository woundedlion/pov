import json
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
PROJECT = GEN.parent / "phantasm.kicad_pro"
sys.path.insert(0, str(GEN))

from constraints import DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS  # noqa: E402


class CommittedProjectConstraintTests(unittest.TestCase):
    def test_project_satisfies_fabrication_constraints(self):
        project = json.loads(PROJECT.read_text(encoding="utf-8"))
        rules = project["board"]["design_settings"]["rules"]
        for field, minimum in RULE_MINIMUMS.items():
            with self.subTest(field=field):
                self.assertGreaterEqual(rules[field], minimum)

        classes = project["net_settings"]["classes"]
        default = next(item for item in classes if item["name"] == "Default")
        for field, minimum in DEFAULT_CLASS_MINIMUMS.items():
            with self.subTest(field=f"Default.{field}"):
                self.assertGreaterEqual(default[field], minimum)


if __name__ == "__main__":
    unittest.main()
