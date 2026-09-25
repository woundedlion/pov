"""Exercise the optional engine-completion notification without network access."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[2]


class EngineDispatchTests(unittest.TestCase):
    def run_notification(self, token="", sha="a" * 40, run_id="1234", status=0):
        source = (ROOT / ".github/workflows/notify-daydream.yml").read_text(encoding="utf-8")
        body = textwrap.dedent(source.split("        run: |\n", 1)[1])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "bin").mkdir()
            (root / "notify.sh").write_text(body, encoding="utf-8", newline="\n")
            (root / "bin/gh").write_text(
                '#!/usr/bin/env bash\nprintf "%s\\n" "$@" > args.txt\n'
                f"exit {status}\n", encoding="utf-8", newline="\n")
            (root / "bin/gh").chmod(0o755)
            env = dict(os.environ, GH_TOKEN=token, ENGINE_SHA=sha,
                       ENGINE_RUN_ID=run_id, RUNNER_TEMP=str(root))
            result = subprocess.run(
                ["bash", "-c", 'export PATH="$PWD/bin:$PATH"; bash -e notify.sh'],
                cwd=root, env=env, capture_output=True, text=True, timeout=15)
            args = (root / "args.txt").read_text().splitlines() if (root / "args.txt").exists() else []
            payload = root / "daydream-dispatch.json"
            return result, args, json.loads(payload.read_text()) if payload.exists() else None

    def test_missing_token_does_not_call_the_api(self):
        result, args, payload = self.run_notification()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args, [])
        self.assertIsNone(payload)
        self.assertIn("scheduled reconciliation", result.stdout)

    def test_success_dispatches_the_exact_completed_run(self):
        result, args, payload = self.run_notification(token="fixture-token")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args[:4], ["api", "--method", "POST", "repos/woundedlion/daydream/dispatches"])
        self.assertEqual(payload, {"event_type": "holosphere-engine-ready", "client_payload": {
            "engine_sha": "a" * 40, "engine_run_id": 1234}})

    def test_api_failure_leaves_scheduled_reconciliation_available(self):
        result, args, _ = self.run_notification(token="fixture-token", status=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(args)
        self.assertIn("::warning::", result.stdout)

    def test_invalid_run_identity_never_reaches_the_api(self):
        for sha, run_id in [("a" * 40 + "-dirty", "1234"), ("a" * 40, "not-a-run")]:
            with self.subTest(sha=sha, run_id=run_id):
                result, args, payload = self.run_notification(token="fixture-token", sha=sha, run_id=run_id)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(args, [])
                self.assertIsNone(payload)


if __name__ == "__main__":
    unittest.main()
