"""Inject the profile image's Git SHA into pullback telemetry."""

# ruff: noqa: F821

import subprocess


Import("env")  # noqa: F821

if env["PIOENV"] in ("profile", "profile_o3"):
    project_dir = env["PROJECT_DIR"]
    try:
        short_sha = subprocess.run(
            ["git", "-C", project_dir, "rev-parse", "--short=12", "HEAD"],
            check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"pullback_profile_build: cannot resolve Git SHA: {error}")
    env.Append(CPPDEFINES=[("HS_PULLBACK_SHORT_SHA", f'\\"{short_sha}\\"')])
