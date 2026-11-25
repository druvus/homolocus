"""Stand-in executables for tests that must not need the real tools."""

import stat
from pathlib import Path


def fake_diamond(tmp_path: Path, script: str) -> str:
    """An executable shell script standing in for DIAMOND; returns its path."""
    path = tmp_path / "fake-diamond"
    path.write_text("#!/bin/sh\n" + script, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)
