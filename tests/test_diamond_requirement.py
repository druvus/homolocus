"""DIAMOND tests skip locally but fail in CI when DIAMOND is missing (spec 15)."""

from pathlib import Path

import pytest

pytest_plugins = ["pytester"]

CONFTEST = Path(__file__).with_name("conftest.py").read_text(encoding="utf-8")


def run_without_diamond(pytester: pytest.Pytester, monkeypatch, require: str | None):
    pytester.makeconftest(CONFTEST)
    pytester.makeini("[pytest]\nmarkers =\n    diamond: needs DIAMOND\n")
    pytester.makepyfile(
        """
        import pytest

        @pytest.mark.diamond
        def test_needs_diamond():
            pass

        def test_plain():
            pass
        """
    )
    monkeypatch.setenv("PATH", "")
    if require is None:
        monkeypatch.delenv("HOMOLOCUS_REQUIRE_DIAMOND", raising=False)
    else:
        monkeypatch.setenv("HOMOLOCUS_REQUIRE_DIAMOND", require)
    return pytester.runpytest("-rs")


def test_missing_diamond_skips_tests_by_default(pytester, monkeypatch) -> None:
    result = run_without_diamond(pytester, monkeypatch, require=None)
    result.assert_outcomes(passed=1, skipped=1)
    result.stdout.fnmatch_lines(["*DIAMOND executable not found on PATH*"])


def test_missing_diamond_fails_when_required(pytester, monkeypatch) -> None:
    result = run_without_diamond(pytester, monkeypatch, require="1")
    assert result.ret != 0
    result.stderr.fnmatch_lines(["*HOMOLOCUS_REQUIRE_DIAMOND=1 but DIAMOND executable not found*"])


@pytest.mark.parametrize("value", ["", "0"])
def test_only_one_enables_the_requirement(pytester, monkeypatch, value: str) -> None:
    run_without_diamond(pytester, monkeypatch, require=value).assert_outcomes(passed=1, skipped=1)
