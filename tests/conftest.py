"""Pytest configuration mirroring the ClipKIT project conventions."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest


_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_PATH = _PROJECT_ROOT / "src"
if str(_SRC_PATH) not in sys.path:
    sys.path.insert(0, str(_SRC_PATH))


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "integration: mark integration tests")
    config.addinivalue_line("markers", "slow: mark slow tests")
