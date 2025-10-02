"""Tests for compatibility helpers in :mod:`ecomp._compat`."""

from __future__ import annotations

import sys

import pytest

from ecomp import _compat


@pytest.mark.skipif(sys.version_info < (3, 10), reason="slots require Python >= 3.10")
def test_dataclass_sets_slots_on_modern_python(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("ecomp._compat.sys.version_info", (3, 11, 0))

    @_compat.dataclass
    class WithSlots:
        value: int

    instance = WithSlots(7)
    assert hasattr(WithSlots, "__slots__")
    with pytest.raises(AttributeError):
        instance.extra = 1


def test_dataclass_strips_slots_on_legacy_python(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("ecomp._compat.sys.version_info", (3, 9, 0))

    @_compat.dataclass(slots=True)
    class WithoutSlots:
        value: int

    instance = WithoutSlots(3)
    assert not hasattr(WithoutSlots, "__slots__")
    instance.extra = 2
    assert instance.extra == 2


@pytest.mark.skipif(sys.version_info < (3, 10), reason="zip(strict=...) requires Python >= 3.10")
def test_zip_strict_modern_python(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("ecomp._compat.sys.version_info", (3, 11, 0))

    result = list(_compat.zip_strict([1, 2], ["a", "b"]))
    assert result == [(1, "a"), (2, "b")]

    with pytest.raises(ValueError):
        list(_compat.zip_strict([1, 2], ["only"]))


def test_zip_strict_legacy_python(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("ecomp._compat.sys.version_info", (3, 8, 0))

    result = list(_compat.zip_strict([1, 2], ["a", "b"]))
    assert result == [(1, "a"), (2, "b")]

    with pytest.raises(ValueError, match="different lengths"):
        list(_compat.zip_strict([1], ["alpha", "beta"]))
