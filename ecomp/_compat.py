"""Compatibility helpers for differing Python releases."""

from __future__ import annotations

import sys
from dataclasses import dataclass as _dataclass
from itertools import zip_longest
from typing import Iterable, Iterator, Tuple, TypeVar


def dataclass(*args, **kwargs):
    """dataclass decorator that uses ``slots`` when supported."""

    if sys.version_info >= (3, 10):
        kwargs.setdefault("slots", True)
    else:
        kwargs.pop("slots", None)
    return _dataclass(*args, **kwargs)


_Sentinel = object()
_T = TypeVar("_T")


def zip_strict(*iterables: Iterable[_T]) -> Iterator[Tuple[_T, ...]]:
    """Backport of ``zip(strict=True)`` for Python < 3.10."""

    if sys.version_info >= (3, 10):
        yield from zip(*iterables, strict=True)
        return

    for values in zip_longest(*iterables, fillvalue=_Sentinel):
        if _Sentinel in values:
            if any(value is not _Sentinel for value in values):  # pragma: no branch - rare error path
                raise ValueError("zip() argument iterables have different lengths")
            break
        yield values


__all__ = ["dataclass", "zip_strict"]
