"""Compatibility helpers for differing Python releases."""

from __future__ import annotations

import sys
from dataclasses import dataclass as _dataclass


def dataclass(*args, **kwargs):
    """dataclass decorator that uses ``slots`` when supported."""

    if sys.version_info >= (3, 10):
        kwargs.setdefault("slots", True)
    else:
        kwargs.pop("slots", None)
    return _dataclass(*args, **kwargs)


__all__ = ["dataclass"]
