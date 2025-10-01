"""Legacy compatibility layer that re-exports the new ``ecomp`` package."""

from __future__ import annotations

import warnings

warnings.warn(
    "The `evolutionary-compression` package has been renamed to `ecomp`. "
    "Install and depend on `ecomp` directly.",
    DeprecationWarning,
    stacklevel=2,
)

from ecomp import *  # noqa: F401,F403
import ecomp as _ecomp

__all__ = getattr(_ecomp, "__all__", [])
__version__ = getattr(_ecomp, "__version__", "0.0.0")
