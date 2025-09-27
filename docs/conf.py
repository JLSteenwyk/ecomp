"""Sphinx configuration for the evolutionary compression project."""

from __future__ import annotations

import datetime
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists() and str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

project = "Evolutionary Compression"
current_year = datetime.date.today().year
copyright = f"{current_year} Evolutionary Compression Team"
author = "Evolutionary Compression Team"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = ".rst"
master_doc = "index"
language = "en"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
}
html_logo = None
html_show_sourcelink = True

pygments_style = "sphinx"

napoleon_google_docstring = True
napoleon_numpy_docstring = True

rst_epilog = "\n.. |project| replace:: Evolutionary Compression\n"
