# -*- coding: utf-8 -*-
"""Sphinx configuration for the eComp documentation site."""

import os
import sys

import sphinx_rtd_theme

# -- Path setup --------------------------------------------------------------
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(BASE_DIR, "..", "src")))

# -- Project information -----------------------------------------------------
project = "eComp"
author = "Evolutionary Compression Developers"
copyright = "2024"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_logo = "_static/img/logo.png"
html_favicon = "_static/img/favicon.ico"
html_static_path = ["_static"]
html_show_sourcelink = False
html_theme_options = {
    "body_max_width": "900px",
    "logo_only": True,
}

html_sidebars = {
    "**": [
        "sidebar-top.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
    ]
}

# -- General content options -------------------------------------------------
pygments_style = "friendly"
language = "en"
