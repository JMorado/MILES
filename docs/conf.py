"""Sphinx configuration for the MILES documentation.

Build with ``sphinx-build -b html docs docs/_build/html``.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from miles import __version__

project = "MILES"
author = "The MILES authors"
copyright = "2026, the MILES authors"
version = release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

numfig = True
exclude_patterns = ["_build"]

autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_default_options = {"members": True, "show-inheritance": True}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

html_theme = "furo"
html_title = f"MILES {version}"
html_logo = "../miles/ui/images/miles.svg"   # the app's own logo, so there is one file to edit
html_favicon = "_static/favicon.png"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False
html_theme_options = {
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#3A7BFF",
        "color-brand-content": "#3A7BFF",
        "color-background-primary": "#FFFFFF",
        "color-background-secondary": "#F6F7F9",
        "font-stack": "-apple-system, BlinkMacSystemFont, 'Inter', 'Segoe UI', Helvetica, Arial, sans-serif",
        "font-stack--monospace": "'SF Mono', 'JetBrains Mono', Menlo, Consolas, monospace",
    },
    "dark_css_variables": {
        "color-brand-primary": "#5B93FF",
        "color-brand-content": "#5B93FF",
        "color-background-primary": "#141619",
        "color-background-secondary": "#1B1E22",
    },
}
