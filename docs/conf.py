# Sphinx configuration for HOTPANTS documentation
import os
import sys

# Add parent directory to path for potential autodoc imports
sys.path.insert(0, os.path.abspath('..'))

# --- Project information ---
project = 'HOTPANTS'
copyright = '2013–2026, Andy Becker'
author = 'Andy Becker (original); Maintained by the HOTPANTS community'
release = '0.0.1'

# --- General configuration ---
extensions = [
    'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# --- Options for Breathe (Doxygen integration) ---
breathe_projects = {
    'HOTPANTS': os.path.join(os.path.dirname(__file__), '_doxygen', 'xml'),
}
breathe_default_project = 'HOTPANTS'
breathe_show_include = False

# --- Options for HTML output ---
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'logo_only': False,
}

# --- MathJax for equation rendering ---
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'

# --- Syntax highlighting ---
pygments_style = 'sphinx'

# --- Suppress warnings for expected Breathe behavior ---
suppress_warnings = ['ref.c:member']
