# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Sphinx configuration for the Stereon manual.

The manual is built by the ``Stereon_Docs_Sphinx`` CMake target
(``-DSTEREON_BUILD_DOCS=ON``) or directly with
``sphinx-build -b html docs <out>``. API pages come from Doxygen XML through
Breathe. When the XML is missing (no Doxygen on the machine) the API page
shows a notice instead, so the rest of the manual still builds.
"""

import os
import re
import shutil
import subprocess

_DOCS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_DOCS_DIR, '..'))


def _read_project_version():
    """Return the VERSION from the top-level CMakeLists.txt project() call."""
    with open(os.path.join(_ROOT_DIR, 'CMakeLists.txt'), encoding='utf-8') as f:
        match = re.search(r'project\s*\(\s*Stereon\s+VERSION\s+([0-9][0-9.]*)', f.read())
    return match.group(1) if match else '0.0.0'


# -- Project information -----------------------------------------------------

project = 'Stereon'
author = 'Onur Tuncer and Stereon contributors'
copyright = '2026, Onur Tuncer and Stereon contributors'
version = _read_project_version()
release = version

# -- General configuration ---------------------------------------------------

extensions = [
    'breathe',
    'sphinx.ext.graphviz',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
    'sphinxcontrib.bibtex',
]

source_encoding = 'utf-8'
templates_path = []
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Placeholder pages are marked with ``.. todo::`` and collected on the index
# page by ``.. todolist::``. Switch this off for a release build.
todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_logo = '../assets/stereon_logo.png'
html_static_path = []  # add 'docs/_static' here once there are files to ship
html_theme_options = {
    'logo_only': False,
    'navigation_depth': 3,
}

# -- Doxygen / Breathe -------------------------------------------------------

# CMake passes the XML location through STEREON_DOXYGEN_XML. The fallback is
# where a manual ``doxygen docs/Doxyfile`` run (or the Read the Docs hook
# below) writes it.
_doxygen_xml = os.environ.get(
    'STEREON_DOXYGEN_XML',
    os.path.join(_DOCS_DIR, '_build', 'doxygen', 'xml'))


def _run_doxygen():
    """Configure Doxyfile.in by hand and run Doxygen (no CMake available)."""
    doxygen = shutil.which('doxygen')
    if doxygen is None:
        print('conf.py: doxygen not found; API pages will be skipped')
        return
    output_dir = os.path.join(_DOCS_DIR, '_build', 'doxygen')
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(_DOCS_DIR, 'Doxyfile.in'), encoding='utf-8') as f:
        doxy_config = f.read()
    for key, value in {
        'CMAKE_SOURCE_DIR': _ROOT_DIR,
        'CMAKE_PROJECT_NAME': project,
        'PROJECT_VERSION': version,
        'PROJECT_DESCRIPTION': 'A modern C++26 B-rep geometric modeling kernel',
        'DOXYGEN_OUTPUT_DIR': output_dir,
    }.items():
        doxy_config = doxy_config.replace(f'@{key}@', value)
    doxyfile = os.path.join(_DOCS_DIR, '_build', 'Doxyfile')
    with open(doxyfile, 'w', encoding='utf-8') as f:
        f.write(doxy_config)
    subprocess.run([doxygen, doxyfile], check=True, cwd=_DOCS_DIR)


# Read the Docs (and anyone setting STEREON_RUN_DOXYGEN=1) has no CMake step,
# so Doxygen runs from here.
if os.environ.get('READTHEDOCS') == 'True' or os.environ.get('STEREON_RUN_DOXYGEN') == '1':
    _run_doxygen()

breathe_projects = {project: _doxygen_xml}
breathe_default_project = project
breathe_default_members = ('members',)
add_function_parentheses = True

# ``api.rst`` wraps its Breathe directives in ``.. only:: have_doxygen``.
# Breathe still parses directives inside a false ``only`` block and fails on
# missing XML, so without XML the block between the ``stereon-api-begin`` and
# ``stereon-api-end`` comments is also cut from the page source.
_HAVE_DOXYGEN = os.path.isfile(os.path.join(_doxygen_xml, 'index.xml'))
if _HAVE_DOXYGEN:
    tags.add('have_doxygen')  # noqa: F821  (``tags`` is injected by Sphinx)
    print('conf.py: Doxygen XML found at', _doxygen_xml)
else:
    print('conf.py: no Doxygen XML at', _doxygen_xml, '- API pages skipped')

_API_BLOCK = re.compile(r'^\.\. stereon-api-begin$.*?^\.\. stereon-api-end$', re.M | re.S)


def _strip_api_block(app, docname, source):
    if not _HAVE_DOXYGEN:
        source[0] = _API_BLOCK.sub('', source[0])


def setup(app):
    app.connect('source-read', _strip_api_block)

# -- MathJax macros ----------------------------------------------------------

mathjax3_config = {
    'tex': {
        'macros': {
            'R': r'\mathbb{R}',
            'norm': [r'\left\lVert #1 \right\rVert', 1],
            'abs': [r'\left\lvert #1 \right\rvert', 1],
            'tol': r'\varepsilon',
            'coloneqq': r'\mathrel{:=}',
        }
    }
}

# -- BibTeX references -------------------------------------------------------

bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'unsrt'
