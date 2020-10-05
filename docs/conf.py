# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'AvaFrame'
copyright = '2020, AvaFrame developers'
author = 'AvaFrame developers'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autosectionlabel',
    'sphinxcontrib.bibtex'
]
autosectionlabel_prefix_document = True
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Set the master document name for readthedocs builds
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'logo_only': True,
    'style_nav_header_background': '#343131'
    # 'display_version': False,
}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/logo.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for LaTeX output ---------------------------------------------
latex_maketitle = r'''
\begin{titlepage}%
    \let\footnotesize\small
    \let\footnoterule\relax
    \rule{\textwidth}{1pt}%
    \ifsphinxpdfoutput
        \begingroup
        % These \defs are required to deal with multi-line authors; it
        % changes \\ to ', ' (comma-space), making it pass muster for
        % generating document info in the PDF file.
        \def\\{, }
        \def\and{and }
        \pdfinfo{
            /Author (\@author)
            /Title (\@title)
            }
        \endgroup
    \fi
    \begin{flushright}%
        \sphinxlogo%
        {\rm\Huge\py@HeaderFamily \@title \par}%
        {\em\LARGE\py@HeaderFamily \py@release\releaseinfo \par}
        \vfill
        {\LARGE\py@HeaderFamily
        \begin{tabular}[t]{c}
            \@author
        \end{tabular}
        \par}
        \vfill\vfill
        {\large
        \@date \par
        \vfill
        \py@authoraddress \par
        }%
    \end{flushright}%\par
    \@thanks
\end{titlepage}%
\cleardoublepage%
\setcounter{footnote}{0}%
\let\thanks\relax\let\maketitle\relax
%\gdef\@thanks{}\gdef\@author{}\gdef\@title{}
\pagenumbering{roman}
\tableofcontents
'''

latex_elements = {
    'sphinxsetup': '',
    'passoptionstopackages': r'\PassOptionsToPackage{table}{xcolor}',
    'maketitle': latex_maketitle,
}


# -- Options for referencing -------------------------------------------
numfig = True
math_numfig = True
math_eqref_format = "Eq.{number}"
# make referencing unique if twice the same section heading apear
autosectionlabel_prefix_document = True

def setup(app):
    app.add_stylesheet('css/custom.css')
