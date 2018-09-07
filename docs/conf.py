# -*- coding: utf-8 -*-

import os

XML_PATH = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), 'doxygen', 'xml')

import parpydtk2


def generate_doxygen():
    import subprocess
    if subprocess.call('cd doxygen; doxygen Doxyfile', shell=True):
        raise OSError('doxygen failed')


generate_doxygen()


project = 'ParPyDTK2'
copyright = '2018, Qiao Chen'
author = 'Qiao Chen'
version = parpydtk2.__version__
release = version

breathe_projects = {
    'parpydtk2': XML_PATH
}
breathe_default_project = 'parpydtk2'
breathe_default_members = ('members', 'protected-members')
autoclass_content = 'both'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.imgmath',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'breathe',
    'sphinx.ext.napoleon',
    'numpydoc'
]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = None
pygments_style = 'sphinx'
html_theme = 'classic'
# html_theme_options = {}
html_static_path = ['_static']
html_logo = '_static/logo.png'
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'ParPyDTK2doc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'extraclassoptions': ',openany,oneside',
}

latex_documents = [
    (master_doc, 'ParPyDTK2.tex', 'ParPyDTK2 Documentation',
     'Qiao Chen', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'parpydtk2', 'ParPyDTK2 Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'ParPyDTK2', 'ParPyDTK2 Documentation',
     author, 'ParPyDTK2', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
todo_include_todos = True

numpydoc_show_class_members = False
