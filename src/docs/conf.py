
# -- Project information -----------------------------------------------------
project = 'iCEA'
author = 'Yohey Ishizu'
release = '0.2.2'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'ja'

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_show_copyright = False
html_show_sourcelink = False

# -- Options for todo extension ----------------------------------------------
todo_include_todos = True

# -- Additional CSS ----------------------------------------------------------
def setup(app):
    app.add_css_file('custom.css')
