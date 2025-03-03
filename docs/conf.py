# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
import tomli
import shutil
sys.path.insert(
    0, os.path.abspath('..'))
# sys.path.append(".")

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with open(os.path.join(os.path.abspath('..'), "pyproject.toml"), "rb") as f:
    toml = tomli.load(f)

project = 'flight-mech'
copyright = '2025, PaulCreusy'
author = 'PaulCreusy'

# version = '1.0.6'
# release = '1.0.6'

pyproject = toml["project"]

version = pyproject["version"]
release = pyproject["version"]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx_rtd_theme',
    'sphinx.ext.napoleon',
    'nbsphinx',
    'sphinx_sitemap'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme'
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/PaulCreusy/flight-mech",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ],
    "logo": {
        "text": "Flight-Mech",
        "image_light": "figures/logo-light.png",
        "image_dark": "figures/logo-dark.png",
    }
}

# html_logo = "figures/logo.png"
# html_title = "Flight-Mech"

html_favicon = "figures/logo-light.png"

# -- Options for todo extension ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/todo.html#configuration

todo_include_todos = True

# -- Options for sitemap -----------------------------------------------------

html_baseurl = 'https://flight-mech.creusy.fr/'
sitemap_url_scheme = "{link}"

# -- Copy example files to _static -------------------------------------------

examples_folder_path = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "examples")
sphinx_examples_folder_path = os.path.join(
    os.path.dirname(__file__), "examples")
if os.path.exists(sphinx_examples_folder_path):
    shutil.rmtree(sphinx_examples_folder_path)
os.mkdir(sphinx_examples_folder_path)
for file in os.listdir(examples_folder_path):
    if file.endswith(".ipynb"):
        shutil.copy(os.path.join(examples_folder_path, file),
                    os.path.join(sphinx_examples_folder_path, file))
