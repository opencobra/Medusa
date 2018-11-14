Medusa - Contribution guide
===================================================================


Installation
~~~~~~~~~~~~

Create a fork of the `development` branch and clone to your computer. We recommend working within a virtualenv, and installing your clone locally by navigating to the root of the project (e.g. /Medusa/) and running::

    python setup.py development

You may need to rerun this command after making changes to medusa for those changes to take effect.


Getting your contribution into medusa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For now, we are following the cobrapy contribution guidelines, which can be found [here](https://github.com/opencobra/cobrapy/blob/devel/.github/CONTRIBUTING.rst).


Documentation
~~~~~~~~~~~~~

Documentation is managed with [Readthedocs](https://docs.readthedocs.io/en/latest/) using [Sphinx](http://www.sphinx-doc.org/en/master/index.html). Readthedocs manages automated deployment as html to https://medusa.readthedocs.io/en/latest/, and we primarily use Jupyter notebooks for each section of the documentation. Sphinx and Readthedocs essentially handle the conversion of our collection of Jupyter notebooks to a nice HTML document.

TODO: add instructions for updating the documentation and changing the build process here.


Deployment
~~~~~~~~~~

Releases of medusa are organized through the python package index (PyPI), which as accessible through `pip`. The current URL for the project is https://pypi.org/project/medusa-cobra/.

TODO: add instructions for updating the release on `pip` here.

