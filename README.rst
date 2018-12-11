Medusa - Analysis of ensembles of metabolic network reconstructions
===================================================================

|Build Status| |PyPI|

Medusa is a tool for constraint-based reconstruction and analysis (COBRA) of ensembles. It builds on the cobrapy package (https://github.com/opencobra/cobrapy) by extending most single-model functionality to efficient ensemble-scale analysis. Additionally, Medusa provides novel functions for the analysis of ensembles.

Medusa is developed openly and we welcome contributions in the form of ideas and code. Please contact us by opening an issue if you are interested in a feature we have not yet implemented.


Installation
~~~~~~~~~~~~

Use pip to install medusa from PyPI (Note: we support Python 3 only; use Python 2 at your own peril)::

    pip install medusa-cobra


Getting Started
~~~~~~~~~~~~~~~

Medusa can be loaded in python with::

    import medusa


How to cite medusa
~~~~~~~~~~~~~~~~~~

When using medusa, please cite cobrapy, the paper describing the method you used that was implemented within cobrapy or medusa, and the medusa software itself. We are preparing a preprint describing the software--check back soon for the formal citation. Until then, In place of the software paper citation, please link to this github repository and note either the release (if installed via pip) or the commit (if working from the development branch between versions) in your publication.


Resources
~~~~~~~~~

We recommend familiarizing yourself with standard COBRA methods in the cobrapy package (https://github.com/opencobra/cobrapy). Some examples of ensemble analyses with genome-scale metabolic models can be found in our `recent preprint <https://doi.org/10.1101/460071>`_ and the accompanying `code repository <https://github.com/gregmedlock/ssl_ensembles>`_.

For reading on ensembles in metabolic modeling and systems biology, see our recent papers:

1. Biggs MB, Papin JA (2017) Managing uncertainty in metabolic network structure and improving predictions using EnsembleFBA. PLoS Comput Biol 13(3): e1005413. https://doi.org/10.1371/journal.pcbi.1005413

2. Medlock GL, Papin JA (2018) Guiding the refinement of biochemical knowledgebases with ensembles of metabolic networks and semi-supervised learning. bioRxiv, 2018. https://doi.org/10.1101/460071


Contact
~~~~~~~

We greatly appreciate feedback from current and prospective users from all backgrounds.
Feel free to post issues via github or contact us directly via email with any questions or thoughts.

Authors:
Greg Medlock (gmedlo [at] gmail [dot] com)

.. |Build Status| image:: https://api.travis-ci.org/gregmedlock/Medusa.svg?branch=master
   :target: https://travis-ci.org/gregmedlock/Medusa/
.. |PyPI| image:: https://badge.fury.io/py/medusa-cobra.svg
   :target: https://pypi.python.org/pypi/medusa-cobra
