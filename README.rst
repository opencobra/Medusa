Medusa - Analysis of ensembles of metabolic network reconstructions
===================================================================

|Build Status| |PyPI|

Medusa is a tool for constraint-based reconstruction and analysis (COBRA) of ensembles. It builds on the cobrapy package (https://github.com/opencobra/cobrapy) by extending most single-model functionality to efficient ensemble-scale analysis. Additionally, Medusa provides novel functions for the analysis of ensembles.

Medusa is developed openly. We are currently porting code from previous/in-progress projects to generalize beyond the specific use cases that project code was originally developed for. See below for general release plans, and please contact us via email or by opening an issue if you are interested in a feature we have not yet implemented.


Installation
~~~~~~~~~~~~

Use pip to install medusa from PyPI::

    pip install medusa-cobra


Getting Started
~~~~~~~~~~~~~~~

Medusa can be loaded in python with::

    import medusa




Resources
~~~~~~~~~

We highly recommend familiarizing yourself with standard COBRA methods in the cobrapy package (https://github.com/opencobra/cobrapy). Some examples of ensemble analyses with genome-scale metabolic models can be found in our `recent preprint <https://doi.org/10.1101/460071>`_ and the accompanying `code repository <https://github.com/gregmedlock/ssl_ensembles>`_.

For reading on ensembles in metabolic modeling and systems biology, see:

1. Biggs MB, Papin JA (2017) Managing uncertainty in metabolic network structure and improving predictions using EnsembleFBA. PLoS Comput Biol 13(3): e1005413. https://doi.org/10.1371/journal.pcbi.1005413

2. Kuepfer, Lars, Matthias Peter, Uwe Sauer, and Jörg Stelling. 2007. “Ensemble Modeling for Analysis of Cell Signaling Dynamics.” Nature Biotechnology 25 (9): 1001–6. https://doi.org/10.1038/nbt1330


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
