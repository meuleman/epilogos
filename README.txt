.. image:: https://raw.githubusercontent.com/meuleman/epilogos/main/data/logo.png
  :width: 840
  :alt: Epilogos

------


Information-theoretic navigation of multi-tissue functional genomic annotations
===============================================================================

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos using Python 3.7. 

------


Installation
============

To install Epilogos simply run the following command

.. code-block:: bash

    $ pip install epilogos

Alternatively, install Epilogos directly from the Git repositority using

.. code-block:: bash

    $ pip install git+https://github.com/meuleman/epilogos


Prerequisites
=============

To compute epilogos, you will need to have the following python libraries installed: `cython <https://cython.org/>`_, `pyranges <https://github.com/biocore-ntnu/pyranges>`_, `statsmodels <https://www.statsmodels.org/stable/index.html>`_, `click <https://click.palletsprojects.com/en/7.x/>`_, `numpy <https://numpy.org/>`_, `scipy <https://www.scipy.org/>`_, `matplotlib <https://matplotlib.org/stable/index.html>`_, and `pandas <https://pandas.pydata.org/>`_. These can be installed with one of the following commands

.. code-block:: bash

    $ pip install cython; pip install click numpy pandas pyranges scipy matplotlib statsmodels

or while in the epilogos directory (we use cat and xargs to ensure installation order as pyranges is dependent on cython)

.. code-block:: bash

    $ cat requirements.txt | xargs -n 1 -L 1 pip install


Additionally, it is recommended that python is updated to version 3.7 or later. In earlier python versions, ``src/scores.py`` may raise an OSError 16. It is worth noting that in our testing this error has not affected the results. 


Running Epilogos
================

To be presented with basic documentation of arguments needed to run epilogos, simply run the command ``epilogos --help`` or ``python -m epilogos --help`` (More in-depth explanation is given on the `github README <https://github.com/meuleman/epilogos)>`_).

By default, Epilogos assumes access to a computational cluster managed by `SLURM <https://slurm.schedmd.com/>`_. 
A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the ``-l`` flag to your command (e.g. ``epilogos -l``).

--------------

For a more extensive set of documentation, please refer to `our github <https://github.com/meuleman/epilogos>`_.
