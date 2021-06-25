.. image:: https://raw.githubusercontent.com/meuleman/epilogos/main/data/logo.png
  :width: 840
  :alt: Epilogos

------


Information-theoretic navigation of multi-tissue functional genomic annotations
===============================================================================

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos using Python 3.7. We provide a proof-of-principle dataset based on chromatin state calls from the EpiMap dataset (`Boix et al., Nature 2021 <https://www.nature.com/articles/s41586-020-03145-z>`_).

Created by: Wouter Meuleman, Jacob Quon, Alex Reynolds, and Eric Rynes

------


Installation
============

Although not required, it is good practice to create a virtual environment in which
specific versions of Python and its libraries are installed.
This can be done using ``conda``, for instance as such:

.. code-block:: bash

    $ conda init bash  ## only needed upon first use of conda. Restart shell after this.
    $ conda create -n epilogos python=3.8
    $ conda activate epilogos

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

    $ pip install 'cython>=0.29.23,<1.0.0'; pip install 'click>=7.1.2,<8.0.0' 'numpy>=1.19.2,<2.0.0' 'pandas>=1.1.3,<2.0.0' 'pyranges>=0.0.97,<1.0.0' 'scipy>=1.5.2,<2.0.0' 'matplotlib>=3.3.2,<4.0.0' 'statsmodels>=0.12.0,<1.0.0'

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
