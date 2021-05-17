<h1 align="center">
  <a href="https://github.com/meuleman/epilogos"><img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/logo.png" width="840"></a>
</h1>

---

<h2 align="center">
    Information-theoretic navigation of multi-tissue functional genomic annotations
</h2>

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos using Python 3.7. 
We provide a proof-of-principle dataset based on chromatin state calls from the EpiMap dataset (<a href="https://www.nature.com/articles/s41586-020-03145-z">Boix et al., Nature 2021</a>).

---

<br>


## Installation

To install Epilogos simply run the following command
```bash
$ pip install epilogos
```

Alternatively, install Epilogos directly from this Git repositoriy using
```bash
$ pip install git+https://github.com/meuleman/epilogos
```


## Prerequisites

To compute epilogos, you will need to have the following python libraries installed: 
[click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/stable/index.html), and [pandas](https://pandas.pydata.org/). 
These can be installed with one of the following commands.
```bash
$ pip install click numpy pandas pyranges scipy matplotlib statsmodels
```
or while in the epilogos directory
```bash
$ pip install -r requirements.txt
```

Additionally, it is recommended that python is updated to version 3.7 or later. 
In earlier python versions, `src/scores.py` may raise an OSError 16. 
It is worth noting that in our testing this error has not affected the results. 


## Running Epilogos

To be presented with basic documentation of arguments needed to run epilogos, simply run the command `epilogos --help` (More in-depth explanation is given [below](#command-line-options)).

By default, Epilogos assumes access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/). 
A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the `-l` flag to your command (e.g. `epilogos -l`).

---

<div align="center">This README is incomplete. For the full README, visit [our github](https://github.com/meuleman/epilogos)</div>


