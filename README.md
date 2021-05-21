[![Python package](https://github.com/meuleman/epilogos/actions/workflows/devel.yml/badge.svg?branch=main)](https://github.com/meuleman/epilogos/actions/workflows/devel.yml)

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

<p align="center">
    Created by: Wouter Meuleman, Jacob Quon, Alex Reynolds, and Eric Rynes
</p>

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#installation">Installation</a> •
    <a href="#prerequisites">Prerequisites</a> •
    <a href="#running-epilogos">Running Epilogos</a> •
    <a href="#slurm-examples">SLURM Examples</a> •
    <a href="#non-slurm-examples">Non-SLURM Examples</a> •
    <a href="#command-line-options">Command Line Options</a> •
    <a href="#pairwise-epilogos">Pairwise Epilogos</a>
  </h3>
</div>

---

<br>

<a name="installation"></a>

## Installation

To install Epilogos simply run the following command
```bash
$ pip install epilogos
```

Alternatively, install Epilogos directly from this Git repositoriy using
```bash
$ pip install git+https://github.com/meuleman/epilogos
```

<a name="prerequisites"></a>

## Prerequisites

To compute epilogos, you will need to have the following python libraries installed: [cython](https://cython.org/), [pyranges](https://github.com/biocore-ntnu/pyranges), [statsmodels](https://www.statsmodels.org/stable/index.html), [click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/stable/index.html), and [pandas](https://pandas.pydata.org/).
These can be installed with one of the following commands.
```bash
$ pip install cython; pip install click numpy pandas pyranges scipy matplotlib statsmodels
```
or while in the epilogos directory (we use cat and xargs to ensure installation order as pyranges is dependent on cython)
```bash
$ cat requirements.txt | xargs -n 1 -L 1 pip install
```

Additionally, it is recommended that python is updated to version 3.7 or later.
In earlier python versions, `src/scores.py` may raise an OSError 16.
It is worth noting that in our testing this error has not affected the results.

<a name="running-epilogos"></a>

## Running Epilogos

To be presented with basic documentation of arguments needed to run epilogos, simply run the command `epilogos --help` or `python -m epilogos --help` (More in-depth explanation is given [below](#command-line-options)).

By default, Epilogos assumes access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/).
A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the `-l` flag to your command (e.g. `epilogos -l`).

<a name="slurm-examples"></a>

## SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>Example data has been provided under <code>data/pyData/male/</code>.
The file, <code>epilogos_matrix_chrX.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome X.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -i data/pyData/male/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>scores_male_s1_epilogos_matrix_chrX.txt.gz</code> and <code>greatestHits_male_s1.txt</code> in <code>OUTPUTDIR</code></p>

<p>To customize your run of epilogos see the <a href="#command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>Before you can run Epilogos on your own data, you will need to complete two steps.</p>

<p>First, you will need to format your data such that Epilogos can parse it.
To assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information, run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome name
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```

<p>Second, you will need to create a state info file.
This is a tab separated file containing various information about each of the states in the chromatin state model.
We have provided some files already for common models in the <code>data/state_metadata/</code> directory.
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -i PATH_TO_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the same number of scores files as in your input directory in <code>OUTPUTDIR</code>.
Each of these files will be named <code>scores_*.txt.gz</code>, where 'scores_' is followed by the input directory name, the saliency metric, and the corresponding input file name (extensions removed).
Additionally, you will find a <code>greatestHits_*.txt</code> file which follows the same naming convention minus the input file name.</p>

<p>If you would like to visualize these results as seen on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<a name="non-slurm-examples"></a>

## Non-SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>Example data has been provided under <code>data/pyData/male/</code>.
The file, <code>epilogos_matrix_chrX.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome X.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -l -i data/pyData/male/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the file <code>scores_male_s1_epilogos_matrix_chrX.txt.gz</code> and <code>greatestHits_male_s1.txt</code> in <code>OUTPUTDIR</code></p>

<p>To customize your run of epilogos see the <a href="#command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>


<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>Before you can run Epilogos on your own data, you will need to complete two steps.</p>

<p>First, you will need to format your data such that Epilogos can parse it.
To assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```

<p>Second, you will need to create a state info file.
This is a tab separated file containing various information about each of the states in the chromatin state model.
We have provided some files already for common models in the <code>data/state_metadata/</code> directory.
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -l -i PATH_TO_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the same number of scores files as in your input directory in <code>OUTPUTDIR</code>.
Each of these files will be named <code>scores_*.txt.gz</code>, where 'scores_' is followed by the input directory name, the saliency metric, and the corresponding input file name (extensions removed).
Additionally, you will find a <code>greatestHits_*.txt</code> file which follows the same naming convention except for the input file name.</p>

<p>If you would like to visualize these results as seen on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>


<a name="command-line-options"></a>

## Command Line Options

<a name="mode"></a>
<details><summary><b> Mode [-m, --mode]</b></summary>
<p></p>
<p>Epilogos supports a single group and a paired group mode.
The single group mode identifies "salient" regions compared to a background of itself.
The paired group mode identifies regions which differ significantly between two groups of datasets.</p>

<p>
The argument to this flag either <code>single</code> or <code>paired</code> as the mode of operation, with <code>single</code> being the default.
</p>

```bash
e.g. $ epilogos -m single
```
</details>

<a name="command-line"></a>
<details><summary><b> Local Run [-l, --local]</b></summary>
<p></p>
<p>By default, Epilogos assumes access to a SLURM cluster. However, if you would like to run Epilogos locally in your terminal enable this flag.</p>

```bash
e.g. $ epilogos -l
```
</details>

<a name="input-directory"></a>
<details><summary><b> Input Directory [-i, --input-directory]</b></summary>
<p></p>
<p>Rather than just read in one input file, Epilogos reads the contents of an entire directory.
This allows the computation to be chunked and parallelized.
Note that the genome files in the directory <strong>MUST</strong> be split by chromosome.</p>

<p>
The argument to this flag is the path to the directory which contains the files to be read in.
Note that <strong>ALL</strong> files in this directory will be read in and errors may occur if other files are present.
</p>

```bash
e.g. $ epilogos -i data/pyData/male/
```

<p>Epilogos input data must be formatted specifically for Epilogos.
In order to help you create your own input data files, we have provided a script to transform chromHMM files into Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```
</details>

<a name="output-directory"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>
The output of Epilogos will vary depending on the number of input files present in the <a href="#input-directory">input directory [-i, --input-directory]</a>.
All resulting score files are gzipped .txt files named <code>scores_*.txt.gz</code>, where 'scores_' is followed by the input directory name, the saliency metric (e.g. S1), and the corresponding input file name (extensions removed).</p>

<p>Additionally, you will find a <code>greatestHits_*.txt</code> file which follows the same naming convention minus the input file name.
This file contains the top 1000 highest scoring regions (after merging directly adjacent high-scoring regions), with each row formatted as follows.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest scoring state
Column 5: Kullback-Leibler score
Column 6: Sign of the Kullback-Leibler score
```

<p>The argument to this flag is the path to the directory to which you would like to output.
Note that this may not be the same as the input directory.</p>

```bash
e.g. $ epilogos -o epilogosOutput/
```

</details>

<a name="state-info"></a>
<details><summary><b> State Info [-n, --state-info]</b></summary>
<p></p>
<p>The argument to this flag is a tab separated file specifying information about the state model being used.
This file must contain a header row with the exact names as shown below, and values should be formatting as shown below as well.
</p>

| zero_index | one_index | short_name | long_name | hex | rgba | color |
|------------|-----------|------------|-----------|-----|------|-------|
| 0 | 1 | TssA | Active TSS | #ff0000 | rgba(255,0,0,1) | Red |

<p>For more detail see <code>epilogos/data/state_metadata/README.md</code> or <code>epilogos/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv</code></p>

```bash
e.g. $ epilogos -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv
```

</details>

<a name="saliency"></a>
<details><summary><b> Saliency Level [-s, --saliency]</b></summary>
<p></p>
<p>Epilogos implements information-theoretic metrics to quantify saliency levels of datasets.
The <code>-s</code> flag to the coordination script allows one to choose one of three possible metrics:</p>

```
1. Metric S1, implementing a standard Kullback-Leibler relative entropy

2. Metric S2, implementing a version of S1 that additionally models label co-occurrence patterns

3. Metric S3, implementing a version of S2 that additionally models between-biosample similarities
```

<p>
Note that each increase in saliency level involves more computation and thus requires more time and computational power.
</p>

<p>
The argument to this flag must be an integer <code>1, 2, or 3</code>.
Note that Epilogos defaults to a saliency level of 1.
</p>

<p>Example:</p>

```bash
Saliency 1: $ epilogos -i data/pyData/male/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR

Saliency 2: $ epilogos -i data/pyData/male/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 2

Saliency 3: $ epilogos -i data/pyData/male/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 3
```

</details>

<a name="number-of-cores"></a>
<details><summary><b> Number of Cores [-c, --num-cores]</b></summary>
<p></p>
<p>Epilogos will always try and parallelize where it can.
Computation done on each input file is parallelized using python's <a href="https://docs.python.org/3/library/multiprocessing.html">multiprocessing library</a>.</p>

<p>The argument to this flag is an integer number of cores you would like to utilize to perform this multiprocessing.
Note that Epilogos defaults to using all available cores (equivalent to <code>-c 0</code>).</p>

```bash
e.g. $ epilogos -c 4
```
</details>

<a name="exit"></a>
<details><summary><b> Exit [-x, --exit]</b></summary>
<p></p>
<p>By default <code>src/computeEpilogosSlurm.py</code> prints progress updates to the console and only exits after it has completed all slurm jobs.
If you would like the program to instead exit when all jobs are submitted (allowing use of the terminal while the jobs are running), enable this flag.</p>

```bash
e.g. $ epilogos -x
```
</details>

<br>
<br>

<a name="pairwise-epilogos"></a>

<h1 align="center">
  <a href="https://github.com/meuleman/epilogos#pairwise-epilogos"><img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/pairwise_logo.png" width="840"></a>
</h1>

Pairwise Epilogos, like single-group Epilogos, is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations.
However, its role is to provide a structure by which to compare these genomic annotations accross different groups.

The software provided in this repository implements the methods underlying Pairwise Epilogos using only python.
We provide a proof-of-principle dataset based on chromatin state calls from the EpiMap dataset.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#running-epilogos-pairwise">Running Pairwise Epilogos</a> •
    <a href="#slurm-examples-pairwise">SLURM Examples</a> •
    <a href="#non-slurm-examples-pairwise">Non-SLURM Examples</a> •
    <a href="#command-line-options-pairwise">Command Line Options</a> •
    <a href="#visual-output">Visual Output</a> •
  </h3>
</div>

---

<br>

<a name="running-epilogos-pairwise"></a>

## Running Pairwise Epilogos

To be presented with minimal documentation of arguments needed to run epilogos, simply run the command `epilogos --help` or `python -m epilogos --help` (More in-depth explanation is given [below](#command-line-options-pairwise))

By default, Epilogos assumes access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/).
A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the `-l` flag to your command (e.g. `epilogos -l`).

<a name="slurm-examples-pairwise"></a>

## SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>Example data has been provided under <code>data/pyData/male/</code> and <code>data/pyData/female/</code>.
The files, both named <code>epilogos_matrix_chrX.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome X.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -a data/pyData/male/ -b data/pyData/female/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chrX.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>significantLoci_male_female_s1.txt</code>, and <code>greatestHits_male_female_s1.txt</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>.
For further explanation of the contents of these outputs see <a href="#output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>Before you can run Epilogos on your own data, you will need to complete two steps.</p>

<p>First, you will need to modify your data such that Epilogos can understand it.
To assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```

<p>Second, you will need to create a state info file.
This is a tab separated file containing various information about each of the states in the chromatin state model.
We have provided some files already for common models in the <code>data/state_metadata/</code> directory.
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -a PATH_TO_FIRST_INPUT_DIR -b PATH_TO_SECOND_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_*.txt.gz</code>, <code>pairwiseMetrics_*.txt.gz</code>, <code>significantLoci_*.txt</code>, and <code>greatestHits_*.txt</code> as well as the directory <code>manhattanPlots_*</code> in <code>OUTPUTDIR</code>.
Each of the wildcards will be replaced by a string containing the name of input directory one, the name of input directory two, the saliency metric, and the corresponding input file name when relevant (extensions removed)</p>

<p>If you would like to visualize these results as on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>


<a name="non-slurm-examples-pairwise"></a>

## Non-SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>Example data has been provided under <code>data/pyData/male/</code> and <code>data/pyData/female/</code>.
The files, both named <code>epilogos_matrix_chrX.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome X.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -l -a data/pyData/male/ -b data/pyData/female/ -n data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chrX.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>significantLoci_male_female_s1.txt</code>, and <code>greatestHits_male_female_s1.txt</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>.
For further explanation of the contents of these outputs see <a href="#output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>Before you can run Epilogos on your own data, you will need to complete two steps.</p>

<p>First, you will need to modify your data such that Epilogos can understand it.
To assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```

<p>Second, you will need to create a state info file.
This is a tab separated file containing various information about each of the states in the chromatin state model.
We have provided some files already for common models in the <code>data/state_metadata/</code> directory.
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -l -a PATH_TO_FIRST_INPUT_DIR -b PATH_TO_SECOND_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_*.txt.gz</code>, <code>pairwiseMetrics_*.txt.gz</code>, <code>significantLoci_*.txt</code>, and <code>greatestHits_*.txt</code> as well as the directory <code>manhattanPlots_*</code> in <code>OUTPUTDIR</code>.
Each of the wildcards will be replaced by a string containing the name of input directory one, the name of input directory two, the saliency metric, and the corresponding input file name when relevant (extensions removed)</p>

<p>If you would like to visualize these results as on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<a name="command-line-options-pairwise"></a>

## Command Line Options

Pairwise Epilogos has additional command line options beyond the options offered for [single group epilogos](#command-line-options).
These are outlined below.

<a name="directories-pairwise"></a>
<details><summary><b> Input Directories One and Two [-a, --directory-one] and [-b, --directory-two]</b></summary>
<p></p>
<p>Rather than just read in one input file, Epilogos reads the contents of an entire directory.
This allows the computation to be chunked and parallelized.
Note that the genome files in the directory <strong>MUST</strong> be split by chromosome.</p>

<p>In the paired group version of epilogos, the user must input two directories (one for each group).
Note that <strong>ALL</strong> files in this directory will be read in and errors may occur if other files are present.
Additionally, the files to compare within the directories must have the same name in both directories (e.g. chrX_male.txt and chrX_female.txt would need to be renamed to chrX.txt and chrX.txt, in their respective directories)</p>

<p>The arguments to these flags are the two directories you would like to read data from.</p>

```bash
e.g. $ epilogos -a data/pyData/male/ -b data/pydata/female/
```

<p>Epilogos input data must be formatted specifically for Epilogos.
To help you create your own input data files, we have provided a script to transform chromHMM files into Epilogos input files.
This can be found at <code>scripts/preprocess_data_ChromHMM.sh</code> (to get usage information run without arguments).
If you would prefer not to use the script, data is to be formatted as follows:</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: State data for epigenome 1
...
Column n: State data for epigenome n-3
```
</details>

<a name="output-directory-pairwise"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>The output of paired group Epilogos will vary depending on the number of input files present in the input directories <a href="#directories-pairwise">[-a, --directory-one]</a> or <a href="#directories-pairwise">[-b, --directory-two]</a>.
All score difference files will be gzipped txt files and of the format <code>pairwiseDelta_*.txt.gz</code> where 'pairwiseDelta_' is followed by the names of input directory one, input directory two, the saliency metric, and the name of the corresponding input file (extensions removed).
All other outputs follow this same name suffix format, with the exception that the corresponding input file is omitted in the case that the relevant file is a summary accross all input files.</p>

<p>The output directory will contain one <code>pairwiseMetrics_*.txt.gz</code> file which contains scores for all inputted data.
Each row is formatted as below.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Signed squared euclidean distance between the scores (sign determined by the higher signal between groups 1 and 2)
Column 6: P-Value of the distance
```

<p>The output directory will contain one <code>significantLoci_*.txt</code> file.
This file contains the all loci deemed to be significant below a threshold of 10% FDR after applying the Benjamini-Hochberg Procedure.
Each row is formatted as below.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
Column 7: P-Value of the distance
Column 8: Stars indicating multiple hypothesis adjusted p-value of distance ('***' at .01, '**' at .05, and '*' at .1)
```

<p>The output directory will contain one <code>greatestHits_*.txt</code> file.
This file contains the all significant loci with adjacent regions merged.
If there are less than 1000 significant loci, it takes the top 1000 highest distance regions and merges those.
Each row is formatted as below.</p> 

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
Column 7: P-Value of the distance
Column 8: Stars indicating multiple hypothesis adjusted p-value of distance ('***' stars at .01, '**' at .05, '*' at .1, '.' if not significant)
```

<p>The output directory will contain one <code>manhattanPlots_*</code> directory.
This directory will contain all the manhattan plots generated by pairwise epilogos.
These plots show the signed squared euclidean distances between groups 1 and 2 as well as the p-values of these distances.
There is one genome-wide plot generated and another plot generated for each chromosome.</p>

<p>Depending on the <a href="#diagnostic-figures">[-d, --diagnostic-figures]</a> flag the output directory may contain one <code>diagnosticFigures_*</code> directory.
This directory will contain figures showing the quality of the fit the null data and comparisons between the null and real data.</p>

<p>The argument to this flag is the path to the directory to which you would like to output.
Note that this <strong>CANNOT</strong> be the same as the input directory.</p>

```bash
e.g. $ epilogos -o epilgosOutput/
```
</details>

<a name="diagnostic-figures"></a>
<details><summary><b> Diagnostic Figures [-d, --diagnostic-figures]</b></summary>
<p></p>
<p>If this flag is enabled, Pairwise Epilogos will output diagnostic figures of the gennorm fit on the null data and comparisons between the null and real data.
These can be found in a sub-directory of the output directory named <code>diagnosticFigures_*</code> directory where 'diagnosticFigures_' is followed by the names of input directory one, input directory two, and the saliency metric.</p>

```bash
e.g. $ epilogos -d
```
</details>

<a name="num-trials"></a>
<details><summary><b> Number of Trials [-t, --num-trials]</b></summary>
<p></p>
<p>In order to save time when fitting in paired group Epilogos, random samplings of the null data are fit -t times with the median fit being used.</p>

<p>The argument to this flag is the number of random samplings fit. Epilogos defaults to 101</P>

```bash
e.g. $ epilogos -t 1001
```
</details>

<a name="sampling-size"></a>
<details><summary><b> Sampling Size [-z, --sampling-size]</b></summary>
<p></p>
<p>In order to save time when fitting in paired group Epilogos, random samplings of the null data are fit -t times with the median fit being used.</p>

<p>The argument to this flag is the size  of random samplings fit. Epilogos defaults to 100,000</P>

```bash
e.g. $ epilogos -t 10000
```
</details>

<a name="quiescent-state"></a>
<details><summary><b> Quiescent State [-q, --quiescent-state]</b></summary>
<p></p>
<p>As a proxy for disregarding unmappable genomic regions, epilogos does not consider regions entirely in a "Quiescent" state in the construction of the empirical null distribution used to determine p-values.</p>

<p>The argument to this flag is the 1-index of Quiescent state.
Note that by default epilogos assumes the Quiescent state is the last state in the provided chromatin state model (in an 18 state model <code>$ epilogos -q 18</code> is equivalent to not passing the <code>-q</code> flag.</p>

```bash
e.g. $ epilogos -q 18
```

<p>If you prefer epilogos not filtering out any data, use 0 as the argument for this flag</p>

```bash
e.g. $ epilogos -q 0
```
</details>

<a name="visual-output"></a>

## Visual Output

Unlike single group Epilogos, pairwise Epilogos has a visual componenet in the form of Manhattan plots.
Located in the <code>manhattanPlots_*</code> output directory, these plots offer users a way to visually locate differences between two groups.
These plots color in any points deemed to be significant above a threshold of 10% FDR according to the colors specifed in the [state info tsv](#state-info) provided by the user.
The background color represents the level of significance.
White means insignificant and the shades of gray mean significant at 10%, 5%, and 1% FDR (with darker being more significant).
Additionally, points have varying opacity determined by the ratio of their distance to the most extreme distance.

The example plot below shows a genome-wide Manhattan plot generated by a pairwise Epilogos run of material from male donors versus female donors in the [EpiMap dataset](https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486).
This view makes it immediately clear that the vast majority of chromatin state differences between male and female biosamples is on the X chromosome.

<h1 align="center">
  <img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/manhattan_male_female_genome.png" width="840">
</h1>

Considering that these differences are primarily present on chromosome X, we now can take a closer look at only chromosome X.
With newfound granularity, we can now begin to identify individual genes where the differences are most pronounced.
We can clearly see differences around 73 Mbp (XIST locus) and around 131 Mbp (FIRRE locus).
Furthermore, the directionality of these points towards female donor biosamples indicates that these differences were driven by higher Epilogos scores in female donor biosamples.
This is consistent with our knowledge of these loci playing crucial roles in the female-specific X-inactivation process.

<h1 align="center">
  <img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/manhattan_male_female_chrX.png" width="840">
</h1>




