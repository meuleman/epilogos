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
    <a href="#plot-region">Plot Region</a> •
    <a href="#pairwise-epilogos">Pairwise Epilogos</a> •
    <a href="#similarity-search">Similarity Search</a>
  </h3>
</div>

---

<br>

<a name="installation"></a>

## Installation

Although not required, it is good practice to create a virtual environment in which
specific versions of Python and its libraries are installed.
This can be done using `conda`, for instance as such:
```bash
$ conda init bash  ## only needed upon first use of conda. Restart shell after this.
$ conda create -n epilogos python=3.9
$ conda activate epilogos
```

[comment]: <> ($ conda install -c anaconda libopenblas)


To install Epilogos simply run the following command
```bash
$ pip install epilogos
```

Alternatively, install Epilogos directly from this Git repository using
```bash
$ pip install git+https://github.com/meuleman/epilogos
```

<a name="prerequisites"></a>

## Prerequisites

To compute epilogos, you will need to have the following python libraries installed: [statsmodels](https://www.statsmodels.org/stable/index.html), [click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](https://matplotlib.org/stable/index.html), and [pandas](https://pandas.pydata.org/).
In case the abovementioned commands not automatically and correctly take care of this, the libraries can be installed with one of the following commands.
```bash
$ pip install 'click>=8.1.3,<9.0.0' 'numpy>=1.19.2,<2.0.0' 'pandas>=1.1.3,<2.0.0' 'scipy>=1.5.2,<2.0.0' 'matplotlib>=3.3.2,<4.0.0' 'statsmodels>=0.12.0,<1.0.0' 'scikit-learn>=1.1.2,<2.0.0' 'pysam>=0.19.1,<1.0.0' 'filter-regions @ git+https://github.com/alexpreynolds/filter-regions@8c2ef14dec35b7a4b6092fb2afe2eac409d58275'
```
or while in the epilogos directory
```bash
$ pip install -r requirements.txt
```

Additionally, it is recommended that python is updated to version 3.7 or later.
In earlier python versions, `epilogos/scores.py` may raise an OSError 16.
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

<p>If you cloned this git repository, example data has been provided under <code>data/pyData/male/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The file, <code>epilogos_matrix_chr1.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -i data/pyData/male/ -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code> and <code>regionsOfInterest_male_s1.txt</code> in <code>OUTPUTDIR</code></p>

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
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-j, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -i PATH/TO/INPUT_DIR -j PATH/TO/STATE_INFO_TSV -o PATH/TO/OUTPUT_DIR
```

<p>Upon completion of the run, you should see the same number of scores files as in your input directory in <code>OUTPUTDIR</code>.
Each of these files will be named <code>scores_*.txt.gz</code>, where 'scores_' is followed by the input directory name, the saliency metric, and the corresponding input file name (extensions removed).
Additionally, you will find a <code>regionsOfInterest_*.txt</code> file which follows the same naming convention minus the input file name.</p>

<p>If you would like to visualize these results as seen on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<a name="non-slurm-examples"></a>

## Non-SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>If you cloned this git repository, example data has been provided under <code>data/pyData/male/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The file, <code>epilogos_matrix_chr1.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -l -i data/pyData/male/ -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the file <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code> and <code>regionsOfInterest_male_s1.txt</code> in <code>OUTPUTDIR</code></p>

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
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-j, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -l -i PATH/TO/INPUT_DIR -j PATH/TO/STATE_INFO_TSV -o PATH/TO/OUTPUT_DIR
```

<p>Upon completion of the run, you should see the same number of scores files as in your input directory in <code>OUTPUTDIR</code>.
Each of these files will be named <code>scores_*.txt.gz</code>, where 'scores_' is followed by the input directory name, the saliency metric, and the corresponding input file name (extensions removed).
Additionally, you will find a <code>regionsOfInterest_*.txt</code> file which follows the same naming convention except for the input file name.</p>

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

<p>Additionally, you will find a <code>regionsOfInterest_*.txt</code> file which follows the same naming convention minus the input file name.
This file contains the top 100 regions of interest. Each row is formatted as follows.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest scoring state
Column 5: Sum of the Kullback-Leibler scores
Column 6: Sign of the sum of Kullback-Leibler scores
```

<p>The argument to this flag is the path to the directory to which you would like to output.
Note that this may not be the same as the input directory.</p>

```bash
e.g. $ epilogos -o epilogosOutput/
```

</details>

<a name="state-info"></a>
<details><summary><b> State Info [-j, --state-info]</b></summary>
<p></p>
<p>The argument to this flag is a tab separated file specifying information about the state model being used.
This file must contain a header row with the exact names as shown below, and values should be formatting as shown below as well.
</p>

| zero_index | one_index | short_name | long_name | hex | rgba | color |
|------------|-----------|------------|-----------|-----|------|-------|
| 0 | 1 | TssA | Active TSS | #ff0000 | rgba(255,0,0,1) | Red |

<p>For more detail see <code>epilogos/data/state_metadata/README.md</code> or <code>epilogos/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv</code></p>

```bash
e.g. $ epilogos -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv
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

<p>The argument to this flag is an integer number of cores you would like to utilize to perform this multiprocessing. Should you like to use all available cores, input 0 to the option (i.e. <code>-c 0</code>). Epilogos defaults to using 1 core.</p>

```bash
e.g. $ epilogos -c 4
```
</details>

<a name="exit"></a>
<details><summary><b> Exit [-x, --exit]</b></summary>
<p></p>
<p>By default <code>epilogos/run.py</code> prints progress updates to the console and only exits after it has completed all slurm jobs.
If you would like the program to instead exit when all jobs are submitted (allowing use of the terminal while the jobs are running), enable this flag.</p>

```bash
e.g. $ epilogos -x
```
</details>

<a name="version"></a>
<details><summary><b> Version [-v, --version]</b></summary>
<p></p>
<p>If this flag is enabled epilogos will print the installed version number and exit</p>

```bash
e.g. $ epilogos -v
```
</details>

<a name="partition"></a>
<details><summary><b> Partition [-p, --partition]</b></summary>
<p></p>
<p>By default <code>epilogos/run.py</code>, uses the default partition designated by the system administrator to submit SLURM jobs.
Use this flag if you would like to specify the partition for SLURM resource allocation.</p>

<p>The argument to the flag is the name of the partition you would like the SLURM jobs to run on</p>

```bash
e.g. $ epilogos -p queue1
```
</details>

<a name="roi-width"></a>
<details><summary><b> Region of Interest Width [-w, --roi-width]</b></summary>
<p></p>
<p>This flag controls the size of the regions of interest in <code>regionsOfInterest_*.txt</code>.</p>

<p>The argument to this flag is size of the regions of interest in bins. Epilogos defaults to 50 for <code>single</code> mode and 125 for <code>paired</code> mode (results in 10kb and 25kb regions respectively with standard 200bp bins).</p>

```bash
e.g. $ epilogos -w 10
```
</details>

<a name="file-tag"></a>
<details><summary><b> File Tag [-f, --file-tag]</b></summary>
<p></p>
<p>This flag controls the string appended to each of the output files</p>

<p>The argument to this flag is string you want appended to each of the output files. Epilogos defaults to using <code>INPUT-DIR_SALIENCY</code> (where INPUT-DIR is replaced with the name of the directory containing the epilogos input files and SALIENCY is replaced with the chosen saliency metric</p>

```bash
e.g. $ epilogos -f male_s2
```
</details>

<a name="exp-freq-mem"></a>
<details><summary><b> Expected Frequency Calculation Memory Allocation [--exp-freq-mem]</b></summary>
<p></p>
<p>This flag controls the amount of memory (in MB) assigned to each of the expected frequency calculation slurm jobs. This flag has no effect when running Epilogos using the <a href="#command-line">[-l, --local]</a> flag.</p>

<p>The argument to this flag is the MB of memory desired for each of the expected frequency calculation slurm jobs. Epilogos defaults to 16000</p>

```bash
e.g. $ epilogos --exp-freq-mem 8000
```
</details>

<a name="exp-comb-mem"></a>
<details><summary><b> Expected Frequency Combination Memory Allocation [--exp-comb-mem]</b></summary>
<p></p>
<p>This flag controls the amount of memory (in MB) assigned to the expected frequency combination slurm job. This flag has no effect when running Epilogos using the <a href="#command-line">[-l, --local]</a> flag.</p>

<p>The argument to this flag is the MB of memory desired for the expected frequency combination slurm job. Epilogos defaults to 8000</p>

```bash
e.g. $ epilogos --exp-comb-mem 4000
```
</details>

<a name="score-mem"></a>
<details><summary><b> Score Calculation Memory Allocation [--score-mem]</b></summary>
<p></p>
<p>This flag controls the amount of memory (in MB) assigned to each of the score calculation slurm jobs. This flag has no effect when running Epilogos using the <a href="#command-line">[-l, --local]</a> flag.</p>

<p>The argument to this flag is the MB of memory desired for each of the score calculation slurm jobs. Epilogos defaults to 16000</p>

```bash
e.g. $ epilogos --score-mem 8000
```
</details>

<a name="roi-mem"></a>
<details><summary><b> Region of Interest Calculation Memory Allocation [--roi-mem]</b></summary>
<p></p>
<p>This flag controls the amount of memory (in MB) assigned to the region of interest calculation slurm job. This flag has no effect when running Epilogos using the <a href="#command-line">[-l, --local]</a> flag.</p>

<p>The argument to this flag is the MB of memory desired for the region of interest calculatino slurm job. Epilogos defaults to 20000 for <code>single</code> mode and 100000 for <code>paired</code> mode.</p>

```bash
e.g. $ epilogos --roi-mem 10000
```
</details>


<a name="plot-region"></a>

<h1 align="center">
  <a href="https://github.com/meuleman/epilogos"><img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/plotregion_logo.png" width="840"></a>
</h1>

Though Epilogos reduces the dimensionality of datasets, its output can still present challenges with regards to exploration and interpretation. To this end, we recommend using higlass to visualize results. This can be an arduous process and thus we include the `plotregion` command-line utility. Plot Region provides a simple way to visualize small amounts of epilogos data in a style like the <a href="https://epilogos.altius.org/">Epilogos web browser</a>.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#running-plot-region">Running Plot Region</a> •
    <a href="#plot-region-examples">Examples</a> •
    <a href="#command-line-options-plot-region">Command Line Options</a>
  </h3>
</div>

---

<br>

<a name="running-plot-region"></a>

## Running Plot Region

To be presented with minimal documentation of arguments needed to run plotregion, simply run the command `plotregion --help` or `python -m plotregion --help` (More in-depth explanation is given [below](#command-line-options-plot-region)).

Plot Region runs locally in the command line. It is lightweight enough to not present computational issues at alow scale. If you are plotting many regions and wish to use a computational cluster, we recommend [SLURM](https://slurm.schedmd.com/)</code>

<a name="plot-region-examples"></a>

## Plot Region Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>
<p>The Plot Region commandline interface requires as input four flags: <a href="regions-plot-region">[-r, --regions]</a> <a href="scores-plot-region">[-s, --scores-file]</a>, <a href="state-info-plot-region">[-j, --state-info]</a>, and <a href="output-directory-plot-region">[-o, --output-directory]</a>.</p>

<p>If you cloned this git repository, example data has been provided under <code>data/plotregion/male/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The file, named <code>scores_male_s1_matrix_chr1.txt.gz</code>, contains epilogos scores for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.

The data consist of epilogos scores calculated on chromosome 1 of the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.

Additionally, a tab-separated bed file of 5 regions has been provided as region input at <code>data/plotregion/male/regionsOfInterest_male_s1_chr1.bed</code></p>

<p>To compute Plot Region results for this sample data run the following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ plotregion -r data/plotregion/male/regionsOfInterest_male_s1_chr1.bed -s data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see five files following the <code>epilogos_region_CHR_START_END.pdf</code> naming convention (one for each region in <code>data/plotregion/male/</code>). These files contain the epilogos scores plotted for each of the input regions. For further explanation of the contents of these outputs see <a href="#output-directory-plot-region">Output Directory [-o, --output-directory]</a></p>

</details>

<details><summary><b> Running Plot Region with your own data</b></summary>
<p></p>

<p>Before you can run Plot Region on your own data, you will first need an epilogos scores file. When epilogos is run, it outputs scores split by chromosome. Because Plot Region can only read in one file, if you want to run similarity search across the whole genome, you will have to combine these files into one singular scores file. This file can have chromosomes sorted by genomic (i.e. chr9 before chr12) or lexicographic (i.e. chr12 before chr9) order. We recommend using the either of following commands:</p>

<p><strong>Genomic:</strong></p>

```bash
$ prefix="PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_matrix"; suffix="txt.gz"; paths=""; for chr in GENOMIC_ORDER; do chr="chr${chr}"; path="${prefix}_${chr}.${suffix}"; paths="${paths} ${path}"; done; cat ${paths} > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt.gz
```
<p>Where is GENOMIC_ORDER is replaced with the names of the relevant chromosomes in order separated by spaces. (e.g. <code>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y</code> or <code>`seq 1 22` X Y</code> for humans)

<p><strong>Lexicographic (requires bedops):</strong></p>

```bash
$ zcat PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_chr*.txt.gz | sort-bed - > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt; gzip PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt
```

<br>

<p>Once you have your desired scores file you can generate Plot Region results with the following command:</p>

```bash
$ plotregion -r CHR:START-END -s PATH/TO/EPILOGOS_SCORES_FILE -j PATH/TO/METADATA_FILE -o PATH/TO/BUILD_OUTPUT_DIR
```

<p>Upon completion of the run, you should see files following the <code>epilogos_region_CHR_START_END.pdf</code> naming convention (one for each region input with -r). These files contain the epilogos scores plotted for each of the input regions. For further explanation of the contents of these outputs see <a href="#output-directory-plot-region">Output Directory [-o, --output-directory]</a></p>

</details>


<a name="command-line-options-plot-region"></a>

## Plot Region Command Line Options

<a name="regions-plot-region"></a>
<details><summary><b> Regions [-r, --regions]</b></summary>
<p></p>
<p>The Plot Region commandline interface requires as input at least one region. The <code>-r</code> flag can handle both single & multi region inputs. If using a single region, the argument to the flag should be the region coordinates formatted as chr:start-end. If using multiple regions, the argument to the flag should be the path to a tab-separated bed file.</p>

```bash
e.g. $ plotregion -r chr4:93305800-93330800
```

or

```bash
e.g. $ plotregion -r /PATH/TO/input_regions.bed
```
</details>

<a name="scores-plot-region"></a>
<details><summary><b> Scores File [-s, --scores-fles]</b></summary>
<p></p>
<p>The Plot Region commandline interface requires as input an Epilogos scores file to read region information from.When epilogos is run, it outputs scores split by chromosome. Because Plot Region can only read in one file, if you want to run similarity search across the whole genome, you will have to combine these files into one singular scores file. This file can have chromosomes sorted by genomic (i.e. chr9 before chr12) or lexicographic (i.e. chr12 before chr9) order. We recommend using the either of following commands:</p>

<p><strong>Genomic:</strong></p>

```bash
$ prefix="PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_matrix"; suffix="txt.gz"; paths=""; for chr in GENOMIC_ORDER; do chr="chr${chr}"; path="${prefix}_${chr}.${suffix}"; paths="${paths} ${path}"; done; cat ${paths} > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt.gz
```
<p>Where is GENOMIC_ORDER is replaced with the names of the relevant chromosomes in order separated by spaces. (e.g. <code>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y</code> or <code>`seq 1 22` X Y</code> for humans)

<p><strong>Lexicographic (requires bedops):</strong></p>

```bash
$ zcat PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_chr*.txt.gz | sort-bed - > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt; gzip PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt
```

<br>

<p>The argument to this flag should be the path to a previously run Epilogos scores file (see <a href="#slurm-examples">Epilogos</a> for details).</p>

```bash
e.g. $ plotregion -s EPILOGOS_OUTPUTDIR/scores.txt.gz
```
</details>

<a name="state-info-plot-region"></a>
<details><summary><b> State Info [-j, --state-info]</b></summary>
<p></p>
<p>The argument to this flag is a tab separated file specifying information about the state model being used.
This file must contain a header row with the exact names as shown below, and values should be formatting as shown below as well.
</p>

| zero_index | one_index | short_name | long_name | hex | rgba | color |
|------------|-----------|------------|-----------|-----|------|-------|
| 0 | 1 | TssA | Active TSS | #ff0000 | rgba(255,0,0,1) | Red |

<p>For more detail see <code>epilogos/data/state_metadata/README.md</code> or <code>epilogos/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv</code></p>

```bash
e.g. $ plotregion -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv
```
</details>

<a name="output-directory-plot-region"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>Plot Region outputs a matplotlib figure for each of the input regions. Each of these figures is saved to a file named <code>epilogos_region_CHR_START_END.pdf</code> in the user-defined output directory. The scores are drawn top to bottom in descending order of score individually for each bin. This allows the user to quickly identify the most prevalent state.</p>

<p>The argument to this flag is the path to the desired output directory.</p>

```bash
e.g. $ plotregion -o PATH/TO/OUTPUT_DIRECTORY
```
</details>

<a name="ylims-plot-region"></a>
<details><summary><b> Individual Y-Limits [-y, --individual-ylims]</b></summary>
<p></p>
<p>By default, Plot Region draws all output figures on the same y-scale. The y-max and y-min are calculated to be the minimum necessary to fit all figures. This flag causes Plot Region to abandon this behavior and instead plot each region on individually sized axes.</p>

<p>When this flag is used, <code>plotregion</code> graphs each region on individually sized axes. This flag only has an effect when <a href="#regions-plot-region">graphing multiple regions</a> (using a tab-separated bed file).</p>

```bash
e.g. $ plotregion -y
```
</details>

<a name="file-format-plot-region"></a>
<details><summary><b> File Format [-f, --file-format]</b></summary>
<p></p>
<p>By default, Plot Region saves all output figures as pdfs. Use this flag to save them in a different format (any format supported by <a href="https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html">matplotlib savefig</a>).</p>

<p>The argument to this flag is the desired file format.</p>

```bash
e.g. $ plotregion -f png
```
</details>

<br>
<br>

<a name="pairwise-epilogos"></a>

<h1 align="center">
  <a href="https://github.com/meuleman/epilogos#pairwise-epilogos"><img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/pairwise_logo.png" width="840"></a>
</h1>

Pairwise Epilogos, like single-group Epilogos, is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations.
However, its role is to provide a structure by which to compare these genomic annotations across different groups.

The software provided in this repository implements the methods underlying Pairwise Epilogos using only python.
We provide a proof-of-principle dataset based on chromatin state calls from the EpiMap dataset.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#running-epilogos-pairwise">Running Pairwise Epilogos</a> •
    <a href="#slurm-examples-pairwise">SLURM Examples</a> •
    <a href="#non-slurm-examples-pairwise">Non-SLURM Examples</a> •
    <a href="#command-line-options-pairwise">Command Line Options</a> •
    <a href="#visual-output">Visual Output</a>
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

<p>If you cloned this git repository, example data has been provided under <code>data/pyData/male/</code> and <code>data/pyData/female/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The files, both named <code>epilogos_matrix_chr1.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -a data/pyData/male/ -b data/pyData/female/ -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>significantLoci_male_female_s1.txt</code>, and <code>regionsOfInterest_male_female_s1.txt</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>.
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
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-j, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -a PATH/TO/FIRST_INPUT_DIR -b PATH/TO/SECOND_INPUT_DIR -j PATH/TO/STATE_INFO_TSV -o PATH/TO/OUTPUT_DIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_*.txt.gz</code>, <code>pairwiseMetrics_*.txt.gz</code>, <code>significantLoci_*.txt</code>, and <code>regionsOfInterest_*.txt</code> as well as the directory <code>manhattanPlots_*</code> in <code>OUTPUTDIR</code>.
Each of the wildcards will be replaced by a string containing the name of input directory one, the name of input directory two, the saliency metric, and the corresponding input file name when relevant (extensions removed)</p>

<p>If you would like to visualize these results as on <a href="https://epilogos.altius.org">epilogos.altius.org</a>, we recommend using higlass.</p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>


<a name="non-slurm-examples-pairwise"></a>

## Non-SLURM Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>If you cloned this git repository, example data has been provided under <code>data/pyData/male/</code> and <code>data/pyData/female/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The files, both named <code>epilogos_matrix_chr1.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.
The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -l -a data/pyData/male/ -b data/pyData/female/ -j data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>significantLoci_male_female_s1.txt</code>, and <code>regionsOfInterest_male_female_s1.txt</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>.
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
For more information on the structure of these files see <code>data/state_metadata/README.txt</code> or <a href="#state-info">State Info [-j, --state-info]</a></p>

<p>Once you have completed these two steps, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -l -a PATH/TO/FIRST_INPUT_DIR -b PATH/TO/SECOND_INPUT_DIR -j PATH/TO/STATE_INFO_TSV -o PATH/TO/OUTPUT_DIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_*.txt.gz</code>, <code>pairwiseMetrics_*.txt.gz</code>, <code>significantLoci_*.txt</code>, and <code>regionsOfInterest_*.txt</code> as well as the directory <code>manhattanPlots_*</code> in <code>OUTPUTDIR</code>.
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
All other outputs follow this same name suffix format, with the exception that the corresponding input file is omitted in the case that the relevant file is a summary across all input files.</p>

<p>The output directory will contain one <code>pairwiseMetrics_*.txt.gz</code> file which contains scores for all inputted data.
Each row is formatted as below.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
```

<p>The output directory will contain one <code>regionsOfInterest_*.txt</code> file.
This file contains the top 100 regions of interest.
Each row is formatted as below.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
Column 7: Z-score of the distance
Column 8: Stars indicating Z-score of distance ('***' stars at >=3, '**' at >=2, '*' at >=1, '.' if <1)
```

<p>The output directory will contain one <code>manhattanPlots_*</code> directory.
This directory will contain all the manhattan plots generated by pairwise epilogos.
These plots show the signed squared euclidean distances between groups 1 and 2 as well as the z-scores of these distances.
There is one genome-wide plot generated and another plot generated for each chromosome.</p>

<p>The argument to this flag is the path to the directory to which you would like to output.
Note that this <strong>CANNOT</strong> be the same as the input directory.</p>

```bash
e.g. $ epilogos -o epilogosOutput/
```

<p>All summary files will have format changes if the <a href="#null-distribution">[-n, --null-distribution]</a> flag is used. Go to the <a href="#null-distribution">[-n, --null-distribution]</a> section to view these changes<p>

</details>

<a name="diagnostic-figures"></a>
<details><summary><b> Diagnostic Figures [-d, --diagnostic-figures]</b></summary>
<p></p>
<p>If this flag is enabled, Pairwise Epilogos will output diagnostic figures of the gennorm fit on the null data and comparisons between the null and real data.
These can be found in a sub-directory of the output directory named <code>diagnosticFigures_*</code> directory where 'diagnosticFigures_' is followed by the names of input directory one, input directory two, and the saliency metric.</p>

<p>NOTE: This flag only caused output in conjunction with <a href="#null-distribution">[-n, --null-distribution]</a></p>

```bash
e.g. $ epilogos -d
```
</details>

<a name="null-distribution"></a>
<details><summary><b> Null Distribution [-n, --null-distribution]</b></summary>
<p></p>
<p>If you would like the Epilogos scores to be given p-values and multiple hypothesis corrected p-values enable this flag. These p-values will be output in <code>pairwiseMetrics_*.txt</code> and will be used for <code>regionsOfInterest_*.txt</code> and <code>significantLoci_*.txt</code> generation. They will also be used for manhattan plot coloring.</p>

```bash
e.g. $ epilogos -n
```

<p>Using [-n, --null-distribution] will change the files in the <a href="#output-directory-pairwise">output directory</a>. The <code>pairwiseMetrics_*.txt.gz</code> file which contains scores for all inputted data will contain p-values (details below).</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
Column 7: P-Value of the distance
Column 8: Benjamini–Hochberg adjusted P-Value of the distance
```

<p>The output directory will now additionally contain one <code>significantLoci_*.txt</code> file.
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
Column 8: Benjamini–Hochberg adjusted P-Value of the distance
Column 9: Stars indicating multiple hypothesis adjusted p-value of distance ('***' at .01, '**' at .05, and '*' at .1)
```

<p>The <code>regionsOfInterest_*.txt</code> file will now contain the top 100 or fewer significant regions.
Each row is formatted as below.</p>

```
Column 1: Chromosome
Column 2: Start coordinate
Column 3: End coordinate
Column 4: Name of the largest difference state
Column 5: Squared euclidean distance between the scores
Column 6: Direction of the distance (sign determined by the higher signal between groups 1 and 2)
Column 7: P-Value of the distance
Column 8: Benjamini–Hochberg adjusted P-Value of the distance
Column 9: Stars indicating multiple hypothesis adjusted p-value of distance ('***' stars at .01, '**' at .05, '*' at .1, '.' if not significant)
```

<p>The manhattan plots in the <code>manhattanPlots_*</code> sub-directory will now display p-values rather than z-scores on the second y-axis.</p>

<p>Depending on the <a href="#diagnostic-figures">[-d, --diagnostic-figures]</a> flag the output directory may contain one <code>diagnosticFigures_*</code> directory.
This directory will contain figures showing the quality of the fit the null data and comparisons between the null and real data.</p>
</details>

<a name="num-trials"></a>
<details><summary><b> Number of Trials [-t, --num-trials]</b></summary>
<p></p>
<p>Only to be used in conjunction with <a href="#null-distribution">[-n, --null-distribution]</a>. In order to save time when fitting the null distribution in paired group Epilogos, random samplings of the null data are fit -t times with the median fit being used.</p>

<p>The argument to this flag is the number of random samplings fit. Epilogos defaults to 101</P>

```bash
e.g. $ epilogos -t 1001
```
</details>

<a name="sampling-size"></a>
<details><summary><b> Sampling Size [-z, --sampling-size]</b></summary>
<p></p>
<p>Only to be used in conjunction with <a href="#null-distribution">[-n, --null-distribution]</a>. In order to save time when fitting the null distribution in paired group Epilogos, random samplings of the null data are fit -t times with the median fit being used.</p>

<p>The argument to this flag is the size  of random samplings fit. Epilogos defaults to 100,000</P>

```bash
e.g. $ epilogos -z 10000
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

<a name="group-size"></a>
<details><summary><b> Group Size [-g, --group-size]</b></summary>
<p></p>
<p>In Pairwise Epilogos, all inputted data is considered when generating scores. If you would only like a certain amount of the epigenomes considered, you can do so with this flag. The epigenomes considered are chosen at random for each bin (that is the epigenomes chosen in bin 1 could be different than those chosen in bin 2). This flag equalizes the size of the two groups to the inputted value.</p>

<p>The argument to this flag is the number of epigenomes per group you would like considered in calculations. Epilogos defaults to all (equivalent to <code>-g -1</code>)</P>

```bash
e.g. $ epilogos -g 30
```
</details>

<a name="visual-output"></a>

## Visual Output

Unlike single group Epilogos, pairwise Epilogos has a visual component in the form of Manhattan-like plots.
Located in the <code>manhattanPlots_*</code> output directory, these plots offer users a way to visually locate differences between two groups of epigenomes.
Points are colored in case their scores have a z-score greater than 1. When using <a href="#null-distribution">[-n, --null-distribution]</a> points are colored in case they are found to be significantly different between groups, exceeding a threshold of 10% FDR.
The colors are as specified in the [state info tsv](#state-info) provided by the user.
Additionally, points have varying opacity determined by the ratio of their distance to the most extreme distance.
The background color represents the z-score.
White indicates less than 10 and the shades of gray indicate greater than 10, 20, and 30 (with darker being a higher z-score).
When using <a href="#null-distribution">[-n, --null-distribution]</a>, the background color represents level of signficicance.
White indicates insignificant and the shades of gray indicate significance at 10%, 5%, and 1% FDR (with darker being more significant).

The example plot below shows a genome-wide Manhattan plot generated by a pairwise Epilogos run of material from male donors versus female donors in the [EpiMap dataset](https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486).
This view makes it immediately clear that the vast majority of chromatin state differences between male and female biosamples is on the X chromosome.

<h1 align="center">
  <img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/manhattan_male_female_genome.png" width="840">
</h1>

Considering that these differences are primarily present on chromosome X, we now can take a closer look at only chromosome X.
With newfound granularity, we can begin to identify individual genes where differences are most pronounced.
We can clearly see differences around 73 Mbp (XIST locus) and around 131 Mbp (FIRRE locus).
Furthermore, the directionality of these points towards female donor biosamples indicates that these differences were driven by higher Epilogos scores in female donor biosamples.
This is consistent with our knowledge of these loci playing crucial roles in the female-specific X-inactivation process.

<h1 align="center">
  <img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/manhattan_male_female_chrX.png" width="840">
</h1>

<br>
<br>

<a name="similarity-search"></a>

<h1 align="center">
  <a href="https://github.com/meuleman/epilogos#similarity-search"><img src="https://raw.githubusercontent.com/meuleman/epilogos/main/data/simsearch_logo.png" width="840"></a>
</h1>

Epilogos similarity search is a method to enhance the browsing experience of genome-wide epilogos datasets. Similarity search reduces an Epilogos dataset into a set non-overlapping regions and performs a nearest neighbor search to identify the highest similarity regions for each region.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#running-similarity-search">Running Similarity Search</a> •
    <a href="#similarity-search-example">Examples</a> •
    <a href="#command-line-options-similarity-search">Command Line Options</a>
  </h3>
</div>

---

<br>

<a name="running-similarity-search"></a>

## Running Similarity Search

To be presented with minimal documentation of arguments needed to run similarity search, simply run the command `simsearch --help` or `python -m simsearch --help` (More in-depth explanation is given [below](#command-line-options-similarity-search))

By default, similarity search runs locally in the command line. Should you desire to use a computational cluster managed by [SLURM](https://slurm.schedmd.com/) (as used in Epilogos score calculation), a example bash script has been provided in <code>epilogos/bin/sim_search_multizoom.sh</code>

<a name="similarity-search-example"></a>

## Similarity Search Examples

<details><summary><b> Minimal example on provided example data</b></summary>
<p></p>

<p>If you cloned this git repository, example data has been provided under <code>data/simsearch/male/</code>. Otherwise it is available for download using the script in <code>bin/download_example_data.sh</code>. The script uses <a href="https://curl.se/">cURL</a> to download neccessary files and places them in a file hierarchy generated within the current directory.
The file, named <code>scores_male_s1_matrix_chr1.txt.gz</code>, contains epilogos scores for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1.

The data consist of epilogos scores calculated on chromosome 1 of the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">EpiMap dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute similarity search results for this sample data run the following command within the <code>epilogos/</code> directory (replacing <code>BUILD_OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz -o BUILD_OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>simsearch_cube.npz</code>, <code>simsearch_knn.npz</code>, <code>simsearch.bed.gz</code>, and <code>simsearch.bed.gz.tbi</code> in <code>BUILD_OUTPUTDIR</code>.
For further explanation of the contents of these outputs see <a href="#output-similarity-search-build">Build Output Directory [-o, --output-directory]</a></p>

<p>To query the Similarity Search results for a given region run the following command (replacing <code>CHR, START, & END</code> with the region coordinates, <code>BUILD_OUTPUTDIR</code> with the build step output directory, and <code>QUERY_OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ simsearch -q CHR:START-END -m BUILD_OUTPUTDIR/simsearch.bed.gz -o QUERY_OUTPUTDIR
```

<p>Upon completion of the run, you should see the file <code>similarity_search_region_CHR_START_END_recs.bed</code> in <code>QUERY_OUTPUTDIR</code>.
For further explanation of the contents of these outputs see <a href="#output-similarity-search-query">Query Output Directory [-o, --output-directory]</a></p>

<p>To customize your run of Similarity Search see the <a href="#command-line-options-similarity-search">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Similarity Search with your own data</b></summary>
<p></p>

<p>Before you can run Similarity Search on your own data, you will first need an epilogos scores file. When epilogos is run, it outputs scores split by chromosome. Because Similarity Search can only read in one file, if you want to run similarity search across the whole genome, you will have to combine these files into one singular scores file. This file can have chromosomes sorted by genomic (i.e. chr9 before chr12) or lexicographic (i.e. chr12 before chr9) order. We recommend using the either of following commands:</p>

<p><strong>Genomic:</strong></p>

```bash
$ prefix="PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_matrix"; suffix="txt.gz"; paths=""; for chr in GENOMIC_ORDER; do chr="chr${chr}"; path="${prefix}_${chr}.${suffix}"; paths="${paths} ${path}"; done; cat ${paths} > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt.gz
```
<p>Where is GENOMIC_ORDER is replaced with the names of the relevant chromosomes in order separated by spaces. (e.g. <code>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y</code> or <code>`seq 1 22` X Y</code> for humans)

<p><strong>Lexicographic (requires bedops):</strong></p>

```bash
$ zcat PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_chr*.txt.gz | sort-bed - > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt; gzip PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt
```

<br>

<p>Once you have your desired scores file you can generate Similarity Search results with the following command:</p>

```bash
$ simsearch -b -s PATH/TO/EPILOGOS_SCORES_FILE -o PATH/TO/BUILD_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the files <code>simsearch_cube.npz</code>, <code>simsearch_knn.npz</code>, <code>simsearch.bed.gz</code>, and <code>simsearch.bed.gz.tbi</code> in <code>PATH/TO/BUILD_OUTPUT_DIR</code>.
For further explanation of the contents of these outputs see <a href="#output-similarity-search-build">Build Output Directory [-o, --output-directory]</a></p>

<p>To query these Similarity Search results for a given region run the following command (replacing <code>CHR, START, & END</code> with the region coordinates, <code>PATH/TO/BUILD_OUTPUT_DIR</code> with the build step output directory, and <code>PATH/TO/BUILD_OUTPUT_DIR</code> with the output directory of your choice).</p>

```bash
$ simsearch -q CHR:START-END -m PATH/TO/BUILD_OUTPUT_DIR/simsearch.bed.gz -o PATH/TO/BUILD_OUTPUT_DIR
```

<p>Upon completion of the run, you should see the file <code>similarity_search_region_CHR_START_END_recs.bed</code> in <code>PATH/TO/BUILD_OUTPUT_DIR</code>.
For further explanation of the contents of these outputs see <a href="#output-similarity-search-query">Query Output Directory [-o, --output-directory]</a></p>

<p>To further customize your run of epilogos see the <a href="#command-line-options-similarity-search">Command Line Options</a> of the <code>README</code></p>

</details>

<a name="command-line-options-similarity-search"></a>

## Build Command Line Options

<a name="build-similarity-search"></a>
<details><summary><b> Build Similarity Search Results [-b, --build]</b></summary>
<p></p>
<p>The Similarity Search commandline interface has two modes: build and <a href="query-similarity-search">query</a>. Build mode takes in a single epilogos scores file, finds regions of a specified size, and finds the best matches for each of these regions.</p>

<p>When <code>-b</code> is enabled, run Similarity Search in build mode</p>

```bash
e.g. $ simsearch -b
```
</details>

<a name="scores-similarity-search"></a>
<details><summary><b> Scores Path [-s, --scores-path]</b></summary>
<p></p>
<p>In <a href="build-similarity-search">[-b, --build]</a> mode, Similarity Search takes as input a single epilogos scores file. When epilogos is run, it outputs scores split by chromosome. Because Similarity Search can only read in one file, if you want to run similarity search across the whole genome, you will have to combine these files into one singular scores file. This file can have chromosomes sorted by genomic (i.e. chr9 before chr12) or lexicographic (i.e. chr12 before chr9) order. We recommend using the either of following commands:</p>

<p><strong>Genomic:</strong></p>

```bash
$ prefix="PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_matrix"; suffix="txt.gz"; paths=""; for chr in GENOMIC_ORDER; do chr="chr${chr}"; path="${prefix}_${chr}.${suffix}"; paths="${paths} ${path}"; done; cat ${paths} > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt.gz
```
<p>Where is GENOMIC_ORDER is replaced with the names of the relevant chromosomes in order separated by spaces. (e.g. <code>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y</code> or <code>`seq 1 22` X Y</code> for humans)

<p><strong>Lexicographic (requires bedops):</strong></p>

```bash
$ zcat PATH/TO/EPILOGOS_OUTPUT_DIR/scores_*_chr*.txt.gz | sort-bed - > PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt; gzip PATH/TO/EPILOGOS_OUTPUT_DIR/scores.txt
```

<br>

<p>The argument to this flag is the path to the file you would like to read data from.</p>

```bash
e.g. $ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz
```
</details>

<a name="output-similarity-search-build"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>In <a href="build-similarity-search">[-b, --build]</a> mode, Similarity Search will output 4 files: <code>simsearch_cube.npz</code>, <code>simsearch_knn.npz</code>, <code>simsearch.bed.gz</code>, and <code>simsearch.bed.gz.tbi</code></p>

<p><code>simsearch_cube.npz</code> is a zipped archive of 2 numpy arrays: <code>scores</code> and <code>coords</code>. <code>scores</code> contains the epilogos scores across each region present in the similarity search output (regions are chosen by the mean-max algorithm as outlined in the <a href="https://github.com/alexpreynolds/filter-regions">filter-regions package</a>). Scores are condensed into 25 blocks regardless of chosen window size. Scores in each of these blocks consist of the epilogos scores of the bin with the highest sum of scores in each block. <code>coords</code> contains the coordinates for each region, and is sorted in the same order as <code>scores</code>.</p>

<p><code>simsearch_knn.npz</code> is a zipped archive of 3 numpy arrays: <code>arr</code>, <code>idx</code>, and <code>dist</code>. <code>arr</code> contains the coordinates of each region followed by its top 100 matches. <code>idx</code> contains the index of each region in the same order as <code>arr</code> followed by the the indices of its top 100 matches. <code>dist</code> contains the distance of each region to itself (0) followed by the distances to each of its top 100 matches.</p>

<p><code>simsearch.bed.gz</code> is a gzipped bed file containing the coordinates of each region followed by the coordinates of its top 100 matches.</p>

<p><code>simsearch.bed.gz.tbi</code> is a <a href="http://www.htslib.org/doc/tabix.html">tabix</a> equivalent of <code>simsearch.bed.gz</code>.

<p>The argument to this flag is the path to the directory to which you would like to output.</p>

```bash
e.g. $ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz -o BUILD_OUTPUTDIR/
```
</details>

<a name="window-similarity-search"></a>
<details><summary><b> Window Size [-w, --window-bp]</b></summary>
<p></p>
<p>Similarity Search uses the mean-max algorithm outlined in the <a href="https://github.com/alexpreynolds/filter-regions">filter-regions package</a> to divide the genome into non-overlapping regions of a user-specified size (prioritizing regions based on information content). Scores within these regions are condensed into 25 blocks regardless of chosen window size. Scores in each of these blocks consist of the epilogos scores of the bin with the highest sum of scores in each block. </p>

<p>Only a set number of window size options are available for Similarity Search. When run on a 200bp bin dataset, these are 5000, 10000, 25000, 50000, 75000, and 100000. When run on a 20bp bin dataset, these are 500, 1000, 2500, 5000, 7500, and 10000</p>

<p>The argument to this flag is the desired window size in bp (default=25000(200bp dataset)/2500(20bp dataset)).</p>

```bash
e.g. $ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz -o OUTPUT_DIR -w 10000
```
</details>

<a name="jobs-similarity-search"></a>
<details><summary><b> Number of Jobs [-j, --num-jobs]</b></summary>
<p></p>
<p>Similarity Search uses <a href="https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.NearestNeighbors.html">sklearn's nearest neighbors function</a> to determine the best recommendations for each region. This flag specifies the number of parallel jobs to run for this nearest neighbors search.</p>

<p>The argument to this flag is the number of parallel jobs to run for the nearest neighbors search (default=8).</p>

```bash
e.g. $ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz -o OUTPUT_DIR -j 4
```
</details>

<a name="num-neighbors-similarity-search"></a>
<details><summary><b> Number of Matches [-n, --num-matches]</b></summary>
<p></p>
<p>Similarity Search uses <a href="https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.NearestNeighbors.html">sklearn's nearest neighbors function</a> to determine the best matches for each region. This flag specifies the number of desired neighbors (i.e. matches) for this nearest neighbors search. NOTE that by definition of sklearn, the first neighbor is always the query region (i.e. if you want 100 matches, the argument must be 101).</p>

<p>The argument to this flag is the number desired matches for each region (note the query region itself is counted as a matches) (default=101).</p>

```bash
e.g. $ simsearch -b -s data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz -o OUTPUT_DIR -n 51
```
</details>

## Query Command Line Options

<a name="query-similarity-search"></a>
<details><summary><b> Query Similarity Search Results [-q, --query]</b></summary>
<p></p>
<p>The Similarity Search commandline interface has two modes: <a href="build-similarity-search">build</a> and query. Query mode takes in a region query(ies) and a previously calculated <code>simsearch.bed.gz</code> matches file and outputs the top 100 matches for the query(ies).</p>

<p>The <code>-q</code> flag can handle both single & multi region queries. If querying a single region, the argument to the flag should be the region coordinates formatted as chr:start-end. If querying multiple regions, the argument to the flag should be the path to a tab-separated bed file.</p>

```bash
e.g. $ simsearch -q chr4:93305800-93330800 -m /PATH/TO/simsearch.bed.gz
```

or

```bash
e.g. $ simsearch -q /PATH/TO/query_regions.bed -m /PATH/TO/simsearch.bed.gz
```
</details>

<a name="matches-similarity-search"></a>
<details><summary><b> Matches File [-m, --matches-file]</b></summary>
<p></p>
<p>In <a href="query-similarity-search">[-q, --query]</a> mode, Similarity Search takes as argument a pre-built <code>simsearch.bed.gz</code> matches file (see <a href="output-similarity-search-build">build output</a> for details). This file is queried to find the top 100 matches for each of the query regions.</a>

```bash
e.g. $ simsearch -q chr4:93305800-93330800 -m BUILD_OUTPUTDIR/simsearch.bed.gz
```
</details>

<a name="output-similarity-search-query"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p></p>
<p>In <a href="query-similarity-search">[-q, --query]</a> mode, Similarity Search will output 1 file for each query region: <code>similarity_search_region_CHR_START_END_recs.bed</code> (where <code>CHR</code>, <code>START</code>, & <code>END</code> are replaced with the query coordinates</p>

<p><code>similarity_search_region_CHR_START_END_recs.bed</code> is a tab-separated bed file containing the coordinates for each of the top 100 matches for the query region in sorted order.</p>

<p>The argument to this flag is the path to the directory to which you would like to output.</p>

```bash
e.g. $ simsearch -q chr4:93305800-93330800 -m BUILD_OUTPUTDIR/simsearch.bed.gz -o QUERY_OUTPUTDIR/
```
</details>

<br>
<br>

## Development

To modify a development version of `epilogos`, first set up a virtual environment via e.g. [`conda`](#installation).
After activating the environment and installing [dependencies](#prerequisites), install `epilogos` in editable mode:

```bash
$ pip install -e .
```

### Big Sur

If you are using Mac OS X 11 (Big Sur) or later, it may be necessary to first install OpenBLAS, before installing Python dependencies that require it (such as `scipy`).

This can be done via [Homebrew](https://brew.sh/) and setting environment variables to point to relevant library and header files:

```bash
$ brew install openblas
$ export SYSTEM_VERSION_COMPAT=1
$ export LAPACK=/usr/local/opt/openblas/lib
$ export LAPACK_SRC=/usr/local/opt/openblas/include
$ export BLAS=/usr/local/opt/openblas/lib
```

Then run `pip install -e .`, as before.


