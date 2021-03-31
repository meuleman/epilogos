<div align="center">
  <pre>
                  d8b 888                                     
                  Y8P 888                                     
                      888                                     
 .d88b.  88888b.  888 888  .d88b.   .d88b.   .d88b.  .d8888b  
d8P  Y8b 888 "88b 888 888 d88""88b d88P"88b d88""88b 88K      
88888888 888  888 888 888 888  888 888  888 888  888 "Y8888b. 
Y8b.     888 d88P 888 888 Y88..88P Y88b 888 Y88..88P      X88 
 "Y8888  88888P"  888 888  "Y88P"   "Y88888  "Y88P"   88888P' 
         888                            888                   
         888                       Y8b d88P                   
         888                        "Y88P"                    
  </pre>
</div>

---

<h2 align="center">
    INSERT TAGLINE HERE
</h2>

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos using only python. We provide a proof-of-principle dataset based on chromatin state calls from the BOIX dataset.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#prerequisites">Prerequisites</a> •
    <a href="#installation">Installation</a> •
    <a href="#running-epilogos">Running Epilogos</a> •
    <a href="#slurm-examples">SLURM Examples</a> •
    <a href="#non-slurm-examples">Non-SLURM Examples</a> •
    <a href="#command-line-options">Command Line Options</a> •
    <a href="#pairwise-epilogos">Pairwise Epilogos</a>
  </h3>
</div>

---

<br>


<a name="prerequisites"></a>

## Prerequisites

To compute epilogos, you will need to have the following python libraries: [click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), and [pandas](https://pandas.pydata.org/). These can be installed with one of the following commands.
```bash
$ pip install click, numpy, pandas
```
or while in the epilogos directory
```bash
$ pip install requirements.txt
```

Additionally, it is recommended that python is updated to version 3.7 or later. In earlier versions, `src/computeEpilogosScores.py` may raise an OSError 16. It is worth noting that in our testing this error has not affected the results. 

<a name="installation"></a>

## Installation

To install Epilogos simply run the following command
```bash
$ pip install epilogos
```

<a name="running-epilogos"></a>

## Running Epilogos

To be presented with minimal documentation of arguments needed to run epilogos, simply run the command `epilogos --help` (More in-depth explanation is given [below](#command-line-options))

By default, Epilogos assumes access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/). A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the `-l` flag to your command (e.g. `epilogos -l`).

<a name="slurm-examples"></a>

## SLURM Examples

<details><summary><b> Basic run</b></summary>
<p></p>

<p>Sample data has been provided under <code>~/epilogos/data/pyData/male/</code>. The file, <code>epilogos_matrix_chr1.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1. The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">BOIX dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>exp_freq_male_s1.npy</code> and <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code> in <code>OUTPUTDIR</code></p>

<p>To customize your run of epilogos see the <a href="command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>


<details><summary><b> Changing the saliency metric</b></summary>
<p></p>

<p>We will use the same sample data as for the basic run above (<code>~/epilogos/data/pyData/male/</code>).</p>

<p>To compute epilogos for this sample data run one of the following commands (depending on the desired saliency metric) within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
Saliency 1: $ epilogos -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR

Saliency 2: $ epilogos -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 2

Saliency 3: $ epilogos -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 3
```

<p>Upon completion of the run, you should see the files <code>exp_freq_male.npy</code> and <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code>, <code>scores_male_s2_epilogos_matrix_chr1.txt.gz</code>, or <code>scores_male_s3_epilogos_matrix_chr1.txt.gz</code> depending on the saliency metric you chose. in <code>OUTPUTDIR</code></p>

<p>To further customize your run of epilogos see the <a href="command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>In order to run Epilogos on your own data, you will need to do two things.</p>

<p>First, you will need to modify your data such that Epilogos can understand it. In order to assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files. This can be found at <code>~/epilogos/scripts/preprocess_data_ChromHMM.sh</code></p>

<p>Second, you will need to create a state info file. This is a tab separated file which tells epilogos various information about each of the states in the state model. We have provided some files already for common models in the <code>~/epilogos/data/state_metadata/</code> directory. For more information on the structure of these files see <code>~/epilogos/data/state_metadata/README.txt</code> or <a href="state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two things, you can run epilogos with the following command:</p>

```bash
$ epilogos -i PATH_TO_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

</details>

<a name="non-slurm-examples"></a>

## Non-SLURM Examples

<details><summary><b> Basic run</b></summary>
<p></p>

<p>Sample data has been provided under <code>~/epilogos/data/pyData/male/</code>. The file, <code>epilogos_matrix_chr1.txt.gz</code>, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1. The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">BOIX dataset</a> and contains only those epigenomes which are tagged <code>Male</code> under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -l -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>exp_freq_male_s1.npy</code> and <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code> in <code>OUTPUTDIR</code></p>

<p>To customize your run of epilogos see the <a href="command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>


<details><summary><b> Changing the saliency metric</b></summary>
<p></p>

<p>We will use the same sample data as for the basic run above (<code>~/epilogos/data/pyData/male/</code>).</p>

<p>To compute epilogos for this sample data run one of the following commands (depending on the desired saliency metric) within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
Saliency 1: $ epilogos -l -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR

Saliency 2: $ epilogos -l -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 2

Saliency 3: $ epilogos -l -i ./data/pyData/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 3
```

<p>Upon completion of the run, you should see the files <code>exp_freq_male.npy</code> and <code>scores_male_s1_epilogos_matrix_chr1.txt.gz</code>, <code>scores_male_s2_epilogos_matrix_chr1.txt.gz</code>, or <code>scores_male_s3_epilogos_matrix_chr1.txt.gz</code> depending on the saliency metric you chose. in <code>OUTPUTDIR</code></p>

<p>To further customize your run of epilogos see the <a href="command-line-options">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>In order to run Epilogos on your own data, you will need to do two things.</p>

<p>First, you will need to modify your data such that Epilogos can understand it. In order to assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files. This can be found at <code>~/epilogos/scripts/preprocess_data_ChromHMM.sh</code></p>

<p>Second, you will need to create a state info file. This is a tab separated file which tells epilogos various information about each of the states in the state model. We have provided some files already for common models in the <code>~/epilogos/data/state_metadata/</code> directory. For more information on the structure of these files see <code>~/epilogos/data/state_metadata/README.txt</code> or <a href="state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two things, you can run epilogos with the following command:</p>

```bash
$ epilogos -l -i PATH_TO_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

</details>


<a name="command-line-options"></a>

## Command Line Options

<a name="mode"></a>
<details><summary><b> Mode [-m, --mode]</b></summary>
<p></p>
<p>Epilogos supports a single group and a paired group mode. The single group mode finds interesting regions compared to a background of itself. Whereas the paired group mode finds regions which differ significantly between the two groups.</p>

<p>
The argument to this flag either <code>single</code> or <code>paired</code> as the mode of operation, with <code>single</code> being the default.
</p>
</details>

<a name="command-line"></a>
<details><summary><b> Command Line [-l, --command-line]</b></summary>
<p></p>
<p>By default, Epilogos assumes use of a SLURM cluster. However, if you would like to run Epilogos directly in the command line enable this flag</p>
</details>

<a name="input-directory"></a>
<details><summary><b> Input Directory [-i, --input-directory]</b></summary>
<p></p>
<p>Rather than just read in one input file, Epilogos reads the contents of an entire directory. This allows the computation to be chunked and parallelized. Note that the genome files in the directory <strong>MUST</strong> be split by chromosome. Input files should be formatted such that the first three columns are the chromosome, bin start, and bin end respectively with the rest of the columns containing state data.</p>

<p>
The argument to this flag is the path to the directory which contains the files to be read in. Note that <strong>ALL</strong> files in this directory will be read in and errors may occur if other files are present.
</p>
</details>

<a name="output-directory"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>
The output of Epilogos will vary depending on the number of input files present in the input directory <a href="input-directory">[-i, --input-directory]</a>. All scores files will be gzipped txt files and of the format <code>scores_{}_s$_[].txt.gz</code> where <code>{}</code> is replaced with the input directory name, <code>$</code> is replaced with the saliency level, and <code>[]</code> is replaced with the name of the corresponding input file (extensions removed).</p>
<p>
The argument to this flag is the path to the directory to which you would like to output. Note that this may not be the same as the input directory.</p>
</details>

<a name="state-info"></a>
<details><summary><b> State Info [-n, --state-info]</b></summary>
<p></p>
<p>The argument to this flag is a tab separated file specifying information about the state model being used. Example below (for more detail see <code>epilogos/data/state_metadata/README.md</code> or <code>epilogos/data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv</code></p>
| zero_index | one_index | short_name | long_name | hex | rgba | color |
|------------|-----------|------------|-----------|-----|------|-------|
| 0 | 1 | TssA | Active TSS | #ff0000 | rgba(255,0,0,1) | Red |
<p>
Note that tsv must contain a header row with the exact names above and that values within the table should follow the same formatting as above.
</p>
</details>

<a name="saliency"></a>
<details><summary><b> Saliency Level [-s, --saliency]</b></summary>
<p></p>
<p>Epilogos implements information-theoretic metrics to quantify saliency levels of datasets. The <code>-l</code> flag to the coordination script allows one to choose one of three possible metrics:</p>

```
1. Metric S1, implementing a standard Kullback-Leibler relative entropy

2. Metric S2, implementing a version of S1 that additionally models label co-occurrence patterns

3. Metric S3, implementing a version of S2 that additionally models between-biosample similarities
```

<p>
Note that each increase in saliency level involves much more computation and thus each increase requires more time and computational power.
</p>

<p>
The arguement to this flag must be an integer <code>1, 2, or 3</code>. Note that Epilogos defaults to a saliency of 1.
</p>
</details>

<a name="number-of-cores"></a>
<details><summary><b> Number of Cores [-c, --num-cores]</b></summary>
<p></p>
<p>Epilogos will always try and parallelize where it can. Computation done on each input file is parallelized using python's <a href="https://docs.python.org/3/library/multiprocessing.html">multiprocessing library</a>.</p>

<p>
The argument to this flag is an integer number of cores you would like to utilize to perform this multiprocessing. Note that Epilogos defaults to using all available cores (equivalent to <code>-c 0</code>).</p>
</details>

<a name="exit"></a>
<details><summary><b> Exit [-x, --exit]</b></summary>
<p></p>
<p>By default <code>src/computeEpilogosSlurm.py</code> only exits after it has completed all slurm jobs and prints progress updates to the console. If you would like the program to instead exit when all jobs are submitted (allowing use of the terminal while the jobs are running), enable this flag.</p>
</details>

<br>
<br>

<a name="pairwise-epilogos"></a>

<h2 align="center">
    Pairwise Epilogos
</h2>

Pairwise Epilogos, like Epilogos, is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations. However, its role is to provide a structure by which to compare these genomic annotations accross different groups.

The software provided in this repository implements the methods underlying Pairwise Epilogos using only python. We provide a proof-of-principle dataset based on chromatin state calls from the BOIX dataset.

---

<div align="center"><a name="menu"></a>
  <h3>
    <a href="#prerequisites-pairwise">Prerequisites</a> •
    <a href="#installation-pairwise">Installation</a> •
    <a href="#running-epilogos-pairwise">Running Pairwise Epilogos</a> •
    <a href="#slurm-examples-pairwise">SLURM Examples</a> •
    <a href="#non-slurm-examples-pairwise">Non-SLURM Examples</a> •
    <a href="#command-line-options-pairwise">Command Line Options</a> •
  </h3>
</div>

---

<br>

<a name="prerequisites-pairwise"></a>

## Prerequisites

To compute epilogos, you will need to have the following python libraries: [click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [matplotlib](#https://matplotlib.org/stable/index.html), and [pandas](https://pandas.pydata.org/). These can be installed with one of the following commands.
```bash
$ pip install click, numpy, pandas, scipy
```
or while in the epilogos directory
```bash
$ pip install requirements.txt
```

Additionally, it is recommended that python is updated to version 3.7 or later. In earlier versions, `src/computeEpilogosScores.py` may raise an OSError 16. It is worth noting that in our testing this error has not affected the results. 

<a name="installation-pairwise"></a>

## Installation

Installation of Pairwise Epilogos follows the same procedure as single group Epilogos:
```bash
$ pip install epilogos
```

<a name="running-epilogos-pairwise"></a>

## Running Pairwise Epilogos

To be presented with minimal documentation of arguments needed to run epilogos, simply run the command `epilogos --help` (More in-depth explanation is given [below](#command-line-options-pairwise))

By default, Epilogos assumes access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/). A version of epilogos has been created for those without access to a SLURM cluster and can be run by using the `-l` flag to your command (e.g. `epilogos -l`).

<a name="slurm-examples-pairwise"></a>

## SLURM Examples

<details><summary><b> Basic run</b></summary>
<p></p>

<p>Sample data has been provided under <code>~/epilogos/data/pyData/male/</code> and <code>~/epilogos/data/pyData/female/</code>. The files, both named <code>epilogos_matrix_chr1.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1. The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">BOIX dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -a ./data/pyData/male/ -b ./data/pyData/female/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>greatestHits_male_female_s1.bed</code>, and <code>exp_freq_male_female_s1.npy</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>. For further explanation of the contents of these outputs see <a href="output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To customize your run of epilogos see the <a href="command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>


<details><summary><b> Changing the saliency metric</b></summary>
<p></p>

<p>We will use the same sample data as for the basic run above (<code>~/epilogos/data/pyData/male/</code> and <code>~/epilogos/data/pyData/female/</code>).</p>

<p>To compute epilogos for this sample data run one of the following commands (depending on the desired saliency metric) within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
Saliency 1: $ epilogos -m paired -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR

Saliency 2: $ epilogos -m paired -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 2

Saliency 3: $ epilogos -m paired -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 3
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s$_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s$.txt.gz</code>, <code>greatestHits_male_female_s$.bed</code>, and <code>exp_freq_male_female_s$.npy</code> as well as the directory <code>manhattanPlots_male_female_s$</code> (where $ is replaced with the saliency metric you chose) in <code>OUTPUTDIR</code>. For further explanation of the contents of these outputs see <a href="output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To further customize your run of epilogos see the <a href="command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>In order to run Epilogos on your own data, you will need to do two things.</p>

<p>First, you will need to modify your data such that Epilogos can understand it. In order to assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files. This can be found at <code>~/epilogos/scripts/preprocess_data_ChromHMM.sh</code></p>

<p>Second, you will need to create a state info file. This is a tab separated file which tells epilogos various information about each of the states in the state model. We have provided some files already for common models in the <code>~/epilogos/data/state_metadata/</code> directory. For more information on the structure of these files see <code>~/epilogos/data/state_metadata/README.txt</code> or <a href="state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two things, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -a PATH_TO_FIRST_INPUT_DIR -b PATH_TO_SECOND_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

</details>


<a name="non-slurm-examples-pairwise"></a>

## Non-SLURM Examples

<details><summary><b> Basic run</b></summary>
<p></p>

<p>Sample data has been provided under <code>~/epilogos/data/pyData/male/</code> and <code>~/epilogos/data/pyData/female/</code>. The files, both named <code>epilogos_matrix_chr1.txt.gz</code>, contain chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1. The data was pulled from the <a href="https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486">BOIX dataset</a> and contains only those epigenomes which are tagged <code>Male</code> or <code>Female</code> respectively under the <code>Sex</code> column.</p>

<p>To compute epilogos (using the S1 saliency metric) for this sample data run following command within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
$ epilogos -m paired -l -a ./data/pyData/male/ -b ./data/pyData/female/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s1_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s1.txt.gz</code>, <code>greatestHits_male_female_s1.bed</code>, and <code>exp_freq_male_female_s1.npy</code> as well as the directory <code>manhattanPlots_male_female_s1</code> in <code>OUTPUTDIR</code>. For further explanation of the contents of these outputs see <a href="output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To customize your run of epilogos see the <a href="command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>


<details><summary><b> Changing the saliency metric</b></summary>
<p></p>

<p>We will use the same sample data as for the basic run above (<code>~/epilogos/data/pyData/male/</code> and <code>~/epilogos/data/pyData/female/</code>).</p>

<p>To compute epilogos for this sample data run one of the following commands (depending on the desired saliency metric) within the <code>~/epilogos/</code> directory (replacing <code>OUTPUTDIR</code> with the output directory of your choice).</p>

```bash
Saliency 1: $ epilogos -m paired -l -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR

Saliency 2: $ epilogos -m paired -l -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 2

Saliency 3: $ epilogos -m paired -l -a ./data/pyData/male/ -b ./data/pydata/male/ -n ./data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv -o OUTPUTDIR -s 3
```

<p>Upon completion of the run, you should see the files <code>pairwiseDelta_male_female_s$_epilogos_matrix_chr1.txt.gz</code>, <code>pairwiseMetrics_male_female_s$.txt.gz</code>, <code>greatestHits_male_female_s$.bed</code>, and <code>exp_freq_male_female_s$.npy</code> as well as the directory <code>manhattanPlots_male_female_s$</code> (where $ is replaced with the saliency metric you chose) in <code>OUTPUTDIR</code>. For further explanation of the contents of these outputs see <a href="output-directory-pairwise">Output Directory [-o, --output-directory]</a></p>

<p>To further customize your run of epilogos see the <a href="command-line-options-pairwise">Command Line Options</a> of the <code>README</code></p>

</details>

<details><summary><b> Running Epilogos with your own data</b></summary>
<p></p>

<p>In order to run Epilogos on your own data, you will need to do two things.</p>

<p>First, you will need to modify your data such that Epilogos can understand it. In order to assist with this, we have provided a bash script which takes ChromHMM files and generates Epilogos input files. This can be found at <code>~/epilogos/scripts/preprocess_data_ChromHMM.sh</code></p>

<p>Second, you will need to create a state info file. This is a tab separated file which tells epilogos various information about each of the states in the state model. We have provided some files already for common models in the <code>~/epilogos/data/state_metadata/</code> directory. For more information on the structure of these files see <code>~/epilogos/data/state_metadata/README.txt</code> or <a href="state-info">State Info [-n, --state-info]</a></p>

<p>Once you have completed these two things, you can run epilogos with the following command:</p>

```bash
$ epilogos -m paired -l -a PATH_TO_FIRST_INPUT_DIR -b PATH_TO_SECOND_INPUT_DIR -n PATH_TO_STATE_INFO_TSV -o PATH_TO_OUTPUT_DIR
```

</details>

<a name="command-line-options-pairwise"></a>

## Command Line Options

<a name="mode-pairwise"></a>
<details><summary><b> Mode [-m, --mode]</b></summary>
<p></p>
<p>Epilogos supports a single group and a paired group mode. The single group mode finds interesting regions compared to a background of itself. Whereas the paired group mode finds regions which differ significantly between the two groups.</p>

<p>
The argument to this flag either <code>single</code> or <code>paired</code> as the mode of operation, with <code>single</code> being the default.
</p>
</details>

<a name="command-line-pairwise"></a>
<details><summary><b> Command Line [-l, --command-line]</b></summary>
<p></p>
<p>By default, Epilogos assumes use of a SLURM cluster. However, if you would like to run Epilogos directly in the command line enable this flag</p>
</details>

<a name="directory-one"></a>
<details><summary><b> Input Directory One [-a, --directory-one]</b></summary>
<p></p>
<p>Rather than just read in one input file, Epilogos reads the contents of an entire directory. This allows the computation to be chunked and parallelized. Note that the genome files in the directory <strong>MUST</strong> be split by chromosome. Input files should be formatted such that the first three columns are the chromosome, bin start, and bin end respectively with the rest of the columns containing state data.</p>

<p>In the paired group version of epilogos, the user must input two directories (one for each group), the argument to this flag is the path to the first directory which contains the files to be read in. Note that <strong>ALL</strong> files in this directory will be read in and errors may occur if other files are present.
</p>
</details>

<a name="directory-two"></a>
<details><summary><b> Input Directory Two [-b, --directory-two]</b></summary>
<p></p>
<p>Rather than just read in one input file, Epilogos reads the contents of an entire directory. This allows the computation to be chunked and parallelized. Note that the genome files in the directory <strong>MUST</strong> be split by chromosome. Input files should be formatted such that the first three columns are the chromosome, bin start, and bin end respectively with the rest of the columns containing state data.</p>

<p>In the paired group version of epilogos, the user must input two directories (one for each group), the argument to this flag is the path to the second directory which contains the files to be read in. Note that <strong>ALL</strong> files in this directory will be read in and errors may occur if other files are present.
</p>
</details>

<a name="output-directory-pairwise"></a>
<details><summary><b> Output Directory [-o, --output-directory]</b></summary>
<p></p>
<p>The output of paired group Epilogos will vary depending on the number of input files present in the input directories <a href="#directory-one">[-a, --directory-one]</a><a href="#directory-two">[-b, --directory-two]</a>. All score difference files will be gzipped txt files and of the format <code>pairwiseDelta_{}_()_s$_[].txt.gz</code> where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively, <code>$</code> is replaced with the saliency level, and <code>[]</code> is replaced with the name of the corresponding input file (extensions removed).</p>

<p>The output directory will contain one <code>pairwiseMetrics_{}_()_s$.txt.gz</code> file where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level. Columns 1-3 contain the locations, column 4 contains the state with the largest difference between the scores, column 5 contains the squared euclidean distance between the scores, and column 6 contains the p-value of this distance.</p>

<p>The output directory will contain one <code>greatestHits_{}_()_s$.bed</code> file where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level. This file contains the state with top 1000 highest distance regions (adjacent regions are merged). Columns 1-3 contain the locations, column 4 contains name of the largest difference states, column 5 contains the squared euclidean distance between the scores, and column 6 contains the direction of this distance (determined by whether group 1 or 2 had higher signal).</p>

<p>The output directory will contain one <code>exp_freq_{}_()_s$.npy</code> file where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level. The file contains the expected frequencies for each of the states.</p>

<p>The output directory will contain one <code>manhattanPlots_{}_()_s$</code> directory where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level. This directory will contain all the manhattan plots generated by pairwise epilogos. These plots show the signed squared euclidean distances between groups 1 and 2 as well as the p-values of these distances. There is one genome-wide plot generated and another plot generate for each chromosome.</p>

<p>Depending on the <a href="#diagnostic-figures">[-d, --diagnostic-figures] flag</a> the output directory may contain one <code>diagnosticPlots_{}_()_s$</code> directory where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level. This directory will contain diagnostic figures of the gennorm fit on the null data and comparisons between the null and real data.</p>

<p>The argument to this flag is the path to the directory to which you would like to output. Note that this <strong>CANNOT</strong> be the same as the input directory.</p>
</details>

<a name="state-info-pairwise"></a>
<details><summary><b> State Info [-n, --state-info]</b></summary>
<p></p>
<p>The argument to this flag is a tab separated file specifying information about the state model being used. Example below (for more detail see <code>epilogos/data/state_metadata/README.md</code> or <code>epilogos/data/state_metadata/human/Adsera_et_al_833_sample/hg19/18/metadata.tsv</code></p>
| zero_index | one_index | short_name | long_name | hex | rgba | color |
|------------|-----------|------------|-----------|-----|------|-------|
| 0 | 1 | TssA | Active TSS | #ff0000 | rgba(255,0,0,1) | Red |
<p>
Note that tsv must contain a header row with the exact names above and that values within the table should follow the same formatting as above.
</p>
</details>


<a name="saliency-pairwise"></a>
<details><summary><b> Saliency Level [-s, --saliency]</b></summary>
<p></p>
<p>Epilogos implements information-theoretic metrics to quantify saliency levels of datasets. The <code>-l</code> flag to the coordination script allows one to choose one of three possible metrics:</p>

```
1. Metric S1, implementing a standard Kullback-Leibler relative entropy

2. Metric S2, implementing a version of S1 that additionally models label co-occurrence patterns

3. Metric S3, implementing a version of S2 that additionally models between-biosample similarities
```

<p>
Note that each increase in saliency level involves much more computation and thus each increase requires more time and computational power.
</p>

<p>
The arguement to this flag must be an integer <code>1, 2, or 3</code>. Note that Epilogos defaults to a saliency of 1.
</p>
</details>

<a name="number-of-cores-pairwise"></a>
<details><summary><b> Number of Cores [-c, --num-cores]</b></summary>
<p></p>
<p>Epilogos will always try and parallelize where it can. Computation done on each input file is parallelized using python's <a href="https://docs.python.org/3/library/multiprocessing.html">multiprocessing library</a>.</p>

<p>
The argument to this flag is an integer number of cores you would like to utilize to perform this multiprocessing. Note that Epilogos defaults to using all available cores (equivalent to <code>-c 0</code>).</p>
</details>

<a name="exit-pairwise"></a>
<details><summary><b> Exit [-x, --exit]</b></summary>
<p></p>
<p>By default <code>src/computeEpilogosSlurm.py</code> only exits after it has completed all slurm jobs and prints progress updates to the console. If you would like the program to instead exit when all jobs are submitted (allowing use of the terminal while the jobs are running), enable this flag.</p>
</details>

<a name="diagnostic-figures"></a>
<details><summary><b> Diagnostic Figures [-d, --diagnostic-figures]</b></summary>
<p></p>
<p>If this flag is enabled, Pairwise Epilogos will output diagnostic figures of the gennorm fit on the null data and comparisons between the null and real data. These can be found in a sub-directory of the output directory named <code>manhattanPlots_{}_()_s$</code> directory where <code>{}</code> and <code>()</code> are replaced with the names of input directories one and two respectively and <code>$</code> is replaced with the saliency level.</p>
</details>

<a name="num-trials"></a>
<details><summary><b> Number of Trials [-t, --num-trials]</b></summary>
<p></p>
<p>In order to save time when fitting in paired group Epilogos, samples of the null data are fit rather than the whole null data and then the median fit is used.</p>

<p>The argument to this flag is the number of times these samples are fit. Epilogos defaults to 101</P>
</details>

<a name="sampling-size"></a>
<details><summary><b> Sampling Size [-z, --sampling-size]</b></summary>
<p></p>
<p>In order to save time when fitting in paired group Epilogos, samples of the null data are fit rather than the whole null data and then the median fit is used.</p>

<p>The argument to this flag is the size of the samples that are fit. Epilogos defaults to 100000</P>
</details>