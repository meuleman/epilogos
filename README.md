# Epilogos

1. [About](#about)
2. [Prerequisites](#prerequisites)
3. [Running epilogos](#running-epilogos)
    * [Modes of operation](#modes-operation)
    * [Output](#output)
    * [Minimal example](#minimal-example)
4. [Visualizing results](#visualizing-results)
5. [Support](#support)

## About

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, 
with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos, using a combination of C++ and bash.
We provide a proof-of-principle dataset based on chromatin state calls from the Roadmap Epigenomics project.

## Prerequisites

To compute epilogos, we recommend installing the external programs `bgzip` and `starch`, for efficient storage and manipulation of output files.
The `bgzip` tool is part of the [htslib](https://github.com/samtools/htslib) kit, 
and `starch` is part of [BEDOPS](https://github.com/bedops/bedops).
Both binaries must be accessible through your `PATH` environment variable.

## Running epilogos

A single bash script, `scripts/computeEpilogos.sh`, coordinates all the processing.
Run the script with no arguments to see documentation on the list of arguments it requires.

The script calls three executables that need to be compiled, by running `make` from this directory.
This will compile the three programs and output the result in `bin`.
To allow the script to find this `bin` directory, add it to your PATH environment variable, e.g. using
```bash
$ export PATH=${PWD}/bin:${PATH}
```

By default, the script assumes the availability of a compute cluster managed by [SLURM](https://slurm.schedmd.com/).
See [below](#minimal-example) for a more elementary implementation of the script for use on a single machine, using data from a single chromosome only.

The first argument to `computeEpilogos.sh` must be the name of the cluster/queue.

The second argument is a text file with paths to *uncompressed* input data files (one per chromosome).
The first three columns of each input file must specify genomic coordinates (`seqname`, `start`, `end`),
and the remaining columns contain labels (e.g. chromatin state calls), representing (chromatin state) annotations -- one column per biosample.

The third argument, `numStates`, specifies the number of distinct labels (chromatin states) provided in the input data,
to prevent having to scan through the full input dataset.

### Modes of operation

Epilogos implements information-theoretic metrics to quantify saliency levels of datasets.
The fourth argument to the coordination script allows one to choose one of three possible metrics:
1. Metric S1, implementing a standard Kullback-Leibler relative entropy 
2. Metric S2, implementing a version of S1 that additionally models label co-occurrence patterns
3. Metric S3, implementing a version of S2 that additionally models between-biosample similarities

The fifth argument, `groupSpec`, specifies which biosamples, i.e. columns in the input data, are to be used in the computation.
This argument can be specified as comma- and/or hyphen-delimited list (e.g. "1,3,6-10,12,13").
For instance, "1,2,3" corresponds to the first 3 samples, which are in columns 4,5,6 of the input data.

Importantly, Epilogos implements the possibility for comparing annotations between two groups of biosamples.
This pairwise comparison mode can be enabled with a sixth argument, denoting a separate set of columns/biosamples to compare the first set to.
The sixth argument has the same format as the fifth argument.

### Output

Regardless of the mode of operation, upon successful completion of execution, 
three output files will be produced: `observations.starch`, `scores.txt.gz`, and `exemplarRegions.txt`.
In case of errors during the exuction of the code, error messages will be collected in separate log files.

`observations.starch` is a compressed file (uncompress it using `unstarch` from the `bedops` tool suite) containing one line of summary statistics for each input site.
Its exact format depends on the specified mode of operation (see below).
Regardless, the value that is likely of most interest to typical users are the overal (S1, S2 or S3) metric scores for each genomic coordinate.

`scores.txt.gz` is a compressed file (compressed using `bgzip`) containing all sites and the signed contributions from each state to the metric at each site.
Columns 1-3 hold the site's coordinates.
Column 4 holds the contribution from state 1, column 5 holds the contribution from state 2, etcetera; if there are *n* states, the file will contain *n* + 3 columns.
This file can be used for downstream analyses or visualization.

`exemplarRegions.txt` is a text file containing single sites or contiguous stretches of sites where the highest saliency scores were observed, sorted in descending order by score.
Its format is identical to the corresponding `observations.starch` file.

#### Single group of biosamples

For saliency metric S1 (standard relative entropy), `observations.starch` will contain 7 columns, e.g.:
``chr1	2468800	     2469000	 13	    2.56266  1	3.79557``
The specifications of the columns are as follows:
1. Chromosome
2. Start coordinate
3. End coordinate
4. Label/state with the largest contribution to the S1 metric score (i.e, largest relative information content)
5. Magnitude of this state's contribution.
6. Constant value of `1` (used in other modes of operation, see below).
7. Overall S1 metric score, i.e. the relative entropy of the specified genomic region.

For saliency metrics S2 and S3, which both deal with co-occurrence patterns of pairs of labels/states, `observations.starch` contains 10 columns, e.g.:
``chr3		     125932600	    125932800 13   3.65048  1	    (2,11)	 2.27667	1    20.2829``
The first 6 columns are specified identical to the S1 metric, and the remainder of the columns are specified as follows:
<ol start="7">
<li>The label/state pair with the largest contribution to the S2 or S3 metric</li>
<li>Magnitude of the contribution of this pair</li>
<li>Constant value of `1` (used in other modes of operation, see below).</li>
<li>Overall S2 or S3 metric score</li>
</ol>

#### Pairwise comparison between groups of biosamples

For saliency metric S1 (standard relative entropy), `observations.starch` will contain 8 columns, 
the first 7 columns are specified similarly to the single-group case, with the following exceptions and additional column:
<ol start="4">
<li>Label/state with the largest difference in information between the two groups of biosamples.</li>
<li>Absolute magnitude of this difference</li>
<li>Sign of the difference, i.e. wether the first (`1`) or second (`-1`) group of biosamples contributes more information.</li>
<li>Overall differential S1 metric score, i.e. the sum of all absolute per-state information differences between the two groups.</li>
<li>Empirical p-value of this differential S1 metric score.</li>
</ol>

P-values are obtained by creating an empirical null distribution of differential S1 metric values by randomly shuffling state and epigenome labels.
These are nominal p-values *not* adjusted for multiple hypothesis testing; an external procedure must be used to estimate false discovery rates (FDRs) from these p-values.

For saliency metrics S2 and S3, which both deal with co-occurrence patterns of pairs of labels/states, `observations.starch` contains 11 columns, 
specified analogous to what is described above for pairwise metric S1 and single-group S2 and S3, with an additional empirical p-value column.

<!--
the appearance of -1 in column 9 means that for state pair (2,11), the contribution from *B* increased DKL\*\* while *A* decreased it, with a next contribution of 2.27667 (column 8) from state pair (2,11).
In this case, it was estimated that p < 6.58697e-08 for observing a DKL\*\* score of 20.2829 or higher due to random chance alone.
-->

### Minimal example

A smaller version of the script, `computeEpilogos_minimal.sh`, has been provided for running epilogos without access to a SLURM cluster, on data from a single chromosome only.

The following sample input data are provided in the `data` directory:
* The file `chr1_127epigenomes_15observedStates.txt.gz` contains ChromHMM chromatin state calls for a 15-state chromatin state model, across 200bp genomic bins spanning human chromosome 1.
* The file `Blood_T-cellGroupSpec.txt` contains the column specifications for a group of blood and T-cell biosamples in the input data (33-34,37-45,47-48,61).
* The file `HSC_B-cellGroupSpec.txt` contains the column specifications for a group of hematopoietic stem cell (HSC) and B-cell samples in the input data (29-32,35-36,46,50-51).

To compute epilogos (using the S1 saliency metric) for blood and T-cell biosamples only, run the following command from the `data` subdirectory:
```bash
$ ../scripts/computeEpilogos_minimal.sh chr1_127epigenomes_15observedStates.txt.gz 1 15 OUTPUTDIR "33-34,37-45,47-48,61"
```
The resulting output files `observations.starch`, `scores.txt.gz`, and `exemplarRegions.txt` in directory `OUTPUTDIR` should match the corresponding files in `data/results_Blood_T-cell/KL`.

To compute differential epilogos (using the differential S1 saliency metric) between two groups of biosamples, run the following command, again from the epilogos `data` subdirectory:
```bash
$ ../scripts/computeEpilogos_minimal.sh chr1_127epigenomes_15observedStates.txt 1 15 OUTPUTDIR2 "29-32,35-36,46,50-51" "33-34,37-45,47-48,61"
```
The resulting output files `observations.starch`, `scores.txt.gz`, and `exemplarRegions.txt` in directory `OUTPUTDIR2` should match the corresponding files in `data/results_HSC_B-cell_vs_Blood_T-cell/DKL`.

## Visualizing results

We recommend using [HiGlass](https://higlass.io) to visualize the per-site per-state results written to the file `scores.txt.gz`.
Further instructions are forthcoming.

## Support

To get additional help or if you have questions about this software, open an [issue ticket](https://github.com/Altius/epilogos/issues).


