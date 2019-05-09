# Epilogos

1. [About](#about)
2. [Prerequisites](#prerequisites)
3. [Running epilogos](#running-epilogos)
    * [Output](#output)
    * [Single-chromosome execution](#single-chromosome-execution)
4. [Visualizing results](#visualizing-results)
5. [Support](#support)

## About

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, 
with an emphasis on chromatin state maps generated with, e.g. ChromHMM or Segway.

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
See [below](#single-chromosome-execution) for a more elementary implementation of the script for use on a single machine, using data from a single chromosome only.
The first argument to `computeEpilogos.sh` must be the name of the cluster/queue.
The second argument is a text file with paths to *uncompressed* input data files (one per chromosome).

### Output

Three output files will be produced: `observations.starch`, `scores.txt.gz`, and `exemplarRegions.txt`.
In case of errors during the exuction of the code, error messages will be collected in separate log files.

`observations.starch` is a compressed file (uncompress it using `unstarch` from the `bedops` tool suite) containing one line of summary statistics for each input site.
Its exact contents will depend on whether you're measuring occurrences of states (metric KL or DKL, specified via `0`), state pairs (metric KL\* or DKL\*, specified via `1`), or state pairs while tracking the epigenome pairs in which they were observed (metric KL\*\* or DKL\*\*, specified via `2`).

For KL, `observations.starch` will contain 7 columns.
An example: `chr1	2468800	     2469000	 13	    2.56266  1	3.79557`.
Columns 1-3 hold the site's coordinates.
Column 4 holds the state that made the largest contribution to the overall score (i.e., to the KL or DKL value).
Column 5 holds the magnitude of this contribution.
When DKL is being computed to compare two groups of epigenomes, the contribution to DKL given in column 5 is the difference between the contributions from that state in the two groups, and column 6 will hold either 1 or -1, indicating whether the second group (-1) or the first group (1) contributed more to DKL; when KL is being computed on a single collection of epigenomes, all entries in column 6 will be 1.
Column 7 holds the overall score (KL or DKL value).
When DKL is being computed, an empirical p-value estimate for DKL (column 7) is given in column 8.
The p-values are obtained by creating an empirical null distribution of DKL values by randomly shuffling state and epigenome labels.
These p-values are *not* adjusted for the large number of measurements made; an external program of the user's choosing must be used to estimate false discovery rates (FDRs) from these p-values.

For KL\* and KL\*\*, `observations.starch` will contain 10 columns; for DKL\* and DKL\*\*, p-values will be appended in an 11th column.
An example using DKL\*\*, comparing two groups of epigenomes labeled *A* and *B* : `chr3		     125932600	    125932800 13   3.65048  1	    (2,11)	 2.27667	-1    20.2829	 6.58696e-08`.
Columns 1-6 have the same meanings as for KL and DKL:  in this example, at chr3:125932600-125932800, state 13 (column 4) made the largest contribution, and because column 6 is +1, we know that the contribution to DKL\*\* from state 13 (3.65048, column 5), which is the difference between the contributions from state 13 from *A* and *B*, was larger in *A*.
At the same time, the co-occurrence of states 2 and 11 (column 7) in one or more pairs of epigenomes had the largest impact on DKL\*\*; the appearance of -1 in column 9 means that for state pair (2,11), the contribution from *B* increased DKL\*\* while *A* decreased it, with a next contribution of 2.27667 (column 8) from state pair (2,11).
In this case, it was estimated that p < 6.58697e-08 for observing a DKL\*\* score of 20.2829 or higher due to random chance alone.

`exemplarRegions.txt` is a text file containing single sites or contiguous stretches of sites where the highest scores were observed, sorted in descending order by score.
Its format is identical to the corresponding `observations.starch` file.

`scores.txt.gz` is a compressed file (compressed using `bgzip`) containing all sites and the signed contributions from each state to the metric at each site.
Columns 1-3 hold the site's coordinates.
Column 4 holds the contribution from state 1, column 5 holds the contribution from state 2, etc.; if there are *N* states, the file will contain *N* + 3 columns.
This file can be used to visualize the data in a browser; see the section below regarding data visualization.

### Single-chromosome execution

A slightly smaller version of the script, `computeEpilogos_singleChromosomeSingleProcessor.sh`, has been provided for running data from a single chromosome on a single processor, as opposed to `computeEpilogos.sh`, which assumes multiple uncompressed input files (one per chromosome) and access to a compute cluster managed by a [SLURM](https://slurm.schedmd.com/) job scheduler.

Sample input data has been provided.
If you set up epilogos correctly, you should be able to use the input data, write results into new directories of your choosing, and then ensure that those results match the results provided alongside the input data.
This input file (`chr1_127epigenomes_15observedStates.txt.gz`) can be provided as-is as the first argument to `computeEpilogos_singleChromosomeSingleProcessor.sh`.
(To use it with `computeEpilogos.sh`, you will need to decompress it using `gunzip` and supply a text file containing the path to the decompressed version of the file.)

The file `Blood_T-cellGroupSpec.txt` contains the column specifications for a group of blood and T-cell samples in the input data (33-34,37-45,47-48,61).
To compute KL from this subset of the input data, cd to the epilogos `data` subdirectory, then run the following command:

```bash
$ ../scripts/computeEpilogos_singleChromosomeSingleProcessor.sh chr1_127epigenomes_15observedStates.txt.gz 0 15 yourOutputDirectory1/KL "33-34,37-45,47-48,61"
```

The files `yourOutputDirectory1/KL/observations.starch`, `yourOutputDirectory1/KL/scores.txt.gz`, and `yourOutputDirectory1/KL/exemplarRegions.txt` that your run produces should match the corresponding files in `data/results_Blood_T-cell/KL`.

The file `HSC_B-cellGroupSpec.txt` contains the column specifications for a group of stem-cell and B-cell samples in the input data (29-32,35-36,46,50-51).
To compute DKL for a comparison of these two subsets of the input data, run the following command, again from the epilogos `data` subdirectory:

```bash
$ ../scripts/computeEpilogos_singleChromosomeSingleProcessor.sh chr1_127epigenomes_15observedStates.txt 0 15 yourOutputDirectory2/DKL "29-32,35-36,46,50-51" "33-34,37-45,47-48,61"
```

The files `yourOutputDirectory2/DKL/observations.starch`, `yourOutputDirectory2/DKL/scores.txt.gz`, and `yourOutputDirectory2/DKL/exemplarRegions.txt` that your run produces should match the corresponding files in `data/results_HSC_B-cell_vs_Blood_T-cell/DKL`.

## Visualizing results

We strongly recommend using [HiGlass](https://higlass.io) to visualize the per-site per-state results written to the file `scores.txt.gz`.
Further instructions are forthcoming.

## Support

To get additional help or if you have questions about this software, open an [issue ticket](https://github.com/Altius/epilogos/issues).
