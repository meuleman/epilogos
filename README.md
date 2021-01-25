# Epilogos

1. [About](#about)
2. [Prerequisites](#prerequisites)
3. [Running Epilogos](#running-epilogos)
    * [Input Directory](#input-directory)
    * [Output Directory](#output-directory)
    * [State Model](#state-model)
    * [Saliency Level](#saliency-level)
    * [Mode of Operation](#mode-of-operation)
    * [Background Directory](#background-directory)
    * [Number of Cores](#number-of-cores)
    * [Exit When Complete](#exit-when-complete)
4. [Minimal example](#minimal-example)

## About

Epilogos is an approach for analyzing, visualizing, and navigating multi-biosample functional genomic annotations, with an emphasis on chromatin state maps generated with e.g. ChromHMM or Segway.

The software provided in this repository implements the methods underlying Epilogos using only python. We provide a proof-of-principle dataset based on chromatin state calls from the BOIX dataset.

## Prerequisites

To compute epilogos, you will need to have the following python libraries: [click](https://click.palletsprojects.com/en/7.x/), [numpy](https://numpy.org/), and [pandas](https://pandas.pydata.org/). These can be installed with the following command.
```bash
$ pip install click, numpy, pandas
```

## Running Epilogos:

A single python script, `src/computeEpilogosSlurm.py`, controls all of the processing. To be presented with minimal documentation of arguments needed to run epilogos, simply run python `src/computeEpilogosSlurm.py --help` (More in-depth explanation is given [below](#input-directory))

The script, `src/computeEpilogosSlurm.py`, requires access to a computational cluster managed by [SLURM](https://slurm.schedmd.com/). A minimal version of epilogos, `src/minimalEpilogos.py`, has been created for those without access to a SLURM cluster. It functions identically to `src/computeEpilogosSlurm.py` but runs everything within one terminal command.

### Input Directory (-f, --file-directory)

Rather than just read in one input file, Epilogos reads the contents of an entire directory. This allows the computation to be chunked and parallelized. Additionally, it allows users to separate data as makes sense to them (e.g. split up the genome by chromosome)

The argument to this flag is the path to the directory which contains the files to be read in. Note that ALL files in this directory will be read in.

### Output Directory (-o, --output-directory)

The output of Epilogos will vary depending on the number of input files present in the input directory (-f, --file-directory). All scores files will be gzipped txt files and of the format `scores_{}_[].txt.gz` where {} is replaced with the input directory name and [] is replaced with the name of the corresponding input file (extensions removed).

The argument to this flag is the path to the directory to which you would like to output. Note that this may not be the same as the input directory.

### State Model (-s, --state-model)

The argument to this flag specifies the number of distinct labels (chromatin states) provided in the input data.

### Saliency Level (-l, --saliency-level)

Epilogos implements information-theoretic metrics to quantify saliency levels of datasets. The -l flag to the coordination script allows one to choose one of three possible metrics:
1. Metric S1, implementing a standard Kullback-Leibler relative entropy
2. Metric S2, implementing a version of S1 that additionally models label co-occurrence patterns
3. Metric S3, implementing a version of S2 that additionally models between-biosample similarities

Note that each increase in saliency level involves much more computation and thus each increase requires more time and computational power.

The arguement to this flag must be an integer 1, 2, or 3. Note that Epilogos defaults to a saliency of 1.

### Mode of Operation (-m, --mode-of-operation)

As epilogos has 2 different types of output files, we allow the user to designate which they would like to receive and thus minimize potentially repeated computation.

The argument to this flag must be one of three strings: `bg`, `s`, `both`. If you would like to calculate only the background frequencies of the chromatin statesm use `bg`. If you already have a file containing the background frequencies and would only like to calculate the per state scores, use `s`. If you would like to calculate both the background frequencies and the scores, use `both`. Note that Epilogos defaults to `both`.

### Background Directory (-b, --background-directory)

In the case that the user chooses `s` as the mode of operation, the argument to this flag is the directory in which the background frequency file resides. Note that the file must maintain the same name as it was given upon original output. The format for this name is `exp_freq_{}.npy` where {} is replace with the name of the input directory. Note that Epilogos defaults to the ouput directory.

In the case that the user chooses either `bg`" or `both` as the mode of operation, the argument to this flag is the directory to which the background frequencies should be written. This is in case you want the background frequency output directory to be different from the score output directory. Note that Epilogos defaults to the ouput directory.

### Number of Cores (-c, --num-cores)

Epilogos will always try and parallelize where it can. Computation done on each input file is parallelized using python's [multiprocessing](https://docs.python.org/3/library/multiprocessing.html) library.

The argument to this flag is an integer number of cores you would like to utilize to perform this multiprocessing. Note that Epilogos defaults to using all available cores (equivalent to `-c 0`).

### Exit Timing (-x, --exit-when-complete)

By default `src/computeEpilogosSlurm.py` exits after it has submitted all slurm jobs. This allows the user to continue use of their terminal while the jobs are running. If you would like the program to instead exit when all jobs are done, enable this flag.

## Minimal Example

Sample data has been provided under `~/epilogos/data/pyData/male/`. The file, `epilogos_matrix_chr1.txt.gz`, contains chromatin state calls for a 18-state chromatin model, across 200bp genomic bins spanning human chromosome 1. The data was pulled from the [BOIX dataset](https://docs.google.com/spreadsheets/d/103XbiwChp9sJhUXDJr9ztYEPL00_MqvJgYPG-KZ7WME/edit#gid=1813267486) and contains only those epigenomes which are tagged `Male` under the `Sex` column

To compute epilogos (using the S1 saliency metric) for this sample data run one of the following commands (depending on if you want to use SLURM or not) within the `~/epilogos/` directory.
```bash
$ python ./src/computeEpilogosSlurm.py -f ./data/pyData/male/ -s 18 -o OUTPUTDIR
```
```bash
$ python ./src/minimalEpilogos.py -f ./data/pyData/male/ -s 18 -o OUTPUTDIR
```
Replacing `OUTPUTDIR` with the output directory of your choice. 

Upon completion of the run, you should see the files `exp_freq_male.npy` and `scores_male_epilogos_matrix_chr1.txt.gz` in `OUTPUTDIR`

To customize your run of epilogos see the [Running Epilogos](#running-epilogos) of the `README`