This is the repository for the software that computes Epilogos components
from observations of states at sites across the genome.

More info will be added to this README file soon.

In the meantime, here's a very brief tutorial:

* To compute epilogos, you must have the external programs `bgzip`
and `starch` (the latter is part of `bedops`) in your PATH
environment variable.
* A single script, `computeEpilogos.sh`, does all the processing.
Before you run the script, you need to make the executable programs
that the script calls.  This is done by running `make` from
this directory.  This will make 3 programs, and place them into a
parallel directory named `bin`.
* After you have made these programs, and before you run the script,
you need to do one more thing:  Enable the script to find the `bin`
directory so it can access and run those programs. To do this,
add the full path to the `bin` directory to your PATH environment
variable.
* Run the script `computeEpilogos.sh` with no arguments to see
the list of arguments it requires and descriptions of them.
* Then you can run the script and supply your arguments and get your
results.  All processing will be done on the cluster, managed by SLURM.
For the time being, the cluster/queue name is hardcoded in the script.
* The key output files are `observations.starch` and `qcat.bed.gz`.
Additional files will be created during execution and deleted when
no longer needed. Various small files that will persist will contain
error messages in the event of an error, and otherwise be empty.


A slightly smaller version of the script, `computeEpilogos_singleChromosomeSingleProcessor.sh`,
has been provided for running data from a single chromosome on a single processor,
as opposed to `computeEpilogos.sh`, which assumes multiple input files (one per chromosome)
and access to a computer cluster managed by SLURM software.


Sample input data has been provided. If you set up epilogos correctly,
you should be able to use the input data, write results into new
directories of your choosing, and then ensure that those results
match the results provided alongside the input data. NOTE: You must
use gunzip to decompress the input data file before using it.


The file `Blood_T-cellGroupSpec.txt` contains the column specifications
for a group of blood and T-cell samples in the input data. To compute KL
from this subset of the input data, cd to the epilogos `data` subdirectory,
decompress the input data as mentioned above, then run the following command:


../scripts/computeEpilogos_singleChromosomeSingleProcessor.sh chr1_127epigenomes_15observedStates.txt 0 15 yourOutputDirectory1/KL "33-34,37-45,47-48,61"


The files `yourOutputDirectory1/KL/observations.starch` and `yourOutputDirectory1/KL/qcat.bed.gz`
that your run produces should match the corresponding files in `data/results_Blood_T-cell/KL`.


The file `HSC_B-cellGroupSpec.txt` contains the column specification
for a group of stem-cell and B-cell samples in the input data. To compute DKL
for a comparison of these two subsets of the input data, run the following command,
again from the epilogos `data` subdirectory:



../scripts/computeEpilogos_singleChromosomeSingleProcessor.sh chr1_127epigenomes_15observedStates.txt 0 15 yourOutputDirectory2/DKL "29-32,35-36,46,50-51" "33-34,37-45,47-48,61"


The files `yourOutputDirectory2/DKL/observations.starch` and `yourOutputDirectory2/DKL/qcat.bed.gz`
that your run produces should match the corresponding files in `data/results_HSC_B-cell_vs_Blood_T-cell/DKL`.


