This is the repository for the software that computes Epilogos components
from observations of states at sites across the genome.

More info will be added to this README file soon.

In the meantime, here's a very brief tutorial:

* To compute epilogos, you must have the external programs `bgzip`
and `starch` (the latter is part of `bedops`) in your PATH
environment variable.
* A single script, `computeEpilogos.sh`, does all the processing.
(But see below for important differences between `computeEpilogos.sh`
and `computeEpilogos_singleChromosomeSingleProcessor.sh`.)
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
results.  All processing will be done on a computer cluster, managed by SLURM.
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


---

## Visualizing a qcat file with the WashU Epigenome Browser

This section describes how to visualize the qcat result file, using the WashU Epigenome Browser. 

In broad terms, this can be done by way of:

1. Setting up a public web server
2. Indexing the qcat result with tabix
3. Creating a JSON-formatted *datahub* file specific to the type of epilogos (single or paired group) and your web server
4. Hosting the tabix and JSON assets on your web server
5. Loading the WashU browser with your custom datahub

### Setting up a web server

A public-facing web server is a basic requirement for visualizing your qcat data with the WashU browser. Your web server makes your qcat-formatted epilogos data available, as well as related metadata about the epilogos track, chromatin states, and other annotation data.

This web server must be able to serve files via port tcp/80 (http). You will need the public-facing static IP address assigned to this web server; in other words, you must be able to access this server through firewalls and outside your local network.

Your institution or research facility may already offer a basic, public web service, and contacting your local IT help may get you information on this option to get you started, such as where you would store your files and, importantly, the web address of the server.

Depending on what is available through your local setup, there are also (fee-based) web hosting services, such as [Dreamhost](https://www.dreamhost.com/) or [Amazon Lightsail](https://lightsail.aws.amazon.com). Management tools for these services are web-based and make configuration easy. 

In the case of Lightsail, for example, you could set up a virtual machine instance running *nginx* (a web server), and then assign a static IP address to that instance. Epilogos-related files are copied to this instance in or under `/opt/bitnami/nginx/html` (or similar). The contents of this directory are available from a web browser, using the static address that Amazon has assigned the instance.

For the purposes of this section, we will assume that you have a public-facing web server that is available at `http://192.168.0.1`, and that you are able to copy files to the required directory on this server, so that they are available through this address.

### Preparing tabix-indexed qcat files

The next step is to preprocess the qcat output from epilogos, so that smaller regions of interest can be quickly retrieved via the web. This is a two-step process: 

1. Compressing with `bgzip`
2. Indexing with `tabix`

Both `bgzip` and `tabix` are programs that are available through the larger `htslib` software package. This can be installed through package managers like [Homebrew](https://brew.sh/) or [Bioconda](https://bioconda.github.io/), or you can download and compile source from [Github](https://github.com/samtools/htslib/blob/develop/INSTALL). You should be able to run `bgzip --help` and `tabix --help` or similar to verify their availability.

#### Compressing

The following creates a `bgzip`-compressed version of the `qcat.bed.gz` output file, which is called `qcat.bed.bgz`:

```bash
$ gunzip -c qcat.bed.gz | bgzip -c > qcat.bed.bgz
```

#### Indexing

The following makes an index file called `qcat.bed.bgz.tbi`

```bash
$ tabix -p bed qcat.bed.bgz
```

Copy both `qcat.bed.bgz` and `qcat.bed.bgz.tbi` to the same directory associated with your web server. You can put these files into a subdirectory, but you must keep both files together in the same directory, at the same directory level.

### Creating a datahub

#### Single group

```json

```

#### Paired group

```json

```

Copy this JSON file to the directory associated with your web server. You can put this file into a subdirectory; the name of a subdirectory simply modifies the web address for the JSON file.

### Loading the WashU browser

Now that you have generated the tabix-indexed qcat and datahub JSON assets and added them to your web server, you are ready to generate a web address that will load the WashU browser and instruct it to load your data.

The basic format of this address is:

> http://epigenomegateway.wustl.edu/browser/?genome=**build_name**&datahub=**datahub_address**

The *build_name* is replaced with one of `hg19`, `hg38`, or `mm10`, depending on how you made your epilogos dataset.

The *datahub_address* is replaced with the web address that points to your datahub JSON file.

Here's a more concrete example, which specifies the `hg19` assembly and points to a hypothetical datahub file available at `http://192.168.0.1/mydatahub.json`:

> http://epigenomegateway.wustl.edu/browser/?genome=hg19&datahub=http://192.168.0.1/mydatahub.json

You can open this in a web browser, test it, modify it, and share it with others. It will persist as long as you have your web server up and running.

Other parameters may be added to this address, which customize the behavior and appearance of the WashU browser. A more complete listing of parameters is available from the WashU browser [wiki](http://wiki.wubrowse.org/URL_parameter).