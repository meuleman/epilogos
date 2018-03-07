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

This section describes how to visualize the qcat result file ("epilogos"), using the [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/). 

In broad terms, this can be done by way of:

1. Setting up a public web server
2. Indexing the qcat result with tabix
3. Creating a JSON-formatted *datahub* file
4. Copying the qcat file, tabix index, and JSON datahub assets to your web server
5. Loading the WashU browser with your custom datahub

### Setting up a web server

A public-facing web server is a basic requirement for visualizing your qcat data with the WashU browser. Your web server will make your qcat-formatted epilogos data available, as well as related metadata about the epilogos track, chromatin states, and optional annotation data.

This web server must be able to serve files via port tcp/80 (http). You will need the public-facing static IP address assigned to this web server; in other words, you must be able to access this server through firewalls and outside your local network.

Your institution or research facility may already offer a public-facing web service, and contacting your local IT help may get you information on this option to get you started, such as where you would store your files and, importantly, the web address of the server.

Depending on what is available locally, there are also (fee-based) web hosting services, such as [Dreamhost](https://www.dreamhost.com/) or [Amazon Lightsail](https://lightsail.aws.amazon.com), for instance. Management tools for these services are web-based and make configuration easy.

For the purposes of this section, we will assume that you have a working, public-facing web server that is available at `http://192.168.0.1` (your actual IP address will be different), and that you are able to copy or upload files to the required directory on this server, so that these files are available through this address.

### Preparing tabix-indexed qcat files

The next step is to preprocess the qcat output from epilogos, so that smaller regions of interest can be quickly retrieved via the web. This is a two-step process: 

1. Compressing with `bgzip`
2. Indexing with `tabix`

Both `bgzip` and `tabix` are programs that are available through the larger `htslib` software package. This can be installed through package managers like [Homebrew](https://brew.sh/) or [Bioconda](https://bioconda.github.io/), or you can download and compile source from [Github](https://github.com/samtools/htslib/blob/develop/INSTALL). You should be able to run `bgzip --help` and `tabix --help` or similar to verify their availability.

#### Compressing

The following command creates a `bgzip`-compressed version of the `qcat.bed.gz` output file, which is called `qcat.bed.bgz`:

```bash
$ gunzip -c qcat.bed.gz | bgzip -c > qcat.bed.bgz
```

#### Indexing

The following command makes an index file called `qcat.bed.bgz.tbi`

```bash
$ tabix -p bed qcat.bed.bgz
```

Copy or upload both `qcat.bed.bgz` and `qcat.bed.bgz.tbi` to the same directory associated with your web server. You can put these files into a subdirectory, but you must keep both files together in the same directory, at the same directory level.

### Creating a datahub

A **datahub** is a JSON-formatted text file that tells the WashU browser where your qcat files are located, and how they should be rendered. The structure of the datahub is the same for single- and paired-group qcat files, which we describe in the template below.

#### Datahub template

```json
[
  {
    "set": {
      "1": {
        "11": [
          "Flanking Bivalent TSS/Enh",
          "#e9967a"
        ],
        "10": [
          "Bivalent/Poised TSS",
          "#cd5c5c"
        ],
        "13": [
          "Repressed PolyComb",
          "#808080"
        ],
        "12": [
          "Bivalent Enhancer",
          "#bdb76b"
        ],
        "15": [
          "Quiescent/Low",
          "#ffffff"
        ],
        "14": [
          "Weak Repressed PolyComb",
          "#c0c0c0"
        ],
        "1": [
          "Active TSS",
          "#ff0000"
        ],
        "3": [
          "Transcr at gene 5' and 3'",
          "#32cd32"
        ],
        "2": [
          "Flanking Active TSS",
          "#ff4500"
        ],
        "5": [
          "Weak transcription",
          "#006400"
        ],
        "4": [
          "Strong transcription",
          "#008000"
        ],
        "7": [
          "Enhancers",
          "#ffff00"
        ],
        "6": [
          "Genic enhancers",
          "#c2e105"
        ],
        "9": [
          "Heterochromatin",
          "#8a91d0"
        ],
        "8": [
          "ZNF genes & repeats",
          "#66cdaa"
        ]
      },
      "2": {
        "11": [
          "Weak Enhancer",
          "#FFFF00"
        ],
        "10": [
          "Active Enhancer 2",
          "#FFC34D"
        ],
        "13": [
          "Heterochromatin",
          "#8A91D0"
        ],
        "12": [
          "ZNF genes & repeats",
          "#66CDAA"
        ],
        "15": [
          "Bivalent Enhancer",
          "#BDB76B"
        ],
        "14": [
          "Bivalent/Poised TSS",
          "#CD5C5C"
        ],
        "17": [
          "Weak Repressed PolyComb",
          "#C0C0C0"
        ],
        "16": [
          "Repressed PolyComb",
          "#808080"
        ],
        "18": [
          "Quiescent/Low",
          "#FFFFFF"
        ],
        "1": [
          "Active TSS",
          "#FF0000"
        ],
        "3": [
          "Flanking TSS Upstream",
          "#FF4500"
        ],
        "2": [
          "Flanking TSS",
          "#FF4500"
        ],
        "5": [
          "Strong transcription",
          "#008000"
        ],
        "4": [
          "Flanking TSS Downstream",
          "#FF4500"
        ],
        "7": [
          "Genic enhancer 1",
          "#C2E105"
        ],
        "6": [
          "Weak transcription",
          "#006400"
        ],
        "9": [
          "Active Enhancer 1",
          "#FFC34D"
        ],
        "8": [
          "Genic enhancer 2",
          "#C2E105"
        ]
      },
      "3": {
        "24": [
          "Repressed PolyComb",
          "#808080"
        ],
        "25": [
          "Quiescent/Low",
          "#ffffff"
        ],
        "20": [
          "ZNF genes & repeats",
          "#66cdaa"
        ],
        "21": [
          "Heterochromatin",
          "#8a91d0"
        ],
        "22": [
          "Poised Promoter",
          "#e6b8b7"
        ],
        "23": [
          "Bivalent Promoter",
          "#7030a0"
        ],
        "1": [
          "Active TSS",
          "#ff0000"
        ],
        "3": [
          "Promoter Downstream TSS with DNase",
          "#ff4500"
        ],
        "2": [
          "Promoter Upstream TSS",
          "#ff4500"
        ],
        "5": [
          "Transcription 5'",
          "#008000"
        ],
        "4": [
          "Promoter Downstream TSS",
          "#ff4500"
        ],
        "7": [
          "Transcription 3'",
          "#008000"
        ],
        "6": [
          "Transcription",
          "#008000"
        ],
        "9": [
          "Transcription Regulatory",
          "#c2e105"
        ],
        "8": [
          "Weak transcription",
          "#009600"
        ],
        "11": [
          "Transcription 3' Enhancer",
          "#c2e105"
        ],
        "10": [
          "Transcription 5' Enhancer",
          "#c2e105"
        ],
        "13": [
          "Active Enhancer 1",
          "#ffc34d"
        ],
        "12": [
          "Transcription Weak Enhancer",
          "#c2e105"
        ],
        "15": [
          "Active Enhancer Flank",
          "#ffc34d"
        ],
        "14": [
          "Active Enhancer 2",
          "#ffc34d"
        ],
        "17": [
          "Weak Enhancer 2",
          "#ffff00"
        ],
        "16": [
          "Weak Enhancer 1",
          "#ffff00"
        ],
        "19": [
          "DNase only",
          "#ffff66"
        ],
        "18": [
          "Enhancer Acetylation Only",
          "#ffff00"
        ]
      },
      "4": {
        "11": [
          "Transcription - Permissive",
          "#deecf7"
        ],
        "10": [
          "Transcription - Strong",
          "#0454a3"
        ],
        "13": [
          "Heterochromatin - Polycomb",
          "#f48c8f"
        ],
        "12": [
          "Transcription - Initiation",
          "#4290cf"
        ],
        "15": [
          "No signal",
          "#ffffff"
        ],
        "14": [
          "Heterochromatin - H3K9me3",
          "#fde2e5"
        ],
        "1": [
          "Promoter - Active",
          "#0e6f37"
        ],
        "3": [
          "Promoter - Bivalent",
          "#cdcdcd"
        ],
        "2": [
          "Promoter - Weak/Inactive",
          "#c7e4c0"
        ],
        "5": [
          "Enhancer - Strong, TSS-distal",
          "#f3eb1a"
        ],
        "4": [
          "Promoter - Flanking",
          "#41ac5e"
        ],
        "7": [
          "Enhancer - Weak, TSS-distal",
          "#faf8c8"
        ],
        "6": [
          "Enhancer - Strong, TSS-proximal",
          "#f3eb1a"
        ],
        "9": [
          "Enhancer - Poised, TSS-proximal",
          "#808080"
        ],
        "8": [
          "Enhancer - Poised, TSS-distal",
          "#808080"
        ]
      }
    },
    "type": "category_set"
  },
  {
    "category_set_index": 1,
    "name": "My sample of interest",
    "url": "http://192.168.0.1/qcat.bed.bgz",
    "height": 225,
    "mode": "show",
    "backgroundcolor": "#000000",
    "type": "quantitativeCategorySeries"
  }
]
```

You can copy this text to your text editor of choice and make modifications, described below.

For this example JSON, there are two main objects. The first object `set` contains chromatin states, names, and colors. You must include this object, unedited, in order to color the epilogos data correctly.

You would modify three variables in the second and last object; specifically: `category_set_index`, `name`, and `url`. Optionally, you can adjust the track height by changing the `height` key. Include the remaining three key-value pairs for `mode`, `backgroundcolor`, and `type`, and leave their values unchanged.

The value of `category_set_index` should be set to one of `1`, `2`, `3`, or `4`. The first three numbers `1`, `2`, and `3` refer to 15-, 18-, and 25-state chromatin models, respectively, for `hg19` and `hg38` genome assemblies. The number `4` refers to the labels and colors used for the 15-state chromatin state model for `mm10`. This setting should match your epilogos analysis.

The `name` variable (set here to "My sample of interest") can be set to a descriptive name of your choice. It is recommended to make this short enough to read within the space provided in a WashU browser track label column -- perhaps no more than 20-25 characters.

The `url` variable is the web address of the bgzip-compressed epilogos result. While the `url` key contains the link to the bgzip file, as mentioned above, the tabix index should be in the same directory as the bzgip file.

Copy this JSON file to the directory associated with your web server. You can put this file into a subdirectory; the name of a subdirectory simply modifies the web address for the JSON file. You may give this file any name that you like; this simply changes the web address used to load the JSON-formatted datahub file.

For purposes of demonstration, we can call this file `mydatahub.json` and place it in the root of the web server directory. This would make the JSON file available at the address `http://192.168.0.1/mydatahub.json`. 

### Loading the datahub in the WashU browser

Now that you have generated the tabix-indexed qcat and datahub JSON assets and added them to your web server, you are ready to generate a web address that will load the WashU browser and instruct the browser to load your data.

The basic format of this address is:

#### &nbsp;&nbsp;http&#8288;://epigenomegateway.wustl.edu/browser/?genome=**build_name**&datahub=**datahub_address**

The *build_name* is replaced with one of `hg19`, `hg38`, or `mm10`, depending on which assembly was used to generate your epilogos dataset.

The *datahub_address* is replaced with the web address that points to your datahub JSON file.

Here is an example that specifies the `hg19` assembly and points to a hypothetical datahub file available at `http://192.168.0.1/mydatahub.json`:

#### &nbsp;&nbsp;http&#8288;://epigenomegateway.wustl.edu/browser/?genome=**hg19**&datahub=**http&#8288;://192.168.0.1/mydatahub.json**

You can open this link in a web browser, test it, modify it, and share it with others. This link will persist as long as you have your web server up and running, serving your datahub and qcat files.

Other parameters may be added to this address, which customize the behavior and appearance of the WashU browser. A more complete listing of track parameters is available from the WashU browser [wiki](http://wiki.wubrowse.org/URL_parameter).