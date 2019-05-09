1. [Visualizing a qcat file](#visualizing-a-qcat-file)
    * [Setting up a web server](#setting-up-a-web-server)
    * [Preparing tabix-indexed qcat files](#preparing-tabix-indexed-qcat-files)
        - [Compressing](#compressing)
        - [Indexing](#indexing)
    * [Creating a datahub](#creating-a-datahub)
        - [Datahub template](#datahub-template)
    * [Loading the datahub in the WashU browser](#loading-the-datahub-in-the-washu-browser)
2. [Qcat file specification](#qcat-file-specification)

## Visualizing a qcat file (not recommended)

This section describes how to visualize a [qcat](#qcat-file-specification) result file ("epilogos"), using the [WashU Epigenome Browser](https://epigenomegateway.wustl.edu/). **We strongly recommend using [HiGlass](https://higlass.io) to visualize results via the output file** `scores.txt.gz`**.**

In broad terms, this can be done by way of:

1. Converting output file `scores.txt.gz` to qcat format
2. Setting up a public web server
3. Indexing the qcat result with tabix
4. Creating a JSON-formatted *datahub* file
5. Copying the qcat file, tabix index, and JSON datahub assets to your web server
6. Loading the WashU browser with your custom datahub

### Converting output file `scores.txt.gz` to qcat format

Use the Python script `qcatCreator.py` in the `scripts` subdirectory to convert `scores.txt.gz` to qcat format:

```bash
$ bgzip -cd yourDirectory/scores.txt.gz | qcatCreator.py | bgzip > yourDirectory.qcat.bed.gz
```

You will need to have python (either 2.X or 3.X) available on your system. If the above command does not work, you might need to add permission to execute qcatCreator.py by executing the command `chmod +x` on `qcatCreator.py`.

### Setting up a web server

A public-facing web server is a basic requirement for visualizing your qcat data with the WashU browser. Your web server will make your qcat-formatted epilogos data available, as well as related metadata about the epilogos track, chromatin states, and optional annotation data.

This web server must be able to serve files via port `tcp/80` (http). You will need the public-facing static IP address assigned to this web server; in other words, you must be able to access this server through firewalls and outside your local network.

Your institution or research facility may already offer a public-facing web service, and contacting your local IT help may get you information on this option to get you started, such as where you would store your files and, importantly, the web address of the server.

Depending on what resources are available locally, there are also (fee-based) web hosting services, such as [Dreamhost](https://www.dreamhost.com/) or [Amazon Lightsail](https://lightsail.aws.amazon.com), for instance. Management tools for these services are web-based and make configuration easy.

For the purposes of this section, we will assume that you have a working, public-facing web server that is available at `http://192.168.0.1` (your actual IP address will be different), and that you are able to copy or upload files to the required directory on this server, so that these files are available through this address.

### Preparing tabix-indexed qcat files

The next step is to index the `bgzip`-compressed qcat file with `tabix`. `tabix`, along with `bgzip`, is available through the larger `htslib` software package. This can be installed through package managers like [Homebrew](https://brew.sh/) or [Bioconda](https://bioconda.github.io/), or you can download and compile source from [Github](https://github.com/samtools/htslib/blob/develop/INSTALL). You should be able to run `bgzip --help` and `tabix --help` or similar to verify their availability. The following command makes an index file called `qcat.bed.bgz.tbi`:

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

## Qcat file specification

After extraction with `gunzip` or `bgzip`, the qcat output file is a BED4 (four-column BED) formatted file.

The first three columns of the BED4 file represent the genomic interval or bin that contains [histone modifications used to generate chromatin state calls](http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html). Epilogos allows bins of any size; however, bins will be generally be 200 nt in size, representing the distance between nucleosomes.

The fourth column contains a JSON-like string in "qcat" or "quantitative category" format, which describes a 15-, 18- or 25-element array of "state: value" pairs. These pairs are the chromatin states ordered by ascending signed per-state contribution to the overall epilogos score.  When running Epilogos on a single group of epigenomes, the overall score (KL, KL\*, or KL\**) is the sum of the individual per-state values; when comparing two groups of epigenomes, the overall score (DKL, DKL\*, or DKL\*\*) is the sum of the absolute values of the per-state values.

A complete description of the quantitative category format is provided at the [WashU Epigenome Browser wiki](http://wiki.wubrowse.org/QuantitativeCategorySeries).

## Support

To get additional help or if you have questions about this software, open an [issue ticket](https://github.com/Altius/epilogos/issues).
