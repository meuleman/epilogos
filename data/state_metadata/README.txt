epilogos-metadata
=================

This directory contains metadata tables describing chromatin states for current datasets.

Data are organized by organism, dataset, genome assembly, and state model (15, 18, or 25, where applicable or available). 

Each directory contains a tab-delimited file called `metadata.tsv`. Columns contain zero- and one-based indices, short and long names, and hexadecimal, RGBA, and HTML color names.

*Note*: Tick marks or single-quote characters in short and long names are replaced with the letter "p", to avoid escaping problems.

```bash
state_metadata
|
+-- README.txt
|
+-- human
|
|  +-- Adsera_et_al_833_sample
|    +-- hg19
|      +-- 18
|    +-- hg38
|      +-- 18
|
|  +-- Roadmap_Consortium_127_sample
|    +-- hg19
|      +-- 15
|      +-- 18
|      +-- 25
|    +-- hg38
|      +-- 15
|      +-- 18
|      +-- 25
|
+-- mouse
|  +-- Gorkin_et_al_65_sample
|      +-- 15
```

## References

Human
-----

For human, 15-, 18- and 25-state model metadata are derived from the Roadmap Consortium portal:

https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
https://egg2.wustl.edu/roadmap/web_portal/imputed.html

Mouse
-----

For mouse, 15-state model metadata are derived from Gorkin, et al.:

https://doi.org/10.1038/s41586-020-2093-3

