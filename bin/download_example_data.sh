#!/bin/bash
# Script to download example data for epilogos

mkdir -p data/pyData/male/
mkdir -p data/pyData/female/
mkdir -p data/state_metadata/human/Boix_et_al_833_sample/hg19/18/

curl -L https://github.com/meuleman/epilogos/blob/main/data/pyData/male/epilogos_matrix_chr1.txt.gz?raw=true --output data/pyData/male/epilogos_matrix_chr1.txt.gz
curl -L https://github.com/meuleman/epilogos/blob/main/data/pyData/female/epilogos_matrix_chr1.txt.gz?raw=true --output data/pyData/female/epilogos_matrix_chr1.txt.gz
curl -L https://github.com/meuleman/epilogos/raw/main/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv --output data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv

echo "Find files at:"
echo "data/pyData/male/epilogos_matrix_chr1.txt.gz"
echo "data/pyData/female/epilogos_matrix_chr1.txt.gz"
echo "data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv"
