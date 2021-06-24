#!/bin/bash
# Script to download example data for epilogos

if [[ ! -f data/pyData/male/epilogos_matrix_chr1.txt.gz ]]
then
    echo "Downloading male data..."
    mkdir -p data/pyData/male/
    curl -L https://github.com/meuleman/epilogos/blob/main/data/pyData/male/epilogos_matrix_chr1.txt.gz?raw=true --output data/pyData/male/epilogos_matrix_chr1.txt.gz
    echo ""
else
    echo "FILE NOT DOWNLOADED: data/pyData/male/epilogos_matrix_chr1.txt.gz already exists"
    echo ""
fi

if [[ ! -f data/pyData/female/epilogos_matrix_chr1.txt.gz ]]
then
    echo "Downloading female data..."
    mkdir -p data/pyData/female/
    curl -L https://github.com/meuleman/epilogos/blob/main/data/pyData/female/epilogos_matrix_chr1.txt.gz?raw=true --output data/pyData/female/epilogos_matrix_chr1.txt.gz
    echo ""
else
    echo "FILE NOT DOWNLOADED: data/pyData/female/epilogos_matrix_chr1.txt.gz already exists"
    echo ""
fi

if [[ ! -f data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv ]]
then
    echo "Downloading state metadata..."
    mkdir -p data/state_metadata/human/Boix_et_al_833_sample/hg19/18/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv --output data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv
    echo
else
    echo "FILE NOT DOWNLOADED: data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv already exists"
    echo ""
fi

echo ""
echo "Find files at:"
echo "data/pyData/male/epilogos_matrix_chr1.txt.gz"
echo "data/pyData/female/epilogos_matrix_chr1.txt.gz"
echo "data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv"
