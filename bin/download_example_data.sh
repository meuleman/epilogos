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
    echo "Downloading Boix state metadata..."
    mkdir -p data/state_metadata/human/Boix_et_al_833_sample/hg19/18/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv --output data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv
    echo
else
    echo "FILE NOT DOWNLOADED: data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv already exists"
    echo ""
fi

if [[ ! -f data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/metadata.tsv ]]
then
    echo "Downloading Roadmap state metadata..."
    mkdir -p data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/metadata.tsv --output data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/metadata.tsv
    echo
else
    echo "FILE NOT DOWNLOADED: data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/metadata.tsv already exists"
    echo ""
fi

if [[ ! -f data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz ]]
then
    echo "Downloading similarity search data..."
    mkdir -p data/simsearch/male/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz --output data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz
    echo
else
    echo "FILE NOT DOWNLOADED: data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz already exists"
    echo ""
fi

if [[ ! -f data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz ]]
then
    echo "Downloading plot region data..."
    mkdir -p data/simsearch/male/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz --output data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz
    echo
else
    echo "FILE NOT DOWNLOADED: data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz already exists"
    echo ""
fi

if [[ ! -f data/plotregion/male/regionsOfInterest_male_s1_chr1.bed ]]
then
    echo "Downloading plot region data..."
    mkdir -p data/simsearch/male/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/plotregion/male/regionsOfInterest_male_s1_chr1.bed --output data/plotregion/male/regionsOfInterest_male_s1_chr1.bed
    echo
else
    echo "FILE NOT DOWNLOADED: data/plotregion/male/regionsOfInterest_male_s1_chr1.bed already exists"
    echo ""
fi

if [[ ! -f data/plotregion/male_vs_female/regionsOfInterest_Male_vs_Female_s1.bed ]]
then
    echo "Downloading plot region data..."
    mkdir -p data/plotregion/male_vs_female/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/plotregion/male_vs_female/regionsOfInterest_Male_vs_Female_s1.bed --output data/plotregion/male_vs_female/regionsOfInterest_Male_vs_Female_s1.bed
    echo
else
    echo "FILE NOT DOWNLOADED: data/plotregion/male_vs_female/regionsOfInterest_Male_vs_Female_s1.bed already exists"
    echo ""
fi

if [[ ! -f data/plotregion/male_vs_female/scores_Male_s1_matrix_chrX.txt.gz ]]
then
    echo "Downloading plot region data..."
    mkdir -p data/plotregion/male_vs_female/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/plotregion/male_vs_female/scores_Male_s1_matrix_chrX.txt.gz --output data/plotregion/male_vs_female/scores_Male_s1_matrix_chrX.txt.gz
    echo
else
    echo "FILE NOT DOWNLOADED: data/plotregion/male_vs_female/scores_Male_s1_matrix_chrX.txt.gz already exists"
    echo ""
fi

if [[ ! -f data/plotregion/male_vs_female/scores_Female_s1_matrix_chrX.txt.gz ]]
then
    echo "Downloading plot region data..."
    mkdir -p data/plotregion/male_vs_female/
    curl -L https://github.com/meuleman/epilogos/raw/main/data/plotregion/male_vs_female/scores_Female_s1_matrix_chrX.txt.gz --output data/plotregion/male_vs_female/scores_Female_s1_matrix_chrX.txt.gz
    echo
else
    echo "FILE NOT DOWNLOADED: data/plotregion/male_vs_female/scores_Female_s1_matrix_chrX.txt.gz already exists"
    echo ""
fi

echo ""
echo "Find files at:"
echo "data/pyData/male/epilogos_matrix_chr1.txt.gz"
echo "data/pyData/female/epilogos_matrix_chr1.txt.gz"
echo "data/state_metadata/human/Boix_et_al_833_sample/hg19/18/metadata.tsv"
echo "data/state_metadata/human/Roadmap_Consortium_127_sample/hg19/18/metadata.tsv"
echo "data/simsearch/male/scores_male_s1_matrix_chr1.txt.gz"
echo "data/plotregion/male/scores_male_s1_matrix_chr1.txt.gz"
echo "data/plotregion/male/regionsOfInterest_male_s1_chr1.bed"
echo "data/plotregion/male_vs_female/regionsOfInterest_Male_vs_Female_s1.bed"
echo "data/plotregion/male_vs_female/scores_Male_s1_matrix_chrX.txt.gz"
echo "data/plotregion/male_vs_female/scores_Female_s1_matrix_chrX.txt.gz"
