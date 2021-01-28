#!/bin/bash

GENOME=$1
METADATA=$2
DIR="~meuleman/work/projects/imputation//WM20191003_curated_data_mixed/public_ChromHMM_released/observed_aux_18_hg19/CALLS/STATEBYLINE/"

for CHR in $(cut -f1 ${GENOME}.genome); do
  echo "Processing ${CHR}..."
  CMD="paste"; 
  for BIOSAMPLE in $(cut -f1 ${METADATA} | tail -n +2); do
    CMD="$CMD <(zcat ${DIR}/*${BIOSAMPLE}*${CHR}_*statebyline.txt.gz)"; 
  done
  eval $CMD | awk 'BEGIN{OFS="\t";chr=""}
    {if (NR==1) { chr=$2 } else { if (NR>2) print chr, (NR-3)*200, (NR-2)*200, $0}}' \
    > matrix_${CHR}.txt
done



