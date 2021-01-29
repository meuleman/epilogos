#!/bin/bash

### The ${GENOME}.genome file contains chromosome sizes, and can be obtained in several ways, including:
### wget -q -O - ftp://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/database/chromInfo.txt.gz | gunzip - | cut -f 1-2
GENOME=$1

### the METADATA file is expected to contain a header line
METADATA=$2

### Directory where all STATEBYLINE ChromHMM output is stored, one file per biosample-chromosome combination.
DIR="/home/meuleman/work/projects/imputation//WM20191003_curated_data_mixed/public_ChromHMM_released/observed_aux_18_hg19/CALLS/STATEBYLINE/"

for CHR in $(cut -f1 ${GENOME}.genome); do
  echo -n "Processing ${CHR}: "
  CMD="paste"; 
  NUMFILES=`ls -1 ${DIR}/*${CHR}_*.txt* 2>/dev/null | wc -l`
  if [ ${NUMFILES} -gt 0 ];
  then
    echo "${NUMFILES} files"
  for BIOSAMPLE in $(cut -f1 ${METADATA} | tail -n +2); do
    NAM=`stat -t -c %n ${DIR}/*${BIOSAMPLE}*${CHR}_*.txt*`
    CMD="$CMD <(zcat -f -- ${NAM})";
  done
  eval $CMD | awk 'BEGIN{OFS="\t";chr=""}
    {if (NR==1) { chr=$2 } else { if (NR>2) print chr, (NR-3)*200, (NR-2)*200, $0}}' \
    > matrix_${CHR}.txt
  else
    echo "no files found"
  fi
done



