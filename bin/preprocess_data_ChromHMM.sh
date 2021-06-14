#!/bin/bash
set -ueo pipefail

if [ $# -lt 3 ]; then
  echo
  echo "Usage: $0 <datadir> <metadata> <chromsizes>"
  echo "Insufficient information -- please provide:"
  echo
  echo "1. <datadir>: Directory with input data, one file per biosample-chromosome combination."
  echo "   These can be obtained from ChromHMM using the '-printstatebyline' parameter."
  echo "   Each line expected to represent a single 200bp genomic bin and its chromatin state."
  echo "   For more information on the formatting of these files, see the ChromHMM manual at:"
  echo "     http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf"
  echo "   Example files can be found in the 'data/ChromHMM' directory."
  echo
  echo "2. <metadata>: File containing metadata information for each input dataset/biosample."
  echo "   This file is expected to contain a header line, and a unique biosample name in column 1,"
  echo "   represented in the input data directory. Example file: 'data/metadata_Boix.txt'."
  echo
  echo "3. <chromsizes>: File containing chromosome sizes in number of basepairs."
  echo "   This information can be obtained from UCSC, e.g. for human genome hg19, run:"
  echo "   $ wget -q -O - ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz | \\"
  echo "       gunzip - | cut -f 1-2 > hg19.genome"
  echo "   Examples for hg19 and hg38 can be found in the 'data' directory."
  echo
  echo "Example run: $0 data/ChromHMM data/metadata_Boix.txt data/hg19.genome"
  exit
fi  

datadir=${1}
metadata=${2}
chromsizes=${3}

for chr in $(cut -f1 ${chromsizes}); do
  echo -n "Processing ${chr}: "
  cmd="paste"; 
  i=0;
  for biosample in $(cut -f1 ${metadata} | tail -n +2); do
    if compgen -G "${datadir}/*${biosample}*${chr}_*.txt*" > /dev/null; then
      nam=`stat -t -c %n ${datadir}/*${biosample}*${chr}_*.txt*`
      cmd="$cmd <(zcat -f -- "${nam}")";
      i=$((i+1))
    fi
  done
  echo -n "${i} files found. "
  if [ ${i} -gt 0 ]; then
    eval $cmd | awk 'BEGIN{OFS="\t";chr=""}
      {if (NR==1) { chr=$2 } else { if (NR>2) print chr, (NR-3)*200, (NR-2)*200, $0}}' \
      > matrix_${chr}.txt
    echo "Done."
  else
    echo "Skipping."
  fi
done



