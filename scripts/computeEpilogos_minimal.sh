#! /bin/bash

# set -e -o pipefail

usage() {
    echo -e "Usage:  $0 singleChromInputFile numStates metric outdir groupSpec [group2spec]"
    echo -e "where"
    echo -e "* singleChromInputFile consists of tab-delimited coordinates (chrom, start, stop)"
    echo -e "  and states observed at those loci for epigenome 1 (in column 4), epigenome 2 (in column 5), etc."
    echo -e "  States are positive integers.  singleChromInputFile must only contain data for a single chromosome."
    echo -e "* numStates is the number of possible states (the positive integers given in columns 4 and following)"
    echo -e "* metric is either 1 (S1), 2 (S2), or 3 (S3)"
    echo -e "* groupSpec specifies epigenomes to use in the analysis;"
    echo -e "  these are epigenome numbers corresponding to the order in which they appear in the input files"
    echo -e "  (e.g. 1, 2, 3 for the first three epigenomes, whose states are in columns 4, 5, 6)."
    echo -e "  groupSpec is a comma- and/or hyphen-delimited list (e.g. \"1,3,6-10,12,13\")."
    echo -e "  By default, the analysis will be performed on a single group of epigenomes."
    echo -e "* If an optional specification for a second group of epigenomes is provided as the last argument (group2spec),"
    echo -e "  then the analysis will be a comparison between the two groups of epigenomes."
    echo -e "* outdir will be created if necessary, and all results files will be written there."
    exit 2
}
    
if [[ $# != 5 ]] && [[ $# != 6 ]]; then # invalid number of arguments
    usage
fi    

# This script requires bgzip and starch (bedops).
BGZIP_EXE=`which bgzip 2> /dev/null`
if [ ! -x "$BGZIP_EXE" ]; then
    echo -e "Error:  Required external program bgzip was not found, or it is not executable."
    echo -e "Make sure this program is in your \$PATH, then try again."
    exit 2
fi
STARCH_EXE=`which starch 2> /dev/null`
if [ ! -x "$STARCH_EXE" ]; then
    echo -e "Error:  Required external program starch (part of bedops) was not found, or it is not executable."
    echo -e "Make sure this program is in your \$PATH, then try again."
    exit 2
fi

singleChromInputFile=$1
numStates=$2
metric=$3
outdir=$4
groupAspec=$5
groupBspec=""
if [ "$6" != "" ]; then
    groupBspec=$6
fi

# metric directs the programs to use S1, S2, or S3
if [ "$metric" != "1" ] && [ "$metric" != "2" ] && [ "$metric" != "3" ]; then
    echo -e "Error:  Invalid \"metric\" (\"$metric\") received."
    usage
fi
S1=0
S2=0
S3=0
if [ "$metric" == "1" ]; then
    S1=1
elif [ "$metric" == "2" ]; then
    S2=1
elif [ "$metric" == "3" ]; then
    S3=1
else
    echo -e "Error:  Mode of operation \"$metric\" not implemented."
    exit 2
fi

if [ ! -s $singleChromInputFile ]; then
    echo -e "Error:  File \"$singleChromInputFile\" was not found, or it is empty."
    exit 2
fi

EXE1=`which computeEpilogosPart1_perChrom 2> /dev/null`
EXE2=`which computeEpilogosPart2_perChrom 2> /dev/null`
EXE3=`which computeEpilogosPart3_perChrom 2> /dev/null`

if [ ! -x "$EXE1" ]; then
    echo -e "Error:  Required executable \"computeEpilogosPart1_perChrom\" was not found, or it is not executable."
    echo -e "Make sure you've run the \"make\" command from within the \"epilogos\" directory"
    echo -e "and added the path to computeEpilogosPart1_perChrom to your \$PATH environment variable."
    echo -e "The command \"which computeEpilogosPart1_perChrom\" (executed from any directory)"
    echo -e "will confirm that computeEpilogosPart1_perChrom is in your \$PATH by showing you its location,"
    echo -e "or return a message stating it could not be found if it has not been successfully added to your \$PATH."
    exit 2
fi
if [ ! -x "$EXE2" ]; then
    echo -e "Error:  Required executable \"computeEpilogosPart2_perChrom\" was not found, or it is not executable."
    echo -e "Make sure you've run the \"make\" command from within the \"epilogos\" directory"
    echo -e "and added the path to computeEpilogosPart2_perChrom to your \$PATH environment variable."
    echo -e "The command \"which computeEpilogosPart2_perChrom\" (executed from any directory)"
    echo -e "will confirm that computeEpilogosPart2_perChrom is in your \$PATH by showing you its location,"
    echo -e "or return a message stating it could not be found if it has not been successfully added to your \$PATH."
    exit 2
fi
if [ ! -x "$EXE3" ]; then
    echo -e "Error:  Required executable \"computeEpilogosPart3_perChrom\" was not found, or it is not executable."
    echo -e "Make sure you've run the \"make\" command from within the \"epilogos\" directory"
    echo -e "and added the path to computeEpilogosPart3_perChrom to your \$PATH environment variable."
    echo -e "The command \"which computeEpilogosPart3_perChrom\" (executed from any directory)"
    echo -e "will confirm that computeEpilogosPart3_perChrom is in your \$PATH by showing you its location,"
    echo -e "or return a message stating it could not be found if it has not been successfully added to your \$PATH."
    exit 2
fi

numStatesIsValid=`echo $numStates | awk '{if($0<1){print "0"; exit 1}else{print}}' | sed 's/[0-9]//g' | awk '{if(""==$0){print "1"}else{print "0"}}'`
if [ "0" == "$numStatesIsValid" ]; then
    echo -e "Error:  Invalid number of states specified (\"$numStates\").  This must be a positive integer."
    usage
fi

# Validate the group specifications
$EXE1 $groupAspec $groupBspec > /dev/null
if [ "$?" != 0 ]; then
    if [ "$groupBspec" == "" ]; then
	echo -e "Error in the group specification (\"$groupAspec\")."
    else
	echo -e "Error in one or both group specifications (\"$groupAspec\" and \"$groupBspec\")."
    fi
    usage
fi

groupBsize=""
groupAsize=`$EXE1 $groupAspec`
if [ "$groupBspec" != "" ]; then
    groupBsize=`$EXE1 $groupBspec`
fi

mkdir -p $outdir

# If all 3 final output files already exist, assume this is a re-run that isn't necessary.
if [ -s ${outdir}/observations.starch ] && [ -s ${outdir}/scores.txt.gz ] && [ -s ${outdir}/exemplarRegions.txt ]; then
    echo -e "All 3 final output files already exist and are non-empty:"
    echo -e "\t${outdir}/observations.starch"
    echo -e "\t${outdir}/scores.txt.gz"
    echo -e "\t${outdir}/exemplarRegions.txt"
    echo -e "Exiting, under the assumption that these files should not be recreated from scratch."
    echo -e "To recreate these files, move/rename them before running this program."
    exit 0
fi

# ----------------------------------
# echo -e "Executing \"part 1a\"..."
# ----------------------------------

INFILE_IS_COMPRESSED=`echo $singleChromInputFile | awk '{len=split($0,x,".");if("gz" == x[len]){print 1}else{print 0}}'`
if [ "$INFILE_IS_COMPRESSED" == "1" ]; then
    ZCAT_EXE=`which zcat 2> /dev/null`
    if [ ! -x "$ZCAT_EXE" ]; then
	echo -e "Error:  Failed to decompress input file $singleChromInputFile using zcat."
	echo -e "Try decompressing it yourself and specifying the decompressed version of it as input."
	exit 2
    fi
    chr=`zcat $singleChromInputFile | head -n 1 | cut -f1`
    infile1=${outdir}/uncompressedInfile_${chr}.txt
    zcat $singleChromInputFile > $infile1
    if [ ! -s $infile1 ]; then
	echo -e "Error:  Failed to decompress input file $singleChromInputFile using zcat."
	echo -e "Try decompressing it yourself and specifying the decompressed version of it as input."
	exit 2
    fi
else
    infile1=$singleChromInputFile
    chr=`head -n 1 $infile1 | cut -f1`
fi

outfileRand="" # used only if 2 groups of epigenomes are being compared
outfileQB=""   # used only if 2 groups of epigenomes are being compared
PfilenameString=""
randFilenameString=""

outfileNsites=${outdir}/${chr}_numSites.txt
if [ $S1 == 1 ]; then
    PfilenameString="_P1numerators.txt"
    outfileQ=${outdir}/${chr}_Q1numerators.txt
    if [ "$groupBspec" != "" ]; then
	randFilenameString="_randP1numerators.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}
	outfileQ=${outdir}/${chr}_Q1Anumerators.txt
	outfileQB=${outdir}/${chr}_Q1Bnumerators.txt
    fi
elif [ $S2 == 1 ]; then
    PfilenameString="_P2numerators.txt"
    outfileQ=${outdir}/${chr}_Q2numerators.txt
    if [ "$groupBspec" != "" ]; then
	randFilenameString="_randP2numerators.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}
	outfileQ=${outdir}/${chr}_Q2Anumerators.txt
	outfileQB=${outdir}/${chr}_Q2Bnumerators.txt
    fi
else
    PfilenameString="_statePairs.txt"
    outfileQ=${outdir}/${chr}_Q3Tallies.txt
    if [ "$groupBspec" != "" ]; then
	randFilenameString="_randomizedStatePairs.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}	    
	outfileQ=${outdir}/${chr}_Q3Atallies.txt
	outfileQB=${outdir}/${chr}_Q3Btallies.txt
    fi
fi
outfileP=${outdir}/${chr}${PfilenameString}



jobName="p1a_$chr"
if [[ ! -s $outfileP || ! -s $outfileQ || ! -s $outfileNsites || ($groupBspec != "" && (! -s $outfileQB || ! -s $outfileRand)) ]]; then
    $EXE1 $infile1 $metric $numStates $outfileP $outfileQ $outfileNsites $groupAspec $groupBspec $outfileRand $outfileQB > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
        exit 2
    fi
fi

# ----------------------------------
# echo -e "Executing \"part 1b\"..."
# ----------------------------------

# Add up the Q1 (or Q2 or Q3) tallies to get the genome-wide Q1 (or Q2 or Q3) matrices,
# or to get the Q1A and Q1A (or Q2A and Q2B, or Q3A and Q3B) matrices
# if two groups of epigenomes are being compared.

PID=$$
QQjobName=QQ_1b_${PID}
outfileQB=""
QBfilenameString=""
if [ $S1 == 1 ]; then
    outfileQ=${outdir}/Q1.txt
    QfilenameString="_Q1numerators.txt"
    if [ "$groupBspec" != "" ]; then
	outfileQ=${outdir}/Q1A.txt
	outfileQB=${outdir}/Q1B.txt
	QfilenameString="_Q1Anumerators.txt"
	QBfilenameString="_Q1Bnumerators.txt"
    fi
elif [ $S2 == 1 ]; then
    outfileQ=${outdir}/Q2.txt
    QfilenameString="_Q2numerators.txt"
    if [ "$groupBspec" != "" ]; then
	outfileQ=${outdir}/Q2A.txt
	outfileQB=${outdir}/Q2B.txt
	QfilenameString="_Q2Anumerators.txt"
	QBfilenameString="_Q2Bnumerators.txt"
    fi
else
    outfileQ=${outdir}/Q3.txt
    QfilenameString="_Q3Tallies.txt"
    if [ "$groupBspec" != "" ]; then
	outfileQ=${outdir}/Q3A.txt
	outfileQB=${outdir}/Q3B.txt
	QfilenameString="_Q3Atallies.txt"
	QBfilenameString="_Q3Btallies.txt"
    fi
fi

if [[ ! -s $outfileQ || ("$groupBspec" != "" && ! -s $outfileQB) ]]; then
    file=${outdir}/${chr}${QfilenameString}
    if [ ! -s $file ]; then
	echo -e "Error:  File "$file" is empty." > ${outdir}/${QQjobName}.stdout
	exit 2
    fi
    mv $file $outfileQ
    if [ "$groupBspec" != "" ]; then
	file=${outdir}/${chr}${QBfilenameString}
	if [ ! -s $file ]; then
	    echo -e "Error:  File "$file" is empty." > ${outdir}/${QQjobName}.stdout
	    exit 2
	fi
	mv $file $outfileQB
    fi
fi

totalNumSites=`cat $outfileNsites`
rm -f $outfileNsites

# ----------------------------------
# echo -e "Executing \"part 2a\"..."
# ----------------------------------

if [ "$INFILE_IS_COMPRESSED" == "1" ]; then
    rm -f $infile1
fi
infile=${outdir}/${chr}${PfilenameString}
infileQ=$outfileQ
infileQB=$outfileQB # empty ("") if only one group of epigenomes was specified
outfileObserved=${outdir}/${chr}_observed.txt
outfileScores=${outdir}/${chr}_scores.txt
outfileNulls=""
if [ "$groupBspec" != "" ]; then
    outfileNulls=${outdir}/${chr}_nulls.txt
fi
jobName="p2_$chr"

if [ ! -s $outfileObserved ]; then
    $EXE2 $infile $metric $totalNumSites $infileQ $outfileObserved $outfileScores $chr $infileQB > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
	exit 2
    else # clean up
	rm -f $infile
    fi
fi

if [ "$groupBspec" != "" ] && [ ! -s $outfileNulls ]; then
    infile=${outdir}/${chr}${randFilenameString}
    jobName="p2r_$chr"
    # The smaller number of arguments in the following call informs $EXE2 that it should only write the metric to $outfileNulls, with no additional info.
    $EXE2 $infile $metric $totalNumSites $infileQ $infileQB $outfileNulls > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
	exit 2
    else # clean up
	rm -f $infile
    fi
fi

# ----------------------------------
# echo -e "Executing \"part 2b\"..."
# ----------------------------------

cat ${outdir}/${chr}_scores.txt \
      | bgzip \
      > ${outdir}/scores.txt.gz
rm -f ${outdir}/${chr}_scores.txt $outfileQ $outfileQB
          
# ----------------------------------
# echo -e "Executing \"part 3\"..."
# ----------------------------------

infile=${chr}_observed.txt
if [ "$outfileNulls" != "" ]; then
    outfile=`echo $infile | sed 's/_observed.txt$/_withPvals.bed/g'`
else
    outfile=`echo $infile | sed 's/txt$/bed/g'`
fi
outfile=${outdir}/$outfile
finalOutfile=${outdir}/observations.starch
jobName="p3_$chr"

if [ ! -s $outfile ]; then
    if [ "$outfileNulls" != "" ]; then
      infile=${outdir}/$infile
      $EXE3 $infile $outfileNulls $outfile > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
      if [ $? == 0 ]; then
          rm -f $infile $outfileNulls
      else
	  exit 2
      fi
   else
      infile=${outdir}/$infile
      mv $infile $outfile
    fi
    starch $outfile > $finalOutfile 2> ${outdir}/${jobName}.stderr
    if [ $? == 0 ]; then
	rm -f $outfile
    fi
fi

exemplarRegions=${outdir}/exemplarRegions.txt
if [ -s $finalOutfile ] && [ ! -s $exemplarRegions ]; then
    STATE_COLUMN=4
    if [ $S1 == 1 ]; then
	SCORE_COLUMN=7
    else
	SCORE_COLUMN=10
    fi
    unstarch $finalOutfile \
	 | awk -v stateCol=$STATE_COLUMN -v scoreCol=$SCORE_COLUMN 'BEGIN {OFS="\t"; chrom=""; state=""; bestScore=""; line=""} \
           { \
              if ( $(stateCol) != state || $1 != chrom ) { \
                 if ( state != "" ) { \
                    print line; \
                 } \
                 chrom = $1; \
                 state = $(stateCol); \
                 bestScore = $(scoreCol); \
                 line = $0; \
              } else { \
                 if ( $(scoreCol) >= bestScore ) { \
                    bestScore = $(scoreCol); \
                    line = $0; \
                 } \
              } \
           } END { \
                      if ( line != "" ) { \
                         print line; \
                      } \
                 }' \
	 | sort -gr -k${SCORE_COLUMN},${SCORE_COLUMN} \
	 > $exemplarRegions
fi

exit 0
