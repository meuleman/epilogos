#! /bin/bash

# set -e -o pipefail

usage() {
    echo -e "Usage:  $0 singleChromInputFile measurementType numStates outdir groupSpec [group2spec]"
    echo -e "where"
    echo -e "* singleChromInputFile consists of tab-delimited coordinates (chrom, start, stop)"
    echo -e "  and states observed at those loci for epigenome 1 (in column 4), epigenome 2 (in column 5), etc."
    echo -e "  States are positive integers.  singleChromInputFile must only contain data for a single chromosome."
    echo -e "* measurementType is either 0 (to use KL), 1 (to use KL*), or 2 (to use KL**)"
    echo -e "* numStates is the number of possible states (the positive integers given in columns 4 and following)"
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
measurementType=$2
numStates=$3
outdir=$4
group1spec=$5
group2spec=""
if [ "$6" != "" ]; then
    group2spec=$6
fi

# measurementType directs the programs to use KL, KL*, or KL**
if [ "$measurementType" != "0" ] && [ "$measurementType" != "1" ] && [ "$measurementType" != "2" ]; then
    echo -e "Error:  Invalid \"measurementType\" (\"$measurementType\") received."
    usage
fi
KL=0
KLs=0
KLss=0
if [ "$measurementType" == "0" ]; then
    KL=1   # KL
elif [ "$measurementType" == "1" ]; then
    KLs=1  # KL*
else
    KLss=1 # KL**
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
$EXE1 $group1spec $group2spec > /dev/null
if [ "$?" != 0 ]; then
    if [ "$group2spec" == "" ]; then
	echo -e "Error in the group specification (\"$group1spec\")."
    else
	echo -e "Error in one or both group specifications (\"$group1spec\" and \"$group2spec\")."
    fi
    usage
fi

group2size=""
group1size=`$EXE1 $group1spec`
if [ "$group2spec" != "" ]; then
    group2size=`$EXE1 $group2spec`
fi

mkdir -p $outdir

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
outfileQ2=""   # used only if 2 groups of epigenomes are being compared
PfilenameString=""
randFilenameString=""

outfileNsites=${outdir}/${chr}_numSites.txt
if [ $KL == 1 ]; then
    PfilenameString="_Pnumerators.txt"
    outfileQ=${outdir}/${chr}_Qnumerators.txt
    if [ "$group2spec" != "" ]; then
	randFilenameString="_randPnumerators.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}
	outfileQ=${outdir}/${chr}_Q1numerators.txt
	outfileQ2=${outdir}/${chr}_Q2numerators.txt
    fi
elif [ $KLs == 1 ]; then
    PfilenameString="_PsNumerators.txt"
    outfileQ=${outdir}/${chr}_QsNumerators.txt
    if [ "$group2spec" != "" ]; then
	randFilenameString="_randPsNumerators.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}
	outfileQ=${outdir}/${chr}_Qs1numerators.txt
	outfileQ2=${outdir}/${chr}_Qs2numerators.txt
    fi
else
    PfilenameString="_statePairs.txt"
    outfileQ=${outdir}/${chr}_QssTallies.txt
    if [ "$group2spec" != "" ]; then
	randFilenameString="_randomizedStatePairs.txt"
	outfileRand=${outdir}/${chr}${randFilenameString}	    
	outfileQ=${outdir}/${chr}_Qss1tallies.txt
	outfileQ2=${outdir}/${chr}_Qss2tallies.txt
    fi
fi
outfileP=${outdir}/${chr}${PfilenameString}



jobName="p1a_$chr"
if [[ ! -s $outfileP || ! -s $outfileQ || ! -s $outfileNsites || ($group2spec != "" && (! -s $outfileQ2 || ! -s $outfileRand)) ]]; then
    $EXE1 $infile1 $measurementType $numStates $outfileP $outfileQ $outfileNsites $group1spec $group2spec $outfileRand $outfileQ2 > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
        exit 2
    fi
fi

# ----------------------------------
# echo -e "Executing \"part 1b\"..."
# ----------------------------------

# Add up the Q (or Q* or Q**) tallies to get the genome-wide Q (or Q* or Q**) matrices,
# or to get the Q1 and Q2 (or Q1* and Q2*, or Q1** and Q2**) matrices
# if two groups of epigenomes are being compared.

PID=$$
QQjobName=QQ_1b_${PID}
outfileQ2=""
Q2filenameString=""
if [ $KL == 1 ]; then
    outfileQ=${outdir}/Q.txt
    QfilenameString="_Qnumerators.txt"
    if [ "$group2spec" != "" ]; then
	outfileQ=${outdir}/Q1.txt
	outfileQ2=${outdir}/Q2.txt
	QfilenameString="_Q1numerators.txt"
	Q2filenameString="_Q2numerators.txt"
    fi
elif [ $KLs == 1 ]; then
    outfileQ=${outdir}/Qs.txt
    QfilenameString="_QsNumerators.txt"
    if [ "$group2spec" != "" ]; then
	outfileQ=${outdir}/Qs1.txt
	outfileQ2=${outdir}/Qs2.txt
	QfilenameString="_Qs1numerators.txt"
	Q2filenameString="_Qs2numerators.txt"
    fi
else
    outfileQ=${outdir}/Qss.txt
    QfilenameString="_QssTallies.txt"
    if [ "$group2spec" != "" ]; then
	outfileQ=${outdir}/Qss1.txt
	outfileQ2=${outdir}/Qss2.txt
	QfilenameString="_Qss1tallies.txt"
	Q2filenameString="_Qss2tallies.txt"
    fi
fi

if [[ ! -s $outfileQ || ("$group2spec" != "" && ! -s $outfileQ2) ]]; then
    file=${outdir}/${chr}${QfilenameString}
    if [ ! -s $file ]; then
	echo -e "Error:  File "$file" is empty." > ${outdir}/${QQjobName}.stdout
	exit 2
    fi
    mv $file $outfileQ
    if [ "$group2spec" != "" ]; then
	file=${outdir}/${chr}${Q2filenameString}
	if [ ! -s $file ]; then
	    echo -e "Error:  File "$file" is empty." > ${outdir}/${QQjobName}.stdout
	    exit 2
	fi
	mv $file $outfileQ2
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
infileQ2=$outfileQ2 # empty ("") if only one group of epigenomes was specified
outfileObserved=${outdir}/${chr}_observed.txt
outfileQcat=${outdir}/${chr}_qcat.txt
outfileNulls=""
if [ "$group2spec" != "" ]; then
    outfileNulls=${outdir}/${chr}_nulls.txt
fi
jobName="p2_$chr"

if [ ! -s $outfileObserved ]; then
    $EXE2 $infile $measurementType $totalNumSites $infileQ $outfileObserved $outfileQcat $chr $infileQ2 > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
	exit 2
    else # clean up
	rm -f $infile
    fi
fi

if [ "$group2spec" != "" ] && [ ! -s $outfileNulls ]; then
    infile=${outdir}/${chr}${randFilenameString}
    jobName="p2r_$chr"
    # The smaller number of arguments in the following call informs $EXE2 that it should only write the metric to $outfileNulls, with no additional info.
    $EXE2 $infile $measurementType $totalNumSites $infileQ $infileQ2 $outfileNulls > ${outdir}/${jobName}.stdout 2> ${outdir}/${jobName}.stderr
    if [ $? != 0 ]; then
	exit 2
    else # clean up
	rm -f $infile
    fi
fi

# ----------------------------------
# echo -e "Executing \"part 2b\"..."
# ----------------------------------

cat ${outdir}/${chr}_qcat.txt \
      | sort -k1b,1 -s \
      | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"id:"NR",qcat:"$4}' \
      | bgzip \
      > ${outdir}/qcat.bed.gz
rm -f ${outdir}/${chr}_qcat.txt
          
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
    if [ $KL == 1 ]; then
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
