#! /bin/bash

# set -e -o pipefail

usage() {
    echo -e "Usage:  $0 fileOfPerChromFilenames measurementType numStates outdir groupSpec [group2spec]"
    echo -e "where"
    echo -e "* fileOfPerChromFilenames contains one input filename per line (including path if needed);"
    echo -e "  each chrom-specific infile consists of tab-delimited coordinates (chrom, start, stop)"
    echo -e "  and states observed at those loci for epigenome 1 (in column 4), epigenome 2 (in column 5), etc."
    echo -e "  States are positive integers."
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

fileOfFilenames=$1
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
    
if [ ! -s "$fileOfFilenames" ]; then
    echo -e "Error:  File \"$fileOfFilenames\" was not found, or it is empty."
    exit 2
fi

EXE1=`which computeEpilogosPart1_perChrom`
EXE2=`which computeEpilogosPart2_perChrom`
EXE3=`which computeEpilogosPart3_perChrom`

if [ ! -x "$EXE1" ]; then
    echo -e "Error:  Required executable \"$EXE1\" not found, or it is not executable."
    exit 2
fi
if [ ! -x "$EXE2" ]; then
    echo -e "Error:  Required executable \"$EXE2\" not found, or it is not executable."
    exit 2
fi
if [ ! -x "$EXE3" ]; then
    echo -e "Error:  Required executable \"$EXE3\" not found, or it is not executable."
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

outfileRand="" # used only if 2 groups of epigenomes are being compared
outfileQ2=""   # used only if 2 groups of epigenomes are being compared
PfilenameString=""
randFilenameString=""
dependencies="afterok" # SLURM syntax
linenum=0
# Important:
# Need to code this as follows (i.e., "while read line do ... done <<< input"
# instead of "input | while read line do ... done"),
# because otherwise the body of the loop gets executed in a subshell (due to the '|'),
# and the updates to variable "dependencies" aren't passed back to the parent shell
# (i.e., they're lost after the "done" statement).
while read line
do
    ((linenum++))
    file=$line
    if [ ! -s $file ]; then
	echo -e "Error:  File \"$file\", on line $linenum of $fileOfFilenames, was not found, or it is empty."
	exit 2
    fi
    chr=`head -n 1 $file | cut -f1`
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
    memSize="50M"
    if [[ ! -s $outfileP || ! -s $outfileQ || ! -s $outfileNsites || ($group2spec != "" && (! -s $outfileQ2 || ! -s $outfileRand)) ]]; then
    thisJobID=$(sbatch --parsable --partition=queue1 --job-name=$jobName --output=${outdir}/${jobName}.o%j --error=${outdir}/${jobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
       $EXE1 $file $measurementType $numStates $outfileP $outfileQ $outfileNsites $group1spec $group2spec $outfileRand $outfileQ2
       if [ \$? != 0 ]; then
          exit 2
       fi
EOF
    )
    dependencies="${dependencies}:${thisJobID}"
    fi

done <<< "$(cat $fileOfFilenames)" # see note above regarding this syntax

# ----------------------------------
# echo -e "Executing \"part 1b\"..."
# ----------------------------------

# Add up the Q (or Q* or Q**) tallies to get the genome-wide Q (or Q* or Q**) matrices,
# or to get the Q1 and Q2 (or Q1* and Q2*, or Q1** and Q2**) matrices
# if two groups of epigenomes are being compared.

PID=$$
QQjobName=QQ_1b_${PID}
memSize="50M" # more than enough
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
totalNumSites=0 # placeholder

tempfile=${outdir}/tempfile_Q_${PID}.txt
TEMP1=${outdir}/tempfile1_Q_${PID}.txt
TEMP2=${outdir}/tempfile2_Q_${PID}.txt
if [ $dependencies == "afterok" ]; then
    dependencyString=""
else
    dependencyString="--dependency="$dependencies
fi

dependencyString2=""
if [[ ! -s $outfileQ || ("$group2spec" != "" && ! -s $outfileQ2) ]]; then
    QQjobID=$(sbatch --parsable --partition=queue1 $dependencyString --job-name=$QQjobName --output=${outdir}/${QQjobName}.o%j --error=${outdir}/${QQjobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
firstFile=\`ls -1 ${outdir}/chr*${QfilenameString} | head -n 1\`
if [ ! -s \$firstFile ]; then
   echo -e "Error:  File "\$firstFile" is empty."
   exit 2
fi
numRows=\`wc -l < \$firstFile\`
awk '{print NF}' \$firstFile | sort -n | uniq > $tempfile
if [ \`wc -l < $tempfile\` != "1" ]; then
   echo -e "Error:  File "\$firstFile" has at least two rows with different numbers of columns."
   exit 2
fi
numCols=\`cat $tempfile\`

cp \$firstFile $TEMP1
for file in \`ls -1 ${outdir}/chr*${QfilenameString} | tail -n +2\`
do
   if [ ! -s \$file ]; then
      echo -e "Error:  File "\$file" is empty."
      exit 2
   fi
   thisNumRows=\`wc -l < \$file\`
   awk '{print NF}' \$file | sort -n | uniq > $tempfile
   if [ \`wc -l < $tempfile\` != "1" ]; then
      echo -e "Error:  File "\$file" has at least two rows with different numbers of columns."
      exit 2
   fi
   thisNumCols=\`cat $tempfile\`
   if [ \$thisNumRows != \$numRows ]; then
      echo -en "Error:  Files "\$firstFile" and "\$file" have different numbers of rows ("
      echo \$numRows" and "\$thisNumRows")."
      exit 2
   fi
   if [ \$thisNumCols != \$numCols ]; then
      echo -en "Error:  Files "\$firstFile" and "\$file" have different numbers of columns ("
      echo \$numCols" and "\$thisNumCols")."
      exit 2
   fi

   paste $TEMP1 \$file \
      | awk -v nCols=\$numCols '{printf("%d", \$1 + \$(nCols + 1)); \
                                 for (i=2; i<= nCols; i++){ \
                                    printf("\t%d", \$i + \$(nCols + i)); \
                                 } \
                                 printf("\n"); \
                                }' \
      > $TEMP2
   mv $TEMP2 $TEMP1
done
cp $TEMP1 $outfileQ
rm -f ${outdir}/chr*${QfilenameString}

if [ "$group2spec" != "" ]; then
   firstFile=\`ls -1 ${outdir}/chr*${Q2filenameString} | head -n 1\`
   if [ ! -s \$firstFile ]; then
      echo -e "Error:  File "\$firstFile" is empty."
      exit 2
   fi
   numRows=\`wc -l < \$firstFile\`
   awk '{print NF}' \$firstFile | sort -n | uniq > $tempfile
   if [ \`wc -l < $tempfile\` != "1" ]; then
      echo -e "Error:  File "\$firstFile" has at least two rows with different numbers of columns."
      exit 2
   fi
   numCols=\`cat $tempfile\`

   cp \$firstFile $TEMP1
   for file in \`ls -1 ${outdir}/chr*${Q2filenameString} | tail -n +2\`
   do
      if [ ! -s \$file ]; then
         echo -e "Error:  File "\$file" is empty."
         exit 2
      fi
      thisNumRows=\`wc -l < \$file\`
      awk '{print NF}' \$file | sort -n | uniq > $tempfile
      if [ \`wc -l < $tempfile\` != "1" ]; then
         echo -e "Error:  File "\$file" has at least two rows with different numbers of columns."
         exit 2
      fi
      thisNumCols=\`cat $tempfile\`
      if [ \$thisNumRows != \$numRows ]; then
         echo -en "Error:  Files "\$firstFile" and "\$file" have different numbers of rows ("
         echo \$numRows" and "\$thisNumRows")."
         exit 2
      fi
      if [ \$thisNumCols != \$numCols ]; then
         echo -en "Error:  Files "\$firstFile" and "\$file" have different numbers of columns ("
         echo \$numCols" and "\$thisNumCols")."
         exit 2
      fi

      paste $TEMP1 \$file \
         | awk -v nCols=\$numCols '{printf("%d", \$1 + \$(nCols + 1)); \
                                    for (i=2; i<= nCols; i++){ \
                                       printf("\t%d", \$i + \$(nCols + i)); \
                                    } \
                                    printf("\n"); \
                                   }' \
         > $TEMP2
      mv $TEMP2 $TEMP1
   done
   cp $TEMP1 $outfileQ2
   rm -f ${outdir}/chr*${Q2filenameString}
fi

rm -f $tempfile $TEMP1 $TEMP2

EOF
	   )
    dependencyString2="--dependency=afterok:"${QQjobID}
fi

# ----------------------------------
# echo -e "Executing \"part 1c\"..."
# ----------------------------------

totalSitesJobName=QQ_1c_${PID}
memSize="5M" # more than enough
totalSitesJobID=$(sbatch --parsable --partition=queue1 $dependencyString2 --job-name=$totalSitesJobName --output=${outdir}/${totalSitesJobName}.o%j --error=${outdir}/${totalSitesJobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   cat ${outdir}/chr*_numSites.txt \
      | awk 'BEGIN{sum=0}{sum += \$0}END{print sum}' \
      > $tempfile
   rm ${outdir}/chr*_numSites.txt
EOF
	       )
dependencyString2="--dependency=afterok:"${totalSitesJobID}

# ----------------------------------
# echo -e "Executing \"part 2a\"..."
# ----------------------------------

dependencies="afterok" # SLURM syntax
linenum=0
# Important:
# Need to code this as follows (i.e., "while read line do ... done <<< input"
# instead of "input | while read line do ... done"),
# because otherwise the body of the loop gets executed in a subshell (due to the '|'),
# and the updates to variable "dependencies" aren't passed back to the parent shell
# (i.e., they're lost after the "done" statement).
while read line
do
    ((linenum++))
    file=$line
    if [ ! -s $file ]; then
	echo -e "Error:  File \"$file\", on line $linenum of $fileOfFilenames, was not found, or it is empty."
	exit 2
    fi
    chr=`head -n 1 $file | cut -f1`
    firstBegPos=`head -n 1 $file | cut -f2`
    firstEndPos=`head -n 1 $file | cut -f3`
    regionWidth=`echo $firstBegPos | awk -v e=$firstEndPos '{print e - $1}'`
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
    memSize="5G"
    if [ ! -s $outfileObserved ]; then
	thisJobID=$(sbatch --parsable --partition=queue1 $dependencyString2 --job-name=$jobName --output=${outdir}/${jobName}.o%j --error=${outdir}/${jobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   totalNumSites=\`cat $tempfile\`
   $EXE2 $infile $measurementType \$totalNumSites $infileQ $outfileObserved $outfileQcat $chr $firstBegPos $regionWidth $infileQ2
   if [ \$? != 0 ]; then
      exit 2
   fi
EOF
		 )
    dependencies="${dependencies}:${thisJobID}"
    fi
    if [ "$group2spec" != "" ] && [ ! -s $outfileNulls ]; then
	infile=${outdir}/${chr}${randFilenameString}
	jobName="p2r_$chr"
	thisJobID=$(sbatch --parsable --partition=queue1 $dependencyString2 --job-name=$jobName --output=${outdir}/${jobName}.o%j --error=${outdir}/${jobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   totalNumSites=\`cat $tempfile\`
   # The smaller number of arguments in the following call informs $EXE2 that it should only write the metric to $outfileNulls, with no additional info.
   $EXE2 $infile $measurementType \$totalNumSites $infileQ $infileQ2 $outfileNulls
   if [ \$? != 0 ]; then
      exit 2
   else # clean up
      rm -f $infile
   fi
EOF
		 )
    dependencies="${dependencies}:${thisJobID}"
    fi
	
done <<< "$(cat $fileOfFilenames)" # see note above regarding this syntax

# ----------------------------------
# echo -e "Executing \"part 2b\"..."
# ----------------------------------

memSize="10G" # is this a good amount???
catJobName=catNulls_${PID}
if [ $dependencies == "afterok" ]; then
    dependencyString=""
else
    dependencyString="--dependency="$dependencies
fi
dependencyString2=""

catJobID=$(sbatch --parsable --partition=queue1 $dependencyString --job-name=$catJobName --output=${outdir}/${catJobName}.o%j --error=${outdir}/${catJobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   module load htslib
   echo -n "Total number of sites is:  "
   cat $tempfile
# rm -f $tempfile
   cat ${outdir}/*_qcat.txt \
      | sort -k1b,1 -s \
      | awk 'BEGIN{FS = OFS="\t"}{print \$1,\$2,\$3,"id:"NR",qcat:"\$4}' \
      | bgzip \
      > ${outdir}/qcat.bed.gz
   rm -f ${outdir}/*_qcat.txt
          
   if [ "$group2spec" != "" ] && [ ! -s ${outdir}/allNullsGenomewide.txt ]; then
      cat ${outdir}/*_nulls.txt > ${outdir}/allNullsGenomewide.txt
      if [ \$? != 0 ]; then
         exit 2
      fi
      rm -f ${outdir}/*_nulls.txt
   fi
   # Clean up some files
   rm -f ${outdir}/*${PfilenameString}
EOF
	)
dependencyString2="--dependency=afterok:"${catJobID}

# ----------------------------------
# echo -e "Executing \"part 3\"..."
# ----------------------------------

dependencies="afterok"

while read line
do
    origFile=$line
    chr=`head -n 1 $origFile | cut -f1`
    begPos=`head -n 1 $origFile | cut -f2`
    endPos=`head -n 1 $origFile | cut -f3`
    infile=${chr}_observed.txt
    if [ -s ${outdir}/allNullsGenomewide.txt ]; then
	outfile=`echo $infile | sed 's/_observed.txt$/_withPvals.bed/g'`
    else
	outfile=`echo $infile | sed 's/txt$/bed/g'`
    fi

    outfile=${outdir}/${outfile}
    jobName="p3_$chr"
    memSize="10G" # is this a good amount???
    if [ ! -s $outfile ]; then
	thisJobID=$(sbatch --parsable --partition=queue1 $dependencyString2 --job-name=$jobName --output=${outdir}/${jobName}.o%j --error=${outdir}/${jobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   if [ -s ${outdir}/allNullsGenomewide.txt ]; then
      outfile=\`echo $infile | sed 's/_observed.txt$/_withPvals.bed/g'\`
      infile=${outdir}/$infile
      outfile=${outdir}/\$outfile
      $EXE3 \$infile ${outdir}/allNullsGenomewide.txt \$outfile
      if [ \$? == 0 ]; then
         rm -f \$infile
      fi
   else
      outfile=\`echo $infile | sed 's/txt$/bed/g'\`
      infile=${outdir}/$infile
      outfile=${outdir}/\$outfile
      mv \$infile \$outfile
   fi
EOF
	 )
    dependencies="${dependencies}:${thisJobID}"
    fi
done <<< "$(cat $fileOfFilenames)"

# --------------------------------------
# echo -e "Collating and cleaning up..."
# --------------------------------------

dependencyString="--dependency=$dependencies"
finalJobName=finalStep
memSize="5G" # is this a good number?

finalJobID=$(sbatch --parsable --partition=queue1 $dependencyString --job-name=$finalJobName --output=${outdir}/${finalJobName}.o%j --error=${outdir}/${finalJobName}.e%j --mem=$memSize <<EOF
#! /bin/bash
   module load bedops
   if [ -s ${outdir}/allNullsGenomewide.txt ] && [ \`ls -1 ${outdir}/*_withPvals.bed | wc -l\` != "0" ]; then
      rm -f ${outdir}/allNullsGenomewide.txt
   fi
   bedops -u ${outdir}/*.bed | starch - > ${outdir}/observations.starch
   rm -f ${outdir}/*.bed
EOF
	   )

exit 0
