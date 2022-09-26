#!/bin/bash

numjobs=8
numneighbors=101

while getopts t:i:o:s:j:n: flag
do
    case "${flag}" in
        t) title=${OPTARG};;
        i) src=${OPTARG};;
        o) out=${OPTARG};;
        s) scores=${OPTARG};;
        j) numjobs=${OPTARG};;
        n) numneighbors=${OPTARG};;
    esac
done
echo "Job Title: $title";
echo "Src Directory: $src";
echo "Output Directory: $out";
echo "Scores File: $scores";
echo "Nearest Neighbor Jobs: $numjobs";
echo "Number Nearest Neighbors: $numneighbors";
echo ""


if [ ! -d "$out"/5k ]
then
    mkdir -p "$out"/5k/
fi

> "$out"/5k/similarity_search.err
> "$out"/5k/similarity_search.out

job5k=$(sbatch --job-name=Sim_Search_"$title"_5k.job --partition=queue0 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/5k/similarity_search.err --output=$out/5k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/5k/ -w 5 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job5k"


if [ ! -d "$out"/10k ]
then
    mkdir -p "$out"/10k/
fi

> "$out"/10k/similarity_search.err
> "$out"/10k/similarity_search.out

job10k=$(sbatch --job-name=Sim_Search_"$title"_10k.job --partition=queue0 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/10k/similarity_search.err --output=$out/10k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/10k/ -w 10 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job10k"


if [ ! -d "$out"/25k ]
then
    mkdir -p "$out"/25k/
fi

> "$out"/25k/similarity_search.err
> "$out"/25k/similarity_search.out

job25k=$(sbatch --job-name=Sim_Search_"$title"_25k.job --partition=queue0 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/25k/similarity_search.err --output=$out/25k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/25k/ -w 25 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job25k"


if [ ! -d "$out"/50k ]
then
    mkdir -p "$out"/50k/
fi

> "$out"/50k/similarity_search.err
> "$out"/50k/similarity_search.out

job50k=$(sbatch --job-name=Sim_Search_"$title"_50k.job --partition=queue1 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/50k/similarity_search.err --output=$out/50k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/50k/ -w 50 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job50k"


if [ ! -d "$out"/75k ]
then
    mkdir -p "$out"/75k/
fi

> "$out"/75k/similarity_search.err
> "$out"/75k/similarity_search.out

job75k=$(sbatch --job-name=Sim_Search_"$title"_75k.job --partition=queue1 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/75k/similarity_search.err --output=$out/75k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/75k/ -w 75 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job75k"


if [ ! -d "$out"/100k ]
then
    mkdir -p "$out"/100k/
fi

> "$out"/100k/similarity_search.err
> "$out"/100k/similarity_search.out

job100k=$(sbatch --job-name=Sim_Search_"$title"_100k.job --partition=queue1 -n 1 --nodes=1 --ntasks-per-node=6 --mem 16000 --error=$out/100k/similarity_search.err --output=$out/100k/similarity_search.out --wrap="python $src/similarity_search.py -s $scores -o $out/100k/ -w 100 -j $numjobs -n $numneighbors" | sed 's/Submitted batch job //')

echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job100k"

echo ""
echo "sacct --format=JobName%30,JobID,Elapsed,State,MaxRSS,ReqMem --jobs $job5k,$job10k,$job25k,$job50k,$job75k,$job100k"