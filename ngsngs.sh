#!/bin/bash -l
# # conda run -n my_env python my_script.py

# create reads and all pipeline from given folders and fasta files

# This script creates a 1 Simulation Pedigree 
## simulation name should be of type simx (x number of )

# cmd should be: screen -S -L source sim1 ./fullSim.sh sim1 no_deam 10
# example cmd:  source screen -S variantCall_childA -L -Logfile sim1/no_deam/cov10/screen_childA.log  ./fullSim.sh sim1 no_deam 10

# screen -S manyKinships -L -Logfile sim1/no_deam/cov10/many_kinships.log -dm bash -c "source ./fullSim2.sh sim1 no_deam 10"

original_path=$PWD

mkdir -p $PWD/$1/$2/cov$3
basedir=$PWD/$1/$2/cov$3 ## something like
helper_files_path=$PWD/helper_files


source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis

## SIMULATE NGSNGS READS - START ##
# # # ancient_fragment_length_distribution=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/Size_dist_sampling.txt
# # # read_quality=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/AccFreqL150R1.txt
mkdir -p $basedir/fq
fq_path=$basedir/fq

ancient_fragment_length_distribution=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/Size_dist_sampling.txt
# note: for no
read_quality_no_deam=$PWD/newqualprofile_R1.txt
read_quality_deam=$PWD/AccFreqL150R1.txt
adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG'
read_quality=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/AccFreqL150R1.txt

start=$(date +%s.%N)


if [ $2 == "no_deam" ]; then
    for fa in $(ls $PWD/final_fasta/*fa)
    do
        ngsngs -i $fa -c $3 -f fq.gz -ld norm,350,30 -a1 $adapter -p G -t 32 -t2 12 -seq SE -q1 $read_quality_no_deam -o $fq_path/$(basename $fa .fa)
        
        fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

        # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
    done
else [ $2 == "deam" ];
    for fa in $(ls $PWD/final_fasta/*fa)
    do
        ngsngs -i $fa -c $3 -f fq.gz -cl 100 -lf $ancient_fragment_length_distribution -a1 $adapter -t 32 -t2 12 -seq SE -q1 $read_quality_deam -m b7,0.024,0.36,0.68,0.0097 -o $fq_path/$(basename $fa .fa)

        fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

        # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
    done
fi


end=$(date +%s.%N)

echo "$(echo "$end - $start" | bc)" > "$basedir/ngsngs.time"