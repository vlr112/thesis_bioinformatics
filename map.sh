
mkdir -p $PWD/$1/$2/cov$3
basedir=$PWD/$1/$2/cov$3 ## something like
# basedir=$1
helper_files_path=$PWD/helper_files

mkdir -p $basedir/fq
fq_path=$basedir/fq

ancient_fragment_length_distribution=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/Size_dist_sampling.txt
# note: for no
read_quality_no_deam=$PWD/newqualprofile_R1.txt
read_quality_deam=$PWD/AccFreqL150R1.txt
adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG'
read_quality=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/AccFreqL150R1.txt

# start=$(date +%s.%N)


# if [ $2 == "no_deam" ]; then
#     for fa in $PWD/final_fasta/childAA.h1.final.fa $PWD/final_fasta/childAA.h2.final.fa
#     do
#         ngsngs -i $fa -c $3 -f fq.gz -ld norm,350,30 -a1 $adapter -p G -t 32 -t2 12 -seq SE -q1 $read_quality_no_deam -o $fq_path/$(basename $fa .fa)
        
#         fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

#         # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
#     done
# else [ $2 == "deam" ];
#     for fa in $PWD/final_fasta/childAA.h1.final.fa $PWD/final_fasta/childAA.h2.final.fa
#     do
#         ngsngs -i $fa -c $3 -f fq.gz -cl 100 -lf $ancient_fragment_length_distribution -a1 $adapter -t 32 -t2 12 -seq SE -q1 $read_quality_deam -m b7,0.024,0.36,0.68,0.0097 -o $fq_path/$(basename $fa .fa)

#         fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

#         # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
#     done
# fi


# ## MAP READS - START ##
mkdir -p $basedir/bam
bam_path=$basedir/bam

# for f in $(ls $fq_path/*.gz|cut -f1 -d.|sort -u) 
# do
# for f in "founder1" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder2" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9"

for f in "founder12"
do
    # echo $f
    # bwa mem fasta/new_reference_genome.fa <(cat $fq_path/$f.h1.final.trimm.fq.gz $fq_path/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$f.bam
    bwa mem fasta/new_reference_genome.fa <(cat $fq_path/$f.h1.final.trimm.fq.gz $fq_path/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$f.bam
done

for f in $(ls $bam_path/*bam)
do
    samtools index $f
done |parallel
## MAP READS - END ##

# "founder1" "founder20" "founder21" "founder22" "founder23" "founder24" "founder25" "founder26" "founder27" "founder28" "founder29" "founder2" "founder3" "founder4" "founder5" "founder6" "founder7" "founder8" "founder9"

# "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI"


# for f in  "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI" "founder0" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19"; do bwa mem $PWD/fasta/new_reference_genome.fa <(cat $PWD/fq/$f.h1.final.trimm.fq.gz $PWD/fq/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $PWD/bam/$f.bam ;done




# "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI" "founder0" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19"