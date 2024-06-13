
original_path=$PWD

# mkdir -p $PWD/cov$1/$2/
# basedir=$PWD/cov$1/$2 ## something like
# # helper_files_path=$PWD/helper_files


mkdir -p $PWD/$1/$2/cov$3
basedir=$PWD/$1/$2/cov$3 ## something like
mkdir -p $basedir/fq/
# fq_path=$basedir/fq

source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis
# rm -r $basedir/results 


# ## make final_fasta files. RUN ONCE!

# for i in $(ls $original_path/variable_sites/*txt)
#     do 
#         python make_invariable.py $i
#     done

# ## SIMULATE NGSNGS READS - START ##
# ancient_fragment_length_distribution=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/Size_dist_sampling.txt
# read_quality_no_deam=$PWD/newqualprofile_R1.txt
# read_quality_deam=$PWD/AccFreqL150R1.txt
# adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG'
# read_quality=/home/vlr112/thesis_vlr112/simulations/NGSNGS/Test_Examples/AccFreqL150R1.txt



# if [ $2 == "no_deam" ]; then
#     for fa in $(ls $PWD/final_fasta/*fa)
#     do
#         ngsngs -i $fa -c $3 -f fq.gz -ld norm,350,30 -a1 $adapter -p G -t 32 -t2 12 -seq SE -q1 $read_quality_no_deam -o $fq_path/$(basename $fa .fa)
        
#         fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

#         # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
#     done
# else [ $2 == "deam" ];
#     for fa in $(ls $PWD/final_fasta/*fa)
#     do
#         ngsngs -i $fa -c $3 -f fq.gz -cl 100 -lf $ancient_fragment_length_distribution -a1 $adapter -t 32 -t2 12 -seq SE -q1 $read_quality_deam -m b7,0.024,0.36,0.68,0.0097 -o $fq_path/$(basename $fa .fa)

#         fastp  -i $fq_path/$(basename $fa .fa).fq.gz -o $fq_path/$(basename $fa .fa).trimm.fq.gz

#         # fastqc  $fq_path/$(basename $fa .fa).trimm.fq.gz
#     done
# fi


# # ## MAP READS - START ##
# rm -rf $basedir/bam
mkdir -p $basedir/bam
bam_path=$basedir/bam
fq_path=$basedir/fq

for f in $(ls $fq_path/*.gz | cut -f 11 -d'/' | sort -u | cut -f1 -d'.'| uniq) 
# for f in $(ls $fq_path/*.gz|cut -f-2 -d.|sort -u)
do
    # echo $fq_path/$f
    bwa mem fasta/new_reference_genome.fa <(cat $fq_path/$f.h1.final.trimm.fq.gz $fq_path/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$f.bam
done

for f in $(ls $bam_path/*bam)
do
    samtools index $f
done |parallel

## VARIANT CALLING -START ##

ref=$PWD/fasta/new_reference_genome.fa
mkdir -p $basedir/results/
dir=$basedir/results/
blist=$basedir/bamlist.txt
ls $basedir/bam/*.bam > $blist
pos=$PWD/helper_files/pos

bcftools mpileup -f  $ref --bam-list  $blist --regions-file $pos --threads 10 | bcftools call -m -Oz --threads 10 - > $basedir/allbams_final.vcf.gz
# ## get plink binary files from bcftools vcf ##
plink --vcf $basedir/allbams_final.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_final

# ## need to substitute the family column in .fam file with familyX
awk '{ $1 = "familyX" }1' $basedir/allbams_final.fam > $basedir/tmp && mv $basedir/tmp $basedir/allbams_final.fam
awk '{split($2,a,"/"); sub(/\.bam$/, "", a[length(a)]); $2=a[length(a)]; print}' $basedir/allbams_final.fam | column -t -s '' > $basedir/tmp2 && mv $basedir/tmp2 $basedir/allbams_final.fam
## VARIANT CALLING -END ##


# ## LCMLKIN - START ##
# # Includes special type of variant Calling for lcmlkin
# ./lcmlkin.sh $basedir
# cd $original_path
# ## LCMLKIN - END ##

# # for curiosity sake, calculate accuracy of genotype calls for bcftools and lcmlkin vcf.
# # # screen -S geno_acc_$1_$2_cov$3 -L -Logfile $basedir/many_kinships.log -dm bash -c "source ./genotype_acc_fast.sh $basedir $original_path"
# # source ./genotype_acc.sh $basedir $original_path

# # #### READ_lcmlkin - START #####
# # # uses vcf from lcmlkin
# # ./read.sh $basedir allbams_lcmlkin READ_lcmlkin
# # cd $original_path
# # #### READ_lcmlkin - END #####

# # ### READ - START #####
# # ## uses vcf from bcftools
# # ./read.sh $basedir allbams_final READ
# # cd $original_path
# # ##### READ_lcmlkin - END #####

# # ##### PLINK_lcmlkin - START #####
# # ## uses vcf from lcmlkin
# # ./plink.sh $basedir allbams_lcmlkin plink_lcmlkin $original_path 
# # ##### PLINK_lcmlkin - START #####

# #### PLINK - START #####
# # uses vcf from bcftools
# ./plink.sh $basedir allbams_final plink $original_path 
# #### PLINK - START #####

# ##### KING_lcmlkin - START #####
# ## uses vcf from lcmlkin
# ./king.sh $basedir allbams_lcmlkin KING_lcmlkin $original_path
# ##### KING_lcmlkin - END #####

# ##### KING - START #####
# ## uses vcf from bcftools
# ./king.sh $basedir allbams_final KING $original_path
# ##### KING - END #####

# ##  CORRECTKIN, KIN AND NGSRELATE TAKE LONG TIME 
# ##### KIN - START #####
# ./kin.sh $basedir
# cd $original_path
# ##### KIN - END #####

# ##### correctKin - START #####
# # Includes special type of variant Calling for lcmlkin
# ./correctKin.sh $basedir
# cd $original_path
# ##### correctKin - END #####

# ##### NGSRELATE -START #####
# ./relate.sh $basedir
# cd $original_path
# ##### NGSRELATE - END #####


# ####### PREP FILES FOR ANALYSIS #######

# # ## RUN ONCE - START ###
# # for i in $(cat $basedir/bamlist.txt); do basename $i | cut -d'.' -f1 ; done > names.txt 
# # for i in $(seq 0 38); do echo $i; done > numbers.txt
# # # read the numbers and names files into arrays
# # mapfile -t numbers < numbers.txt
# # mapfile -t names < names.txt
# # # use paste to join the lines from the two files
# # # and use awk to print the key-value pairs
# # paste -d '\t' numbers.txt names.txt | awk '{print $1 "\t" $2}' > my_dict.txt
# # ## RUN ONCE - END ###

# infile=$basedir/results/relate/data
# tmpfile=$basedir/results/relate/tmp
# outfile=$basedir/results/relate/data_out
# tmp=$basedir/results/relate/tmp

# # substitute values in column 2 of the input file
# awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' my_dict.txt $infile > $tmpfile
# awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' my_dict.txt $tmpfile > $outfile

# rm $tmpfile

# awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $outfile > $tmp && mv $tmp $outfile
# awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $outfile > $tmp && mv $tmp $outfile

# paste $outfile correction_ngsrelate.txt > $tmp && mv $tmp $outfile


