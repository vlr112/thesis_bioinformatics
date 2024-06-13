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
# mkdir -p $basedir/bam
# bam_path=$basedir/bam

# for f in $(ls $fq_path/*.gz|cut -f1 -d.|sort -u) 
# do
#     bwa mem fasta/new_reference_genome.fa <(cat $f.h1.final.trimm.fq.gz $f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$(basename $f .fq.gz).bam
# done

# for f in $(ls $bam_path/*bam)
# do
#     samtools index $f
# done |parallel
# ## MAP READS - END ##


# # ## VARIANT CALLING -START ##

# ref=$PWD/fasta/new_reference_genome.fa
# mkdir -p $basedir/results/
# dir=$basedir/results/
# blist=$basedir/bamlist.txt
# ls $basedir/bam/*.bam > $blist

# pos=$PWD/helper_files/pos

# bcftools mpileup -f  $ref --bam-list  $blist --regions-file $pos --threads 15 | bcftools call -m -Oz --threads 15 - > $basedir/allbams_final.vcf.gz


# # ## get plink binary files from bcftools vcf ##
# plink --vcf $basedir/allbams_final.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_final

# # ## need to substitute the family column in .fam file with familyX
# awk '{ $1 = "familyX" }1' $basedir/allbams_final.fam > $basedir/tmp && mv $basedir/tmp $basedir/allbams_final.fam
# awk '{split($2,a,"/"); sub(/\.bam$/, "", a[length(a)]); $2=a[length(a)]; print}' $basedir/allbams_final.fam | column -t -s '' > $basedir/tmp2 && mv $basedir/tmp2 $basedir/allbams_final.fam
# # # ## VARIANT CALLING -END ##


# ## LCMLKIN - START ##
# # Includes special type of variant Calling for lcmlkin
# ./lcmlkin.sh $basedir
# cd $original_path
# ## LCMLKIN - END ##

# # for curiosity sake, calculate accuracy of genotype calls for bcftools and lcmlkin vcf.
# # # screen -S geno_acc_$1_$2_cov$3 -L -Logfile $basedir/many_kinships.log -dm bash -c "source ./genotype_acc_fast.sh $basedir $original_path"
# source ./genotype_acc.sh $basedir $original_path

# #### READ_lcmlkin - START #####
# # uses vcf from lcmlkin
# ./read.sh $basedir allbams_lcmlkin READ_lcmlkin
# cd $original_path
# #### READ_lcmlkin - END #####

# ### READ - START #####
# ## uses vcf from bcftools
# ./read.sh $basedir allbams_final READ
# cd $original_path
# ##### READ_lcmlkin - END #####

# ##### PLINK_lcmlkin - START #####
# ## uses vcf from lcmlkin
# ./plink.sh $basedir allbams_lcmlkin plink_lcmlkin $original_path 
# ##### PLINK_lcmlkin - START #####

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
# ./king.sh $basedir allbams_final_KING_bf KING $original_path
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

# ## ngsrelate alternative - start ##
# # ngsrelate -h $basedir/allbams_final.vcf.gz -O $basedir/results/relate/vcf.res
# ## ngsrelate alternative - end ##


# # ####### PREP FILES FOR ANALYSIS #######

# ## RUN ONCE - START ###
# for i in $(cat $basedir/bamlist.txt); do basename $i | cut -d'.' -f1 ; done > names.txt 
# for i in $(seq 0 38); do echo $i; done > numbers.txt
# # read the numbers and names files into arrays
# mapfile -t numbers < numbers.txt
# mapfile -t names < names.txt
# # use paste to join the lines from the two files
# # and use awk to print the key-value pairs
# paste -d '\t' numbers.txt names.txt | awk '{print $1 "\t" $2}' > my_dict.txt
# ## RUN ONCE - END ###


# ### NGSRELATE FIX - START ####
# # infile=$basedir/results/relate/data
# infile1=$basedir/results/relate/real_data
# infile2=$basedir/results/relate/angsd_data
# infile3=$basedir/results/relate/randomFreq_data
# tmpfile=$basedir/results/relate/tmp
# # outfile=$basedir/results/relate/data_out
# # outfile=$basedir/results/relate/ouf
# tmp=$basedir/results/relate/tmp

# for infile in "$infile3" "$infile1" "$infile2"
# do

#     outfile="$infile"_out


#     # substitute values in column 2 of the input file
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' my_dict.txt "$infile" > "$tmpfile"
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' my_dict.txt "$tmpfile" > "$outfile"

#     rm "$tmpfile"

#     awk -v col1="$3" 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' "$outfile" > "$tmp" && mv "$tmp" "$outfile"
#     awk -v col1="$2" 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' "$outfile" > "$tmp" && mv "$tmp" "$outfile"

#     paste "$outfile" correction_ngsrelate_twin_fix.txt > "$tmp" && mv "$tmp" "$outfile"
#     if [[ "$infile" == "$infile1" ]]
#     then
#         # add column with program name
#         # # # # awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="ngsrelate"}1' $outfile | sed 's/[[:space:]]\+/\t/g' > $tmp && mv $tmp $outfile
#         awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="ngsrelate_trueFreq"}1' "$outfile" | sed 's/[[:space:]]\+/\t/g' > "$tmp" && mv "$tmp" "$outfile"
#     elif [[ "$infile" == "$infile2" ]]
#     then
#         awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="ngsrelate_angsdFreq"}1' "$outfile" | sed 's/[[:space:]]\+/\t/g' > "$tmp" && mv "$tmp" "$outfile"
#     else
#         awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="ngsrelate_randomFreq"}1' "$outfile" | sed 's/[[:space:]]\+/\t/g' > "$tmp" && mv "$tmp" "$outfile"
#     fi
# done

# ## NGSRELATE FIX - END ####


# ### PLINK FIX - START ####
# infile=$basedir/results/plink/plink.genome
# tmpfile=$basedir/results/plink/tmp
# outfile=$basedir/results/plink/plink_out.genome
# tmp=$basedir/results/plink/tmp

# # infile=$basedir/results/plink_lcmlkin/plink.genome
# # tmpfile=$basedir/results/plink_lcmlkin/tmp
# # outfile=$basedir/results/plink_lcmlkin/plink_lcmlkin_out.genome
# # tmp=$basedir/results/plink_lcmlkin/tmp

# # substitute values in column 2 of the input file
# awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' my_dict.txt $infile > $tmpfile
# awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' my_dict.txt $tmpfile > $outfile

# rm $tmpfile

# awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $outfile > $tmp && mv $tmp $outfile
# awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $outfile > $tmp && mv $tmp $outfile

# paste $outfile correction_ngsrelate_twin_fix.txt  > $tmp && mv $tmp $outfile

# # add column with program name
# awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="plink"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f2-  > $tmp && mv $tmp $outfile
# # awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="plink_Lcmlkin"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f2-  > $tmp && mv $tmp $outfile

# ## PLINK FIX - END ####


# ### LCMLKIN FIX - START ####
# # infile=$basedir/results/lcmlkin/lcmlkin_results
# infile=$basedir/results/lcmlkin/lcmlkin_results_update
# tmpfile=$basedir/results/lcmlkin/tmpfile
# # outfile=$basedir/results/lcmlkin/lcmlkin_results_out
# outfile=$basedir/results/lcmlkin/lcmlkin_results_update_out
# tmp=$basedir/results/lcmlkin/tmp

# # rm $outfile

# awk '$2 == "cov0" { next } { print }' $infile >  $outfile
# awk '$1 == "cov0" { next } { print }' $outfile >  $tmp



# # order rows in lcmlkin according to relate/data_out order of first 2 columns (samples id columns)
# awk 'FNR == NR {lineno[$1,$2] = NR; next} {print lineno[$1,$2], $0;}' $basedir/results/relate/real_data_out $tmp | tail -n +2 | cut -f1- | sort -n -k1,1  > $tmpfile
# sed  '1i sort Ind1 Ind2 k0_hat k1_hat k2_hat pi_HAT nbSNP' $tmpfile | sed 's/[[:space:]]\+/\t/g' > $outfile

# # awk 'FNR == NR {lineno[$1,$2] = NR; next} {print lineno[$1,$2], $0;}' $PWD/sim1/deam/cov1/results/relate/real_data_out $tmp > $tmpfile

# awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $outfile > $tmp && mv $tmp $outfile
# awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $outfile > $tmp && mv $tmp $outfile

# paste $outfile correction_ngsrelate_twin_fix.txt > $tmp && mv $tmp $outfile
# rm $tmpfile
# ## add column with program name
# awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="lcmlkin"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f2-  > $tmp && mv $tmp $outfile
# ### LCMLKIN FIX - END ####




# ### KING FIX - START ####
# # infile=$basedir/results/lcmlkin/lcmlkin_results
# # infile=$basedir/results/KING/king.kin0
# infile=$basedir/results/KING/king.kin

# tmpfile=$basedir/results/KING/tmp
# # outfile=$basedir/results/lcmlkin/lcmlkin_results_out
# # outfile=$basedir/results/KING/king_out.kin0
# outfile=$basedir/results/KING/king_out.kin

# tmp=$basedir/results/KING/tmp

# # rm $outfile
# # order rows in lcmlkin according to relate/data_out order of first 2 columns (samples id columns)
# # (awk 'NR==1 {header=$0; print; next} $2=="childA" && $3=="childAA" {print; next} {print | "sort -k1,1n"}' $infile) > $outfile
# awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $infile > $tmp && mv $tmp $outfile
# awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $outfile > $tmp && mv $tmp $outfile

# # awk 'BEGIN{ FS = OFS = "\t" } { print (NR==1? "cov" : col1), $0 }' $outfile > $tmp && mv $tmp $outfile

# paste $outfile correction_ngsrelate_twin_fix.txt > $tmp && mv $tmp $outfile

# ## add column with program name
# # awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="king0"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f1-  > $tmp && mv $tmp $outfile
# awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="king"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f1-  > $tmp && mv $tmp $outfile

# # ##remove rows that were not outputed by KING program
# awk 'NF>5' $outfile > $tmp && mv $tmp $outfile
# ### KING FIX - END ####


# ### CORRECTKIN FIX - START ####
# infile=$basedir/results/correctKin/angsdHaploCall.relatives.tsv
# outfile=$basedir/results/correctKin/out_angsdHaploCall.relatives.tsv
# tmpfile=$basedir/results/correctKin/tmp
# tmpfile2=$basedir/results/correctKin/tmp2

# awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $infile > $tmpfile 
# awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $tmpfile > $tmpfile2
# python helper_files/fix_correctKin.py $basedir/results/correctKin tmp2 

# awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="correctKin"}1' $outfile | sed 's/[[:space:]]\+/\t/g'  > $tmp && mv $tmp $outfile
# rm $tmp $tmp2
# ### CORRECTKIN FIX - END ####



## READ -START ##

infile=$basedir/results/READ/Read_intermediate_output
outfile=$basedir/results/READ/Read_intermediate_output_out
tmpfile=$basedir/results/READ/tmp
tmpfile2=$basedir/results/READ/tmp2

sed -i '1s/^/PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow\tIBS2\tIBS0\tP1\tP0\tMissing\tTMP\n/' $infile


awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $infile > $tmpfile 
awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $tmpfile > $tmpfile2

awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="READ"}1' $tmpfile2 | sed 's/[[:space:]]\+/\t/g'  > $tmpfile && mv $tmpfile $outfile

# sed -i '1s/^/PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow\tIBS2\tIBS0\tP1\tP0\tMissing\tTMP\n/' $outfile


rm $tmpfile $tmpfile2
## READ -END ##



# # ## after previous corrections, create another file where we calculate RMSE
# # ## according to true IBD proportions (from fasta files)

# # ### NGSRELATE - TRUE IBD ###
# # infile=real_data_out
# # inpath=$basedir/results/relate/
# # outfile=$basedir/results/relate/rmse_true

# # python $helper_files_path/rmse_true.py $inpath $infile $outfile

# # ### PLINK - TRUE IBD ###
# # infile=plink_out.genome
# # inpath=$basedir/results/plink/
# # outfile=$basedir/results/plink/rmse_true

# # python $helper_files_path/rmse_true.py $inpath $infile $outfile

# # ### PLINK - TRUE IBD ###
# # infile=lcmlkin_results_out
# # inpath=$basedir/results/lcmlkin/
# # outfile=$basedir/results/lcmlkin/rmse_true

# # python $helper_files_path/rmse_true.py $inpath $infile $outfile