origin_bam=/projects/korneliussen/data/vikings2020/bam
mkdir -p $PWD/faroe/bam
faroe=$PWD/faroe/bam
# REF=/maps/projects/skywalker/scratch/abi/ref/hs37d5.fa
REF=/projects/korneliussen/data/hs37d5/hs37d5.fa
helper_files=$PWD/helper_files


# # for i in $(cat faroe/samples.txt)
# # do
# #   cp $origin_bam/$i.final.bam $faroe/$i.bam
# # done

# ## index bam files
# for bam in $(ls $faroe/*.bam)
# do
#   samtools index $bam
# done


basedir=$PWD/faroe
# blist=$basedir/bamlist.txt
# ls $basedir/bam/*.bam > $blist
# mkdir -p $basedir/results/mt/vcf
# mt_vcf=$basedir/results/mt/vcf
# mkdir -p $basedir/results/mt/haplo
# mt_haplo=$basedir/results/mt/haplo

# # mkdir -p $basedir/results/yhaplo
# # y=$basedir/results/yhaplo

# # # ####  MT chromosome - START ####
# # for bam in $(ls $faroe/*.bam)
# # do
# #   bcftools mpileup  $bam --ignore-RG -Ou --fasta-ref $REF -I -r MT --threads 15| bcftools call --threads 15 -m -Oz -o $mt_vcf/$(basename $bam .bam).vcf.gz
# #   bcftools index $mt_vcf/$(basename $bam .bam).vcf.gz
# #   java -jar /projects/korneliussen/apps/haplogrep/haplogrep-2.1.25.jar --format vcf --in $mt_vcf/$(basename $bam .bam).vcf.gz --out $mt_haplo/$(basename $bam).haplo
# # done
# # # ####  MT chromosome - END ####



# mkdir -p $basedir/results/y/vcf
# y_vcf=$basedir/results/y/vcf
# mkdir -p $basedir/results/y/haplo
# y_haplo=$basedir/results/y/haplo

# ####  Y chromosome - START ####
# for bam in $(ls $faroe/*.bam)
# do
#   bcftools mpileup $bam --ignore-RG -f $REF  -Q 13 -q 30 -r Y  | bcftools call -m -Oz -o $y_vcf/$(basename $bam).chrY.vcf.gz
#   yhaplo -i $y_vcf/$(basename $bam).chrY.vcf.gz -o $y_haplo/$(basename $bam).haplo
# done
# ####  Y chromosome - END ####

###### extract chromosome 22 and do variant calling ####


# ####  select only chr22 
# mkdir -p $PWD/faroe/bam22
# bam22=$PWD/faroe/bam22
# for bam in $(ls $faroe/*.bam)
# do
#   samtools view -h $bam 22 > $bam22/$(basename $bam .bam).chr22.bam
#   samtools sort $bam22/$(basename $bam .bam).chr22.bam > $bam22/$(basename $bam .bam).chr22.sort.bam
#   samtools index $bam22/$(basename $bam .bam).chr22.sort.bam
# done

# blist22=$basedir/bamlist22.txt
# ls $basedir/bam22/*.sort.bam > $blist22
# # do variant calling on chr22 of samples
# bcftools mpileup -f $REF --bam-list  $blist22 -r 22 --threads 15 | bcftools call -m -Oz --threads 15 - > $basedir/allbams_final.vcf.gz

# plink --vcf $basedir/allbams_final.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_final

# ## KING- START ##
mkdir -p $basedir/results/KING/
dir=$basedir/results/KING
cd $dir
## would make more sense if it is allbams_clean.bed, since it's the one used in lcmlkin
# king -b $basedir/allbams_final.bed --kinship 
king -b $basedir/allbams_clean.bed --kinship 
cd $basedir
# ## KING- END ##


## PLINK - START ##
mkdir -p $basedir/results/plink/
dir=$basedir/results/plink
cd $dir
tmp=$dir/tmp
plink --bfile $basedir/allbams_clean --genome --out plink.clean 
sed 's/[[:space:]]\+/\t/g'  $dir/plink.clean.genome  | cut -f2- > $dir/tmp && mv $dir/tmp $dir/plink.clean.genome
cd $basedir
## PLINK - END ##

# ## READ - START ##
# ./read.sh
# ## READ - END ##


# ## LMCLKIN - START ###
# ./lcmlkin_realData.sh
# ## LMCLKIN - END ###


# ## RELATE - START ##
# ./relate_realData.sh
# ## RELATE - END ##

# ## CORRECTKIN - START ##
# ./correctKin_realData.sh
# ## CORRECTKIN - END ##



# ## FIX ##
# for i in $(seq 0 16); do echo $i; done > $basedir/numbers.txt
# paste -d '\t' $basedir/numbers.txt $basedir/samples.txt | awk '{print $1 "\t" $2}' > $basedir/my_dict.txt

# # ### NGSRELATE FIX - START ####
# dir=$basedir/results/relate
# tmpfile=$basedir/results/relate/tmp

# for infile in $dir/data_AFR_AF $dir/data_EUR_AF $dir/data_AMR_AF $dir/data_EAS_AF $dir/data_SAS_AF $dir/data_angsd
# do

#     outfile="$infile"_out


#     # substitute values in column 2 of the input file
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' $basedir/my_dict.txt "$infile" > "$tmpfile"
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' $basedir/my_dict.txt "$tmpfile" > "$outfile"

#     rm "$tmpfile"
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

# paste $outfile correction_ngsrelate_twin.txt  > $tmp && mv $tmp $outfile

# # add column with program name
# awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="plink"}1' $outfile | sed 's/[[:space:]]\+/\t/g' | cut -f2-  > $tmp && mv $tmp $outfile

# ## PLINK FIX - END ####
