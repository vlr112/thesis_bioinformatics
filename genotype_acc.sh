# get path to fasta files with haplotypes

# screen -S manyKinships -L -Logfile sim1/no_deam/cov10/many_kinships.log -dm bash -c "source ./genotype_acc.sh $PWD/sim1/no_deam/cov10 $PWD"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis



# mkdir -p $PWD/$1/$2/cov$3
# basedir=$PWD/$1/$2/cov$3 ## something like
basedir=$1
original_path=$2
fasta_files_path=$2/final_fasta
fastalist=$basedir/fastalist.txt
ls $fasta_files_path/*.fa > $fastalist

mkdir -p $basedir/accuracy
acc=$basedir/accuracy

# #### ACCURCY INDIVIDUAL ANALYSIS - START ####

# vcf=$basedir/allbams_lcmlkin.vcf
# vcf_no_header=$acc/allbams_lcmlkin.single.no_header.vcf
# grep -v "^#" $vcf > $vcf_no_header
# python $original_path/helper_files/genotype_acc.py $vcf_no_header $fastalist lcmlkin $acc/ > $acc/genotype_accuracy_lcmlkin.single.txt

# # for bcftools vcf file:
# bcftools view -H $basedir/allbams_final.vcf.gz > $acc/allbams_final.single.no_header.vcf

# python $original_path/helper_files/genotype_acc.py $acc/allbams_final.no_header.vcf $fastalist bcftools $acc/ > $acc/genotype_accuracy_bcftools.single.txt
# #### ACCURCY INDIVIDUAL ANALYSIS - END ####


#### ACCURCY JOINT ANALYSIS - START ####
## in both vcf, keep only shares SNP positions
gunzip -c $basedir/allbams_final.vcf.gz | bgzip -c > $acc/allbams_final.vcf.bgz
bgzip -c $basedir/allbams_lcmlkin.vcf > $acc/allbams_lcmlkin.vcf.bgz

bcftools index -f $acc/allbams_final.vcf.bgz > $acc/allbams_final.vcf.bgz.csi
bcftools index -f $acc/allbams_lcmlkin.vcf.bgz > $acc/allbams_lcmlkin.vcf.bgz.csi

bcftools isec $acc/allbams_final.vcf.bgz $acc/allbams_lcmlkin.vcf.bgz -p $acc

## remove headers now
# for bcftools vcf file:
bcftools view -H $acc/0002.vcf > $acc/allbams_final.common.no_header.vcf
bcftools view -H $acc/0003.vcf > $acc/allbams_lcmlkin.common.no_header.vcf

vcf_bcftools=$$acc/allbams_final.common.no_header.vcf
vcf_lcmlkin=$acc/allbams_lcmlkin.common.no_header.vcf

python $original_path/helper_files/genotype_acc.changed.py $vcf_bcftools $vcf_lcmlkin $fastalist $acc/> $acc/joint_accuracy.txt

#### ACCURCY JOINT ANALYSIS - END ####









############################################# trash


# # # # # get lists of POS for each vcf
# # # grep "^[^#]" $vcf | awk '{print $2}' > $basedir/lcmlkin_list
# # # zgrep "^[^#]" $basedir/allbams_final.vcf.gz  | awk '{print $2}' > $basedir/allbams_final_list 

# # # # get commom POS 
# # # comm -12  <(sort $basedir/lcmlkin_list) <(sort $basedir/allbams_final_list) | wc -l

# get percentage of missing genotypes:
# grep -v "^#" allbams_final.vcf.gz  |cut -f 10- | tr "\t" "\n" | cut -d ':' -f 1 | awk '/^\.\/\./ {NC++;} END{printf("%f\n",NC/(1.0*NR))}' 

## comparison test ##

# python $original_path/helper_files/genotype_acc.changed.py $vcf $basedir/allbams_final.vcf.gz  