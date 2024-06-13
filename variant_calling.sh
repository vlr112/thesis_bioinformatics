
mkdir -p $PWD/$1/$2/cov$3
basedir=$PWD/$1/$2/cov$3 ## something like
# basedir=$1
helper_files_path=$PWD/helper_files


# ## VARIANT CALLING -START ##

ref=$PWD/fasta/new_reference_genome.fa
mkdir -p $basedir/results/
dir=$basedir/results/
blist=$basedir/bamlist.txt
ls $basedir/bam/*.bam > $blist

pos=$PWD/helper_files/pos

start=$(date +%s.%N)

bcftools mpileup -f  $ref --bam-list  $blist --regions-file $pos --threads 10 | bcftools call -m -Oz --threads 10 - > $basedir/allbams_final.vcf.gz


# ## get plink binary files from bcftools vcf ##
plink --vcf $basedir/allbams_final.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_final

# ## need to substitute the family column in .fam file with familyX
awk '{ $1 = "familyX" }1' $basedir/allbams_final.fam > $basedir/tmp && mv $basedir/tmp $basedir/allbams_final.fam
awk '{split($2,a,"/"); sub(/\.bam$/, "", a[length(a)]); $2=a[length(a)]; print}' $basedir/allbams_final.fam | column -t -s '' > $basedir/tmp2 && mv $basedir/tmp2 $basedir/allbams_final.fam
# # ## VARIANT CALLING -END ##

end=$(date +%s.%N)
echo "$(echo "$end - $start" | bc)" > "$basedir/vcf.time"


