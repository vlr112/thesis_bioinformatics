# get list with CEU individuals

basedir=$PWD/faroe
dir=$basedir/trash

cd $dir

# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# grep CEU $dir/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > $dir/CEU.samples.list


# vcf-subset -c $dir/CEU.samples.list $basedir/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
#  fill-an-ac | \
#  bgzip -c > $dir/CEU.chr22.phase3.vcf.gz


## restrict the chr22 vcf file to be biallelic, the ancient 1240k panel positions, only SNPs ... 
# vcftools --gzvcf $dir/CEU.chr22.phase3.vcf.gz --positions $basedir/v42.4.1240K_autosomal.sites --remove-indels --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | bgzip -c > $dir/CEU.chr22.biallelic.vcf.gz
# bcftools index -f $dir/CEU.chr22.biallelic.vcf.gz

# bcftools +fill-tags $dir/CEU.chr22.biallelic.vcf.gz  -- -t AF > $dir/CEU.AF.chr22.biallelic.vcf.gz 

# bcftools query -f '%POS %INFO/AF\n' $dir/CEU.AF.chr22.biallelic.vcf.gz | sed 's/[[:space:]]\+/\t/g' > $dir/CEU.AF.chr22.biallelic.freq

awk 'FNR==NR{positions[$1]; next} $1 in positions' $basedir/results/relate/freq  $dir/CEU.AF.chr22.biallelic.freq | cut -f2 | sed 1d > $dir/clean.CEU.AF.chr22.biallelic.freq

