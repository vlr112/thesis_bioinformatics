##### my code ###

basedir=$PWD/faroe

REF=/projects/korneliussen/data/hs37d5/hs37d5.fa
blist22=$basedir/bamlist22.txt
mkdir -p $basedir/results/lcmlkin/
dir=$basedir/results/lcmlkin


# ## get chr22 vcf file
# wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi}

# ## get snplist with CHRM and POS from 1240k panel (from ancient individuals)
# wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V42/V42.4/SHARE/public.dir/v42.4.1240K.snp

# # for kinship analysis restrict only to autosomal sites and chr22 only
# sed -r 's/^[ ]+//;s/[ ]+/\t/g' v42.4.1240K.snp | cut -f 2,4 | awk '{if($1>=1 && $1 =22){print}}'| sed 's/[[:space:]]\+/\t/g' > v42.4.1240K_autosomal.sites


# ## restrict the chr22 vcf file to be biallelic, the ancient 1240k panel positions, only SNPs ... 
# #doesn't work
# # bcftools view -c 1 -c 1:nonmajor --min-alleles 2 --max-alleles 2 --exclude-types indels --regions-file $basedir/v42.4.1240K_autosomal.sites --threads 20 $basedir/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > $basedir/chr22.biallelic.vcf.gz  
# vcftools --gzvcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --positions v42.4.1240K_autosomal.sites --remove-indels --min-alleles 2 --max-alleles 2 --recode-INFO-all --recode --stdout | bgzip -c > chr22.biallelic.vcf.gz 
# bcftools index -f $basedir/chr22.biallelic.vcf.gz  


# bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $basedir/chr22.biallelic.vcf.gz | bgzip -c > $basedir/chr22.biallelic.tsv.gz 
# tabix -s1 -b2 -e2 $basedir/chr22.biallelic.tsv.gz 

# # do variant calling on 17 individuals with restrictions 

# bcftools mpileup -f $REF -B -q20 -Q5 -I -a 'FORMAT/DP' -T $basedir/chr22.biallelic.vcf.gz  -b $blist22 -Ou | bcftools call -Am -C alleles -f GP,GQ -T $basedir/chr22.biallelic.tsv.gz  -Oz -o $basedir/allbams_clean.vcf.gz 
# plink --vcf $basedir/allbams_clean.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_clean
# lcmlkin -i $basedir/allbams_clean.vcf -o $dir/lcmlkin_results -g all -l phred

lcmlkin -i $basedir/allbams_clean.vcf -u $basedir/samples_unrelated.txt -o $dir/lcmlkin_results_U -g all -l phred 

# bcftools mpileup -f $REF  -q20 -Q5 -T $basedir/chr22.biallelic.vcf.gz  -b $blist22 -Ou | bcftools call -m -T $basedir/chr22.biallelic.tsv.gz  -Oz -o $basedir/allbams_test.vcf.gz 
# plink --vcf $basedir/allbams_test.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_test
