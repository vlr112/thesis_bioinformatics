


basedir=$PWD/faroe
blist22=$basedir/bamlist22.txt
ls $basedir/bam22/*.sort.bam > $blist22

# mkdir -p $basedir/results/relate
dir0=$basedir/results/relate

# # ## fetch and index the positions of interest that exist on the 1000g vcf and 1240K ancient panel
# # bcftools query -f '%CHROM %POS\n' chr22.biallelic.vcf.gz | sed 's/[[:space:]]\+/\t/g' > $basedir/positions
# # angsd sites index $basedir/positions

# # angsd -checkBamHeaders 0 -b $blist22 -gl 1 -doMajorMinor 1 -doMaf 1 -minMapQ 30 -minQ 20 -minMaf 0.05 -doGlf 3 -nThreads 10 -sites $basedir/positions -out $dir/data 
# # zcat $dir/data.mafs.gz | cut -f2,5 > $dir/freq

# # awk 'FNR==NR{positions[$1]; next} $1 in positions' $basedir/chr22_freq $dir/freq | cut -f2 | sed 1d >  $dir/real_freq_angsd
# # ngsrelate  -g $dir/data.glf.gz -n 17 -f $dir/real_freq_angsd  -O $dir/data_angsd



# # ### get frequencies for supepopulations
# mkdir -p $basedir/superPop

# for  AF in 'EUR_AF' 'EAS_AF' 'AMR_AF' 'AFR_AF' 'SAS_AF'
# # for  AF in 'EUR_AF' 
# do
#     # get allele frequencies of each super population
#     bcftools query -f '%POS %INFO/'"$AF"'\n' $basedir/chr22.biallelic.vcf.gz | sed 's/[[:space:]]\+/\t/g' > $basedir/superPop/chr22_"$AF"_freq
#     # filter for only positions that also exist from angsd command 
#     awk 'FNR==NR{positions[$1]; next} $1 in positions' $dir/freq  $basedir/superPop/chr22_"$AF"_freq | cut -f2 | sed 1d > $basedir/superPop/real_freq_$AF
#     ngsrelate  -g $dir/data.glf.gz -l 0.1 -n 17 -f $basedir/superPop/real_freq_$AF  -O $dir/data_"$AF"_tmp

#     tmpfile=$basedir/results/relate/tmp

#     # substitute the numbers with actual samples names in column 2 of the input file
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' $basedir/my_dict.txt $dir/data_$AF_tmp > $tmpfile
#     # awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' $basedir/my_dict.txt "$tmpfile" > $dir/data_"$AF"_tmp_out

#     rm "$tmpfile"

#     # awk 'FNR==NR{positions[$1]; next} $1 in positions' $dir/freq  $basedir/superPop/chr22_"$AF"_freq > $basedir/superPop/real_freq_tmp_$AF

# done
 

mkdir -p $basedir/results/relate2
dir=$basedir/results/relate2

# ## fetch and index the positions of interest that exist on the 1000g vcf and 1240K ancient panel
# bcftools query -f '%CHROM %POS\n' $basedir/chr22.biallelic.vcf.gz | sed 's/[[:space:]]\+/\t/g' > $basedir/positions
# angsd sites index $basedir/positions

# angsd -checkBamHeaders 0 -b $blist22 -gl 1 -doMajorMinor 1 -doMaf 1 -minMapQ 30 -minQ 20 -minMaf 0.05 -doGlf 3 -nThreads 10 -sites $basedir/positions -out $dir/data 
# zcat $dir/data.mafs.gz | cut -f2,5 > $dir/freq

# awk 'FNR==NR{positions[$1]; next} $1 in positions' $basedir/chr22_freq $dir/freq | cut -f2 | sed 1d >  $dir/real_freq_angsd
# ngsrelate  -g $dir/data.glf.gz -n 17 -f $dir/real_freq_angsd  -O $dir/data_angsd

tmpfile=$basedir/results/relate2/tmp

awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' $basedir/my_dict.txt $dir/data_angsd > $tmpfile
awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' $basedir/my_dict.txt $tmpfile > $dir/data_angsd_out


# # ### get frequencies for supepopulations
# mkdir -p $basedir/superPop2

# for  AF in 'EUR_AF' 'EAS_AF' 'AMR_AF' 'AFR_AF' 'SAS_AF'
# # for  AF in 'EUR_AF' 
# do
#     # get allele frequencies of each super population
#     # bcftools query -f '%POS %INFO/'"$AF"'\n' $basedir/chr22.biallelic.vcf.gz | sed 's/[[:space:]]\+/\t/g' > $basedir/superPop2/chr22_"$AF"_freq
#     # # filter for only positions that also exist from angsd command 
#     # awk 'FNR==NR{positions[$1]; next} $1 in positions' $dir/freq  $basedir/superPop2/chr22_"$AF"_freq | cut -f2 > $basedir/superPop2/real_freq_$AF
#     # ngsrelate  -g $dir/data.glf.gz -n 17 -f $basedir/superPop2/real_freq_$AF  -O $dir/data_"$AF"

#     tmpfile=$basedir/results/relate2/tmp

#     # substitute the numbers with actual samples names in column 2 of the input file
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next}{if ($2 in a) $2=a[$2]}1' $basedir/my_dict.txt $dir/data_$AF > $tmpfile
#     awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]=$2;next} FNR==1{print;next} $1 in a{$1=a[$1]}1' $basedir/my_dict.txt $tmpfile > $dir/data_"$AF"_out

#     # rm $tmpfile

# done
 


