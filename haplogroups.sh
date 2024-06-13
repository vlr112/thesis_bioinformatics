# # join results Y haplogroup

# find faroe/haplo_results/y/haplo/ -name "haplogroups.*"  -exec cat {} + | sed 's/[[:space:]]\+/\t/g' > tmp && \
# cut -f1 tmp | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp2  && \
# awk -v i=1  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp2 tmp | sed 's/[[:space:]]\+/\t/g' > faroe/haplo_results/y/results.txt

# # join results MT haplogroup

# find faroe/haplo_results/mt/haplo/ -name "*.bam.haplo"  > faroe/haplo_results/mt/filelist.txt

find faroe/haplo_results/mt/haplo/ -name "*.bam.haplo" -exec cat {} + |  sed -n '0~2p' > tmp  && \
cut -f1 tmp | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp2 && \
awk -v i=1  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp2 tmp | \
sed '1 i\SampleID      Range Haplogroup    Rank  Quality'| sed 's/[[:space:]]\+/\t/g' > faroe/haplo_results/mt/results.txt

rm tmp tmp2  
# header.haploMT