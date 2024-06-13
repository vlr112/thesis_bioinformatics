
mkdir -p $PWD/$1/$2/cov$3
basedir=$PWD/$1/$2/cov$3 ## something like
# basedir=$1
helper_files_path=$PWD/helper_files

fq_path=$basedir/fq

# ## MAP READS - START ##
mkdir -p $basedir/bam
bam_path=$basedir/bam

# # for f in $(ls $fq_path/*.gz|cut -f1 -d.|sort -u) 
# # do
for f in "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI" "founder0" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19" 
do
    bwa mem fasta/new_reference_genome.fa <(cat $fq_path/$f.h1.final.trimm.fq.gz $fq_path/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $bam_path/$f.bam
done

for f in $(ls $bam_path/*bam)
do
    samtools index $f
done |parallel
## MAP READS - END ##


# "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI"
# for f in  "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI" "founder0" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19"; do bwa mem $PWD/fasta/new_reference_genome.fa <(cat $PWD/fq/$f.h1.final.trimm.fq.gz $PWD/fq/$f.h2.final.trimm.fq.gz) -t 16 |samtools sort -@4 -m4G - > $PWD/bam/$f.bam ;done
# "childAA" "childA" "childB" "childC" "childD" "childE" "childF" "childG" "childH" "childI" "founder0" "founder10" "founder11" "founder12" "founder13" "founder14" "founder15" "founder16" "founder17" "founder18" "founder19" "founder1" "founder20" "founder21" "founder22" "founder23" "founder24"