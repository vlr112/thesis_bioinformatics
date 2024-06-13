

basedir=$PWD/$1/$2/cov$3 ## something like
bam_path=$basedir/bam

# for f in $(ls $bam_path/*bam)
for f in "childAA.bam" "childA.bam" "childB.bam" "founder10.bam"
do
    samtools index $f
done |parallel