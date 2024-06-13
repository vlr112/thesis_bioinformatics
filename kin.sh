# ## clean previous results
# rm -r $PWD/$basedir/results/kin/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate kin-3.1.3

basedir=$PWD/$1
# rm -r $basedir/results/kin/
original_path=$PWD
helper_files_path=$PWD/helper_files
unrelated_ind=$PWD/unrelated.txt

mkdir -p $basedir/results/kin/
mkdir -p $basedir/results/kin/bam
dir=$basedir/results/kin
names=$PWD/names_copy.txt

cp $helper_files_path/snp.txt $dir/bedfile.bed

mkdir -p $dir/bam/unrelated/bam
unrelated=$dir/bam/unrelated/bam
cp $helper_files_path/snp.txt $unrelated/bedfile.bed

# ln -sf $basedir/bam/founder0.bam $unrelated/founder0.bam
# ln -sf $basedir/bam/founder0.bam.bai $unrelated/founder0.bam.bai
# ln -sf $basedir/bam/founder1.bam $unrelated/founder1.bam
# ln -sf $basedir/bam/founder1.bam.bai $unrelated/founder1.bam.bai

# cd $unrelated

# # KINgaroo -bam $PWD -bed $dir/bedfile.bed -T $unrelated_ind -cnt 0 -c 8 -r 1 -N 1 -s 0
# KINgaroo -bam $PWD -bed $unrelated/bedfile.bed -T $unrelated_ind -i 1000000 -cnt 0 -c 8 -N 1 -s 0

cd $dir/bam/

while read f; do
    ln -sf "$basedir/bam/$f.bam" "$dir/bam/$f.bam"
    ln -sf "$basedir/bam/$f.bam.bai" "$dir/bam/$f.bam.bai"
done < $names

p_0=$(cat $unrelated/hmm_parameters/p_0.txt)
filtered_windows=$(cat $unrelated/filtered_windows.txt)

# KINgaroo -bam $PWD -bed $dir/bedfile.bed -T $unrelated_ind -p $p_0 -cnt 0 -c 8 -r 1 -N 1 -s 0
# KIN -p $p_0 -I $PWD -O $PWD/KIN_results/ -r $PWD 
KINgaroo -bam $PWD -bed $dir/bedfile.bed -T $names -i 1000000 -p $p_0 -cnt 0 -c 10 -N 1 -s 0 #-n $filtered_windows
KIN -I $PWD -O $PWD/KIN_results/ -p $p_0 #-i 1000000

# # # cd $dir
# cd $original_path

# conda deactivate









# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate kin-3.1.3
# # source ./kin.sh cov10

# start=$(date +%s.%N)

# basedir=$1
# helper_files_path=$PWD/helper_files
# names=$PWD/names.txt

# # ## clean previous results
# # rm -r $basedir/results/kin/


# mkdir -p $basedir/results/kin/
# mkdir -p $basedir/results/kin/bam
# dir=$basedir/results/kin


# blist=$dir/bamlist.txt
# ls $basedir/bam/*.bam > $blist

# # ./make_bed.sh -> run once in helper_files directory

# cp $helper_files_path/snp.txt $dir/bedfile.bed
# # for i in $(ls $basedir/bam/*); do ln -s  $dir/bam/

# while read f; do
#     ln -sf "$basedir/bam/$f.bam" "$dir/bam/$f.bam"
#     ln -sf "$basedir/bam/$f.bam.bai" "$dir/bam/$f.bam.bai"
# done < $names

# cd $dir/bam/

# KINgaroo -bam $PWD -bed $dir/bedfile.bed -T $names -cnt 0 -c 30 -r 1 -N 1

# KIN -I $PWD -O $PWD/KIN_results/ -r $PWD

# end=$(date +%s.%N)

# cd $dir

# echo "$(echo "$end - $start" | bc)" > "KIN.time"

# conda deactivate