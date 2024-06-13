#conda env thesis

source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis

# run eg: ./lcmlkin.sh sim1/no_deam/cov10
start=$(date +%s.%N)
# bams to vcf
basedir=$PWD/$1
helper_files_path=$PWD/helper_files
base_path=/projects/korneliussen/people/vlr112/final_simulations

# ## clean previous results
# rm -r $basedir/results/lcmlkin/
# rm $basedir/lcmlkin_mp.log

mkdir -p $basedir/results/lcmlkin/
dir=$basedir/results/lcmlkin
blist=$dir/bamlist.txt
ls $basedir/bam/*.bam > $blist

# #run it once
# # ./helper_files/prepare_lcmlkin.sh 
cp $helper_files_path/target_SNP_file $dir/target_SNP_file

# cp $helper_files_path/SNPbam2vcf.py $dir/SNPbam2vcf.py
cp $helper_files_path/SNPbam2vcf.francisca.changed.py $dir/SNPbam2vcf.francisca.changed.py

cd $dir

# # # # This works fine, but takes over 5h!!!
# # # # ./SNPbam2vcf.py $blist  $dir/allbams.vcf $dir/target_SNP_file

## So lets try with python multiprocessing Pool thing. This took a little less than 13 mins!!! Hurray!
./SNPbam2vcf.francisca.changed.py $blist  $dir/allbams_mp.vcf $dir/target_SNP_file

python $helper_files_path/updateAF_lcmlkinVCF.py $helper_files_path/freq_tmp.txt $dir/allbams_mp.vcf $dir/allbams_mp_update.vcf

lcmlkin -i $dir/allbams_mp_update.vcf -o lcmlkin_results_update -g all

# lcmlkin -i $basedir/allbams_final.vcf -o $dir/lcmlkin_results -g all -l phred


# #since we're here, obtain the binary plink files
# plink --vcf $dir/allbams_mp.vcf  --recode --make-bed --out $basedir/allbams_lcmlkin

# cp $dir/allbams_mp.vcf  $basedir/allbams_lcmlkin.vcf

# ## need to substitute the family column in .fam file with familyX
# awk '{ $1 = "familyX" }1' $basedir/allbams_lcmlkin.fam > $basedir/tmp && mv $basedir/tmp $basedir/allbams_lcmlkin.fam
end=$(date +%s.%N)

echo "$(echo "$end - $start" | bc)" > "LCMLKIN.time"

cd $base_path

conda deactivate



# # # rm $basedir/lcmlkin_mp.log

# # mkdir -p $basedir/results/lcmlkin/
# # dir=$basedir/results/lcmlkin
# # blist=$dir/bamlist.txt
# # ls $basedir/bam/*.bam > $blist

# # #run it once
# # # ./helper_files/prepare_lcmlkin.sh 
# # cp $helper_files_path/target_SNP_file $dir/target_SNP_file

# # # cp $helper_files_path/SNPbam2vcf.py $dir/SNPbam2vcf.py
# # cp $helper_files_path/SNPbam2vcf.francisca.changed.py $dir/SNPbam2vcf.francisca.changed.py

# # cd $dir

# # # # # # This works fine, but takes over 5h!!!
# # # # # # ./SNPbam2vcf.py $blist  $dir/allbams.vcf $dir/target_SNP_file

# # ## So lets try with python multiprocessing Pool thing. This took a little less than 13 mins!!! Hurray!
# # ./SNPbam2vcf.francisca.changed.py $blist  $dir/allbams_mp.vcf $dir/target_SNP_file

# # lcmlkin -i $dir/allbams_mp.vcf -o lcmlkin_results -g all

# # #since we're here, obtain the binary plink files
# # plink --vcf $dir/allbams_mp.vcf  --recode --make-bed --out $basedir/allbams_lcmlkin

# # cp $dir/allbams_mp.vcf  $basedir/allbams_lcmlkin.vcf

# # ## need to substitute the family column in .fam file with familyX
# # awk '{ $1 = "familyX" }1' $basedir/allbams_lcmlkin.fam > $basedir/tmp && mv $basedir/tmp $basedir/allbams_lcmlkin.fam
# # end=$(date +%s.%N)

# # echo "$(echo "$end - $start" | bc)" > "LCMLKIN.time"

# # cd $base_path
