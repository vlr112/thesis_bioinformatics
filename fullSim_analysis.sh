## takes results obtained by fullSim.sh and conducts the analysis for various programs

### NGSRELATE - START ###
# original_path=$PWD

# mkdir -p $PWD/$1/$2/cov$3
# basedir=$PWD/$1/$2/cov$3 ## something like
# helper_files_path=$PWD/helper_files


source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis

# # aggregate all data_out and so on from all simulations
# # rm -r results/relate/*.txt


# ### NGSRELATE alternative - START ###
# mkdir -p results/relateREAL
# dir1=$PWD/results/relateREAL
# find . -name "real_data_out" > $dir1/relateREAL_list.txt
# python helper_files/join_sameProgram.py $dir1/relateREAL_list.txt $dir1/relateREAL.txt 
# # find . -name "real_data_out" -exec cat {} + > $dir1/ngsrelateREAL.txt
# ### NGSRELATE alternative - END ###

# ### PLINK - START ###
# mkdir -p results/plink
# dir2=$PWD/results/plink
# find . -name "plink_out.genome" > $dir2/plink_list.txt
# python helper_files/join_sameProgram.py $dir2/plink_list.txt $dir2/plink.txt
# ### PLINK - START ###

# ### LCMLKIN - START ###
# mkdir -p results/lcmlkin
# dir3=$PWD/results/lcmlkin
# find . -name "lcmlkin_results_update_out" > $dir3/lcmlkin_list.txt
# python helper_files/join_sameProgram.py $dir3/lcmlkin_list.txt $dir3/lcmlkin.txt
# ## LCMLKIN - END ###

# ### CORRECTKIN - START ###
# mkdir -p results/correctKin
# dir4=$PWD/results/correctKin
# find . -name "out_angsdHaploCall.relatives.tsv" > $dir4/correctKin_list.txt
# python helper_files/join_sameProgram.py $dir4/correctKin_list.txt $dir4/correctKin.txt
# ## CORRECTKIN - END ###

# dir_list=("plink" "lcmlkin" "relateREAL")
# python merger.py "${dir_list[@]}" $PWD/results/merged_ALL.txt

# ##### KING - START #####
# mkdir -p results/king
# dir4=$PWD/results/king
# find . -name "king_out.kin" > $dir4/king_list.txt
# python helper_files/join_sameProgram.py $dir4/king_list.txt $dir4/king.txt
# ##### KING - END #####

# ### NGSRELATE True & Angsd & random - START ###
# mkdir -p results/relateREAL
# dir1=$PWD/results/relateREAL
# find . -name "real_data_out" -o -name "angsd_data_out" -o -name "randomFreq_data_out"  > $dir1/compareRelate_list.txt
# python helper_files/join_sameProgram.py $dir1/compareRelate_list.txt $dir1/compareRelate.txt 
# ### NGSRELATE True & Angsd - END ###



### READ - START ###
mkdir -p results/READ
dir4=$PWD/results/READ
find . -name "meansP0_AncientDNA_normalized_out" > $dir4/READ_list.txt
python helper_files/join_sameProgram.py $dir4/READ_list.txt $dir4/READ.txt
## READ - END ###


# ### KING within and between family - START ###
# mkdir -p results/king
# dir1=$PWD/results/king
# # find . -name "king_out.kin" -o -name "king_out.kin0"  > $dir1/compareKing_list.txt
# python helper_files/join_sameProgram.py $dir1/compareKing_list.txt $dir1/compareKing.txt 
# ### KING within and between family - END ###

# # python helper_files/fix.ngsrelate.plink.lcmlkin.py $dir1/ngsrelate.txt $dir1/ &
# # python helper_files/fix.ngsrelate.plink.lcmlkin.py $dir2/plink.txt $dir2/ &
# # python helper_files/fix.ngsrelate.plink.lcmlkin.py $dir3/lcmlkin.txt $dir3/ &


# # python helper_files/fix.ngsrelate.plink.lcmlkin.py $dir1/ngsrelateREAL.txt $dir1/ &
# # python helper_files/fix.ngsrelate.plink.lcmlkin.py $dir2/plink_lcmlkin.txt $dir2/ &

# # wait

# # rm $dir/kinship_degree.txt
# # cd results/
# # # loop through subdirectories
# # this will loop over the 3 subdicts (plink, lcmlkin and ngsrelate) in dict results and merge all files 
# # with same relatedness, fx. all FS files will be merged together. 

# # dir_list=("plink" "plink_lcmlkin")
# # cmd='files=($(find $PWD/results/{plink, plink_lcmlkin} -mindepth 1 -type f -name "${prefix}_*"))'






# # dir_str=$(IFS=','; echo "${dir_list[*]}")

# # for prefix in A_2 FC_3 FS_1 GP_2 HC_4 HS_2 PO_1 Unrelated_Unrelated; do
# #     python merger.py "$prefix" "$dir_str"
# # done


# ### lets look into rmse for true ibd  (from fasta files)

# mkdir -p results/true_rmse
# dir1=$PWD/results/true_rmse
# # # find . -name "rmse_true" -exec cat {} + > $dir1/rmse_true.txt
# # find $PWD -name "rmse_true"  > $dir1/rmse_true.txt


# python helper_files/fix.rmse_true.py $dir1/rmse_true.txt 


# # # ## PLINK lcmlkin - START ###
# # # mkdir -p results/plink_lcmlkin
# # # dir2=$PWD/results/plink_lcmlkin
# # # find . -name "plink_lcmlkin_out.genome" -exec cat {} + > $dir2/plink_lcmlkin.txt
# # # ## PLINK lcmlkin - START ###

# # # ### NGSRELATE - START ###
# # mkdir -p results/relate
# # dir1=$PWD/results/relate
# # # find . -name "data_out" -exec cat {} + > $dir1/ngsrelate.txt
# # # ### NGSRELATE - END ###