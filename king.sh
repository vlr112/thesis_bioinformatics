source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis

# start=$(date +%s.%N)

basedir=$1
binary_plink=$2
folder_name=$3
original_path=$4
helper_files_path=$PWD/helper_files

FOLDER_NAME=$(echo "$folder_name" | tr '[:lower:]' '[:upper:]')

# ## clean previous results
# rm -r $basedir/results/READ_lcmlkin/

mkdir -p $basedir/results/$folder_name/
dir=$basedir/results/$folder_name
# cd $dir
# # # king -b $basedir/$binary_plink.bed --related --degree 3 
# # # king -b $basedir/$binary_plink.bed --ibs
# # # king -b $basedir/$binary_plink.bed --ibdseg --degree 3
# king -b $basedir/$binary_plink.bed --related 

# python $helper_files_path/fix_KING.py $dir $basedir/results/relate/real_data_out  king.kin
# mv king.kin2 king.kin
# cd $original_path
# cut -f1 king.kin | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp
# cut -f2 king.kin | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp2
# cut -f3 king.kin | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp3
# cut -f4 king.kin | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp4

# awk -v i=1  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp king.kin > tm && mv tm king.kin0


# # # cut -f2 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp2
# awk -v i=2 'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp2 king.kin > tm && mv tm king.kin0


# # # cut -f3 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp
# awk -v i=3  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp3 king.kin > tm && mv tm king.kin0


# # # cut -f4 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp
# awk -v i=4  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp4 king.kin > tm && mv tm king.kin0

# sed -i 's/[[:space:]]\+/\t/g' king.kin
# rm tmp tmp2 tmp3 tmp4

# end=$(date +%s.%N)

# echo "$(echo "$end - $start" | bc)" > "$FOLDER_NAME.time"

# cd $original_path


#### now for between family inference ####

# plink --vcf $basedir/allbams_final.vcf.gz --recode --double-id --make-bed --out $basedir/allbams_final_KING_bf
cd $dir
king -b $basedir/$binary_plink.bed --kinship 

cut -f1 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp
cut -f2 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp2
cut -f3 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp3
cut -f4 king.kin0 | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' > tmp4

awk -v i=1  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp king.kin0 > tm && mv tm king.kin0
awk -v i=2 'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp2 king.kin0 > tm && mv tm king.kin0
awk -v i=3  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp3 king.kin0 > tm && mv tm king.kin0
awk -v i=4  'FNR==NR{a[NR]=$1;next}{$i=a[FNR]}1' tmp4 king.kin0 > tm && mv tm king.kin0
sed -i 's/[[:space:]]\+/\t/g' king.kin0
rm tmp tmp2 tmp3 tmp4

# python $helper_files_path/fix_KING.py $dir $basedir/results/relate/real_data_out  king.kin0

python $helper_files_path/fix_KING.py $dir $PWD/sim20/no_deam/cov20/results/relate/real_data_out  king.kin0


mv king.kin02 king.kin0
cd $original_path

