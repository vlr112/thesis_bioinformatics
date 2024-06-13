source ~/miniconda3/etc/profile.d/conda.sh
conda activate thesis

start=$(date +%s.%N)

basedir=$1
binary_plink=$2
folder_name=$3
original_path=$4

FOLDER_NAME=$(echo "$folder_name" | tr '[:lower:]' '[:upper:]')

# ## clean previous results
# rm -r $basedir/results/READ_lcmlkin/

mkdir -p $basedir/results/$folder_name/
dir=$basedir/results/$folder_name
cd $dir

plink --bfile $basedir/$binary_plink --genome
end=$(date +%s.%N)

echo "$(echo "$end - $start" | bc)" > "$FOLDER_NAME.time"

cd $original_path
