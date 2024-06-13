# start=$(date +%s.%N)

basedir=$PWD/faroe
binary_plink=allbams_test
# folder_name=
original_path=$PWD

# FOLDER_NAME=$(echo "$folder_name" | tr '[:lower:]' '[:upper:]')

# ## clean previous results
# rm -r $basedir/results/READ_lcmlkin/

mkdir -p $basedir/results/plink/
dir=$basedir/results/plink
cd $dir

plink --bfile $basedir/$binary_plink --genome
# end=$(date +%s.%N)

# echo "$(echo "$end - $start" | bc)" > "$FOLDER_NAME.time"

cd $original_path

