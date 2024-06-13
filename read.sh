# Must run on conda env py2 ##
# cmd line: source ./read.sh xx 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py2

# start=$(date +%s.%N)

# basedir=$PWD/$1
# binary_plink=$2
# folder_name=$3

# FOLDER_NAME=$(echo "$folder_name" | tr '[:lower:]' '[:upper:]')

# # ## clean previous results
# # rm -r $basedir/results/READ/

# mkdir -p $basedir/results/$3/
# dir=$basedir/results/$3

# # cd $dir

# # paths to READ
# READ=/projects/korneliussen/people/vlr112/bin/read

# # ln -sf $READ/READ.py $dir/READ.py
# # ln -sf $READ/READscript.R $dir/READscript.R

# # Run READ for unrelated individuals first!
# mkdir -p $dir/unrelated/
# unrelated_dir=$dir/unrelated

# # ln -sf $READ/READ.py $unrelated_dir/READ.py
# # ln -sf $READ/READscript.R $unrelated_dir/READscript.R

# unrelatedlist=$unrelated_dir/unrelatedInds.txt

# # # actually only need 2  unrelated individuals for this part
# # printf "familyX founder0\nfamilyX founder1" > $unrelatedlist

# # # Obtain tped and tfam for the two unrelated individuals 
# # plink --bfile $basedir/$2 --keep $unrelatedlist --aec --recode transpose --out $unrelated_dir/$2.unrelated

# # cd $unrelated_dir

# # # # Run READ for  unrelated individuals
# # python $dir/READ.py $2.unrelated
# # cd ..

# # #### All individuals ####

# # # Obtain tped and tfam for all individuals 
# # plink --bfile $basedir/$2 --aec --recode transpose --out $2.all

# ## normalization_value:
# normalization_value=$(cut -d' ' -f4 $unrelated_dir/meansP0_AncientDNA_normalized | sed -n 2p)

# # # # Run for all individuals
# # python READ.py $2.all value $normalization_value

# # # end=$(date +%s.%N)

# # # echo "$(echo "$end - $start" | bc)" > "$FOLDER_NAME.time"

# # conda deactivate

# # ## plink relatedness
# # plink --bfile $dir/angsdHaploCall.all --recode --transpose --out $dir/angsdHaploCall.out


###### fix ###

## READ -START ##

basedir=$PWD/$1/$2/cov$3
dir=$basedir/results/READ

unrelated_dir=$dir/unrelated


## normalization_value:
normalization_value=$(cut -d' ' -f4 $unrelated_dir/meansP0_AncientDNA_normalized | sed -n 2p)

infile=$basedir/results/READ/meansP0_AncientDNA_normalized
outfile=$basedir/results/READ/meansP0_AncientDNA_normalized_out
tmpfile=$basedir/results/READ/tmp
tmpfile2=$basedir/results/READ/tmp2

## divide by unrelated non-normalized P0
awk -F' ' -v OFS=' ' -v normalization_value="$normalization_value" 'NR>1 { $4 = $4/ normalization_value } 1' $infile > $tmpfile


# sed '1s/^/PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow\tIBS2\tIBS0\tP1\tP0\tMissing\tTMP\n/' $infile > $tmpfile2


awk -v col1=$3 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "cov" : col1) }' $tmpfile > $tmpfile2 
awk -v col1=$2 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "deamination" : col1) }' $tmpfile2 > $tmpfile

awk 'NR==1{$(NF+1)="program"} NR>1{$(NF+1)="READ"}1' $tmpfile | sed 's/[[:space:]]\+/\t/g'  > $tmpfile2 && mv $tmpfile2 $outfile

# sed -i '1s/^/PairIndividuals Normalized2AlleleDifference StandardError NonNormalizedP0 NonNormalizedStandardError
# PairIndividuals\tChromosome\tWindowIndex\tSNVperWindow\tIBS2\tIBS0\tP1\tP0\tMissing\tTMP\n/' $outfile




rm $tmpfile 


