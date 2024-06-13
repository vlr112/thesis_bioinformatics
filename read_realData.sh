source ~/miniconda3/etc/profile.d/conda.sh
conda activate py2

# ## REAL FAROESE DATA ##
# basedir=$PWD/faroe
# mkdir -p $basedir/results/READ/
# dir=$basedir/results/READ

# ln -sf $READ/READ.py $dir/READ.py
# ln -sf $READ/READscript.R $dir/READscript.R
# cd $dir

# # # Obtain tped and tfam for all individuals 
# plink --bfile $basedir/allbams_final --aec --recode transpose --out allbams_final.all

# ## normalization_value:
# # normalization_value=$(cut -d' ' -f4 $unrelated_dir/meansP0_AncientDNA_normalized | sed -n 2p)

# # # Run for all individuals
# python READ.py allbams_final.all 

# conda deactivate

# ##################################



basedir=$PWD/faroe
mkdir -p $basedir/results/READ2/
dir=$basedir/results/READ2

cd $dir

# paths to READ
READ=/projects/korneliussen/people/vlr112/bin/read

ln -sf $READ/READ.py $dir/READ.py
ln -sf $READ/READscript.R $dir/READscript.R

# Run READ for unrelated individuals first!
mkdir -p $dir/unrelated/
unrelated_dir=$dir/unrelated

ln -sf $READ/READ.py $unrelated_dir/READ.py
ln -sf $READ/READscript.R $unrelated_dir/READscript.R

unrelatedlist=$unrelated_dir/unrelatedInds.txt

# actually only need 2  unrelated individuals for this part
printf "VK24 VK24\nVK239 VK239" > $unrelatedlist

# Obtain tped and tfam for the two unrelated individuals 
plink --bfile $basedir/allbams_final --keep $unrelatedlist --aec --recode transpose --out $unrelated_dir/allbams_final.unrelated

cd $unrelated_dir

# # Run READ for  unrelated individuals
python $dir/READ.py allbams_final.unrelated
cd ..

#### All individuals ####

# Obtain tped and tfam for all individuals 
plink --bfile $basedir/allbams_final --aec --recode transpose --out allbams_final.all

## normalization_value:
normalization_value=$(cut -d' ' -f4 $unrelated_dir/meansP0_AncientDNA_normalized | sed -n 2p)

# # Run for all individuals
python READ.py allbams_final.all value $normalization_value

conda deactivate




