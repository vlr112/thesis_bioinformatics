


basedir=$1
dir=$basedir/results/relate
python helper_files/change_randomAF.py $dir/freq $dir/random_freq
ngsrelate  -g $dir/data.glf.gz -n 40 -f $dir/random_freq  -O $dir/randomFreq_data
