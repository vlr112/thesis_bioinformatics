start=$(date +%s.%N)

pcangsdPath=/projects/korneliussen/people/vlr112/bin/pcangsd-v.0.99
correctKinPath=/projects/korneliussen/people/vlr112/bin/correctKin/bin

# haploToPlink=/home/vlr112/prog/angsd/misc/haploToPlink
# path007=/projects/korneliussen/people/vlr112/final_simulations/
path0=/projects/korneliussen/people/vlr112/final_simulations/helper_files 
# # path=/projects/korneliussen/people/vlr112/simulations5

basedir=$PWD/$1
# mkdir -p $basedir/results/correctKin/
dir=$basedir/results/correctKin
# blist22=$basedir/bamlist22.txt


# # # ./make_snp.sh  
# cp $path0/file.snp $dir/angsdHaploCall.snp
# cp $path0/file.fam $dir/angsdHaploCall.fam

# ## ALREADY RAN IN HELPER_FILES/MAKE_BED.SH
# # python scripts/variable_pos_list.py > pos
# refpos=$path0/pos
# # /home/vlr112/prog/angsd/angsd sites index $refpos

# conda activate thesis

# cd $dir
# # random allele calling with ANGSD using the prepared sites file
# angsd -doHaploCall 1 -sites $refpos -doCounts 1 -bam $blist -nThreads 16 -out $dir/angsdHaploCall

# # importing ANGSD haploid call output into PLINK binary data set
# $correctKinPath/importHaploCall $dir/angsdHaploCall.snp $dir/angsdHaploCall.bed $dir/angsdHaploCall.haplo.gz

# # put all .log in separate folder
# mkdir -p $basedir/results/correctKin/individuals/
# ind=$basedir/results/correctKin/individuals
# mv *log $ind/

# cd $path007

# cp $path0/go.mod $dir/
# cp $path0/go.sum $dir/

# ##NOTE: cannot install in thesis conda env.
# source ~/miniconda3/etc/profile.d/conda.sh 
# conda deactivate
# # performing kinship analysis with PCangsd v0.99
# python $pcangsdPath/pcangsd.py -plink $dir/angsdHaploCall -o $dir/angsdHaploCall -inbreed 1 -kinship

# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate thesis
# # # calculating the pairwise marker overlap fraction matrix of samples
# $correctKinPath/markerOverlap $dir/angsdHaploCall.bed

# # ln -s $correctKinPath/correctKin/cmd/filterRelates/main $dir/
# # ln -s $correctKinPath/correctKin/cmd/filterRelates/main.go $dir/


ln -sf /projects/korneliussen/people/vlr112/bin/correctKin/cmd/filterRelates/main $dir/main
ln -sf /projects/korneliussen/people/vlr112/bin/correctKin/cmd/filterRelates/main.go $dir/main.go



# # using the estimated kinship coefficient matrix and marker overlap fraction matrix 
# #to calculate the corrected kinship coefficient matrix, statistical analysis 
# #of corrected kinship estimates, and filter relatives
# # $correctKinPath/filterRelates  $dir/angsdHaploCall.overlap $dir/angsdHaploCall.kinship.npy  
# # > $dir/angsdHaploCall.relatives.tsv
$dir/main -sigma 0 $dir/angsdHaploCall.overlap $dir/angsdHaploCall.kinship.npy > $dir/angsdHaploCall.relatives.tsv

# end=$(date +%s.%N)

# echo "$(echo "$end - $start" | bc)" > "$basedir/results/correctKin/CORRECTKIN.time"

# cd $path007
# ## NOW RUN READ PROGRAM WITH THE INFORMATION FROM CORRECTKIN

# # $haploToPlink $dir/angsdHaploCall.haplo.gz $dir/coisa

#########

# basedir=$PWD/faroe
# mkdir -p $basedir/results/correctKin/
# dir=$basedir/results/correctKin
# blist22=$basedir/bamlist22.txt

# cd $dir
# angsd sites index $basedir/v42.4.1240K_autosomal.sites

# # random allele calling with ANGSD using the prepared sites file
# angsd -doHaploCall 1 -sites $basedir/v42.4.1240K_autosomal.sites -doCounts 1 -bam $blist22 -nThreads 16 -out $dir/GT_1240

# # subset for chr22
# awk '$2 == 22' $basedir/v42.4.1240K.snp > $basedir/v42.4.1240K.chr22.snp

# # importing ANGSD haploid call output into PLINK binary data set
# $correctKinPath/importHaploCall $basedir/v42.4.1240K.chr22.snp GT_1240_haploid.bed GT_1240.haplo.gz

# # put all .log in separate folder
# mkdir -p $basedir/results/correctKin/individuals/
# ind=$basedir/results/correctKin/individuals
# mv *log $ind/

# awk 'BEGIN{FS=OFS=" "} {print $1, $1, 0, 0, -9, -9}' $basedir/samples.txt  > $dir/GT_1240_haploid.fam


# ##NOTE: cannot install in thesis conda env.
# source ~/miniconda3/etc/profile.d/conda.sh 
# conda deactivate
# # performing kinship analysis with PCangsd v0.99
# python $pcangsdPath/pcangsd.py -plink $dir/GT_1240_haploid -o $dir/GT_1240_haploid -inbreed 1 -kinship

# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate thesis
# cp $path0/go.mod $dir/
# cp $path0/go.sum $dir/

# ln -sf /projects/korneliussen/people/vlr112/bin/correctKin/cmd/filterRelates/main $dir/main
# ln -sf /projects/korneliussen/people/vlr112/bin/correctKin/cmd/filterRelates/main.go $dir/main.go


# $correctKinPath/markerOverlap GT_1240_haploid.bed

# $correctKinPath/main  GT_1240_haploid.overlap GT_1240_haploid.kinship.npy > GT_1240_haploid.relatives.tsv



