
# replace base N with ancestral bases.
path=$PWD
output_path=$path/helper_files 

# python $output_path/correct_fasta_genome.py $path/fasta/reference_genome.fa $path/fasta/ancestral.fa $path/fasta/new_reference_genome.fa

# cd $path/fasta

# bwa index $path/fasta/new_reference_genome.fa

# cd $path


# ## preprare quality profile for modern DNA with need at least 350 length cycle reads.
# cat AccFreqL150R1.txt | head -2 > header.txt
# cat AccFreqL150R1.txt | tail -n+3|awk '{for(i=0;i<3;i++)print}' > qualprob.txt
# cat header.txt qualprob.txt > newqualprofile_R1.txt

########################################################################################################################################################3

# ancestral=$path/fasta/ancestral.fa
# derived=$path/fasta/derived.fa
# python $output_path/variable_pos_list.py > $output_path/pos
# refpos=$output_path/pos

# grep -v ">" $ancestral | tr -d '\n' | fold -w 1 | nl | cut -f2 > $output_path/ancestral.txt
# grep -v ">" $derived | tr -d '\n' | fold -w 1 | nl | cut -f2 > $output_path/derived.txt

# cut $refpos -f2 | sed 's/[[:space:]]\+/\t/g' | cut -f2 | sed 's/^/rs/' > $output_path/id 
# cut $refpos -f2 > $output_path/position


# #### MAKE REAL FREQUENCY FILE: POSITION FREQ - START####
# paste  $output_path/position af.txt |  sed 's/[[:space:]]\+/\t/g' > $output_path/freq_tmp.txt
# sed -i '1iposition\tknownEM' $output_path/real_freq.txt
# #### MAKE REAL FREQUENCY FILE: POSITION FREQ - END ####


# #### MAKE EIGENSTRAT.SNP - START  ####
# paste $output_path/id $output_path/position $output_path/derived.txt $output_path/ancestral.txt >  $output_path/snp_tmp.txt
# awk 'BEGIN{FS=OFS="\t"} {print $1, 1, 0, $2, $3, $4}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/ /g' | cut -f1- > $output_path/file.snp
#### MAKE EIGENSTRAT.SNP - END  ####

#### MAKE .FAM FILE - START  ####
for i in $(ls $path/final_fasta/*.fa); do basename $i | cut -d'.' -f1 ; done | uniq > $output_path/names.txt
awk 'BEGIN{FS=OFS=" "} {print "familyX", $1, 0, 0, -9, -9}' $output_path/names.txt  > $output_path/file.fam 
#### MAKE .FAM FILE - END  ####

# # index pos file
# angsd sites index $refpos

# required for KIN programm
for i in $(ls $path/final_fasta/*.fa); do basename $i | cut -d'.' -f1 ; done | uniq > $path/names.txt


# ### MAKE file for KIN programm - START ###
# paste $output_path/position $output_path/ancestral.txt $output_path/derived.txt  >  $output_path/snp_tmp.txt
# awk 'BEGIN{FS=OFS="\t"} {print "1", $1, $1+1, $2, $3}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/snp.txt
# ### MAKE file for KIN programm - END ###


# ###### prepare_lcmlkin
# #### MAKE target_SNP_file - START  ####
# paste $output_path/position $output_path/ancestral.txt $output_path/derived.txt  >  $output_path/snp_tmp.txt
# awk 'BEGIN{FS=OFS="\t"} {print 1, $1, $2, $3, $4}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/target_SNP_file 
# #### MAKE target_SNP_file - END  ####


# # #real BED file
# # awk 'BEGIN{FS=OFS="\t"} {print "1", $1, $1+1}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/snp.bed
# # awk 'BEGIN{FS=OFS="\t"} {print "1", $1-1, $1}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/file.bed

rm $output_path/snp_tmp.txt $output_path/derived.txt $output_path/ancestral.txt $output_path/position $output_path/id


#### MAKE FILE WITH IBD0, IBD1, IBD2 FROM ORIGINAL FASTA FILES ####
python $output_path/IBD_true.py $path/names.txt IBD_twin
