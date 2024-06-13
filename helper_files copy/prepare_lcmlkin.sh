# conda activate thesis

#### MAKE target_SNP_file - START  ####
path=/projects/korneliussen/people/vlr112/simulations5
output_path=/projects/korneliussen/people/vlr112/simulations5/helper_files 

ancestral=$path/fasta/ancestral.fa
derived=$path/fasta/derived.fa
python variable_pos_list.py > pos
refpos=pos

grep -v ">" $ancestral | tr -d '\n' | fold -w 1 | nl | cut -f2 > ancestral.txt
grep -v ">" $derived | tr -d '\n' | fold -w 1 | nl | cut -f2 > derived.txt

cut $refpos -f2 | sed 's/[[:space:]]\+/\t/g' | cut -f2 | sed 's/^/rs/' > id 
cut $refpos -f2 > position

paste position ancestral.txt derived.txt  >  $output_path/snp_tmp.txt
awk 'BEGIN{FS=OFS="\t"} {print 1, $1, $2, $3, $4}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/target_SNP_file 
rm $output_path/snp_tmp.txt $output_path/derived.txt $output_path/ancestral.txt $output_path/position $output_path/id
#### MAKE target_SNP_file - END  ####

# #### MAKE .FAM FILE - START  ####
# awk 'BEGIN{FS=OFS=" "} {print "familyX", $1, 0, 0, -9, -9}' ../names.txt  > $output_path/file.fam 
# #### MAKE .FAM FILE - END  ####

# # index pos file
# angsd sites index $refpos