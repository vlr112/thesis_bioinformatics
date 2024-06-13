## join various files to make snp file -> use in correctKin importHaploCall 

#### MAKE EIGENSTRAT.SNP - START  ####
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

paste id position derived.txt ancestral.txt >  $output_path/snp_tmp.txt
awk 'BEGIN{FS=OFS="\t"} {print $1, 1, 0, $2, $3, $4}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/ /g' | cut -f1- > $output_path/file.snp
rm $output_path/snp_tmp.txt $output_path/derived.txt $output_path/ancestral.txt $output_path/position $output_path/id
#### MAKE EIGENSTRAT.SNP - END  ####

#### MAKE .FAM FILE - START  ####
for i in $(ls $path/final_fasta/*.fa); do basename $i | cut -d'.' -f1 ; done | uniq > $output_path/names.txt
awk 'BEGIN{FS=OFS=" "} {print "familyX", $1, 0, 0, -9, -9}' $output_path/names.txt  > $output_path/file.fam 
#### MAKE .FAM FILE - END  ####

# index pos file
angsd sites index $refpos
