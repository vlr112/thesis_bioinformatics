# program KIN requires .bed file

## join various files to make snp file -> use in correctKin importHaploCall 
path=/projects/korneliussen/people/vlr112/simulations5
output_path=/projects/korneliussen/people/vlr112/simulations5/helper_files 

ancestral=$path/fasta/ancestral.fa
derived=$path/fasta/derived.fa
python variable_pos_list.py > pos
refpos=pos

grep -v ">" $ancestral | tr -d '\n' | fold -w 1 | nl | cut -f2 > ancestral.txt
grep -v ">" $derived | tr -d '\n' | fold -w 1 | nl | cut -f2 > derived.txt

cut $refpos -f2 > position

paste position ancestral.txt derived.txt  >  $output_path/snp_tmp.txt

awk 'BEGIN{FS=OFS="\t"} {print "1", $1, $1+1, $2, $3}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/snp.txt

#real BED file
awk 'BEGIN{FS=OFS="\t"} {print "1", $1, $1+1}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/snp.bed

awk 'BEGIN{FS=OFS="\t"} {print "1", $1-1, $1}' $output_path/snp_tmp.txt | sed 's/[[:space:]]\+/\t/g' | cut -f1- > $output_path/file.bed

rm $output_path/snp_tmp.txt $output_path/derived.txt $output_path/ancestral.txt $output_path/position 
