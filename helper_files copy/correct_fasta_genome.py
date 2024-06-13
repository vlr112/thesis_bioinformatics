from Bio import SeqIO
# from Bio import fastaSequence
from sys import argv


###Input arguments
original_genome=argv[1] 
ancestral_positions=argv[2] 
outfile=argv[3]

# Open the reference genome file and read the sequence
with open(original_genome) as f:
    ref_seq = f.read().strip().split('\n')[1]  # assuming the first line is a header

# Open the substitution file and read the sequence
with open(ancestral_positions) as f:
    subs_seq = f.read().strip().split('\n')[1]  # assuming the first line is a header

# Replace each N base in the reference sequence with the corresponding base in the substitution sequence
new_seq = ''
i = 0
for base in ref_seq:
    if base == 'N':
        new_seq += subs_seq[i]
        i += 1
    else:
        new_seq += base

# Write the new sequence to a new file
with open(outfile, 'w') as f:
    f.write('>1\n' + new_seq)


