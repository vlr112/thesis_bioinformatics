#!/usr/bin/python
import numpy as np
import random
import time
import glob
import sys

# Takes in a "variable sites" file
geno_in = sys.argv[1]

nucs = ["A", "C", "G", "T"]
last_time = time.time()

# Nice print
def msg(m):
    global last_time
    this_time = time.time()
    elapsed = this_time - last_time
    print("%s - took %.5f secs" % (m, elapsed))
    last_time = this_time

# Read in variable sites to a list
def read_geno(geno_in):
    variable_sites = []
    with open (geno_in, 'r') as file:
        for i in file.read():
            nuc = i.rstrip('\n')
            if nuc != "":
                variable_sites.append(nuc)
    return variable_sites


variable_sites = read_geno(geno_in)
msg("Read variable sites")

# Read the reference genome
with open("fasta/reference_genome.fa") as f:
    header = f.readline().rstrip(">").rstrip('\n')
    seq = f.readline().rstrip('\n')
msg("Read in reference")

# Loop through the reference
final_seq = []
c = 0
for nuc in list(seq):
    # Whenever the reference allele is an N, this means the site is variable
    # so we use the site from our variable sites list instead
    if nuc == "N":
        final_seq.append(variable_sites[c])
        c += 1
    else:
        final_seq.append(nuc)
msg("Made final sequence")

# Now we save the whole thing as a fasta file
outname = geno_in.split("/")[-1].rstrip(".txt")
out = "final_fasta/" + outname + ".final.fa"
with open (out, 'w') as file:
    header = ">" + outname + "---final\n"
    file.write(header)
    file.write("".join(final_seq))
    file.write("\n")
