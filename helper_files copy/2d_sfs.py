#!/usr/bin/python
import numpy as np
import random
import time
import glob
import sys

ind1 = sys.argv[1]
ind2 = sys.argv[2]
last_time = time.time()

# Read which nucleotide is ancestral
def get_ancestral():
    with open("fasta/ancestral.fa") as f:
        header = f.readline().rstrip(">").rstrip('\n')
        seq = f.readline().rstrip('\n')
        return list(seq)

# Read in genotypes
def read_geno(geno_in):
    variable_sites = []
    with open (geno_in, 'r') as file:
        for i in file.read():
            nuc = i.rstrip('\n')
            if nuc != "":
                variable_sites.append(nuc)
    return variable_sites

# Allele genotype to 0,1,2 form
def get_geno(g, a):
    if g[0] == a and g[1] == a:
        return 0
    if g[0] == g[1]:
        return 2
    if g[0] != g[1]:
        return 1

v1_h1_geno = read_geno((ind1 + ".h1.txt"))
v2_h1_geno = read_geno((ind2 + ".h1.txt"))
v1_h2_geno = read_geno((ind1 + ".h2.txt"))
v2_h2_geno = read_geno((ind2 + ".h2.txt"))
ancestral = get_ancestral()

sfs = np.zeros(shape=(3,3))
# Compare genotypes at each pos to get sfs
for v1_h1, v2_h1, v1_h2, v2_h2, a in zip(v1_h1_geno, v2_h1_geno, v1_h2_geno, v2_h2_geno, ancestral):
    g1 = [v1_h1, v1_h2]
    g2 = [v2_h1, v2_h2]
    ix1 = get_geno(g1,a)
    ix2 = get_geno(g2,a)
    sfs[ix1,ix2] += 1
print(sfs)
