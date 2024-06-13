#!/usr/bin/python
import numpy as np
import random
import time
import glob
import sys
from tabulate import tabulate
from multiprocessing import Pool

# last_time = time.time()

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

def compare_individuals(pair):
    # pair is a tuple containing two individuals to compare
    ind1, ind2 = pair
    # perform pairwise comparison of individual1 and individual2
    # return the result

    v1_h1_geno = read_geno(('variable_sites/' + ind1 + ".h1.txt"))
    v1_h2_geno = read_geno(('variable_sites/' + ind1 + ".h2.txt"))
    v2_h1_geno = read_geno(('variable_sites/' + ind2 + ".h1.txt"))
    v2_h2_geno = read_geno(('variable_sites/' + ind2 + ".h2.txt"))
    ancestral = get_ancestral()

    ibd0=0
    ibd1=0
    ibd2=0
    # sfs = np.zeros(shape=(3,3))
    # Compare genotypes at each pos to get sfs

    i=0
    for v1_h1, v1_h2, v2_h1, v2_h2, a in zip(v1_h1_geno, v1_h2_geno,  v2_h1_geno, v2_h2_geno, ancestral):

        g1 = [v1_h1, v1_h2]
        g2 = [v2_h1, v2_h2]

        g1.sort()
        g2.sort()
        if g1 == g2:
            ibd2+=1
        elif g1[0] == g2[0] or g1[1] == g2[1]:
            ibd1 +=1
        else:
            ibd0+=1

    total= ibd2 +ibd1 + ibd0

    IBD0=ibd0/total
    IBD1=ibd1/total
    IBD2=ibd2/total

    out = str(ind1) + '\t' + str(ind2) + '\t' + str(IBD0) + '\t' + str(IBD1) + '\t' + str(IBD2)
    out = out + '\n'

    # return ind1, ind2, IBD0, IBD1, IBD2
    return out

if __name__ == '__main__':

    names = sys.argv[1]
    filenameout = sys.argv[2]
    outfile=open(filenameout,'w')
    out='a\tb\tJ9\tJ8\tJ7'
    outfile.write(out)   
    outfile.write('\n')

    with open(names, "r") as file:
        contents = file.read()

    individuals = contents.split("\n")[:-1]
    # individuals = ['individual1', 'individual2', 'individual3', 'individual4']
    pairs = []
    for i in range(len(individuals)):
        for j in range(i+1, len(individuals)):
            pairs.append((individuals[i], individuals[j]))

    # create the process pool
    with Pool() as pool:
        # # call the same function with different data in parallel
        for result in pool.imap(compare_individuals, pairs):
        #     # report the value to show progress
            outfile.write(result)
