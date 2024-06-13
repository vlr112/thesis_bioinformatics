#!/usr/bin/python

## sctipt's goal: calculate IBD probabilities based on IBS counts. 
## for this, count ibs states in pairwise comparison of fasta files between samples
## then use bias adjustment from plink to obtain the conditional (E00, E10 and so on)
## question: not really sure, this adjustment is necessary, since we're using 
## the true allele frequency. however, without it the results are nonsense

import numpy as np
import random
import time
import glob
import sys
from tabulate import tabulate
from multiprocessing import Pool
import pandas as pd

# last_time = time.time()


def count():
    """The TrueIBD has the IBS counts. Now we need to get the IBD probabilities.
    We can get that by following what is done in PLINK, minus the adjustment for the bias 
    (the bias is considering a finite size population we calculating the AF. In this case,
    we have the true AF values)  """

    af=pd.DataFrame(pd.read_csv('af.txt', delimiter="\t", header=None))
    af = af.rename(columns={0: 'freq'})
    af['freq'] = af['freq'].apply(lambda x: x if x >= 0.5 else 1 - x )

    E00=E10=E20=E01=E11=E21=E02=E12=E22=0
    cnt=0

    for t in range(len(af['freq'])):

        p = af['freq'][t]
        q = 1 - af['freq'][t]
        Na = 78
        x = p*Na
        y = q*Na

        if q != 0:

            a00 = 2 * p * p * q * q * ((x - 1) / x * (y - 1) / y * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3)))

            a01 = 4 * p * p * p * q * ((x - 1) / x * (x - 2) / x * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3))) \
                + 4 * p * q * q * q * ((y - 1) / y * (y - 2) / y * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3)))

            a02 = q * q * q * q * ((y - 1) / y * (y - 2) / y * (y - 3) / y * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3))) \
                + p * p * p * p * ((x - 1) / x * (x - 2) / x * (x - 3) / x * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3))) \
                + 4 * p * p * q * q * ((x - 1) / x * (y - 1) / y * (Na / (Na - 1)) * (Na / (Na - 2)) * (Na / (Na - 3)))

            a11 = 2 * p * p * q * ((x - 1) / x * Na / (Na - 1) * Na / (Na - 2)) \
                + 2 * p * q * q * ((y - 1) / y * Na / (Na - 1) * Na / (Na - 2))

            a12 = p * p * p * ((x - 1) / x * (x - 2) / x * Na / (Na - 1) * Na / (Na - 2)) \
                + q * q * q * ((y - 1) / y * (y - 2) / y * Na / (Na - 1) * Na / (Na - 2)) \
                + p * p * q * ((x - 1) / x * Na / (Na - 1) * Na / (Na - 2)) \
                + p * q * q * ((y - 1) / y * Na / (Na - 1) * Na / (Na - 2))
        # if q != 0:
        #     a00 = 2 * p * p * q * q 

        #     a01 = 4 * p * p * p * q  + 4 * p * q * q * q 

        #     a02 = q * q * q * q \
        #         + p * p * p * p * \
        #         + 4 * p * p * q * q 

        #     a11 = 2 * p * p * q \
        #         + 2 * p * q * q 

        #     a12 = p * p * p \
        #         + q * q * q \
        #         + p * p * q \
        #         + p * q * q 

            E00 += a00
            E01 += a01
            E02 += a02
            E11 += a11
            E12 += a12
            cnt+=1

    E00 = E00/cnt
    E01 = E01/cnt
    E02 = E02/cnt
    E11 = E11/cnt
    E12 = E12/cnt

    print('count', cnt)

    return E00, E01, E02, E11, E12

E00, E01, E02, E11, E12 = count()


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

    ibs0=0
    ibs1=0
    ibs2=0

    i=0
    for v1_h1, v1_h2, v2_h1, v2_h2, a in zip(v1_h1_geno, v1_h2_geno,  v2_h1_geno, v2_h2_geno, ancestral):

        g1 = [v1_h1, v1_h2]
        g2 = [v2_h1, v2_h2]

        g1.sort()
        g2.sort()
        if g1 == g2:
            ibs2+=1
        elif g1[0] == g2[0] or g1[1] == g2[1]:
            ibs1 +=1
        else:
            ibs0+=1

    # total allele count count (should be 100 000 * 2)
    S = ibs0 + ibs1 + ibs2

    e00 = E00*S 
    e01 = E01*S 
    e02 = E02*S  
    e11 = E11*S 
    e12 = E12*S 
    e22 = 1*S

    z0 =  ibs0 / e00
    z1 = (ibs1 - z0*e01) / e11
    z2 = (ibs2 - z0*e02 - z1*e12) / e22
    kin = 0.5*z2 + 0.25*z1

    out = str(ind1) + '\t' + str(ind2) + '\t' + str(z0) + '\t' + str(z1) + '\t' + str(z2) + '\t' + str(kin)
    out = out + '\n'
    return out

if __name__ == '__main__':

    names = sys.argv[1]
    filenameout = sys.argv[2]
    outfile=open(filenameout,'w')
    out='a\tb\tJ9\tJ8\tJ7\tKinship'
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
