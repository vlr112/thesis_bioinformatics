import sys
from pyfaidx import Fasta
###import libraries
import string
import numpy as np
import pysam
import gzip
import math
import copy
from sys import argv
import time
from multiprocessing import Pool 
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import multiprocessing
from pysam import VariantFile
from tabulate import tabulate


# Get input arguments
vcf_file = sys.argv[1] # VCF file path
fasta_list_file = sys.argv[2] # Text file path containing list of fasta files
outfile_suffix = sys.argv[3]
outfile__path=sys.argv[4]

# Read list of fasta files
with open(fasta_list_file, 'r') as f:
    fasta_files = f.read().splitlines()

def task(g):
    # Initialize counters
    total_correct = 0
    total_called = 0
    total_missed = 0

    fields = f[g].split('\t')
    pos = int(fields[1])
    ref_allele = fields[3]
    # alt_allele = fields[4]
    # print(alt_allele)
    genotypes = [x.split(':')[0] for x in fields[9:]] # Extract genotype calls
    # print(genotypes)
    count=0
    # print('len(genotypes)', len(genotypes))

    # Iterate over each haplotype fasta file and check genotypes. check 2 haplotypes at a time. so same individual
    for i in range(0, len(genotypes)):
        # print(i)
        fasta1 = Fasta(fasta_files[2*i])
        fasta2 = Fasta(fasta_files[2*i+1])
        # name of fasta sequence header
        key1 = (list(fasta1.keys())[0])
        key2 = (list(fasta2.keys())[0])

        allele1_true = fasta1[key1][pos -1]
        allele2_true = fasta2[key2][pos -1]

        allele1, allele2 = genotypes[i].split('/')

        # not called
        if allele1 == '.' or allele2 == '.':
            total_missed += 1

        # homozygous reference
        elif allele1 == '0' and allele2 == '0':
            if (allele1_true == ref_allele) and (allele2_true == ref_allele):
                total_correct += 1
            total_called += 1

        # homozygous alternative
        elif allele1 == '1' and allele2 == '1':
            if (allele1_true != ref_allele) and (allele2_true != ref_allele):
                total_correct += 1
            total_called += 1

        # heterozygous
        elif (allele1 == '1' and allele2 == '0'): 
            if ((allele1_true != ref_allele) and (allele2_true == ref_allele)):
                total_correct += 1
            total_called += 1
        
        #heterozygous
        elif (allele1 == '0' and allele2 == '1'):
            if ((allele1_true == ref_allele) and (allele2_true != ref_allele)):
                total_correct += 1
            total_called += 1

    return total_correct, total_called, total_missed


if __name__ == '__main__':

    start_time = time.perf_counter()

    f = open(vcf_file, 'r').readlines()

    total_correct=0
    mask_correct=[]
    mask_called=[]
    mask_missed=[]
    # create the process pool
    with Pool() as pool:

        # # call the same function with different data in parallel
        for total_correct, total_called, total_missed in pool.imap(task, range(len(f)), chunksize=100):
            mask_correct.append(total_correct)
            mask_called.append(total_called)
            mask_missed.append(total_missed)

        # print('final total correct', total_correct)
    # print('mask_correct', mask_correct)
    # print('mask_called', mask_called)
    total= sum(mask_called) + sum(mask_missed)
    total_missed = sum(mask_missed)
    total_called = sum(mask_called)
    total_correct = sum(mask_correct)

    table = [[' ','MISSED CALL', 'CORRECT CALL', 'WRONG CALL'], [outfile_suffix, f"{total_missed} , ({total_missed/total*100:.3f} %) ",  f"{total_correct} , ({total_correct/total*100:.3f} %) ", f"{total - total_correct} , ({(total - total_correct)/total*100:.3f} %) " ]]

    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))

    with open(outfile__path + f"missed.{outfile_suffix}.individual", 'w') as f:
        f.write(str(total_missed))
        f.close()

    with open(outfile__path + f"correct.{outfile_suffix}.individual", 'w') as f:
        f.write(str(total_correct))
        f.close()

    with open(outfile__path + f"wrong.{outfile_suffix}.individual", 'w') as f:
        f.write(str(total - total_correct))    
        f.close()

    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")

    with open(outfile__path + f"{outfile_suffix}.individual.TIME", 'w') as f:
        f.write(str(finish_time-start_time)) 
        f.close()  
 
# print('Finished. Laterz')


