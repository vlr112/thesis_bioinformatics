"""adaptation of genotype_acc.py. compares to vcf files at same time, only on SNPs both share. return:
 
Find out with lcmlkin how many times that lcmlkin, 1) does not call a genotype 2) call genotypes correctly 3) calls genotypes incorrectly
Find out with bcftools how many times bcftools, 1) does not call a genotype 2) calls genotypes correctly 3) call genotypes incorrectly
Now compare the genotypecalls between lcmlkind and bcftools. 1) How many times is either of the called genotypes missing 2) how many times do they call the same genotype 3) how many times do they call different genotypes.
"""

# cmd line: ./genotype_acc.changed.py allbams_mp.vcf allbams_final.vcf.gz

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
import vcf
import gzip
import subprocess
from prettytable import PrettyTable
from tabulate import tabulate


# from cyvcf2 import VCF


# Get input arguments. Note: this files are the header free files
vcf_bcftools = sys.argv[1] # VCF file path
vcf_lcmlkin= sys.argv[2] # VCF file path
# common_list=sys.argv[3] # vcf positions shared by both vcf files
fasta_list_file = sys.argv[3] # Text file path containing list of fasta files
outfile__path=sys.argv[4]



# Read list of fasta files
with open(fasta_list_file, 'r') as f:
    fasta_files = f.read().splitlines()

def call_genotypes(sample, genotype, pos, ref_allele, alt_allele):

    correct = 0 
    missed = 0
    called = 0 
    called_allele1= 0
    called_allele2 = 0
    fasta1 = Fasta(fasta_files[2*sample])
    fasta2 = Fasta(fasta_files[2*sample+1])
    # name of fasta sequence header
    key1 = (list(fasta1.keys())[0])
    key2 = (list(fasta2.keys())[0])

    # becasue positions are the same in both vcf files, it doesn't matter from which i get the position value to fetch the real fasta base.         
    allele1_true = fasta1[key1][pos -1]
    allele2_true = fasta2[key2][pos -1]

    allele1, allele2 = genotype[sample].split('/')  # 0/0, 0/1, 1/1

    # homozygous reference for bcftools vcf
    if (allele1 == '0' and allele2 == '0'):
        called_allele1 = ref_allele
        called_allele2 = ref_allele
        if (allele1_true == ref_allele) and (allele2_true == ref_allele):
            correct += 1
        called += 1

    # homozygous alternative
    elif allele1 == '1' and allele2 == '1':
        called_allele1 = alt_allele
        called_allele2 = alt_allele
        if (allele1_true != ref_allele) and (allele2_true != ref_allele):
            correct += 1
        called += 1

    # heterozygous
    elif (allele1 == '1' and allele2 == '0') or (allele1 == '0' and allele2 == '1'):
        called_allele1 = ref_allele
        called_allele2 = alt_allele
        if ((allele1_true != ref_allele) and (allele2_true == ref_allele)) or ((allele1_true == ref_allele) and (allele2_true != ref_allele)):
            correct += 1
        called += 1

    # miss call genotype
    elif allele1 == '.' or allele2 == '.':
        # total_called = 1
        missed += 1
        called += 1

    return called_allele1, called_allele2, correct, called, missed


def task(g):
    """ Only SNP positions on the common list are considered"""

    # print('g', g)

    # Initialize counters
    total_correct_b = 0
    total_correct_l = 0

    total_missed_b = 0
    total_missed_l = 0

    total_called_b = 0
    total_called_l = 0

    total_b = 0
    total_l = 0

    same_geno_called=0  # 2) how many times do they call the same genotype 
    diff_geno_called=0  # 2) how many times do they call the different genotype 
    same_correct_geno_called=0
    same_wrong_geno_called=0

    
    fields_bcftools = bcftools[g].split('\t')
    pos_bcftools = int(fields_bcftools[1])
    ref_allele_bcftools = fields_bcftools[3]
    alt_allele_bcftools = fields_bcftools[4]

    fields_lcmlkin = lcmlkin[g].split('\t')
    pos_lcmlkin = int(fields_lcmlkin[1])
    ref_allele_lcmlkin = fields_lcmlkin[3]
    alt_allele_lcmlkin = fields_lcmlkin[4]


    genotypes_b = [x.split(':')[0] for x in fields_bcftools[9:]] # Extract genotype calls
    genotypes_l = [x.split(':')[0] for x in fields_lcmlkin[9:]] # Extract genotype calls

    assert pos_bcftools == pos_lcmlkin
    assert len(genotypes_b) == len(genotypes_l)

    # Iterate over each haplotype fasta file and check genotypes. check 2 haplotypes at a time. so same individual
    for i in range(0, len(genotypes_b)):

        called_allele1_b, called_allele2_b, correct_b, called_b, missed_b = call_genotypes(i, genotypes_b, pos_bcftools, ref_allele_bcftools, alt_allele_bcftools )
        called_allele1_l, called_allele2_l, correct_l, called_l, missed_l = call_genotypes(i, genotypes_l, pos_lcmlkin, ref_allele_lcmlkin, alt_allele_lcmlkin )

        ## check if they called same genotype
        if ((called_allele1_b == called_allele1_l) and (called_allele2_b == called_allele2_l)) or ((called_allele1_b == called_allele2_l) and (called_allele2_b == called_allele1_l)):
            same_geno_called +=1
            if correct_b==correct_l==1:
                same_correct_geno_called += 1
            else:
                same_wrong_geno_called+=1

        # they each call different genotypes
        elif (missed_b==missed_l==0) :
            diff_geno_called += 1


        # check if individually they correctly called genotypes or not.
        total_called_b += called_b
        total_called_l += called_l

        total_correct_b += correct_b
        total_correct_l += correct_l

        total_missed_b += missed_b
        total_missed_l += missed_l     

    total_b = total_missed_b +  total_called_b
    total_l = total_missed_l +  total_called_l


    return total_called_b, total_called_l, total_correct_b, total_correct_l, total_missed_b, total_missed_l, same_correct_geno_called, same_wrong_geno_called, diff_geno_called

# total_missed_b, total_called_b



if __name__ == '__main__':

    start_time = time.perf_counter()

    lcmlkin = open(vcf_lcmlkin, 'r').readlines()
    bcftools = open(vcf_bcftools, 'r').readlines()

    mask_total_called_b = []
    mask_total_called_l =[]
    mask_total_correct_b = []
    mask_total_correct_l= [] 
    mask_total_missed_b = []
    mask_total_missed_l = []
    mask_same_correct_geno_called = []
    mask_same_wrong_geno_called = []
    mask_diff_geno_called = []
    # total_correct=0
    # mask_correct=[]
    # mask_called=[]

    # # # total_b, total_l= task(1)
    # # # print(total_b)
    # create the process pool
    with Pool() as pool:

        # # call the same function with different data in parallel
        for total_called_b, total_called_l, total_correct_b, total_correct_l, total_missed_b, total_missed_l, same_correct_geno_called, same_wrong_geno_called, diff_geno_called in pool.imap(task, range(len(lcmlkin)), chunksize=100):
            mask_total_called_b.append(total_called_b)
            mask_total_called_l.append(total_called_l)
            mask_total_correct_b.append(total_correct_b)
            mask_total_correct_l.append(total_correct_l)
            mask_total_missed_b.append(total_missed_b)
            mask_total_missed_l.append(total_missed_l)
            mask_same_correct_geno_called.append(same_correct_geno_called)
            mask_same_wrong_geno_called.append(same_wrong_geno_called)
            mask_diff_geno_called.append(diff_geno_called)

    correct_call_b = sum(mask_total_correct_b)
    correct_call_l = sum(mask_total_correct_l)

    total_call_b = sum(mask_total_called_b)
    total_call_l = sum(mask_total_called_l)

    total_missed_b = sum(mask_total_missed_b)
    total_missed_l = sum(mask_total_missed_l)

    total_b= total_call_b + total_missed_b
    total_l= total_call_l + total_missed_l

    
    # table = [[' ','MISSED CALL', 'CORRECT CALL', 'WRONG CALL'], ['bcftools',f"{total_missed_b} , ({total_missed_b/total_b*100:.3f} %) ",  f"{correct_call_b} , ({correct_call_b/total_b*100:.3f} %) ", total_call_b - correct_call_b ], ['lcmlkin',total_missed_l, correct_call_l, total_l - correct_call_l ]]
    
    table = [[' ','MISSED CALL', 'CORRECT CALL', 'WRONG CALL'], ['bcftools', f"{total_missed_b} , ({total_missed_b/total_b*100:.3f} %) ",  f"{correct_call_b} , ({correct_call_b/total_b*100:.3f} %) ", f"{total_call_b - correct_call_b} , ({(total_call_b - correct_call_b)/total_b*100:.3f} %) " ], ['lcmlkin' , f"{total_missed_l} , ({total_missed_l/total_l*100:.3f} %) ",  f"{correct_call_l} , ({correct_call_l/total_l*100:.3f} %) ", f"{total_call_l - correct_call_l} , ({(total_call_l - correct_call_l)/total_l*100:.3f} %) " ]]


    print(
    """
    MISSED CALL: Number of genotypes that are ./.
    CORRECT CALL: called genotype corresponds to the truth  
    WRONG CALL: called genotype does not correspond to the truth

    """
          )

    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
    print('Number of times the same genotype was correctly called: ', sum(mask_same_correct_geno_called))
    print('Number of times the same genotype was wrongly called: ', sum(mask_same_wrong_geno_called))
    print('Number of times each program called a different genotype: ', sum(mask_diff_geno_called))
    print('Total calls + missed: ', total_b)

    with open(outfile__path + f"missed.bcftools.joint", 'w') as f:
        f.write(str(total_missed_b))
        f.close()

    with open(outfile__path + f"correct.bcftools.joint", 'w') as f:
        f.write(str(correct_call_b))
        f.close()

    with open(outfile__path + f"wrong.bcftools.joint", 'w') as f:
        f.write(str(total_call_b - correct_call_b))    
        f.close()

    with open(outfile__path + f"missed.lcmlkin.joint", 'w') as f:
        f.write(str(total_missed_l))
        f.close()

    with open(outfile__path + f"correct.lcmlkin.joint", 'w') as f:
        f.write(str(correct_call_l))
        f.close()

    with open(outfile__path + f"wrong.lcmlkin.joint", 'w') as f:
        f.write(str(total_call_l - correct_call_l))    
        f.close()
    # tab = PrettyTable(table[0])
    # print(tab)


    #     # print('final total correct', total_correct)
    # print('mask_correct', mask_correct)
    # print('mask_called', mask_called)
    # # percent_correct = 100.0 * sum(mask_correct) / sum(mask_called)
    # print("Percentage of correctly called genotypes over all individuals and all SNPs: {:.2f}%".format(percent_correct))

    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
 
# print('Finished. Laterz')



# if __name__ == '__main__':

#     start_time = time.perf_counter()

#     lcmlkin = open(vcf_lcmlkin, 'r').readlines()
#     bcftools = open(vcf_bcftools, 'r').readlines()

#     total_correct=0
#     mask_correct=[]
#     mask_called=[]
#     # create the process pool
#     with Pool() as pool:

#         # # call the same function with different data in parallel
#         for total_correct, total_called in pool.imap(task, range(len(f)), chunksize=100):
#             mask_correct.append(total_correct)
#             mask_called.append(total_called)

#         # print('final total correct', total_correct)
#     # print('mask_correct', mask_correct)
#     # print('mask_called', mask_called)
#     percent_correct = 100.0 * sum(mask_correct) / sum(mask_called)
#     print("Percentage of correctly called genotypes over all individuals and all SNPs: {:.2f}%".format(percent_correct))

#     finish_time = time.perf_counter()
#     print(f"Program finished in {finish_time-start_time} seconds")
 
# print('Finished. Laterz')












    # # heterozygous
        # elif (allele1 == '1' and allele2 == '0'): 
        #     if ((allele1_true != ref_allele) and (allele2_true == ref_allele)):
        #         total_correct += 1
        #     total_called += 1
        
        # #heterozygous
        # elif (allele1 == '0' and allele2 == '1'):
        #     if ((allele1_true == ref_allele) and (allele2_true != ref_allele)):
        #         total_correct += 1
        #     total_called += 1

########################################################################


    #     fasta1 = Fasta(fasta_files[2*i])
    #     fasta2 = Fasta(fasta_files[2*i+1])
    #     # name of fasta sequence header
    #     key1 = (list(fasta1.keys())[0])
    #     key2 = (list(fasta2.keys())[0])


    #     # becasue positions are the same in both vcf files, it doesn't matter from which i get the position value to fetch the real fasta base.         
    #     allele1_true = fasta1[key1][pos_bcftools -1]
    #     allele2_true = fasta2[key2][pos_bcftools -1]

    #     # allele1_true_l = fasta1[key1][pos_lcmlkin -1]
    #     # allele2_true_l = fasta2[key2][pos_lcmlkin -1]

    #     if i==3:
    #         print('allele1_true', allele1_true)
    #         print('allele2_true', allele2_true)

    #     allele1_b, allele2_b = genotypes_b[i].split('/')  # 0/0, 0/1, 1/1
    #     allele1_l, allele2_l = genotypes_l[i].split('/')


    #     # homozygous reference for bcftools vcf
    #     if allele1_b == '0' and allele2_b == '0':
    #         called_allele1_b = ref_allele_bcftools
    #         called_allele2_b = ref_allele_bcftools
    #         if (allele1_true == ref_allele_bcftools) and (allele2_true == ref_allele_bcftools):
    #             total_correct += 1
    #             print('allele1_true', allele1_true)
    #             print('allele2_true', allele2_true)
    #         total_called += 1

    #     # homozygous alternative
    #     elif allele1 == '1' and allele2 == '1':
    #         if (allele1_true != ref_allele) and (allele2_true != ref_allele):
    #             # print('allele1_true', allele1_true)
    #             # print('allele2_true', allele2_true)
    #             total_correct += 1
    #         total_called += 1

    #     # heterozygous
    #     elif (allele1 == '1' and allele2 == '0') or (allele1 == '0' and allele2 == '1'):
    #         if ((allele1_true != ref_allele) and (allele2_true == ref_allele)) or ((allele1_true == ref_allele) and (allele2_true != ref_allele)):
    #             total_correct += 1
    #             # print('allele1_true', allele1_true)
    #             # print('allele2_true', allele2_true)
    #         total_called += 1

    #  # not called
    #     elif allele1 == '.' and allele2 == '.':
    #         total_called += 1
    #         # count += 1

    #     else:
    #         total_called += 1
    #         count += 1
    #         print('oar bolas!!!!')
    # print(count)
    # print(total_called)

    # return total_correct, total_called