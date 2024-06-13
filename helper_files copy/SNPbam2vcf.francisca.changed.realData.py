#!/usr/bin/env python
# -*- coding: ASCII -*-
 
#####This python script takes a list of bamfiles and performs variant calling such that the output is in the format required for lcMLkin to calculate pairwise relatedness.
#####Please note it is slow, and is only designed for calling genotypes for specific SNPs in smallish numbers (tens to hundreds) of individuals. I do not recommend using this as a general variant caller.
#####We are working on a better version (this was an experimental one developed to test lcMLkin), but for now, if using 100 samples with 100K loci with 20x mean coverage takes about 4 hours on our machine. Lower coverage datasets, which is what lcMLkin is designed for, will be quicker though, maybe just over an hour.
#####Pysam and Numpy must be installed in the version of python used.
#####Bam files must be indexed.
#####The SNP file must have the tab seperated fields in the following order: chromosome, position (one-based), reference allele, alternate allele
#####If you want to change things like mapping and base quality threshold, edit the python code under the section "Input arguments"
#####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)

##### New:multithreading/multiproccessing in the script. As such, some functions were adapted.  
#####  Changed so it is recognized by bcftools and can operate with it
##### Change2: allow for beyesian caller to run. 
    # before it would skip individuals whith genotype below quality threshold. what was happening was
    # at some SNPs (the ones with at least one individual in condition mentioned above),
    # the len(GLs) would be smaller than len samples. And that's a problem! Anyway, solved now :)
 
 
#####usage is ./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file> 
 
 
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

 
###Input arguments
BAMin=argv[1] #flat text file consisiting of list of bams (must be indexed), one per line
filenameout=argv[2] #creates a vcf
SNPfile=argv[3] #must have the tab seperated fields in the following order chromosome, position (one-based), reference allele, alternate allele
BC=1 #choose whether to use the bayesian caller for individual genotypes (0=no, 1=yes). I don't recommend turning this on the current  version, as it's quite slow, and allele frequencies shouldn't really be biased if the genotypes are incorreclty called, especially as the genotype calling we use is not reference aware
MQ_t=20 #mapping quality threshold
BQ_t=5 #base_qualitythreshold
GQ_t=0.1 #GQ threshold
min_sams=1 #minimum samples with reads at locus needed to attempt bayesian calling of genotypes
exc_fix=0  #exclude fixed differences (doesn't make much sense to turn on if only looking at one sample
 
###converts phred score to probability
def phred2prob(x):
    return 10.0**(-x/10.0)
 
###converts probability to phred score
def prob2phred(x):
    return -10*math.log10(x)
 
###diploid caller assuming all alternative alleles are possible (error is divided by three for now, could add a more complex model)
def geno_caller_10GT(X):
 
    GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T
 
    all_dic={}
    all_dic['A']=0
    all_dic['C']=1
    all_dic['G']=2
    all_dic['T']=3
 
    count=0
    for g in range(len(X)):
        if X[g][0] not in all_dic:
            continue
        hap[g]=phred2prob(X[g][1])*(1.0/3.0)
        hap[g][all_dic[X[g][0]]]=1-phred2prob(X[g][1])
 
        GL[0]=GL[0]+math.log10(hap[g][0])
        GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)
        GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
        GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)
 
        GL[4]=GL[4]+math.log10(hap[g][1])
        GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
        GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)
 
        GL[7]=GL[7]+math.log10(hap[g][2])
        GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)
 
        GL[9]=GL[9]+math.log10(hap[g][3])
        count+=1
 
    if count==0:
        GL.fill(-9.0)
    return GL
 
all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3
 
def geno_caller_3GT(X,ref,alt,all_dic):
    #diploid caller assuming that only assesses likelihood for three possible genotypes (ref/ref,ref/alt,alt/alt)
    GL=[0.0,0.0,0.0]
 
    count=0
    for g in range(len(X)):
        if X[g][0] not in all_dic:
            # print('X[g][0]', X[g][0])
            continue
        err=phred2prob(X[g][1])*(1.0/3.0)
        tru=1-phred2prob(X[g][1])
 
        if X[g][0]==ref:
            GL[0]=GL[0]+math.log10(tru)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(err)
        elif X[g][0]==alt:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10((tru+err)/2)
            GL[2]=GL[2]+math.log10(tru)
        else:
            GL[0]=GL[0]+math.log10(err)
            GL[1]=GL[1]+math.log10(err)
            GL[2]=GL[2]+math.log10(err)
        count+=1
    
    if count==0:
        GL=[-9.0,-9.0,-9.0]
        # print('count',count)

    # print('len(GL)',GL)
    return GL
 
 
###extract reads for a give position in a bam
def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
    var_list=[]
    for pileupcolumn in samfile.pileup(chromo,pos-1,pos,truncate=True,stepper='all'):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):
                        var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
    return var_list
 
#####extract reads for a give position in a bam (original version of above,not clear which is faster)
##def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
##    var_list=[]
##    for pileupcolumn in samfile.pileup(chromo,pos-1,pos):
##        if pileupcolumn.pos==pos-1:
##            for pileupread in pileupcolumn.pileups:
##                if not pileupread.is_del and not pileupread.is_refskip:
##                    if pileupread.alignment.mapping_quality>=MQ:
##                        if pileupread.alignment.is_duplicate==False:
##                            if ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ:
##                                var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
##    return var_list
 
 
###Version of Depristo's bayesian caller
def bayes_caller(X):
    het=0.001   # heterozygosity starting point
    nb_ind=len(X) ## number of individuals
 
    # print('len(X)',len(X))
    # print('X', X)
    thetas=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix with theta priors
    # print(thetas)
    for g in range(1,len(thetas)):
        thetas[g]=het/float(g)
        # print(thetas)
    thetas[0]=1-sum(thetas[1:])
    # print(thetas)
 
    best_genos=np.zeros(((nb_ind*2)+1,nb_ind),dtype='int64')  #matrix for storing best genotypes for a given allele frquency
    # print(best_genos)
    best_genos_ll=np.zeros(((nb_ind*2)+1,nb_ind),dtype='float32')  #matrix for storing best genotypes for a given allele frquency
    # print(best_genos_ll)

    p_best_genos=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix for storing the likelihood for each p
 
    for g in range(len(X)):   ### p = 0
        # print(best_genos_ll[0][g])
        best_genos_ll[0][g]=best_genos_ll[0][g]+X[g][0] # for individual g get GT of AA
        # print(X[g][0])
        # print(best_genos_ll)
    p_best_genos[0]=np.sum(best_genos_ll[0]) ## sum of genotype likelihoods of all individuals with GT=0/0
    # print(p_best_genos)
 
    ##perform greedy algorithm to find best likelihoods for each possible configuration as described in Depristo
    for g in range(1,len(p_best_genos)):
    # for g in range(1,3):

        p_best_genos[g]=-1000000000000
        # print(p_best_genos[g])
        genos=np.copy(best_genos[g-1]) # zeros 1D matrix (1,number ind.)
        # print(genos)
        # print('len(genos)',len(genos))
        genos_ll=np.copy(best_genos_ll[g-1])
        # print('genos_ll', genos_ll)
        for gg in range(len(genos)): # iterate each individual 
            if genos[gg]<2: # 2 because its biallelic. 0: homo ref, 1:het, 2:homo alt
                # print('genos[gg]', genos[gg])

                genos_mod=np.copy(genos)
                genos_mod[gg]+=1
                # print('genos_mod', genos_mod)
                genos_ll_mod=np.copy(genos_ll)
                genos_ll_mod[gg]=X[gg][genos[gg]+1] # genos[gg] = 0 (always, because matrix is full with zeros). X[gg] = [AA, AT, TT] for individual X. So X[gg][genos[gg]+1]= X[gg][1] = AA for ind. gg
                # print('[genos[gg]+1]', [genos[gg]+1])
                # print('X[gg]',X[gg])
                # print('genos_ll_mod', genos_ll_mod)
                ll=np.sum(genos_ll_mod)
                # print('ll', ll)
                # print('p_best_genos[g]', p_best_genos[g])
                if ll>p_best_genos[g]:
                    p_best_genos[g]=ll
                    # print(p_best_genos)
                    # print('(p_best_genos[g]', p_best_genos[g])
                    best_genos[g]=np.copy(genos_mod)
                    best_genos_ll[g]=np.copy(genos_ll_mod)

        # print('###########')
        # print('')
        
        # # print('X[1]',X[1])
        # print('genos_mod_end', genos_mod)
    # print('best_genos_ll[g]', best_genos_ll[g])
    # print('p_best_genos', len(p_best_genos))
    # print(X)
    # print(X[0])
    # print('best_genos', best_genos)
 
    ##Calculate variant quality
    #denom=np.max(p_best_genos)+np.log10(np.sum(10**(p_best_genos-np.max(p_best_genos))))
    #phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)
    marg_comp=np.log10(thetas)+p_best_genos
    # print('thetas', thetas)
    # print('marg_comp', marg_comp)
    denom=np.max(marg_comp)+np.log10(np.sum(10**(marg_comp-np.max(marg_comp))))
    phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)
 
    VQ=int(round(phred))
 
    geno_comb=best_genos[np.argmax(p_best_genos)] # geno_comb=[1 1 1 1 1 0 1 1 1 1 1 1 1 1]
    # print('geno_comb', len(geno_comb))
 
    if np.sum(geno_comb)==0:
        stat='REF'
 
    else:
        stat='ALT'
    return geno_comb,VQ,stat
 
###dictionary for GL order
geno_ord={}
geno_ord['AA']=0
geno_ord['AC']=1
geno_ord['AG']=2
geno_ord['AT']=3
geno_ord['CA']=1
geno_ord['CC']=4
geno_ord['CG']=5
geno_ord['CT']=6
geno_ord['GA']=2
geno_ord['GC']=5
geno_ord['GG']=7
geno_ord['GT']=8
geno_ord['TA']=3
geno_ord['TC']=6
geno_ord['TG']=8
geno_ord['TT']=9
 
###dictionary for vcf genotype codes
vcf_gt={}
vcf_gt[0]='0/0'
vcf_gt[1]='0/1'
vcf_gt[2]='1/1'
 
 
###Read bamlist
file=open(BAMin,'r')
samples=file.read()
samples=str.split(samples,'\n')
if samples[-1]=='':
    del(samples[-1])

# Open file with samples names
###Get sample name for each bam. Each bam must only contain one sample
child_names = []
with open(BAMin, "r") as f:
    for line in f:
        file_path = line.strip()
        child_name = os.path.basename(file_path).split(".")[0]
        child_names.append(child_name)
        
 
###Read SNPlist
file=open(SNPfile,'r')
SNPs=file.read()
SNPs=str.split(SNPs,'\n')
if SNPs[-1]=='':
    del(SNPs[-1])
 
# sample_name=[]
# samfiles={}
# ###Get sample name for each bam. Each bam must only contain one sample, otherwise the program will take the sample in first readgroup in the header
# for g in range(len(samples)):
#     samfiles[g] = pysam.AlignmentFile(samples[g], "rb")
#     sample_name.append(samfiles[g])#.header['RG'][0]['SM'])
 
 
###start writing outfile with header
outfile=open(filenameout,'w')
out="""##fileformat=VCFv4.2
##bcftoolsVersion=1.9+htslib-1.9
##FILTER=<ID=.,Description="Missing FILTER">
##contig=<ID=1,length=100100000>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=AF,Number=A,Type=Float,Description="allele frequencies (from multisampling or single-base calling)">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="genotype likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""
# out='##Custom variant caller for_lcMLkin\n##Bayesian_caller='+str(BC)+',MappingQuality_filter='+str(MQ_t)+',BaseQuality_filter='+str(BQ_t)+',GenotypeQuality_filter='+str(GQ_t)+',bamlist='+BAMin+',filenameout='+filenameout+',SNPfile='+SNPfile+'\n'
out=out+'##Time created: '+(time.strftime("%H:%M:%S"))+' '+(time.strftime("%d/%m/%Y"))+'\n'
out=out+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
for x in child_names:
    out=out+'\t'+ str(x)
 
out=out+'\n'
 
outfile.write(out)

def task(g):
 
    ###Iterate through SNPs
    # for g in range(len(SNPs)):
    k=str.split(SNPs[g])
    # return k
    chromo=k[0]
    pos=int(k[1])
    ref=k[2]
    alt=k[3]

    ###iterate through samples
    GLs=[]
    GLs_nomiss=[]
    GLs_miss=[]
    GQs=[]
    GTs=[]

    for gg in range(len(samples)):

        samfile = pysam.AlignmentFile(samples[gg], "rb")
        # print(samples[gg])
        ###extract reads for a position
        read_list=extract_bam_SNP(samfile,chromo,pos,BQ_t,MQ_t)

        # print('read_list', len(read_list))

        ###get the genotype likelihood for the 3 possible genotypes based on the ref and alt alleles provided
        GL=geno_caller_3GT(read_list,ref,alt,all_dic)

        ###calcule genotype quality
        GL_sort=copy.deepcopy(GL)
        GL_sort.sort()
        GQ = -10 * (GL_sort[1] - GL_sort[2])
        GQs.append(GQ)

        ###assign genotype using highest likelihood
        if BC==1:
            if GQ>GQ_t:
                GLs_nomiss.append(GL)
            else:
                # GLs_nomiss.append([-9.0,-9.0,-9.0])
                # print('GQ',GQ)
                # print('GQ_t',GQ_t)
                # print(samples[gg])
                # print(gg)
                # print(GL)
                GLs_miss.append(gg) # save line number of sample that will not be included in bayesian caller.

        elif BC==0:
            if GQ>GQ_t:
                GTs.append(vcf_gt[np.argmax(np.asarray(GL))])
            else:
                GTs.append('./.')
                # GLs
        GLs.append(GL)
    # print('GLs', len(GLs))
    skip=0

    ###estimate allele frequencies for multisample calling
    if BC==1:
        if len(GLs_nomiss)>min_sams:  #must have genotype data for at least 20 individuals to consider SNP
            genos_b,VQ,s=bayes_caller(GLs_nomiss)
            # print('GLs_nomiss', GLs_nomiss)
            # print('len(GLs_nomiss)', len(GLs_nomiss))
            # print(GQ_t)
            # print('GLs', GLs)

            AF=(np.sum(genos_b)/(float(len(genos_b))*2))
            it=0
            for ggg in range(0,len(GLs)):
                # print(ggg)
                # for excluded in range(len())
                if -9.0 in GLs[ggg]:
                    GTs.append('./.')

                elif GLs_miss.count(ggg): # basically if GQ fails quality then the genotype IS NOT called
                    GTs.append('./.')
                    # print(GTs)
                    # print(GLs[ggg])

                else:
                    GTs.append(vcf_gt[genos_b[it]])
                    # print(GTs)
                    # print(len(GTs))
                    it+=1

            # print('bang',len(GTs))
                

        # else:
            # skip=1
            # GTs.append('./.')

    ###estimate allele frequencies for single-base calling
    elif BC==0:
        if (float((len(GTs)-GTs.count('./.')))*2) != 0.0:
            AF=1-((GTs.count('0/1')+GTs.count('1/1')*2) / (float((len(GTs)-GTs.count('./.')))*2))
            VQ=100
            s='.'
        else:
            skip=1

    ###only write to file if site is variable (otherwise useless for kinship estimation)
    if skip==0:
        out=k[0]+'\t'+k[1]+'\t.\t'+k[2]+'\t'+k[3]+'\t'+str(VQ)+'\t'+s+'\tAF='+str(AF)+'\tGT:GL:GQ'
        if exc_fix==0: ###Don't throw out sites that are fixed. If examining only a few or one individual, and population allele frequencies in lcMLkin come from elsewhere, this setting makes sense to use.
            for ggg in range(len(GTs)):
                if GTs[ggg]=='./.':
                    out=out+'\t'+GTs[ggg]+':'+str(GLs[ggg][0])+','+str(GLs[ggg][1])+','+str(GLs[ggg][2])+':'+str(int(round(GQs[ggg])))
                else:
                    out=out+'\t'+GTs[ggg]+':'+str(10**GLs[ggg][0])+','+str(10**GLs[ggg][1])+','+str(10**GLs[ggg][2])+':'+str(int(round(GQs[ggg])))

        elif 0.0<AF<1.0:###Only include SNPs that are variable in the final VCF
            for ggg in range(len(GTs)):
                if GTs[ggg]=='./.':
                    out=out+'\t'+GTs[ggg]+':'+str(GLs[ggg][0])+','+str(GLs[ggg][1])+','+str(GLs[ggg][2])+':'+str(int(round(GQs[ggg])))
                else:
                    out=out+'\t'+GTs[ggg]+':'+str(10**GLs[ggg][0])+','+str(10**GLs[ggg][1])+','+str(10**GLs[ggg][2])+':'+str(int(round(GQs[ggg])))

        out=out+'\n'
        # return out
        # return outfile.write(out)
        return out

    # if g%100==0:
    #     print('Processed '+BAMin+' up to ' +chromo+':'+str(pos)+', SNP nb = '+str(g))    

# entry point for the program
if __name__ == '__main__':

    start_time = time.perf_counter()
    # create the process pool
    with Pool(processes=15) as pool:
        # # call the same function with different data in parallel
        for result in pool.imap(task, range(len(SNPs))):
        #     # report the value to show progress
            outfile.write(result)

    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")

    outfile.close()
 
print('Finished. Laterz')



#############################################################################################################

# #!/usr/bin/env python
# # -*- coding: ASCII -*-
 
# #####This python script takes a list of bamfiles and performs variant calling such that the output is in the format required for lcMLkin to calculate pairwise relatedness.
# #####Please note it is slow, and is only designed for calling genotypes for specific SNPs in smallish numbers (tens to hundreds) of individuals. I do not recommend using this as a general variant caller.
# #####We are working on a better version (this was an experimental one developed to test lcMLkin), but for now, if using 100 samples with 100K loci with 20x mean coverage takes about 4 hours on our machine. Lower coverage datasets, which is what lcMLkin is designed for, will be quicker though, maybe just over an hour.
# #####Pysam and Numpy must be installed in the version of python used.
# #####Bam files must be indexed.
# #####The SNP file must have the tab seperated fields in the following order: chromosome, position (one-based), reference allele, alternate allele
# #####If you want to change things like mapping and base quality threshold, edit the python code under the section "Input arguments"
# #####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)

# ##### New:multithreading/multiproccessing in the script. As such, some functions were adapted.  
# #####  Changed so it is recognized by bcftools and can operate with it
 
# #####usage is ./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file> 
 
 
# ###import libraries
# import string
# import numpy as np
# import pysam
# import gzip
# import math
# import copy
# from sys import argv
# import time
# from multiprocessing import Pool 
# import os
# from concurrent.futures import ProcessPoolExecutor, as_completed
# import time
# import multiprocessing

 
# ###Input arguments
# BAMin=argv[1] #flat text file consisiting of list of bams (must be indexed), one per line
# filenameout=argv[2] #creates a vcf
# SNPfile=argv[3] #must have the tab seperated fields in the following order chromosome, position (one-based), reference allele, alternate allele
# BC=1 #choose whether to use the bayesian caller for individual genotypes (0=no, 1=yes). I don't recommend turning this on the current  version, as it's quite slow, and allele frequencies shouldn't really be biased if the genotypes are incorreclty called, especially as the genotype calling we use is not reference aware
# MQ_t=20 #mapping quality threshold
# BQ_t=5 #base_qualitythreshold
# GQ_t=0.1 #GQ threshold
# min_sams=1 #minimum samples with reads at locus needed to attempt bayesian calling of genotypes
# exc_fix=0  #exclude fixed differences (doesn't make much sense to turn on if only looking at one sample
 
# ###converts phred score to probability
# def phred2prob(x):
#     return 10.0**(-x/10.0)
 
# ###converts probability to phred score
# def prob2phred(x):
#     return -10*math.log10(x)
 
# ###diploid caller assuming all alternative alleles are possible (error is divided by three for now, could add a more complex model)
# def geno_caller_10GT(X):
 
#     GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
#     hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T
 
#     all_dic={}
#     all_dic['A']=0
#     all_dic['C']=1
#     all_dic['G']=2
#     all_dic['T']=3
 
#     count=0
#     for g in range(len(X)):
#         if X[g][0] not in all_dic:
#             continue
#         hap[g]=phred2prob(X[g][1])*(1.0/3.0)
#         hap[g][all_dic[X[g][0]]]=1-phred2prob(X[g][1])
 
#         GL[0]=GL[0]+math.log10(hap[g][0])
#         GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)
#         GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
#         GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)
 
#         GL[4]=GL[4]+math.log10(hap[g][1])
#         GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
#         GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)
 
#         GL[7]=GL[7]+math.log10(hap[g][2])
#         GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)
 
#         GL[9]=GL[9]+math.log10(hap[g][3])
#         count+=1
 
#     if count==0:
#         GL.fill(-9.0)
#     return GL
 
# all_dic={}
# all_dic['A']=0
# all_dic['C']=1
# all_dic['G']=2
# all_dic['T']=3
 
# def geno_caller_3GT(X,ref,alt,all_dic):
#     #diploid caller assuming that only assesses likelihood for three possible genotypes (ref/ref,ref/alt,alt/alt)
#     GL=[0.0,0.0,0.0]
 
#     count=0
#     for g in range(len(X)):
#         if X[g][0] not in all_dic:
#             continue
#         err=phred2prob(X[g][1])*(1.0/3.0)
#         tru=1-phred2prob(X[g][1])
 
#         if X[g][0]==ref:
#             GL[0]=GL[0]+math.log10(tru)
#             GL[1]=GL[1]+math.log10((tru+err)/2)
#             GL[2]=GL[2]+math.log10(err)
#         elif X[g][0]==alt:
#             GL[0]=GL[0]+math.log10(err)
#             GL[1]=GL[1]+math.log10((tru+err)/2)
#             GL[2]=GL[2]+math.log10(tru)
#         else:
#             GL[0]=GL[0]+math.log10(err)
#             GL[1]=GL[1]+math.log10(err)
#             GL[2]=GL[2]+math.log10(err)
#         count+=1
 
#     if count==0:
#         GL=[-9.0,-9.0,-9.0]
#     return GL
 
 
# ###extract reads for a give position in a bam
# def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
#     var_list=[]
#     for pileupcolumn in samfile.pileup(chromo,pos-1,pos,truncate=True,stepper='all'):
#         for pileupread in pileupcolumn.pileups:
#             if not pileupread.is_del and not pileupread.is_refskip:
#                 if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):
#                         var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
#     return var_list
 
# #####extract reads for a give position in a bam (original version of above,not clear which is faster)
# ##def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
# ##    var_list=[]
# ##    for pileupcolumn in samfile.pileup(chromo,pos-1,pos):
# ##        if pileupcolumn.pos==pos-1:
# ##            for pileupread in pileupcolumn.pileups:
# ##                if not pileupread.is_del and not pileupread.is_refskip:
# ##                    if pileupread.alignment.mapping_quality>=MQ:
# ##                        if pileupread.alignment.is_duplicate==False:
# ##                            if ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ:
# ##                                var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
# ##    return var_list
 
 
# ###Version of Depristo's bayesian caller
# def bayes_caller(X):
#     het=0.001
#     nb_ind=len(X)
 
 
#     thetas=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix with theta priors
#     for g in range(1,len(thetas)):
#         thetas[g]=het/float(g)
#     thetas[0]=1-sum(thetas[1:])
 
#     best_genos=np.zeros(((nb_ind*2)+1,nb_ind),dtype='int64')  #matrix for storing best genotypes for a given allele frquency
#     best_genos_ll=np.zeros(((nb_ind*2)+1,nb_ind),dtype='float32')  #matrix for storing best genotypes for a given allele frquency
 
#     p_best_genos=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix for storing the likelihood for each p
 
#     for g in range(len(X)):   ### p = 0
#         best_genos_ll[0][g]=best_genos_ll[0][g]+X[g][0]
 
#     p_best_genos[0]=np.sum(best_genos_ll[0])
 
#     ##perform greedy algorithm to find best likelihoods for each possible configuration as described in Depristo
#     for g in range(1,len(p_best_genos)):
#         p_best_genos[g]=-1000000000000
#         genos=np.copy(best_genos[g-1])
#         genos_ll=np.copy(best_genos_ll[g-1])
#         for gg in range(len(genos)):
#             if genos[gg]<2:
#                 genos_mod=np.copy(genos)
#                 genos_mod[gg]+=1
#                 genos_ll_mod=np.copy(genos_ll)
#                 genos_ll_mod[gg]=X[gg][genos[gg]+1]
#                 ll=np.sum(genos_ll_mod)
#                 if ll>p_best_genos[g]:
#                     p_best_genos[g]=ll
#                     best_genos[g]=np.copy(genos_mod)
#                     best_genos_ll[g]=np.copy(genos_ll_mod)
 
 
#     ##Calculate variant quality
#     #denom=np.max(p_best_genos)+np.log10(np.sum(10**(p_best_genos-np.max(p_best_genos))))
#     #phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)
#     marg_comp=np.log10(thetas)+p_best_genos
#     denom=np.max(marg_comp)+np.log10(np.sum(10**(marg_comp-np.max(marg_comp))))
#     phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)
 
#     VQ=int(round(phred))
 
#     geno_comb=best_genos[np.argmax(p_best_genos)]
 
#     if np.sum(geno_comb)==0:
#         stat='REF'
 
#     else:
#         stat='ALT'
#     return geno_comb,VQ,stat
 
# ###dictionary for GL order
# geno_ord={}
# geno_ord['AA']=0
# geno_ord['AC']=1
# geno_ord['AG']=2
# geno_ord['AT']=3
# geno_ord['CA']=1
# geno_ord['CC']=4
# geno_ord['CG']=5
# geno_ord['CT']=6
# geno_ord['GA']=2
# geno_ord['GC']=5
# geno_ord['GG']=7
# geno_ord['GT']=8
# geno_ord['TA']=3
# geno_ord['TC']=6
# geno_ord['TG']=8
# geno_ord['TT']=9
 
# ###dictionary for vcf genotype codes
# vcf_gt={}
# vcf_gt[0]='0/0'
# vcf_gt[1]='0/1'
# vcf_gt[2]='1/1'
 
 
# ###Read bamlist
# file=open(BAMin,'r')
# samples=file.read()
# samples=str.split(samples,'\n')
# if samples[-1]=='':
#     del(samples[-1])

# # Open file with samples names
# ###Get sample name for each bam. Each bam must only contain one sample
# child_names = []
# with open(BAMin, "r") as f:
#     for line in f:
#         file_path = line.strip()
#         child_name = os.path.basename(file_path).split(".")[0]
#         child_names.append(child_name)
        
 
# ###Read SNPlist
# file=open(SNPfile,'r')
# SNPs=file.read()
# SNPs=str.split(SNPs,'\n')
# if SNPs[-1]=='':
#     del(SNPs[-1])
 
# # sample_name=[]
# # samfiles={}
# # ###Get sample name for each bam. Each bam must only contain one sample, otherwise the program will take the sample in first readgroup in the header
# # for g in range(len(samples)):
# #     samfiles[g] = pysam.AlignmentFile(samples[g], "rb")
# #     sample_name.append(samfiles[g])#.header['RG'][0]['SM'])
 
 
# ###start writing outfile with header
# outfile=open(filenameout,'w')
# out="""##fileformat=VCFv4.2
# ##bcftoolsVersion=1.9+htslib-1.9
# ##FILTER=<ID=.,Description="Missing FILTER">
# ##contig=<ID=1,length=100100000>
# ##ALT=<ID=*,Description="Represents allele(s) other than observed.">
# ##INFO=<ID=AF,Number=A,Type=Float,Description="allele frequencies (from multisampling or single-base calling)">
# ##FORMAT=<ID=GL,Number=G,Type=Float,Description="genotype likelihoods">
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""
# # out='##Custom variant caller for_lcMLkin\n##Bayesian_caller='+str(BC)+',MappingQuality_filter='+str(MQ_t)+',BaseQuality_filter='+str(BQ_t)+',GenotypeQuality_filter='+str(GQ_t)+',bamlist='+BAMin+',filenameout='+filenameout+',SNPfile='+SNPfile+'\n'
# out=out+'##Time created: '+(time.strftime("%H:%M:%S"))+' '+(time.strftime("%d/%m/%Y"))+'\n'
# out=out+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
# for x in child_names:
#     out=out+'\t'+ str(x)
 
# out=out+'\n'
 
# outfile.write(out)

# def task(g):
 
#     ###Iterate through SNPs
#     # for g in range(len(SNPs)):
#     k=str.split(SNPs[g])
#     # return k
#     chromo=k[0]
#     pos=int(k[1])
#     ref=k[2]
#     alt=k[3]

#     ###iterate through samples
#     GLs=[]
#     GLs_nomiss=[]
#     GQs=[]
#     GTs=[]

#     for gg in range(len(samples)):

#         samfile = pysam.AlignmentFile(samples[gg], "rb")
#         ###extract reads for a position
#         read_list=extract_bam_SNP(samfile,chromo,pos,BQ_t,MQ_t)


#         ###get the genotype likelihood for the 3 possible genotypes based on the ref and alt alleles provided
#         GL=geno_caller_3GT(read_list,ref,alt,all_dic)

#         ###calcule genotype quality
#         GL_sort=copy.deepcopy(GL)
#         GL_sort.sort()
#         GQ = -10 * (GL_sort[1] - GL_sort[2])
#         GQs.append(GQ)

#         ###assign genotype using highest likelihood
#         if BC==1:
#             if GQ>GQ_t:
#                 GLs_nomiss.append(GL)

#         elif BC==0:
#             if GQ>GQ_t:
#                 GTs.append(vcf_gt[np.argmax(np.asarray(GL))])
#             else:
#                 GTs.append('./.')
#                 GLs
#         GLs.append(GL)

#     skip=0

#     ###estimate allele frequencies for multisample calling
#     if BC==1:
#         if len(GLs_nomiss)>min_sams:  #must have genotype data for at least 20 individuals to consider SNP
#             genos_b,VQ,s=bayes_caller(GLs_nomiss)
#             AF=(np.sum(genos_b)/(float(len(genos_b))*2))
#             it=0
#             for ggg in range(len(GLs)):
#                 if -9.0 in GLs[ggg]:
#                     GTs.append('./.')
#                 else:
#                     GTs.append(vcf_gt[genos_b[it]])
#                     it+=1

#         else:
#             skip=1

#     ###estimate allele frequencies for single-base calling
#     elif BC==0:
#         if (float((len(GTs)-GTs.count('./.')))*2) != 0.0:
#             AF=1-((GTs.count('0/1')+GTs.count('1/1')*2) / (float((len(GTs)-GTs.count('./.')))*2))
#             VQ=100
#             s='.'
#         else:
#             skip=1

#     ###only write to file if site is variable (otherwise useless for kinship estimation)
#     if skip==0:
#         out=k[0]+'\t'+k[1]+'\t.\t'+k[2]+'\t'+k[3]+'\t'+str(VQ)+'\t'+s+'\tAF='+str(AF)+'\tGT:GL:GQ'
#         if exc_fix==0: ###Don't throw out sites that are fixed. If examining only a few or one individual, and population allele frequencies in lcMLkin come from elsewhere, this setting makes sense to use.
#             for ggg in range(len(GTs)):
#                 if GTs[ggg]=='./.':
#                     out=out+'\t'+GTs[ggg]+':'+str(GLs[ggg][0])+','+str(GLs[ggg][1])+','+str(GLs[ggg][2])+':'+str(int(round(GQs[ggg])))
#                 else:
#                     out=out+'\t'+GTs[ggg]+':'+str(10**GLs[ggg][0])+','+str(10**GLs[ggg][1])+','+str(10**GLs[ggg][2])+':'+str(int(round(GQs[ggg])))

#         elif 0.0<AF<1.0:###Only include SNPs that are variable in the final VCF
#             for ggg in range(len(GTs)):
#                 if GTs[ggg]=='./.':
#                     out=out+'\t'+GTs[ggg]+':'+str(GLs[ggg][0])+','+str(GLs[ggg][1])+','+str(GLs[ggg][2])+':'+str(int(round(GQs[ggg])))
#                 else:
#                     out=out+'\t'+GTs[ggg]+':'+str(10**GLs[ggg][0])+','+str(10**GLs[ggg][1])+','+str(10**GLs[ggg][2])+':'+str(int(round(GQs[ggg])))

#         out=out+'\n'
#         # return out
#         # return outfile.write(out)
#         return out

#     # if g%100==0:
#     #     print('Processed '+BAMin+' up to ' +chromo+':'+str(pos)+', SNP nb = '+str(g))    

# # entry point for the program
# if __name__ == '__main__':

#     start_time = time.perf_counter()
#     # create the process pool
#     with Pool() as pool:
#         # # call the same function with different data in parallel
#         for result in pool.imap(task, range(len(SNPs))):
#         #     # report the value to show progress
#             outfile.write(result)

#     finish_time = time.perf_counter()
#     print(f"Program finished in {finish_time-start_time} seconds")

#     outfile.close()
 
# print('Finished. Laterz')
