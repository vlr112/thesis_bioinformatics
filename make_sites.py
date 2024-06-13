#!/usr/bin/python
import numpy as np
import random
import copy
import time
import os
import glob

# We want 100k variable sites with 1000 sites in between each one
sites = 100000
invariable_sites = 1000

# Start with four founders
founders = 30

# Basic DNA
geno = [0, 1]
nucs = ["A", "C", "G", "T"]
ploidy = len(geno)

# Start time
last_time = time.time()

# Hack for sampling random nucleotides
def restore():
    return copy.deepcopy(nucs)

# Nice print statements with time elapsed
def msg(m):
    global last_time
    this_time = time.time()
    elapsed = this_time - last_time
    print("%s - took %.5f secs" % (m, elapsed))
    last_time = this_time

# Make an ancestral genome
# This is just a random choice from the four nucleotides
# One for each variable site
def get_ancestral():
    ancestral = random.choices(nucs, k=sites)
    save_genome(ancestral, "ancestral")
    return ancestral

# Derived genome
# I'm sure there is a nicer way of doing this
def get_derived(ancestral):
    derived = []
    nucs_der = restore()
    # For each site in ancestral genome
    for a in ancestral:
        # Choose a random nucleotide that is not the ancestral
        nucs_der.remove(a)
        derived.append(random.choices(nucs_der)[0])
        nucs_der = restore()
    save_genome(derived, "derived")
    return derived

# This makes random genotypes for the founder individuals based on the allele frequencies
# Returns a dictionary because dictionaries are my fave
def gen_genotypes():
    ind_genotypes = {}
    for ind in range(0, founders):
        # Each founder is a list of {0,1,2} for each site
        ind_genotypes[ind] = []
        for f in range(0, sites):
            freq = allele_freqs[f] # Get the allele freq for this site
            # This part samples a random value [0 or 1] with weight given by allele freq
            # It does this twice and sums them together
            # For example freq=0.7, this is the chance of the allele at this site being ancestral,
            # So we sample it from [0,1] with the weights [0.7,0.3]
            # Let's say once it was 0 and once 1, so the genotype is 0+1 = heterozygous
            ind_genotypes[ind].append(sum(random.choices(geno, weights=[freq,1-freq], k=2)))
    return ind_genotypes

# Two individuals fall in love
# This function takes two lists of genotypes
def mate(ind1, ind2):
    # Child genotype will also be a list
    child_geno = []
    # For each site
    for g in range(0, sites):
        # Empty child genotype
        this_child_geno = 0
        # For each parent
        for ind in [ind1, ind2]:
            # Get genotype
            geno = ind[g]
            # If it is a heterozygous parent
            if geno == 1:
                # Sample 0 or 1 (anc or der) and add to child's genotype
                r = random.randint(0, 1)
                this_child_geno += r
            # If homozgous anc, it gives anc to child
            if geno == 0:
                this_child_geno += 0
            # If homozgous der, it gives der to child
            if geno == 2:
                this_child_geno += 1
        # Since we looped through both parents we now have a list of genotypes 0,1,2 for the child
        child_geno.append(this_child_geno)
    return child_geno

# Here we actually get the nucleotides behind the 0,1,2 lists
def geno_out(ind):
    # Two lists per individual, one per haplotype
    f1, f2 = [], []
    for g in range(0, sites):
        geno = ind[g]
        # Get the nucleotide at this site for the ancestral and the derived
        anc = ancestral[g]
        der = derived[g]
        if geno == 0:
            f1.append(anc)
            f2.append(anc)
        if geno == 1:
            f1.append(anc)
            f2.append(der)
        if geno == 2:
            f1.append(der)
            f2.append(der)
    return f1, f2

# Write the variable sites as a list
def save_sites(f, outname):
    out = "variable_sites/" + outname + ".txt"
    with open (out, 'w') as file:
        for i in f:
            o = i + "\n"
            file.write(o)

# Makes the inveriable reference genome
def make_ref():
    ref = []
    # For each of the variable sites
    for i in range(0,sites):
        # Choose 1000 (or whatever) random nucleotides
        random_sites = random.choices(nucs, k=invariable_sites)
        for s in random_sites:
            # and add to reference
            ref.append(s)
        # then add our one base for the variable site
        ref.append("N")
    save_genome(ref, "reference_genome")

# Prints variable sites
def print_sites(f, id):
    f1, f2 = geno_out(f)
    save_sites(f1, (id + ".h1"))
    save_sites(f2, (id + ".h2"))

# Saves fasta genomes
def save_genome(fa, outname):
    out = "fasta/" + outname + ".fa"
    with open (out, 'w') as file:
        header = ">" + outname + "\n"
        file.write(header)
        file.write("".join(fa))
        file.write("\n")

# Saves the allele freq vector
def save_freq(allele_freqs):
    out = "af.txt"
    with open (out, 'w') as file:
        for i in allele_freqs:
            o = str(i) + "\n"
            file.write(o)

# CAREFUL: This removes all files in variable sites dir, just to avoid having old ones
def remove_old():
    files = glob.glob('variable_sites/*')
    for f in files:
        os.remove(f)


# Generate Nsites allele frequencies
allele_freqs = np.random.uniform(low=0.1, high=0.9, size=(sites,))
msg("Generated allele frequencies")

# Saves them as a file
save_freq(allele_freqs)
msg("Saved frequencies")

# Generate genotypes for founders
founders = gen_genotypes()
msg("Created founders")

# Gets the ancestral and derived variable sites
ancestral = get_ancestral()
msg("Generated ancestral")
derived = get_derived(ancestral)
msg("Generated derived")

children = {}
children["A"] = mate(founders[0], founders[1])
children["B"] = mate(founders[0], founders[1])

# half family
children["C"] = mate(founders[1], founders[2])
children["G"] = mate(children["C"], founders[5])
###

children["D"] = mate(founders[3], children["A"])
children["E"] = mate(founders[4], children["B"])
children["F"] = mate(founders[4], children["B"])

## INBREED ##
children["H"] = mate(children["A"], children["E"])
children["I"] = mate(children["F"], children["E"])


# Makes the reference genome, with Ns in the variable positions
make_ref()
msg("Made invariable reference")

# Removes old stuff, delete if scared
remove_old()
msg("Removed old stuff in site dir")

# Saves the allele at the variable sites for all of the population
for f in founders:
    msg("Saving founder " + str(f))
    id = "founder%i" % f
    print_sites(founders[f], id)
for c in children:
    msg("Saving child " + str(c))
    id = "child%s" % c
    print_sites(children[c], id)
