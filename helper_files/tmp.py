# python script to test temporary ideas

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

## check what is wrong with the new_reference_genome.fa. 

fasta = Fasta(sys.argv[1]) # VCF file path

print(fasta.keys())

print(len(fasta['1']))