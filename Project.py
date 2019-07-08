import numpy as np
from random import *
from scipy.optimize import fsolve
from scipy.optimize import least_squares
from math import log
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from motif_helper import get_motif_column


BASES = ('A', 'C', 'T', 'G')

def gen_sequence(SL)->int:
    sequences = ''
    for i in range(SL):
        sequences += np.random.choice(BASES, p=(0.25, 0.25, 0.25, 0.25))
    return sequences

def f(c, *data):
    a, b, ICPC = data
    tmp = a*log(4*a,2) + b*log(4*b,2)
    return tmp + c*log(4*c,2)+(1-a-b-c)*log(4*(1-a-b-c),2) - ICPC

# Assume, a:P(A), b:P(G), c:P(C), d:P(T)
# We first fix a and b, then use the following two equations to solve for c and d
# a + b + c + d = 1
# a*log(4a) + b*log(4b) + c*log(4c) + d*log(4d) = ICPC
# --->
# a*log(4a) + b*log(4b) + c*log(4c) + (1-a-b-c)*log(4(1-a-b-c)) = ICPC
def gen_motif(ICPC, ML):

    motif = []
    for i in range(ML):
        motif.append(get_motif_column(ICPC))
    print(np.array(motif))
    return np.array(motif)


def gen_subsequence(motif):
    subseq=''
    for prob in motif:
        subseq += np.random.choice(BASES, p=prob)
    return subseq

def substitue(sequences, subseq):
    n = randint(0, len(sequences) - len(subseq))
    new_seq = sequences[0:n] + subseq + sequences[n+len(subseq):]
    return n,new_seq

def build_benchmark(ICPC, ML, SL, SC):
    for set_id in range(10):
        motif = gen_motif(ICPC, ML)
        newseqs, sites = [], []
        for seq_id in range(SC):
            seq = gen_sequence(SL)
            subseq = gen_subsequence(motif)
            n, new_seq = substitue(seq, subseq)
            newseqs.append(new_seq)
            sites.append(n)
        print(newseqs)

        path = "./data/%s-%s-%s-%s/dataset_%s/" %(ICPC,ML,SL,SC,set_id+1)
        # 1. motif.txt
        with open(path+"motif.txt", 'w') as fout:
            fout.write(">MOTIF"+str(set_id)+"\t"+str(ML)+"\n")
            for w in motif:
                fout.write("%s %s %s %s\n" %(tuple(w)))
            fout.write("<")

        # # 2. sequences.fa
        with open(path+"sequences.fa", 'w') as fout:
            for idx, sequence in enumerate(newseqs):
                record = SeqRecord(Seq(sequence, IUPAC.protein), id=str(idx), description="")
                SeqIO.write(record, fout, "fasta")
        # 3. sites.txt
        with open(path + "sites.txt", 'w') as fout:
            for site in sites:
                fout.write(str(site)+"\n")
        # 4. motiflength.txt
        with open(path + "motiflength.txt", 'w') as fout:
            fout.write(str(ML))

def step_one():
    build_benchmark(1, 8, 500, 10)
    build_benchmark(1.5, 8, 500, 10)
    # build_benchmark(2, 8, 500, 10)
    # build_benchmark(2, 6, 500, 10)
    # build_benchmark(2, 7, 500, 10)
    # build_benchmark(2, 8, 500, 5)
    # build_benchmark(2, 8, 500, 20)

step_one()
