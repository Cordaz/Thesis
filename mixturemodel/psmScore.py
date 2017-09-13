'''
@author Andrea Corneo
@year 2017

Implementation of a JASPAR/region best score calculator, given a FASTA file.

'''

#!/usr/bin/python

import sys
import numpy as np
import re
import math
# sys.path.append("../utilities") # Not necessary anymore. Use for import custom libs

# return the index associated to each base/pairs
def get_base_index(b):
    if b == 'A':
        return 0
    if b == 'C':
        return 1
    if b == 'G':
        return 2
    if b == 'T':
        return 3
    return -1

# complements a base, that is A <-> T, C <-> G
def complement_base(b):
    if b == 'A':
        return 'T'
    if b == 'T':
        return 'A'
    if b == 'C':
        return 'G'
    if b == 'G':
        return 'C'
    return 'N'

# complements a sequence
# return complementary as list
def complement(seq):
    c_seq = []
    for b in reversed(list(seq)):
        c_seq.append(complement_base(b))
    return c_seq

argv = sys.argv[1:]

if not argv:
    print "usage: python psmScore.py psm_file fasta_file score_output_file"
    sys.exit(1)
psm_file = argv[0]
del argv[0]
if not argv:
    print "usage: python psmScore.py psm_file fasta_file score_output_file"
    sys.exit(1)
fa_file = argv[0]
del argv[0]
if not argv:
    print "usage: python psmScore.py psm_file fasta_file score_output_file"
    sys.exit(1)
score_file = argv[0]
del argv[0]
if argv:
    print "usage: python psmScore.py psm_file fasta_file score_output_file"
    sys.exit(1)


with open(psm_file, 'r') as psm_fp:
    line = psm_fp.readline()
    # determine length to alloc matrix
    values = line.split(' ')[:-1] # remove last empty

    d = len(values)

    C = np.zeros( (4, d), dtype=np.double ) # counts matrix

    for i in xrange(d):
        C[0,i] = float(values[i])

    # compute other lines
    for j in xrange(1,4):
        line = psm_fp.readline()
        values = line.split(' ')[:-1]
        for i in xrange(d):
            C[j,i] = float(values[i])

# print C

# Normalize by column

# np version
col_sums = C.sum(axis=0) # axis 0 is column
F = C / col_sums # frequencies matrix

'''
# extended version (equivalent to np, just for double check)
F = np.zeros( (4,d), dtype=np.double ) # frequencies matrix
for i in xrange(d):
    col_sum = 0.0
    for j in xrange(4):
        col_sum += C[j,i]
    for j in xrange(4):
        F[j,i] = C[j,i] / col_sum
'''

# print F

# Add 0.01 and renormalize

# np version
F = F + 0.01
col_sums = F.sum(axis=0)
PSM = F / col_sums

'''
# extended version (equivalent to np, just for double check)
PSM = np.zeros( (4,d), dtype=np.double )
for i in xrange(d):
    col_sum = 0.0
    for j in xrange(4):
        F[j,i] += 0.01
        col_sum += F[j,i]
    for j in xrange(4):
        PSM[j,i] = F[j,i] / col_sum
'''

# compute normalization factor
p_max = 1.0
p_min = 1.0
for i in xrange(d):
    max_of_col = -2.0 # ensure out of bound
    min_of_col = 2.0
    for j in xrange(4):
        if PSM[j,i] > max_of_col:
            max_of_col = PSM[j,i]
        if PSM[j,i] < min_of_col:
            min_of_col = PSM[j,i]
    p_max = p_max * max_of_col
    p_min = p_min * min_of_col

p_max = math.log(p_max, 2)
p_min = math.log(p_min, 2)
norm_denom = p_max - p_min

# print PSM
# print p_max
# print p_min
# print norm_denom

### PSM generated, can now compute on fa file

already_processed = {}

with open(fa_file, 'r') as fa_fp:
    with open(score_file, 'w+') as sc_fp:
        sc_fp.write("CHROM\tREG_START\tREG_END\tSCORE\tSITE\n")
        l = fa_fp.readline().strip()
        while True:
            if not l:
                break
            if l.startswith(">"):
                header = l
                match = re.search("(chr\w+):(\d+)[-]?(?:\d+)", header)
                # print header
                # print match.group(0)
                # check if is a duplicate
                if (match.group(1), int(match.group(2))) in already_processed:
                    # is a duplicate, just read till new sequence and restart the loop
                    l = fa_fp.readline().strip()
                    if not l:
                        break # EOF
                    while not l.startswith('>'):
                        l = fa_fp.readline().strip()
                        if not l:
                            break # EOF
                    continue
                else:
                    already_processed[(match.group(1), int(match.group(2)))] = True
                # Get coordinates
                geco = ( match.group(1), int(match.group(2)) ) #GEnomic COordinates
                # print geco
            l = fa_fp.readline().strip()
            seq = []

            while not l.startswith(">"):
                for c in list(l):
                    seq.append(c)

                l = fa_fp.readline().strip()
                if not l:
                    break

            # print "".join(seq) + "\n" + str(geco)

            max_score = p_min
            for i in xrange(len(seq) - d + 1):
                score = 1.0
                flag = True
                for j in xrange(d): # sliding window size of matrix
                    bi = get_base_index(seq[i+j])
                    if bi < 0:
                        flag = False
                        break
                    pre_score = score
                    score *= PSM[bi,j]
                    # print "Base: " + seq[i+j] + "\tPre_score= " + str(pre_score) + "\tPSM[" + str(bi) + "," + str(j) + "]=" + str(PSM[bi,j]) + "\tscore= " + str(score)
                if flag and (score > max_score):
                    max_score = score
                    max_offset = i

            # complementary
            c_seq = complement(seq)
            # print "".join(seq)
            # print "".join(c_seq)

            for i in xrange(len(c_seq) - d + 1): # c_seq is already reversed, slide positive
                score = 1.0
                flag = True
                for j in xrange(d): # sliding window size of matrix
                    bi = get_base_index(c_seq[i+j])
                    if bi < 0:
                        flag = False
                        break
                    score *= PSM[bi,d-j-1] # consider matrix as reversed
                if flag and (score > max_score):
                    max_score = score
                    max_offset = i + d - 1

            chrom, pos = geco
            score = ( math.log(max_score,2) - p_min ) / norm_denom # normalize score
            sc_fp.write(chrom + "\t" + str(pos + max_offset) + "\t" + str(pos + max_offset + d) + "\t" + str(score) + "\t" + "".join(seq[max_offset : max_offset + d]) + "\n")
