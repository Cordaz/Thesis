#!/usr/bin/python

import sys
import numpy as np
import re

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
col_sums = C.sum(axis=0) # axis 0 is column
F = C / col_sums # frequencies matrix

# print F

# Add 0.01 and renormalize
F = F + 0.01
col_sums = F.sum(axis=0)
PSM = F / col_sums

# compute normalization factor
p_max = 1.0
p_min = 1.0
for i in xrange(d):
    max_of_col = 0
    min_of_col = 1
    for j in xrange(4):
        if PSM[j,i] > max_of_col:
            max_of_col = PSM[j,i]
        if PSM[j,i] < min_of_col:
            min_of_col = PSM[j,i]
    p_max = p_max * max_of_col
    p_min = p_min * min_of_col

norm_denom = p_max - p_min

# print PSM

## Create specular matrix for complementary kmer
c_PSM = np.zeros( (4, d), dtype=np.double )
for i in xrange(d):
    c_PSM[0,i] = PSM[3,d-i-1] # A row becomes T row right-left
    c_PSM[3,i] = PSM[0,d-i-1] # viceversa
    c_PSM[1,i] = PSM[2,d-i-1] # same as above for C/G
    c_PSM[2,i] = PSM[1,d-i-1]

# print c_PSM

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

            max_score = 0.0
            for i in xrange(len(seq) - d + 1):
                score = 1.0
                for j in xrange(d):
                    bi = get_base_index(seq[i+j])
                    if bi < 0:
                        score = 0.0
                        break
                    score *= PSM[bi,j]
                if score > max_score:
                    max_score = score
                    max_offset = i

            c_max_score = 0.0
            for i in xrange(len(seq) - d + 1):
                score = 1.0
                for j in xrange(d):
                    bi = get_base_index(seq[i+j])
                    if bi < 0:
                        score = 0.0
                        break
                    score *= c_PSM[bi,j]
                if score > c_max_score:
                    c_max_score = score
                    c_max_offset = i

            if c_max_score > max_score: # doing this lost info about strand
                max_score = c_max_score
                max_offset = c_max_offset

            if max_score > 0.0:
                chrom, pos = geco
                # print pos
                # print max_offset
                # print pos + max_offset
                score = ( max_score - p_min ) / norm_denom # normalize score
                sc_fp.write(chrom + "\t" + str(pos + max_offset) + "\t" + str(pos + max_offset + d) + "\t" + str(score) + "\t" + "".join(seq[max_offset : max_offset + d + 1]) + "\n")
