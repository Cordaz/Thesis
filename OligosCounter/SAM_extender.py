#!/usr/bin/python

import sys

#Constant
dataset_path = "../../dataset/"
genome_path = dataset_path + "genome/"
result_path = dataset_path + "extended/"

pos_strand = 0
neg_strand = 16
not_mapped = 4

ext_length = 100
read_length = 34

#Global vars
chrom_in_use = ""
chrom_seq = ""
start = 0
chrom_length = {}

#Functions

def to_fasta(path):
    filename = ""
    for c in path:
        filename = filename + c
        if c == ".":
            filename = filename + "fa"
            break
    return filename


def get_value(string):
    return string.split(":")[1]

def get_pair_chr_length(line):
    args = line.strip().split("\t")
    return (get_value(args[1]), get_value(args[2]))

def get_sam_info(line):
    args = line.strip().split("\t")
    return (args[0], int(args[1]), args[2], int(args[3])) #(>id, strand, chrom, pos)

def get_chromosome_seq(chrom):
    seq = ""
    with open(genome_path + chrom + ".fa", "r") as chrom_file:
        next(chrom_file)
        for line in chrom_file:
            seq = seq + line.strip().upper()
    return seq


def check_seq(s):
    return 'N' not in s


def complementary(seq):
   comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
   c_seq = ""
   for char in reversed(seq):
      c_seq = c_seq + comp_bases[char]

   return c_seq

def extend_read(chrom_index, seq, strand, chrom, pos):
    if strand == not_mapped:
        return "", chrom_index

    while chrom_index < chrom_length[chrom]:
        if chrom_index == pos:
            break
        chrom_index = chrom_index + 1

    if strand == pos_strand:
        return extend_pos(seq, chrom_index), chrom_index
    if strand == neg_strand:
        return extend_neg(seq, chrom_index), chrom_index
    return "", chrom_index

def extend_pos(chrom_seq, chrom_index):
    seq = chrom_seq[chrom_index -1 : chrom_index + ext_length - 1]
    if check_seq(seq):
        return seq
    return ""

def extend_neg(chrom_seq, chrom_index):
    seq = chrom_seq[chrom_index - (ext_length - read_length) - 1 : chrom_index + read_length - 1 - 1]
    if check_seq(seq):
        return complementary(seq)
    return ""

def file_len(fname):
    with open(fname, "r") as f:
        for i, l in enumerate(f):
            pass
    return i+1


#Main


args = sys.argv[1:]
if len(args) < 1:
    print "usage: input_file"
    sys.exit(1)

input_file_path = dataset_path + args[0]

input_length = file_len(input_file_path)

with open(input_file_path, "r") as input_file:
    output_file = open(result_path + to_fasta(input_file_path), "w+")
    line = input_file.readline()
    while line[0] == '@':
        (chrom, length) = get_pair_chr_length(line)
        chrom_length[chrom] = length
        line = input_file.readline()
    #Process first line before loop
    (id_seq, strand, chrom, pos) = get_sam_info(line)
    chrom_in_use = chrom
    chrom_seq = get_chromosome_seq(chrom)
    #Seaarch for read from starting point
    seq, start = extend_read(start, chrom_seq, strand, chrom, pos)
    progress = 1
    if seq:
        output_file.write(id_seq + "\n" + seq)

    #Process rest of file
    for line in input_file:
        (id_seq, strand, chrom, pos) = get_sam_info(line)
        if chrom_in_use != chrom:
            chrom_seq = get_chromosome_seq(chrom)
            chrom_in_use = chrom
            start = 0

        seq, start = extend_read(start, chrom_seq, strand, chrom, pos)
        if seq:
            output_file.write("\n" + id_seq + "\n" + seq)

        progress = progress + 1
        if (float(progress) / input_length) * 100 % 5 == 0:
            print "*"

    output_file.close()
