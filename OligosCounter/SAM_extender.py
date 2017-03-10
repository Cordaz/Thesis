#!/usr/bin/python

'''

Works with nohup output: please be careful to remove already existent 'nohup.out'.

'''

import sys
import os

home = os.getenv("HOME")

#Constant
dataset_path = home + "/dataset/"
genome_path = dataset_path + "genome/"

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
            filename = filename + "." + str(ext_length) + ".txt"
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


def complementary(seq):
   comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
   c_seq = ""
   for char in reversed(seq):
      c_seq = c_seq + comp_bases[char]

   return c_seq

def extend_read(chrom_index, seq, strand, chrom, pos):
    chrom_len = int(chrom_length[chrom])
    if strand == not_mapped:
        return "", chrom_index

    while chrom_index < chrom_len:
        if chrom_index == pos:
            break
        chrom_index = chrom_index + 1

    if strand == pos_strand:
        if (chrom_len - chrom_index) <= ext_length:
            return "", chrom_index
        return extend_pos(seq, chrom_index), chrom_index
    if strand == neg_strand:
        if (chrom_len - chrom_index) < read_length:
            return "", chrom_index
        return extend_neg(seq, chrom_index), chrom_index
    return "", chrom_index

def extend_pos(chrom_seq, chrom_index):
    return chrom_seq[chrom_index - 1 : chrom_index + ext_length - 1]

def extend_neg(chrom_seq, chrom_index):
    seq = chrom_seq[chrom_index - (ext_length - read_length) - 1 : chrom_index + read_length - 1]
    return complementary(seq)


#Main

args = sys.argv[1:]
if len(args) < 1:
    print "usage: input_file"
    sys.exit(1)

input_file_path = dataset_path + args[0]

with open(input_file_path, "r") as input_file:
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
    if seq:
        print ">" + id_seq + "\n" + seq

    #Process rest of file
    for line in input_file:
        (id_seq, strand, chrom, pos) = get_sam_info(line)
        if chrom_in_use != chrom:
            chrom_seq = get_chromosome_seq(chrom)
            chrom_in_use = chrom
            start = 0

        seq, start = extend_read(start, chrom_seq, strand, chrom, pos)
        if seq:
            print ">" + id_seq + "\n" + seq

'''
output_file = open(result_path + to_fasta(input_file_path), "w+")
for id_seq in reads:
    output_file.write(id_seq + "\n" + reads[id_seq] + "\n")
output_file.close()
'''
