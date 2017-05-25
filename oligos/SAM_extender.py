#!/usr/bin/python

import sys
import os

home = os.getenv("HOME")

'''

Works with nohup output: please be careful to remove already existent 'nohup.out'.

input file must be in ~/dataset/
genome files ('chr_.fa') must be in ~/dataset/genome/

to change:
'''
dataset_path = home + "/dataset/"
genome_path = dataset_path + "genome/"

#Extension
ext_length = 34


#Constant
pos_strand = 0
neg_strand = 16
not_mapped = 4

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
    return (args[0], int(args[1]), args[2], int(args[3]), len(args[9])) #(>id, strand, chrom, pos, read_len)

def get_chromosome_seq(chrom):
    seq = ""
    with open(genome_path + chrom + ".fa", "r") as chrom_file:
        next(chrom_file)
        for line in chrom_file:
            seq = seq + line.strip().upper()
    return seq


def extend_read(chrom_index, seq, strand, chrom, pos, read_len):
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
        if (chrom_len - chrom_index) < read_len:
            return "", chrom_index
        return extend_neg(seq, chrom_index, read_len), chrom_index
    return "", chrom_index

def extend_pos(chrom_seq, chrom_index):
    return chrom_seq[chrom_index - 1 : chrom_index + ext_length - 1]

def extend_neg(chrom_seq, chrom_index, read_len):
    seq = chrom_seq[chrom_index - (ext_length - read_len) - 1 : chrom_index + read_len - 1]
    return seq


#Main

args = sys.argv[1:]
if len(args) < 2:
    print "usage: input_file extension"
    sys.exit(1)

input_file_path = args[0]
ext_length = int(args[1])
extend = True
if ext_length == 0:
    extend = False

last_pos = -1

with open(input_file_path, "r") as input_file:
    line = input_file.readline()
    while line[0] == '@':
        if line[0:3] == '@SQ':
            (chrom, length) = get_pair_chr_length(line)
            chrom_length[chrom] = length
        line = input_file.readline()
    #Process first line before loop
    (id_seq, strand, chrom, pos, read_len) = get_sam_info(line)
    chrom_in_use = chrom
    chrom_seq = get_chromosome_seq(chrom)
    #Seaarch for read from starting point
    if not extend:
        ext_length = read_len
    seq, start = extend_read(start, chrom_seq, strand, chrom, pos, read_len)
    if seq:
        print ">" + id_seq + "\n" + seq
	last_pos = pos

    #Process rest of file
    for line in input_file:
        (id_seq, strand, chrom, pos, read_len) = get_sam_info(line)
        if chrom != '*':
            if chrom_in_use != chrom:
                chrom_seq = get_chromosome_seq(chrom)
                chrom_in_use = chrom
                start = 0
                last_pos = -1
            if pos > last_pos:
                if not extend:
                    ext_length = read_len
                seq, start = extend_read(start, chrom_seq, strand, chrom, pos, read_len)
                last_pos = pos
                if seq:
                    print ">" + id_seq + "\n" + seq

'''
output_file = open(result_path + to_fasta(input_file_path), "w+")
for id_seq in reads:
    output_file.write(id_seq + "\n" + reads[id_seq] + "\n")
output_file.close()
'''
