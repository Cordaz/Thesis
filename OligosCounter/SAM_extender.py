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
chrom_index = 0
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
    with open(genome_path + chrom + ".fa", "r") as chrom_file:
        next(chrom_file)
        for line in chrom_file:
            global chrom_seq
            chrom_seq = chrom_seq + line.strip()
    chrom_index = 0

def check_chrom(chrom):
    if chrom_in_use != chrom:
        chrom_in_use = chrom
        get_chromosome_seq(chrom)


def extend_read(strand, chrom, pos):
    if strand == not_mapped:
        return ""

    while chrom_index < chrom_length[chrom]:
        if chrom_index == pos:
            break
        chrom_index = chrom_index + 1

    if strand == pos_strand:
        return extend_pos()
    if strand == neg_strand:
        return extend_neg()
    return ""

def extend_pos():
    return chrom_seq[chrom_index : chrom_index + ext_length]

def extend_neg():
    temp = chrom_seq[chrom_index - (ext_length - read_length) : chrom_index + read_length - 1 ]
    return temp[::-1]




#Main


args = sys.argv[1:]
if len(args) < 1:
    print "usage: input_file"
    sys.exit(1)

input_file_path = args[0]


with open(dataset_path + input_file_path, "r") as input_file:
    output_file = open(result_path + to_fasta(input_file_path), "w+")
    line = input_file.readline()
    while line[0] == '@':
        (chrom, length) = get_pair_chr_length(line)
        chrom_length[chrom] = length
        line = input_file.readline()
    #Process first line before loop
    (id_seq, strand, chrom, pos) = get_sam_info(line)
    chrom_in_use = chrom
    get_chromosome_seq(chrom)
    seq = extend_read(strand, chrom, pos)
    if seq:
        output_file.write(id_seq + "\n" + seq)

    #Process rest of file
    for line in input_file:
        (id_seq, strand, chrom, pos) = get_sam_info(line)
        check_chrom(chrom)
        seq = extend_read(strand, chrom, pos)
        if seq:
            output_file.write("\n" + id_seq + "\n" + seq)
    output_file.close()
