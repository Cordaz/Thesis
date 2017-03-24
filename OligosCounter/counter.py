#!/usr/bin/python

import sys

k = 8 #k-mer length

def drop_newline(s):
   return s[0:len(s)-1]


def get_reads_fq(f):
   reads = []
   i = 0
   for line in f:
      if i == 4:
         i = 0
      i = i + 1
      if i != 2:
         continue
      reads.append(drop_newline(line))

   return reads

def get_reads_fa(f):
   reads = []
   i = 0
   for line in f:
      if i == 2:
          i = 0
      i = i + 1
      if i == 1:
          continue
      reads.append(drop_newline(line))

   return reads


def get_kmers(read):
   kmers = []
   for i in range(0, len(read)-k+1):
      kmers.append(read[i:i+k])

   return kmers


def is_valid_kmer(kmer):
   return 'N' not in kmer


def complementary_kmer(kmer):
   comp_bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
   c_kmer = ""
   for char in reversed(kmer):
      c_kmer = c_kmer + comp_bases[char]

   return c_kmer

def is_palyndrome(kmer, c_kmer):
   return kmer == c_kmer


args = sys.argv[1:]
if not args or len(args) != 3:
   print "usage: [--fasta|--fastq] input_file output_file";
   sys.exit(1)

fasta = False
fastq = False
if args[0] == "--fasta":
    fasta = True
    del args[0]
elif args[0] == "--fastq":
    fastq = True
    del args[0]

if not (fasta or fastq):
    print "select input type (--fasta or --fastq)"
    sys.exit(1)

input_file_path = args[0]
del args[0]
output_file_path = args[0]
del args[0]

kmers_hash = {}
kmer_total = 0

with open(input_file_path, "r") as input_file:
    if fasta:
        reads = get_reads_fa(input_file)
    elif fastq:
        reads = get_reads_fq(input_file)

    for read in reads:
       kmers = get_kmers(read)

       for kmer in kmers:
          if is_valid_kmer(kmer):
             if kmer in kmers_hash:
                kmers_hash[kmer] = kmers_hash[kmer] + 1
             else:
                kmers_hash[kmer] = 2
                kmer_total = kmer_total + 1
             kmer_total = kmer_total + 1
             c_kmer = complementary_kmer(kmer)
             if not is_palyndrome(kmer, c_kmer):
                if c_kmer in kmers_hash:
                   kmers_hash[c_kmer] = kmers_hash[c_kmer] + 1
                else:
                   kmers_hash[c_kmer] = 2
                   kmer_total = kmer_total + 1
                kmer_total = kmer_total + 1

output_file = open(output_file_path, "w")
output_file.write("k-mer\tcount\tfreq\n")
output_file.write("@total\t" + str(kmer_total) + "\n")
for kmer in sorted(kmers_hash):
   output_file.write(str(kmer) + "\t" + str(kmers_hash[kmer]) + "\t" + str(float(kmers_hash[kmer]) / kmer_total) + "\n")
output_file.close()
