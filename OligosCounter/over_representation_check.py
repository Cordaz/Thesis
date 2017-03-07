#!/usr/bin/python

import sys

args = sys.argv[1:]
if not args:
   print "usage: [-t threshold (default 0.1 (10%))] input_file control_file output_file"
   sys.exit(1)

threshold = 0.1
if args[0] == "-t":
   threshold = args[1]
del[0:2]

if len(args) < 3:
   print "usage: [-t threshold (default 0.1 (10%))] input_file control_file output_file"

input_file_path = args[0]
del args[0]
control_file_path = args[0]
del args[0]
output_file_path = args[0]
del args[0]


control_file = open(control_file_path, "r")
#Create lookup table
lookup = {}
for line in control_file:
   (kmer, count, freq) = line.strip().split("\t")
   lookup[kmer] = (count, freq)

over_represented = {}

input_file = open(input_file_path, "r")
for line in input_file:
   (kmer, count, freq) = line.strip().split("\t")
   if lookup.has_key(kmer):
      (ctrl_count, ctrl_freq) = lookup[kmer]
   else:
      (ctrl_count, ctrl_freq) = (0,0)
   
   if freq >= (1 + threshold) * ctrl_freq:
      over_represented[kmer] = (count, freq, ctrl_count, ctrl_freq)

output_file = open(output_file_path, "w")

output_file.write("k-mer\tcount\tfreq\tcontrol count\tcontrol freq\n")
for kmer in sorted(over_represented):
   (count, freq, ctrl_count, ctrl_freq) = over_represented[kmer]
   output_file.write(str(kmer) +  "\t" + str(count) + "\t" + str(freq) + "\t" + str(ctrl_count) + "\t" + str(ctrl_freq) + "\n")
