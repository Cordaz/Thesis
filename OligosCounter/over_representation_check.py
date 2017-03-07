#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np

threshold = 0.5

args = sys.argv[1:]
if len(args) < 3:
   print "usage: input_file control_file output_file"
   sys.exit(1)

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
differential_kmer = {}
kmer_total = 0
diff_total = 0.0

input_file = open(input_file_path, "r")
for line in input_file:
   (kmer, count, freq) = line.strip().split("\t")
   if lookup.has_key(kmer):
      (ctrl_count, ctrl_freq) = lookup[kmer]
      diff = (float(freq) - float(ctrl_freq)) / float(ctrl_freq)
   else:
      (ctrl_count, ctrl_freq) = (0,0)
      diff = 1
   
   kmer_total += 1
   
   diff_total += diff
   differential_kmer[kmer] = diff
   
   if abs(diff) >= threshold:
      over_represented[kmer] = (count, freq, ctrl_count, ctrl_freq)

diff_mean = diff_total / kmer_total

fig, ax = plt.subplots()
ax.bar(range(len(differential_kmer)), differential_kmer.values(), align='center', color='r', edgecolor='r')
#Mean
ax.plot((0,len(differential_kmer)), (diff_mean,diff_mean), color='b')
ax.axes.get_xaxis().set_visible(False)

fig.savefig(output_file_path + ".pdf", format='pdf')

output_file = open(output_file_path + ".txt", "w")

output_file.write("k-mer\tcount\tfreq\tcontrol count\tcontrol freq\tdifference\n")
for kmer in sorted(over_represented):
   (count, freq, ctrl_count, ctrl_freq) = over_represented[kmer]
   output_file.write(str(kmer) + "\t" + str(count) + "\t" + str(freq) + "\t" + str(ctrl_count) + "\t" + str(ctrl_freq) + "\t" + str(differential_kmer[kmer]) + "\n")

