#!/usr/bin/python

import sys
import matplotlib.pyplot as plt
import math

args = sys.argv[1:]
if len(args) < 3:
   print "usage: [--plot] input_file control_file output_file"
   sys.exit(1)

plot = False
if args[0] == "--plot":
   plot = True
   del args[0]


input_file_path = args[0]
del args[0]
control_file_path = args[0]
del args[0]
output_file_path = args[0]
del args[0]


control_file = open(control_file_path, "r")
next(control_file) #skip header
#Create lookup table
lookup = {}
for line in control_file:
   (kmer, count, freq) = line.strip().split("\t")
   lookup[kmer] = (count, freq)

differential_kmer = {}
differential_kmer_log = {}
kmer_total = 0
diff_total = 0.0

input_file = open(input_file_path, "r")
next(input_file) #skip header
for line in input_file:
   (kmer, count, freq) = line.strip().split("\t")
   if not lookup.has_key(kmer): #skip if not present
      continue

   (ctrl_count, ctrl_freq) = lookup[kmer]
   diff = float(freq) / float(ctrl_freq)

   kmer_total = kmer_total + 1
   diff_total = diff_total + math.log(diff,2)

   differential_kmer[kmer] = (count, freq, ctrl_count, ctrl_freq, diff)
   differential_kmer_log[kmer] = math.log(diff, 2)

diff_mean = diff_total / kmer_total

if plot:
   fig, ax = plt.subplots()
   ax.bar(range(len(differential_kmer_log)), differential_kmer_log.values(), align='center', color='#cc3300', edgecolor='#cc3300')
   #Mean
   ax.plot((0,len(differential_kmer_log)), (diff_mean,diff_mean), color='#3366ff')
   ax.plot((0,len(differential_kmer_log)), (0,0), color='white')
   ax.axes.get_xaxis().set_visible(False)

   fig.savefig(output_file_path + ".pdf", format='pdf')


output_file = open(output_file_path + ".txt", "w")

output_file.write("k-mer\tcount\tfreq\tcontrol count\tcontrol freq\tdifference\tdifference log2\n")
for kmer, log_diff  in sorted(differential_kmer_log.items(), key=lambda x: x[1], reverse = True):
   (count, freq, ctrl_count, ctrl_freq, diff) = differential_kmer[kmer]
   output_file.write(str(kmer) + "\t" + str(count) + "\t" + str(freq) + "\t" + str(ctrl_count) + "\t" + str(ctrl_freq) + "\t" + str(diff) + "\t" + str(log_diff) + "\n")
