#!/usr/bin/python

import sys

args = sys.argv[1:]
if len(args) < 1:
   print "usage: input_file"
   sys.exit(1)


input_file_path = args[0]
output_file_path = input_file_path + ".srt"
del args[0]

edges_count = {}
edges = {}

with open(input_file_path, "r") as inf:
    with open(output_file_path, "w+") as outf:
        l = inf.readline()
        outf.write(l)
        #header written
        for l in inf:
            args = l.strip().split('\t')
            eid = int(args[0])
            edges_count[eid] = int(args[1])
            edges[eid] = l
        for e, count in sorted(edges_count.items(), key=lambda x: x[1], reverse=True):
            outf.write(edges[e])
