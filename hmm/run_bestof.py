#!/usr/bin/python

import sys
import subprocess

args = sys.argv[1:]

if len(args) < 4:
    print "usage: <file> n d ntry"
    sys.exit(1)

path = args[0]
n = int(args[1])
d = int(args[2])
ntry = int(args[3])

maxll = 0.0
plotcommandline = ""
for i in xrange(ntry):
    commandline = "./kmeans.R " + str(n) +  " " +  str(d) + " " + path + " | python hmm_1.py " + path + " " + str(n) + " " + str(d)
    train = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE)
    train.wait()
    # read from subprocess stdout the loglikelihood and the commandline to use
    ll_line = train.stdout.readline().strip()
    ll = float(ll_line)
    if ll > maxll:
        maxll = ll
        plotcommandline = train.stdout.readline().strip()
    sys.stdout.write("try: " + str(i) + "\tll: " + str(ll) + "\tmaxll: " + str(maxll) + "\n")

sys.stdout.write("tried " +  str(i) + " times\nmaxll: " +  str(maxll) + "\nrunning plot script\n")
plot = subprocess.Popen(plotcommandline, shell=True)
plot.wait()
