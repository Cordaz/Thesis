import lib.mvGHMM as hmm
import sys
import subprocess
import numpy as np
import random as rand

# Randomize observation vector

def randomizeobserv(observ):
    obs = observ.copy()

    rand.shuffle(obs)

    return obs


### args
args = sys.argv[1:]
path = args[0]
n = int(args[1])
d = int(args[2])
m = 2 #num of mixture components

### read observ
obs = []
with open(path, "r") as in_f:
    header = in_f.readline() # read and store header
    for l in in_f:
        s = l.strip().split('\t')
        obs.append( [ float(s[3]) , float(s[6]) ] )

observ = np.asarray(obs) #convert in a [TxD] numpy array
obs = randomizeobserv(observ)

### setting up model
## parameters

#define precision
precision = np.double

#transition matrix
A = np.ones((n,n), precision) / float(n) # uniform probability for each state

# initial probability
pi = np.ones(n, precision) / float(n) # uniform probability for each state


# means (coming from kmeans.R via stdin)
means = np.zeros( (n,d), dtype=precision )
l = raw_input()
s = l.split(" ")
for x in xrange(n):
    for i in xrange(d):
        means[x][i] = float(s[0])
        del s[0]

# covariances (coming from kmeans.R via stdin)

covars = [ np.matrix(np.zeros((d,d), precision)) for i in xrange(n)]
l = raw_input()
s = l.split(" ")
for x in xrange(n):
    for i in xrange(d):
        for j in xrange(d):
            while not s[0]:
                del s[0] # remove double white spaces
            covars[x][i,j] = float(s[0])
            del s[0]
'''


# coming from mixture
means = np.zeros( (n,d), dtype=precision )
means[0][0] = 0.7577889
means[0][1] = 0.7685576
means[1][0] = 0.9225547
means[1][1] = 0.9717965
means[2][0] = 0.7577889
means[2][1] = 0.9717965
means[3][0] = 0.9225547
means[3][1] = 0.7685576

covars = [ np.matrix(np.zeros((d,d), precision)) for i in xrange(n)]
covars[0][0,0] = 0.04501044
covars[0][1,1] = 0.06785575
covars[1][0,0] = 0.03627211
covars[1][1,1] = 0.02172594
covars[2][0,0] = 0.04501044
covars[2][1,1] = 0.02172594
covars[3][0,0] = 0.03627211
covars[3][1,1] = 0.06785575
'''

## build model
model = hmm.mvGHMM(n,d,A,means,covars,pi,precision,fixedA=True,fixedPi=True,verbose=False)


### training
model.train(observ, iterations=1000, threshold=-0.000001)
'''
# using chunks
chunks = [observ[i:i + 100] for i in xrange(0, len(observ), 100)]
tot = len(observ) / 100 + 1
i = 1
for chunk in chunks:
    model.train(chunk, iterations=1000, threshold=-0.000001)
    sys.stdout.write("chunk %d of %d\n" % (i, tot))
    i += 1
'''

# for c in model.covars:
#    print c

'''
### outputs
# print means
for i in xrange(n-1):
    for j in xrange(d):
        sys.stdout.write(str(model.means[i][j]))
        sys.stdout.write(",")
for j in xrange(d-1):
    sys.stdout.write(str(model.means[n-1][j]))
    sys.stdout.write(",")
sys.stdout.write(str(model.means[n-1][d-1]))
sys.stdout.write("\n")


# print covars
for i in xrange(n-1):
    for j in xrange(d):
        sys.stdout.write(str(model.covars[i].item(j,j)))
        sys.stdout.write(",")
for j in xrange(d-1):
    sys.stdout.write(str(model.covars[n-1].item(j,j)))
    sys.stdout.write(",")
sys.stdout.write(str(model.covars[n-1].item(d-1,d-1)))
sys.stdout.write("\n")

# print predicted
predicted = model.viterbi(observ)
for t in xrange(len(predicted)-1):
    sys.stdout.write(str(predicted[t]))
    sys.stdout.write(",")
sys.stdout.write(str(predicted[len(predicted)-1]))
sys.stdout.write("\n")
'''

### plot
# command line
commandline = "./plot.R " + path + " " + str(n) + " "
# add means
for i in xrange(n-1):
    for j in xrange(d):
        commandline += str(model.means[i][j]) + ","
for j in xrange(d-1):
    commandline += str(model.means[n-1][j]) + ","
commandline += str(model.means[n-1][d-1]) + " "

# add covars
for i in xrange(n-1):
    for j in xrange(d):
        commandline += str(model.covars[i].item(j,j)) + ","
for j in xrange(d-1):
    commandline += str(model.covars[n-1].item(j,j)) + ","
commandline += str(model.covars[n-1].item(d-1,d-1)) + " "

# add predicted states
predicted = model.viterbi(observ)
for t in xrange(len(predicted)-1):
    commandline += str(predicted[t]) + ","
commandline += str(predicted[len(predicted)-1]) + " "

print model.forwardbackward(observ)
print commandline

# plot_proc = subprocess.Popen(commandline, shell=True)
# plot_proc.wait()
