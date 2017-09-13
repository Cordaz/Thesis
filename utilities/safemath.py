import math

'''

#author Andrea Corneo

This library extend the standard math lib to include the management of ln(0)

Based on:
 - Numerically Stable Hidden Markov Model Implementation,
   Tobias P. Mann, February 21, 2006
   http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf

'''

# define a known value for ln(0)
LOGZERO = float('nan')

# return the exponential of x or 0 is x is LOGZERO
def safeexp(x):
    if math.isnan(x):
        return 0
    return math.exp(x)

# return ln(x) if x>0, LOGZERO if x=0, raise exception if x<0
def safeln(x):
    if x < 0:
        raise ValueError('Cannot compute logarithm of negative number')
    if x == 0:
        return LOGZERO
    return math.log(x)

# compute the sum of two safe-logarithm
def safelnsum(safelnx, safelny):
    if math.isnan(safelnx) and math.isnan(safelny):
        return LOGZERO
    if math.isnan(safelnx):
        return safelny
    if math.isnan(safelny):
        return safelnx
    if safelnx > safelny:
        return safelnx + safeln(1 + math.exp(safelny-safelnx))
    return safelny + safeln(1 + math.exp(safelnx-safelny))

# compute the product of two safe-logarithm
def safelnprod(safelnx, safelny):
    if math.isnan(safelnx) or math.isnan(safelny):
        return LOGZERO
    return safelnx + safelny

# compare two safeln: greater than
def safeisgt(safelnx, safelny):
    if math.isnan(safelnx) and math.isnan(safelny):
        return False
    if math.isnan(safelnx):
        return False
    if math.isnan(safelny):
        return True
    return safelnx > safelny

# compare two safeln: less than
def safeislt(safelnx, safelny):
    if math.isnan(safelnx) and math.isnan(safelny):
        return False
    if math.isnan(safelnx):
        return True
    if math.isnan(safelny):
        return False
    return safelnx < safelny

# compare two safeln: equal
def safeiseq(safelnx, safelny):
    if math.isnan(safelnx) and math.isnan(safelny):
        return True
    if math.isnan(safelnx) or math.isnan(safelny):
        return False
    return safelnx == safelny
