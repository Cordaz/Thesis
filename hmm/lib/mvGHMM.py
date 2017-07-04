import numpy as np
import lib.safemath as safemath
import random as rand
from scipy.stats import multivariate_normal

'''

@author: Andrea Corneo

This code is mostly based on GuyZ HMM implementation
(https://bitbucket.org/GuyZ/hmm/)

Changes are:
 - a single class implementing only MultiVariate Gaussian HMM
 - Introduced use of safemath lib to avoid vanishing alpha

'''

class mvGHMM(object):
    '''
    MultiVariate Gaussian HMM

    Attributes:
    n           number of hidden states
    d           dimensionality
    A           transition probability matrix, [NxN] numpy array
    means       initial means of each gaussian, [NxD] numpy array
    covars      initial covariance matrix of each state, [Nx1] array of [DxD] numpy matrix
    pi          initial probability, [Nx1] numpy array

    Additional attributes:
    precision   numpy element denoting the precision
    fixedA      mantain initial A at each training iteration [boolean]
    fixedPi     mantain initial pi at each training iteration [boolean]
    verbose     print status after each training iteration [boolean]
    '''

    def __init__(self,n,d,A,means,covars,pi,precision=np.double,fixedA=False,fixedPi=False,verbose=False):
        self.n = n
        self.precision = precision
        self.verbose = verbose
        self.d = d
        self.A = A
        self.fixedA = fixedA
        self.pi = pi
        self.fixedPi = fixedPi
        self.means = means
        self.covars = covars


    # Used to compute once per iteration the vecotr B_j(t),
    # that is the probability of observing value O_t from state j at time t

    def calcB(self,observ):
        self.B = np.zeros( (self.n, len(observ)), dtype=self.precision )

        for j in xrange(self.n):
            for t in xrange(len(observ)):
                self.B[j][t] = self.pdf(observ[t], self.means[j], self.covars[j])


    # Update current model attributes according to a new model

    def updatemodel(self,new_model):
        if not self.fixedA:
            self.A = new_model['A']
        if not self.fixedPi:
            self.pi = new_model['pi']
        self.means = new_model['means']
        self.covars = new_model['covars']


    # Compute all variables used for training
    # Perform the E-step of the Baum-Welch

    def calcvariables(self, observ):
        variables = {}

        variables['alpha'] = self.calcalpha(observ)
        variables['beta'] = self.calcbeta(observ)
        variables['xi'] = self.calcxi(observ,variables['alpha'],variables['beta'])
        variables['gamma'] = self.calcgamma(variables['alpha'],variables['beta'],len(observ))

        return variables


    # Helper method of calcvariable: compute the (logarithm) value of alpha (forward variable).

    def calcalpha(self, observ):
        alpha = np.ones((len(observ),self.n),dtype=self.precision) # allocation

        # init stage - alpha_1(x) = pi(x)b_x(O1)
        for x in xrange(self.n):
            alpha[0][x] = safemath.safelnprod(safemath.safeln(self.pi[x]), safemath.safeln(self.B[x][0]))

        # induction
        for t in xrange(1,len(observ)):
            for j in xrange(self.n):
                logalpha = safemath.LOGZERO
                for i in xrange(self.n):
                    logalpha = safemath.safelnsum( logalpha, safemath.safelnprod(alpha[t-1][i], safemath.safeln(self.A[i][j])) )
                alpha[t][j] = safemath.safelnprod(logalpha, safemath.safeln(self.B[j][t]))

        return alpha


    # Helper method of calcvariable: compute the (logarithm) value of beta (backward variable).

    def calcbeta(self, observ):
        beta = np.zeros((len(observ),self.n),dtype=self.precision)

        # induction
        for t in xrange(len(observ)-2,-1,-1):
            for i in xrange(self.n):
                logbeta = safemath.LOGZERO
                for j in xrange(self.n):
                    logbeta = safemath.safelnsum( logbeta,  safemath.safelnprod(
                        safemath.safeln(self.A[i][j]), safemath.safelnprod(
                            safemath.safeln(self.B[j][t+1]), beta[t+1][j]
                            )
                        )
                    )
                beta[t][i] = logbeta

        return beta


    # Helper method of calcvariable: compute the (logarithm) value of xi (Baum-Welch parameters update variable).

    def calcxi(self,observ,alpha=None,beta=None):
        if alpha is None:
            alpha = self._calcalpha(observ)
        if beta is None:
            beta = self._calcbeta(observ)
        xi = np.zeros((len(observ),self.n,self.n),dtype=self.precision)

        for t in xrange(len(observ)-1):
            normalizer = safemath.LOGZERO
            for i in xrange(self.n):
                for j in xrange(self.n):
                    xi[t][i][j] = safemath.safelnprod( alpha[t][i],
                        safemath.safelnprod( safemath.safeln(self.A[i][j]),
                            safemath.safelnprod( safemath.safeln(self.B[j][t+1]), beta[t+1][j] )
                        )
                    )
                    normalizer = safemath.safelnsum(normalizer, xi[t][i][j])

            for i in xrange(self.n):
                for j in xrange(self.n):
                    xi[t][i][j] = safemath.safelnprod(xi[t][i][j], -1*normalizer)

        return xi


    # Helper method of calcvariable: compute the (logarithm) value of gamma (state probabilities variable).

    def calcgamma(self, alpha, beta, seqlen):
        gamma = np.zeros((seqlen,self.n),dtype=self.precision)

        for t in xrange(seqlen):
            normalizer = safemath.LOGZERO
            for i in xrange(self.n):
                gamma[t][i] = safemath.safelnprod(alpha[t][i], beta[t][i])
                normalizer = safemath.safelnsum(normalizer, gamma[t][i])
            for i in xrange(self.n):
                gamma[t][i] = safemath.safelnprod(gamma[t][i], -1*normalizer)

        return gamma


    # Performs the M-step of the Baum-Welch algorithm

    def reestimate(self,variables,observ):
        new_model = {}

        if not self.fixedPi:
            new_model['pi'] = self.reestimatePi(variables['gamma'][0])
        if not self.fixedA:
            new_model['A'] = self.reestimateA(observ,variables['xi'],variables['gamma'])
        new_means, new_covars = self.reestimateGaussian(observ,variables['gamma'])
        new_model['means'] = new_means
        new_model['covars'] = new_covars

        return new_model


    # Helper method of reestimate: reestimate initial probability pi
    # as probability of being in state x at time t=0

    def reestimatePi(self, gamma0):
        new_pi = np.zeros( self.n, dtype=self.precision )

        for x in xrange(self.n):
            new_pi[x] = safemath.safeexp( gamma0[x] )

        return new_pi


    # Helper method of reestimate: reestimate matrix A

    def reestimateA(self,observ,xi,gamma):
        new_A = np.zeros((self.n, self.n), dtype=self.precision)
        for i in xrange(self.n):
            for j in xrange(self.n):
                numer = safemath.LOGZERO
                denom = safemath.LOGZERO
                for t in xrange(len(observ)-1):
                    numer = safemath.safelnsum( numer, xi[t][i][j] )
                    denom = safemath.safelnsum( denom, gamma[t][i] )
                new_A[i][j] = safemath.safeexp( safemath.safelnprod(numer, -1*denom) )

        return new_A


    # Helper method of reestimate: reestimate means and covars

    def reestimateGaussian(self,observ,gamma):
        new_means = np.zeros( (self.n,self.d), dtype=self.precision )
        new_covars = [ np.matrix(np.zeros((self.d,self.d), dtype=self.precision)) for i in xrange(self.n) ]

        for j in xrange(self.n):
            numer = [0.0] * self.d
            denom = 0.0
            for t in xrange(len(observ)):
                numer += safemath.safeexp(gamma[t][j]) * observ[t]
                denom += safemath.safeexp(gamma[t][j])
            new_means[j] = numer/denom

        #covars_prior = [ np.matrix(0.001*np.eye(self.d, dtype=self.precision)) for i in xrange(self.n) ] # in case of underflowing
        for j in xrange(self.n):
            numer = np.matrix(np.zeros( (self.d,self.d), dtype=self.precision))
            denom = np.matrix(np.zeros( (self.d,self.d), dtype=self.precision))
            for t in xrange(len(observ)):
                vector_as_mat = np.matrix( (observ[t]-self.means[j]), dtype=self.precision )
                numer += safemath.safeexp(gamma[t][j]) * np.dot(vector_as_mat.T, vector_as_mat)
                denom += safemath.safeexp(gamma[t][j])
            new_covars[j] = numer/denom
            #new_covars[j] += new_covars[j] + covars_prior[j] # in case of underflowing

        return new_means, new_covars


    # Compute the PDF of value x given N(mean,covar) (multivariate gaussian)

    def pdf(self,x,mean,covar):
        dist = multivariate_normal(mean=mean, cov=covar)
        return dist.pdf(x)


    # Compute the probability of the observation given the model,
    # return the log of the probability (that is the loglikelihood)

    def forwardbackward(self,observ,cache=False):
        if not cache:
            self.calcB(observ)

        alpha = self.calcalpha(observ)

        loglikelihood = safemath.LOGZERO
        for i in xrange(self.n):
            loglikelihood = safemath.safelnsum( loglikelihood, alpha[-1][i] )
        return loglikelihood


    # Find the path that maximizes the loglikelihood.
    # Implementation of Viterbi algorithm.
    #
    # delta[t][i] = max(P[q1..qt=i,O1...Ot|model] - the path ending in Si and until time t,
    # that generates the highest probability.
    #
    # psi[t][i] = argmax(delta[t-1][i]*aij) - the index of the maximizing state in time (t-1),
    # i.e: the previous state.

    def viterbi(self,observ):
        self.calcB(observ)

        delta = np.zeros((len(observ),self.n), dtype=self.precision)
        psi = np.zeros((len(observ),self.n), dtype=np.int32)

        for t in xrange(len(observ)):
            for j in xrange(self.n):
                delta[t][j] = safemath.LOGZERO
                psi[t][j] = 0

        # init
        for x in xrange(self.n):
            delta[0][x] = safemath.safelnprod(safemath.safeln(self.pi[x]), safemath.safeln(self.B[x][0]))

        # induction
        for t in xrange(1,len(observ)):
            for j in xrange(self.n):
                for i in xrange(self.n):
                    if safemath.safeislt( delta[t][j], ( safemath.safelnprod(delta[t-1][i], safemath.safeln(self.A[i][j])) ) ):
                        delta[t][j] = safemath.safelnprod(delta[t-1][i], safemath.safeln(self.A[i][j]))
                        psi[t][j] = i
                delta[t][j] = safemath.safelnprod(delta[t][j], safemath.safeln(self.B[j][t]))

        # find max prob in time T
        p_max = safemath.LOGZERO
        path = np.zeros(len(observ), dtype=np.int32)
        for i in xrange(self.n):
            if safemath.safeislt( p_max, delta[len(observ)-1][i] ):
                p_max = delta[len(observ)-1][i]
                path[len(observ)-1] = i

        # backtracking
        for i in xrange(1, len(observ)):
            path[len(observ)-i-1] = psi[len(observ)-i][path[len(observ)-i]]

        return path


    # Updates the HMMs parameters given a new set of observed sequences.

    def train(self,observ,iterations=1,eps=0.0001,threshold=-0.001):
        self.calcB(observ)
        for i in xrange(iterations):
            old_p, new_p = self.performiter(observ)
            if self.verbose:
                print "iter: ", i, ", L(model|O)= ", old_p, ", L(new_model|O)= ", new_p, ", converging= ", ( new_p-old_p > threshold )

            if abs(new_p-old_p) < eps:
                break # converged


    # A single iteration of an EM algorithm, which given the current HMM,
    # computes new model parameters and internally replaces the old model
    # with the new one.
    #
    # Returns the log likelihood of the old model (before the update),
    # and the one for the new model.

    def performiter(self,observ):
        # start EM algorithm
        new_model = self.baumwelch(observ)

        # calculate old loglikelihood
        old_p = self.forwardbackward(observ, cache=True)

        # update model
        self.updatemodel(new_model)

        # calculate new loglikelihood
        new_p = self.forwardbackward(observ, cache=False)

        return old_p, new_p


    # An EM(expectation-maximization) algorithm devised by Baum-Welch. Finds a local maximum
    # that outputs the model that produces the highest probability, given a set of observations.
    #
    # Returns the new maximized model parameters

    def baumwelch(self, observ):
        #E-step
        variables = self.calcvariables(observ)

        #M-step
        return self.reestimate(variables,observ)
