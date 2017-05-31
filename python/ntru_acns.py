#!/usr/bin/python
# written for python 2.7, broken on python 3.x at this point

# uses mpmath
# http://docs.sympy.org/0.7.6/modules/mpmath/setup.html
# https://github.com/fredrik-johansson/mpmath/releases
# using mpmath-0.19, download the .tgz, untar it, and:
# cd ~/Downloads/mpmath-0.19
# sudo python setup.py install
# test by..
#import mpmath
#mpmath.runtests()


# implementation based on IEEE P1363.1/D12, October 2008

# A special note:
#
# 12 The attack thus has two stages: the lattice reduction stage and the combinatorial stage. The total time for the
# 13 attack is the sum of the time for these stages. This standard requires that for a security level k, both of these
# 14 stages shall take more than k bits of work.

import unittest
import math 		#turns to floats and limits to double precision
import sys
import getopt
import json

from time import gmtime, strftime

# used for multiprocessing (process is distinct from thread in python..)
#from Queue import Queue
import multiprocessing #from multiprocessing import *  #Pool
from multiprocessing import Process, Queue, freeze_support
#from multiprocessing import *

import operator as op

#from multiprocessing import freeze_support
# mpmath for arbitrary precision
from mpmath import *
#mp.prec = 80 		# bit precision, default works fine

# for testing if a file exists, to skip work if it's already been done.
import os.path

# use basic floats instead of mpf, unless absolutely necessary
def mpf_necessary(a):
	return mpf(a)

def mpf_nice(a):
	return float(a)

# output file prefix, put the final output jsons here with our naming convention
results_output_file_prefix = "output/"

# output detailed logs to this place. I was unable to use absolute paths in my environment..
# TODO the code currently will always output these detailed logs. This must point to a valid
# directory or else the code will probably crash for you..
detailed_log_file_prefix = "../../../Desktop/detailed_ntru_logs/"

# default bit security
k = 112

# set q to 2048 as per the paper
constant_q = 2048

# regulate debug output
debug = 1

# enable or disable caching
caching_enabled = True

# how fast we move forward through the alpha
# determined through heuristics
alpha_step = 15

# number of concurrent processes to run
single_threaded = False
#single_threaded = True

# Number of cores to use for processing
# I recommend total cores (threads) on your system - 1 or - 2
# Using all threads on your system will probably make your system unresponsive
wmitm_parallelism = 15

# for the binomial coefficient caching..
binomial_cache_hits = 0
binomial_cache_misses = 0
binomial_cache_dict = {}

# for the multinomial coefficient caching..
multinomial_cache_hits = 0
multinomial_cache_misses = 0
multinomial_cache_dict = {}

# cache logs
LOG = [mpf_nice(0.0)]
# update LOG array
for i in range (1, 2049):
	LOG.append(mpf_nice(math.log(mpf_nice(i), mpf_nice(2))))



########################################
# calculations for Wmitm


def calculateSigma0(N, d):
	return mpf_nice(math.sqrt(mpf_nice(8*d)/mpf_nice(3*N)))


def calculateSigma1(N, d1, d2, c1, c2, y2):
	# sigma
	dg = N // 3
	# This functions correctly.
	sigma = math.sqrt(mpf_nice(2 * dg + (d1 + d2) - (c1 + c2))/y2)
	return mpf_nice(sigma)


def calculatePsFunc(D, sigma):
	# Where: f(D,sigma) =
	# 	 erfc(D / sigma*sqrt(2)) - |(sigma*sqrt(2))/(D*sqrt(pi))|(e^(-1 * (D^2)/(2*sigma^2)) - 1)
	arg = (mpf_nice(D) / (mpf_nice(sigma) * mpf_nice(math.sqrt(2.0))))
	left_term = mpf_nice(math.erfc(arg))
	right_term = mpf_nice(mpf_nice(math.exp(-arg*arg)) - mpf_nice(1.0))/(mpf_nice(arg) * mpf_nice(math.sqrt(math.pi)))

	return mpf_nice(left_term - right_term)

def calculatePs(q, N, y2, alpha, sigma):
	# calculating Ps without y1, 2nd equation from the 2008/2009 paper
	# ps =
	# (1 - (2 / 3q)) ^ (((2N - y2(1 + alpha))) / (1-alpha))    times..
	# Pi: i = 0 to (2(y2 - N) / 1-alpha) for:
	#		(1 - f|(q^(2*alpha(y2-N)+i(1-alpha)^2)), sigma| )
	# Where: f(D,sigma) =
	#	 erfc(D / sigma*sqrt(2)) - |(sigma*sqrt(2))/(D*sqrt(pi))|(e^(-1 * (D^2)/(2*sigma^2)) - 1)
	first_term = (mpf_nice(1) - (mpf_nice(2.0) / mpf_nice(3.0 * mpf_nice(q))))
	first_term_exp = mpf_nice((2.0 * N) - (y2 * (1.0 + alpha))) / mpf_nice(1.0 - alpha)
	complete_first_term = pow(first_term, first_term_exp)

	loop_range_start = 0
	loop_range_finish = mpf_nice(2.0 * (y2 - N)) // mpf_nice(1.0 - alpha)

	running_total = mpf_nice(1.0)

	for i in xrange(loop_range_start, int(loop_range_finish) + 1):
		# f(q^((alpha(y2-y1)+i(1-alpha)) / y2-y1)
		# f(q^((2alpha(y2-N)+i(1-alpha)^2) / 2(y2-N))
		func_term = pow(mpf_nice(q), (mpf_nice(2.0 * alpha * (y2 - N) + i * mpf_nice(1.0 - alpha) * mpf_nice(1.0 - alpha)) / mpf_nice(2.0 * (y2 - N))))
		term = mpf_nice(1.0) - (calculatePsFunc(func_term, sigma))
		running_total *= term

	result = complete_first_term * running_total
	return result

# Wsearch
def calculateWsearch(K, c1, c2):
	# Wsearch = [note, optimized from previous version w/sqrt on denominator due to bday atk]
	# |    K   || K - c1 / 2 |
	# | c1 / 2 ||   c2 / 2   |
	# ------------------------
	# |   c1   ||   c2   |
	# | c1 / 2 || c2 / 2 |

	topLeft = nchooser(K, c1 // 2)
	topRight = nchooser(K - c1 // 2, c2 // 2)
	bottomLeft = nchooser(c1, c1 // 2)
	bottomRight = nchooser(c2, c2 // 2)

	if (topLeft == None) or (topRight == None) or (bottomLeft == None) or (bottomRight == None):
		return None

	top = topLeft + topRight
	bottom = bottomLeft + bottomRight

	return mpf_nice(pow(mpf_nice(2), mpf_nice(top) - mpf_nice(bottom)))

# the 'conservative' version
def MITMTrialsBeforeCollision(q, N, c, y2, alpha, sigma):
	# this is aka 'Wsearch'
	# N' meaning conservative version is: N0 / Ps
	# N0:
	# [ 2N-y2   ] ([    c    ] )^{-1}
	# [ c/2,c/2 ] ([ c/2, c/2] )
	# (multinomial coefficient)
	#K = 2 * N - y2
	y_term = calculateWsearch((2 * N) - y2, c, c)
	inner_ps_term = mpf_nice(calculatePs(q, N, y2, alpha, sigma))
	if (inner_ps_term == None):
		return None
	ps_term = pow(mpf_nice(2), -1.0 * math.log(inner_ps_term, 2))
#	print "Y: ", math.log(y_term, 2), ", ps: ", math.log(ps_term, 2)

	if(y_term == None) or (ps_term == None):
		return None

	result = mpf_nice(y_term) * mpf_nice(ps_term)
#	print "MITMTrials alpha: ", alpha, " is: ", result, " ; log2: ", math.log(result, 2), ", ps: ", ps_term, ", ps log2: ", math.log(ps_term, 2)

	return result


def calculatePsplit(N, K, d1, d2, c1, c2, whichPsplitIteration):
	# Robert: Verified 4/26 that this is equivalent to 2009 paper Psplit
	# Note, whichPsplitIteration is the ",1" or ",N".
	# Psplit,1 =
	# |N  - K ||N  - K  - (d1-c1)| * | K  || K - c1 |
	# |d1 - c1||     d2 - c2     |   | c1 ||   c2   |
	# -----------------------------------------------
	#		| N  || N - d1 |
	#    		| d1 ||   d2   |
	# (remember, x choose y = x!/y!(x-y)!)  (math.factorial(x))

	if (debug >= 5):
		print "calculatePsplit: N:", N, ", K:, ", K, ", d1:", d1, ", d2:", d2, ", c1:", c1, ", c2:", c2, ", whichIt:", whichPsplitIteration

	topLeftmost = nchooser(N - K, d1 - c1)
	topLeft = nchooser(N - K - (d1 - c1), d2 - c2)
	topRight = nchooser(K, c1)
	topRightmost = nchooser(K - c1, c2)

	if (topLeftmost == None) or (topLeft == None) or (topRight == None) or (topRightmost == None):
		return None

	if (debug >= 5):
		print "topLeftmost: ", topLeftmost, ", topLeft: ", topLeft, ", topRight: ", topRight, ", topRightmost: ", topRightmost

	top = topLeftmost + topLeft + topRight + topRightmost

	bottomLeft = nchooser(N, d1)
	bottomRight = nchooser(N - d1, d2)

	if (bottomLeft == None) or (bottomRight == None):
		return None

	bottom = bottomLeft + bottomRight

	PsplitOne = mpf_nice(top) - mpf_nice(bottom)

	ret = mpf_nice(pow(mpf_nice(2), PsplitOne))

	if whichPsplitIteration == 1:
		return ret

	# otherwise..
	# Psplit,N = 1 - (1 - Psplit,1)^N
	# Note that N is prime so it's always going to be odd
	c = mpf_nice(1)
	p = mpf_nice(-1)
	PsplitN = mpf_nice(0)
	# this is the only place where we actually need high precision numbers, because this requires precision or becomes nan
	for i in range(0,N/2):
		c *= mpf_necessary(N - i)/mpf_necessary(i + 1)
		p *= -ret
		PsplitN += c*p

	return PsplitN


# Wreduction
def calculateWreduction(y2):
	# Wreduction = N / 2^1.06 [NOTE: the value N^2/2^1.06 is also given, but 'optimized']
	return mpf_nice(y2) / mpf_nice(pow(2.0, 1.06))



def MITMRunningTime(q, y2, d, c, N, K, alpha, sigma):
	# = tN/Psplit   where the Psplit appears to use N instead of 1
	# t is Wreduction
	t_term = calculateWreduction(y2)
	n_term = MITMTrialsBeforeCollision(q, N, c, y2, alpha, sigma)
	# n_term is returned as 2^values..
	# using , N is so that we are using the conservative (not current) estimate
	psplit_term = calculatePsplit(N, K, d, d, c, c, N)
#	print "psplit at alpha: ", alpha, " is: ", psplit_term
#	print "t = ", math.log(t_term, 2), ", psplit: ", math.log(psplit_term, 2)

	# Invalid part of the profile
	if(t_term == None) or (n_term == None) or (psplit_term == None):
		return None

#	print "about to do log, t_term: ", t_term, ", n_term: ", n_term, ", psplit_term: ", psplit_term
	if (psplit_term <= 0) or (psplit_term > 1):
		return None

	return math.log((t_term * n_term) / psplit_term, 2)


def LatticeRunningTime(N, y2, alpha):
	# Assumption #4 from Hirschorn 2009 paper sets c and m, conservative not current estimates
	local_c_term = -50
	local_m_term = 0.2
	# the other published terms, reflecting 'current' lattice attacks.
	# These don't match the paper's published results, since those are based on
	# the assumption that there are better attacks that exist
#	local_c_term = -110
#	local_m_term = 0.45

	w_left_term = (2.0 * local_m_term * (y2 - N)) / ((1 - alpha)*(1 - alpha))
	w_mid_term = 3 * log((2 * (y2 - N)) / (1 - alpha))
	w_right_term = local_c_term

	w = w_left_term + w_mid_term + w_right_term

	# Running time is 2^w where w is as above
	return math.log(pow(2, w), 2)



########################################
# Decryption failure
#
# this equation comes from A hybrid lattice-reduction and meet-in-the-middle attack against NTRU
# Nick Howgrave-Graham, pg16...
# alpha = 2N - m - r / m - r, we thought that m and r were y2 and y1
def calculateDecryptFailure(q, N, d):

	sigma = calculateSigma0(N, d)
	p = mpf_nice(N) * math.erfc( mpf_nice(q - 2)/mpf_nice((6.0*math.sqrt(2.0*N)*sigma)) )

	if(p > 0):
		return -int(math.log(p, 2))
	else:
		return 0


# Algorithm from 2009 ACNS paper, 'Algorithm 2'
def calculateForGivenParameters(q, N, d):
	aLog = open_unique_detailed_logfile(q, N, d)

	if(debug >= 1):
		print "calculateForGivenParameters(N=", N, ", q=", q, ", d=", d, ")"

	if (results_file_exists(q, N, d)):
		cachedResults = results_from_file(q, N, d)
		if(debug >= 1):
			print "Using cached results for q:", q, ", N:", N, ", d:", d, " = ", json.dumps(cachedResults, indent = 4)
		aLog.write("Using cached results for q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + "\n")
		aLog.close()
		return cachedResults

	# best w calc for a given q, N, d
	bestWork = None

	# Calculate decryption faiure
	decryption_failure = calculateDecryptFailure(q, N, d)
	if(debug >= 2):
		print "        For N: ", N, ", d: ", d, "df: ", decryption_failure
	aLog.write("decryption failure calculated as: " + str(decryption_failure) + "\n")

	# Go through the y2 range
	for y2 in range(N + 1, 2 * N + 1, 1):
		for c in range(1, d + 1, 1):
			alpha = findAlpha(q, N, d, c, y2)
			if (alpha == None):
				if (debug >= 5):
					print "Skipping q: ", q, ", N:", N, ", d:", d, ", y2:", y2, ", c:", c, " because of invalid profile"
				aLog.write("Skipping q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + ", y2: " + str(y2) + ", c: " + str(c) + " because of invalid profile\n");
				continue # skip this q, N, d, y2, c because it's an invalid profile
			sigma = mpf_nice(calculateSigma1(N, d, d, c, c, y2))
			mitm = MITMRunningTime(q, y2, d, c, N, 2 * N - y2, alpha, sigma)
			lattice = LatticeRunningTime(N, y2, alpha)
			if (mitm == None) or (lattice == None): # NOTE, probably wouldn't get here, would catch above..
				if (debug >= 5):
					print "Skipping q: ", q, ", N:", N, ", d:", d, ", y2:", y2, ", c:", c, " because of invalid profile"
				aLog.write("Skipping q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + ", y2: " + str(y2) + ", c: " + str(c) + " because of invalid profile\n");
				continue # skip this q, N, d, y2, c because it's an invalid profile
			workCandidate = max(lattice, mitm) # take the largest. they should be nearly equal
			if (debug >= 4):
				print "Processed q: ", q, ", N:", N, ", d:", d, ", y2:", y2, ", c:", c, ", result:", workCandidate, ", alpha:", alpha
			aLog.write("Processed q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + ", y2: " + str(y2) + ", c: " + str(c) + ", result: " + str(workCandidate) + ", alpha: " + str(alpha) + "\n")
			if((bestWork == None) or (bestWork > workCandidate)):
				bestWork = workCandidate
				bestParameters = { 'securitylevel':bestWork, 'N':N, 'd':d, 'q':q, 'c':c, 'y2':y2, 'alpha':alpha, 'decryption_failure':decryption_failure }
				if(debug >= 3): # output it
					print "New best Candidate found for q: ", q, ", N: ", N, ", d: ", d
					print(json.dumps(bestParameters, indent = 4))
				aLog.write("New best Candidate found for q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + "\n")
				aLog.write(json.dumps(bestParameters, indent = 4))


	if(debug >= 1):
		print "Final Candidate found for q: ", q, ", N: ", N, ", d: ", d
		print(json.dumps(bestParameters, indent = 4))

	aLog.write("Final Candidate found for q: " + str(q) + ", N: " + str(N) + ", d: " + str(d) + "\n")
	aLog.write(json.dumps(bestParameters, indent = 4))

	# output to a file
	write_results_to_file(q, N, d, bestParameters)

	aLog.write("Wrote results to file. Finished")
	aLog.close()

	return bestParameters



# Work through the entire set of N, all calculations
# Note that we output results as we obtain them, and then
# those completed computations are skipped instead of re-computed
# so running this program repeatedly is idempotent and also will
# eventually get to a completed state even if it is interrupted a lot
def GenerateNtruParameters(q):
	######
	# TODO: make sure the output directory exists ...
	######
	if(debug >= 1):
		print "Generating NTRU Parameters for q: ", q

	if q <= 1:
		print "Error, invalid q: ", q, "!"
		sys.exit(-1)

	N = 2
	N = calculateNextGoodN(N, q)

	while ( N != 0 ):

		if(debug >= 1):
			print "Generating NTRU Parameters for N: ", N, " q: ", q

		#########################################################################
		## Note, THIS is where we will parallelize, each d should be a parallel
		## calculation and independently write the outputs
		#########################################################################
		if __name__ == "__main__":
			freeze_support()
			pool = None
			if(single_threaded == False):
				pool = multiprocessing.Pool(wmitm_parallelism)

			d = 1 # start with d=1, inside of the loop it becomes 2..
			# d=1 gives exceptional results and shouldn't be evaluated.
			# It is not really a reasonable parameter choice.

			while d < (N/3) - 1:
				d = d + 1
				if(debug >= 2):
					print "Dispatching job for q:", q, ", N:", N, ", d:", d

				# calculate and output for this N, q, d
				# This causes results to be output to files, and also detailed logs are generated.
				if(single_threaded == True):
					calculateForGivenParameters(q, N, d)
				else:
					pool.apply_async(calculateForGivenParameters, (q, N, d,))

			if(single_threaded == False):
				pool.close() # tells the pool we aren't accepting new work
				pool.join() # blocks until all workers are done

		if(debug >= 1):
			print "Moving on to the next N."

		# Gets the next N or 0
		N = calculateNextGoodN(N, q)

	print "All Jobs Completed."


def findAlpha(q, N, d, c, y2):
	sigma = mpf_nice(calculateSigma1(N, d, d, c, c, y2))
	smallest_difference = None
	smallest_difference_alpha = None

	continue_alpha = None

	# walk forward in large steps until we cross the root
	# Note: alpha_step is how many we skip ahead, this has been
	# balanced experimentally. If it's too large then we do get to
	# skip a lot in the beginning but we end up having to do a large
	# amount of 'recheck' work
	for alpha in range(0, 1000, alpha_step):
		this_alpha = float(alpha) / float(1000.0)
		# K = 2 * N - y2
		mitm = MITMRunningTime(q, y2, d, c, N, 2 * N - y2, this_alpha, sigma)
		lattice = LatticeRunningTime(N, y2, this_alpha)
		if (mitm == None) or (lattice == None):
			return None # skip this q, N, d, y2, c because it's an invalid profile
		abs_diff = (mitm - lattice)
		if abs_diff < 0:
			abs_diff = abs_diff * -1.0
		if(smallest_difference == None) or (smallest_difference > abs_diff):
			smallest_difference = abs_diff
			smallest_difference_alpha = this_alpha
#			print "this_alpha: ", this_alpha, " has become smallest with diff: ", smallest_difference
		else:
			# after the difference goes back up, we're done
			break

	# redo the same element, it's the new 'worst' one
	# and we step through each one to make sure we've found the smallest.
	continue_alpha = alpha
	# reset these
	smallest_difference = None
	smallest_difference_alpha = None

	# then walk backwards slowly until we cross it
	for alpha in range(continue_alpha, 0, -1):
		this_alpha = float(alpha) / float(1000.0)
#		print "at this_alpha: ", this_alpha
		mitm = MITMRunningTime(q, y2, d, c, N, 2 * N - y2, this_alpha, sigma)
		lattice = LatticeRunningTime(N, y2, this_alpha)
		abs_diff = (mitm - lattice)
		if abs_diff < 0:
			abs_diff = abs_diff * -1.0
	#	print "diff is: ", abs_diff
		if(smallest_difference == None) or (smallest_difference > abs_diff):
			smallest_difference = abs_diff
			smallest_difference_alpha = this_alpha
	#		print "this_alpha: ", this_alpha, " has become smallest"
		else:
			# after the difference goes back up, we're done
			break

#		print this_alpha, ", ", mitm, ", ", mitm_log2, ", ", lattice, ", ", lattice_log2, ", ", abs_diff
	return smallest_difference_alpha


# Test wrapper
def checkFindingAlpha(q, N, d, c, y2):
	print "Checking if we get the published alpha within the q, N, d, c, y2 loop"
	smallest_difference_alpha = findAlpha(q, N, d, c, y2)

	print "Found best alpha of: ", smallest_difference_alpha
	return smallest_difference_alpha


class TestOverallValidation(unittest.TestCase):

	def testReproduceResults(self):
		print "I) Running through 12 parameter sets given by original paper authors, calculating Wmitm " \
			"with these parameters to ensure that we obtain the correct values."
		print ""

		# Uncomment the following 3 blocks to do full 'qNd' searches.
		# Note that these take hours (1, 2) or days (12) to complete

#		print "1) k=112, N=401, d=113, c=27, y=693, alpha=0.095"
#		bestParameters = calculateForGivenParameters(2048, 401, 113)
#		self.assertGreaterEqual(bestParameters['securitylevel'], 112)
#		self.assertLess(bestParameters['securitylevel'], 113)
#		self.assertEqual(bestParameters['y2'], 693)
#		self.assertEqual(bestParameters['c'], 27)
#		self.assertAlmostEqual(bestParameters['alpha'], 0.095, delta=0.002)

#		print "2) k=128, N=449, d=134, c=35, y2=770, alpha=0.100"
#		bestParameters = calculateForGivenParameters(2048, 449, 134)
#		self.assertGreaterEqual(bestParameters['securitylevel'], 128)
#		self.assertLess(bestParameters['securitylevel'], 129)
#		self.assertEqual(bestParameters['y2'], 770)
#		self.assertEqual(bestParameters['c'], 35)
#		self.assertAlmostEqual(bestParameters['alpha'], 0.100, delta=0.002)

#		print "12) k=256, N=1499, d=79, c=29, y2=1984, alpha=0.174"
#		bestParameters = calculateForGivenParameters(2048, 1499, 79)
#		self.assertGreaterEqual(bestParameters['securitylevel'], 256)
#		self.assertLess(bestParameters['securitylevel'], 257)
#		self.assertEqual(bestParameters['y2'], 1984)
#		self.assertEqual(bestParameters['c'], 29)
#		self.assertAlmostEqual(bestParameters['alpha'], 0.174, delta=0.002)

		print "1) k=112, N=401, d=113, c=27, y=693, alpha=0.095"
#		self.assertAlmostEqual(checkWmitm(2048, 401, 48, 693, 113, 113, 27, 27), 112.204883082, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 401, 113, 27, 693), 0.095, delta=0.002)
		print "2) k=128, N=449, d=134, c=35, y2=770, alpha=0.100"
#		self.assertAlmostEqual(checkWmitm(2048, 449, 57, 770, 134, 134, 35, 35), 128.068860169, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 449, 134, 35, 770), 0.100, delta=0.002)
		print "3) k=192, N=677, d=157, c=45, y2=1129, alpha=0.096"
#		self.assertAlmostEqual(checkWmitm(2048, 677, 129, 1129, 157, 157, 45, 45), 192.387537617, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 677, 157, 45, 1129), 0.096, delta=0.002)
		print "4) k=256, N=1087, d=120, c=39, y2=1630, alpha=0.127"
#		self.assertAlmostEqual(checkWmitm(2048, 1087, 386, 1630, 120, 120, 39, 39), 256.690210632, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 1087, 120, 39, 1630), 0.127, delta=0.002)
		print "5) k=112, N=541, d=49, c=15, y2=800, alpha=0.149"
#		self.assertAlmostEqual(checkWmitm(2048, 541, 192, 800, 49, 49, 15, 15), 112.271625998, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 541, 49, 15, 800), 0.149, delta=0.002)
		print "6) k=128, N=613, d=55, c=17, y2=905, alpha=0.142"
#		self.assertAlmostEqual(checkWmitm(2048, 613, 225, 905, 55, 55, 17, 17), 128.281639552, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 613, 55, 17, 905), 0.142, delta=0.002)
		print "7) k=192, N=887, d=81, c=27, y2=1294, alpha=0.143"
#		self.assertAlmostEqual(checkWmitm(2048, 887, 345, 1294, 81, 81, 27, 27), 192.166342224, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 887, 81, 27, 1294), 0.143, delta=0.002)
		print "8) k=256, N=1171, d=106, c=37, y2=1693, alpha=0.144"
#		self.assertAlmostEqual(checkWmitm(2048, 1171, 474, 1693, 106, 106, 37, 37), 256.294181758, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 1171, 106, 37, 1693), 0.144, delta=0.002)
		print "9) k=112, N=659, d=38, c=13, y2=902, alpha=0.175"
#		self.assertAlmostEqual(checkWmitm(2048, 659, 313, 902, 38, 38, 13, 13), 112.12377268, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 659, 38, 13, 902), 0.175, delta=0.002)
		print "10) k=128, N=761, d=42, c=15, y2=1026, alpha=0.183"
#		self.assertAlmostEqual(checkWmitm(2048, 761, 378, 1026, 42, 42, 15, 15), 128.190480858, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 761, 42, 15, 1026), 0.183, delta=0.002)
		print "11) k=192, N=1087, d=63, c=23, y2=1464, alpha=0.175"
#		self.assertAlmostEqual(checkWmitm(2048, 1087, 550, 1464, 63, 63, 23, 23), 192.385016473, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 1087, 63, 23, 1464), 0.175, delta=0.002)
		print "12) k=256, N=1499, d=79, c=29, y2=1984, alpha=0.174"
#		self.assertAlmostEqual(checkWmitm(2048, 1499, 809, 1984, 79, 79, 29, 29), 256.161994836, places=2)
		self.assertAlmostEqual(checkFindingAlpha(2048, 1499, 79, 29, 1984), 0.174, delta=0.002)



	# highest level tests
	def testNtruParameterSelection(self):
		print "testNtruParameterSelection... TODO NO TEST HERE YET"
#		self.assertEqual(NtruParameterSelection("size=Nxlog2q", 112, constant_q), {'N':401, 'd':113, 'q':constant_q})



########################################
# Main routine
#
def main(argv):

	global k
	global debug
	global constant_q

	unit_test = False
	validation_test = False

	# parse options
	try:
		opts, args = getopt.getopt(argv, "u:k:q:d:v", ["unit-test","k-bit-security=","q-large-mod=","debug=", "validation-test"])
	except getopt.GetoptError:
		print "Options not understood"
		sys.exit()

	for i in range (0, len(opts)):
		opt, arg = opts[i]
		if opt in ("-u", "--unit-test"):
			unit_test = True
		elif opt in ("-v", "--validation-test"):
			validation_test = True
		elif opt in ("-k", "--k-bit-security"):
			k = int(arg)
		elif opt in ("-q", "--q-large-mod"):
			constant_q = int(arg)
		elif opt in ("-d", "--debug"):
			debug = int(arg)

	sys.argv = [ sys.argv[0] ]

	if(unit_test == True):
		print "Running unit tests..."
		unittest.main()
	if(validation_test == True):
		print "Running validation testing to ensure that your environment is calculating as expected"

		suite = unittest.TestLoader().loadTestsFromTestCase(TestOverallValidation)
		unittest.TextTestRunner(verbosity=2).run(suite)

	else:
		# Do all the work.
		GenerateNtruParameters(constant_q)


if __name__ == '__main__':
	freeze_support()
	# main(sys.argv[1:])
	q = 2048
	y2 = 402
	d = 134
	c = 134
	N = 401
	alpha = 0
	sigma = mpf_nice(calculateSigma1(N, d, d, c, c, y2))
	mitm = MITMRunningTime(q, y2, d, c, N, 2 * N - y2, alpha, sigma)





############################################
## File helpers
############################################

# True or False..
def results_file_exists(q, N, d):
	return os.path.isfile(results_filename(q, N, d))

# Returns hash of results or None
def results_from_file(q, N, d):
	try:
		json_data_file = open(results_filename(q, N, d))
		return_data = json.load(json_data_file)
		json_data_file.close()
		return return_data

	except ValueError:
		print "Unable to read or decode JSON in results file " + str(results_filename(q, N, d))
		return None

def write_results_to_file(q, N, d, results):
	# NOTE, i'm probably missing some exception handling here..
	with open(results_filename(q, N, d), 'w') as outfile:
		json.dump(results, outfile)
	outfile.close()

def results_filename(q, N, d):
	return (str(results_output_file_prefix) + "results_" + str(q) + "_" + str(N) + "_" + str(d) + ".json")


def open_unique_detailed_logfile(q, N, d):
	# Unique filename based on timestamp..
	fd = open((str(detailed_log_file_prefix) + strftime("%Y%m%d_%H%M%S_", gmtime()) + str(q) + "_" + str(N) + "_" + str(q) + ".out"), "w")
	return fd



###################################
# Unused functions, calculating the 'current' best known attacks
# We don't calculate that, we calculate more conservatively,
# meaning that we assume there are more effective attacks than we know today.


# this is the 'current' version, should not be using this for our
# calculations because we are calculating 'conservative' estimates
def MITMTrialsBeforeCollision_Current(q, N, c, y2, alpha, sigma):
# 'N' (trials before collision) is:
# [ 2N-y2   ] ([  [    c    ] [ Ps ] )^{-1/2}
# [ c/2,c/2 ] ([  [ c/2, c/2] [    ] )
# multinomial coefficient

	left_term = nchoosers(2*N - y2, c // 2, c // 2)
	right_inner_term = nchoosers(c, c // 2, c // 2)
	ps_term = calculatePs(q, N, y2, alpha, sigma)
	right_term = (right_inner_term * ps_term) ** (-0.5)
	# TODO; check for None
	result = left_term * right_term
	return result


########################################
# prime helpers to avoid installing extra deps like gmpy2..
# performance of this software doesn't matter as much as ease of use
def isPrime(candidate):

	if(candidate <= 1):
		return False

	if(candidate <= 12):
		if(candidate == 2 or candidate == 3 or candidate == 5 or candidate == 7 or candidate == 11):
			return True
		else:
			return False

	if((pow(mpf_nice(2), candidate - 1))%candidate != 1):
		return False

	for x in range(7, candidate, 6):
		if((candidate % x) == 0):
			return False

	for x in range(11, candidate, 6):
		if((candidate % x) == 0):
			return False

	return True


# Note that prime numbers > 6 have to be equiv +/-1 mod 6, to limit
# required number of trials
def nextIncrement(startingpoint):

	rem = (startingpoint % 6)

	if(rem == 1):
		return 4
	if(rem == 5 or rem == 3):
		return 2
	if(rem == 2):
		return 3
	return 1


# return the next prime. Input doesn't need to be prime.
def nextPrime(startingpoint):

	if(startingpoint < 2):
		return 2
	if(startingpoint == 2):
		return 3
	if(startingpoint <= 4):
		return 5

	startingpoint = startingpoint + nextIncrement(startingpoint)

	while(isPrime(startingpoint) == False):
		startingpoint = startingpoint + nextIncrement(startingpoint)

	return startingpoint


# next prime satisfying the NTRU requirements
def calculateNextGoodN(N, q):
	# to calculate the order of 2 mod N and see if it's N-1 or (N-1)/2
	# go from 2^1 upward, (2^2, 2^3) and determine when the result is = 1
	# if the exponent = N-1 or (N-1)/2 then we can proceed
	# check for order of 2 mod N. Iterate the primes until we get one
	# if N > q then return 0
	global LOG

	while 1:
		N = nextPrime(N)

		if(N > q):
			return 0

		answer = 0

		for i in range(1, N):
			if (pow(2, i) % N) == 1:
				answer = i
				break

		# TODO are we being too flexible with N by including /2?
		if(answer == N-1) or (answer == (N-1)/2):
			return N

# verify if we have a good N or not
# Not used for our main loop but used in verification of results..
def isGoodN(N, q):
	if(N <= 0):
		return False

	if(N > q):
		return False

	if(isPrime(N) == False):
		return False

	answer = 0

	for i in range(1, N):
		if (pow(2, i) % N) == 1:
			answer = i
			break

	if(answer == N-1) or (answer == (N-1)/2):
		return True

	return False


def outputCacheStatistics():
	if (caching_disabled):
		print "Caching Disabled. No statistics to output."
		return

	print "Cache Statistics"
	print "================"
	print "binomial hits: ", binomial_cache_hits, ", misses: ", binomial_cache_misses, \
	", hit%%: ~", float(binomial_cache_misses) / float(binomial_cache_misses + binomial_cache_hits + 1), \
	", cache entries: ", len(binomial_cache_dict)

	print "multinomial hits: ", multinomial_cache_hits, ", misses: ", multinomial_cache_misses, \
	", hit%%: ~", float(multinomial_cache_misses) / float(multinomial_cache_misses + multinomial_cache_hits + 1), \
	", cache entries: ", len(multinomial_cache_dict)


# binomial coefficient with caching implemented
def nchooser(n, r):
	key = str(n) + ":" + str(r)
	global binomial_cache_hits
	global binomial_cache_misses
	global binomial_cache_dict
	if key in binomial_cache_dict:
		# cache hit
		binomial_cache_hits += 1
		if (caching_enabled == True):
			return binomial_cache_dict[key]

	if n < r:
		return None
	if r < 0:
		return None

	r = min(r, n - r)

	if r == 0:
		return 0
	if n == 0:
		return 0

	res = mpf_nice(0)

	for i in range(0, r):
		res += (LOG[n - i] - LOG[r - i])

	result = mpf_nice(res)

	binomial_cache_dict[key] = result
	binomial_cache_misses += 1

	return result

def Nchooser(n, r):
	nchoose_result = nchooser(n, r)
	if(nchoose_result == None):
		return None

	return mpf_nice(pow( mpf_nice(2), nchoose_result))


# multinomial coefficient
def nchoosers(n, r, s):
	key = str(n) + ":" + str(r) + ":" + str(s)
	global multinomial_cache_hits
	global multinomial_cache_misses
	global multinomial_cache_dict
	if key in multinomial_cache_dict:
# (  n  )
# ( r, s)
		# cache hit
		multinomial_cache_hits += 1
		if (caching_enabled == True):
			return multinomial_cache_dict[key]

	if n < r + s:
		return None
	if r < 0:
		return None
	if s < 0:
		return None
	t = n - r - s

	if(t < r):
		a = t
		t = r
		r = a
	if(t < s):
		a = t
		t = s
		s = a
	elif(s < r):
		a = s
		s = r
		r = a

	if n == 0:
		return 0
	if r == 0:
		return nchooser(n, s)
	if s == 0 or t == 0:
		return nchooser(n, r)

	res = mpf_nice(0)

	for i in range(1, r + 1):
		res += mpf_nice(LOG[t + i] - LOG[i])
	for i in range(1, s + 1):
		res += mpf_nice(LOG[t + r + i] - LOG[i])

	result = mpf_nice(pow(mpf_nice(2), res))

	multinomial_cache_dict[key] = result
	multinomial_cache_misses += 1

	return result


def Nchoosers(n, r, s):
	return mpf_nice(pow(mpf_nice(2), nchoosers(n, r, s)))

