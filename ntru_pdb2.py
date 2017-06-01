#!/usr/bin/python

# implement the improved NTRU parameter finding algorithm
# from the paper "Choosing Parameters for NTRUEncrypt" by
# Jeff Hoffstein, Jill Pipher, John M. Schanck, Joseph H.
# Silverman, William Whyte, and Zhenfei Zhang

import unittest
import math
import sys
import getopt
import gmpy2
import sympy
import scipy

from gmpy2 import mpz
from mpmath import *
from scipy.misc import comb
from sympy.ntheory.generate import nextprime
from scipy.special import gammaln
from numpy import *

mp.prec = 128
res = mpf(1)/mpf(pow(2, 107))

# security level k and debug level debug
security_levels = [112, 128, 192, 256]
k = 112
debug = 0

# ---------------- debug functions ------------------#

def DEB_out(*args):

    global debug
    if(len(args) < 2): return
    deb_level = int(args[0])
    if(debug >= deb_level):
        for i in range(1, deb_level):
            sys.stdout.write('\t')
        for i in range(1,len(args)):
            sys.stdout.write(str(args[i]))
            sys.stdout.write(' ')
        sys.stdout.write('\n')
        sys.stdout.flush()

    return


# ---------------- helper functions -----------------#

def array_to_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]

def binomial(n,k):

    if(n < 0): return 0
    if(k < 0): return 0
    if(n < k): return 0

    return mpz(scipy.misc.comb(n, k, exact = 1))

# BKZ 2.0 simulator
def bkz(N, beta, target, abort = 50):
    if beta < 50:
        return None
    r = [0.4809337322749968, 0.4889068929146757, 0.4629910732303647, 0.4384921120061095, 0.4198271756529734,
    0.3940124751357192, 0.3793579556691379, 0.3552017168415738, 0.3375032857978846, 0.3229996676156046,
    0.3103169826524305, 0.2978627511364960, 0.2828082600293407, 0.2685092222965025, 0.2470246073218571,
    0.2345601366183950, 0.2236298423327614, 0.2026125221670087, 0.1833511717333619, 0.1635239915325074,
    0.1460909610754462, 0.1239402813211751, 0.1033442833745716, 0.08072183355489210, 0.05747352858422083,
    0.03615285314640355, 0.009734731674006085, -0.01314276679308946, -0.03859536413875225, -0.06166664730992491,
    -0.08732858253410711, -0.1159733213895935, -0.1395873057069733, -0.1685959449423031, -0.2009987452911466,
    -0.2272943479144534, -0.2548892487960738, -0.2845907037340676, -0.3130406180111631, -0.3439519155213564,
    -0.3729166620199606, -0.4037203497626708, -0.4279121623225402, -0.4591242077605871, -0.4851668230787535,
    -0.5069333755274962, -0.5312523582495852, -0.5480002333962808, -0.5470408985906416, -0.5201614648988958]

    c = zeros(beta)
    for d in range(51,size(c)+1):
        c[d-1] = (mpf(1.0)/mpf(d))*mpf(gammaln(d/2.0 + 1)) - mpf(1)/mpf(2)*mpf(log(pi))

    ll = zeros(N)
    for i in range(0, N): ll[i] = mpf((N-2*i))*mpf(log(1.01263))
    vs = mpf(sum(ll))/mpf(N)
    for i in range(0, N): ll[i] = mpf(ll[i]) - mpf(vs)
    llp = zeros(N)

    R = 0
    while(exp(ll[0]/N) > target and R < abort):
        phi = 1
        for k in range(1, N - 50 + 1):
            d = min(beta, N - k + 1)
            f = min(k + beta, N)
            logV = mpf(sum(ll[0:f])) - mpf(sum(llp[0:(k-1)]))
            if phi:
                if mpf(logV)/mpf(d) + mpf(c[d-1]) < mpf(ll[k-1]):
                    llp[k-1] = mpf(logV)/mpf(d) + mpf(c[d-1])
                    phi = 0
            else:
                llp[k-1] = mpf(logV)/mpf(d) + mpf(c[d-1])
        logV = mpf(sum(ll)) - mpf(sum(llp[0:(N-50)]))
        for k in range(0, 50):
            llp[N-50+k] = mpf(logV)/mpf(50.0) + mpf(r[k])
        ll = copy(llp)
        R += 1
        if phi:
            return "failure", exp(ll[0]/N), R

    if R >= abort:
        return "failure", exp(ll[0]/N), R
    return "success", exp(ll[0]/N), R

# ---------------- NTRU parameter functions ---------#

# increase N up to the next prime such that 2 has order N-1 mod N
def NTRU_next_Prime_Table4(N):

    # Algorithm allows only for primes > 100
    if(N < 100):
        N = 100

    # Search for next prime until we find one where 2 has maximal order
    while(True):
        N = sympy.ntheory.generate.nextprime(N)

	x = (mpz(2)**mpz((N - 1)/2)) % N;
	print x, " "
	if(x == N-1): return N;
	if(x != 1): continue;

        n = (N - 1)/2
        p = 3
        ok = True
        while(n > 1 and ok):
            if( (mpz(2)**mpz((N - 1)/p) - 1) % N == 0 ):
                ok = False
                break

            while(n % p == 0): n = n // p
            if(n == 1): break
            while(n % p != 0): p = sympy.ntheory.generate.nextprime(p)

        if(ok == True): break

    return N

# increase N up to the next prime such that 2 has order (N-1)/2 mod N
def NTRU_next_Prime_Table5(N):

    # Algorithm allows only for primes > 100
    if(N < 100):
        N = 100

    # Search for next prime until we find one where 2 has maximal order
    while(True):

        N = sympy.ntheory.generate.nextprime(N)

        n = (N - 1)/2
	if(n%2 == 0):
            if( (mpz(2)**mpz((N - 1)/4) - 1) % N == 0 ):
            	continue
	    while(n%2 == 0): n = n//2
		
        p = 3
        ok = True
        while(n > 1 and ok):
            if( (mpz(2)**mpz((N - 1)/p) - 1) % N == 0 ):
                ok = False
                break

            while(n % p == 0): n = n // p
            if(n == 1): break
            while(n % p != 0): p = sympy.ntheory.generate.nextprime(p)

        if(ok == True): break

    return N

def NTRU_next_Prime(N):
    if(N < 100):
        N = 100

    while(True):
        N = sympy.ntheory.generate.nextprime(N)
        n = (N - 1)/2
        p = 2
        ok = True
        isTable5 = False
        if( (mpz(2)**mpz((N - 1)/2) - 1) % N == 0 ):
            isTable5 = True
            p = 3

        if (isTable5):
            if( (mpz(2)**mpz((N - 1)/4) - 1) % N == 0 ):
                continue

            while(n % 2 == 0): n = n // 2

        while(n > 1 and ok):
            if( (mpz(2)**mpz((N - 1)/p) - 1) % N == 0 ):
                ok = False
                break

            while(n % p == 0): n = n // p
            if(n == 1): break
            while(n % p != 0): p = sympy.ntheory.generate.nextprime(p)

        if(ok == True): break

    return N

# ---------------- functions defined in the paper -----------------#

# size of the space of ternary polynomials P_N(d1, d2, d3)
# see definition on p.5
def NTRU_P_N(N, d1, d2, d3):

    t1 = mpz(binomial(N, d1))
    t2 = mpz(binomial(N - d1, d1))
    t3 = mpz(binomial(N, d2))
    t4 = mpz(binomial(N - d2, d2))
    t5 = mpz(binomial(N, d3))
    t6 = mpz(binomial(N - d3, d3))

    return mpz(t1*t2*t3*t4*t5*t6)


# The cost of a non-hybrid MITM search on f
def NTRU_product_form_search(N, d1, d2, d3):

    ret = math.sqrt(mpf(NTRU_P_N(N, d1, d2, d3))/mpf(N))
    DEB_out(4, "NTRU_product_form_search =", ret)
    return ret


#####################################################################

# experimenal delta* coefficient
# see definition on p.10
def NTRU_delta_star(Lambda):

    if(Lambda <= 60):
        return 1.009
    if(Lambda <= 80):
        return 1.008
    if(Lambda <= 128):
        return 1.007
    if(Lambda <= 256):
        return 1.005
    return 1.0


#####################################################################

# expression p(v_{a,b})
# see p.9 - needed for eq (4)
# for a polynomial with d 1's and (d+1) -1's, what is the probability
# that the last K coordiantes contain a +1's and b -1's
def NTRU_p(N, K, d, a, b):

    t1 = mpz(binomial(N - K, d - a))
    #t2 = mpz(binomial(N - K - d + a, d - b))	# paper uses e = d, but should be e = d + 1
    t2 = mpz(binomial(N - K - d + a, d + 1 - b))
    t3 = mpz(binomial(N, d))
    #t4 = mpz(binomial(N - d, d))			# paper uses e = d, but should be e = d + 1
    t4 = mpz(binomial(N - d, d + 1))

    return ( mpf(t1 * t2) / mpf(t3 * t4) )


# number of choices of distinct v_{a,b}
# how many was are there of selecting a +1's and b -1's from K coordinates
def NTRU_S_ab(K, a, b):

    t1 = mpz(binomial(K, a))
    t2 = mpz(binomial(K - a, b))

    return ( mpz(t1 * t2) )


# Shannon entropy -
# This is the expression just below eq (4) on p.9
# contribution to the total entropy from the (a,b) component
def NTRU_shannon_p_ab(N, K, d, a, b):

    p = NTRU_p(N, K, d, a, b)
    n = NTRU_S_ab(K, a, b)
    u = n * p
    # ensure we don't crash for small p - eq below ensures we
    # discard small contributions that don't matter anyway
    if (u < res): return mpf(0.0)
    return  u * math.log(mpf(p), 2)


# Difficulty of the MITM search for KEY with the hybrid improvement
# eq (4) p.9
def NTRU_shannon_p(N, K, d):

    ret = mpf(0)

    for a in range(0, d + 1):
        for b in range(0,  d + 2):				# really should be d + 2
            ret -= NTRU_shannon_p_ab(N, K, d, a, b)

    ret = ( mpf(0.5)*mpf(ret - math.log(N, 2)) )

    DEB_out(4, "NTRU_shannon_p =", ret)
    return ret



# the q-term in eq (7) on p.9. Note the definition of admissible
# index pairs (i,j) in I(dm) together with the condition i - j = y
# imply that dm <= j = i - y => i >= dm + y
def q(e1, e2, y, N, dm):

    low = dm + y				# the lowest possible value for i so that i >= dm, and j = (i - y) >= dm
    high1 = (N - dm + y - 1) // 2 + 1	# i < high1 ensures j = (i - y) < N - dm - i <=> 2i < N - dm + y
    high2 = N - 2*dm			# i < high2 is required by definition of I(dm)
    high = min(high1, high2)		# final upper bound
    if(low >= high): return mpf(0)
    t = mpf(0)
    for i in range(low, high):
        t += binomial(N, i) * binomial(N - i, i - y)

    return -mpf(math.log( NTRU_S_ab(N, e1, e2), 2 ) - math.log( t ))


# to be used for message attacks
# probability that a message with e1 +1's and e2 -1's projects onto
# a +1's and b -1's in the last K coordinates
def NTRU_p_msg(N, K, e1, e2, a, b):

    t1 = mpz(binomial(N - K, e1 - a))
    t2 = mpz(binomial(N - K - e1 + a, e2 - b))
    t3 = mpz(binomial(N, e1))
    t4 = mpz(binomial(N - e1, e2))

    return ( mpf(t1 * t2) / mpf(t3 * t4) )


# Size of the message space M of N-vectors satisfying the dm constraints,
# i.e., there are >= dm +1's and >= dm -1's
def NTRU_S_msg(N, dm):
    # size of the message space satisfying the dm constraint (>= dm +1's and -1's)
    M = mpz(0)
    for e1 in range(dm, N - dm + 1):
        for e2 in range(dm, N - dm - e1 + 1):
            M += NTRU_S_ab(N, e1, e2)

    return M


# contribution of the (a,b) component to the total entropy
def NTRU_shannon_p_ab_msg(N, K, e1, e2, a, b):

    p = NTRU_p_msg(N, K, e1, e2, a, b)
    n = NTRU_S_ab(K, a, b)
    u = n * p
    # ensure we don't crash for small p - eq below ensures we
    # discard small contributions that don't matter anyway
    if (u < res): return mpf(0.0)
    return  u * math.log(mpf(p), 2)


# Difficulty of the MITM search for MESSAGE with the hybrid improvement
def NTRU_shannon_p_msg(N, K, dm, e1, e2):

    # size of message space M satisfying dm constraint
    M = NTRU_S_msg(N, dm)

    # number of elements from T_N(e1, e2) in M
    m = NTRU_S_ab(N, e1, e2)

    ret = mpf(0)
    for a in range(0, e1 + 1):
        for b in range(0, e2 + 1):
            ret -= NTRU_shannon_p_ab_msg(N, K, e1, e2, a, b)

    DEB_out(4, "NTRU_shannon_p_msg", ret)
    return mpf(0.5) * mpf(ret)


def NTRU_shannon_p_msg_min(N, K, dm, k):

    e1 = int(round(float(N)/float(3)))
    e2 = e1

    min1 = NTRU_shannon_p_msg(N, K, dm, e1, e2) - q(e1, e2, e1 - e2, N, dm)
    if(min1 < k):
        return min1

    diff = N - 3*dm
    for e1 in range (dm, N - diff):
        e2 = e1 + diff
        min2 = NTRU_shannon_p_msg(N, K, dm, e1, e2) - q(e1, e2, -diff, N, dm)
        if(min2 < min1):
            min1 = min2
            DEB_out(4, "NTRU_shannon_p_msg_min", min1)
            if(min1 < k):
                return min1

    for e1 in range (dm + diff, N):
        e2 = e1 - diff
        min2 = NTRU_shannon_p_msg(N, K, dm, e1, e2) - q(e1, e2, diff, N, dm)
        if(min2 < min1):
            min1 = min2
            DEB_out(4, "NTRU_shannon_p_msg_min", min1)
            if(min1 < k):
                return min1

    DEB_out(3, "NTRU_shannon_p_msg_min", min1)
    return min1

#####################################################################

# compute log_2(eta) as defined by the equation on p.10
# just above eq (6)
# ensures we can recover the closest vector in T
def NTRU_log2_eta(N, log2_q, r1, K):

    denom = mpf(2*N - (K + r1))	# a.k.a. y2 - y1
    eta = mpf( 1 + mpf((N - r1)*log2_q)/(denom*denom) - mpf(1.0)/denom )
    DEB_out(4, "eta:", round(eta, 4))
    return eta

# compute eta as defined by the equation on p.10
# just above eq (6)
def NTRU_eta(N, log2_q, r1, K):
    denom = mpf(2*N - (K + r1))	# a.k.a. y2 - y1
    eta = mpf(mpf((N - r1)*log2_q)/(denom*denom) - mpf(1.0)/denom)
    return mpz(2)**mpf(eta)

#####################################################################


# on top of p.9 there is a condition dm must satisfy
def NTRU_dm_valid(N, dm):
    z = mpz(0)
    for i in range(dm, N - 2*dm + 1):
        for j in range(dm, N - dm - i + 1):
            z += binomial(N, i)*binomial(N - i, j)

    if(z >= (mpf(3)**N) * mpf(1023) / mpf(1024)):
        return True
    else:
        return False

# Paper says in step 8 that dm must be the biggest value satisfying eq (5),
# however, eq (5) is not abiyt dm at all. I think they mean the biggest value
# satisfying the top equation on p11
def NTRU_dm(N):
    r1 = 0
    r2 = N
    while(r1 < r2 - 1):
        r3 = (r1 + r2) // 2
        if(NTRU_dm_valid(N, r3)):
            r1 = r3
        else:
            r2 = r3
    DEB_out(4, "dm =", r1, "(N =", N, ")")
    return r1


#####################################################################


# cost of hybrid search
def NTRU_wmitm(N, log2_q, dg, dm, k):

    r1 = k
    for r2 in range(N + 1, 2*N + 1):

        K = 2*N - r2

        # We consider only K's for which the cost of lattice reduction is > 2^k (r1 == k)
        eta = NTRU_eta(N, log2_q, r1, K)
        if(eta > NTRU_delta_star(k) or eta <= 1):
            DEB_out(4, "NTRU_log2_eta(N =", N, ", log2(q) =", log2_q, ", y1 =", r1,
                   ", y2 =", r2, ", K =", K, ") =", round(eta, 4))
            continue

        # Given any such a K, we need to ensure the MITM attack costs at least 2^k (r1 == k)
        k2 = NTRU_shannon_p(N, K, dg)

        if(k2 < k):
            # k2 gets smaller as K decreases, so no point continuing
            DEB_out(3, "NTRU_shannon_p(", N, ",", K, ",", dg, ") =", round(k2,4))
            break

        return k2, K

        # doesn't work
        k2m = NTRU_shannon_p_msg_min(N, K, dm, k)
        if(k2m < k):
            # k2 gets smaller as K decreases, so no point continuing
            DEB_out(3, "NTRU_shannon_p(", N, ",", K, ",", dg, ") =", round(k2,4))
            DEB_out(3, "NTRU_shannon_p_msg(", N, ",", K, ",", dm, ") =", round(k2m,4))
            break

        DEB_out(2, "NTRU_log2_eta(N =", N, ", log2(q) =", log2_q, ", y1 =", r1, ", y2 =", r2,
               ", K =", K, ") =", round(eta, 4))
        DEB_out(2, "NTRU_shannon_p(", N, ",", K, ",", dg, ") =", round(k2,4))
        DEB_out(2, "NTRU_shannon_p(", N, ",", K, ",", dm, ") =", round(k2m,4))

        return k2, K

    return 0.0, 0


#####################################################################


# erfc term - note p = 3 is implicit (factor 6 in denominator is 3p)
def NTRU_erfc(N, q, sigma):

    return mpf(N * math.erfc((q - 2)/(6*sqrt(2)*sigma)))


#####################################################################


# ---------------- main algorithm -----------------#

def NTRU_parameters(k):

    # initialise local variables
    N = 0
    d1 = 0
    d2 = 0
    d3 = 0
    dg = 0
    dm = 0
    q = 0
    k1 = 0
    k2 = 0
    pdf = 0
    K = 0

    while(True):

        # The inner loop finds a (N, d1, d2, d3, dg, dm) tuple that satisfies the search-hardness
        while(True):

            #------------ Steps 1-3 -----------#
            N = NTRU_next_Prime(N)

            #------------ Step 4 --------------#
            dg = int(math.floor(float(N)/float(3) + 0.5))

            #------------ Step 5 --------------#
            d1 = int(math.ceil(0.25*(math.sqrt(1 + float(8*N)/float(3)) - 1)))

            #------------ Step 6 --------------#
            d2 = int(math.ceil((float(N)/float(3) - d1)/float(2*d1)))

            #------------ Step 7 --------------#
            d3 = int(max(math.ceil(float(d1)*0.5 + 1), math.ceil(float(N)/float(3) - 2*d1*d2)))

            #------------ Step 8 --------------#
            dm = NTRU_dm(N)

            #------------ Step 9 --------------#
            k1 = mpz(0.5*math.log(mpf(NTRU_P_N(N, d1, d2, d3))/mpf(N), mpf(2)))

            DEB_out(1, "Step  9:", "N =", N, "d1 =", d1, "d2 =", d2, "d3 =",\
                d3, "dg =", dg, "dm =", dm, "k1 =", k1)

            #------------ Steps 10-13 ---------#
            # Cost of combinatorial search is k1
            if(k1 >= k):
                break

        DEB_out(1, "Step 13:", "Combo search ok ( k1 =", round(k1, 4), ")")

        #------------ Step 14 --------------#
        sigma = math.sqrt( (4*d1*d2 + 2*d3) * (float(N - dm + 2*dg + 1)/float(N)) )

        #------------ Step 15 --------------#
        # This step makes sure the risk of decryption failure is < 2^{-k}
        # by making q large enough - q not used so far
        q = mpz(2)
        log2_q = 1
        twopowk1 = mpf(pow(2, k1))
        while(True):

            pdf = NTRU_erfc(N, q, sigma)
            if( (pdf * twopowk1) < 1 ):
                break
            q = mpz(2*q)
            log2_q += 1

        DEB_out(4, "Working with sigma =", round(sigma, 4))
        DEB_out(1, "Step 15:", "Decryption failure ok for q =", q, "( log2(p) =", round(log(pdf, 2), 4), ")" )

        #------------ Step 16 --------------#
        # Cost of hybrid search is k2
        k2, K = NTRU_wmitm(N, log2_q, dg, dm, k)
        if(k2 > 0):
            DEB_out(1, "Step 16:", "NTRU_wmitm =", round(k2, 4))
        else:
            DEB_out(1, "Step 16:", "NTRU_wmitm failed")

        #------------ Step 17-20 -----------#
        cont = True
        while(cont):

            if(k2 < k):
                # The parameter is not strong enough for the WMITM attack
                # breaking here while cont == True means we will continue the outer while(True) loop
                break

            DEB_out(1, "Step 17:", "Hybrid search ok ( k2 =", round(k2, 4), ")")

            #------------ Step 21 --------------#
            # Check if we might be able to use half of q
            qprime = q // 2

            #------------ Step 22 --------------#
            pdfb = NTRU_erfc(N, qprime, sigma)
            if( pdfb * mpf(pow(2, k)) < 1 ):

                #------------ Step 23 --------------#
                q = qprime
                log2_q -= 1
                DEB_out(1, "Step 23:", "Decryption failure still ok ( p =", round(log(pdf, 2), 4), ")" )

                #------------ Step 24 --------------#
                k2b, Kb = NTRU_wmitm(N, log2_q, dg, dm, k)

                #------------ Step 25 --------------#
                if(k2b >= k):

                    # ok now - we use this one and quit
                    # just cache the new parameters
                    if(k2b > k2):
                        k2 = k2b
                        K = Kb
                        pdf = pdfb

                    cont = False

                #------------ Step 26 --------------#
                DEB_out(1, "Step 24:", "Hybrid search security ( k2 =", round(k2, 4), ")")

            else:

                # just keep the original q
                cont = False

        # if cont == True, k2 was too small - go back to the beginning (continue the while(True) loop)
        if(not cont):
            break

    print "N =", N, "q =", q, "dg =", dg, "d1 =", d1, "d2 =", d2, "d3 =", d3, "dm =", dm
    print "K =", K, "Cost (k2):",  round(k2, 4)
    print "Product form search cost (k1):", k1
    print "Decryption fail prob. (p_dec):", round(math.log(pdf, 2), 2)


# ---------------- unit tests ---------------------#

class TestNTRU(unittest.TestCase):

    def test_binomial(self):
        self.assertEqual(binomial(4,2), 6)
        self.assertEqual(binomial(16,8), 12870)
        self.assertEqual(binomial(123,65), 626530645220345271243285139240656129)
        self.assertEqual(binomial(12,0), 1)
        self.assertEqual(binomial(4,12), 0)
        self.assertEqual(binomial(23,15), 490314)
        self.assertEqual(binomial(233,155), 1786364572172824703996858613353160486834994843862378906323850064)

    def test_NTRU_next_Prime_Table4(self):
        self.assertEqual(NTRU_next_Prime_Table4(100), 101)
        self.assertEqual(NTRU_next_Prime_Table4(101), 107)
        self.assertEqual(NTRU_next_Prime_Table4(107), 131)
        self.assertEqual(NTRU_next_Prime_Table4(131), 139)
        self.assertEqual(NTRU_next_Prime_Table4(139), 149)
        self.assertEqual(NTRU_next_Prime_Table4(149), 163)
        self.assertEqual(NTRU_next_Prime_Table4(163), 173)
        self.assertEqual(NTRU_next_Prime_Table4(173), 179)
        self.assertEqual(NTRU_next_Prime_Table4(179), 181)
        self.assertEqual(NTRU_next_Prime_Table4(181), 197)
        self.assertEqual(NTRU_next_Prime_Table4(197), 211)
        self.assertEqual(NTRU_next_Prime_Table4(211), 227)
        self.assertEqual(NTRU_next_Prime_Table4(227), 269)
        self.assertEqual(NTRU_next_Prime_Table4(269), 293)
        self.assertEqual(NTRU_next_Prime_Table4(293), 317)
        self.assertEqual(NTRU_next_Prime_Table4(317), 347)
        self.assertEqual(NTRU_next_Prime_Table4(347), 349)
        self.assertEqual(NTRU_next_Prime_Table4(349), 373)
        self.assertEqual(NTRU_next_Prime_Table4(373), 379)
        self.assertEqual(NTRU_next_Prime_Table4(379), 389)
        self.assertEqual(NTRU_next_Prime_Table4(389), 419)
        self.assertEqual(NTRU_next_Prime_Table4(419), 421)
        self.assertEqual(NTRU_next_Prime_Table4(421), 443)
        self.assertEqual(NTRU_next_Prime_Table4(443), 461)
        self.assertEqual(NTRU_next_Prime_Table4(461), 467)
        self.assertEqual(NTRU_next_Prime_Table4(467), 491)
        self.assertEqual(NTRU_next_Prime_Table4(491), 509)
        self.assertEqual(NTRU_next_Prime_Table4(509), 523)
        self.assertEqual(NTRU_next_Prime_Table4(523), 541)
        self.assertEqual(NTRU_next_Prime_Table4(541), 547)
        self.assertEqual(NTRU_next_Prime_Table4(547), 557)
        self.assertEqual(NTRU_next_Prime_Table4(557), 563)
        self.assertEqual(NTRU_next_Prime_Table4(563), 587)
        self.assertEqual(NTRU_next_Prime_Table4(587), 613)
        self.assertEqual(NTRU_next_Prime_Table4(613), 619)
        self.assertEqual(NTRU_next_Prime_Table4(619), 653)
        self.assertEqual(NTRU_next_Prime_Table4(653), 659)
        self.assertEqual(NTRU_next_Prime_Table4(659), 661)
        self.assertEqual(NTRU_next_Prime_Table4(661), 677)
        self.assertEqual(NTRU_next_Prime_Table4(677), 701)
        self.assertEqual(NTRU_next_Prime_Table4(701), 709)
        self.assertEqual(NTRU_next_Prime_Table4(709), 757)
        self.assertEqual(NTRU_next_Prime_Table4(757), 773)
        self.assertEqual(NTRU_next_Prime_Table4(773), 787)
        self.assertEqual(NTRU_next_Prime_Table4(787), 797)
        self.assertEqual(NTRU_next_Prime_Table4(797), 821)
        self.assertEqual(NTRU_next_Prime_Table4(821), 827)
        self.assertEqual(NTRU_next_Prime_Table4(827), 829)
        self.assertEqual(NTRU_next_Prime_Table4(829), 853)
        self.assertEqual(NTRU_next_Prime_Table4(853), 859)
        self.assertEqual(NTRU_next_Prime_Table4(859), 877)
        self.assertEqual(NTRU_next_Prime_Table4(877), 883)
        self.assertEqual(NTRU_next_Prime_Table4(883), 907)
        self.assertEqual(NTRU_next_Prime_Table4(907), 941)
        self.assertEqual(NTRU_next_Prime_Table4(941), 947)
        self.assertEqual(NTRU_next_Prime_Table4(947), 1019)
        self.assertEqual(NTRU_next_Prime_Table4(1019), 1061)
        self.assertEqual(NTRU_next_Prime_Table4(1061), 1091)
        self.assertEqual(NTRU_next_Prime_Table4(1091), 1109)
        self.assertEqual(NTRU_next_Prime_Table4(1109), 1117)
        self.assertEqual(NTRU_next_Prime_Table4(1117), 1123)
        self.assertEqual(NTRU_next_Prime_Table4(1123), 1171)
        self.assertEqual(NTRU_next_Prime_Table4(1171), 1187)
        self.assertEqual(NTRU_next_Prime_Table4(1187), 1213)
        self.assertEqual(NTRU_next_Prime_Table4(1213), 1229)
        self.assertEqual(NTRU_next_Prime_Table4(1229), 1237)
        self.assertEqual(NTRU_next_Prime_Table4(1237), 1259)
        self.assertEqual(NTRU_next_Prime_Table4(1259), 1277)
        self.assertEqual(NTRU_next_Prime_Table4(1277), 1283)
        self.assertEqual(NTRU_next_Prime_Table4(1283), 1291)
        self.assertEqual(NTRU_next_Prime_Table4(1291), 1301)
        self.assertEqual(NTRU_next_Prime_Table4(1301), 1307)
        self.assertEqual(NTRU_next_Prime_Table4(1307), 1373)
        self.assertEqual(NTRU_next_Prime_Table4(1373), 1381)
        self.assertEqual(NTRU_next_Prime_Table4(1381), 1427)
        self.assertEqual(NTRU_next_Prime_Table4(1427), 1451)
        self.assertEqual(NTRU_next_Prime_Table4(1451), 1453)
        self.assertEqual(NTRU_next_Prime_Table4(1453), 1483)
        self.assertEqual(NTRU_next_Prime_Table4(1483), 1493)
        self.assertEqual(NTRU_next_Prime_Table4(1493), 1499)
        self.assertEqual(NTRU_next_Prime_Table4(1499), 1523)
        self.assertEqual(NTRU_next_Prime_Table4(1523), 1531)
        self.assertEqual(NTRU_next_Prime_Table4(1531), 1549)
        self.assertEqual(NTRU_next_Prime_Table4(1549), 1571)
        self.assertEqual(NTRU_next_Prime_Table4(1571), 1619)
        self.assertEqual(NTRU_next_Prime_Table4(1619), 1621)
        self.assertEqual(NTRU_next_Prime_Table4(1621), 1637)
        self.assertEqual(NTRU_next_Prime_Table4(1637), 1667)
        self.assertEqual(NTRU_next_Prime_Table4(1667), 1669)
        self.assertEqual(NTRU_next_Prime_Table4(1669), 1693)
        self.assertEqual(NTRU_next_Prime_Table4(1693), 1733)
        self.assertEqual(NTRU_next_Prime_Table4(1733), 1741)
        self.assertEqual(NTRU_next_Prime_Table4(1741), 1747)
        self.assertEqual(NTRU_next_Prime_Table4(1747), 1787)
        self.assertEqual(NTRU_next_Prime_Table4(1787), 1861)
        self.assertEqual(NTRU_next_Prime_Table4(1861), 1867)
        self.assertEqual(NTRU_next_Prime_Table4(1867), 1877)
        self.assertEqual(NTRU_next_Prime_Table4(1877), 1901)
        self.assertEqual(NTRU_next_Prime_Table4(1901), 1907)
        self.assertEqual(NTRU_next_Prime_Table4(1907), 1931)
        self.assertEqual(NTRU_next_Prime_Table4(87863764), 87863813)

    def test_NTRU_next_Prime_Table5(self):
        self.assertEqual(NTRU_next_Prime_Table5(100),103)
        self.assertEqual(NTRU_next_Prime_Table5(103),137)
        self.assertEqual(NTRU_next_Prime_Table5(137),167)
        self.assertEqual(NTRU_next_Prime_Table5(167),191)
        self.assertEqual(NTRU_next_Prime_Table5(191),193)
        self.assertEqual(NTRU_next_Prime_Table5(193),199)
        self.assertEqual(NTRU_next_Prime_Table5(199),239)
        self.assertEqual(NTRU_next_Prime_Table5(239),263)
        self.assertEqual(NTRU_next_Prime_Table5(263),271)
        self.assertEqual(NTRU_next_Prime_Table5(271),311)
        self.assertEqual(NTRU_next_Prime_Table5(311),313)
        self.assertEqual(NTRU_next_Prime_Table5(313),359)
        self.assertEqual(NTRU_next_Prime_Table5(359),367)
        self.assertEqual(NTRU_next_Prime_Table5(367),383)
        self.assertEqual(NTRU_next_Prime_Table5(383),401)
        self.assertEqual(NTRU_next_Prime_Table5(401),409)
        self.assertEqual(NTRU_next_Prime_Table5(409),449)
        self.assertEqual(NTRU_next_Prime_Table5(449),463)
        self.assertEqual(NTRU_next_Prime_Table5(463),479)
        self.assertEqual(NTRU_next_Prime_Table5(479),487)
        self.assertEqual(NTRU_next_Prime_Table5(487),503)
        self.assertEqual(NTRU_next_Prime_Table5(503),521)
        self.assertEqual(NTRU_next_Prime_Table5(521),569)
        self.assertEqual(NTRU_next_Prime_Table5(569),599)
        self.assertEqual(NTRU_next_Prime_Table5(599),607)
        self.assertEqual(NTRU_next_Prime_Table5(607),647)
        self.assertEqual(NTRU_next_Prime_Table5(647),719)
        self.assertEqual(NTRU_next_Prime_Table5(719),743)
        self.assertEqual(NTRU_next_Prime_Table5(743),751)
        self.assertEqual(NTRU_next_Prime_Table5(751),761)
        self.assertEqual(NTRU_next_Prime_Table5(761),769)
        self.assertEqual(NTRU_next_Prime_Table5(769),809)
        self.assertEqual(NTRU_next_Prime_Table5(809),823)
        self.assertEqual(NTRU_next_Prime_Table5(823),839)
        self.assertEqual(NTRU_next_Prime_Table5(839),857)
        self.assertEqual(NTRU_next_Prime_Table5(857),863)
        self.assertEqual(NTRU_next_Prime_Table5(863),887)
        self.assertEqual(NTRU_next_Prime_Table5(887),929)
        self.assertEqual(NTRU_next_Prime_Table5(929),967)
        self.assertEqual(NTRU_next_Prime_Table5(967),977)
        self.assertEqual(NTRU_next_Prime_Table5(977),983)
        self.assertEqual(NTRU_next_Prime_Table5(983),991)
        self.assertEqual(NTRU_next_Prime_Table5(991),1009)
        self.assertEqual(NTRU_next_Prime_Table5(1009),1031)
        self.assertEqual(NTRU_next_Prime_Table5(1031),1039)
        self.assertEqual(NTRU_next_Prime_Table5(1039),1063)
        self.assertEqual(NTRU_next_Prime_Table5(1063),1087)
        self.assertEqual(NTRU_next_Prime_Table5(1087),1129)
        self.assertEqual(NTRU_next_Prime_Table5(1129),1151)
        self.assertEqual(NTRU_next_Prime_Table5(1151),1223)
        self.assertEqual(NTRU_next_Prime_Table5(1223),1231)
        self.assertEqual(NTRU_next_Prime_Table5(1231),1279)
        self.assertEqual(NTRU_next_Prime_Table5(1279),1297)
        self.assertEqual(NTRU_next_Prime_Table5(1297),1303)
        self.assertEqual(NTRU_next_Prime_Table5(1303),1319)
        self.assertEqual(NTRU_next_Prime_Table5(1319),1361)
        self.assertEqual(NTRU_next_Prime_Table5(1361),1367)
        self.assertEqual(NTRU_next_Prime_Table5(1367),1409)
        self.assertEqual(NTRU_next_Prime_Table5(1409),1439)
        self.assertEqual(NTRU_next_Prime_Table5(1439),1447)
        self.assertEqual(NTRU_next_Prime_Table5(1447),1487)
        self.assertEqual(NTRU_next_Prime_Table5(1487),1489)
        self.assertEqual(NTRU_next_Prime_Table5(1489),1511)
        self.assertEqual(NTRU_next_Prime_Table5(1511),1543)
        self.assertEqual(NTRU_next_Prime_Table5(1543),1559)
        self.assertEqual(NTRU_next_Prime_Table5(1559),1567)
        self.assertEqual(NTRU_next_Prime_Table5(1567),1583)
        self.assertEqual(NTRU_next_Prime_Table5(1583),1607)
        self.assertEqual(NTRU_next_Prime_Table5(1607),1663)
        self.assertEqual(NTRU_next_Prime_Table5(1663),1697)
        self.assertEqual(NTRU_next_Prime_Table5(1697),1759)
        self.assertEqual(NTRU_next_Prime_Table5(1759),1783)
        self.assertEqual(NTRU_next_Prime_Table5(1783),1823)
        self.assertEqual(NTRU_next_Prime_Table5(1823),1847)
        self.assertEqual(NTRU_next_Prime_Table5(1847),1871)
        self.assertEqual(NTRU_next_Prime_Table5(1871),1873)
        self.assertEqual(NTRU_next_Prime_Table5(1873),1879)
        self.assertEqual(NTRU_next_Prime_Table5(1879),1951)
        self.assertEqual(NTRU_next_Prime_Table5(1951),1993)
        self.assertEqual(NTRU_next_Prime_Table5(1993),2039)

    def test_NTRU_next_Prime(self):
        self.assertEqual(NTRU_next_Prime(100), 101)
        self.assertEqual(NTRU_next_Prime(101), 103)
        self.assertEqual(NTRU_next_Prime(103), 107)
        self.assertEqual(NTRU_next_Prime(107), 131)
        self.assertEqual(NTRU_next_Prime(131), 137)
        self.assertEqual(NTRU_next_Prime(137), 139)
        self.assertEqual(NTRU_next_Prime(139), 149)
        self.assertEqual(NTRU_next_Prime(149), 163)
        self.assertEqual(NTRU_next_Prime(163), 167)
        self.assertEqual(NTRU_next_Prime(167), 173)
        self.assertEqual(NTRU_next_Prime(173), 179)
        self.assertEqual(NTRU_next_Prime(179), 181)
        self.assertEqual(NTRU_next_Prime(181), 191)
        self.assertEqual(NTRU_next_Prime(191), 193)
        self.assertEqual(NTRU_next_Prime(193), 197)
        self.assertEqual(NTRU_next_Prime(197), 199)
        self.assertEqual(NTRU_next_Prime(199), 211)
        self.assertEqual(NTRU_next_Prime(211), 227)
        self.assertEqual(NTRU_next_Prime(227), 239)
        self.assertEqual(NTRU_next_Prime(239), 263)
        self.assertEqual(NTRU_next_Prime(263), 269)
        self.assertEqual(NTRU_next_Prime(269), 271)
        self.assertEqual(NTRU_next_Prime(271), 293)
        self.assertEqual(NTRU_next_Prime(293), 311)
        self.assertEqual(NTRU_next_Prime(311), 313)
        self.assertEqual(NTRU_next_Prime(313), 317)
        self.assertEqual(NTRU_next_Prime(317), 347)
        self.assertEqual(NTRU_next_Prime(347), 349)
        self.assertEqual(NTRU_next_Prime(349), 359)
        self.assertEqual(NTRU_next_Prime(359), 367)
        self.assertEqual(NTRU_next_Prime(367), 373)
        self.assertEqual(NTRU_next_Prime(373), 379)
        self.assertEqual(NTRU_next_Prime(379), 383)
        self.assertEqual(NTRU_next_Prime(383), 389)
        self.assertEqual(NTRU_next_Prime(389), 401)
        self.assertEqual(NTRU_next_Prime(401), 409)
        self.assertEqual(NTRU_next_Prime(409), 419)
        self.assertEqual(NTRU_next_Prime(419), 421)
        self.assertEqual(NTRU_next_Prime(421), 443)
        self.assertEqual(NTRU_next_Prime(443), 449)
        self.assertEqual(NTRU_next_Prime(449), 461)
        self.assertEqual(NTRU_next_Prime(461), 463)
        self.assertEqual(NTRU_next_Prime(463), 467)
        self.assertEqual(NTRU_next_Prime(467), 479)
        self.assertEqual(NTRU_next_Prime(479), 487)
        self.assertEqual(NTRU_next_Prime(487), 491)
        self.assertEqual(NTRU_next_Prime(491), 503)
        self.assertEqual(NTRU_next_Prime(503), 509)
        self.assertEqual(NTRU_next_Prime(509), 521)
        self.assertEqual(NTRU_next_Prime(521), 523)
        self.assertEqual(NTRU_next_Prime(523), 541)
        self.assertEqual(NTRU_next_Prime(541), 547)
        self.assertEqual(NTRU_next_Prime(547), 557)
        self.assertEqual(NTRU_next_Prime(557), 563)
        self.assertEqual(NTRU_next_Prime(563), 569)
        self.assertEqual(NTRU_next_Prime(569), 587)
        self.assertEqual(NTRU_next_Prime(587), 599)
        self.assertEqual(NTRU_next_Prime(599), 607)
        self.assertEqual(NTRU_next_Prime(607), 613)
        self.assertEqual(NTRU_next_Prime(613), 619)
        self.assertEqual(NTRU_next_Prime(619), 647)
        self.assertEqual(NTRU_next_Prime(647), 653)
        self.assertEqual(NTRU_next_Prime(653), 659)
        self.assertEqual(NTRU_next_Prime(659), 661)
        self.assertEqual(NTRU_next_Prime(661), 677)
        self.assertEqual(NTRU_next_Prime(677), 701)
        self.assertEqual(NTRU_next_Prime(701), 709)
        self.assertEqual(NTRU_next_Prime(709), 719)
        self.assertEqual(NTRU_next_Prime(719), 743)
        self.assertEqual(NTRU_next_Prime(743), 751)
        self.assertEqual(NTRU_next_Prime(751), 757)
        self.assertEqual(NTRU_next_Prime(757), 761)
        self.assertEqual(NTRU_next_Prime(761), 769)
        self.assertEqual(NTRU_next_Prime(769), 773)
        self.assertEqual(NTRU_next_Prime(773), 787)
        self.assertEqual(NTRU_next_Prime(787), 797)
        self.assertEqual(NTRU_next_Prime(797), 809)
        self.assertEqual(NTRU_next_Prime(809), 821)
        self.assertEqual(NTRU_next_Prime(821), 823)
        self.assertEqual(NTRU_next_Prime(823), 827)
        self.assertEqual(NTRU_next_Prime(827), 829)
        self.assertEqual(NTRU_next_Prime(829), 839)
        self.assertEqual(NTRU_next_Prime(839), 853)
        self.assertEqual(NTRU_next_Prime(853), 857)
        self.assertEqual(NTRU_next_Prime(857), 859)
        self.assertEqual(NTRU_next_Prime(859), 863)
        self.assertEqual(NTRU_next_Prime(863), 877)
        self.assertEqual(NTRU_next_Prime(877), 883)
        self.assertEqual(NTRU_next_Prime(883), 887)
        self.assertEqual(NTRU_next_Prime(887), 907)
        self.assertEqual(NTRU_next_Prime(907), 929)
        self.assertEqual(NTRU_next_Prime(929), 941)
        self.assertEqual(NTRU_next_Prime(941), 947)
        self.assertEqual(NTRU_next_Prime(947), 967)
        self.assertEqual(NTRU_next_Prime(967), 977)
        self.assertEqual(NTRU_next_Prime(977), 983)
        self.assertEqual(NTRU_next_Prime(983), 991)
        self.assertEqual(NTRU_next_Prime(991), 1009)
        self.assertEqual(NTRU_next_Prime(1009), 1019)
        self.assertEqual(NTRU_next_Prime(1019), 1031)
        self.assertEqual(NTRU_next_Prime(1031), 1039)
        self.assertEqual(NTRU_next_Prime(1039), 1061)
        self.assertEqual(NTRU_next_Prime(1061), 1063)
        self.assertEqual(NTRU_next_Prime(1063), 1087)
        self.assertEqual(NTRU_next_Prime(1087), 1091)
        self.assertEqual(NTRU_next_Prime(1091), 1109)
        self.assertEqual(NTRU_next_Prime(1109), 1117)
        self.assertEqual(NTRU_next_Prime(1117), 1123)
        self.assertEqual(NTRU_next_Prime(1123), 1129)
        self.assertEqual(NTRU_next_Prime(1129), 1151)
        self.assertEqual(NTRU_next_Prime(1151), 1171)
        self.assertEqual(NTRU_next_Prime(1171), 1187)
        self.assertEqual(NTRU_next_Prime(1187), 1213)
        self.assertEqual(NTRU_next_Prime(1213), 1223)
        self.assertEqual(NTRU_next_Prime(1223), 1229)
        self.assertEqual(NTRU_next_Prime(1229), 1231)
        self.assertEqual(NTRU_next_Prime(1231), 1237)
        self.assertEqual(NTRU_next_Prime(1237), 1259)
        self.assertEqual(NTRU_next_Prime(1259), 1277)
        self.assertEqual(NTRU_next_Prime(1277), 1279)
        self.assertEqual(NTRU_next_Prime(1279), 1283)
        self.assertEqual(NTRU_next_Prime(1283), 1291)
        self.assertEqual(NTRU_next_Prime(1291), 1297)
        self.assertEqual(NTRU_next_Prime(1297), 1301)
        self.assertEqual(NTRU_next_Prime(1301), 1303)
        self.assertEqual(NTRU_next_Prime(1303), 1307)
        self.assertEqual(NTRU_next_Prime(1307), 1319)
        self.assertEqual(NTRU_next_Prime(1319), 1361)
        self.assertEqual(NTRU_next_Prime(1361), 1367)
        self.assertEqual(NTRU_next_Prime(1367), 1373)
        self.assertEqual(NTRU_next_Prime(1373), 1381)
        self.assertEqual(NTRU_next_Prime(1381), 1409)
        self.assertEqual(NTRU_next_Prime(1409), 1427)
        self.assertEqual(NTRU_next_Prime(1427), 1439)
        self.assertEqual(NTRU_next_Prime(1439), 1447)
        self.assertEqual(NTRU_next_Prime(1447), 1451)
        self.assertEqual(NTRU_next_Prime(1451), 1453)
        self.assertEqual(NTRU_next_Prime(1453), 1483)
        self.assertEqual(NTRU_next_Prime(1483), 1487)
        self.assertEqual(NTRU_next_Prime(1487), 1489)
        self.assertEqual(NTRU_next_Prime(1489), 1493)
        self.assertEqual(NTRU_next_Prime(1493), 1499)
        self.assertEqual(NTRU_next_Prime(1499), 1511)
        self.assertEqual(NTRU_next_Prime(1511), 1523)
        self.assertEqual(NTRU_next_Prime(1523), 1531)
        self.assertEqual(NTRU_next_Prime(1531), 1543)
        self.assertEqual(NTRU_next_Prime(1543), 1549)
        self.assertEqual(NTRU_next_Prime(1549), 1559)
        self.assertEqual(NTRU_next_Prime(1559), 1567)
        self.assertEqual(NTRU_next_Prime(1567), 1571)
        self.assertEqual(NTRU_next_Prime(1571), 1583)
        self.assertEqual(NTRU_next_Prime(1583), 1607)
        self.assertEqual(NTRU_next_Prime(1607), 1619)
        self.assertEqual(NTRU_next_Prime(1619), 1621)
        self.assertEqual(NTRU_next_Prime(1621), 1637)
        self.assertEqual(NTRU_next_Prime(1637), 1663)
        self.assertEqual(NTRU_next_Prime(1663), 1667)
        self.assertEqual(NTRU_next_Prime(1667), 1669)
        self.assertEqual(NTRU_next_Prime(1669), 1693)
        self.assertEqual(NTRU_next_Prime(1693), 1697)
        self.assertEqual(NTRU_next_Prime(1697), 1733)
        self.assertEqual(NTRU_next_Prime(1733), 1741)
        self.assertEqual(NTRU_next_Prime(1741), 1747)
        self.assertEqual(NTRU_next_Prime(1747), 1759)
        self.assertEqual(NTRU_next_Prime(1759), 1783)
        self.assertEqual(NTRU_next_Prime(1783), 1787)
        self.assertEqual(NTRU_next_Prime(1787), 1823)
        self.assertEqual(NTRU_next_Prime(1823), 1847)
        self.assertEqual(NTRU_next_Prime(1847), 1861)
        self.assertEqual(NTRU_next_Prime(1861), 1867)
        self.assertEqual(NTRU_next_Prime(1867), 1871)
        self.assertEqual(NTRU_next_Prime(1871), 1873)
        self.assertEqual(NTRU_next_Prime(1873), 1877)
        self.assertEqual(NTRU_next_Prime(1877), 1879)
        self.assertEqual(NTRU_next_Prime(1879), 1901)
        self.assertEqual(NTRU_next_Prime(1901), 1907)
        self.assertEqual(NTRU_next_Prime(1907), 1931)
        self.assertEqual(NTRU_next_Prime(1931), 1949)
        self.assertEqual(NTRU_next_Prime(1949), 1951)

    def testNTRU_P_N(self):
        self.assertEqual(NTRU_P_N( 20 , 3 , 3 , 3 ), 465844843008000000 )
        self.assertEqual(NTRU_P_N( 20 , 3 , 3 , 4 ), 5298985089216000000 )
        self.assertEqual(NTRU_P_N( 20 , 3 , 4 , 3 ), 5298985089216000000 )
        self.assertEqual(NTRU_P_N( 20 , 4 , 3 , 3 ), 5298985089216000000 )
        self.assertEqual(NTRU_P_N( 20 , 4 , 3 , 4 ), 60275955389832000000 )
        self.assertEqual(NTRU_P_N( 20 , 4 , 3 , 5 ), 318257044458312960000 )
        self.assertEqual(NTRU_P_N( 20 , 5 , 5 , 4 ), 19114518090166276377600 )
        self.assertEqual(NTRU_P_N( 20 , 5 , 5 , 5 ), 100924655516077939273728 )
        self.assertEqual(NTRU_P_N( 100 , 3 , 3 , 3 ), 13551146061118253102592000000000 )
        self.assertEqual(NTRU_P_N( 100 , 3 , 5 , 4 ), 1354711683550829297666917306602240000000 )
        self.assertEqual(NTRU_P_N( 100 , 4 , 5 , 5 ), 247871678021996249643873719604050248934400000 )
        self.assertEqual(NTRU_P_N( 349 , 3 , 3 , 3 ), 111085108640464020924562731096033494528000 )
        self.assertEqual(NTRU_P_N( 349 , 4 , 4 , 4 ), 43778066278602474284200956502132121375470866421263000 )
        self.assertEqual(NTRU_P_N( 349 , 5 , 4 , 3 ), 27691700021979688642069166913786224212351789693996800 )
        self.assertEqual(NTRU_P_N( 349 , 5 , 5 , 4 ), 941549482821304560792127066402597686039082201278480101738880 )
        self.assertEqual(NTRU_P_N( 349 , 5 , 5 , 5 ),
                 4366529881532082031129568483148687028774847616649079319824229888 )

    def test_NTRU_shannon_p(self):

        self.assertAlmostEqual(NTRU_shannon_p(401, 131, 113),  97.377083840317757911093337789895373407334)
        self.assertAlmostEqual(NTRU_shannon_p(401, 154, 133), 117.64544671566375467936245484379700328014)
        self.assertAlmostEqual(NTRU_shannon_p(439, 174, 146), 133.42409866220694)
        self.assertAlmostEqual(NTRU_shannon_p(539, 261, 197), 200.24555924629416704288801661072368736311)
        #self.assertAlmostEqual(NTRU_shannon_p(743, 350, 247), 117.643379324493)	# too big

    def test_NTRU_log2_eta(self):
        self.assertAlmostEqual(NTRU_log2_eta(401, 11, 113, 131), 1.0083824719620765407)

    def testNTRU_wmitm(self):
        self.assertAlmostEqual(NTRU_wmitm(541, 11, 180, 144, 128)[0], 166.55629436501071594517972168760542593151)


def NTRU_test(N, K, q, Lambda):

    dg = int(math.floor(float(N)/float(3) + 0.5))
    d1 = int(math.ceil(0.25*(math.sqrt(1 + float(8*N)/float(3)) - 1)))
    d2 = int(math.ceil( ( float(N)/float(3) - d1 ) / float(2*d1) ))
    d3 = int(max(math.ceil(float(d1)*0.5 + 1), math.ceil(float(N)/float(3) - 2*d1*d2)))
    dm = NTRU_dm(N)

    print 4*d1*d2 + 2*d3, int(2.0*float(N)/3.0)

    if(d2 == 8): d3 = 5
    if(d2 == 10): d3 = 8
    if(d2 == 11): d3 = 15

    p = NTRU_shannon_p(N, K, dg)
    eta = NTRU_log2_eta(N, int(math.log(q, 2)), Lambda, K)
    sigma = math.sqrt( (4*d1*d2 + 2*d3) * (float(N - dm + 2*dg + 1)/float(N)) )
    pdf = math.log(NTRU_erfc(N, q, sigma), 2)
    ps = math.log(NTRU_product_form_search(N, d1, d2, d3), 2)

    print "N = ", N, "K =", K, "q = ", q, "lambda =", Lambda
    print "d1 =", d1, "d2 =", d2, "d3 =", d3, "dg =", dg, "dm =", dm, "Cost =", round(p, 2), "eta =", \
        round(eta, 4), "pdf =", round(pdf,2), "product =", round(ps, 2)
    print 4*d1*d2 + 2*d3, int(2.0*float(N)/3.0)
    print ""


def main(argv):

    # parse options
    try:
        opts, args = getopt.getopt(argv, "uk:d:", ["unit-test","k-bit-security=","debug="])
    except getopt.GetoptError:
        print "Options not understood"
        sys.exit()

    # initialise globals
    global k
    global debug
    unit = False
    for i in range (0, len(opts)):
        opt, arg = opts[i]
        if opt in ("-k", "--k-bit-security"):
            if(int(arg) in security_levels):
                k = int(arg)
            else:
                print "Security", int(arg), "is not supported"
                sys.exit()
        elif opt in ("-d", "--debug"):
            debug = int(arg)
            if(debug > 0): print "( debug level", debug, ")"
        elif opt in ("-u", "--unit-test"):
            unit = True
        else:
            print "Option", opt, "not understood"
            sys.exit()

    # delete parsed options so as not to confuse unit tests
    sys.argv = [ sys.argv[0] ]

    # run calculation
    if(unit):
        unittest.main()
    else:
        print "Security level", k
        NTRU_parameters(k)

if __name__ == '__main__':
    main(sys.argv[1:])
    bkz(601, 331, 1, 50) # to compare with the GP implementation


