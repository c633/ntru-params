read("bkzsim.gp");

log2(x) = {if(x == 0, 0, log(x)/log(2))};

binsearchLT(f, val, low, high) = {
  my(mp);
  if(low >= high-1, return(low));
  mp = ceil((low+high)/2);
  if(f(mp) <= val,
    binsearchLT(f, val, mp, high),
    binsearchLT(f, val, low, mp));
}
addhelp(binsearchLT,\
"binsearchLT(f, val, low, high): "\
"Binary search for largest x satisfying f(x) < val in range [low, high]");


dmRejectProb(N,dm) = {
  my(a);
  a = 0;
  for(p1s=dm, N-2*dm, \
  for(m1s=dm, N-dm-p1s, \
      a += binomial(N, p1s)*binomial(N-p1s,m1s)));
  log2(1 - a/3^N);
}
addhelp(dmRejectProb,\
"dmRejectProb(N, dm): "\
"Probability that the number of 1s, -1s, or 0s is less than dm in a uniform "\
"random trinary polynomial of degree <= N");


hybridMITM(N, K, dg1, dg2) = {
    my(Nc, tot, H, awork, d1, d2, d1a, d2a);
    Nc = N-K;
    tot = (binomial(N,dg1)*binomial(N-dg1,dg2));
    H = vector(dg1+1, d1a,
               vector(dg2+1, d2a,
               d1 = d1a - 1;
               d2 = d2a - 1;
               p = binomial(Nc, dg1-d1) * binomial(Nc - dg1 + d1, dg2-d2) / tot;
               -(binomial(K, d1) * binomial(K-d1, d2)) * p * log2(p)));
    awork = .5*(sum(i=0, dg1, sum(j=0, dg2, H[i+1][j+1])) - log2(N));
    \\printf("HybridMITM: K: %d Work: %d\n", K, awork);
    awork;
}

minHybridMITM(N, K, dm) = {
    my(est,tot,t);
    c = N-3*dm;
    t = hybridMITM(N, K, dm, dm+c) + .5*log2(N);
}
addhelp(minHybridMITM,\
"minHybridMITM(N, K, dm): Calculate cost of performing hybrid attack on"\
"a message with maximal information leakage via m(1), i.e. m(1) = N-3dm");


hybridHermiteRequirement(N,q,L,K) = {
  ld = (N - L)*log2(q);
  ld /= (4*N^2 - 4*N*(K+L) + (K^2 + 2*K*L + L^2));
  ld -= 1/(2*N - (K+L));
  2^ld
}
addhelp(hybridHermiteRequirement,\
"hybridHermiteRequirement(N,q,L,K): "\
"Root Hermite factor required to prepare an NTRU basis for hybrid "\
"meet-in-the-middle attack. [L, N-K] is column index range of block "\
"to be reduced. L should be taken equal to the security parameter");


deltaStar(k) = {
  \\if(k <= 50, return(1.01),
  \\if(k <= 60, return(1.009),
  \\if(k <= 80, return(1.008),
  \\if(k <= 128, return(1.007),
  \\if(k <= 192, return(1.006),
  \\if(k <= 256, return(1.005),
  \\return(1) ))))));
  return(1.91157e-8*k^2 - 2.32633e-5*k + 1.00972);
  \\return(9.09036e-8*k^2 - 5.7995*k + 1.01419);
}
addhelp(deltaStar, \
"deltaStar(k): Conjectured root Hermite factor reachable by BKZ "\
"with 2^k operations.");


decFailSig(pm, pg, d1,d2,d3) = {
  sig = 3 * sqrt((4*d1*d2 + 2*d3)*pm + (4*d1*d2 + 2*d3)*pg);
}
addhelp(decFailSig, \
"decFailSig(pm,pg,d1,d2,d3): returns expected standard deviation "\
"of a coefficient of 3*(r*g + m*F) + m.")

numProdForm(N, a, b, c) = {
  my(S);
  S = binomial(N,a)*binomial(N-a,a) * binomial(N,b)*binomial(N-b,b) * binomial(N,c)*binomial(N-c,c);
  return(S);
}

blockSize(N, q) = {
  h = sqrt(q)^(1/(2*N));
  binsearchLT((x)->(-(simulate(2*N, x, h)[2])), -h, 60, N);
}

bkzCost(dim, bs, iter) = {
  my(logNodes);
  \\ Quad fit to published BKZ-2.0 paper table 3 row 1
  \\logNodes1 = 0.00405892*bs^2 - 0.337913*bs + 34.9018;
  \\ Quad fit to Full BKZ-2.0 paper table 4
  logNodes2 = 0.000784314*bs^2 + 0.366078*bs - 6.125;
  round(logNodes2 + log2(dim*iter) + 7);
}

cn11est(dim, hreq) = {
  my(bs, iter, cost);

  bs = binsearchLT((x)->(-simulate(dim, x, hreq)[2]), -hreq, 60, dim) + 1;
  iter = simulate(dim, bs, hreq)[3];
  cost = bkzCost(dim, bs, iter);
  [iter, bs, cost];
}

getDs(N) = {
  my(d1,d2,d3);
  d1 = d2 = d3 = ceil(vecmax(real(polroots(2*x^2 + x - N/3))));
  d2 = ceil((N/3 - d1)/(2*d1));
  d3 = max(ceil((d1 + 1)/2), ceil(N/3 - 2*d1*d2));
  [d1,d2,d3]
}

genParams(N, verbose=0) = {
  my(lambda, directMITM, dm, d1, d2, d3, dg, sig,\
  q, q2, decFail, decFail2, Kh, Kc, Kb, LL, LM, high, low, out);

  /* Standard choices for dg, d1, d2, and d3 */
  dg = round(N/3);
  [d1,d2,d3] = getDs(N);

  /* Pick initial dm based on rejection probability below 2^-10 */
  dm = binsearchLT((x)->(dmRejectProb(N,x)), -10, 0, floor(N/3));
  mRej = dmRejectProb(N,dm);

  /* Use direct product-form MITM for upper bound on security */
  directMITM = floor(0.5 * log2(numProdForm(N,d1,d2,d3)/N));

  /* Choose q as smallest power of 2 admitting negligible
     decryption failure probability */
  sig = decFailSig((1-dm/N), (2*dg + 1)/N, d1, d2, d3);
  q = binsearchLT((x)->(-log2(N*erfc(x/(sig*sqrt(2))))), directMITM, 2^4, 2^16);
  q = 2^ceil(log2(2*q + 2));
  decFail = round(log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));

  /* Kh is smallest K that can be prepared via lattice reduction in O(2^\lambda) time. */
  [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm);

  /* Redo search for q with new security estimate */
  /* TODO: This favors estimate 1 */
  q2 = binsearchLT((x)->(-log2(N*erfc(x/(sig*sqrt(2))))), lambda, 2^4, 2^16);
  q2 = 2^ceil(log2(2*q2 + 2));
  decFail2 = round(log2(N*erfc(((q2-2)/2)/(sig*sqrt(2)))));
  if(q2 != q, /* If we can lower q, rederive security */
    q = q2;
    decFail = decFail2;
    [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm));

  /* Update security estimate. Either keep direct MITM cost or
   * choose larger of the two costs from hybrid attack (which should
   * be roughly equal anyway). */
  \\lambda = min(lambda, max(LM,LL));

  if(verbose,
    checkParams([N, q, d1, d2, d3, dg, dm, lambda, K]));

  out = [N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM];
  printf("[N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM]\n%s\n", out);

  out
}

optimalK(N,q,d1,d2,d3,dg,dm) = {
  my(lambda, Kb, Llr, Lmitm, high, low, diff);
  lambda = floor(0.5 * log2(numProdForm(N,d1,d2,d3)/N));

  low = 0;
  high = N;
  Kb = ceil((low + high)/2);
  Llr = cn11est(2*N - Llr - Kb, hybridHermiteRequirement(N, q, Llr, Kb))[3];
  Lmitm = floor(hybridMITM(N, Kb, dg+1, dg));
  diff = abs(Llr - Lmitm);
  if(Llr < Lmitm,
     high = Kb;
     low = Kb - 2*floor(diff/2),
     low = Kb;
     high = low + 2*floor(diff/2));

  /* Search for a K between Kc and Kh that balances the two costs */
  Kb = ceil((low + high)/2);
  while(high-low > 1, /* Binary search for balanced costs */
    Llr = cn11est(2*N - Llr - Kb, hybridHermiteRequirement(N, q, Llr, Kb))[3];
    Lmitm = floor(hybridMITM(N, Kb, dg+1, dg));
    diff = abs(Llr - Lmitm);
    if(Llr < Lmitm,
       high = Kb;
       low = Kb - 2*floor(diff/2),
       low = Kb;
       high = low + 2*floor(diff/2));
    Kb = ceil((low + high)/2));

  Lmsg = floor(minHybridMITM(N, Kb, dm));

  /* Update security estimate. Either keep direct MITM cost or
   * choose smaller (message or key recovery) of the hybrid attack costs.  */
  lambda = min(lambda, max(min(Lmitm, Lmsg),Llr));

  [lambda, Kb];
}



/* Extra parameters needed for reference implementation */

probKUniq(N, K, M) = {
  my(Z, T, R, S, RHS14, JT2);
  Z = 1.0*(K/N);
  T = 1.0*(M/N);

  /* Equation 3 of reference */
  R = solve(X=0.01,(1/Z)-(1e-9),-1/X*log(1-X*Z)-T);
  /* Equation 12 */
  S = sqrt(Z/(1-R*Z) - T);
  /* RHS of Equation 14 */
  RHS14 = (2*Pi*S^2)^(-1/2) * (R/(R-1)) * sqrt((1 - R*Z)/(1-Z));
  /* Equation 2 */
  JT2 = (T - Z)*log(R) + (1-Z)*log(1-Z) - (1-R*Z)/(R)*log(1-R*Z);
  /* Probability from equation 14 */
  -log2(RHS14 / (exp(N*JT2)*sqrt(N)));
}
addhelp(probKUniq, \
"probKUniq(N,K,M): Probability that a set of M integers uniformly " \
"distributed in [1,N] contains K unique values. Expressed as -log2(prob). "\
"Estimate from: Dupuis, Zhang, Whiting, \"Refined Large Deviation Asymptotics "\
"for the Classical Occupancy Problem.\" 2006");


logBinomialCDF(k, n, p) = {
  my(i);
  log2(sum(i=0,k,binomial(n,i)*(p)^i * (1-p)^(n-i)));
}
addhelp(logBinomialCDF, \
"logBinomialCDF(k,n,p): log2(Pr(X <= k)), X binomial with parameters n and p.")


altOccupancy(C,N,n,L,d) = {
  global(grid);
  if(type(grid) != "t_VEC" || grid[1] != [C,N,n] || grid[2] < L,
  /* Initialize storage for memoization */
  grid = [[C,N,n],L,matrix(L+2,L+2,iL,id,if(iL < id, 0, if(id == 1, (1.0 - n*N/C)^(iL-1), -1)))]);

  if(grid[3][L+1, d+1] != -1,
    /* Reuse previously calculated entry */
    return(grid[3][L+1,d+1]));

  /* Calculate a new entry */
  grid[3][L+1,d+1] = altOccupancy(C,N,n,L-1,d-1) * (1.0*n*(N-d+1)/C) \
                   + altOccupancy(C,N,n,L-1,d)   * (1 - n*(N-d)/C);
  return(grid[3][L+1,d+1]);
}
addhelp(altOccupancy, \
"altOccupancy(C,N,n,L,d): probability that a list of L random numbers chosen"\
"uniformly in [0,C) contains exactly d values in [0, nN] that are distinct"\
"modulo N. As in Silverman and Whyte 2007.");

altOccupancy2(C,N,n,L,d) = {
  global(grid);
  if(type(grid) != "t_VEC" || grid[1] != [C,N,n] || grid[2] < L,
  /* Initialize storage for memoization */
    grid = [[C,N,n], L, matrix(L+2,L+2,iL,id,
            if(iL < id, 0,
            if(id == 1, (1.0 - n*N/C)^(iL-1),
                        -1)))
           ]);

  if(grid[3][L+1, d+1] != -1,
    /* Reuse previously calculated entry */
    return(grid[3][L+1,d+1]));

  /* Calculate a new entry */
  grid[3][L+1,d+1] = altOccupancy(C,N,n,L-1,d-1) * (1.0*n*(N-d+1)/C) \
                   + altOccupancy(C,N,n,L-1,d)   * (1 - n*(N-d)/C);
  return(grid[3][L+1,d+1]);
}

minIndexSamples(N, c, need, lambda) = {
  my(i, cm, minSamp, err, minRand);
  minSamp = binsearchLT((x)->probKUniq(N, need, x), lambda, need+1, 10*need);
  cm = 2^c - lift(Mod(2^c, N));
  err = 1.0 * cm/2^c;
  minRand = binsearchLT((s)->(-logBinomialCDF(minSamp, s, err)), lambda, minSamp, floor(N*log(N)));

  /* Use upper bound to initialize memoization array for exact computation */
  altOccupancy(2^c, N, 2^c\N, minRand, need);
  minRand = need;
  while(-log2(sum(i=0, need-1, altOccupancy(2^c, N, 2^c\N, minRand, i))) < lambda, minRand+=1);

  minRand;
}

minMGFCalls(N, hashLen, lambda) = {
  my(n5, minMGF);
  n5 = ceil(N/5);
  minMGF = binsearchLT((x)->(-logBinomialCDF(n5, x, 243/256)), lambda, n5, 2*n5);
  minMGF = ceil(minMGF/hashLen);
}


formatParams(genoutput) = {
  my(lambda, c, cs, cm, lLen, secOct, mLen, hashLen, minSamp, minRand, err);
  [N, q, d1, d2, d3, dg, dm, lambda, K] = genoutput;
  c = ceil(log2(N));
  for(cs=c, 13, if(lift(Mod(2^cs, N))/N < lift(Mod(2^c, N))/N, c = cs));
  cm = 2^c - lift(Mod(2^c, N));

  /* Upper bound on security */
  \\lambda = round(0.5 * log2(numProdForm(N,d1,d2,d3)/N));

  secOct = min(32, floor(lambda/8));
  /* max message length (bytes) using 3 bit to 2 trit encoding
     assumes mLen will be < 256, assumes b is secOct bytes */
  mLen = floor(((N-1)/2)*3/8) - min(32, 2*secOct);
  lLen = 1 + floor(log(mLen)/log(256));
  mLen -= lLen;

  /* hashLen = 20; if(secOct > 20, hashLen = 32);*/
  hashLen = 32;

  /* Calculate minIGFCalls */
  need = 2*(d1+d2+d3);
  minRand = minIndexSamples(N, c, need, lambda);
  minCalls = ceil(minRand*c/(8*hashLen));

  /* minMGF calls assuming 5 trits are extracted from 1 byte with probability 243/256 */
  minMGF = minMGFCalls(N, hashLen, lambda);

printf( \
"    {\n"
"        NTRU_EES%dEPX,\t\t\t/* parameter-set id */\n" \
"        \"ees%depX\",\t\t\t/* human readable param set name */\n" \
"        {0xFF, 0xFF, 0xFF},\t\t/* OID */\n" \
"        0xFF,\t\t\t\t/* DER id */\n" \
"        %d,\t\t\t\t/* no. of bits in N (i.e., in an index) */\n" \
"        %d,\t\t\t\t/* N */\n" \
"        %d,\t\t\t\t/* security strength in octets */\n" \
"        %d,\t\t\t\t/* no. of octets for random string b */\n"
"        %d,\t\t\t\t/* q */\n" \
"        %d,\t\t\t\t/* no. of bits in q (i.e., in a coeff) */\n" \
"        TRUE,\t\t\t\t/* product form */\n" \
"        %d + (%d << 8) + (%d << 16),\t/* df, dr */\n" \
"        %d,\t\t\t\t/* dg */\n" \
"        %d,\t\t\t\t/* maxMsgLenBytes */\n" \
"        %d,\t\t\t\t/* dm0 */\n" \
"        %d,\t\t\t\t/* 2^c - (2^c mod N) */\n" \
"        %d,\t\t\t\t/* c */\n" \
"        1,\t\t\t\t/* lLen */\n" \
"        %d,\t\t\t\t/* min. no. of hash calls for IGF-2 */\n" \
"        %d,\t\t\t\t/* min. no. of hash calls for MGF-TP-1 */\n" \
"        NTRU_CRYPTO_HASH_ALGID_SHA256,\t/* hash function for MGF-TP-1, HMAC-DRBG, etc. */\n"
"    },\n", N,N,ceil(log2(N)), N, secOct, min(32, 2*secOct), q, ceil(log2(q)), \
d1, d2, d3, dg, mLen, dm, cm, c, minCalls, minMGF);
}

checkParams(genoutput) = {
  my(goodNQ, hB, L, K, bCombSec, bMSec, sig, decFail, mRej);
  [N, q, d1, d2, d3, dg, dm, lambda, K] = genoutput;

  goodNQ = if(isprime(N),
              pAbove2 = (N-1)/znorder(Mod(2,N));
              if(isprimepower(q) && Mod(q,2) == 0,
              if(pAbove2 < 3, "Yes.",
                Strprintf("No. %d primes above 2.", pAbove2)),
                "No. q is not a power of 2."),
                "No. Composite N.");

  /* Upper bound on security */
  directMITM = floor(0.5 * log2(numProdForm(N,d1,d2,d3)/N));

  hB = hybridHermiteRequirement(N,q,lambda,K);

  bCombSec = hybridMITM(N, K, dg+1, dg);
  bMSec = minHybridMITM(N, K, dm);

  preB = cn11est(2*N - lambda - K, hB);

  sig = decFailSig((1-dm/N), 2/3, d1, d2, d3);
  decFail = round(log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));
  mRej = dmRejectProb(N,dm);

  printf("[N, q, d1, d2, d3, dg, dm] = %s\n", [N,q,d1,d2,d3,dg,dm]);
  printf("Valid N and q? %s\n", goodNQ);
  printf("Decryption failure prob. = %.1f\n", decFail);
  printf("Message rejection prob. = %.1f\n\n", mRej);
  printf("Security Estimates\n\n");

  printf("Direct MITM search cost = %d\n", directMITM);

  printf("Hybrid attack\n\tSearch K = %d coordinates\n\tMust reach root Hermite factor = %.4f\n", K, hB);
  printf("\tCN11 estimates %d rounds of BKZ-%d. Total cost = %d\n", preB[1], preB[2], preB[3]);
  printf("\tHybrid MITM cost = %d\n\tHybrid MITM cost [msg with max m(1)] = %d\n\n", bCombSec, bMSec);
}

P = [107, 113, 131, 139, 149, 163, 173, 181, 191, 199, 211, 227, 239, 251, 263, 271, 281, 293, 307, 317, 331, 347, 359, 367, 379, 389, 401, 439, 593, 743]

