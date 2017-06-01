\\ BKZ simulator.
\\ usage: simulate( dimension, blocksize, target hermite factor )
\\  simulates BKZ-blocksize on a unit volume lattice to which BKZ-20
\\  has been applied as preprocessing.
\\ returns ["success"/"failure", achieved hermite factor, # of rounds required ]

simdata(really=0) = {
    \\ average log(gram schmidt norm) over 100 HKZ reduced bases of dim 50 lattices
    [0.4809337322749968, 0.4889068929146757, 0.4629910732303647, 0.4384921120061095, 0.4198271756529734,
    0.3940124751357192, 0.3793579556691379, 0.3552017168415738, 0.3375032857978846, 0.3229996676156046,
    0.3103169826524305, 0.2978627511364960, 0.2828082600293407, 0.2685092222965025, 0.2470246073218571,
    0.2345601366183950, 0.2236298423327614, 0.2026125221670087, 0.1833511717333619, 0.1635239915325074,
    0.1460909610754462, 0.1239402813211751, 0.1033442833745716, 0.08072183355489210, 0.05747352858422083,
    0.03615285314640355, 0.009734731674006085, -0.01314276679308946, -0.03859536413875225, -0.06166664730992491,
    -0.08732858253410711, -0.1159733213895935, -0.1395873057069733, -0.1685959449423031, -0.2009987452911466,
    -0.2272943479144534, -0.2548892487960738, -0.2845907037340676, -0.3130406180111631, -0.3439519155213564,
    -0.3729166620199606, -0.4037203497626708, -0.4279121623225402, -0.4591242077605871, -0.4851668230787535,
    -0.5069333755274962, -0.5312523582495852, -0.5480002333962808, -0.5470408985906416, -0.5201614648988958]
}

simulate(N, beta, target, abort=50) = {
  my(r, c, ll, vs, llp, R, \
    phi, d, f, logV);
  if(beta < 50, return(0));
  r = simdata(0);
  c = vector(beta);
  for(d=51, #c, c[d] = (1/d)*lngamma(d/2 + 1) - 1/2*log(Pi));

  \\ Start with BKZ-20 preprocessing.
  ll  = vector(N, i, (N-2*(i-1))*log(1.01263));
  vs = vecsum(ll)/N;
  for(i=1,N,ll[i] -= vs);
  llp = vector(N);

  R = 0;
  while(exp(ll[1]/N) > target && R < abort,
    phi = 1;
    for(k=1, N-50,
      d = min(beta,N-k+1);
      f = min(k+beta, N);
      logV = sum(i=1,f,ll[i]) - sum(i=1,k-1,llp[i]);
      if(phi,
          if(logV/d + c[d] < ll[k],
            llp[k] = logV/d + c[d];
            phi = 0),
          llp[k] = logV/d + c[d]);
      );
    logV = vecsum(ll) - vecsum(llp[1..(N-50)]);
    for(k=1, 50, llp[N-50+k] = logV/50 + r[k]);
    ll = llp;
    R += 1;
    if(phi, return(["failure", exp(ll[1]/N), R])));
  if(R >= abort, return(["failure", exp(ll[1]/N), R]));
  ["success", exp(ll[1]/N), R];
}

ntru_log2(x) = {if(x == 0, 0, log(x)/log(2))};

ntruNextPrime(N) = {
  my(np, e2, ord);
  if(N<100, np=100, np=N);
  while(1,
    np = nextprime(np+1);
    e2 = Mod(2, np);
    ord = znorder(e2);
    if(ord == (np-1)/2 || ord == (np-1), break;);
  );
  return(np);
}

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
  ntru_log2(1 - a/3^N);
}
addhelp(dmRejectProb,\
"dmRejectProb(N, dm): "\
"Probability that the number of 1s, -1s, or 0s is less than dm in a uniform "\
"random trinary polynomial of degree <= N");


hybridMITM(N, K, dg1, dg2) = {
    my(Nc, tot, H, awork, d1, d2, p);
    Nc = N-K;
    tot = (binomial(N,dg1)*binomial(N-dg1,dg2));
    H = vector(dg1+1, d1a,
               vector(dg2+1, d2a,
               d1 = d1a - 1;
               d2 = d2a - 1;
               p = binomial(Nc, dg1-d1) * binomial(Nc - dg1 + d1, dg2-d2) / tot;
               -(binomial(K, d1) * binomial(K-d1, d2)) * p * ntru_log2(p)));
    awork = .5*(sum(i=0, dg1, sum(j=0, dg2, H[i+1][j+1])) - ntru_log2(N));
    \\printf("HybridMITM: K: %d Work: %d\n", K, awork);
    awork;
}

minHybridMITM(N, K, dm) = {
    my(t, c);
    c = N-3*dm;
    t = hybridMITM(N, K, dm, dm+c) + .5*ntru_log2(N);
}
addhelp(minHybridMITM,\
"minHybridMITM(N, K, dm): Calculate cost of performing hybrid attack on"\
"a message with maximal information leakage via m(1), i.e. m(1) = N-3dm");


hybridHermiteRequirement(N,q,L,K) = {
  my(ld);
  ld = (N - L)*ntru_log2(q);
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
  return(3 * sqrt((4*d1*d2 + 2*d3)*pm + (4*d1*d2 + 2*d3)*pg));
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
  my(h);
  h = sqrt(q)^(1/(2*N));
  binsearchLT((x)->(-(simulate(2*N, x, h)[2])), -h, 60, N);
}

bkzCost(dim, bs, iter) = {
  my(logNodes1, logNodes2);
  \\ Quad fit to published BKZ-2.0 paper table 3 row 1
  logNodes1 = 0.00405892*bs^2 - 0.337913*bs + 34.9018;
  \\ Quad fit to Full BKZ-2.0 paper table 4
  logNodes2 = 0.000784314*bs^2 + 0.366078*bs - 6.125;
  [round(logNodes1 + ntru_log2(dim*iter) + 7), round(logNodes2 + ntru_log2(dim*iter) + 7)];
}

cn11est(dim, hreq) = {
  my(bs, iter, cost);

  bs = binsearchLT((x)->(-simulate(dim, x, hreq)[2]), -hreq, 60, dim) + 1;
  iter = simulate(dim, bs, hreq)[3];
  cost = bkzCost(dim, bs, iter);
  [iter, bs, cost[1], cost[2]];
}

getDs(N) = {
  my(d1,d2,d3);
  d1 = d2 = d3 = ceil(vecmax(real(polroots(2*x^2 + x - N/3))));
  d2 = ceil((N/3 - d1)/(2*d1));
  d3 = max(ceil((d1 + 1)/2), ceil(N/3 - 2*d1*d2));
  [d1,d2,d3]
}

/* default is estimate 2 */
genParams(N, estimate=2, verbose=0) = {
  my(lambda, directMITM, dm, d1, d2, d3, dg, sig,\
  q, q2, decFail, decFail2, out, mRej, K);

  /* Standard choices for dg, d1, d2, and d3 */
  dg = round(N/3);
  [d1,d2,d3] = getDs(N);

  /* Pick initial dm based on rejection probability below 2^-10 */
  dm = binsearchLT((x)->(dmRejectProb(N,x)), -10, 0, floor(N/3));
  mRej = dmRejectProb(N,dm);

  /* Use direct product-form MITM for upper bound on security */
  directMITM = floor(0.5 * ntru_log2(numProdForm(N,d1,d2,d3)/N));

  /* Choose q as smallest power of 2 admitting negligible
     decryption failure probability */
  sig = decFailSig((1-dm/N), (2*dg + 1)/N, d1, d2, d3);
  q = binsearchLT((x)->(-ntru_log2(N*erfc(x/(sig*sqrt(2))))), directMITM, 2^4, 2^16);
  q = 2^ceil(ntru_log2(2*q + 2));
  decFail = round(ntru_log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));

  /* Kh is smallest K that can be prepared via lattice reduction in O(2^\lambda) time. */
  [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm,estimate);

  /* Redo search for q with new security estimate */
  q2 = binsearchLT((x)->(-ntru_log2(N*erfc(x/(sig*sqrt(2))))), lambda, 2^4, 2^16);
  q2 = 2^ceil(ntru_log2(2*q2 + 2));
  decFail2 = round(ntru_log2(N*erfc(((q2-2)/2)/(sig*sqrt(2)))));
  if(q2 != q, /* If we can lower q, rederive security */
    q = q2;
    decFail = decFail2;
    [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm,estimate));

  /* Update security estimate. Either keep direct MITM cost or
   * choose larger of the two costs from hybrid attack (which should
   * be roughly equal anyway). */
  \\lambda = min(lambda, max(LM,LL));

  if(verbose,
    checkParams([N, q, d1, d2, d3, dg, dm, lambda, K], estimate));

  out = [N, q, d1, d2, d3, dg, dm, K, lambda, directMITM, decFail];
  printf("[N, q, d1, d2, d3, dg, dm, K, cost, directMITM, decFail]\n%s\n", out);
}

optimalK(N,q,d1,d2,d3,dg,dm,estimate=2) = {
  my(lambda, Kb, Llr, Lmitm, Lmsg, high, low, diff);
  lambda = floor(0.5 * ntru_log2(numProdForm(N,d1,d2,d3)/N));

  low = 0;
  high = N;
  Kb = ceil((low + high)/2);
  Llr = cn11est(2*N - Llr - Kb, hybridHermiteRequirement(N, q, Llr, Kb))[2+estimate];
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
    Llr = cn11est(2*N - Llr - Kb, hybridHermiteRequirement(N, q, Llr, Kb))[2+estimate];
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
  -ntru_log2(RHS14 / (exp(N*JT2)*sqrt(N)));
}
addhelp(probKUniq, \
"probKUniq(N,K,M): Probability that a set of M integers uniformly " \
"distributed in [1,N] contains K unique values. Expressed as -ntru_log2(prob). "\
"Estimate from: Dupuis, Zhang, Whiting, \"Refined Large Deviation Asymptotics "\
"for the Classical Occupancy Problem.\" 2006");


logBinomialCDF(k, n, p) = {
  ntru_log2(sum(i=0,k,binomial(n,i)*(p)^i * (1-p)^(n-i)));
}
addhelp(logBinomialCDF, \
"logBinomialCDF(k,n,p): ntru_log2(Pr(X <= k)), X binomial with parameters n and p.")

global(grid);
altOccupancy(C,N,n,L,d) = {
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
  my(cm, minSamp, err, minRand);
  minSamp = binsearchLT((x)->probKUniq(N, need, x), lambda, need+1, 10*need);
  cm = 2^c - lift(Mod(2^c, N));
  err = 1.0 * cm/2^c;
  minRand = binsearchLT((s)->(-logBinomialCDF(minSamp, s, err)), lambda, minSamp, floor(N*log(N)));

  /* Use upper bound to initialize memoization array for exact computation */
  altOccupancy(2^c, N, 2^c\N, minRand, need);
  minRand = need;
  while(-ntru_log2(sum(i=0, need-1, altOccupancy(2^c, N, 2^c\N, minRand, i))) < lambda, minRand+=1);

  minRand;
}

minMGFCalls(N, hashLen, lambda) = {
  my(n5, minMGF);
  n5 = ceil(N/5);
  minMGF = binsearchLT((x)->(-logBinomialCDF(n5, x, 243/256)), lambda, n5, 2*n5);
  minMGF = ceil(minMGF/hashLen);
}


formatParams(genoutput) = {
  my(lambda, c, cm, lLen, secOct, mLen, hashLen, minCalls, minRand, minMGF,\
  N, q, d1, d2, d3, dg, dm, K, need);
  [N, q, d1, d2, d3, dg, dm, lambda, K] = genoutput;
  c = ceil(ntru_log2(N));
  for(cs=c, 13, if(lift(Mod(2^cs, N))/N < lift(Mod(2^c, N))/N, c = cs));
  cm = 2^c - lift(Mod(2^c, N));

  /* Upper bound on security */
  \\lambda = round(0.5 * ntru_log2(numProdForm(N,d1,d2,d3)/N));

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
"    },\n", N,N,ceil(ntru_log2(N)), N, secOct, min(32, 2*secOct), q, ceil(ntru_log2(q)), \
d1, d2, d3, dg, mLen, dm, cm, c, minCalls, minMGF);
}

checkParams(genoutput, estimate=2) = {
  my(goodNQ, hB, K, bCombSec, bMSec, sig, decFail, mRej, \
  N, q, d1, d2, d3, dg, dm, lambda, \
  pAbove2, directMITM, preB);
  [N, q, d1, d2, d3, dg, dm, lambda, K] = genoutput;

  goodNQ = if(isprime(N),
              pAbove2 = (N-1)/znorder(Mod(2,N));
              if(isprimepower(q) && Mod(q,2) == 0,
              if(pAbove2 < 3, "Yes.",
                Strprintf("No. %d primes above 2.", pAbove2)),
                "No. q is not a power of 2."),
                "No. Composite N.");

  /* Upper bound on security */
  directMITM = floor(0.5 * ntru_log2(numProdForm(N,d1,d2,d3)/N));

  hB = hybridHermiteRequirement(N,q,lambda,K);

  bCombSec = hybridMITM(N, K, dg+1, dg);
  bMSec = minHybridMITM(N, K, dm);

  preB = cn11est(2*N - lambda - K, hB);

  sig = decFailSig((1-dm/N), 2/3, d1, d2, d3);
  decFail = round(ntru_log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));
  mRej = dmRejectProb(N,dm);

  printf("[N, q, d1, d2, d3, dg, dm] = %s\n", [N,q,d1,d2,d3,dg,dm]);
  printf("Valid N and q? %s\n", goodNQ);
  printf("Decryption failure prob. = %.1f\n", decFail);
  printf("Message rejection prob. = %.1f\n\n", mRej);
  printf("Security Estimates (using estimate %d)\n\n", estimate);

  printf("Direct MITM search cost = %d\n", directMITM);

  printf("Hybrid attack\n\tSearch K = %d coordinates\n\tMust reach root Hermite factor = %.4f\n", K, hB);
  printf("\tCN11 estimates %d rounds of BKZ-%d. Total cost = %d\n", preB[1], preB[2], preB[2+estimate]);
  printf("\tHybrid MITM cost = %d\n\tHybrid MITM cost [msg with max m(1)] = %d\n\n", bCombSec, bMSec);
}

/* Algorithm 4. per ePrint 2015/708 paper */
algo4(level=128,verbose=0) = {
  my(P,r1,r2);
  P = 100;
  while(1,
    P = ntruNextPrime(P);
    r1=genParamsInternal(level,P,1,verbose);
    r2=genParamsInternal(level,P,2,verbose);
    if(r1==0||r2==0,,
      printf("[N, q, d1, d2, d3, dg, dm]\n%s\n", r1);
      break;
    );
  );
}

/* default is estimate 2 */
genParamsInternal(level = 128, N, estimate=2, verbose=0) = {
  my(lambda, directMITM, dm, d1, d2, d3, dg, sig,\
  q, q2, decFail, decFail2, out, mRej, K);

  /* Standard choices for dg, d1, d2, and d3 */
  dg = round(N/3);
  [d1,d2,d3] = getDs(N);

  /* Pick initial dm based on rejection probability below 2^-10 */
  dm = binsearchLT((x)->(dmRejectProb(N,x)), -10, 0, floor(N/3));
  mRej = dmRejectProb(N,dm);

  /* Use direct product-form MITM for upper bound on security */
  directMITM = floor(0.5 * ntru_log2(numProdForm(N,d1,d2,d3)/N));
  if(directMITM < level,return(0));

  /* Choose q as smallest power of 2 admitting negligible
     decryption failure probability */
  sig = decFailSig((1-dm/N), (2*dg + 1)/N, d1, d2, d3);
  q = binsearchLT((x)->(-ntru_log2(N*erfc(x/(sig*sqrt(2))))), directMITM, 2^4, 2^16);
  q = 2^ceil(ntru_log2(2*q + 2));
  decFail = round(ntru_log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));

  /* {Estimate security} */
  /* Kh is smallest K that can be prepared via lattice reduction in O(2^\lambda) time. */
  [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm,estimate);
  if(level > min(directMITM,lambda), return(0));

  /* Redo search for q with new security estimate */
  q2 = binsearchLT((x)->(-ntru_log2(N*erfc(x/(sig*sqrt(2))))), lambda, 2^4, 2^16);
  q2 = 2^ceil(ntru_log2(2*q2 + 2));
  decFail2 = round(ntru_log2(N*erfc(((q2-2)/2)/(sig*sqrt(2)))));
  if(q2 != q, /* If we can lower q, rederive security */
    q = q2;
    decFail = decFail2;
    [lambda, K] = optimalK(N,q,d1,d2,d3,dg,dm,estimate));

  if(level > min(directMITM,lambda), return(0));

  if(verbose,
    checkParams([N, q, d1, d2, d3, dg, dm, lambda, K], estimate));

  out = [N, q, d1, d2, d3, dg, dm];
  out
}
