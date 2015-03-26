binsearchLT(f, val, low, high) = {
  my(mp);
\\  printf("%s, %d, %d\n", val, low, high);
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
}

minHybridMITM(N, K, dm) = {
    my(est,tot,t);
    c = N-3*dm;
    t = hybridMITM(N, K, dm, dm+c) + .5*log2(N);
}


hermiteRequirement(N,q,L,K) = {
  ld = (N - L)*log2(q);
  ld /= (4*N^2 - 4*N*(K+L) + (K^2 + 2*K*L + L^2));
  ld -= 1/(2*N - (K+L));
  2^ld
}
addhelp(hermiteRequirement,\
"hermiteRequirement(N,q,L,K): "\
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
}
addhelp(deltaStar, \
"deltaStar(k): Conjectured root Hermite factor reachable by BKZ "\
"with 2^k operations");


decFail(N, q, d1, d2, d3) = {
  my(B,sig,Ds,e1,e2);
  B = (q-2)/6;
  sig1 = sqrt(16*d1*d2/(3*N^2) + 8*d3/(3*N^3));
  sig2 = sqrt(2/3 * (4*d1*d2 + 2*d3));
  e1 = log2(N*erfc(B/(sig*sqrt(2))));
  e2 = log2(N*erfc(B/(sig*sqrt(2))));
  [e1, e2];
}

decFailSig(pm, pg, d1,d2,d3) = {
  sig = 3 * sqrt((4*d1*d2 + 2*d3)*pm + (4*d1*d2 + 2*d3)*pg);
}

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
  \\ Quad fit to Full BKZ-2.0 paper table 4
  \\logNodes = 0.000784314*bs^2 + 0.366078*bs - 6.125;
  \\ Quad fit to published BKZ-2.0 paper table 3 row 1
  logNodes = 0.00405892*bs^2 - 0.337913*bs + 34.9018;
  round(logNodes + log2(dim*iter) + 7);
}

bkzCost2(dim, bs, iter) = {
  my(logNodes);
  \\ Quad fit to Full BKZ2.0 paper table 4
  logNodes = 0.000784314*bs^2 + 0.366078*bs - 6.125;
  \\ Quad fit to published BKZ2.0 paper table 3 row 1
  \\logNodes = 0.00405892*bs^2 - 0.337913*bs + 34.9018;
  round(logNodes + log2(dim*iter) + 7);
}

cn11est(dim, hreq) = {
  my(bs, iter, logNodes);
  printf("\n----%d %f\n",dim, hreq);

  bs = binsearchLT((x)->(-simulate(dim, x, hreq)[2]), -hreq, 60, dim) + 1;
  iter = simulate(dim, bs, hreq)[3];

  [iter, bs, bkzCost(dim, bs, iter)];
}

genParams(N, lambda, verbose=0) = {
  my(dm, d1, d2, d3, dg, sig, q, Kh, Kc, dmSec);

  /* Standard choices for dg, d1, d2, and d3 */
  dg = round(N/3);
  d1 = d2 = d3 = ceil(vecmax(real(polroots(2*x^2 + x - N/3))));
  while((2*d1*(d2-1) + d3) >= round(N/3), d2 -= 1);
  while((2*d1*d2 + (d3-1)) >= round(N/3) && (d3-1) > (d1/2), d3 -= 1);

  /* Increment d3 until product form search space is sufficiently large */
  L = round(0.5 * log2(numProdForm(N,d1,d2,d3)/N));
  while(L < lambda, d3 += 1; L = round(0.5 * log2(numProdForm(N,d1,d2,d3)/N)));

  /* Pick initial dm based on rejection probability below 2^-10 */
  dm = binsearchLT((x)->(dmRejectProb(N,x)), -10, 0, floor(N/3));
  mRej = dmRejectProb(N,dm);

  /* Choose q as smallest power of 2 admitting negligible
     decryption failure probability */
  sig = decFailSig((1-dm/N), 2/3, d1, d2, d3);
  q = binsearchLT((x)->(-log2(N*erfc(x/(sig*sqrt(2))))), lambda, 2^4, 2^16);
  q = 2^ceil(log2(2*q + 2));
  decFail = round(log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));

  /* Kh is largest K that can be prepared via lattice reduction */
  Kh = binsearchLT((x)->(hermiteRequirement(N,q,lambda,x)), deltaStar(lambda), 0, N);
  /* Kc is smallest K that cannot be searched via meet-in-the-middle technique */
  Kc = binsearchLT((x)->(hybridMITM(N, x, dg, dg)), lambda, 0, Kh) + 1;

  if(verbose,
    checkParams([lambda, N, q, d1, d2, d3, dg, dm, Kc, Kh]));

  [lambda, N, q, d1, d2, d3, dg, dm, Kc, Kh]
}

blahblah(gen) = {
  [lambda, N, q, d1, d2, d3, dg, dm, Kc] = gen;
  [N, q, d1, d2, d3, dg, dm, Kc]
}



/* Extra parameters needed for reference implementation */

probKUniq(N, K, M) = {
  my(Z, T, R, S, RHS13, JT2);
  Z = 1.0*(K/N);
  T = 1.0*(M/N);

  /* Equation 3 of reference */
  R = solve(X=0.01,(1/Z)-(1e-9),-1/X*log(1-X*Z)-T);
  /* Equation 12 */
  S = sqrt(Z/(1-R*Z) - T);
  /* RHS of Equation 13 */
  RHS13 = (2*Pi*S^2)^(-1/2) * (R/(R-1)) * sqrt((1 - R*Z)/(1-Z));
  /* Equation 2 */
  JT2 = (T - Z)*log(R) + (1-Z)*log(1-Z) - (1-R*Z)/(R)*log(1-R*Z);
  /* Probability from equation 13 */
  -log2(RHS13 / (exp(N*JT2)*sqrt(N)));
}
addhelp(probKUniq, \
"probKUniq(N,K,M): Probability that a set of M integers uniformly " \
"distributed in [1,N] contains K unique values. Expressed as -log2(prob). "\
"Estimate from: Dupuis, Zhang, Whiting, \"Refined Large Deviation Asymptotics "\
"for the Classical Occupancy Problem.\" 2006");


logBinomialCDF(k, n, p) = {
  log2(sum(i=1,k,binomial(n,i)*(p)^i * (1-p)^(n-i)));
}
addhelp(logBinomialCDF, \
"logBinomialCDF(k,n,p): log2(Pr(X <= k)), X binomial with parameters n and p.")


formatParams(genoutput) = {
  my(c, cs, cm, lLen, secOct, mLen, hashLen, minSamp, minRand, err);
  [lambda, N, q, d1, d2, d3, dg, dm, K] = genoutput;
  c = ceil(log2(N));
  for(cs=c, 13, if(lift(Mod(2^cs, N))/N < lift(Mod(2^c, N))/N, c = cs));
  cm = 2^c - lift(Mod(2^c, N));
  secOct = floor(lambda/8);
  /* max message length (bytes) using 3 bit to 2 trit encoding
     assumes mLen will be < 256, assumes b is secOct bytes */
  mLen = floor(((N-1)/2)*3/8) - secOct;
  lLen = 1 + floor(log(mLen)/log(256));
  mLen -= lLen;

  hashLen = 20; if(secOct > 20, hashLen = 32);

  /* Calculate minIGFCalls */
  need = 2*(d1+d2+d3);
  minSamp = binsearchLT((x)->probKUniq(N, need, x), lambda, need+1, 10*need);
  err = 1.0 * cm/2^c;
  minRand = binsearchLT((s)->(-logBinomialCDF(minSamp, s, err)), lambda, minSamp, floor(N*log(N)));
  minCalls = ceil(ceil(minRand*c/8)/hashLen);

  /* minMGF calls assuming 5 trits are extracted from 1 byte with probability 243/256 */
  n5 = ceil(N/5);
  minMGF = binsearchLT((x)->(-logBinomialCDF(n5, x, 243/256)), lambda, n5, 2*n5);
  minMGF = ceil(minMGF/hashLen);

printf( \
"    {\n"
"        NTRU_EES%dEPX,\t\t\t/* parameter-set id */\n" \
"        \"ees%depX\",\t\t\t/* human readable param set name */\n" \
"        {0xFF, 0xFF, 0xFF},\t\t/* OID */\n" \
"        0xFF,\t\t\t\t/* DER id */\n" \
"        %d,\t\t\t\t/* no. of bits in N (i.e., in an index) */\n" \
"        %d,\t\t\t\t/* N */\n" \
"        %d,\t\t\t\t/* security strength in octets */\n" \
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
"    },\n", N,N,ceil(log2(N)), N, secOct, q, ceil(log2(q)), \
d1, d2, d3, dg, mLen, dm, cm, c, minCalls, minMGF);
}

checkParams(genoutput) = {
  my(L, sig, decFail, hermite, combSec, dmSec);
  [lambda, N, q, d1, d2, d3, dg, dm, Kc, Kh] = genoutput;

  L = round(0.5 * log2(numProdForm(N,d1,d2,d3)/N));

  hMITM = hermiteRequirement(N,q,lambda,Kc);
  hLR = hermiteRequirement(N,q,lambda,Kh);

  curCombSec = hybridMITM(N, Kh, dg, dg);
  bestCombSec = hybridMITM(N, Kc, dg, dg);
  curMSec = minHybridMITM(N, Kh, dm);
  bestMSec = minHybridMITM(N, Kc, dm);

  goodNQ = if(isprime(N),
              if(isprimepower(q) && Mod(q,2) == 0,
              if(#factormod(polcyclo(N), 2, 1)[,1] < 3, "Yes.",
                "No. Primes above 2 have \"small\" inertia degree."),
                "No. q is not a power of 2."),
                "No. Composite N.");

  sig = decFailSig((1-dm/N), 2/3, d1, d2, d3);
  decFail = round(log2(N*erfc(((q-2)/2)/(sig*sqrt(2)))));
  mRej = dmRejectProb(N,dm);

  preMITM = cn11est(2*N - lambda - Kc, hMITM);
  preLR = cn11est(2*N - lambda - Kh, hLR);

  printf("[N, q, d1, d2, d3, dg, dm] = %s\n", [N,q,d1,d2,d3,dg,dm]);
  printf("Safe N and q? %s\n", goodNQ);
  printf("Direct MITM cost = %d\n", L);
  printf("Decryption failure prob. = %.1f\n", decFail);
  printf("Message rejection prob. = %.1f\n\n", mRej);
  printf("Attack characteristics.\nMITM: Assume search cost estimate is tight.\nLR: Assume lattice reduction cost estimate is tight\n\n");
  printf("MITM K = %d : requires root Hermite factor = %.4f\n", Kc, hMITM);
  printf("\tCN11 estimates %d rounds of BKZ-%d. Total cost = %d \n", preMITM[1], preMITM[2], preMITM[3]);
  printf("\tHybrid MITM cost = %d\n\tHybrid MITM cost [msg with max m(1)] = %d\n\n", bestCombSec, bestMSec);
  printf("LR K = %d : requires root Hermite factor = %.4f\n", Kh, hLR);
  printf("\tCN11 estimates %d rounds of BKZ-%d. Total cost = %d \n", preLR[1], preLR[2], preLR[3]);
  printf("\tHybrid MITM cost = %d\n\tHybrid MITM cost [msg with max m(1)] = %d\n\n", curCombSec, curMSec);

  /* Ensure lambda <= log |product form search space| */
  \\L = round(0.5 * log2(numProdForm(N, d1, d2, d3)/N));

  \\hermite = hermiteRequirement(N,q,lambda,K);
  \\combSec = hybridMITM(N, K, dg, dg);
  \\dmSec = minHybridMITM(N, K, dm);

  \\printf("PF comb, hermite, mitm-f, mitm-m, dec fail\n");
  \\printf("need: %s\n", [lambda, deltaStar(lambda), lambda, lambda, -lambda]);
  \\printf("have: %s\n", [L, hermite, combSec, dmSec, decFail]);
  \\[L >= lambda, hermite < deltaStar(lambda), combSec >= lambda, dmSec >= lambda, decFail < -lambda]
}

P = [107, 113, 131, 139, 149, 163, 173, 181, 191, 199, 211, 227, 239, 251, 263, 271, 281, 293, 307, 317, 331, 347, 359, 367, 379, 389, 401]
sec = [34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 100, 103, 106, 109, 112]

