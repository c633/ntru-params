/* Setting for conservative security case */
global(q, p);
q = 2048;
p = 3;

ntruNextPrime(N) = {
  my(np, e2, ord);
  if(N<100, np=100, np=N);
  while(1,
    np = nextprime(np+1);
    e2 = Mod(2, np);
    ord = znorder(e2);
    if(ord == (np-1)/2 || ord == (np-1), break;);
  );
  return (np);
}

ntru_log2(x) = {if(x == 0, 0, log(x)/log(2))};

trinomial(n,r1,r2) = {
    return (binomial(n,r1) * binomial(n-r1,r2));
}

/* Lattice Security */
/* Eq. (3) per paper */
latticeRunningTime(N, y2, alpha) = {
    my(w, m, c, t1, t2);
    m = 0.2;
    c = -50;
    t1 = y2 - N;
    t2 = 1 - alpha;
    w = (2*m*t1) / (t2^2) + \
        3*log( (2*t1)/t2 ) + \
        c;
    return (2^w);
}

/* MITM Security */
f(D, sig) = {
    my(t);
    t = D/(sig*sqrt(2));
    return (erfc(t) - (exp(-(t^2))-1)/(t*sqrt(Pi)));
}

calculatePs(N, df, y2, c, alpha) = {
    my(out, dg, sig, t1, t2);
    dg = round(N/3);
    sig = sqrt((2*dg + 2*(df-c))/y2);
    t1 = 1 - alpha;
    t2 = y2 - N;
    out = (1 - 2/(3*q)) ^ ((2*N-y2*(1+alpha)) / t1) * \
        prod(i=0, 2*t2/t1, \
            (1 - f(q^((2*alpha*t2+i*t1^2) / (2*t2)), sig)) \
        );
    return (out);
}

calculateNu(N, df, y2, c, alpha) = {
    my(Nu0, c_2, ps);
    ps = calculatePs(N, df, y2, c, alpha);
    c_2 = c / 2;
    Nu0 = trinomial(2*N-y2,c_2,c_2) / trinomial(c,c_2,c_2);
    \\ Nu = trinomial(2*N-y2,c/2,c/2) * ((trinomial(c,c/2,c/2)*ps)^(-1/2));
    return (Nu0 / ps);
}

calculateT(y2) = {
    return (y2/(2^1.06));
}

calculatePSplit(N, df, y2, c) = {
    my(pSplit, ppSplit);
    pSplit = trinomial(y2-N,df-c,df-c)*trinomial(2*N-y2,c,c)/trinomial(N,df,df);
    ppSplit = 1 - (1-pSplit)^N;
    return (ppSplit);
}

MITMRunningTime(N, df, y2, c, alpha) = {
    my(t, Nu, pSplit);
    pSplit = calculatePSplit(N, df, y2, c);
    if(pSplit == 0, return (0));
    t = calculateT(y2);
    Nu = calculateNu(N, df, y2, c, alpha);
    return (t * Nu / pSplit);
}

findRoot(N, df, y2, c) = {
    my(alpha, a, b, step);
    step = 0.1;
    /* Should f(0) > 0? */
    if (latticeRunningTime(N, y2, 0) - MITMRunningTime(N, df, y2, c, 0) < 0,
        printf("Something went wrong. Expect f(0) > 0 with %s\n", [N, df, y2, c]);
        exit(1);
    );
    for(range = 0, 1,
        range += step;
        if(latticeRunningTime(N, y2, range) - MITMRunningTime(N, df, y2, c, range) < 0,
            a = range - step;
            b = range;
            break;
        );
    );
    alpha = solve(root = a, b, latticeRunningTime(N, y2, root) - MITMRunningTime(N, df, y2, c, root));
    return (alpha);
}

/* Hybrid Security */
hybridSecurityEstimate(N, df) = {
    my(alpha, w, \
        wStar, y2Star, cStar, alphaStar);
    wStar = -1; /* w* is undefined */
    for(y2 = N + 1, 2*N, /* the equations are undefined if y2 = N */
        for(c = 0, df,
            alpha = findRoot(N, df, y2, c);
            w = max(latticeRunningTime(N, y2, alpha), MITMRunningTime(N, df, y2, c, alpha));
            if(wStar == -1 || w < wStar, [wStar, y2Star, cStar, alphaStar] = [w, y2, c, alpha]);
        );
    );
    return ([wStar, y2Star, cStar, alphaStar]);
}

decryptionFailureProb(N, df) = {
    return (N*erfc(sqrt(3)*(q-2) / (24*sqrt(d))));
    \\ my(sig, dr);
    \\ dr = df;
    \\ sig = sqrt(4*(dr+df) / (3*N));
    \\ return (N*erfc(((q-2)/(2*p)) / (sig*sqrt(2*N))));
}

/* line 4 - 14 per paper */
genParamsInternal(N, k=128) = {
    my(df, k1, k2, NStar, dfStar);
    df = round(N/3);
    NStar = 0;

    while(1,
        k1 = hybridSecurityEstimate(N, df);
        k2 = ntru_log2(decryptionFailureProb(N, df));
        if(k1 >= k && k2 <= -k, [NStar, dfStar] = [N, df]);
        df = df - 1;
        if(NStar > 0 || df < 1, break);
    );

    return ([NStar, dfStar]);
}

costSpace(N, df) = {
    return (N*ntru_log2(q));
}

costSpeed(N, df) = {
    return (N*df);
}

costTradeOff(N, df) = {
    return ((costSpace(N,df)^2)*costSpeed(N,df));
}

/* default security level is 128 */
genParams(k = 128) = {
    /**
     * NStar: the first acceptable prime P which can achieve the required security
     **/
    my(NStar, dfStar, cStar, N, df, k1, c);

    NStar = 0;
    while(NStar == 0,
        N = ntruNextPrime(400);
        [NStar, dfStar] = genParamsInternal(N, k);
    );

    /* TODO: Run the algorithm with the first round then terminate */
    return;
    cStar = costSpace(NStar, dfStar);

    while( 1 /* an increase in N can potentially lower the cost */,
        df = dfStar;
        while(df >=0,
            k1 = hybridSecurityEstimate(N, df);
            c = costSpace(N, df);
            if(k1 >= k && c <= cStar, [cStar,NStar,dfStar] = [c,N,df]);
            df = df - 1;
        );
        N = ntruNextPrime(Nstar);
    );
    printf("[Pi*, df*] = [%d, %d]\n", NStar, dfStar);
}