### Usage
- Run the GP scripts
```
$ gp
[snip]
? \r ntru_params.gp
[snip]
? #
   timer = 1 (on)
? algo4(112)
[N, q, d1, d2, d3, dg, dm]
[401, 2048, 8, 8, 6, 134, 102]
time = 12min, 5,857 ms.
```

- Run the compiled code
```
$ gp2c-run -g ntru_params.gp
[snip]
? default(parisizemax,4*10^8)
  ***   Warning: new maximum stack size = 400003072 (381.473 Mbytes).
? #
   timer = 1 (on)
? algo4(112)
  *** algo4: Warning: increasing stack size to 16000000.
  *** algo4: Warning: increasing stack size to 32000000.
[N, q, d1, d2, d3, dg, dm]
[401, 2048, 8, 8, 6, 134, 102]
time = 6min, 36,532 ms.
```

### Timing
* GP language

| Security Level    | Parameter Sets | Time          |
| :---------------: | :------------: | :-----------: |
| 112 | [401, 2048, 8, 8, 6, 134, 102] | 12min, 5,857 ms. |
| 128 | [443, 2048, 9, 8, 5, 148, 115] | 13min, 26,144 ms. |
| 192 | - | - |
| 256 | - | - |

* C compiled code

| Security Level    | Parameter Sets | Time          |
| :---------------: | :------------: | :-----------: |
| 112 | [401, 2048, 8, 8, 6, 134, 102] | 6min, 36,532 ms. |
| 128 | [443, 2048, 9, 8, 5, 148, 115] | 7min, 29,363 ms. |
| 192 | [677, 4096, 11, 10, 6, 226, 185] | 4min, 37,590 ms. |
| 256 | - | - |


### Sample output
```
[N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM]
[401, 2048, 8, 8, 6, 134, 102, 116, 154, -217, 145]
    {
        NTRU_EES401EPX,			/* parameter-set id */
        "ees401epX",			/* human readable param set name */
        {0xFF, 0xFF, 0xFF},		/* OID */
        0xFF,				/* DER id */
        9,				/* no. of bits in N (i.e., in an index) */
        401,				/* N */
        14,				/* security strength in octets */
        28,				/* no. of octets for random string b */
        2048,				/* q */
        11,				/* no. of bits in q (i.e., in a coeff) */
        TRUE,				/* product form */
        8 + (8 << 8) + (6 << 16),	/* df, dr */
        134,				/* dg */
        46,				/* maxMsgLenBytes */
        102,				/* dm0 */
        2005,				/* 2^c - (2^c mod N) */
        11,				/* c */
        1,				/* lLen */
        5,				/* min. no. of hash calls for IGF-2 */
        5,				/* min. no. of hash calls for MGF-TP-1 */
        NTRU_CRYPTO_HASH_ALGID_SHA256,	/* hash function for MGF-TP-1, HMAC-DRBG, etc. */
    },
[N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM]
[439, 2048, 9, 8, 5, 146, 113, 133, 174, -195, 147]
    {
        NTRU_EES439EPX,			/* parameter-set id */
        "ees439epX",			/* human readable param set name */
        {0xFF, 0xFF, 0xFF},		/* OID */
        0xFF,				/* DER id */
        9,				/* no. of bits in N (i.e., in an index) */
        439,				/* N */
        16,				/* security strength in octets */
        32,				/* no. of octets for random string b */
        2048,				/* q */
        11,				/* no. of bits in q (i.e., in a coeff) */
        TRUE,				/* product form */
        9 + (8 << 8) + (5 << 16),	/* df, dr */
        146,				/* dg */
        49,				/* maxMsgLenBytes */
        113,				/* dm0 */
        439,				/* 2^c - (2^c mod N) */
        9,				/* c */
        1,				/* lLen */
        6,				/* min. no. of hash calls for IGF-2 */
        5,				/* min. no. of hash calls for MGF-TP-1 */
        NTRU_CRYPTO_HASH_ALGID_SHA256,	/* hash function for MGF-TP-1, HMAC-DRBG, etc. */
    },
[N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM]
[593, 4096, 10, 10, 6, 198, 159, 181, 243, -578, 181]
    {
        NTRU_EES593EPX,			/* parameter-set id */
        "ees593epX",			/* human readable param set name */
        {0xFF, 0xFF, 0xFF},		/* OID */
        0xFF,				/* DER id */
        10,				/* no. of bits in N (i.e., in an index) */
        593,				/* N */
        22,				/* security strength in octets */
        32,				/* no. of octets for random string b */
        4096,				/* q */
        12,				/* no. of bits in q (i.e., in a coeff) */
        TRUE,				/* product form */
        10 + (10 << 8) + (6 << 16),	/* df, dr */
        198,				/* dg */
        78,				/* maxMsgLenBytes */
        159,				/* dm0 */
        1779,				/* 2^c - (2^c mod N) */
        11,				/* c */
        1,				/* lLen */
        8,				/* min. no. of hash calls for IGF-2 */
        7,				/* min. no. of hash calls for MGF-TP-1 */
        NTRU_CRYPTO_HASH_ALGID_SHA256,	/* hash function for MGF-TP-1, HMAC-DRBG, etc. */
    },
[N, q, d1, d2, d3, dg, dm, lambda, K, decFail, directMITM]
[743, 4096, 11, 11, 6, 248, 205, 201, 332, -482, 201]
    {
        NTRU_EES743EPX,			/* parameter-set id */
        "ees743epX",			/* human readable param set name */
        {0xFF, 0xFF, 0xFF},		/* OID */
        0xFF,				/* DER id */
        10,				/* no. of bits in N (i.e., in an index) */
        743,				/* N */
        25,				/* security strength in octets */
        32,				/* no. of octets for random string b */
        4096,				/* q */
        12,				/* no. of bits in q (i.e., in a coeff) */
        TRUE,				/* product form */
        11 + (11 << 8) + (6 << 16),	/* df, dr */
        248,				/* dg */
        106,				/* maxMsgLenBytes */
        205,				/* dm0 */
        8173,				/* 2^c - (2^c mod N) */
        13,				/* c */
        1,				/* lLen */
        7,				/* min. no. of hash calls for IGF-2 */
        8,				/* min. no. of hash calls for MGF-TP-1 */
        NTRU_CRYPTO_HASH_ALGID_SHA256,	/* hash function for MGF-TP-1, HMAC-DRBG, etc. */
    },
```