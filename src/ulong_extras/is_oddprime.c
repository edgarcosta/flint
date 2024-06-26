/*
    Copyright (C) 2009 Tom Boothby
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

#if FLINT64

ulong FLINT_ODD_PRIME_LOOKUP[] =
{
    0x816d129a64b4cb6eUL, UWORD(0x2196820d864a4c32), UWORD(0xa48961205a0434c9),
    UWORD(0x4a2882d129861144), UWORD(0x834992132424030),  0x148a48844225064bUL,
    UWORD(0xb40b4086c304205),  UWORD(0x65048928125108a0), UWORD(0x80124496804c3098),
    UWORD(0xc02104c941124221), UWORD(0x804490000982d32),  UWORD(0x220825b082689681),
    UWORD(0x9004265940a28948), UWORD(0x6900924430434006), UWORD(0x12410da408088210),
    UWORD(0x86122d22400c060),  UWORD(0x110d301821b0484),  UWORD(0x14916022c044a002),
    0x92094d204a6400cUL,  UWORD(0x4ca2100800522094), UWORD(0xa48b081051018200),
    UWORD(0x34c108144309a25),  UWORD(0x2084490880522502), UWORD(0x241140a218003250),
    UWORD(0xa41a00101840128),  UWORD(0x2926000836004512), UWORD(0x10100480c0618283),
    0xc20c26584822006dUL, UWORD(0x4520582024894810), UWORD(0x10c0250219002488),
    UWORD(0x802832ca01140868), UWORD(0x60901300264b0400)
};

#else

ulong FLINT_ODD_PRIME_LOOKUP[] =
{
    0x64b4cb6eUL, 0x816d129aUL, UWORD(0x864a4c32), 0x2196820dUL,
    UWORD(0x5a0434c9), UWORD(0xa4896120), UWORD(0x29861144), UWORD(0x4a2882d1),
    UWORD(0x32424030), UWORD(0x8349921),  0x4225064bUL, UWORD(0x148a4884),
    UWORD(0x6c304205), UWORD(0xb40b408),  UWORD(0x125108a0), UWORD(0x65048928),
    UWORD(0x804c3098), UWORD(0x80124496), UWORD(0x41124221), UWORD(0xc02104c9),
    UWORD(0x982d32),   UWORD(0x8044900),  UWORD(0x82689681), UWORD(0x220825b0),
    UWORD(0x40a28948), UWORD(0x90042659), UWORD(0x30434006), UWORD(0x69009244),
    UWORD(0x8088210),  UWORD(0x12410da4), UWORD(0x2400c060), UWORD(0x86122d2),
    UWORD(0x821b0484), UWORD(0x110d301),  UWORD(0xc044a002), UWORD(0x14916022),
    0x4a6400cUL,  UWORD(0x92094d2),  UWORD(0x522094),   UWORD(0x4ca21008),
    UWORD(0x51018200), UWORD(0xa48b0810), UWORD(0x44309a25), UWORD(0x34c1081),
    UWORD(0x80522502), UWORD(0x20844908), UWORD(0x18003250), UWORD(0x241140a2),
    UWORD(0x1840128),  UWORD(0xa41a001),  UWORD(0x36004512), UWORD(0x29260008),
    UWORD(0xc0618283), UWORD(0x10100480), 0x4822006dUL, UWORD(0xc20c2658),
    UWORD(0x24894810), UWORD(0x45205820), UWORD(0x19002488), UWORD(0x10c02502),
    UWORD(0x1140868),  0x802832caUL, UWORD(0x264b0400), UWORD(0x60901300)
};

#endif

int n_is_oddprime_small(ulong n)
{
    ulong q = n / 2;
    ulong x = (q & (FLINT_BITS - UWORD(1)));
    return (FLINT_ODD_PRIME_LOOKUP[q / FLINT_BITS] & (UWORD(1) << x)) >> x;
}

int n_is_oddprime_binary(ulong n)
{
    ulong diff, prime_lo, prime_hi;
    const ulong * primes;

    n_prime_pi_bounds(&prime_lo, &prime_hi, n);
    primes = n_primes_arr_readonly(prime_hi + 1);

    prime_hi--; /* convert to indices of primes in table */
    prime_lo--;

    if (n == primes[prime_hi]) return 1;
    if (n > primes[prime_hi]) return 0;

    diff = (prime_hi - prime_lo + 1) / 2;

    while (1)
    {
        ulong diff2;
        if (primes[prime_lo + diff] <= n) prime_lo += diff;
        if (diff <= UWORD(1)) break;
        diff  = (diff + 1)/2;
        diff2 = (prime_hi - prime_lo + 1)/2;
        if (diff > diff2) diff = diff2;
    }

    return (n == primes[prime_lo]);
}
